use clap::{Arg, ArgAction, Command};
use needletail::{parse_fastx_file, Sequence};
use needletail::kmer::Kmers;
use needletail::sequence::canonical;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crossbeam_channel::bounded;
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::fs;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::thread;
use std::path::Path;
use hyperminhash::Sketch;
use serde_json::{to_writer_pretty};
use serde_json::json;

mod hypermash;
use hypermash::Hypermash; 
use ultraloglog::{Estimator, MaximumLikelihoodEstimator, OptimalFGRAEstimator, UltraLogLog};
mod ultraloglog_utils;
use crate::ultraloglog_utils::{load, save, Estimation};



fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    env_logger::Builder::from_default_env().init();
    // Set up the command-line arguments
    let matches = Command::new("Genome Sketching via HyperMinhash")
        .version("0.1.0")
        .about("Fast and Memory Efficient Genome/Metagenome Sketching via HyperMinhash")
        .subcommand(
            Command::new("sketch")
            .about("Creates vectors and serializes them")
            .arg(
                Arg::new("sketch")
                .short('s')
                .long("sketch")
                .help("One file containing list of FASTA files (.gz supported). File must be UTF-8.")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("output")
                .short('o')
                .long("output")
                .help("Input a prefix/name for your output files")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("kmer_length")
                .short('k')
                .long("kmer")
                .help("Length of the kmer")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("algorithm")
                .short('a')
                .long("algorithm")
                .help("Which algorithm to use: hyperminhash (hmh) or ultraloglog (ull)")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("precision")
                .short('p')
                .long("precision")
                .help("Specifiy precision, for ull only")
                .required(false)
                .requires_if("ull", "algorithm")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("estimator")
                .short('e')
                .long("estimator")
                .help("Specify estimator (optimal or ml), for ull only")
                .required(false)
                .requires_if("ull", "algorithm")
                .action(ArgAction::Set)
            )
        )
        .subcommand(
            Command::new("distance")
            .about("Computes distance between sketches")
            .arg(
                Arg::new("reference")
                .short('r')
                .long("reference")
                .help("Prefix to search for reference genome files")
                .required(true)
                .num_args(2)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("query")
                .short('q')
                .long("query")
                .help("Prefix to search for query genome files")
                .required(true)
                .num_args(2)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("output_file")
                .short('o')
                .long("output_file")
                .help("Name of output file to write results")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set)
            )
            
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sketch", s_matches)) => {
            // organize the inputs
            let sketch_file_name = s_matches.get_one::<String>("sketch").expect("required");
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").expect("required");
            let threads: usize = *s_matches.get_one::<usize>("threads").unwrap_or(&num_cpus::get());

            let output_name = s_matches.get_one::<String>("output").expect("required");
            let alg = s_matches.get_one::<String>("algorithm").expect("required");

            if alg == "hmh" {
                // create hypermash object and sketch
                let hm = Hypermash::new(sketch_file_name.clone(), kmer_length, threads, output_name.clone());
                hm.sketch();
            }
            else if alg == "ull" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").expect("required") as u32;
                let estimator = s_matches.get_one::<String>("estimator").expect("required");
                // create ull object and sketch
                let sketch_file = File::open(sketch_file_name)?;
                let sketch_reader = BufReader::new(sketch_file);
                let files: Vec<String> = sketch_reader
                    .lines()
                    .filter_map(|line| line.ok())
                    .filter(|line| !line.trim().is_empty())
                    .collect();

                // sketch the genomes per file
                files.par_iter().for_each(|file| {
                    let hashes = load::<u64>(file);
                    let mut ull = UltraLogLog::new(precision).expect("failed to create ull");
                    let estimations: Vec<Estimation> = hashes
                        .iter()
                        .enumerate()
                        .map(|(i, num)| {
                            ull.add(*num);
                            let est = match estimator.as_str() {
                                "optimal" => OptimalFGRAEstimator.estimate(&ull),
                                "ml" => MaximumLikelihoodEstimator.estimate(&ull),
                                _ => unreachable!(),
                            };

                            Estimation((i + 1) as u64, est as u64, estimator.to_string())
                        })
                        .collect();

                    let basename = Path::new(file).file_name().unwrap().to_str().unwrap();

                    // write sketches to file
                    let filename = format!("est-p{}-{}-{}", precision, estimator, basename);
                    save(&estimations, filename.as_str(), output_name);
                });

            // write all file names to a json file
            to_writer_pretty(
                &File::create(format!("{}_names.json", output_name))?,
                &files
            )?;

            // save a json of parameters used
            let data = json!({
                "k": kmer_length.to_string(),
                "algorithm": "ull", 
                "precision": precision.to_string(),
                "estimator": estimator
            });
        
            let json_str = serde_json::to_string_pretty(&data).unwrap();
            let mut param_file = File::create(format!("{}_parameters.json", 
                output_name))?;
            param_file.write_all(json_str.as_bytes())?;

            }
            Ok(())
        }
        Some(("distance", s_matches)) => {
            let ref_prefix = s_matches.get_one::<String>("reference").expect("required");
            let query_prefix = s_matches.get_one::<String>("query").expect("required");
            
            fn find_files(prefix: &String) -> std::io::Result<HashMap<&str, String>> {
                let mut files: Vec<String> = Vec::new();
                let dir = "./"; // use curent directory

                for entry in fs::read_dir(dir)? {
                    let entry = entry?;
                    let path = entry.path();

                    // if file starts with prefix, add it to vector of possible reference/ sketch files
                    if path.is_file() {
                        if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                            if filename.starts_with(prefix) {
                                files.push(filename.to_string());
                            }
                        }
                    }
                }

                if files.len() != 3 {
                    panic!("There should be 3 files starting with {} but {} were found instead", 
                    prefix, 
                    files.len());
                }
                let mut file_map: HashMap<&str, String> = HashMap::new();
                for file in &files {
                    if file.ends_with("parameters.json") {
                        file_map.insert("params",file.clone());
                    }
                    else if file.ends_with("names.json") {
                        file_map.insert("names",file.clone());
                    }
                    else { // .bin file for sketches
                        file_map.insert("sketches",file.clone());
                    }
                }
                Ok(file_map)
            }
            
            let output_file: &String = s_matches.get_one::<String>("output_file").expect("required");
            let threads: usize = *s_matches.get_one::<usize>("threads").unwrap();

            // go through the files needed, find name file, sketch file, and param file
            let ref_files = find_files(ref_prefix).unwrap();
            let query_files = find_files(query_prefix).unwrap();
            let ref_param_file = ref_files["params"].clone();
            let query_param_file = query_files["params"].clone();
            
            let query_map: HashMap<String, String> = serde_json::from_str(&query_param_file).expect("Invalid JSON");
            let ref_map: HashMap<String, String> = serde_json::from_str(&ref_param_file).expect("Invalid JSON");
            
            if ref_map["k"] != query_map["k"] {
                panic!("Genomes were not sketched with the same k");
            }
            if ref_map["algorithm"] != query_map["algorithm"] {
                panic!("Algorithms do not match in query and sketch genomes")
            }
            if ref_map["algorithm"] == "ull" {
                if ref_map["precision"] != query_map["precision"] {
                    panic!("ull was not sketched with same precision")
                }
                if ref_map["estimator"] != query_map["estimator"] {
                    panic!("ull was not sketched with same estimator")
                }
            }
            // assign kmer length once k matches
            let kmer_length:usize = ref_map["k"].parse().expect("Not a valid usize");

            ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();

            // function to read in names of the genomes into a vector of names
            fn read_names(file_name: &str) -> std::io::Result<Vec<String>> {
                let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
                let reader = BufReader::new(file);
                let names: Vec<String> = serde_json::from_reader(reader)?;
                Ok(names)
            }

            // function to read in sketches and put them together with their names into a hashmap
            fn read_sketches<'a>(file_name: &'a str, names: &'a Vec<String>) -> std::io::Result<Vec<Sketch>> {
                let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
                let mut reader = BufReader::new(file);
                let mut sketches = Vec::new();
                for name in names {
                    let sketch = Sketch::load(&mut reader)?;
                    sketches.push(sketch);
                }
                Ok(sketches)
            }

            //create query sketch hashmap
            let query_namefile = query_files["names"].clone();
            let query_sketch_file = query_files["sketches"].clone();
            let query_names: Vec<String> = read_names(&query_namefile)
                            .expect(&format!("Error with reading from {}", query_namefile));
            let q_sketch_vec: Vec<Sketch> = read_sketches(&query_sketch_file, &query_names)
                            .expect(&format!("Error with reading from {}", query_sketch_file));
            
            let mut index = 0;
            let mut query_sketches: HashMap<&String, &Sketch> = HashMap::new();
            for sketch in &q_sketch_vec {
                query_sketches.insert(&query_names[index], sketch);
                index += 1;
            }
            
            // create reference sketch hashmap
            let ref_namefile = ref_files["names"].clone();
            let ref_sketch_file = ref_files["sketches"].clone();
            let reference_names: Vec<String> = read_names(&ref_namefile)
                            .expect(&format!("Error with reading from {}", ref_namefile));
            let r_sketch_vec = read_sketches(&ref_sketch_file, &reference_names)
                            .expect(&format!("Error with reading from {}", ref_sketch_file));
            let mut reference_sketches = HashMap::new();
            index = 0;
            for sketch in &r_sketch_vec {
                reference_sketches.insert(&reference_names[index], sketch);
                index += 1;
            }

            //Generate all pairs of query and reference sketches
            let pairs: Vec<(&String, &Sketch, &String, &Sketch)> = query_sketches
            .iter()
            .flat_map(|(q_name, q_sketch)| {
                reference_sketches.iter().map(move |(r_name, r_sketch)| {
                    (*q_name, *q_sketch, *r_name, *r_sketch)
                })
            })
            .collect();


            // Compute similarities and distances in parallel
            let results: Vec<(String, String, f64)> = pairs
                .par_iter()
                .map(|&(query_name, query_sketch, reference_name, reference_sketch)| {
                    let similarity = query_sketch.similarity(reference_sketch);

                    // Avoid division by zero and log of zero
                    let adjusted_similarity = if similarity <= 0.0 {
                        std::f64::EPSILON // Small positive number to avoid log(0)
                    } else {
                        similarity
                    };

                    // Calculate distance using the provided formula
                    let numerator = 2.0 * adjusted_similarity;
                    let denominator = 1.0 + adjusted_similarity;
                    let fraction: f64 = numerator / denominator;
                    let distance = -fraction.ln() / (kmer_length as f64);

                    (query_name.clone(), reference_name.clone(), distance)
                })
                .collect();

            // Open the output file for writing
            let mut output = File::create(output_file)?;

            // Write header line
            writeln!(output, "Query\tReference\tDistance")?;

            // Write the results with file path normalization
            for (query_name, reference_name, distance) in &results {
                // Extract file names from the paths
                let query_basename = Path::new(&query_name)
                    .file_name()
                    .and_then(|os_str| os_str.to_str())
                    .unwrap_or(&query_name);

                let reference_basename = Path::new(&reference_name)
                    .file_name()
                    .and_then(|os_str| os_str.to_str())
                    .unwrap_or(&reference_name);

                let distance = if query_basename == reference_basename {
                    0.0
                } else {
                    *distance
                };
                writeln!(output, "{}\t{}\t{:.6}", query_name, reference_name, distance)?;
            }
            println!("Distances computed.");
            Ok(())
        }
        _ => {
            Ok(())
        }
    }
}
    

