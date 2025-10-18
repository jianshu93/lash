use clap::{Arg, ArgAction, Command};
// use needletail::kmer::Kmers;
// use needletail::sequence::canonical;
use hashbrown::HashMap;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::error::Error;
//use xxhash_rust::xxh3::Xxh3Builder;
use std::fs;
use std::fs::File;
use std::io::{BufReader, Write};
use std::path::Path;
use zstd::stream::Decoder;
mod hasher;
use hasher::Xxh3Builder;
use log::info;
mod utils;
use crate::utils::{hll_distance, hll_sketch, hmh_distance, hmh_sketch, ull_sketch};
use ultraloglog::{Estimator, MaximumLikelihoodEstimator, UltraLogLog};

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    env_logger::Builder::from_default_env().init();
    // Set up the command-line arguments
    let matches = Command::new("Genome Sketching via HyperLogLog, HyperMinhash and UltraLogLog")
        .version("0.1.3")
        .about("Fast and Memory Efficient (Meta)genome Sketching via HyperLogLog, HyperMinhash and UltraLogLog")
        .subcommand(
            Command::new("sketch")
            .about("Sketches genomes and serializes them, sketches are compressed")
            .arg(
                Arg::new("file")
                .short('f')
                .long("file")
                .help("One file containing list of FASTA/FASTQ files (.gz/.bz2/.zstd supported), one per line. File must be UTF-8.")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("output")
                .short('o')
                .long("output")
                .help("Input a prefix/name for your output files")
                .required(false)
                .default_value("sketch")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("kmer_length")
                .short('k')
                .long("kmer")
                .help("Length of the kmer")
                .required(false)
                .default_value("16")
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
                .default_value("1")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("algorithm")
                .short('a')
                .long("algorithm")
                .help("Which algorithm to use: HyperMinHash (hmh), UltraLogLog (ull), or HyperLogLog (hll)")
                .required(false)
                .default_value("hmh")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("precision")
                .short('p')
                .long("precision")
                .help("Specifiy precision, for ull and hll only.")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .default_value("10")
                .action(ArgAction::Set)
            )
        )
        .subcommand(
            Command::new("dist")
            .about("Computes distance between sketches")
            .arg(
                Arg::new("query")
                .short('q')
                .long("query")
                .help("Prefix to search for query genome files")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("reference")
                .short('r')
                .long("reference")
                .help("Prefix to search for reference genome files")
                .required(true)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("output_file")
                .short('o')
                .long("output_file")
                .help("Name of output file to write results")
                .required(false)
                .default_value("dist.txt")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use")
                .required(false)
                .value_parser(clap::value_parser!(usize))
                .default_value("1")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("estimator")
                .short('e')
                .long("estimator")
                .help("Specify estimator (fgra or ml), for ull only")
                .required(false)
                .default_value("fgra")
                .action(ArgAction::Set)
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sketch", s_matches)) => {
            // organize the inputs
            let sketch_file_name = s_matches.get_one::<String>("file").expect("required");
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").expect("required");
            let threads: usize = *s_matches
                .get_one::<usize>("threads")
                .unwrap_or(&num_cpus::get());

            let output_name = s_matches.get_one::<String>("output").expect("required");
            let alg = s_matches.get_one::<String>("algorithm").expect("required");

            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            let result: Result<(), Box<dyn Error>>;
            if alg == "hmh" {
                // create hypermash object and sketch
                result = hmh_sketch(
                    sketch_file_name.clone(),
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                );
            } else if alg == "hll" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = hll_sketch(
                    precision,
                    sketch_file_name.clone(),
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                );
            } else if alg == "ull" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = ull_sketch(
                    precision,
                    sketch_file_name.clone(),
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                );
            } else {
                // input for alg is not hmh, ull, or hll
                panic!("Algorithm must be either hmh, ull, or hll");
            }
            result
        }
        Some(("dist", s_matches)) => {
            let ref_prefix = s_matches.get_one::<String>("reference").expect("required");
            let query_prefix = s_matches.get_one::<String>("query").expect("required");

            fn find_files(prefix: &String) -> std::io::Result<HashMap<&str, String>> {
                let mut files: Vec<String> = Vec::new();
                let dir = "./"; // use curent directory

                let norm_prefix = {
                    let p = Path::new(prefix)
                        .file_name()
                        .and_then(|os_str| os_str.to_str())
                        .unwrap_or(prefix);
                    // .strip_prefix("./")
                    // .unwrap_or(prefix)

                    p.strip_prefix("./").unwrap_or(p)
                };

                for entry in fs::read_dir(dir)? {
                    let entry = entry?;
                    let path = entry.path();

                    // if file starts with prefix, add it to vector of possible reference/ sketch files
                    if path.is_file() {
                        if let Some(filename) = path.file_name().and_then(|n| n.to_str()) {
                            if filename.starts_with(norm_prefix) {
                                files.push(filename.to_string());
                            }
                        }
                    }
                }

                // if files.len() != 3 {
                //     panic!("There should be 3 files starting with {} but {} were found instead",
                //     prefix,
                //     files.len());
                // }
                let mut file_map: HashMap<&str, String> = HashMap::new();
                for file in files {
                    if file.ends_with("parameters.json") {
                        file_map.insert("params", file);
                    } else if file.ends_with("files.json") {
                        file_map.insert("files", file);
                    } else if file.ends_with(".bin") {
                        // .bin file for sketches
                        file_map.insert("sketches", file);
                    }
                }
                if file_map.keys().len() != 3 {
                    panic!(
                        "There should be 3 files starting with {} but {} were found instead",
                        norm_prefix,
                        file_map.keys().len()
                    );
                }
                Ok(file_map)
            }

            let output_file: &String = s_matches
                .get_one::<String>("output_file")
                .expect("required");
            let threads: usize = *s_matches
                .get_one::<usize>("threads")
                .unwrap_or(&num_cpus::get());

            // go through the files needed, find name file, sketch file, and param file
            let ref_files = find_files(ref_prefix).unwrap();
            let query_files = find_files(query_prefix).unwrap();
            let ref_param_file = ref_files["params"].clone();
            let query_param_file = query_files["params"].clone();

            // println!("{}", ref_param_file);
            // println!("{}", query_param_file);

            // read in parameter json files into hashmaps
            let contents_ref = fs::read_to_string(ref_param_file).expect("Unable to read file");
            let contents_query = fs::read_to_string(query_param_file).expect("Unable to read file");
            let ref_map: HashMap<String, String> =
                serde_json::from_str(&contents_ref).expect("Invalid JSON");
            let query_map: HashMap<String, String> =
                serde_json::from_str(&contents_query).expect("Invalid JSON");

            // check that parameters match between ref and query genomes
            if ref_map["k"] != query_map["k"] {
                panic!("Genomes were not sketched with the same k");
            }
            if ref_map["algorithm"] != query_map["algorithm"] {
                panic!("Algorithms do not match in query and sketch genomes")
            }
            if ref_map["algorithm"] == "ull" || ref_map["algorithm"] == "hll" {
                if ref_map["precision"] != query_map["precision"] {
                    panic!("{} was not sketched with same precision btwn genomes", ref_map["algorithm"]);
                }
            }
            // assign kmer length once k matches
            let kmer_length: usize = ref_map["k"].parse().expect("Not a valid usize");

            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            // function to read in names of the genomes, outputs a vector of names
            fn read_names(file_name: &str) -> std::io::Result<Vec<String>> {
                let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
                let reader = BufReader::new(file);
                let names: Vec<String> = serde_json::from_reader(reader)?;
                Ok(names)
            }

            //create query sketch hashmap
            let query_namefile = query_files["files"].clone();
            let query_sketch_file_name = query_files["sketches"].clone();
            let query_names: Vec<String> = read_names(&query_namefile)
                .expect(&format!("Error with reading from {}", query_namefile));

            // create reference sketch hashmap
            let ref_namefile = ref_files["files"].clone();
            let ref_sketch_file_name = ref_files["sketches"].clone();
            let reference_names: Vec<String> = read_names(&ref_namefile)
                .expect(&format!("Error with reading from {}", ref_namefile));

            // let mut results: Vec<(String, String, f64)> = Vec::new();
            // --- compute all pairwise distances into a single results Vec
            let results: Vec<(String, String, f64)> = if ref_map["algorithm"] == "hmh" {
                hmh_distance(
                    reference_names,
                    ref_sketch_file_name,
                    kmer_length,
                    query_names,
                    query_sketch_file_name,
                ).unwrap()
            } else if ref_map["algorithm"] == "ull" {
                // ULL: build maps, pairs, and compute
                let ref_sketch_file = File::open(ref_sketch_file_name).expect("Failed to open file");
                let query_sketch_file = File::open(query_sketch_file_name).expect("Failed to open file");
                let estimator = s_matches
                    .get_one::<String>("estimator")
                    .cloned()
                    .unwrap_or_else(|| "fgra".to_string());

                fn create_ull_map(
                    sketch_file: File,
                    names: &Vec<String>,
                    estimator: &String,
                ) -> Result<HashMap<String, (UltraLogLog, f64), Xxh3Builder>, std::io::Error> {
                    let hasher = Xxh3Builder { seed: 93 };
                    let mut sketches = HashMap::with_hasher(hasher);
                    let reader = BufReader::new(sketch_file);
                    let mut decoder = Decoder::new(reader).expect("failed to create decompressor");
                    for file in names {
                        let ull = UltraLogLog::load(&mut decoder)?;
                        let c: f64 = match estimator.as_str() {
                            "fgra" => ull.get_distinct_count_estimate(),
                            "ml"   => MaximumLikelihoodEstimator.estimate(&ull),
                            _      => panic!("estimator needs to be either fgra or ml"),
                        };
                        sketches.insert(file.clone(), (ull, c));
                    }
                    Ok(sketches)
                }

                let ref_map = create_ull_map(ref_sketch_file, &reference_names, &estimator).unwrap();
                let query_map = create_ull_map(query_sketch_file, &query_names, &estimator).unwrap();

                let pairs: Vec<(&str, &str)> = ref_map
                    .keys()
                    .flat_map(|k1| query_map.keys().map(move |k2| (k1.as_str(), k2.as_str())))
                    .collect();

                pairs.par_iter()
                    .map(|&(ref_name, qry_name)| {
                        let a: f64 = ref_map[ref_name].1;
                        let b: f64 = query_map[qry_name].1;

                        let union_ull = UltraLogLog::merge(&ref_map[ref_name].0, &query_map[qry_name].0)
                            .expect("failed to merge sketches");
                        let union_count = union_ull.get_distinct_count_estimate();

                        info!("Union: {}, a: {}, b: {}", union_count, a, b);

                        let similarity = (a + b - union_count) / union_count;
                        let s = if similarity <= 0.0 { std::f64::EPSILON } else { similarity };
                        let distance = -( (2.0*s) / (1.0 + s) ).ln() / (kmer_length as f64);

                        (ref_name.to_string(), qry_name.to_string(), distance)
                    })
                    .collect()
            } else {
                // HLL
                hll_distance(
                    reference_names,
                    ref_sketch_file_name,
                    kmer_length,
                    query_names,
                    query_sketch_file_name,
                ).unwrap()
            };

            // --- write output ---
            let mut output = File::create(output_file)?; // don't append ".txt" again
            writeln!(output, "Query\tReference\tDistance")?;

            // results are (reference, query, distance)
            for (reference_name, query_name, distance) in &results {
                let query_basename = Path::new(query_name.as_str())
                    .file_name()
                    .and_then(|os_str| os_str.to_str())
                    .unwrap_or(query_name.as_str());

                let reference_basename = Path::new(reference_name.as_str())
                    .file_name()
                    .and_then(|os_str| os_str.to_str())
                    .unwrap_or(reference_name.as_str());

                let d = if query_basename == reference_basename { 0.0 } else { *distance };

                writeln!(output, "{}\t{}\t{:.6}", query_name, reference_name, d)?;
            }
            println!("Distances computed.");
            Ok(())
        }
        _ => Ok(()),
    }
}
