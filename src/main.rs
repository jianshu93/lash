use clap::{Arg, ArgAction, Command};
// use needletail::kmer::Kmers;
// use needletail::sequence::canonical;
use hashbrown::HashMap;
use std::error::Error;
//use xxhash_rust::xxh3::Xxh3Builder;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use streaming_algorithms::HyperLogLog;
use hyperminhash::Sketch;
use ultraloglog::UltraLogLog;
mod hasher;
use serde_json::json;
mod utils;
use crate::utils::{hll_distance, hmh_distance, ull_distance, sketch_files};
use num_traits::Float;
use std::sync::{Arc, Mutex};

fn main() -> Result<(), Box<dyn Error>> {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    env_logger::Builder::from_default_env().init();
    // Set up the command-line arguments
    let matches = Command::new("Genome Sketching via HyperLogLog, HyperMinhash and UltraLogLog")
        .version("0.1.4")
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
                .help("Number of threads to use, default to all logical cores")
                .required(false)
                .value_parser(clap::value_parser!(usize))
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
            .arg(
                Arg::new("seed")
                .short('s')
                .long("seed")
                .help("Random seed")
                .required(false)
                .value_parser(clap::value_parser!(u64))
                .default_value("42")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("aa")
                .long("aa")
                .help("Amino acid sketching")
                .required(false)
                .action(clap::ArgAction::SetTrue)
                .num_args(0)
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
                .default_value("dist")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of threads to use, default to all logical cores")
                .required(false)
                .value_parser(clap::value_parser!(usize))
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
            .arg(
                Arg::new("model")
                .short('m')
                .long("model")
                .help("Equation used to calculate distance: 1 for poisson model or 0 for binomial model")
                .required(false)
                .value_parser(clap::value_parser!(u64))
                .default_value("1")
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("fp32")
                .long("fp32")
                .help("Distance output in float 32 instead of 64")
                .action(clap::ArgAction::SetTrue)
                .num_args(0)
            )
            .arg(
                Arg::new("dm")
                .long("dm")
                .help("Prints distance matrix")
                .action(clap::ArgAction::SetTrue)
                .num_args(0)
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sketch", s_matches)) => {
            // organize the inputs
            let sketch_file_name = s_matches.get_one::<String>("file").expect("required");
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").expect("required");
            let threads = s_matches
                .get_one::<usize>("threads")
                .copied()
                .unwrap_or_else(num_cpus::get);

            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build_global()
                .unwrap();

            let output_name = s_matches.get_one::<String>("output").expect("required");
            let alg = s_matches.get_one::<String>("algorithm").expect("required");
            let seed: u64 = *s_matches.get_one::<u64>("seed").expect("required");

            let aa = s_matches.get_flag("aa");

            let files: Vec<String> = {
                let f = File::open(sketch_file_name)?;
                BufReader::new(f)
                    .lines()
                    .filter_map(|l| l.ok())
                    .filter(|l| !l.trim().is_empty())
                    .collect()
            };

            let result: Result<(), Box<dyn Error>>;
            if alg == "hmh" {
                // create hypermash object and sketch
                result = sketch_files::<Sketch> (
                    None,
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                    aa
                );
            } else if alg == "hll" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = sketch_files::<HyperLogLog<i64>>(
                    Some(precision),
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                    aa
                );
            } else if alg == "ull" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = sketch_files::<UltraLogLog>(
                    Some(precision),
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                    aa
                );
            } else {
                // input for alg is not hmh, ull, or hll
                panic!("Algorithm must be either hmh, ull, or hll");
            }

            let molecule_param = if aa {
                "amino_acid".to_string()
            } else {
                "nucleotide".to_string()
            };

            // parameter JSONs
            let params;
            if alg == "ull" || alg == "hll" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                params = json!({
                    "k": kmer_length.to_string(),
                    "algorithm": alg,
                    "precision": precision.to_string(),
                    "seed": seed.to_string(),
                    "molecule": molecule_param
                });
            } else {
                params = json!({
                    "k": kmer_length.to_string(),
                    "algorithm": alg,
                    "seed": seed.to_string(),
                    "molecule": molecule_param
                });
            }

            // writing out
            File::create(format!("{}_parameters.json", output_name))?
                .write_all(serde_json::to_string_pretty(&params)?.as_bytes())?;

            result
        }
        Some(("dist", s_matches)) => {
            let ref_prefix = s_matches.get_one::<String>("reference").expect("required");
            let query_prefix = s_matches.get_one::<String>("query").expect("required");

            fn find_files(prefix: &str) -> std::io::Result<HashMap<&'static str, String>> {
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
            let threads = s_matches
                .get_one::<usize>("threads")
                .copied()
                .unwrap_or_else(num_cpus::get);

            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build_global()
                .unwrap();

            // go through the files needed, find name file, sketch file, and param file
            let ref_files = find_files(ref_prefix)?;
            let query_files = find_files(query_prefix)?;
            let ref_param_file = ref_files["params"].clone();
            let query_param_file = query_files["params"].clone();

            // println!("{}", ref_param_file);
            // println!("{}", query_param_file);

            // read in parameter json files into hashmaps
            let contents_ref = fs::read_to_string(ref_param_file)?;
            let contents_query = fs::read_to_string(query_param_file)?;
            let ref_map: HashMap<String, String> = serde_json::from_str(&contents_ref)?;
            let query_map: HashMap<String, String> = serde_json::from_str(&contents_query)?;

            // check that parameters match between ref and query genomes
            if ref_map["k"] != query_map["k"] {
                panic!("Genomes were not sketched with the same k");
            }
            if ref_map["algorithm"] != query_map["algorithm"] {
                panic!("Algorithms do not match in query and sketch genomes")
            }
            if ref_map["algorithm"] == "ull" || ref_map["algorithm"] == "hll" {
                if ref_map["precision"] != query_map["precision"] {
                    panic!(
                        "{} was not sketched with same precision btwn genomes",
                        ref_map["algorithm"]
                    );
                }
            }
            // assign kmer length once k matches
            let kmer_length: usize = ref_map["k"].parse()?;

            // function to read in names of the genomes, outputs a vector of names
            fn read_names(file_name: &str) -> std::io::Result<Vec<String>> {
                let file = File::open(file_name)?;
                let reader = BufReader::new(file);
                let names: Vec<String> = serde_json::from_reader(reader)?;
                Ok(names)
            }

            //create query sketch hashmap
            let query_namefile = query_files["files"].clone();
            let query_sketch_file_name = query_files["sketches"].clone();
            let query_names: Vec<String> = read_names(&query_namefile)?;

            // create reference sketch hashmap
            let ref_namefile = ref_files["files"].clone();
            let ref_sketch_file_name = ref_files["sketches"].clone();
            let reference_names: Vec<String> = read_names(&ref_namefile)?;

            let create_matrix = s_matches.get_flag("dm");
            let same_files = query_namefile == ref_namefile;
            let output = Arc::new(Mutex::new(File::create(output_file)?));
            let equation = *s_matches.get_one::<u64>("model").expect("required");
            let fp32 = s_matches.get_flag("fp32");

            if !create_matrix {
                let mut file = output.lock().unwrap();
                writeln!(file, "Reference\tQuery\tDistance")?;
            }
            
            // function to compute distance from fraction
            fn compute_distance<F: Float>(frac: F, kmer_length: usize, equation: u8) -> F {
                let k = F::from(kmer_length).unwrap();

                match equation {
                    1 => (-frac.ln() / k).min(F::one()),
                    0 => F::one() - frac.powf(F::one() / k),
                    _ => panic!("model needs to be 0 or 1"),
                }
            }

            // for triangular matrix printing 
            let file_idx = Arc::new(Mutex::new(HashMap::<String, usize>::new()));

            // callback and emit function
            fn print_dist<T: Float + std::fmt::Display>(
                distance_list: Vec<(&String, &String, T)>, 
                output: &Arc<Mutex<File>>, 
                create_matrix: bool, 
                same_files: bool, 
                file_idx: &Arc<Mutex<HashMap<String, usize>>>, 
                kmer_length: usize, 
                equation: u64) {
                // printing columns for matrix output using the query list
                let mut file = output.lock().unwrap();
                if create_matrix && !distance_list.is_empty() && distance_list[0].0.is_empty() {
                    for (i, col) in distance_list.iter().enumerate() {
                        write!(file, "\t{}", col.1).expect("Error writing columns for matrix output");
                        if same_files {
                            let mut idx = file_idx.lock().unwrap();
                            idx.insert(col.1.clone(), i);
                        }
                    }
                }
                else {
                    for (i, row) in distance_list.iter().enumerate() {
                        let r_name = row.0;
                        let q_name = row.1;
                        let d: T = if q_name == r_name {
                            T::zero()
                        } else {
                            compute_distance::<T>(row.2, kmer_length, equation as u8)
                        };
                        
                        if !create_matrix {
                            writeln!(file, "{}\t{}\t{:.6}", r_name, q_name, d)
                                .expect("Error writing to file");
                        } else {
                            if i == 0 { 
                                write!(file, "\n{}", r_name).expect("Error writing to file"); 
                            }
                            write!(file, "\t{:.6}", d).expect("Error writing to file");            
                        }
                    
                    }
                }
                //Ok(())
            }

            // for each algorithm, use a different generic depending on if user wants F32 or F64
            let result = if ref_map["algorithm"] == "hmh" {
                if fp32 {
                    let emit = move |rows: Vec<(&String, &String, f32)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    hmh_distance::<_, f32>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        create_matrix,
                        same_files,
                        emit,
                    )?
                } else {
                    let emit = move |rows: Vec<(&String, &String, f64)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    hmh_distance::<_, f64>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        create_matrix,
                        same_files,
                        emit
                    )?
                }
            } else if ref_map["algorithm"] == "ull" {
                let estimator = s_matches
                    .get_one::<String>("estimator")
                    .cloned()
                    .unwrap_or_else(|| "fgra".to_string());
                if fp32 {
                    let emit = move |rows: Vec<(&String, &String, f32)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    ull_distance::<_, f32>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        estimator, 
                        create_matrix,
                        same_files,
                        emit
                    )?
                } else {
                    let emit = move |rows: Vec<(&String, &String, f64)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    ull_distance::<_, f64>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        estimator, 
                        create_matrix,
                        same_files,
                        emit
                    )?
                }
            } else {
                // HLL
                if fp32 {
                    let emit = move |rows: Vec<(&String, &String, f32)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    hll_distance::<_, f32>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        create_matrix,
                        same_files,
                        emit
                    )?
                } else {
                    let emit = move |rows: Vec<(&String, &String, f64)>| {
                        print_dist(
                            rows, 
                            &output, 
                            create_matrix, 
                            same_files, 
                            &file_idx, 
                            kmer_length, 
                            equation
                        );
                    };
                    hll_distance::<_, f64>(
                        reference_names,
                        ref_sketch_file_name,
                        query_names,
                        query_sketch_file_name,
                        create_matrix,
                        same_files,
                        emit
                    )?
                }
            };

            println!("Distances computed.");
            Ok(result)
        }
        _ => Ok(()),
    }
}