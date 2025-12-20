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
mod hasher;
use serde_json::json;
mod utils;
use crate::utils::{hll_distance, hll_sketch, hmh_distance, hmh_sketch, ull_sketch, ull_distance};
use num_traits::Float;

pub enum Matrix {
    F32(Vec<Vec<f32>>),
    F64(Vec<Vec<f64>>),
}

impl Matrix {
    pub fn set<F: Float>(&mut self, i: usize, j: usize, value: F) {
        match self {
            Matrix::F32(m) => m[i][j] = F::to_f32(&value).unwrap(),
            Matrix::F64(m) => m[i][j] = F::to_f64(&value).unwrap(),
        }
    }
}

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
                .unwrap_or_else(|| num_cpus::get());

            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build_global()
                .unwrap();

            let output_name = s_matches.get_one::<String>("output").expect("required");
            let alg = s_matches.get_one::<String>("algorithm").expect("required");
            let seed: u64 = *s_matches.get_one::<u64>("seed").expect("required");

            let files: Vec<String> = {
                let f = File::open(&sketch_file_name)?;
                BufReader::new(f)
                    .lines()
                    .filter_map(|l| l.ok())
                    .filter(|l| !l.trim().is_empty())
                    .collect()
            };

            let result: Result<(), Box<dyn Error>>;
            if alg == "hmh" {
                // create hypermash object and sketch
                result = hmh_sketch(
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                );
            } else if alg == "hll" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = hll_sketch(
                    precision,
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                );
            } else if alg == "ull" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                result = ull_sketch(
                    precision,
                    files,
                    kmer_length,
                    output_name.clone(),
                    threads as u32,
                    seed,
                );
            } else {
                // input for alg is not hmh, ull, or hll
                panic!("Algorithm must be either hmh, ull, or hll");
            }

            // parameter JSONs
            let params;
            if alg == "ull" || alg == "hll" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;
                params = json!({
                    "k": kmer_length.to_string(),
                    "algorithm": &alg,
                    "precision": precision.to_string(),
                    "seed": seed.to_string()
                });
            } else {
                params = json!({
                    "k": kmer_length.to_string(),
                    "algorithm": &alg,
                    "seed": seed.to_string()
                });
            }

            // writing out
            File::create(format!("{}_parameters.json", &output_name))?
                .write_all(serde_json::to_string_pretty(&params)?.as_bytes())?;

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
            let threads = s_matches
                .get_one::<usize>("threads")
                .copied()
                .unwrap_or_else(|| num_cpus::get());

            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build_global()
                .unwrap();

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
                    panic!(
                        "{} was not sketched with same precision btwn genomes",
                        ref_map["algorithm"]
                    );
                }
            }
            // assign kmer length once k matches
            let kmer_length: usize = ref_map["k"].parse().expect("Not a valid usize");

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

            let create_matrix = s_matches.get_flag("dm");

            // let mut matrix = if fp32 {
            //     Matrix::F32(Vec::new())
            // } else {
            //     Matrix::F64(Vec::new())
            // };
            // let mut query_idx = HashMap::new();
            // let mut ref_idx = HashMap::new();
            
            // if create_matrix {
            //     if fp32 {
            //         matrix = Matrix::F32(vec![vec![0.0f32; reference_names.len()]; query_names.len()]);
            //     } else {
            //         matrix = Matrix::F64(vec![vec![0.0f64; reference_names.len()]; query_names.len()]);
            //     }
            //     for (i, name) in query_names.iter().enumerate() {
            //         query_idx.insert(name.clone(), i);
            //     }
            //     for (i, name) in reference_names.iter().enumerate() {
            //         ref_idx.insert(name.clone(), i);
            //     }
            // }

            let mut output = File::create(output_file)?; // don't append ".txt"
            let equation = *s_matches.get_one::<u64>("model").expect("required");
            let fp32 = s_matches.get_flag("fp32");

            if !create_matrix {
                writeln!(output, "Reference\tQuery\tDistance")?;
            }
            
            // function to compute distance from fraction
            fn compute_distance<F: Float>(frac: f64, kmer_length: usize, equation: u8) -> F {
                match equation {
                    1 => {
                        // -ln(frac)/kmer_length, clamped to 1.0
                        let frac_f = F::from(frac).unwrap();
                        let k = F::from(kmer_length).unwrap();
                        (-frac_f.ln() / k).min(F::one())
                    }
                    0 => {
                        // 1 - frac^(1/kmer_length)
                        let frac_f = F::from(frac).unwrap();
                        let k = F::from(kmer_length).unwrap();
                        F::one() - frac_f.powf(F::one() / k)
                    }
                    _ => panic!("model needs to be 0 or 1"),
                }
            }
            

            // for knowing when to use "-" in matrix
            let same_files = ref_namefile == query_namefile;
            let mut row = 0;
            let mut col = 0;

            // callback and emit function
            let print_list = | r_name: &String, q_name: &String, frac: f64| {

                // empty r_name means printing matrix column names (q_name)
                if r_name == "" {
                    write!(output, "\t{}", q_name)?;
                }
                // empty q_name means printing a new row (r_name) for matrix
                else if q_name == "" { 
                    writeln!(output, "")?;
                    write!(output, "{}", r_name)?;
                    col = 0;
                    row += 1;
                }

                // print a distance value
                else {
                    // for float32 output
                    if fp32 {
                        let d: f32 = if q_name == r_name {
                            0.0
                        } else {
                            compute_distance::<f32>(frac, kmer_length, equation as u8)
                        };

                        if !create_matrix {
                            writeln!(output, "{}\t{}\t{:.6}", r_name, q_name, d).expect("Error writing to file");
                        }
                        else {
                            col += 1;
                            if row > col && same_files {
                                write!(output, "\t-").expect("Error writing to file"); // redundant distance
                            } else {
                                write!(output, "\t{:.6}", d).expect("Error writing to file");
                            }
                        }
                    } else { // for float64 (default) output
                        let d: f64 = if q_name == r_name {
                            0.0
                        } else {
                            compute_distance::<f64>(frac, kmer_length, equation as u8)
                        };

                        if !create_matrix {
                            writeln!(output, "{}\t{}\t{:.6}", r_name, q_name, d).expect("Error writing to file");
                        }
                        else {
                            col += 1;
                            if row > col && same_files {
                                write!(output, "\t-").expect("Error writing to file"); // redundant distance
                            } else {
                                write!(output, "\t{:.6}", d).expect("Error writing to file");
                            }
                        }
                    }
                }
                Ok(())
            };

            // let mut results: Vec<(String, String, f64)> = Vec::new();
            // compute all pairwise distances into a single results Vec
            let result = if ref_map["algorithm"] == "hmh" {
                hmh_distance(
                    reference_names,
                    ref_sketch_file_name,
                    query_names,
                    query_sketch_file_name,
                    create_matrix,
                    print_list
                )
                .unwrap()
            } else if ref_map["algorithm"] == "ull" {

                let estimator = s_matches
                    .get_one::<String>("estimator")
                    .cloned()
                    .unwrap_or_else(|| "fgra".to_string());

                ull_distance(
                    reference_names,
                    ref_sketch_file_name,
                    query_names,
                    query_sketch_file_name,
                    estimator, 
                    create_matrix,
                    print_list
                ).unwrap()
            } else {
                // HLL
                hll_distance(
                    reference_names,
                    ref_sketch_file_name,
                    query_names,
                    query_sketch_file_name,
                    create_matrix,
                    print_list
                ).unwrap()
            };

            println!("Distances computed.");
            Ok(result)
        }
        _ => Ok(()),
    }
}
