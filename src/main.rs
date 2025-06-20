use clap::{Arg, ArgAction, Command};
use needletail::{parse_fastx_file, Sequence};
// use needletail::kmer::Kmers;
// use needletail::sequence::canonical;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
//use crossbeam_channel::bounded;
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::fs;
use std::io::{BufRead, BufReader, Write, BufWriter};
//use std::thread;
use std::path::Path;
use serde_json::{to_writer_pretty};
use serde_json::json;
use xxhash_rust::xxh3::xxh3_64;

mod hypermash;
use hypermash::Hypermash; 
use crate::hypermash::{hmh_distance, mask_bits, filter_out_n};
use ultraloglog::{UltraLogLog, Estimator, MaximumLikelihoodEstimator};

use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};
use kmerutils::base::KmerT;


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
                Arg::new("file")
                .short('f')
                .long("file")
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
                .help("Specifiy precision, for ull only. Default 10.")
                .required(false)
                .value_parser(clap::value_parser!(usize))
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
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("query")
                .short('q')
                .long("query")
                .help("Prefix to search for query genome files")
                .required(true)
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
            .arg(
                Arg::new("estimator")
                .short('e')
                .long("estimator")
                .help("Specify estimator (fgra or ml), for ull only. Default fgra. ")
                .required(false)
                .action(ArgAction::Set)
            )
            
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sketch", s_matches)) => {
            // organize the inputs
            let sketch_file_name = s_matches.get_one::<String>("file").expect("required");
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").expect("required");
            let threads: usize = *s_matches.get_one::<usize>("threads").unwrap_or(&num_cpus::get());

            let output_name = s_matches.get_one::<String>("output").expect("required");
            let alg = s_matches.get_one::<String>("algorithm").expect("required");

            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();

            if alg == "hmh" {
                // create hypermash object and sketch
                let hm = Hypermash::new(sketch_file_name.clone(), kmer_length, output_name.clone());
                hm.sketch();
            }
            else if alg == "ull" {
                let precision: u32 = *s_matches.get_one::<usize>("precision").unwrap_or(&10) as u32;

                // create vector of files to be sketched
                let sketch_file = File::open(sketch_file_name)?;
                let sketch_reader = BufReader::new(sketch_file);
                let files: Vec<String> = sketch_reader
                    .lines()
                    .filter_map(|line| line.ok())
                    .filter(|line| !line.trim().is_empty())
                    .collect();

                // set up file and writer for sketching
                let sketch_output = File::create(format!("{}_sketches.bin", &output_name)).expect("Failed to create file");
                let mut writer = BufWriter::new(sketch_output);

                // loop through each file, create a ull for each one, add kmers 
                for file in &files {
                    let mut ull = UltraLogLog::new(precision).expect("failed to create ull");                    // add kmers into ull
                    let mut reader = parse_fastx_file(file).expect("Invalid input file");

                    while let Some(record) = reader.next() {
                        let seqrec = record.expect("Error reading record");
                        let seq_vec: Vec<u8> = filter_out_n(&seqrec.seq());
                        let kseq = KSeq::new(&seq_vec, 2);

                        if kmer_length <= 14 {
                            let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                                let hash64 = xxh3_64(&masked.to_le_bytes());
                                ull.add(hash64);
                            }
                        }
                        else if kmer_length == 16 {
                            let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                                let hash64 = xxh3_64(&masked.to_le_bytes());
                                ull.add(hash64);
                            }
                        } else if kmer_length <= 32 {
                            let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                                let hash64 = xxh3_64(&masked.to_le_bytes());
                                ull.add(hash64);
                            }
                        } else {
                            panic!("k-mer length must be 1–32, k=15 is not supported");
                        }
                    }

                    // serialize sketches
                    ull.save(&mut writer).expect("Failed to save UltraLogLog");
                }

                // write all file names to a json file
                to_writer_pretty(
                    &File::create(format!("{}_files.json", &output_name))?,
                    &files
                )?;

                // save a json of parameters used
                let data = json!({
                    "k": kmer_length.to_string(),
                    "algorithm": "ull", 
                    "precision": precision.to_string(),
                });
                let json_str = serde_json::to_string_pretty(&data).unwrap();
                let mut param_file = File::create(format!("{}_parameters.json", 
                    output_name))?;
                param_file.write_all(json_str.as_bytes())?;

                println!("Serialized sketches with ull");
            }
            else { // input for alg is not hmh or ull
                panic!("Algorithm must be either hmh or ull");
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

                // if files.len() != 3 {
                //     panic!("There should be 3 files starting with {} but {} were found instead", 
                //     prefix, 
                //     files.len());
                // }
                let mut file_map: HashMap<&str, String> = HashMap::new();
                for file in files {
                    if file.ends_with("parameters.json") {
                        file_map.insert("params",file);
                    }
                    else if file.ends_with("files.json") {
                        file_map.insert("files",file);
                    }
                    else if file.ends_with(".bin") { // .bin file for sketches
                        file_map.insert("sketches",file);
                    }
                }
                if file_map.keys().len() != 3 {
                    panic!("There should be 3 files starting with {} but {} were found instead", 
                    prefix, 
                    file_map.keys().len());
                }
                Ok(file_map)
            }
            
            let output_file: &String = s_matches.get_one::<String>("output_file").expect("required");
            let threads: usize = *s_matches.get_one::<usize>("threads").unwrap_or(&num_cpus::get());

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
            let ref_map: HashMap<String, String> = serde_json::from_str(&contents_ref).expect("Invalid JSON");
            let query_map: HashMap<String, String> = serde_json::from_str(&contents_query).expect("Invalid JSON");
            
            // check that parameters match between ref and query genomes
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
                
            }
            // assign kmer length once k matches
            let kmer_length:usize = ref_map["k"].parse().expect("Not a valid usize");

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
            
            let mut results: Vec<(String, String, f64)> = Vec::new();
            if ref_map["algorithm"] == "hmh" {
                results = hmh_distance(reference_names, ref_sketch_file_name, kmer_length, query_names, query_sketch_file_name).unwrap();
            }
            else if ref_map["algorithm"] == "ull" {
                
                let ref_sketch_file = File::open(ref_sketch_file_name).expect("Failed to open file");
                let query_sketch_file = File::open(query_sketch_file_name).expect("Failed to open file");
                let estimator = s_matches.get_one::<String>("estimator").cloned().unwrap_or_else(|| "fgra".to_string());

                // function that creates a hashmap holding name of genome and the ull and cardinality for it
                fn create_ull_map(sketch_file: File, names: &Vec<String>, estimator:&String ) -> Result<HashMap<String, (UltraLogLog, f64)>, std::io::Error> {
                    let mut sketches: HashMap<String, (UltraLogLog, f64)> = HashMap::new();
                    let mut reader = BufReader::new(sketch_file);
                    for file in names {
                        let ref_ull = UltraLogLog::load(&mut reader)?;
                        let c: f64 = match estimator.as_str() {
                            "fgra" => ref_ull.get_distinct_count_estimate(),
                            "ml" => MaximumLikelihoodEstimator.estimate(&ref_ull),
                            _ => panic!("estimator needs to be either fgra or ml"),
                        };
                        sketches.insert(file.clone(), (ref_ull, c));
                    }
                    Ok(sketches)
                }
                let ref_map = create_ull_map(ref_sketch_file, &reference_names, &estimator).unwrap();
                let query_map = create_ull_map(query_sketch_file, &query_names, &estimator).unwrap();
                
                // create all pairs between reference and query
                let pairs: Vec<(&str, &str)> = ref_map.keys()
                    .flat_map(|k1| {
                        query_map.keys().map(move |k2| (k1.as_str(), k2.as_str()))
                    })
                    .collect();
                
                results = pairs
                    .par_iter()
                    .map(|&(reference_name, query_name)| {
                        let a: f64 = ref_map[reference_name].1; // cardinality of reference
                        //println!("{}", a);
                        let b: f64 = query_map[query_name].1; // cardinality of query
                        //println!("{}", b);
                        let union_ull = UltraLogLog::merge(
                            &ref_map[reference_name].0, 
                            &query_map[query_name].0,)
                            .expect("failed to merge sketches");

                        let intersection = union_ull.get_distinct_count_estimate();
                        let similarity = (a + b - intersection) / intersection;
                        let adjusted_similarity = if similarity <= 0.0 {
                            std::f64::EPSILON // Small positive number to avoid log(0)
                        } else {
                            similarity
                        };
                        let numerator: f64 = 2.0 * adjusted_similarity;
                        let denominator: f64 = 1.0 + adjusted_similarity;
                        let fraction: f64 = numerator / denominator;
                        let distance: f64 = -fraction.ln() / (kmer_length as f64);
                        (reference_name.to_string(), query_name.to_string(), distance)
                    })
                    .collect();
            } 

            // Open the output file for writing
            let mut output = File::create(format!("{}.txt", output_file))?;

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
    