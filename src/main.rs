use clap::{Arg, ArgAction, Command};
use needletail::{parse_fastx_file, Sequence};
use needletail::kmer::Kmers;
use needletail::sequence::canonical;
use rayon::prelude::*;
use crossbeam_channel::bounded;
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::thread;
use std::path::Path;
use hyperminhash::Sketch;

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
                .help("One file containing list of FASTA files (.gz supported)")
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
        )
        .subcommand(
            Command::new("distance")
            .about("Computes distance between sketches")
            .arg(
                Arg::new("reference")
                .short('r')
                .long("reference")
                .help("File for list of names and file for sketches for the reference genome, respectively, 
                separated by space")
                .required(true)
                .num_args(2)
                .action(ArgAction::Set)
            )
            .arg(
                Arg::new("query")
                .short('q')
                .long("query")
                .help("File for list of names and file for sketches for the query genome, respectively, 
                separated by space")
                .required(true)
                .num_args(2)
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
                Arg::new("output_file")
                .short('o')
                .long("output_file")
                .help("Name of output file to write results")
                .required(true)
                .action(ArgAction::Set)
            )
        )
        .get_matches();

    match matches.subcommand() {
        Some(("sketch", s_matches)) => {
            let sketch_file_name = s_matches.get_one::<String>("sketch").expect("required");
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").unwrap();
            
            // Read the lists of query and reference files
            let sketch_file = File::open(sketch_file_name)?;
            let sketch_reader = BufReader::new(sketch_file);
            let sketch_files: Vec<String> = sketch_reader
                .lines()
                .filter_map(|line| line.ok())
                .filter(|line| !line.trim().is_empty())
                .collect();

            // Process files and create sketches
            let sketches: HashMap<String, Sketch> = sketch_files
            .par_iter()
            .map(|file_name| {
                // Clone the file name for ownership in the closure
                let file_name = file_name.clone();

                // Process the file and create a sketch
                let mut reader = parse_fastx_file(&file_name).expect("Failed to parse file");
                let (sender, receiver) = bounded(10); // Channel with capacity 10

                // Spawn a thread to read sequences and send batches
                let reader_thread = thread::spawn(move || {
                    let mut batch = Vec::new();
                    while let Some(result) = reader.next() {
                        if let Ok(seqrec) = result {
                            batch.push(seqrec.seq().to_vec());
                            if batch.len() == 5000 {
                                if sender.send(batch.clone()).is_err() {
                                    break; // Receiver has hung up
                                }
                                batch.clear();
                            }
                        }
                    }
                    if !batch.is_empty() {
                        let _ = sender.send(batch); // Send remaining batch
                    }
                });

                // Initialize an empty Sketch allocated on the heap
                let mut global_sketch = Box::new(Sketch::default());

                // Process the batches
                for batch in receiver {
                    // Process the batch in parallel
                    let local_sketch = batch
                        .par_iter()
                        .map(|seq| {
                            // Allocate sketch on the heap
                            let mut sketch = Box::new(Sketch::default());
                            let kmer_length_u8 = kmer_length as u8;
                            for kmer in Kmers::new(seq, kmer_length_u8) {
                                let kmer_bytes = canonical(kmer); // Use canonical function
                                sketch.add_bytes(&kmer_bytes);
                            }
                            sketch
                        })
                        .reduce(
                            || Box::new(Sketch::default()),
                            |mut a, b| {
                                a.union(&b);
                                a
                            },
                        );

                    // Merge the local sketch into the global sketch
                    global_sketch.union(&local_sketch);
                }

                // Wait for the reader thread to finish
                reader_thread.join().expect("Reader thread panicked");

                (file_name, *global_sketch)
            })
            .collect();

            let name_file: String = format!("{}_names.bin", sketch_file_name.clone());
            let sketch_list = File::create(name_file)?;
            let mut list_writer = BufWriter::new(sketch_list);
            let sketch_file_name: String = format!("{}_sketch.bin", sketch_file_name.clone());
            let sketch = File::create(sketch_file_name)?;
            let mut sketch_writer = BufWriter::new(sketch);

            // serialize query sketches
            for (string, sketch) in sketches {
                writeln!(list_writer, "{}", string)?;
                sketch.save(&mut sketch_writer)?;
            }
            println!("Serialized sketches");

            Ok(())
        }
        Some(("distance", s_matches)) => {
            let query_files: Vec<&str> = s_matches
                .get_many::<String>("query")
                .expect("required")
                .map(|s| s.as_str())
                .collect();
            let query_namefile = query_files[0];
            let query_sketch_file = query_files[1];

            let ref_files: Vec<&str> = s_matches
                .get_many::<String>("reference")
                .expect("required")
                .map(|s| s.as_str())
                .collect();
            let ref_namefile = ref_files[0];
            let ref_sketch_file = ref_files[1];
            let kmer_length: usize = *s_matches.get_one::<usize>("kmer_length").unwrap();

            let output_file: &String = s_matches.get_one::<String>("output_file").expect("required");

            // function to read in names of the genomes into a vector of names
            fn read_names(file_name: &str) -> std::io::Result<Vec<String>> {
                let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
                let reader = BufReader::new(file);
                let mut names = Vec::new();
                for line in reader.lines() {
                    let name = line?;
                    names.push(name);
                }
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
            let query_names: Vec<String> = read_names(&query_namefile)
                            .expect(&format!("Error with reading from {}", query_namefile));
            let q_sketch_vec: Vec<Sketch> = read_sketches(query_sketch_file, &query_names)
                            .expect(&format!("Error with reading from {}", query_sketch_file));
            
            let mut index = 0;
            let mut query_sketches: HashMap<&String, &Sketch> = HashMap::new();
            for sketch in &q_sketch_vec {
                query_sketches.insert(&query_names[index], sketch);
                index += 1;
            }
            
            // create reference sketch hashmap
            let reference_names: Vec<String> = read_names(&ref_namefile)
                            .expect(&format!("Error with reading from {}", ref_namefile));
            let r_sketch_vec = read_sketches(ref_sketch_file, &reference_names)
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
        
            Ok(())
        }
        _ => {
            Ok(())
        }
    }
}
    

