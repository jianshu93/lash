use needletail::{parse_fastx_file, Sequence};
// use needletail::kmer::Kmers;
// use needletail::sequence::canonical;
use rayon::prelude::*;
// use rayon::ThreadPoolBuilder;
use crossbeam_channel::bounded;
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::thread;
// use std::path::Path;
use hyperminhash::Sketch;
use serde_json::{to_writer_pretty};
use serde_json::json;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};
use kmerutils::base::KmerT;

pub fn filter_out_n(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if matches!(c, b'A' | b'C' | b'T' | b'G') {
            out.push(c);
        }
    }
    out
}

pub fn mask_bits(v: u64, k: usize) -> u64 {
    let b = 2 * k as u32;
    if b == 64 { v } else { v & ((1u64 << b) - 1) }
}


// function for distance
pub fn hmh_distance(reference_names: Vec<String>, ref_sketch_file: String, kmer_length: usize,
    query_names: Vec<String>, query_sketch_file: String) 
    -> Result<Vec<(String, String, f64)>, Box<dyn Error>> {

    // function to read sketches in using hypermash load
    fn read_sketches(file_name: &str, names: &Vec<String>) -> std::io::Result<Vec<Sketch>> {
        let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
        let mut reader = BufReader::new(file);
        let mut sketches = Vec::new();
        for name in names {
            let sketch = Sketch::load(&mut reader)?;
            sketches.push(sketch);
        }
        Ok(sketches)
    }

    let q_sketch_vec: Vec<Sketch> = read_sketches(&query_sketch_file, &query_names)
                        .expect(&format!("Error with reading from {}", query_sketch_file));
    let mut index = 0;
    let mut query_sketches: HashMap<&String, &Sketch> = HashMap::new();
    for sketch in &q_sketch_vec {
        query_sketches.insert(&query_names[index], sketch);
        index += 1;
    }

    let r_sketch_vec = read_sketches(&ref_sketch_file, &reference_names)
                        .expect(&format!("Error with reading from {}", ref_sketch_file));
    let mut reference_sketches = HashMap::new();
    index = 0;
    for sketch in &r_sketch_vec {
        reference_sketches.insert(&reference_names[index], sketch);
        index += 1;
    }
    //Generate all pairs of query and reference sketches
    let pairs: Vec<(&String, &Sketch, &String, &Sketch)> = reference_sketches
    .iter()
    .flat_map(|(r_name, r_sketch)| {
        query_sketches.iter().map(move |(q_name, q_sketch)| {
            (*r_name, *r_sketch, *q_name, *q_sketch)
        })
    })
    .collect();


    // Compute similarities and distances in parallel
    let results: Vec<(String, String, f64)> = pairs
        .par_iter()
        .map(|&(reference_name, reference_sketch, query_name, query_sketch)| {
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

            (reference_name.clone(), query_name.clone(), distance)
        })
        .collect();
    Ok(results)
}

pub struct Hypermash {
    file: String, 
    kmer_length: usize, 
    output_name: String, 
}

impl Hypermash {
    // create a new hypermash object
    pub fn new(file_name: String, k: usize, option: String) -> Hypermash {
        Hypermash {
            file: file_name, 
            kmer_length: k, 
            output_name: option
        }
    }

    // sketches a list of genome file names
    pub fn sketch(&self) -> Result<(), Box<dyn Error>> {
        let sketch_file = File::open(&self.file)?;
        let sketch_reader = BufReader::new(sketch_file);
        let sketch_files: Vec<String> = sketch_reader
            .lines()
            .filter_map(|line| line.ok())
            .filter(|line| !line.trim().is_empty())
            .collect();


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
                let local_sketch = batch.par_iter().map(|seq| {
                    let seq_vec = filter_out_n(seq);
                    if seq_vec.is_empty() { return Box::new(Sketch::default()); }
                    let mut sketch = Box::new(Sketch::default());
                    let kseq = KSeq::new(&seq_vec, 2);
                    if self.kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(self.kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, self.kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if self.kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, self.kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if self.kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(self.kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), self.kmer_length);
                            sketch.add_bytes(&masked.to_le_bytes());
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }
                    sketch
                })
                .reduce(|| Box::new(Sketch::default()), |mut a, b| { a.union(&b); a });

                global_sketch.union(&local_sketch);
            }

            // Wait for the reader thread to finish
            reader_thread.join().expect("Reader thread panicked");

            (file_name, *global_sketch)
        })
        .collect();

        // put names into a json file
        let mut names: Vec<String> = sketches.keys().cloned().collect();
        to_writer_pretty(
            &File::create(format!("{}_files.json", self.output_name))?,
            &names
        )?;

        // serialize sketches
        let file_name = format!("{}_sketches.bin",  self.output_name);
        let mut sketch_writer = BufWriter::new(File::create(file_name)?);

        for name in &names {
            let sketch = &sketches[name];
            sketch.save(&mut sketch_writer)?;
        }

        // write to a parameter json file
        let data = json!({
            "k": self.kmer_length.to_string(),
            "algorithm": "hmh"
        });
    
        let json_str = serde_json::to_string_pretty(&data).unwrap();
        let mut param_file = File::create(format!("{}_parameters.json", 
            self.output_name))?;
        param_file.write_all(json_str.as_bytes())?;

        println!("Serialized sketches with hmh");
        Ok(())
    }

}