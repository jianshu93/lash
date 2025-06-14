use needletail::{parse_fastx_file, Sequence};
use needletail::kmer::Kmers;
use needletail::sequence::canonical;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use crossbeam_channel::bounded;
use std::error::Error;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::thread;
use std::path::Path;
use hyperminhash::Sketch;
use serde_json::{to_writer_pretty};
use serde_json::json;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};
use kmerutils::base::KmerT;

fn filter_out_n(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if matches!(c, b'A' | b'C' | b'T' | b'G') {
            out.push(c);
        }
    }
    out
}

fn mask_bits(v: u64, k: usize) -> u64 {
    let b = 2 * k as u32;
    if b == 64 { v } else { v & ((1u64 << b) - 1) }
}
pub struct Hypermash {
    file: String, 
    kmer_length: usize, 
    threads: usize, 
    output_name: String, 
}

impl Hypermash {
    // create a new hypermash object
    pub fn new(file_name: String, k: usize, t: usize, option: String) -> Hypermash {
        Hypermash {
            file: file_name, 
            kmer_length: k, 
            threads: t, 
            output_name: option
        }
    }

    // sketches a list of genome file names
    pub fn sketch(&self) -> Result<(), Box<dyn Error>> {
        ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build_global()
            .unwrap();

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
            &File::create(format!("{}_names.json", self.output_name))?,
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

        println!("Serialized sketches.");
        Ok(())
    }

}