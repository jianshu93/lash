use crossbeam_channel::bounded;
use hashbrown::HashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::error::Error;
//use xxhash_rust::xxh3::Xxh3Builder;
use crate::hasher::Xxh3Builder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::{Arc, Mutex};
use std::thread;
use ultraloglog::UltraLogLog;
use xxhash_rust::xxh3::xxh3_64;
use zstd::stream::{Decoder, Encoder};
// use std::path::Path;
use hyperminhash::Sketch;
use kmerutils::base::KmerT;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};
use log::info;
use serde_json::json;
use serde_json::to_writer_pretty;

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
    if b == 64 {
        v
    } else {
        v & ((1u64 << b) - 1)
    }
}

// function for distance
pub fn hmh_distance(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    kmer_length: usize,
    query_names: Vec<String>,
    query_sketch_file: String,
) -> Result<Vec<(String, String, f64)>, Box<dyn Error>> {
    // function to read sketches in using hypermash load
    fn read_sketches(file_name: &str, names: &Vec<String>) -> std::io::Result<Vec<Sketch>> {
        let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
        let reader = BufReader::new(file);

        // decompress files
        let mut decoder = Decoder::new(reader).expect("failed to create decompress");
        let mut sketches = Vec::new();
        for _name in names {
            let sketch = Sketch::load(&mut decoder)?;
            sketches.push(sketch);
        }
        Ok(sketches)
    }

    let q_sketch_vec: Vec<Sketch> = read_sketches(&query_sketch_file, &query_names)
        .expect(&format!("Error with reading from {}", query_sketch_file));
    let mut index = 0;

    let hasher = Xxh3Builder { seed: 93 }; // make sure ref/query pairs are in same order each time
    let mut query_sketches = HashMap::with_hasher(hasher.clone());
    for sketch in &q_sketch_vec {
        query_sketches.insert(&query_names[index], sketch);
        index += 1;
    }

    let r_sketch_vec = read_sketches(&ref_sketch_file, &reference_names)
        .expect(&format!("Error with reading from {}", ref_sketch_file));

    let mut reference_sketches = HashMap::with_hasher(hasher.clone());
    index = 0;
    for sketch in &r_sketch_vec {
        reference_sketches.insert(&reference_names[index], sketch);
        index += 1;
    }
    //Generate all pairs of query and reference sketches
    let pairs: Vec<(&String, &Sketch, &String, &Sketch)> = reference_sketches
        .iter()
        .flat_map(|(r_name, r_sketch)| {
            query_sketches
                .iter()
                .map(move |(q_name, q_sketch)| (*r_name, *r_sketch, *q_name, *q_sketch))
        })
        .collect();

    // Compute similarities and distances in parallel
    let results: Vec<(String, String, f64)> = pairs
        .par_iter()
        .map(
            |&(reference_name, reference_sketch, query_name, query_sketch)| {
                let similarity = query_sketch.similarity(reference_sketch);

                // Avoid division by zero and log of zero
                let adjusted_similarity = if similarity <= 0.0 {
                    std::f64::EPSILON // Small positive number to avoid log(0)
                } else {
                    similarity
                };

                // for debugging
                info!(
                    "Union: {}, a: {}, b: {}",
                    reference_sketch.clone().union(query_sketch).cardinality(),
                    reference_sketch.cardinality(),
                    query_sketch.cardinality()
                );

                // Calculate distance using the provided formula
                let numerator = 2.0 * adjusted_similarity;
                let denominator = 1.0 + adjusted_similarity;
                let fraction: f64 = numerator / denominator;
                let distance = -fraction.ln() / (kmer_length as f64);

                (reference_name.clone(), query_name.clone(), distance)
            },
        )
        .collect();
    Ok(results)
}

pub fn hmh_sketch(
    file: String,
    kmer_length: usize,
    output_name: String,
    threads: u32,
) -> Result<(), Box<dyn Error>> {
    //read list of genome paths
    let list_file = File::open(&file)?;
    let list_reader = BufReader::new(list_file);
    let sketch_files: Vec<String> = list_reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.trim().is_empty())
        .collect();

    // process each file in parallel
    let sketches: HashMap<String, Sketch> = sketch_files
        .par_iter()
        .map(|fname| {
            let file_name = fname.clone();

            // streaming reader for this file
            let mut reader = parse_fastx_file(&file_name).expect("Failed to parse file");

            // bounded channel for batches of sequences
            let (sender, receiver) = bounded::<Vec<Vec<u8>>>(64);

            // spawn producer thread
            let producer = thread::spawn(move || {
                let mut batch: Vec<Vec<u8>> = Vec::with_capacity(5000);
                while let Some(rec_res) = reader.next() {
                    if let Ok(rec) = rec_res {
                        batch.push(rec.seq().to_vec());
                        if batch.len() == 5000 {
                            if sender.send(batch).is_err() {
                                return;
                            }
                            batch = Vec::with_capacity(5000);
                        }
                    }
                }
                if !batch.is_empty() {
                    let _ = sender.send(batch);
                }
            });

            // one global sketch per file behind a mutex
            let global = Arc::new(Mutex::new(Sketch::default()));

            // consume batches
            for batch in receiver {
                // process each sequence in batch in parallel
                batch.par_iter().for_each(|seq| {
                    let seq_vec = filter_out_n(seq);
                    if seq_vec.is_empty() {
                        return;
                    }
                    let kseq = KSeq::new(&seq_vec, 2);
                    let mut local = Sketch::default();

                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            local.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            local.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            local.add_bytes(&masked.to_le_bytes());
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }

                    // merge local into global
                    {
                        let mut g = global.lock().unwrap();
                        g.union(&local);
                    }
                });
            }

            producer.join().expect("Producer thread panicked");

            // unwrap final sketch
            let final_sketch = Arc::try_unwrap(global)
                .expect("Arc still has multiple refs")
                .into_inner()
                .unwrap();

            (file_name, final_sketch)
        })
        .collect();

    // write names list
    let names: Vec<String> = sketches.keys().cloned().collect();
    to_writer_pretty(
        &File::create(format!("{}_files.json", output_name))?,
        &names,
    )?;

    // serialize sketches
    let out_bin = format!("{}_sketches.bin", output_name);
    let bin_file = File::create(&out_bin)?;
    let mut writer = BufWriter::new(bin_file);
    let mut encoder = Encoder::new(&mut writer, 3).expect("failed to create compression");
    encoder
        .multithread(threads)
        .expect("failed to multithread compressor");
    for name in &names {
        sketches[name].save(&mut encoder)?;
    }
    encoder.finish().expect("failed to compress");

    // parameter JSON
    let params = json!({
        "k": kmer_length.to_string(),
        "algorithm": "hmh"
    });
    let mut param_file = File::create(format!("{}_parameters.json", output_name))?;
    param_file.write_all(serde_json::to_string_pretty(&params).unwrap().as_bytes())?;

    println!("Serialized sketches with hmh");
    Ok(())
}

pub fn ull_sketch(
    precision: u32,
    sketch_file_name: String,
    kmer_length: usize,
    output_name: String,
    threads: u32,
) -> Result<(), Box<dyn Error>> {
    let sketch_file = File::open(sketch_file_name)?;
    let sketch_reader = BufReader::new(sketch_file);
    let files: Vec<String> = sketch_reader
        .lines()
        .filter_map(|line| line.ok())
        .filter(|line| !line.trim().is_empty())
        .collect();

    // loop through each file, create a ull for each one, add kmers
    let ull_vec: Vec<UltraLogLog> = files
        .par_iter()
        .map(|file| {
            let ull = Mutex::new(UltraLogLog::new(precision).expect("failed to create ull"));
            let mut reader = parse_fastx_file(file).expect("Invalid input file");
            let (sender, receiver) = crossbeam_channel::bounded(64); // Channel with capacity 10

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
            for batch in receiver {
                batch.par_iter().for_each(|seq| {
                    //let seqrec = record.expect("Error reading record");
                    let seq_vec: Vec<u8> = filter_out_n(&seq);
                    let kseq = KSeq::new(&seq_vec, 2);
                    let mut ull_guard = ull.lock().unwrap();
                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            let hash64 = xxh3_64(&masked.to_le_bytes());
                            ull_guard.add(hash64);
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            let hash64 = xxh3_64(&masked.to_le_bytes());
                            ull_guard.add(hash64);
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            let hash64 = xxh3_64(&masked.to_le_bytes());
                            ull_guard.add(hash64);
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }
                })
            }
            reader_thread.join().expect("Reader thread panicked");
            let ull_final = Mutex::into_inner(ull).unwrap();
            ull_final
        })
        .collect();

    // set up file and writer for sketching
    let sketch_output =
        File::create(format!("{}_sketches.bin", &output_name)).expect("Failed to create file");
    let writer = BufWriter::new(sketch_output);

    // compress sketches
    let mut encoder = Encoder::new(writer, 3).expect("failed to create compression");
    encoder
        .multithread(threads)
        .expect("failed to multithread compressor");
    for ull in &ull_vec {
        ull.save(&mut encoder).expect("Failed to save UltraLogLog");
    }
    encoder.finish().expect("failed to compress");

    // write all file names to a json file
    to_writer_pretty(
        &File::create(format!("{}_files.json", &output_name))?,
        &files,
    )?;

    // save a json of parameters used
    let data = json!({
        "k": kmer_length.to_string(),
        "algorithm": "ull",
        "precision": precision.to_string(),
    });
    let json_str = serde_json::to_string_pretty(&data).unwrap();
    let mut param_file = File::create(format!("{}_parameters.json", output_name))?;
    param_file.write_all(json_str.as_bytes())?;

    println!("Serialized sketches with ull");
    Ok(())
}
