// src/utils.rs

use hashbrown::HashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::error::Error;

use crate::hasher::Xxh3Builder;

use std::fs::File;
use std::io::{BufReader, BufWriter};

use ultraloglog::UltraLogLog;
use xxhash_rust::xxh3::xxh3_64_with_seed;
// use xxhash_rust::xxh3::xxh3_64;

use zstd::stream::{Decoder, Encoder};

use hyperminhash::Sketch;
use kmerutils::base::KmerT;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};

use log::info;
use serde_json::to_writer_pretty;
use streaming_algorithms::HyperLogLog;

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

// distances
pub fn hmh_distance(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    kmer_length: usize,
    query_names: Vec<String>,
    query_sketch_file: String,
) -> Result<Vec<(String, String, f64)>, Box<dyn Error>> {
    fn read_sketches(file_name: &str, names: &Vec<String>) -> std::io::Result<Vec<Sketch>> {
        let file = File::open(file_name).expect(&format!("Error opening {}", file_name));
        let reader = BufReader::new(file);
        let mut decoder = Decoder::new(reader).expect("failed to create decompress");
        let mut sketches = Vec::with_capacity(names.len());
        for _name in names {
            let sketch = Sketch::load(&mut decoder)?;
            sketches.push(sketch);
        }
        Ok(sketches)
    }

    let q_sketch_vec: Vec<Sketch> = read_sketches(&query_sketch_file, &query_names)
        .expect(&format!("Error with reading from {}", query_sketch_file));
    let mut index = 0;

    // stable hasher to keep key order deterministic
    let hasher = Xxh3Builder { seed: 93 };
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

    // Cartesian product of (reference, query)
    let pairs: Vec<(&String, &Sketch, &String, &Sketch)> = reference_sketches
        .iter()
        .flat_map(|(r_name, r_sketch)| {
            query_sketches
                .iter()
                .map(move |(q_name, q_sketch)| (*r_name, *r_sketch, *q_name, *q_sketch))
        })
        .collect();

    let results: Vec<(String, String, f64)> = pairs
        .par_iter()
        .map(
            |&(reference_name, reference_sketch, query_name, query_sketch)| {
                let similarity = query_sketch.similarity(reference_sketch).max(0.0);

                // for debugging
                info!(
                    "Union: {}, a: {}, b: {}",
                    reference_sketch.clone().union(query_sketch).cardinality(),
                    reference_sketch.cardinality(),
                    query_sketch.cardinality()
                );

                let numerator = 2.0 * similarity;
                let denominator = 1.0 + similarity;
                let fraction: f64 = numerator / denominator;
                let distance = 1.0f64 - fraction.powf(1.0 / kmer_length as f64);

                (reference_name.clone(), query_name.clone(), distance)
            },
        )
        .collect();
    Ok(results)
}

pub fn hll_distance(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    kmer_length: usize,
    query_names: Vec<String>,
    query_sketch_file: String,
) -> Result<Vec<(String, String, f64)>, Box<dyn Error>> {
    let ref_sketch_file = File::open(ref_sketch_file).expect("Failed to open file");
    let query_sketch_file = File::open(query_sketch_file).expect("Failed to open file");

    fn create_ull_map(
        sketch_file: File,
        names: &Vec<String>,
    ) -> Result<HashMap<String, (HyperLogLog<i64>, f64), Xxh3Builder>, std::io::Error> {
        let hasher = Xxh3Builder { seed: 93 };
        let mut sketches = HashMap::with_hasher(hasher);
        let reader = BufReader::new(sketch_file);

        // decompress sketches
        let mut decoder = Decoder::new(reader).expect("failed to create decompressor");
        for file in names {
            let hll = HyperLogLog::load(&mut decoder)?;
            let count = hll.len();
            sketches.insert(file.clone(), (hll, count));
        }
        Ok(sketches)
    }

    let ref_map = create_ull_map(ref_sketch_file, &reference_names).unwrap();
    let query_map = create_ull_map(query_sketch_file, &query_names).unwrap();

    // all pairs
    let pairs: Vec<(&str, &str)> = ref_map
        .keys()
        .flat_map(|k1| query_map.keys().map(move |k2| (k1.as_str(), k2.as_str())))
        .collect();

    let results = pairs
        .par_iter()
        .map(|&(reference_name, query_name)| {
            let a: f64 = ref_map[reference_name].1; // reference cardinality
            let b: f64 = query_map[query_name].1; // query cardinality
            let mut ref_hll = ref_map[reference_name].0.clone();
            let q_hll = &query_map[query_name].0;
            ref_hll.union(q_hll);
            let union_count = ref_hll.len();

            info!("Union: {}, a: {}, b: {}", union_count, a, b);

            let s = ((a + b - union_count) / union_count).max(0.0);
            // let numerator: f64 = 2.0 * s;
            // let denominator: f64 = 1.0 + s;
            // let fraction: f64 = numerator / denominator;
            // let distance: f64 = -fraction.ln() / (kmer_length as f64);
            let frac = 2.0 * s / (1.0 + s);
            let distance = 1.0f64 - frac.powf(1.0 / kmer_length as f64);
            (reference_name.to_string(), query_name.to_string(), distance)
        })
        .collect();

    Ok(results)
}

// Each file is processed to completion in its own task (no inner parallelism / no channels).
// This is simple, avoids stack overflows, and matches the “parallel by sample” request.

pub fn hmh_sketch(
    files: Vec<String>,
    kmer_length: usize,
    output_name: String,
    threads: u32,
    seed: u64,
) -> Result<(), Box<dyn Error>> {
    // build a Sketch per file in parallel
    let sketches: HashMap<String, Sketch> = files
        .par_iter()
        .map(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Failed to parse file");
            let mut global = Sketch::default();

            while let Some(rec_res) = reader.next() {
                if let Ok(rec) = rec_res {
                    let seq = filter_out_n(rec.seq().as_ref());
                    if seq.len() < kmer_length {
                        continue;
                    }
                    let kseq = KSeq::new(&seq, 2);

                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            global.add_bytes_with_seed(&(masked as u32).to_le_bytes(), seed);
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            global.add_bytes_with_seed(&(masked as u32).to_le_bytes(), seed);
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            global.add_bytes_with_seed(&masked.to_le_bytes(), seed);
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }
                }
            }

            (file_name.clone(), global)
        })
        .collect();

    // write names list
    let names: Vec<String> = sketches.keys().cloned().collect();
    to_writer_pretty(
        &File::create(format!("{}_files.json", output_name))?,
        &names,
    )?;

    // serialize sketches (compressed)
    let out_bin = format!("{}_sketches.bin", output_name);
    let mut writer = BufWriter::new(File::create(&out_bin)?);
    let mut encoder = Encoder::new(&mut writer, 3).expect("failed to create compression");
    encoder
        .multithread(threads)
        .expect("failed to multithread compressor");
    for name in &names {
        sketches[name].save(&mut encoder)?;
    }
    encoder.finish().expect("failed to compress");

    println!("Serialized sketches with hmh");
    Ok(())
}

pub fn hll_sketch(
    precision: u32,
    files: Vec<String>,
    kmer_length: usize,
    output_name: String,
    threads: u32,
    seed: u64,
) -> Result<(), Box<dyn Error>> {
    // build an HLL per file in parallel
    let hll_vec: Vec<HyperLogLog<i64>> = files
        .par_iter()
        .map(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
            let mut hll = HyperLogLog::<i64>::with_p(precision as u8);

            while let Some(res) = reader.next() {
                if let Ok(seqrec) = res {
                    let seq = filter_out_n(seqrec.seq().as_ref());
                    if seq.len() < kmer_length {
                        continue;
                    }
                    let kseq = KSeq::new(&seq, 2);

                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            hll.push_hash64(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            hll.push_hash64(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            hll.push_hash64(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }
                }
            }

            hll
        })
        .collect();

    // write compressed sketches
    let sketch_output =
        File::create(format!("{}_sketches.bin", &output_name)).expect("Failed to create file");
    let writer = BufWriter::new(sketch_output);

    let mut encoder = Encoder::new(writer, 3).expect("failed to create compression");
    encoder
        .multithread(threads)
        .expect("failed to multithread compressor");
    for hll in &hll_vec {
        hll.clone().save(&mut encoder).expect("Failed to save HLL");
    }
    encoder.finish().expect("failed to compress");

    // names & params
    to_writer_pretty(
        &File::create(format!("{}_files.json", &output_name))?,
        &files,
    )?;

    println!("Serialized sketches with hll");
    Ok(())
}

pub fn ull_sketch(
    precision: u32,
    files: Vec<String>,
    kmer_length: usize,
    output_name: String,
    threads: u32,
    seed: u64,
) -> Result<(), Box<dyn Error>> {
    // build a ULL per file in parallel
    let ull_vec: Vec<UltraLogLog> = files
        .par_iter()
        .map(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
            let mut ull = UltraLogLog::new(precision).expect("failed to create ULL");

            while let Some(result) = reader.next() {
                if let Ok(seqrec) = result {
                    let seq = filter_out_n(seqrec.seq().as_ref());
                    if seq.len() < kmer_length {
                        continue;
                    }
                    let kseq = KSeq::new(&seq, 2);

                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            ull.add(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked =
                                mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            ull.add(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            ull.add(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
                        }
                    } else {
                        panic!("k-mer length must be 1–32, k=15 is not supported");
                    }
                }
            }

            ull
        })
        .collect();

    // write compressed sketches
    let sketch_output =
        File::create(format!("{}_sketches.bin", &output_name)).expect("Failed to create file");
    let writer = BufWriter::new(sketch_output);

    let mut encoder = Encoder::new(writer, 3).expect("failed to create compression");
    encoder
        .multithread(threads)
        .expect("failed to multithread compressor");
    for ull in &ull_vec {
        ull.save(&mut encoder).expect("Failed to save UltraLogLog");
    }
    encoder.finish().expect("failed to compress");

    // names
    to_writer_pretty(
        &File::create(format!("{}_files.json", &output_name))?,
        &files,
    )?;

    println!("Serialized sketches with ull");
    Ok(())
}
