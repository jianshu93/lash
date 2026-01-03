// src/utils.rs

use hashbrown::HashMap;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::error::Error;
use crate::hasher::Xxh3Builder;
use num_traits::{Float};
use std::fs::File;
use std::io::{BufReader, BufWriter};

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
use kmerutils::aautils::kmeraa::{KmerAA32bit, KmerAA64bit};
use ultraloglog::{Estimator, MaximumLikelihoodEstimator, UltraLogLog};

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
pub fn hmh_distance<F, T: Float>(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    query_names: Vec<String>,
    query_sketch_file: String,
    create_matrix: bool,
    same_files: bool,
    emit: F,
) -> std::io::Result<()>
where F: Fn(Vec<(&String, &String, T)>) + Send + Sync
{
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

    // send column names if printing matrix
    let mut file_idx: HashMap<&String, usize> = HashMap::new();
    
    if same_files || create_matrix {
        let mut columns: Vec<(&String, &String, T)> = Vec::new();
        let blank_str = &"".to_string();
        for (i, q_name) in query_sketches.keys().enumerate() {
            // empty r_name string signals printing columns
            if create_matrix {
                columns.push((blank_str, q_name, T::one()));
            }
            if same_files {
                file_idx.insert(q_name, i);
            }
        }
        if create_matrix {
            emit(columns);
        }
    }

    // loop through reference sketches (i)
    reference_sketches.par_iter().for_each(|(ref_name, _)| {
        let ref_sketch = reference_sketches[ref_name];
        let mut ref_row: Vec<(&String, &String, T)> = Vec::new();

        // loop through query sketches (j)
        for q_name in query_sketches.keys() {

            // for triangular matrix purposes
            if same_files && file_idx[q_name] > file_idx[ref_name] {
                continue;
            }
            let q_sketch = query_sketches[q_name];

            // calculate distance (i, j, d)
            let similarity = q_sketch.similarity(ref_sketch).max(0.0);
            let numerator = 2.0 * similarity;
            let denominator = 1.0 + similarity;
            let fraction = numerator / denominator;

            info!(
                "Union: {}, a: {}, b: {}",
                ref_sketch.clone().union(q_sketch).cardinality(),
                ref_sketch.cardinality(),
                q_sketch.cardinality()
            );

            let frac_t: T = T::from(fraction).expect("failed to convert to f64 or f32");
            ref_row.push((ref_name, q_name, frac_t));
        }
        emit(ref_row);
    });

    Ok(())
    
}

pub fn ull_distance <F, T: Float>(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    query_names: Vec<String>,
    query_sketch_file: String,
    estimator: String,
    create_matrix: bool,
    same_files: bool,
    emit: F,
)-> std::io::Result<()>
where 
F: Fn(Vec<(&String, &String, T)>) + Send + Sync {

    let ref_sketch_file = File::open(ref_sketch_file).expect("Failed to open file");
    let query_sketch_file = File::open(query_sketch_file).expect("Failed to open file");

    fn create_ull_map(
        sketch_file: File,
        names: &Vec<String>,
        estimator: &String,
    ) -> Result<HashMap<String, (UltraLogLog, f64), Xxh3Builder>, std::io::Error>
    {
        let hasher = Xxh3Builder { seed: 93 };
        let mut sketches = HashMap::with_hasher(hasher);
        let reader = BufReader::new(sketch_file);
        let mut decoder = Decoder::new(reader).expect("failed to create decompressor");
        for file in names {
            let ull = UltraLogLog::load(&mut decoder)?;
            let c: f64 = match estimator.as_str() {
                "fgra" => ull.get_distinct_count_estimate(),
                "ml" => MaximumLikelihoodEstimator.estimate(&ull),
                _ => panic!("estimator needs to be either fgra or ml"),
            };
            sketches.insert(file.clone(), (ull, c));
        }
        Ok(sketches)
    }

    let ref_map =
        create_ull_map(ref_sketch_file, &reference_names, &estimator).unwrap();
    let query_map =
        create_ull_map(query_sketch_file, &query_names, &estimator).unwrap();

    let mut file_idx: HashMap<&String, usize> = HashMap::new();
    if same_files || create_matrix {
        let mut columns: Vec<(&String, &String, T)> = Vec::new();
        let blank = "".to_string();
        for (i, q_name) in query_map.keys().enumerate() {
            // empty r_name string signals printing columns
            if create_matrix {
                columns.push((&blank, q_name, T::one()));
            }
            // used for redundant distances
            if same_files {
                file_idx.insert(q_name, i);
            }
        }
        if create_matrix {
            emit(columns);
        }
    }

    ref_map.par_iter().for_each(|(ref_name, _)| {
        // print ref name on the new line and on the left if matrix
        let a: f64 = ref_map[ref_name].1;

        let mut ref_list: Vec<(&String, &String, T)> = Vec::new();
        // loop through query sketches (j)
        for qry_name in query_map.keys() {
            // for redundant distances
            if same_files && file_idx[qry_name] > file_idx[ref_name] {
                continue;
            }
            let b: f64 = query_map[qry_name].1;
            let union_ull =
                UltraLogLog::merge(&ref_map[ref_name].0, &query_map[qry_name].0)
                    .expect("failed to merge sketches");

            //let union_count = union_ull.get_distinct_count_estimate();
            let union_count: f64 = match estimator.as_str() {
                "fgra" => union_ull.get_distinct_count_estimate(),
                "ml" => MaximumLikelihoodEstimator.estimate(&union_ull),
                _ => panic!("estimator needs to be either fgra or ml"),
            };

            info!("Union: {}, a: {}, b: {}", union_count, a, b);

            let similarity = (a + b - union_count) / union_count;
            let s = if similarity < 0.0 { 0.0 } else { similarity };
            let frac_f64 = 2.0 * s / (1.0 + s);

            let frac_t: T = T::from(frac_f64)
                .expect("failed to convert f64 to T");

            ref_list.push((ref_name, qry_name, frac_t));
        }

        // emit a vector with same ref file
        emit(ref_list);
    });
    
    Ok(())
}

pub fn hll_distance<F, T: Float>(
    reference_names: Vec<String>,
    ref_sketch_file: String,
    query_names: Vec<String>,
    query_sketch_file: String,
    create_matrix: bool,
    same_files: bool,
    emit: F,
) -> std::io::Result<()>
where F: Fn(Vec<(&String, &String, T)>) + Send + Sync {
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

    let mut file_idx: HashMap<&String, usize> = HashMap::new();
    if same_files || create_matrix {
        let mut columns: Vec<(&String, &String, T)> = Vec::new();
        let blank = &"".to_string();
        for (i, q_name) in query_map.keys().enumerate() {
            // empty r_name string signals printing columns
            if create_matrix {
                columns.push((blank, q_name, T::one()));
            }
            if same_files {
                file_idx.insert(q_name, i);
            }
        }
        if create_matrix {
            emit(columns);
        }
    }

    ref_map.par_iter().for_each(|(ref_name, _) | {
        let a: f64 = ref_map[ref_name].1;
        let mut row: Vec<(&String, &String, T)> = Vec::new();

        // loop through query sketches (j)
        for qry_name in query_map.keys() {

            // for triangular matrix
            if same_files && file_idx[qry_name] > file_idx[ref_name] {
                continue;
            }
             // reference cardinality
            let b: f64 = query_map[qry_name].1; // query cardinality
            let mut ref_hll = ref_map[ref_name].0.clone();
            let q_hll = &query_map[qry_name].0;
            ref_hll.union(q_hll);
            let union_count = ref_hll.len();

            info!("Union: {}, a: {}, b: {}", union_count, a, b);

            let s = ((a + b - union_count) / union_count).max(0.0);
            let frac = 2.0 * s / (1.0 + s);
            let frac_t: T = T::from(frac).expect("failed to convert f64 to T");
            // let distance = 1.0f64 - frac.powf(1.0 / kmer_length as f64);
            
            row.push((ref_name, qry_name, frac_t));
        }
        emit(row);
    });
    
    Ok(())
}


// sketch trait shared by HMH, ULL, and HLL
pub trait KmerSketch: Send {
    /// Create a new sketch
    fn new(precision: Option<u32>) -> Self;

    /// Add a canonical masked k-mer
    fn add_kmer(&mut self, masked: u64, seed: u64);

    /// Serialize
    fn save<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>>;
}


// HMH sketching
impl KmerSketch for Sketch {
    fn new(_: Option<u32>) -> Self {
        Sketch::default()
    }

    fn add_kmer(&mut self, masked: u64, seed: u64) {
        // HMH uses bytes + seed
        self.add_bytes_with_seed(&(masked as u32).to_le_bytes(), seed);
    }

    fn save<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        Ok(self.save(writer)?)
    }
}

// sketching for HyperLogLog
impl KmerSketch for HyperLogLog<i64> {
    fn new(precision: Option<u32>) -> Self {
        HyperLogLog::with_p(precision.expect("HLL needs precision") as u8)
    }

    fn add_kmer(&mut self, masked: u64, seed: u64) {
        self.push_hash64(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
    }

    fn save<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        Ok(self.save(writer)?)
    }
}

// sketching for UltraLogLog
impl KmerSketch for UltraLogLog {
    fn new(precision: Option<u32>) -> Self {
        UltraLogLog::new(precision.expect("ULL needs precision"))
            .expect("failed to create ULL")
    }

    fn add_kmer(&mut self, masked: u64, seed: u64) {
        self.add(xxh3_64_with_seed(&masked.to_le_bytes(), seed));
    }

    fn save<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Box<dyn Error>> {
        Ok(self.save(writer)?)
    }
}

// general sketching function
// Each file is processed to completion in its own task (no inner parallelism / no channels).
// This is simple, avoids stack overflows, and matches the “parallel by sample” request.
pub fn sketch_files <S: KmerSketch> (
    precision: Option<u32>,
    files: Vec<String>,
    kmer_length: usize,
    output_name: String,
    threads: u32,
    seed: u64,
    aa: bool
) -> Result<(), Box<dyn Error>> {

    let sketches: Vec<S> = if ! aa { // genome sketching
        files
        .par_iter()
        .map(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
            let mut sketch = S::new(precision);

            while let Some(res) = reader.next() {
                if let Ok(seqrec) = res {
                    let seq = filter_out_n(seqrec.seq().as_ref());
                    if seq.len() < kmer_length {
                        continue;
                    }

                    let kseq = KSeq::new(&seq, 2);

                    if kmer_length <= 14 {
                            let mut it =
                                KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked = mask_bits(
                                    canon.get_compressed_value() as u64,
                                    kmer_length,
                                );
                                sketch.add_kmer(masked, seed);
                            }
                        }
                    else if kmer_length == 16 {
                            let mut it =
                                KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked = mask_bits(
                                    canon.get_compressed_value() as u64,
                                    kmer_length,
                                );
                                sketch.add_kmer(masked, seed);
                            }
                        }
                    else if kmer_length <= 32 {
                            let mut it =
                                KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                            while let Some(km) = it.next() {
                                let canon = km.min(km.reverse_complement());
                                let masked =
                                    mask_bits(canon.get_compressed_value(), kmer_length);
                                sketch.add_kmer(masked, seed);
                            }
                        }
                    else {
                        panic!("k-mer length must be 1–32");
                    } 
                    }
                
            }

            sketch
        })
        .collect()
    }
    else { // amino acid sketching
        files
        .par_iter()
        .map(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
            let mut sketch = S::new(precision);

            while let Some(res) = reader.next() {
                if let Ok(seqrec) = res {
                    let seq = filter_out_n(seqrec.seq().as_ref());
                    if seq.len() < kmer_length {
                        continue;
                    }

                    let kseq = KSeq::new(&seq, 5); // 5 bases for aa

                    if kmer_length <= 6 {
                        let mut it =
                            KmerSeqIterator::<KmerAA32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            //let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(
                                km.get_compressed_value() as u64,
                                kmer_length,
                            );
                            sketch.add_kmer(masked, seed);
                        }
                    }
                    else if kmer_length <= 12 {
                        let mut it =
                            KmerSeqIterator::<KmerAA64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(
                                canon.get_compressed_value() as u64,
                                kmer_length,
                            );
                            sketch.add_kmer(masked, seed);
                        }
                    }
                    else {
                        panic!("k-mer length for amino acid must be 1–12");
                    } 
                }
                
            }

            sketch
        })
        .collect()
    };
    

    // write sketches
    let writer = BufWriter::new(File::create(format!("{}_sketches.bin", output_name))?);
    let mut encoder = Encoder::new(writer, 3)?;
    encoder.multithread(threads)?;

    for sketch in sketches {
        sketch.save(&mut encoder)?;
    }
    encoder.finish()?;

    // write names
    to_writer_pretty(
        &File::create(format!("{}_files.json", output_name))?,
        &files,
    )?;

    Ok(())
}
