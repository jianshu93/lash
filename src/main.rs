use clap::{Arg, ArgAction, Command};
use crossbeam_channel::bounded;
use hyperminhash::Sketch;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
    thread,
};

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

fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_default_env().init();

    let matches = Command::new("Genome Sketching via HyperMinhash")
        .version("0.1.0")
        .about("Fast and Memory Efficient Genome/Metagenome Sketching via HyperMinhash")
        .arg(
            Arg::new("query_files")
                .short('q')
                .long("query_file")
                .help("List of query genome files, one per line with .gz support")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("reference_files")
                .short('r')
                .long("ref_file")
                .help("List of reference genome files, one per line with .gz support")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("kmer_length")
                .short('k')
                .long("kmer")
                .help("Kmer length to use for sketching")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads to use in parallel")
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output_file")
                .short('o')
                .long("output")
                .required(true)
                .action(ArgAction::Set),
        )
        .get_matches();

    let query_files_list = matches.get_one::<String>("query_files").unwrap();
    let reference_files_list = matches.get_one::<String>("reference_files").unwrap();
    let kmer_length: usize = *matches.get_one::<usize>("kmer_length").unwrap();
    let output_file = matches.get_one::<String>("output_file").unwrap();
    let threads = *matches.get_one::<usize>("threads").unwrap();
    // Build global threadpool
    ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    
    let read_paths = |p: &str| -> Vec<String> {
        BufReader::new(File::open(p).unwrap())
            .lines()
            .filter_map(Result::ok)
            .filter(|l| !l.trim().is_empty())
            .collect()
    };
    let query_files = read_paths(query_files_list);
    let reference_files = read_paths(reference_files_list);

    let query_sketches: HashMap<String, Sketch> = query_files
        .par_iter()
        .map(|file_name| {
            let file_name = file_name.clone();
            let mut reader = parse_fastx_file(&file_name).expect("Failed to parse file");
            let (sender, receiver) = bounded(10);

            let reader_thread = thread::spawn(move || {
                let mut batch = Vec::new();
                while let Some(result) = reader.next() {
                    if let Ok(seqrec) = result {
                        batch.push(seqrec.seq().to_vec());
                        if batch.len() == 5000 {
                            if sender.send(batch.clone()).is_err() { break }
                            batch.clear();
                        }
                    }
                }
                if !batch.is_empty() { let _ = sender.send(batch); }
            });

            let mut global_sketch = Box::new(Sketch::default());

            for batch in receiver {
                let local_sketch = batch.par_iter().map(|seq| {
                    let seq_vec = filter_out_n(seq);
                    if seq_vec.is_empty() { return Box::new(Sketch::default()); }
                    let mut sketch = Box::new(Sketch::default());
                    let kseq = KSeq::new(&seq_vec, 2);
                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
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
            reader_thread.join().unwrap();
            (file_name, *global_sketch)
        })
        .collect();

    let reference_sketches: HashMap<String, Sketch> = reference_files
        .par_iter()
        .map(|file_name| {
            let file_name = file_name.clone();
            let mut reader = parse_fastx_file(&file_name).expect("Failed to parse file");
            let (sender, receiver) = bounded(10);
            let reader_thread = thread::spawn(move || {
                let mut batch = Vec::new();
                while let Some(result) = reader.next() {
                    if let Ok(seqrec) = result {
                        batch.push(seqrec.seq().to_vec());
                        if batch.len() == 5000 {
                            if sender.send(batch.clone()).is_err() { break }
                            batch.clear();
                        }
                    }
                }
                if !batch.is_empty() { let _ = sender.send(batch); }
            });

            let mut global_sketch = Box::new(Sketch::default());
            for batch in receiver {
                let local_sketch = batch.par_iter().map(|seq| {
                    let seq_vec = filter_out_n(seq);
                    if seq_vec.is_empty() { return Box::new(Sketch::default()); }
                    let mut sketch = Box::new(Sketch::default());
                    let kseq = KSeq::new(&seq_vec, 2);
                    if kmer_length <= 14 {
                        let mut it = KmerSeqIterator::<Kmer32bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length == 16 {
                        let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value() as u64, kmer_length);
                            sketch.add_bytes(&(masked as u32).to_le_bytes());
                        }
                    } else if kmer_length <= 32 {
                        let mut it = KmerSeqIterator::<Kmer64bit>::new(kmer_length as u8, &kseq);
                        while let Some(km) = it.next() {
                            let canon = km.min(km.reverse_complement());
                            let masked = mask_bits(canon.get_compressed_value(), kmer_length);
                            sketch.add_bytes(&masked.to_le_bytes());
                        }
                    } else {
                        panic!("k-mer length must be 1–32");
                    }
                    sketch
                })
                .reduce(|| Box::new(Sketch::default()), |mut a, b| { a.union(&b); a });
                global_sketch.union(&local_sketch);
            }
            reader_thread.join().unwrap();
            (file_name, *global_sketch)
        })
        .collect();

    let pairs: Vec<(&String,&Sketch,&String,&Sketch)> = query_sketches
        .iter()
        .flat_map(|(qn,qs)| reference_sketches.iter().map(move |(rn,rs)| (qn,qs,rn,rs)))
        .collect();

    let results: Vec<(String,String,f64)> = pairs
        .par_iter()
        .map(|(qn,qs,rn,rs)| {
            let sim  = qs.similarity(rs).max(f64::EPSILON);
            let dist = -(2.0*sim/(1.0+sim)).ln() / (kmer_length as f64);
            (qn.to_string(), rn.to_string(), dist)   // ← own the strings
        })
        .collect();

    let mut out = File::create(output_file)?;
    writeln!(out,"Query\tReference\tDistance")?;
    for (q,r,d) in results {
        let qb = Path::new(&q).file_name().unwrap().to_string_lossy();
        let rb = Path::new(&r).file_name().unwrap().to_string_lossy();
        writeln!(out,"{q}\t{r}\t{:.6}", if qb==rb { 0.0 } else { d })?;
    }
    Ok(())
}
