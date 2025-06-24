//! Genome/Metagenome sketching with **HyperMinHash**

use clap::{Arg, ArgAction, Command};
use crossbeam_channel::bounded;
use hyperminhash::Sketch;
use needletail::parse_fastx_file;
use rayon::iter::ParallelBridge;
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
use std::iter;

use chrono::Local;

// k-mer helpers ­──────────────────────────────────────────────────────────────────────
use kmerutils::base::{
    alphabet::Alphabet2b,
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as SequenceStruct,
    CompressedKmerT, Kmer16b32bit, Kmer32bit, Kmer64bit,
};
use kmerutils::base::KmerT;

/// Convert ASCII sequence → `SequenceStruct` (2-bit alphabet)
fn ascii_to_seq(bases: &[u8]) -> SequenceStruct {
    let alpha = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alpha);
    seq
}

/// Keep only the lowest 2 k bits (canonicalised integer already)
fn mask_bits(v: u64, k: usize) -> u64 {
    let b = 2 * k as u32;
    if b == 64 { v } else { v & ((1u64 << b) - 1) }
}

fn sketch_one_file_hmh(path: &str, k: usize) -> Sketch {
    use std::iter;

    const BATCH: usize = 5_000;

    // 1.  Move the FASTX reader into a generator-style iterator
    let batch_iter = {
        let p = path.to_owned();
        let mut rdr = parse_fastx_file(&p).expect("open FASTX");

        iter::from_fn(move || {
            let mut batch = Vec::with_capacity(BATCH);
            while let Some(rec) = rdr.next() {
                if let Ok(seqrec) = rec {
                    if seqrec.seq().len() > k {
                        batch.push(seqrec.seq().to_vec());
                        if batch.len() == BATCH { break }
                    }
                }
            }
            if batch.is_empty() { None } else { Some(batch) }
        })
    };

    // 2.  Stream those batches in parallel with `par_bridge`
    batch_iter
        .par_bridge()
        .map(|batch| {
            let mut sk = Sketch::default();
            for seq in batch {
                let seq_vec = ascii_to_seq(&seq);

                if k <= 14 {
                    let mut it = KmerSeqIterator::<Kmer32bit>::new(k as u8, &seq_vec);
                    while let Some(km) = it.next() {
                        let masked = mask_bits(
                            km.min(km.reverse_complement()).get_compressed_value() as u64, k);
                        sk.add_bytes(&(masked as u32).to_le_bytes());
                    }
                } else if k == 16 {
                    let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &seq_vec);
                    while let Some(km) = it.next() {
                        let masked = mask_bits(
                            km.min(km.reverse_complement()).get_compressed_value() as u64, k);
                        sk.add_bytes(&(masked as u32).to_le_bytes());
                    }
                } else {
                    let mut it = KmerSeqIterator::<Kmer64bit>::new(k as u8, &seq_vec);
                    while let Some(km) = it.next() {
                        let masked = mask_bits(
                            km.min(km.reverse_complement()).get_compressed_value(), k);
                        sk.add_bytes(&masked.to_le_bytes());
                    }
                }
            }
            sk
        })
        .reduce(|| Sketch::default(), |mut a, b| { a.union(&b); a })
}

// ─────────────────────────────────────────────────────────────────────────────────────
// Main program
// ------------------------------------------------------------------------------------
fn main() -> Result<(), Box<dyn Error>> {
    env_logger::Builder::from_default_env().init();

    // ---------- CLI ----------
    let matches = Command::new("Genome Sketching via HyperMinHash")
        .version("0.2.0-bridge")
        .about("Fast and memory-efficient genome / metagenome sketching (streaming, Rayon)")
        .arg(
            Arg::new("query_files")
                .short('q')
                .long("query_file")
                .help("Text file: one query FASTA/Q path per line")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("reference_files")
                .short('r')
                .long("ref_file")
                .help("Text file: one reference FASTA/Q path per line")
                .required(true)
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("kmer_length")
                .short('k')
                .long("kmer")
                .help("k-mer size (1–32, except 15)")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Rayon thread count")
                .default_value("1")
                .value_parser(clap::value_parser!(usize))
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output_file")
                .short('o')
                .long("output")
                .help("TSV output path")
                .required(true)
                .action(ArgAction::Set),
        )
        .get_matches();

    let query_list   = matches.get_one::<String>("query_files").unwrap();
    let reference_list = matches.get_one::<String>("reference_files").unwrap();
    let kmer_len: usize = *matches.get_one::<usize>("kmer_length").unwrap();
    let out_path    = matches.get_one::<String>("output_file").unwrap();
    let threads     = *matches.get_one::<usize>("threads").unwrap();

    ThreadPoolBuilder::new().num_threads(threads).build_global()?;

    // ---------- read path lists ----------
    let load_paths = |p: &str| -> Vec<String> {
        BufReader::new(File::open(p).unwrap())
            .lines()
            .filter_map(Result::ok)
            .filter(|l| !l.trim().is_empty())
            .collect()
    };
    let query_files     = load_paths(query_list);
    let reference_files = load_paths(reference_list);

    println!("{}  sketching queries …", Local::now());
    let query_sketches: HashMap<String, Sketch> = query_files
        .par_iter()
        .map(|p| (p.clone(), sketch_one_file_hmh(p, kmer_len)))
        .collect();

    println!("{}  sketching references …", Local::now());
    let reference_sketches: HashMap<String, Sketch> = reference_files
        .par_iter()
        .map(|p| (p.clone(), sketch_one_file_hmh(p, kmer_len)))
        .collect();

    // ---------- pairwise distances ----------
    let pairs: Vec<(&String,&Sketch,&String,&Sketch)> = query_sketches
        .iter()
        .flat_map(|(qn,qs)| reference_sketches
            .iter()
            .map(move |(rn,rs)| (qn,qs,rn,rs)))
        .collect();

    let results: Vec<(String,String,f64)> = pairs
        .par_iter()
        .map(|(qn,qs,rn,rs)| {
            let sim  = qs.similarity(rs);
            let dist = (2.0 * sim / (1.0 + sim)).powf(1.0 / kmer_len as f64);
            (qn.to_string(), rn.to_string(), dist)
        })
        .collect();

    // ---------- write output ----------
    let mut out = File::create(out_path)?;
    writeln!(out, "Query\tReference\tDistance")?;
    for (q,r,d) in results {
        let qb = Path::new(&q).file_name().unwrap().to_string_lossy();
        let rb = Path::new(&r).file_name().unwrap().to_string_lossy();
        writeln!(out, "{q}\t{r}\t{:.6}", if qb==rb { 0.0 } else { d })?;
    }
    println!("{}  done.", Local::now());
    Ok(())
}