
use std::env;
use kmerutils::base::KmerT;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit};

use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};

pub fn mask_bits(v: u64, k: usize) -> u64 {
    let b = 2 * k as u32;
    if b == 64 {
        v
    } else {
        v & ((1u64 << b) - 1)
    }
}

pub fn filter_out_n(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &c in seq {
        if matches!(c, b'A' | b'C' | b'T' | b'G') {
            out.push(c);
        }
    }
    out
}


fn main() -> Result<(), Box<dyn Error>> {
    // input file path here
    let args: Vec<String> = env::args().collect();

    if args.len() > 2 {
        eprintln!("Use only one argument (filepath to file that lists genomes)")
    }

    let input_file = &args[1];
    let threads = num_cpus::get();

    rayon::ThreadPoolBuilder::new()
        .num_threads(threads.max(1))
        .build_global()
        .unwrap();
    
    // create vector of files to sketch
    let files: Vec<String> = {
        let f = File::open(input_file).expect("error opening file");
        BufReader::new(f)
            .lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.trim().is_empty())
            .collect()
    };

    // 2^32 possible 16-mers / 64 bits per word = 2^26 words (512 MB)
    // ground truth cardinality
    let seen: Vec<AtomicU64> = (0..(1usize << 26)).map(|_| AtomicU64::new(0)).collect();

    // sketch hll, ell, and ull into one big sketch
    files
        .par_iter()
        .for_each(|file_name| {
            let mut reader = parse_fastx_file(file_name).expect("Invalid input file");

            // looping through each sequence in file. k=16
            while let Some(res) = reader.next() {
                if let Ok(seqrec) = res {
                    let seq = filter_out_n(seqrec.seq().as_ref());
                    let kseq = KSeq::new(&seq, 2);
                    
                    let mut it =
                        KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                    while let Some(km) = it.next() {
                        let canon = km.min(km.reverse_complement());
                        let masked = mask_bits(
                            canon.get_compressed_value() as u64,
                            16,
                        );

                        let m = masked as usize; // must be < 2^32
                        seen[m >> 6].fetch_or(1u64 << (m & 63), Ordering::Relaxed);
                    }
                }
                
            }
        });

    
    let true_unique: u64 = seen
        .iter()
        .map(|w| w.load(Ordering::Relaxed).count_ones() as u64)
        .sum();

    println!("true unique k-mers: {true_unique}");
    Ok(())
}