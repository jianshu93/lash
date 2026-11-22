use clap::{Arg, ArgAction, Command};
use exaloglog::ExaLogLog;
use streaming_algorithms::HyperLogLog;
use kmerutils::base::KmerT;
use kmerutils::base::{
    kmergenerator::{KmerSeqIterator, KmerSeqIteratorT},
    sequence::Sequence as KSeq,
    CompressedKmerT, Kmer16b32bit};
use rand::seq::index;
use std::error::Error;
use std::sync::{Mutex};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use ultraloglog::{Estimator, MaximumLikelihoodEstimator, UltraLogLog};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use xxhash_rust::xxh3::xxh3_64_with_seed;
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

    let matches = Command::new("MVP on genomes")
        .arg(
            Arg::new("file")
            .short('f')
            .long("file")
            .help("One file containing list of FASTA/FASTQ files, one per line")
            .required(true)
            .action(ArgAction::Set),
        )
        .arg(
            Arg::new("seeds")
            .short('s')
            .long("seeds")
            .help("Number of random seeds to simulate")
            .default_value("10")
            .value_parser(clap::value_parser!(usize))
            .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output")
            .short('o')
            .long("output")
            .help("output file name")
            .default_value("log")
            .action(ArgAction::Set),
        )
        .arg(
            Arg::new("p")
            .short('p')
            .long("p")
            .help("precision")
            .value_parser(clap::value_parser!(usize))
            .action(ArgAction::Set),
        ).get_matches();

    let input_file = matches.get_one::<String>("file").expect("required");
    let output_file = matches.get_one::<String>("output").expect("required");
    let num_seeds: usize = *matches.get_one::<usize>("seeds").expect("required");
    let p = *matches.get_one::<usize>("p").expect("required");

    if num_seeds > 100 {
        panic!("Max 100 seeds each time")
    }

    let mut rng = rand::rng();
    let seeds: Vec<usize> = index::sample(&mut rng, 101, num_seeds).into_vec();

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

    let mut writer = BufWriter::new(File::create(output_file)?);
    

    // keeps track of true unique count
    let seen: Vec<AtomicU64> = (0..(1usize << 26)).map(|_| AtomicU64::new(0)).collect();

    // vectors for each sketchs
    let mut hlls: Vec<Mutex<HyperLogLog<i64>>> = Vec::with_capacity(num_seeds);
    let mut ulls: Vec<Mutex<UltraLogLog>> = Vec::with_capacity(num_seeds);
    let mut ells: Vec<Mutex<ExaLogLog>> = Vec::with_capacity(num_seeds);

    // fill the vectors: one sketch of each type per seed
    for _ in 0..num_seeds {
        hlls.push(Mutex::new(HyperLogLog::<i64>::with_p(p as u8)));
        ulls.push(Mutex::new(UltraLogLog::new(p as u32).expect("error creating ull")));
        ells.push(Mutex::new(ExaLogLog::new(2, 20, p as u32).expect("error creating ell")));
    }

    const CAPACITY: usize = 1 << 20;

    let flush = |buf: &mut Vec<u64>| {
        for (j, seed) in seeds.iter().enumerate() {
            // hash outside the locks
            let hashes: Vec<u64> = buf
                .iter()
                .map(|&m| xxh3_64_with_seed(&m.to_le_bytes(), *seed as u64))
                .collect();

            let mut h = hlls[j].lock().unwrap();
            for &x in &hashes { h.push_hash64(x); }
        
            let mut u = ulls[j].lock().unwrap();
            for &x in &hashes { u.add(x); }
        
            let mut e = ells[j].lock().unwrap();
            for &x in &hashes { e.add_hash(x); }
        }
        buf.clear();
    };

    files.par_iter().for_each(|file_name| {
        let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
        let mut buf: Vec<u64> = Vec::with_capacity(CAPACITY);
        while let Some(res) = reader.next() {
            if let Ok(seqrec) = res {
                let seq = filter_out_n(seqrec.seq().as_ref());
                if seq.len() < 16 {
                    continue;
                }

                let kseq = KSeq::new(&seq, 2);
                let mut it = KmerSeqIterator::<Kmer16b32bit>::new(16, &kseq);
                while let Some(km) = it.next() {
                    let canon = km.min(km.reverse_complement());
                    let masked = mask_bits(
                        canon.get_compressed_value() as u64,
                        16,
                    );
                    let m = masked as usize;
                    seen[m >> 6].fetch_or(1u64 << (m & 63), Ordering::Relaxed);
                    buf.push(masked);
                    if buf.len() >= CAPACITY { flush(&mut buf); }
                }
            }
        }
        flush(&mut buf); // once per file
    });
    
    // write files
    let mut sketch_writer = BufWriter::new(File::create("hll_sketch.bin")?);
    hlls[0].lock().unwrap().save(&mut sketch_writer)?;
    sketch_writer.flush()?;

    let mut sketch_writer = BufWriter::new(File::create("ull_sketch.bin")?);
    ulls[0].lock().unwrap().save(&mut sketch_writer)?;
    sketch_writer.flush()?;

    let mut sketch_writer = BufWriter::new(File::create("ell_sketch.bin")?);
    ells[0].lock().unwrap().save(&mut sketch_writer)?;
    sketch_writer.flush()?;

    let mut hll_cards = Vec::<f64>::new();
    let mut ull_cards = Vec::<f64>::new();
    let mut ell_cards = Vec::<f64>::new();

    for i in 0..num_seeds {
        writeln!(writer, "-------------- seed number {} --------------", i).expect("error writing to output");
        let hll_card = hlls[i].lock().unwrap().len();
        writeln!(writer, "hll cardinality: {}", hll_card).expect("error writing hll result");
        hll_cards.push(hll_card);

        let ull_card = MaximumLikelihoodEstimator.estimate(&ulls[i].lock().unwrap());
        writeln!(writer, "ull cardinality: {}", ull_card).expect("error writing ull result");
        ull_cards.push(ull_card);

        let ell_card = ells[i].lock().unwrap().estimate();
        writeln!(writer, "ell cardinality: {}", ell_card).expect("error writing ell result");
        ell_cards.push(ell_card);
    }
    
    // calculate true cardinality
    let true_unique: u64 = seen
    .iter()
    .map(|w| w.load(Ordering::Relaxed).count_ones() as u64)
    .sum();

    writeln!(writer, "**************** True cardinality: {} ****************", true_unique).expect("error writing to file");

    fn mse(estimates: Vec<f64>, truth: f64) -> f64 {
        let n = estimates.len() as f64;
        estimates.iter().map(|e| ((e - truth) / truth).powi(2)).sum::<f64>() / n
    }

    let var_hll = mse(hll_cards, true_unique as f64);
    let var_ull = mse(ull_cards, true_unique as f64);
    let var_ell = mse(ell_cards, true_unique as f64);

    writeln!(writer, "**************** MSE results ****************").expect("error writing to log file");
    writeln!(writer, "hll mse: {}", var_hll).expect("error writing hll mse");
    writeln!(writer, "ull mse: {}", var_ull).expect("error writing ull mse");
    writeln!(writer, "ell mse: {}", var_ell).expect("error writing ell mse");

    writeln!(writer, "**************** MVP results ****************").expect("error writing mvp");
    writeln!(writer, "hll mvp: {}", var_hll * 2i32.pow(p as u32) as f64 * 6.0).expect("error writing hll mvp");
    writeln!(writer, "ull mvp: {}", var_ull * 2i32.pow(p as u32) as f64* 8.0).expect("error writing ull mvp");
    writeln!(writer, "ell mvp: {}", var_ell * 2i32.pow(p as u32) as f64 * 28.0).expect("error writing ell mvp");
    Ok(())
}