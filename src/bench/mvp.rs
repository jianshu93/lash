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
            .default_value("20")
            .value_parser(clap::value_parser!(usize))
            .action(ArgAction::Set),
        )
        .arg(
            Arg::new("cardinality")
            .short('c')
            .long("cardinality")
            .help("True cardinality of genome set")
            .required(true)
            .value_parser(clap::value_parser!(f64))
            .action(ArgAction::Set),
        )
        .arg(
            Arg::new("output")
            .short('o')
            .long("output")
            .help("output file name")
            .default_value("log")
            .action(ArgAction::Set),
        ).get_matches();

    let input_file = matches.get_one::<String>("file").expect("required");
    let output_file = matches.get_one::<String>("output").expect("required");
    let cardinality = *matches.get_one::<f64>("cardinality").expect("required");
    let num_seeds: usize = *matches.get_one::<usize>("seeds").expect("required");

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
    let mut hll_errors = Vec::<f64>::with_capacity(num_seeds);
    let mut ull_errors = Vec::<f64>::with_capacity(num_seeds);
    let mut ell_errors = Vec::<f64>::with_capacity(num_seeds);

    let mut i = 0; // keep track of seed index, save sketches of the last seed only

    // calculate cardinality and squared error for each seed
    for seed in seeds {
        writeln!(writer, "---------------- seed {} ----------------", seed).expect("error writing to file");
        let hll = Mutex::new(HyperLogLog::<i64>::with_p(10));
        let ull = Mutex::new(UltraLogLog::new(10).expect("error creating ull"));
        let ell = Mutex::new(ExaLogLog::new(2, 20, 10).expect("error creating ell"));


        const CAPACITY: usize = 1 << 20; // max capacity of vector holding hashed values

        // adds all entries in the vector to hll, ull, and ell. clears vector
        let flush = |buf: &mut Vec<u64>| {
            let mut h = hll.lock().unwrap();
            let mut u = ull.lock().unwrap();
            let mut e = ell.lock().unwrap();
            for &x in buf.iter() {
                h.push_hash64(x);
                u.add(x);
                e.add_hash(x);
            }
            buf.clear();
        };

        // sketch hll, ell, and ull into one big sketch
        files
            .par_iter()
            .for_each(|file_name| {
                let mut reader = parse_fastx_file(file_name).expect("Invalid input file");
                let mut buf: Vec<u64> = Vec::with_capacity(CAPACITY);

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

                            buf.push(xxh3_64_with_seed(&masked.to_le_bytes(), seed as u64));

                            if buf.len() >= CAPACITY {
                                flush(&mut buf);
                            }
                        }
                        flush(&mut buf);
                    }
                    
                }
            });
        
        let hll_card = hll.lock().unwrap().len();
        let ell_card = ell.lock().unwrap().estimate();
        let ull_card = MaximumLikelihoodEstimator.estimate(&ull.lock().unwrap());

        // relative squared error
        fn squared_error(est: f64, truth: f64) -> f64 {
            ((est - truth) / truth).powi(2)
        }
        
        let hll_err = squared_error(hll_card, cardinality);
        let ull_err = squared_error(ull_card, cardinality);
        let ell_err = squared_error(ell_card, cardinality);
        
        writeln!(writer, "hll cardinality: {}, error: {}", hll_card, hll_err).expect("error writing hll result");
        writeln!(writer, "ull cardinality: {}, error: {}", ull_card, ull_err).expect("error writing ull result");
        writeln!(writer, "ell cardinality: {}, error: {}", ell_card, ell_err).expect("error writing ull result");

        hll_errors.push(hll_err);
        ull_errors.push(ull_err);
        ell_errors.push(ell_err);
        
        i = i + 1;

        // serialize sketches. they are not compressed here. 
        if i == num_seeds { // if we are on the last seed, save it so we know the file size
            let mut writer = BufWriter::new(File::create("hll_sketch.bin")?);
            hll.lock().unwrap().save(&mut writer)?;
            writer.flush()?;

            let mut writer = BufWriter::new(File::create("ull_sketch.bin")?);
            ull.lock().unwrap().save(&mut writer)?;
            writer.flush()?;

            let mut writer = BufWriter::new(File::create("ell_sketch.bin")?);
            ell.lock().unwrap().save(&mut writer)?;
            writer.flush()?;
        }
    }
    

    fn variance(errors: Vec<f64>) -> f64 {
        let sum = errors.iter().sum::<f64>() / errors.len() as f64;
        sum
    }

    writeln!(writer, "**************** RMSE results ****************").expect("error writing to log file");
    writeln!(writer, "hll variance: {}", variance(hll_errors)).expect("error writing hll rmse");
    writeln!(writer, "ull variance: {}", variance(ull_errors)).expect("error writing ull rmse");
    writeln!(writer, "ell variance: {}", variance(ell_errors)).expect("error writing ell rmse");

    Ok(())
}