[![Latest Version](https://img.shields.io/crates/v/hypermash?style=for-the-badge&color=mediumpurple&logo=rust)](https://crates.io/crates/hypermash)

# Fast and Memory Efficient Genome/Metagenome Sketching via HyperMinHash

Genome sketching can be extremely accurate but requires a huge amount of memory for MinHash-like algorithms. Recently, a new algorithm combining MinHash and HyperLogLog, called HyerMinHash was invented (1), which can perform MinHash in loglog space, a significant decrease in space/memory requirement. Together with [lukaslueg](https://github.com/lukaslueg), we first create a Rust library [hyperminhash](https://github.com/lukaslueg/hyperminhash) and then combine rolling hashing with HyperMinHash for extremely fast processing of genomic sequences. Xxhash3 was used as the underlying hashing technique. 

More recently, an algorithm named Ultraloglog was invented (2). It is similar to Hyperloglog but with 28% more space efficiency due to a faster estimator. Ultraloglog also has better compaction when using compressing algorithms. Ultraloglog was implemented with [waynexia] (https://github.com/waynexia). Both HyperMinHash and Ultraloglog are options available for use on our tool. 

We employed a simple producer-consumer model to also reduce memory requirement for large files, e.g., metagenomic files. Both sketching and distance computation are parallelized to make full use of all CPU threads/cores. 

There are two main subcommands, sketch and dist. Sketch is the sketching command and outputs 3 files; one file containing the sketches of the genomes, one file containing the genome files used, and one file containing parameters used for the command. Dist is the command that "reads" the sketch files and outputs a file containing the distances between the query and reference genomes, which is specified by the user. More details on these commands are under "Usage". 

We hope that you find this tool helpful in your scientific endeavors!

## Quick install
```bash
### pre-compiled binary for Linux
wget https://github.com/jianshu93/hypermash/releases/download/v0.1.0/hypermash_Linux_x86-64_v0.1.0.zip
unzip hypermash_Linux_x86-64_v0.1.0.zip
chomd a+x ./hypermash
./hypermash -h

### Install from cargo, install cargo first here: https://rustup.rs, cargo will be installed by default
cargo install hypermash

### compiling from source
git clone https://github.com/jianshu93/hypermash
cd hypermash
cargo build --release
./target/release/hypermash -h

```

## Usage
```bash

 ************** initializing logger *****************

Fast and Memory Efficient Genome Sketching via HyperMinhash

Subcommand 1: Sketching
Sketching: hypermash sketch --file <file> --output <output_prefix> --kmer <kmer_length> --threads <num_threads>-algorighm <algorithm> -precision <precision_ull>

Options:
  -f, --file <file>                 File containing list of FASTA files
  -o, --output <output_prefix>      Prefix you would like your output file names to start with
  -k, --kmer <kmer_length>          Length of k-mers
  -t, --threads <num_threads>       Number of threads you would like to use. Default to the number of cores on your device
  -a, --algorithm <algorithm>       Algorithm of choice. Either hmh for hyperminhash, or ull for ultraloglog  
  -p, --precision <precision_ull>   Precision to use, only for Ultraloglog. Default to 10. 
  -v, --version                     Print version


Subcommand 2: Distance
Distance: hypermash dist --query <query__prefix> --reference <ref_prefix> --output <output_prefix>--threads <num_threads> --estimator <estimator_ull>
Options:
  -q, --query <query_prefix>        Prefix to search for your query genome files. Should match what you put as "output" from sketch. 
  -r, --reference <ref_prefix>      Prefix to search for your reference genome files. Should match what you put as "output" from sketch. 
  -o, --output <output_prefix>      Prefix you would like your output file names to start with
  -t, --threads <num_threads>       Number of threads you would like to use. Default to the number of cores on your device
  -e, --estimator <estimator>       Estimator to use, only for Ultraloglog sketches. Either "fgra" for Fast Graph-based Rank Aggregation or "ml" for maximum likelihood estimator, default to "ml".  
  -v, --version                     Print version


```



```bash
ls ./data/*.fasta > query_list_strep.txt
ls ./data/*.fasta > ref_list_strep.txt
hypermash --query_file ./query_list_strep.txt -r ref_list_strep.txt -k 16 -o dist.txt
```

## Output

Output format is the same with Mash, first column query, second column reference name， third column Mash distance

## References
1. Yu YW, Weber GM. Hyperminhash: Minhash in loglog space. IEEE Transactions on Knowledge and Data Engineering. 2020 Mar 17;34(1):328-39.
2. Ertl O. UltraLogLog: A Practical and More Space-Efficient Alternative to HyperLogLog for Approximate Distinct Counting. Proceedings of the VLDB Endowment. 2024 March 1;17(7):1655-1668. 
