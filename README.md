# Fast and Memory Efficient Genome/Metagenome Sketching via HyperMinHash

Genome sketching can be extremely accurate but requires a huge amount of memory for MinHash-like algorithms. Recently, a new algorithm combining MinHash and HyperLogLog, called HyerMinHash was invented (1), which can perform MinHash in loglog space, a significant decrease in space/memory requirement. We implemented this algorihm and apply it for genome sketching. Together with [lukaslueg](https://github.com/lukaslueg), we first creata a Rust library [hyperminhash](https://github.com/lukaslueg/hyperminhash) and then combine rolling hashing with HyperMinHash for extremely fast processing of genomic sequences. Xxhash3 was used as the underlying hashing technique. 

We employed a simple producer-consumer model to also reduce memory requirement for large files, e.g., metagenomic files. Both sketching and distance computation are parallelized to make full use of all CPU threads/cores. 

## Quick install
```bash
### pre-compiled binary for Linux
wget https://github.com/jianshu93/hypermash/releases/download/v0.1.0/hypermash_Linux_x86-64_v0.1.0.zip
unzip hypermash_Linux_x86-64_v0.1.0.zip
chomd a+x ./hypermash
./hypermash -h


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

Usage: hypermash --query_file <query_files> --ref_file <reference_files> --kmer <kmer_length> --output <output_file>

Options:
  -q, --query_file <query_files>    File containing list of query FASTA files
  -r, --ref_file <reference_files>  File containing list of reference FASTA files
  -k, --kmer <kmer_length>          Length of k-mers
  -o, --output <output_file>        Output file to write results
  -h, --help                        Print help
  -V, --version                     Print version

```


```bash
ls ./data/*.fasta > query_list_strep.txt
ls ./data/*.fasta > ref_list_strep.txt
hypermash --query_file ./query_list_strep.txt -r ref_list_strep.txt -k 16 -o dist.txt
```


## References
1. Yu YW, Weber GM. Hyperminhash: Minhash in loglog space. IEEE Transactions on Knowledge and Data Engineering. 2020 Mar 17;34(1):328-39.
