# Fast and Memory Efficient Genome/Metagenome Sketching via HyperMinhash

Genome sketching can be extremely accurate but requires a huge amount of memory for MinHash-like algorithms. Recently, a new algorithm combining MinHash and HyperLogLog, called HyerMinHash was invented (1), which can perform MinHash in loglog space, a significant decrease in space/memory requirement. We implemented this algorihm and apply it for genome sketching. We first create a library with [lukaslueg](https://github.com/lukaslueg), and then combining rolling hashing with HyperMinHash for extremely fast processing of genomic sequences. 

We use simple producer-consumer model to also reduce memory requirement for large files, e.g., metagenomic files. Both sketching and distance computation are parallelized to make use of all CPU threads/cores. 

## Quick install
```bash
git clone https://github.com/jianshu93/hypermash
cd hypermash
cargo build --release
./target/release/hypermash -h
```

## Usage
```bash
hypermash --query_file ./query_list_strep.txt -r ref_list_strep.txt -k 16 -o dist.txt
```


## References
1. Yu YW, Weber GM. Hyperminhash: Minhash in loglog space. IEEE Transactions on Knowledge and Data Engineering. 2020 Mar 17;34(1):328-39.
