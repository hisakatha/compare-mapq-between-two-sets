This program takes an alignment BAM file and two sets of chromosome (reference) names,
and outputs the highest MapQ scores in multiple alignment hits of one read for each set of chromosomes.
An alignment BAM file MUST be sorted by name beforehand.

An output for a read will be as follows:

`read_name,top_mapq_in_set1,count_in_set1,top_mapq_in_set2,count_in_set2,abs_diff_in_mapq`

For example, given the input alignment has scores as follows:

```
read_name, reference, MapQ
read1, chr_a1, 10
read1, chr_a2, 30
read1, chr_b1, 20
read1, chr_b2, 100
read2, chr_a2, 0
read2, chr_b1, 0
read2, chr_b2, 200
```
and chromosome sets are difined as `chromosome set1 = {chr_a1, chr_a2}`, and `chromosome set2 = {chr_b1, chr_b2}`,
then the output will be as follows:

```
read1,30,2,100,2,70
read2,0,1,200,2,200
```

# How to install
Please set up a Rust environment, and then

- `cargo build --release` and find the executable at `target/release`
- or `cargo install --path .`
