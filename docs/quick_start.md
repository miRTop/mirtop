# Quick Start

## From Bam files to GFF3

You can use the example data. Here the reads have been mapped to the precursor sequences.

```
git clone mirtop
cd mirtop/data/annotate
mirtop gff -sps hsa --hairpin hairpin.fa --gtf hsa.gff3 -o test_out sim_isomir.bam
```
## From `seqbuster::miraligner` files to GFF3

miRNA annotation generated from [miraligner](https://github.com/lpantano/seqbuster) tool:

```
cd mirtop/data
mirtop gff --format seqbuster --sps hsa --hairpin annotate/hairpin.fa --gtf annotate/hsa.gff3 -o test_out seqbuster/reads.mirna
```

## From `sRNAbench` files to GFF3


miRNA annotation generated from [sRNAbench](http://bioinfo2.ugr.es:8080/ceUGR/srnabench/) tool:

```
cd mirtop/data
mirtop gff -sps hsa --hairpin annotate/hairpin.fa --gtf annotate/hsa.gff3 -o test_out srnabench
```

## Get statistics from GFF

Get number of isomiRs and miRNAs annotated in the GFF file by isomiR category.

```
cd mirtop/data
mirtop stats -o test_out example/gff/correct_file.gff
```

## Compare GFF file with reference

Compare the sequences from two or more GFF files. The first one will be used as the reference data.

```
cd mirtop/data
mirtop stats -o test_out example/gff/correct_file.gff example/gff/alternative.gff
```

