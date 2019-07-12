# Quick Start

## Importer

### From Bam files to GFF3

```
git clone mirtop
cd mirtop/data
```

You can use the example data. Here the reads have been mapped to the precursor sequences.

```
mirtop gff -sps hsa --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 -o test_out sim_isomir.bam
```

### From `seqbuster::miraligner` files to GFF3

miRNA annotation generated from [miraligner](https://github.com/lpantano/seqbuster) tool:

```
mirtop gff --format seqbuster --sps hsa --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 -o test_out examples/seqbuster/reads.mirna
```

### From `sRNAbench` files to GFF3

miRNA annotation generated from [sRNAbench](http://bioinfo2.ugr.es:8080/ceUGR/srnabench/) tool:

```
mirtop gff --format sranbench -sps hsa --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 -o test_out srnabench examples/srnabench
```

### From `PROST!` files to GFF3

miRNA annotation generated from [PROST!]() tool. Export isomiRs tab from excel file to a tabular text format file.

```
mirtop gff --format prost -sps hsa --hairpin  examples/annotate/hairpin.fa --gtf  examples/annotate/hsa.gff3 -o test_out examples/prost/prost.example.txt
```

### From `isomiR-SEA` files to GFF3

miRNA annotation generated from [isomiR-SEA]() tool.

```
mirtop validate examples/gff/correct_file.gff
```

## Operations

### Validator

To validate your mirGFF3 file and make sure if follows the current format:

```
mirtop gff --format isomirsea -sps hsa --hairpin  examples/annotate/hairpin.fa --gtf  examples/annotate/hsa.gff3 -o  test_out examples/isomir-sea/tagMir-all.gff
```

### Get statistics from GFF

Get number of isomiRs and miRNAs annotated in the GFF file by isomiR category.

```
cd mirtop/data
mirtop stats -o test_out example/gff/correct_file.gff
```

### Compare GFF file with reference

Compare the sequences from two or more GFF files. The first one will be used as the reference data.

```
cd mirtop/data
mirtop compare -o test_out example/gff/correct_file.gff example/gff/alternative.gff
```

### Updates mirGFF3

Updates older versions with the most current one.

```
cd mirtop/data
mirtop update -o test_out_mirs examples/versions/version1.0.gff
```

## Export

### Export file to isomiRs format

To be compatible with [isomiRs](https://bioconductor.org/packages/release/bioc/html/isomiRs.html) bioconductor package use:

```
cd mirtop/data
mirtop export -o test_out_mirs --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 examples/gff correct_file.gff                                   
```

### Export file to FASTA format

```
cd mirtop/data
mirtop export -o test_out_mirs --format fasta -d -vd --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 examples/gff/correct_file.gff
```

### Export file to VCF format

```
cd mirtop/data
mirtop export -o test_out_mirs --format vcf --hairpin examples/annotate/hairpin.fa --gtf examples/a
nnotate/hsa.gff3 examples/gff/correct_file.gff
```

### Get count file

This file it is useful to load into R as a matrix. It contains the minimal information about each sequence and the count data in columns for each samples.

```
cd mirtop/data
mirtop counts -o test_out_mirs --hairpin examples/annotate/hairpin.fa --gtf examples/annotate/hsa.gff3 examples/synthetic/let7a-5p.gtf                              
```
