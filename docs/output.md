# Output

## GFF command

The `mirtop gff` generates the GFF3 adapter format to capture miRNA variations. The output is explained [here](https://github.com/miRTop/incubator/blob/master/format/definition.md).

## Stats command

The `mirtop stats` generates a table with different statistics for each type of isomiRs:

* total counts
* average counts
* total sequences

It generates as well a JSON file with the same information to be integrated easily with QC tools like [MultiQC](https://multiqc.info/).

## Compare command

The `mirtop compare` generates a tabular file with information about the difference and similarities. The first file in the command line will be considered the reference and the following files will be compared to the reference. Each line of the output has the following information for each file:

* sample
* idu 
* seq 
* tag: `E` if not in reference, `D` detected in both, `M` missing in target file 
* same_mirna: if the sequence map to the same miRNA in the reference and target file
* one column for each isomiR type with the following tags: `FP` (variation not in reference), `TP` (variation in both), `FN` (variation not in target file) 

## Counts command

The `mirtop counts` generates a tabular file with the following columns:

* unique identifier
* read sequence
* miRNA name
* Variant attribute from GFF3 column
* One column for each isomiR type showing the exact variation 
* One column for each sample with the counts for that sequence

## Export command

The `mirtop export` generates different files from a mirGFF3 file:

* [isomiRs](https://bioconductor.org/packages/release/bioc/html/isomiRs.html) compatible files
* [FASTA files](https://en.wikipedia.org/wiki/FASTA_format)
* [VCF files](https://samtools.github.io/hts-specs/VCFv4.2.pdf)