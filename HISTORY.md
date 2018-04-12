- 0.3.*a

 * Fix bug for duplicated isomiRs tags. Thanks @AlisR.
 * Fix bug in order of merged gff file. Thanks @AlisR.
 * Add module to read GFF/GTF line in body.py
 * Add version line to stats output
 * Improve PROST! importer
 * Fix output for isomiRs package

- 0.2.*

 * Make GTF default output
 * Add function to get SNPs from Variant attribute
 * Improve PROST with last version output
 * Add isomiR-SEA compatibility
 * Fix sRNAbench exact match to NA in GFF
 * Change stats to use only 1 level isomiR classification
 * Add GFF to count matrix
 * Add read_attributes function
 * Improve isomiR reading from srnabench tool
 * Add PROST to supported tools

- 0.1.7
 
 * Remove deletion from addition isomiRs
 * Support for srnabench output
 * Fix bug mixing up source column
 * Support Seqbuster output
 * Functin to guess database used from GTF file through --mirna parameter
 * Adapt output format to https://github.com/miRTop/incubator/blob/master/format/definition.md

- 0.1.5
 
 * add function to check correct annotation
 * add test data for SAM parsing
 * add script to simulate isomiRs
 * parse indels from bam file
 
- 0.1.4

 * fix index BAM file command line
 * add function to accept indels and test unit
 * change header from subs -> mism to be compatible with isomiRs
