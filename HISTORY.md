- 0.4.19a

- 0.4.18

 * Cast map object to list to avoid errors in py3.
 * Support Manatee output.
 * Support chunk reading for genomic BAM files.
 * Support chunk reading for seqbuster files.
 * Support chunk reading for BAM files.
 * Normalize functions to support different databases.
 * Support miRgeneDB.
 * Export to VCF. Thanks to Roderic Espin.
 * Support isomiRs that go beyond 5p end
 * Support genomic coordinates.
 * Fix missing reads when using --keep-read in the final mirtop.gff file.
 * Allow longer truncation and addition events.
 * Accept seqbuster input without frequency column.
 * Allow keep name of the sequence.
 * Accept indels in snv category.
 * Additions are only last nucleotides that are mismatches.
 * Adapt mintplate license.
 * Revert sign in iso_5p, replace snp by snv.
 * Skip lines that contain malformed UID.
 * Add FASTA as an exporter from GFF.
 * Fix BAM parsing to new GFF rules.
 * Add the possibility to work with spikeins to detect random variability.
 * Fixing UID attribute for tools that don't use our cypher system
 * Add class to parse GFF line as a first move toward isolation
 * Add JSON log for stats command.

- 0.3.17
 * Normalize the read of the tool outputs.
 * Add docs with autodoc plugin.
 * Validator by @Vbarrera.
 * Improve examples commands and test coverage.
 * Only counts sequences with Filter == Pass during stats.
 * Counts cmd add nucleotide information when --add-extra option is on.
 * Fix error in stats that open the file in addition mode.
 * Importer for sRNAbench just convert lines from input to GFF format.
 * Skip lines with non-valid UID or miRNAs not in reference at counts cmd.
 * Fix separators in counts cmd.
 * Make --sps optional.
 * Add synthetic data with known isomiRs to data set.
 * Allow extra columns when converting to counts TSV file.
 * Allow extra attributes for isomir-sea as well.
 * Allow extra attributes to show the nts
   that change in each isomiR type.
 * Fix Expression attrb when join gff files.Thanks @AlisR.
 * Print help when no files are giving to any subcommand.
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
