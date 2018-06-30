# Structure of the code

* mirtop/bam
  * __bam.py__ 
    * `read_bam`: reads BAM files with pysamtools and store in a key - value object
  * __filter.py__
    * `tune`: if option `--clean` is on, filter according generic rules
    * `clean_hits`: get the top hits
* mirtop/gff
  * __init.py__ wraps the convertion process to GFF3
  * __body.py__ `create` will create the line according GFF format established.
    * `read_gff_line`: Inside a for loop to read line of the file. It'll return and structure key:value dictionary for each column.
  * __header.py__ generate header and read header section.
  * __check.py__ checks header and single lines to be valid according GFF format  (NOT IMPLEMENTED)
  * __stats.py__ GFF stats counting number of isomiR, their total and average expression
  * __query.py__ accept SQlite queries after option -q ""
  * __convert.py__
    * `create_counts` table of counts
    * allow filtering by attribute
    * allow collapse by miRNA/isomiR type
  * __filter.py__, parse from query (NOT IMPLEMENTED)
* mirtop/mirna
  * __fasta.py__: 
    * `read_precursor` fasta file: key - value
  * __realign.py__:
    * `hits`: class that defines hits
    * `isomir`: class that defines each sequence
    * `cigar_correction`: function that use CIGAR to make sequence to miRNA alignemt
    * `read_id` and `make_id`: shorter ID for sequences
    * `make_cigar`: giving an alignment return the CIGAR of it
    * `reverse_complement`: return the reverse complement of a sequence
    * `align`: uses biopython to align two sequences of the same size
    * `expand_cigar`: from a 12M to MMMMMMMMMMMM
    * `cigar2snp`: from CIGAR code to list of changes with position and reference and target nts
  * __mapper.py__: 
    * `read_gtf` file: map genomic miRNA position to precursos position, then it needs genomic position for the miRNA and the precursor. Return would be like {mirna: [start, end]}
  * __annotate.py__:
    * `annotate`: read isomiRs and populate all attributes related to isomiRs
 * mirtop/importer:
    * seqbuster.py
    * prost.py
    * srnabench.py
    * isomirsea.py
 * mirtop/exporter:
    * isomirs.py: export file to match [isomiRs BioC package](https://github.com/lpantano/isomiRs).
 * data/examples/
   * check gff files: example of correct, invalid, warning GFF files
   * check BAM file
   * check mapping from genome position to precursor position, example of +/- strand. Using `mirtop/mirna/map.read_gtf`.
   * check clean option: sequence mapping to multiple precursors/mirna, get the best score. Using `mirtop/bam/filter.clean_hits`.

To add new sub-commands, modify the following:

* mirtop/lib/parse.py
  * query: TODO
  * transform: TODO
  * create: TODO
  * check: TODO 
