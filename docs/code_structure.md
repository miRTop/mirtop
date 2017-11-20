# Structure of the code

* mirtop/bam
  * __bam.py__ 
    * `read_bam`: read BAM file with pysamtools and store in a key - value object
  * __filter.py__
    * `tune`: if option --clean is on, filter according generic rules
    * `clean_hits`: get the top hits
* mirtop/gff
  * __init.py__ wraps the convertion process to GFF3
  * __body.py__ `create` will create the line according GFF format established.
  * __header.py__ generate header
  * __check.py__ check header and single lines to be valid according GFF format   
  * __query.py__ accept SQlite queries after option -q ""
  * __convert.py__
    * `create_counts` table of counts
    * allow filtering by attribute
    * allow collapse by miRNA/isomiR type
  * __filter.py__, parse from query
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
  * __mapper.py__: 
    * `read_gtf` file: map genomic miRNA position to precursos position, then it needs genomic position for the miRNA and the precursor. Return would be like {mirna: [start, end]}
  * __annotate.py__:
    * `annotate`: read isomiRs and populate all attributes related to isomiRs
 * mirtop/importer:
    * seqbuster.py
    * prost.py
    * srnabench.py
    * isomirsea.py
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
