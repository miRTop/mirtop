* mirtop/bam
 * bam.py: 
   * `read_bam`: read BAM file with pysamtools and store in a key - value object
 * filter.py
   * `cleaner`: if option --clean is on, filter according generic rules
   * `annotate`: get isomiRs and CIGAR option
   * `clean_hits`: get the top socred hits
* mirtop/gff
 * body.py create entris with gff python package
 * header.py generate header
 * check.py check header and single lines to be valid according GFF format
 * query.py accept SQlite queries after option -q ""
 * convert.py
  * `create_counts` table of counts
  * allow filtering by attribute
  * allow collapse by miRNA/isomiR type
 * filter.py, parse from query
* mirtop/mirna: the input is. fasta file and a GFF file with position that needs to be mapped to precursors coordinates.
 * fasta.py: 
   * `read_precursor` fasta file: key - value
 * realign.py:
   * `hits`: class that defines hits
   * `isomir`: class that defines each sequence
   * `cigar_correction`: function that use CIGAR to make sequence to miRNA alignemt
 * map.py: 
   * `read_gtf` file: map genomic miRNA position to precursos position, then it needs genomic position for the miRNA and the precursor. Return would be like {mirna: [start, end]}
* data/examples/: example data would be in `data/examples`
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