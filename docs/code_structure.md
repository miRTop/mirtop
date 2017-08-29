* mirtop/bam
 * bam.py read BAM file with pysamtools and store in a key - value object
 * filterBam.py, if option --clean is on, filter according generic rules
* mirtop/gff
 * gff.py, create entris with gff python package
 * header.py, generate header
 * check.py, check header and single lines to be valid according GFF format
 * query.py accept SQlite queries after option -q ""
 * transform.py
  * gff -> table of counts
  * allow filtering by attribute
  * allow collapse by miRNA/isomiR type
 * filterGFF.py, parse from query
* mirtop/mirna: the input is. fasta file and a GFF file with position that needs to be mapped to precursors coordinates.
 * fasta.py: read precursor fasta file: key - value
 * mirna.py: read mature file: key - value: key is precursor name, value is another dict: {mirna: [start, end]}
 * map.py
* tests/: example data would be in `data/examples`
 * check gff files: example of correct, invalid, warning GFF files
 * check BAM -> FILE
 * check mapping from genome position to precursor position, example of +/- strand
 * check clean option: sequence mapping to multiple precursors/mirna, get the best score

Add comands:
* mirtop/lib/parse.py
 * query
 * transform
 * create
 * check