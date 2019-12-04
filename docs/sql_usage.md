# SQLite Command Options

The `mirtop sql` package allows users to create and query a SQLite database using GFF3 file format as input. This option creates a database with two tables such as summary and data\_sets. The summary table contains the header portion of the gff file and data\_sets contains the body of the gff containing rest of the information related to miRNAs. 

Use `mirtop sql -h` to display help information

```
mirtop sql -h
```
SQL arguments:
```
['sql', '-h']
usage: mirtop sql [-h] [--db] (-c | -q) [--gff] [-o] [-t] [-txto] [-col] [-n]
                  [-miR] [-var] [-f] [-l] [-e] [-d] [-vd]

optional arguments:
  -h, --help          show this help message and exit
  --db                SQL Database name. (default: mirtop.db)
  -c, --create        Creates a SQLite database from GFF
  -q, --query         Query from a SQLite database
  -d, --debug         max verbosity mode
  -vd, --print_debug  print debug messages on terminal

SQL create usage mode:
  --gff               GFF file with precursor and mature position to genome
  -o , --out          Directory of output files

SQL query usage mode:
  -t , --table        Specify table name to use
  -txto , --txtout    Writes the output of the query to a file speficied. Format (-fmt) is a tab-delimited text file by default
  -col , --columns    Select specific columns from the table to display (Default: all columns), or use with -n option to return n-counts. For information of the available columns see 'show-schema' or 'show-columns'. NOTE: option -e select must be applied!
  -n , --count        Returns 'n' counts for the query. Options 'T' for True, if not 'F' (Default: -n F). NOTE: option -e select must be applied and accepts only one column from -col option.
  -miR , --miRNA      Specify the miRNA names to query. For multiple miRNAs use comma(,) as separator; or supply a text file (.txt) separated with a new line character.
  -pm , --miRNA_prefix (3 digit species code, ex. hsa)
                      Specify the prefix name for miRNAs to query. Example: -pm hsa -miR let-7a-5p results into querying hsa-let-7a-5p.
  -var , --variant    Specify one or more types of variants to query. Use comma(,) as separator
                          The following choices are supported:
                              iso_5p                  - indicates the shift at the reference 5' miRNA
                              iso_3p                  - indicates the shift at the reference 3' miRNA
                              iso_add3p               - Number of non-template nucleotides added at 3p
                              iso_add5p               - Number of non-template nucleotides added at 5p
                              iso_snv_seed            - when affected nucleotides are between [2-7]
                              iso_snv_central_offset  - when affected nucleotide is at position [8]
                              iso_snv_central         - when affected nucleotides are between [9-12]
                              iso_snv_central_supp    - when affected nucleotides are between [13-17]
                              iso_snv                 - anything else

  -f , --filter       Specify Filter tag attribute. Options: Pass, Reject. (Default: None)
  -l , --limit        Specify the number of rows to output. (Example: --limit 30, to limit the first 30 rows)
  -e , --expr         Expression is the query that you want to run; (-e "<statement>")
                          Choices supports the following:
                             show-tables              - Displays tables in the database (default: mirtop.db)
                             show-schema              - Displays the table schema (requires -t)
                             show-columns             - Displays available columns in the table
                             describe-gff             - Prints out the header information from the GFF file
                             isomirs-per-mirna        - Displays the count of isomiRs for miRNA (requires -miR)
                             select                   - Allows specific query construction.
                                                        Example: mirtop sql --db tmp_mirtop/SRR333680_revised2.db -qe select -var iso_5p,iso_3p -miR hsa-let-7a-5p,hsa-let-7d-5p -l 30
                                                        The above expression evaluates to selecting miRNAs in -miR with variants in -var and prints out the first 30 rows in --limit.


```
## Creating Database 

The `mirtop sql -c` takes `--gff` gff3 file and creates a database with name `--db` 

### From GFF3 to SQLite database

```
git clone mirtop
cd mirtop/data
```

You can use the example data to create a database file from GFF3. 

```
mirtop sql -c --gff examples/annotate/SQL_sample.gff -o examples/annotate/ --db SQL_sample.db
```
NOTE-1: If you are re-creating the database with the same name, make sure to delete the existing database in the working directory. As the database from the same set of samples will append the gff rows to the existing database and different set of samples throws an error as shown below.

`sqlite3.OperationalError: table data_sets has xy columns but yz values were supplied`
Where, xy is the number of columns present in database and yz are the number of column-values sent to append to the database.


## Querying a Database

Understanding the database structure helps to quickly query and fetch the results.  The following content lists the useful commands to query based on the availabe options. To start with the available options type: `mirtop sql -h`. 

NOTE-2: To query from a database, the argument `mirtop sql -q --db <db_name.db>` must always be specified with the correct name of the database along with absolute path (if required). The query will begin with an expression (-e) that is specified by the User. Passing a suitable query is mandatory such as: `mirtop sql -q --db <db_name.db> -e <expression>`. The following expressions allows users the ability to query from the database:

* show-tables              - Displays tables in the database (default: mirtop.db)
* show-schema              - Displays the table schema (requires -t)
* show-columns             - Displays available columns in the table
* describe-gff             - Prints out the header information from the GFF file
* isomirs-per-mirna        - Displays the count of isomiRs for miRNA (requires -miR)
* select                   - Allows specific query construction 

### show-tables:
Display the contents of the database. The data or information from a gff file is loaded as a table, and the collection of tables make it a database. To see the contents of the database, we should see the available tables in the database.  

```
cd mirtop/data
```

You can use the example data (query\_sample.db) to query the contents of the database. 

```
mirtop sql -q --db examples/annotate/query_sample.db -e show-tables
```

OUTPUT: 

```
11/29/2019 01:03:49 INFO Run Convert GFF.
 +------------------------- +
 | Tables                   |
 +------------------------- +
 | data_sets                |
 | summary                  |
 +------------------------- +
11/29/2019 01:03:49 INFO It took 0.000 minutes 
```


### show-schema: 
Display the schema of a table. The table schema represents the data type of each column and lists available columns such as samples names. Generally this information is essential for developers and can aid in debugging any errors during the query operation. show-schema essentially requires one to specify a table name for which the schema is intended to be displayed. 

```
mirtop sql -q --db examples/annotate/query_sample.db -e show-schema -t summary
```

OUTPUT:
```
11/29/2019 01:13:55 INFO Run Convert GFF.
 +---------------------------------------------------------+
 | Sl | Field                         | Type | NULL | Key  |
 +---------------------------------------------------------+
 |  0 | version                       | text |  NO  |      |
 |  1 | source                        | text |  NO  |      |
 |  2 | data_sets                     | text |  NO  |      |
 |  3 | tools                         | text |  NO  |      |
 |  4 | commands_exec                 | text |  NO  |      |
 |  5 | filter_tags                   | text |  NO  |      |
 |  6 | citation                      | text |  NO  |      |
 |  7 | records                       | real |  NO  |      |
 |  8 | date_stamp                    | text |  NO  |      |
 +---------------------------------------------------------+
11/29/2019 01:13:55 INFO It took 0.000 minutes
```

### show-columns: Query the columns of the table

This is similar to show-schema, however, only columns from the table are listed excluding the data type and other information. This is required for select query options, explained later. 
NOTE-3: This parameter is only available for table data\_sets. 

```
mirtop sql -q --db examples/annotate/query_sample.db -e show-columns 
```
OUTPUT:
```
11/29/2019 01:18:55 INFO Run Convert GFF.

Serial  Column names
  1     seqID
  2     source_file
  3     type
  4     start
  5     end
  6     score
  7     strand
  8     phase
  9     UID
  10    Read
  11    Name
  12    Parent
  13    Variant
  14    iso_5p
  15    iso_3p
  16    iso_add3p
  17    iso_add5p
  18    iso_snv
  19    iso_snv_seed
  20    iso_snv_central
  21    iso_snv_central_offset
  22    iso_snv_central_supp
  23    source
  24    cigar
  25    hits
  26    alias
  27    genomic_pos
  28    filter
  29    seed_fam
  30    SRR333680_1

11/29/2019 01:18:55 INFO It took 0.000 minutes
```

### describe-gff: 
This option prints out the header information of the gff3 (table: summary) 

```
mirtop sql -q --db examples/annotate/query_sample.db -e describe-gff 
```
OUTPUT:
```
11/29/2019 01:21:18 INFO Run Convert GFF.

Serial  Column names    Description
  1     version "--"
  2     source  "miRBase22"
  3     data_sets       "SRR333680_1"
  4     tools   "--"
  5     commands_exec   "--"
  6     filter_tags     "--"
  7     citation        "--"
  8     records "11672.0"
  9     date_stamp      "November 22, 2019 16:41:50"

11/29/2019 01:21:18 INFO It took 0.000 minutes
```

### isomirs-per-mirna: 
This expression provides a summary of the number of the isomiRs for the query miRNA (`-miR`)

Query miRNA (`-miR` or `--miRNA`) can be a particular miRNA of interest, or list of a few miRNAs (separated by comma) or a file of miRNAs (.txt file). The output can be redirected to a text document with `-txto <text_file.txt>` argument. The users can also chooses the prefix for the miRNAs with `-pm`. For example: A query `-pm hsa -miR let-7a-5p,let-7d-5p` will result into querying the database as `-miR hsa-let-7a-5p,hsa-let-7d-5p`.   

```
mirtop sql -q --db examples/annotate/query_sample.db -e isomirs-per-mirna -miR hsa-let-7a-5p
			(-- OR --)
mirtop sql -q --db examples/annotate/query_sample.db -e isomirs-per-mirna -miR hsa-let-7a-5p,hsa-let-7d-5p 
			(-- OR --)
mirtop sql -q --db examples/annotate/query_sample.db -e isomirs-per-mirna -miR let-7a-5p,let-7d-5p -pm hsa 
			(-- OR --)
mirtop sql -q --db examples/annotate/query_sample.db -e isomirs-per-mirna -miR examples/annotate/miRNA_sample_list.txt -txto examples/annotate/queryOutput_isomirs.txt
```
OUTPUT:

```
11/29/2019 02:26:44 INFO Run Convert GFF.

OUTPUT:
1. isomiRs for miRNA hsa-let-7a-5p: 397
2. isomiRs for miRNA hsa-let-7d-5p: 59

11/29/2019 02:26:44 INFO It took 0.000 minutes
```

### select: 
This expression represents the SELECT statement used in MySQL database. It offers a lot of query options and can be combined with one or more optional arguments mentioned below. 

* limit (`-l or --limit`): Specify the number of rows to output

The following command uses the conventional "SELECT * FROM data\_sets" option to display the contents of the table data\_sets from database query\_sample.db. However using `-l 2 or --limit 2` limits the output to display only the first two rows. If -l is not provided, it prints out all of the rows on the terminal (or) prints to a file if `-txto fileName.txt` is provided. 
``` 
mirtop sql -q --db examples/annotate/query_sample.db -e select --limit 2
```
OUTPUT:
```
11/29/2019 03:04:28 INFO Run Convert GFF.
seqID   source_file     type    start   end     score   strand  phase   UID     Read    Name    Parent  Variant iso_5p  iso_3p  iso_add3p       iso_add5p       iso_snv iso_snv_seed    iso_snv_central      iso_snv_central_offset  iso_snv_central_supp    source  cigar   hits    alias   genomic_pos     filter  seed_fam        SRR333680_1
hsa-miR-342-3p  miRBase22       isomiR  61.0    85.0    .       +       .       0t@TeV#5v2      TCTCACAAAGAAATCGCACCCGTTC       hsa-miR-342-3p  hsa-mir-342     iso_snp,iso_add:+2      None    None None    None    1.0     0.0     0.0     0.0     0.0     miRBase22       7MA15MTC        None    None    None    Pass    None    1
hsa-let-7i-5p   miRBase22       isomiR  6.0     28.0    .       +       .       7AwwRIB71       TGAGGTAGTAGTTTGTGCTGTTG hsa-let-7i-5p   hsa-let-7i      iso_3p:+1       None    1.0     None    None 0.0     0.0     0.0     0.0     0.0     miRBase22       23M     None    None    None    Pass    None    1
11/29/2019 03:04:28 INFO It took 0.000 minutes
```

NOTE-4: In the above query option, we could limit the number of rows to be printed. What about columns? Can we limit them as well? and the answer is Yes. 

* columns: Select specific columns from the table to display

The following command use the conventional "SELECT seqID,UID,Read,iso\_5p,iso\_3p,start,end FROM data\_sets LIMIT 2" option to display the contents of the table data\_sets from database query\_sample.db for specific columns. It prints out all the rows on the terminal (or) prints to a file if `-txto fileName.txt` is provided. 
``` 
mirtop sql -q --db examples/annotate/query_sample.db -e select -l 2 -col seqID,UID,Read,iso_5p,iso_3p,start,end
```
OUTPUT:
```
11/29/2019 03:28:40 INFO Run Convert GFF.
seqID   UID     Read    iso_5p  iso_3p  start   end
hsa-miR-342-3p  0t@TeV#5v2      TCTCACAAAGAAATCGCACCCGTTC       None    None    61.0    85.0
hsa-let-7i-5p   7AwwRIB71       TGAGGTAGTAGTTTGTGCTGTTG None    1.0     6.0     28.0
11/29/2019 03:28:40 INFO It took 0.000 minutes
``` 
* miRNA: Specify one or more miRNA to query. 

The user can specify miRNAs to query from the database. Use comma (miRNA-1,NO-SPACES,miRNA-n) to separate miRNAs while passing as an argument. For large set of query miRNAs use a text-file as input, separated by new line character. If short names are preffered over including the species name every time, then please refer to argument `-pm` or `--miRNA_prefix` to prefix the sepcies name along with the names of miRNAs. 
``` 
mirtop sql -q --db examples/annotate/query_sample.db -e select -l 4 -col seqID,UID,Read,iso_5p,iso_3p,start,end -miR hsa-let-7i-5p
```
OUTPUT:
```
seqID   UID     Read    iso_5p  iso_3p  start   end
hsa-let-7i-5p   7AwwRIB71       TGAGGTAGTAGTTTGTGCTGTTG None    1.0     6.0     28.0
hsa-let-7i-5p   7AwwRIBL1       TGAGGTAGTAGTTTGTGCTGTTA None    None    6.0     28.0
hsa-let-7i-5p   7AwwRIBU1       TGAGGTAGTAGTTTGTGCTGTTC None    None    6.0     28.0
hsa-let-7i-5p   7AwwRIBQ1       TGAGGTAGTAGTTTGTGCTGTTT None    None    6.0     28.0
11/29/2019 03:43:21 INFO It took 0.000 minutes
```
* variant: Specify the query with one or more variant types. The following variant types can be queried using `-var` argument. 

  * iso_5p                  - indicates the shift at the reference 5' miRNA
  * iso_3p                  - indicates the shift at the reference 3' miRNA
  * iso_add3p               - Number of non-template nucleotides added at 3p
  * iso_add5p               - Number of non-template nucleotides added at 5p
  * iso_snv_seed            - when affected nucleotides are between [2-7]
  * iso_snv_central_offset  - when affected nucleotide is at position [8]
  * iso_snv_central         - when affected nucleotides are between [9-12]
  * iso_snv_central_supp    - when affected nucleotides are between [13-17]
  * iso_snv                 - anything else

The conventional query for selecting rows with TRUE values of iso_5p, iso_3p and iso_snv_central_offset would be as "SELECT * FROM data_sets WHERE iso_5p!="None" AND iso_3p!="None" AND iso_snv_central_offset!=0".  In `mirtop sql` we can specifiy the same as shown in the example. 

```
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_3p,iso_snv_central_offset -l 5
```

OUTPUT:
```
12/02/2019 11:52:39 INFO Run Convert GFF.
seqID   source_file     type    start   end     score   strand  phase   UID     Read    Name    Parent  Variant iso_5p  iso_3p  iso_add3p       iso_add5p       iso_snv iso_snv_seed    iso_snv_central      iso_snv_central_offset  iso_snv_central_supp    source  cigar   hits    alias   genomic_pos     filter  seed_fam        SRR333680_1
hsa-miR-6727-5p miRBase22       isomiR  8.0     26.0    .       +       .       uMcGq6v2        TGGGGCAAGCGGCTGGCTC     hsa-miR-6727-5p hsa-mir-6727    iso_snv_central_offset,iso_5p:-2,iso_3p:-2   -2.0    -2.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       T6MA8MCTC       None    None    None    Pass    None    1
hsa-miR-6809-5p miRBase22       isomiR  8.0     25.0    .       +       .       xh@L$4  CCAAGGAAATAAGGGGAG      hsa-miR-6809-5p hsa-mir-6809    iso_snv_central_offset,iso_5p:-2,iso_3p:-2      -2.0 -2.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       C8MT3MG3MG      None    None    None    Pass    None    1
hsa-miR-6727-5p miRBase22       isomiR  8.0     24.0    .       +       .       hMGGDx1 AGGGGCCGGCGGCAGCC       hsa-miR-6727-5p hsa-mir-6727    iso_snv_central_offset,iso_5p:-2,iso_3p:-4      -2.0 -4.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       A5MC6MAMCC      None    None    None    Pass    None    1
hsa-let-7f-5p   miRBase22       isomiR  8.0     27.0    .       +       .       .       NAGGTAGTAGATTGTATAGT    hsa-let-7f-5p   hsa-let-7f-1    iso_snv_central_offset,iso_5p:-1,iso_3p:-1      -1.0 -1.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       N19M    None    None    None    Pass    None    1
12/02/2019 11:52:39 INFO It took 0.001 minutes
```
* filter: Filter attribute lets a user choose data query such that the attribute is 'Pass' for the reads from the miRNA sequencing is OK; 'Reject' if the reads are false positive; 'Reject lowcount'where the miRNA is rejected due to a low count in data. This filter decision is made by the aligner tools and is supplied to the GFF3. 


```
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_3p,iso_snv_central_offset -l 5 -f Pass
```
OUTPUT:
```
seqID   source_file     type    start   end     score   strand  phase   UID     Read    Name    Parent  Variant iso_5p  iso_3p  iso_add3p       iso_add5p       iso_snv iso_snv_seed    iso_snv_central      iso_snv_central_offset  iso_snv_central_supp    source  cigar   hits    alias   genomic_pos     filter  seed_fam        SRR333680_1
hsa-miR-6727-5p miRBase22       isomiR  8.0     26.0    .       +       .       uMcGq6v2        TGGGGCAAGCGGCTGGCTC     hsa-miR-6727-5p hsa-mir-6727    iso_snv_central_offset,iso_5p:-2,iso_3p:-2   -2.0    -2.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       T6MA8MCTC       None    None    None    Pass    None    1
hsa-miR-6809-5p miRBase22       isomiR  8.0     25.0    .       +       .       xh@L$4  CCAAGGAAATAAGGGGAG      hsa-miR-6809-5p hsa-mir-6809    iso_snv_central_offset,iso_5p:-2,iso_3p:-2      -2.0 -2.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       C8MT3MG3MG      None    None    None    Pass    None    1
hsa-miR-6727-5p miRBase22       isomiR  8.0     24.0    .       +       .       hMGGDx1 AGGGGCCGGCGGCAGCC       hsa-miR-6727-5p hsa-mir-6727    iso_snv_central_offset,iso_5p:-2,iso_3p:-4      -2.0 -4.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       A5MC6MAMCC      None    None    None    Pass    None    1
hsa-let-7f-5p   miRBase22       isomiR  8.0     27.0    .       +       .       .       NAGGTAGTAGATTGTATAGT    hsa-let-7f-5p   hsa-let-7f-1    iso_snv_central_offset,iso_5p:-1,iso_3p:-1      -1.0 -1.0    None    None    0.0     0.0     0.0     1.0     0.0     miRBase22       N19M    None    None    None    Pass    None    1
12/02/2019 12:33:28 INFO It took 0.001 minutes
```
* count: Count is used to retrieve a summary for a custom query. Generally, `select * from data_sets` or `select columns from data_sets` is used to retrieve the information and a WHERE clause is optionally supplied such as `--variants \| --miRNA` etc. This can be leveraged to get the count for a specific query, such that the query `select count(*) from data_sets` or `select count(columns) from data_sets` will be executed along with any optional WHERE clause. Argument `-n T` where T is True else False (Default False).

For example: 
1) How many miRNA isoforms exists for a variant type iso_snv_central_offset and iso_5p?
2) How many miRNA isoforms exists for a variant type iso_5p and iso_3p?
3) How many miRNA isoforms exists for a variant type iso_5p and iso_3p for miRNA hsa-miR-142-5p and hsa-miR-372-3p?

Query in that order of example is below: 
```
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_snv_central_offset -f Pass -n T
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_3p -f Pass -n T
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_3p -f Pass -l 30 -miR hsa-miR-142-5p,hsa-miR-372-3p -n T
```

OUTPUT: 

Example 1: 

```
12/02/2019 01:00:03 INFO Run Convert GFF.
COUNT(*)
Unique counts for all rows is: 10
12/02/2019 01:00:03 INFO It took 0.000 minutes

```

Example 2: 

```
12/02/2019 01:01:59 INFO Run Convert GFF.
COUNT(*)
Unique counts for all rows is: 751
12/02/2019 01:01:59 INFO It took 0.000 minutes
```

Example 3: 

```
12/02/2019 01:02:40 INFO Run Convert GFF.
COUNT(*)
1. hsa-miR-142-5p:      28
2. hsa-miR-372-3p:      46
12/02/2019 01:02:40 INFO It took 0.000 minutes
```

* txtout: Specify the query to redirect the output to a file instead of printing it on the screen. The name of the file must have an extension of (.txt). 

Extending from our previous example (3): How many miRNA isoforms exists for a variant type iso_5p and iso_3p for miRNA hsa-miR-142-5p and hsa-miR-372-3p and redirect the output to sample_count.txt
```
mirtop sql -q --db examples/annotate/query_sample.db -e select -var iso_5p,iso_3p -f Pass -l 30 -miR hsa-miR-142-5p,hsa-miR-372-3p -n T -txto sample_count.txt
```
OUTPUT:
```
12/02/2019 01:10:06 INFO Run Convert GFF.

Writing data to file: sample_count.txt

12/02/2019 01:10:06 INFO It took 0.000 minutes
```

