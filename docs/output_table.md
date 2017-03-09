---
title: "example_mirTop_output"
author: "Lorena Pantano"
date: "3/9/2017"
output: html_document
---





```r
library(readr)
library(dplyr)
read_tsv("~/repos/mirtop/test/test_automated_output/test_out_mirs_fasta/sim_isomir.mirna") %>% DT::datatable()
```

```
## Parsed with column specification:
## cols(
##   name = col_character(),
##   seq = col_character(),
##   freq = col_integer(),
##   chrom = col_character(),
##   start = col_character(),
##   end = col_character(),
##   subs = col_character(),
##   add = col_character(),
##   t5 = col_character(),
##   t3 = col_character(),
##   s5 = col_character(),
##   s3 = col_character(),
##   DB = col_character(),
##   precursor = col_character(),
##   hits = col_integer(),
##   Name = col_character()
## )
```

```
## Error in loadNamespace(name): there is no package called 'webshot'
```

