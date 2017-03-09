---
title: "example_mirTop_output"
author: "Lorena Pantano"
date: "3/9/2017"
output: html_document
---





```r
library(readr)
library(dplyr)
read_tsv("~/repos/mirtop/test/test_automated_output/test_out_mirs_fasta/sim_isomir.mirna") %>% knitr::kable()
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



|name                       |seq                                                     | freq|chrom         |start |end |subs |add |t5   |t3 |s5 |s3 |DB    |precursor    | hits|Name                        |
|:--------------------------|:-------------------------------------------------------|----:|:-------------|:-----|:---|:----|:---|:----|:--|:--|:--|:-----|:------------|----:|:---------------------------|
|ATGAGGTAGAAGGTTGTATAGT     |hsa-let-7a-1_hsa-let-7a-5p_5:26_-1:-1_mut:10A_add:null  |    0|hsa-let-7a-5p |NA    |NA  |9AT  |0   |A    |t  |NA |NA |miRNA |hsa-let-7a-1 |    1|hsa-let-7a-5p.sA.T9A.t      |
|GAGGTAGTAGGTTGTATGGT       |hsa-let-7a-1_hsa-let-7a-5p_7:26_1:-1_mut:18G_add:null   |    0|hsa-let-7a-5p |NA    |NA  |17GA |0   |t    |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.st.A17G.t     |
|ATGAGGTAGTAGGTTGTATAGTT    |hsa-let-7a-1_hsa-let-7a-5p_5:27_-1:0_mut:null_add:null  |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |A    |0  |NA |NA |miRNA |hsa-let-7a-1 |    1|hsa-let-7a-5p.sA..NA        |
|TGAGGTAGTAGGTTGTATAGT      |hsa-let-7a-3_hsa-let-7a-5p_4:24_0:-1_mut:null_add:null  |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |0    |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sNA..t        |
|GGTAGTAGGTTGTATAGT         |hsa-let-7a-1_hsa-let-7a-5p_9:26_3:-1_mut:null_add:null  |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |tga  |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stga..t       |
|GGTGAGGTAGTAGGTTGTATAGTTTT |hsa-let-7a-3_hsa-let-7a-5p_2:24_-2:-1_mut:null_add:TTT  |    0|hsa-let-7a-5p |NA    |NA  |0    |T   |GG   |T  |NA |NA |miRNA |hsa-let-7a-3 |    1|hsa-let-7a-5p.sGG..T.eT     |
|AGAGGTAGTAGGTTGTATAGT      |hsa-let-7a-1_hsa-let-7a-5p_6:26_0:-1_mut:1A_add:null    |    0|hsa-let-7a-5p |NA    |NA  |0AT  |0   |0    |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sNA.T0A.t     |
|GAGGTAGTAGGTTGTATAGTATA    |hsa-let-7a-1_hsa-let-7a-5p_7:26_1:-1_mut:null_add:ATA   |    0|hsa-let-7a-5p |NA    |NA  |20AT |0   |t    |TA |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.st.T20A.TA    |
|TTGAGGTAGTAGGTTGTATAGTT    |hsa-let-7a-2_hsa-let-7a-5p_4:26_-1:0_mut:null_add:null  |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |T    |0  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sT..NA        |
|GTTGAGGTAGTAGGTTGTATAGT    |hsa-let-7a-2_hsa-let-7a-5p_3:25_-2:-1_mut:null_add:null |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |GT   |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sGT..t        |
|GTAGTAGGTTGTATAGTTAA       |hsa-let-7a-3_hsa-let-7a-5p_8:24_4:-1_mut:null_add:TAA   |    0|hsa-let-7a-5p |NA    |NA  |0    |AA  |tgag |0  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stgag..NA.eAA |
|GGTAGTAGGTTGTATAGTATT      |hsa-let-7a-2_hsa-let-7a-5p_8:25_3:-1_mut:null_add:ATT   |    0|hsa-let-7a-5p |NA    |NA  |0    |ATT |tga  |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stga..t.eATT  |
|GAGGTAGTAGGTTGTATAGTTAT    |hsa-let-7a-2_hsa-let-7a-5p_6:25_1:-1_mut:null_add:TAT   |    0|hsa-let-7a-5p |NA    |NA  |0    |AT  |t    |0  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.st..NA.eAT    |
|GAGGTAGTAGGTTGTATAGTT      |hsa-let-7a-3_hsa-let-7a-5p_5:24_1:-1_mut:null_add:T     |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |t    |0  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.st..NA        |
|GGTAGTAGGTTGTATAGTTTT      |hsa-let-7a-1_hsa-let-7a-5p_9:26_3:-1_mut:null_add:TTT   |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |tga  |TT |NA |NA |miRNA |hsa-let-7a-1 |    1|hsa-let-7a-5p.stga..TT      |
|AGGTAGTAGGTTGTATAGTT       |hsa-let-7a-2_hsa-let-7a-5p_7:26_2:0_mut:null_add:null   |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |tg   |0  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stg..NA       |
|GGTGAGGTAGTAGGTTGTATAGT    |hsa-let-7a-3_hsa-let-7a-5p_2:24_-2:-1_mut:null_add:null |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |GG   |t  |NA |NA |miRNA |hsa-let-7a-3 |    1|hsa-let-7a-5p.sGG..t        |
|TGAGGTAGTAGGTTGTATAG       |hsa-let-7a-2_hsa-let-7a-5p_5:24_0:-2_mut:null_add:null  |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |0    |tt |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sNA..tt       |
|AGGTAGTAGGTTGGATAGTTTA     |hsa-let-7a-2_hsa-let-7a-5p_7:25_2:-1_mut:14G_add:TTA    |    0|hsa-let-7a-5p |NA    |NA  |13GT |0   |tg   |TA |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stg.T13G.TA   |
|GTGAGGTAGTAGGTTGTATAGT     |hsa-let-7a-3_hsa-let-7a-5p_3:24_-1:-1_mut:null_add:null |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |G    |t  |NA |NA |miRNA |hsa-let-7a-3 |    1|hsa-let-7a-5p.sG..t         |
|TGAGGTAGTAGGTTGTATAGAT     |hsa-let-7a-2_hsa-let-7a-5p_5:24_0:-2_mut:null_add:AT    |    0|hsa-let-7a-5p |NA    |NA  |0    |AT  |0    |tt |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.sNA..tt.eAT   |
|AGGTAGTAGGTTGTATAGT        |hsa-let-7a-1_hsa-let-7a-5p_8:25_2:-2_mut:null_add:T     |    0|hsa-let-7a-5p |NA    |NA  |0    |0   |tg   |t  |NA |NA |miRNA |hsa-let-7a-2 |    1|hsa-let-7a-5p.stg..t        |
|GATGAGGTAGTAGGTTGTATAGTA   |hsa-let-7a-1_hsa-let-7a-5p_4:26_-2:-1_mut:null_add:A    |    0|hsa-let-7a-5p |NA    |NA  |0    |A   |GA   |t  |NA |NA |miRNA |hsa-let-7a-1 |    1|hsa-let-7a-5p.sGA..t.eA     |

