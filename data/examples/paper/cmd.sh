java -jar ~/repos/seqbuster/modules/miraligner/miraligner.jar -sub 1 -trim 3 -add 3 -s hsa -i hsa-let-7a-5p.fa -db ../annotat
e -o hsa-let-7a-5p

 mirtop gff --hairpin ../annotate/hairpin.fa --gtf ../annotate/hsa.gff3 --out . --format seqbuster hsa-let-7a-5p.mirna
