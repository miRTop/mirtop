current = "1.1"
GFFv = {'1.0': ["iso_5p", "iso_3p", "iso_add", "iso_snp_seed",
                "iso_snp_central_offset", "iso_snp_central",
                "iso_snp_central_supp", "iso_snp", "NA"],
        '1.1': ["iso_5p", "iso_3p", "iso_add3p", "iso_snv_seed",
                "iso_snv_central_offset", "iso_snv_central",
                "iso_snv_central_supp", "iso_snv", "NA"]}
ATTRv = {'1.0': set(('IDU', 'Name', 'Alias', 'Parent',
                     'Read', 'Variant', 'Changes',
                     'Cigar', 'Hits', 'Genomic',
                     'Expression', 'Filter', 'Cigar',
                     'Seed_fam')),
         '1.1': set(('IDU', 'Name', 'Alias', 'Parent',
                     'Read', 'Variant', 'Changes',
                     'Cigar', 'Hits', 'Genomic',
                     'Expression', 'Filter', 'Cigar',
                     'Seed_fam'))}
