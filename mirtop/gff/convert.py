import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam

from mirtop.mirna import fasta, mapper
from mirtop.gff.body import read_gff_line, variant_with_nt
from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def convert_gff_counts(args):
	""" Reads a GFF file, produces output file containing Expression counts of
	the format : UID miRNA Variant Sample1 Sample2 ... Sample N

	"""
        sep = "\t"
        variant_header = sep.join(['iso_5p', 'iso_3p',
                                    'iso_add', 'iso_snp'])
        if args.add_extra:
            database = mapper.guess_database(args.gtf)
            precursors = fasta.read_precursor(args.hairpin, args.sps)
            matures = mapper.read_gtf_to_precursor(args.gtf)
            variant_header = sep.join([variant_header,
                                        'iso_5p_nt', 'iso_3p_nt',
                                        'iso_add_nt', 'iso_snp_nt'])

	logger.info("INFO Reading GFF file %s", args.gff)
	logger.info("INFO Writing TSV file to directory %s", args.out)

	gff_file = open(args.gff, 'r')
	out_file = op.join(args.out, "expression_counts.tsv")
        with open(out_file, 'w') as outh:

            for samples_line in gff_file:
                if samples_line.startswith("## COLDATA:"):
                    samples = sep.join(samples_line.strip().split("COLDATA:")[1].strip().split(","))
                    header = sep.join(['UID', 'Read', 'miRNA', 'Variant', variant_header, samples])
                    print >>outh, header
                    break

            for mirna_line in gff_file:
                mirna_values = read_gff_line(mirna_line)
                Read = mirna_values["attrb"]["Read"]
                UID = mirna_values["attrb"]["UID"]
                mirna = mirna_values["attrb"]["Name"]
                variant = mirna_values["attrb"]["Variant"]
                expression = sep.join(mirna_values["attrb"]["Expression"].strip().split(","))
                cols_variants = sep.join(_expand(variant))
                if args.add_extra:
                    extra = variant_with_nt(mirna_line, precursors, matures)
                    logger.debug("COUNTS::EXTRA:%s" % extra)
                    cols_variants = sep.join([cols_variants] + _expand(extra, True))
                summary = sep.join([UID, Read,  mirna, variant, cols_variants, expression])
                logger.debug(summary)
                print >>outh, summary

	gff_file.close()
        logger.info("Output file is at %s" % out_file)


def _expand(variant, nts = False):
    """Expand Variant field into list for iso_5p, iso_3p, iso_add, iso_snp"""
    list_variant = []
    isomir = {}
    snp_var = []
    for v in variant.split(","):
        if v.find(":") > 0:
            isomir[v.split(":")[0]] = v.split(":")[1]
        elif v.find("snp") > 0:
            snp_var.append(1)

    if "iso_5p" in isomir:
        list_variant.append(isomir["iso_5p"])
    else:
        list_variant.append(0)
    if "iso_3p" in isomir:
        list_variant.append(isomir["iso_3p"])
    else:
        list_variant.append(0)
    if "iso_add" in isomir:
        list_variant.append(isomir["iso_add"])
    else:
        list_variant.append(0)
    if nts:
        if "iso_snp" in isomir:
            list_variant.append(isomir["iso_snp"])
        else:
            list_variant.append(0)
    else:
        snp = sum(snp_var)
        list_variant.append(snp)
    return map(str, list_variant)


