""" Read GFF files and output isomiRs compatible format"""
from __future__ import print_function

import os

import mirtop.libs.logger as mylog
from mirtop.mirna import fasta, mapper
from mirtop.gff.classgff import feature
from mirtop.gff.header import read_samples
from mirtop.mirna.realign import get_mature_sequence, align_from_variants
from mirtop.mirna.realign import read_id, variant_to_5p, \
                                 variant_to_3p, variant_to_add
from mirtop.gff.body import variant_with_nt

logger = mylog.getLogger(__name__)


def convert(args):
    """
    Main function to convert from GFF3 to isomiRs Bioc Package.

    Reads a GFF file to produces output file containing Expression counts

    Args:
        *args(namedtuple)*: arguments parsed from command line with
            *mirtop.libs.parse.add_subparser_counts()*.

    Returns:
        *file (file)*: with columns like:
            UID miRNA Variant Sample1 Sample2 ... Sample N
    """
    logger.info("INFO Writing TSV file to directory %s", args.out)
    for gff in args.files:
        logger.info("INFO Reading GFF file %s", gff)
        _convert_file(gff, args)


def _convert_file(gff, args):
    sep = "\t"
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    variant_header = sep.join(['mism', 'add', 't5', 't3'])

    gff_file = open(gff, 'r')
    out_file = os.path.join(args.out, "%s_rawData.tsv" % os.path.splitext(os.path.basename(gff))[0])
    missing_parent = 0
    missing_mirna = 0
    unvalid_uid = 0
    with open(out_file, 'w') as outh:

        for samples_line in gff_file:
            if samples_line.startswith("## COLDATA:"):
                samples = sep.join(samples_line.strip().split("COLDATA:")[1].strip().split(","))
                header = sep.join(['seq', 'mir',
                                   variant_header, samples])
                print(header, file=outh)
                break

        for mirna_line in gff_file:
            gff = feature(mirna_line)
            attr = gff.attributes
            UID = attr["UID"]
            Read = attr["Read"]
            mirna = attr["Name"]
            parent = attr["Parent"]
            variant = attr["Variant"]
            try:
                Read = read_id(UID)
            except KeyError:
                unvalid_uid += 1
                continue

            expression = sep.join(attr["Expression"].strip().split(","))
            cols_variants = sep.join(_expand(variant))
            logger.debug("COUNTS::Read:%s" % Read)
            logger.debug("COUNTS::EXTRA:%s" % variant)
            if parent not in precursors:
                missing_parent += 1
                continue
            if mirna not in matures[parent]:
                missing_mirna += 1
                continue
            extra = variant_with_nt(mirna_line, precursors, matures)
            if extra == "Invalid":
                continue
            logger.debug("COUNTS::EXTRA:%s" % extra)
            cols_variants = sep.join(_expand(extra, True))
            summary = sep.join([Read,  mirna,
                                cols_variants, expression])
            logger.debug(summary)
            print(summary, file=outh)

    gff_file.close()
    logger.info("Missing Parents in hairpin file: %s" % missing_parent)
    logger.info("Missing MiRNAs in GFF file: %s" % missing_mirna)
    logger.info("Non valid UID: %s" % unvalid_uid)
    logger.info("Output file is at %s" % out_file)


def _expand(variant, nts=False):
    """Expand Variant field into list for iso_5p, iso_3p, iso_add3p, iso_snv"""
    list_variant = []
    isomir = {}
    snp_var = []
    for v in variant.split(","):
        if v.find(":") > 0:
            isomir[v.split(":")[0]] = v.split(":")[1]
        elif v.find("snv") > 0:
            snp_var.append(1)
    if nts:
        if "iso_snv" in isomir:
            list_variant.append(isomir["iso_snv"])
        else:
            list_variant.append(0)
    else:
        snp = sum(snp_var)
        list_variant.append(snp)
    if "iso_add3p" in isomir:
        list_variant.append(isomir["iso_add3p"])
    else:
        list_variant.append(0)
    if "iso_5p" in isomir:
        list_variant.append(isomir["iso_5p"])
    else:
        list_variant.append(0)
    if "iso_3p" in isomir:
        list_variant.append(isomir["iso_3p"])
    else:
        list_variant.append(0)

    return list(map(str, list_variant))
