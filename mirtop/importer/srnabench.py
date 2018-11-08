""" Read sRNAbench files"""

import os
from collections import defaultdict

import mirtop.libs.logger as mylog
from mirtop.gff.body import paste_columns, variant_with_nt
from mirtop.mirna.realign import make_cigar, make_id
from mirtop.gff.classgff import feature

logger = mylog.getLogger(__name__)


def read_file(folder, args):
    """
    Read sRNAbench file and convert to mirtop GFF format.

    Args:
        *fn(str)*: file name with sRNAbench output information.

        *database(str)*: database name.

        *args(namedtuple)*: arguments from command line.
            See *mirtop.libs.parse.add_subparser_gff()*.

    Returns:
        *reads (nested dicts)*:gff_list has the format as
            defined in *mirtop.gff.body.read()*.

    """
    reads_anno = os.path.join(folder, "reads.annotation")
    reads_iso = os.path.join(folder, "microRNAannotation.txt")
    sep = " " if args.out_format == "gtf" else "="
    sample = os.path.basename(folder)
    database = args.database
    precursors = args.precursors
    matures = args.matures

    n_out = 0
    n_in = 0
    n_ns = 0
    n_notassign = 0
    n_notindb = 0
    reads = defaultdict(dict)
    seen = set()

    source_iso = _read_iso(reads_iso)
    logger.info("Reads with isomiR information %s" % len(source_iso))
    with open(reads_anno) as handle:
        for sequence in handle:
            cols = sequence.strip().split("\t")
            query_name = cols[0]
            query_sequence = cols[0]
            if query_name not in reads and not query_sequence:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                n_ns += 1
                continue
            if cols[3].find("mature") == -1:
                n_in += 1
                continue

            counts = int(cols[1])

            hit = len(set([mirna.split("#")[1] for mirna in cols[4].split("$")]))

            for nhit in cols[4].split("$"):
                logger.debug("SRNABENCH::line hit: %s" % nhit)
                hit_info = nhit.split("#")
                pos_info = hit_info[3].split(",")
                start = int(pos_info[1]) - 1
                end = start + len(query_sequence)  # int(pos_info[2]) - 1
                chrom = pos_info[0]
                mirName = hit_info[1]
                if chrom not in precursors or chrom not in matures:
                    n_notindb += 1
                if mirName not in matures[chrom]:
                    n_notindb += 1
                if (query_sequence, mirName) in seen:
                    continue

                seen.add((query_sequence, mirName))

                if (query_sequence, mirName) not in source_iso:
                    continue

                isoformat = source_iso[(query_sequence, mirName)]

                if isoformat == "mv":
                    n_notassign += 1
                    continue

                source = "isomiR" if isoformat != "NA" else "ref_miRNA"

                logger.debug("SRNABENCH::query: {query_sequence}\n"
                             "  precursor {chrom}\n"
                             "  name:  {query_name}\n"
                             "  start: {start}\n"
                             "  external: {isoformat}\n"
                             "  hit: {hit}".format(**locals()))
                logger.debug("SRNABENCH:: start %s end %s" % (start, end))
                if len(precursors[chrom]) < start + len(query_sequence):
                    n_out += 1
                    continue

                Filter = "Pass"
                cigar = make_cigar(query_sequence,
                                   precursors[chrom][start:end])
                preName = chrom
                score = "."
                strand = "+"
                idu = make_id(query_sequence)
                attrb = ("Read {query_sequence}; UID {idu}; Name {mirName};"
                         " Parent {preName}; Variant {isoformat};"
                         " Cigar {cigar}; Expression {counts};"
                         " Filter {Filter}; Hits {hit};").format(**locals())
                line = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t"
                        "{score}\t{strand}\t.\t{attrb}").format(**locals())
                if args.add_extra:
                    extra = variant_with_nt(line, args.precursors,
                                            args.matures)
                    line = "%s Changes %s;" % (line, extra)

                line = paste_columns(feature(line), sep=sep)
                if start not in reads[chrom]:
                    reads[chrom][start] = []
                if Filter == "Pass":
                    n_in += 1
                    reads[chrom][start].append([idu, chrom, counts,
                                                sample, line])

    logger.info("Loaded %s reads with %s hits" % (len(reads), n_in))
    logger.info("Reads without precursor information: %s" % n_notindb)
    logger.info("Reads with MV as variant definition,"
                " not supported by GFF: %s" % n_notassign)
    logger.info("Hit Filtered by having > 3 changes: %s" % n_out)

    return reads


def _read_iso(fn):
    """
    Read definitions of isomiRs by srnabench
    """
    iso = dict()
    with open(fn) as inh:
        inh.readline()
        for line in inh:
            cols = line.strip().split("\t")
            label = cols[3].split("$")
            mirnas = cols[1].split("$")
            if len(mirnas) == 1 and len(label) > 1:
                label = [cols[3].split("$")[0]]
            if len(mirnas) != len(label):
                label = label * (len(mirnas) - len(label))
            anno = dict(zip(mirnas, label))
            logger.debug("TRANSLATE::%s with %s" % (mirnas, label))
            for m in anno:
                iso[(cols[0], m)] = _translate(anno[m], cols[4])
                logger.debug("TRANSLATE::code %s" % iso[(cols[0], m)])
    return iso


def _translate(isomirs, description):
    # TODO: make unit test for this
    iso = []
    labels = isomirs.split("@")
    logger.debug("TRANSLATE::label:%s" % isomirs)
    for label in labels:
        logger.debug("TRANSLATE::label:%s" % label)
        if label == "exact":
            return "NA"
        if label.find("mv") > -1:
            return "mv"
        number_nts = label.split("|")[-1].split("#")[-1]
        if number_nts.find("-") < 0:
            number_nts = "+%s" % number_nts
        if label.find("lv3p") > -1:
            iso.append("iso_3p:%s" % number_nts)
        if label.find("lv5p") > -1:
            if number_nts.startswith("+"):
                number_nts = number_nts.replace("+", "-")
            else:
                number_nts = number_nts.replace("-", "+")
            iso.append("iso_5p:%s" % number_nts)
        if label.find("nta") > -1:
            number_nts = label.split("|")[1].split("#")[-1]
            if number_nts.find("-") < 0:
                number_nts = "%s" % number_nts
            iso.append("iso_add3p:%s" % number_nts)
        if label.find("NucVar") > -1:
            for nt in description.split(","):
                logger.debug("TRANSLATE::change:%s" % description)
                if nt == "-" or nt == "NA":
                    return "notsure"
                iso.extend(_iso_snp(int(nt.split(":")[0])))
        logger.debug("TRANSLATE::iso:%s" % iso)
    return ",".join(iso)


def _iso_snp(pos):
    iso = []
    if pos > 1 and pos < 8:
        iso.append("iso_snv_seed")
    elif pos == 8:
        iso.append("iso_snv_central_offset")
    elif pos > 8 and pos < 13:
        iso.append("iso_snv_central")
    elif pos > 12 and pos < 18:
        iso.append("iso_snv_central_supp")
    else:
        iso.append("iso_snv")
    return iso
