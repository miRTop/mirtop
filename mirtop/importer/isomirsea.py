""" Read isomiR GFF files"""

import os
from collections import defaultdict, Counter

import mirtop.libs.logger as mylog
from mirtop.mirna import mapper
from mirtop.mirna.realign import expand_cigar, make_id
from mirtop.gff.body import paste_columns, read_attributes
from mirtop.gff.body import variant_with_nt
from mirtop.gff.classgff import feature

logger = mylog.getLogger(__name__)


def header(fn):
    """
    Custom header for isomiR-SEA importer.

    Args:
        *fn (str)*: file name with isomiR-SEA GFF output

    Returns:
        *(str)*: isomiR-SEA header string.
    """
    h = ""
    return h


def read_file(fn, args):
    """
    Read isomiR-SEA file and convert to mirtop GFF format.

    Args:
        *fn(str)*: file name with isomiR-SEA output information.

        *database(str)*: database name.

        *args(namedtuple)*: arguments from command line.
            See *mirtop.libs.parse.add_subparser_gff()*.

    Returns:
        *reads (nested dicts)*:gff_list has the format as
            defined in *mirtop.gff.body.read()*.

    """
    database = args.database
    gtf = args.gtf
    sep = " " if args.out_format == "gtf" else "="
    map_mir = mapper.read_gtf_to_mirna(gtf)
    reads = defaultdict(dict)
    reads_in = 0
    sample = os.path.splitext(os.path.basename(fn))[0]
    hits = _get_hits(fn)
    logger.debug("ISOMIRSEA::SAMPLE::%s" % sample)
    with open(fn) as handle:
        for line in handle:
            cols = line.strip().split("\t")
            attr = read_attributes(line, "=")
            query_name = attr['TS']
            query_sequence = attr['TS'].replace("U", "T")
            start = int(cols[3])
            end = int(cols[4])
            isomirseq_iso = attr['ISO']
            if query_name not in reads and query_sequence == None:
                continue
            if query_sequence and query_sequence.find("N") > -1:
                continue
            counts = attr["TC"]
            chrom = cols[0]
            # logger.debug("SEQBUSTER:: cigar {cigar}".format(**locals()))
            cigar = attr['CI'].replace("U", "T")
            idu = make_id(query_sequence)
            isoformat = cigar2variants(cigar, query_sequence, attr['ISO'])
            logger.debug("\nISOMIRSEA::NEW::query: {query_sequence}\n"
                         "  precursor {chrom}\n"
                         "  name: {query_name}\n"
                         "  idu: {idu}\n"
                         "  start: {start}\n"
                         "  cigar: {cigar}\n"
                         "  iso: {isoformat}\n"
                         "  variant: {isoformat}".format(**locals()))
            source = "isomiR" if isoformat != "NA" else "ref_miRNA"
            strand = "+"
            database = cols[1]
            mirName = attr['MIN'].split()[0]
            preName = attr['PIN'].split()[0]
            score = "."
            Filter = attr['FILTER']
            isotag = attr['ISO']
            tchrom, tstart = _genomic2transcript(map_mir[mirName],
                                                 chrom, start)
            start = start if not tstart else tstart
            chrom = chrom if not tstart else tchrom
            end = start + len(query_sequence)
            hit = hits[idu]
            attrb = ("Read {query_sequence}; UID {idu}; Name {mirName};"
                     " Parent {preName}; Variant {isoformat};"
                     " Isocode {isotag}; Cigar {cigar}; Expression {counts};"
                     " Filter {Filter}; Hits {hit};").format(**locals())
            line = ("{chrom}\t{database}\t{source}\t{start}\t{end}\t"
                    "{score}\t{strand}\t.\t{attrb}").format(**locals())
            if args.add_extra:
                extra = variant_with_nt(line, args.precursors, args.matures)
                line = "%s Changes %s;" % (line, extra)

            line = paste_columns(feature(line), sep=sep)
            if start not in reads[chrom]:
                reads[chrom][start] = []
            if Filter == "Pass":
                reads_in += 1
                reads[chrom][start].append([idu, chrom, counts, sample, line])

    logger.info("Hits: %s" % reads_in)
    return reads


def _get_hits(fn):
    hits = Counter()
    with open(fn) as handle:
        for line in handle:
            attr = read_attributes(line, "=")
            query_sequence = attr['TS'].replace("U", "T")
            if query_sequence and query_sequence.find("N") > -1:
                continue
            idu = make_id(query_sequence)
            hits[idu] += 1
    return hits


def cigar2variants(cigar, sequence, tag):
    """From cigar to Variants in GFF format"""
    pos = 0
    iso5p = 0
    logger.debug("\nISOMIRSEA:: expanded: %s" % expand_cigar(cigar))
    for l in expand_cigar(cigar):
        if l == "I":
            iso5p -= 1
        elif l == "D":
            iso5p += 1
        else:
            break
    iso3p = 0
    for l in reversed(expand_cigar(cigar)):
        if l == "I":
            iso3p += 1
        elif l == "D":
            iso3p -= 1
        else:
            break
    isosnp = []
    for l in expand_cigar(cigar):
        if l in ['A', 'T', 'C', 'G' ]:
            isosnp.append([pos, sequence[pos], l])
        if l in ['D']:
            continue
        pos += 1
    iso5p = "iso_5p:%s" % _fix(iso5p) if iso5p else ""
    if tag[-1] == "T" or iso3p < 0:
        iso3p = "iso_3p:%s" % _fix(iso3p) if iso3p else ""
    else:
        iso3p = "iso_add3p:%s" % iso3p if iso3p else ""

    variant = ""
    for iso in [iso5p, iso3p, _define_snp(isosnp)]:
        if iso:
            variant += "%s," % iso

    variant = "NA;" if not variant else variant
    return variant[:-1]


def _define_snp(subs):
    value = ""
    logger.debug("\nISOMIRSEA:: subs %s" % subs)
    for sub in subs:
        if sub:
            if sub[0] > 1 and sub[0] < 8:
                value += "iso_snv_seed,"
            elif sub[0] == 8:
                value += "iso_snv_central_offset,"
            elif sub[0] > 8 and sub[0] < 13:
                value += "iso_snv_central,"
            elif sub[0] > 12 and sub[0] < 18:
                value += "iso_snv_central_supp,"
            else:
                value += "iso_snv,"
    return value[:-1]


def _fix(n):
    if n > 0:
        return "+%s" % n
    return n


def _genomic2transcript(code, chrom, pos):
    for ref in code:
        if _is_chrom(chrom, code[ref][0]):
            if _is_inside(pos, code[ref][1:3]):
                return [ref, _transcript(pos, code[ref][1:4])]
    return [None, None]


def _is_chrom(chrom, annotated):
    logger.debug("TRANSCRIPT::CHROM::read position %s and db position %s" % (chrom, annotated))
    if chrom == annotated:
        return True
    if chrom == annotated.replace("chr", ""):
        return True
    return False


def _is_inside(pos, annotated):
    logger.debug("TRANSCRIPT::INSIDE::read position %s and db position %s" % (pos, annotated))
    if pos > annotated[0] and pos < annotated[1]:
        return True
    return False


def _transcript(pos, annotated):
    logger.debug("TRANSCRIPT::TRANSCRIPT::read position %s and db position %s" % (pos, annotated))
    if annotated[2] == "+":
        return pos - annotated[0]
    elif annotated[2] == "-":
        return annotated[1] - pos
    raise ValueError("Strand information is incorrect %s" % annotated[3])
