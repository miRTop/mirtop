import re
from Bio import pairwise2
from Bio.Seq import Seq
from collections import defaultdict

from mirtop.mirna.mintplates import convert
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


class hits:
    """"Class with alignment information."""

    def __init__(self):
        self.sequence = ""
        self.idseq = ""
        self.precursors = defaultdict(isomir)
        self.score = []
        self.best_hits = []  # maybe sam object?
        self.counts = 0

    def set_sequence(self, seq):
        self.sequence = seq
        self.idseq = make_id(seq)

    def set_precursor(self, precursor, isomir):
        self.precursors[(precursor, isomir.start)] = isomir

    def remove_precursor(self, precursor):
        del self.precursors[precursor]


class isomir:
    """
    Class to represent isomiRs information.
    """

    def __init__(self):
        self.t5 = []
        self.t3 = []
        self.add = []
        self.subs = []
        self.external = ""
        self.align = None
        self.cigar = None
        self.filter = "Pass"
        self.map_score = 0
        self.end = None
        self.start = None
        self.mirna = None
        self.strand = "+"

    def set_pos(self, start, l, strand="+"):
        """Set end position"""
        if start < 0:
            # l = l + start
            start = 0
        self.strand = strand
        self.start = start
        self.end = start + l - 1
        if strand == "-":
            self.start = start + l - 1
            self.end = start

    def formatGFF(self):
        """Create Variant attribute."""
        value = []
        subs = self.subs
        if self.external != "notsure" and self.external != "":
            return self.external
        for sub in subs:
            if sub:
                if sub[0] > 1 and sub[0] < 8:
                    value.append("iso_snv_seed")
                elif sub[0] == 8:
                    value.append("iso_snv_central_offset")
                elif sub[0] > 8 and sub[0] < 13:
                    value.append("iso_snv_central")
                elif sub[0] > 12 and sub[0] < 18:
                    value.append("iso_snv_central_supp")
                else:
                    value.append("iso_snv")

        if self.add:
            value.append("iso_add3p:%s" % len(self.add))
        if self.t5:
            size = len(self.t5)
            direction = "-" if self.t5.isupper() else "+"
            value.append("iso_5p:%s%s" % (direction, size))
        if self.t3:
            size = len(self.t3)
            direction = "+" if self.t3.isupper() else "-"
            value.append("iso_3p:%s%s" % (direction, size))
        if not value:
            value = ["NA"]
        return ",".join(list(set(value)))

    def format(self, sep="\t"):
        """Create tabular line from variant fields."""
        subs = "".join(["".join(map(str, mism)) for mism in self.subs])
        if not subs:
            subs = "0"
        add = "0" if not self.add else self.add
        return "%s%s%s%s%s%s%s" % (subs, sep, add, sep,
                                   self.t5, sep, self.t3)

    def format_id(self, sep="\t"):
        """Create simple identifier from variant fields."""
        subs = ["".join(["".join([c[2], str(c[0]), c[1]]) for c in self.subs])]
        if not subs:
            subs = []
        add = [] if not self.add else ["e%s" % self.add]
        t5 = ["s%s" % self.t5] if self.t5 and self.t5 != "NA" else []
        t3 = ["%s" % self.t3] if self.t3 and self.t3 != "NA" else []
        full = t5 + subs + t3 + add
        return sep.join([f for f in full if f])

    def get_score(self, sc):
        """Get score from variant fields."""
        for a in self.add:
            if a in ['A', 'T']:
                sc -= 0.25
            else:
                sc -= 0.75
        for e in self.subs:
            sc -= 1
        return sc

    def is_iso(self):
        """Define whether element is isomiR or not."""
        if self.external == "NA":
            return False
        if self.t5 or self.t3 or self.add or self.subs or self.external != "":
            return True
        return False


def read_id(idu):
    """
    Read a unique identifier for the sequence and
    convert it to the nucleotides,
    replacing an unique code for 5 nts.

    It uses the code from *mirtop.mirna.keys()*.

    Inspired by MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates

    Args:
        *idu(str)*: unique identifier for the sequence.

    Returns:
        *seq(str)*: nucleotides sequences.
    """
    try:
        seq = convert(idu, False, 'iso')
    except KeyError:
        logger.error("UID is not valid " + idu)
        return False

    return seq


def make_id(seq):
    """
    Create a unique identifier for the sequence from the nucleotides,
    replacing 5 nts for a unique sequence.

    It uses the code from *mirtop.mirna.keys()*.

    Inspired by MINTplate: https://cm.jefferson.edu/MINTbase
    https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates

    Args:
        *seq(str)*: nucleotides sequences.

    Returns:
        *idName(str)*: unique identifier for the sequence.
    """
    try:
        idu = convert(seq, True, 'iso')
    except KeyError as error:
        logger.error("Sequence is not valid " + seq)
        raise

    # If you wanted to add "iso-" into the license plate as the prefix
    # idu = convert(seq, True, "iso")
    return idu


def is_sequence(seq):
    """
    This function check whether the sequence is valid or not.

    Args:
        *seq(str)*: string acting as a sequence.

    Returns:
        *boolean*: whether is or not a valid nucleotide sequence.
    """
    alphabet = re.compile('^[ACTG]*$', re.IGNORECASE)
    return alphabet.match(seq)


def align(x, y, local=False):
    """
    Pairwise alignments between two sequenes.
    https://medium.com/towards-data-science/pairwise-sequence-alignment-using-biopython-d1a9d0ba861f

    Args:
        *x(str)*: short sequence.

        *y(str)*: long sequence.

        *local(boolean)*: local or global alignment.

    Returns:
        *aligned_x(hit)*: alignment information, socre and positions.
    """
    if local:
        aligned_x = pairwise2.align.localxx(x, y)[0]
    else:
        aligned_x = pairwise2.align.globalms(x, y, 1, -1, -1, -0.5)[0]
    aligned_x = list(aligned_x)
    n_x = aligned_x[0]
    if "N" in n_x:
        N_indices = [i for i, ltr in enumerate(n_x) if ltr == 'N']
        n_x = list(n_x)
        for N_index in N_indices:
            n_x[N_index] = y[N_index]
        n_x = ''.join(n_x)

    aligned_x[0] = n_x
    return tuple(aligned_x)


def _add_cigar_char(counter, cache):
    if counter == 1:
        return cache
    else:
        return str(counter) + cache


def make_cigar(seq, mature):
    """
    Function that will create CIGAR string from aligment
    between read and reference sequence.

    Args:
        *seq(str)*: read sequence.

        *mature(str)*: short sequence.

    Return:
        *short(str)*: CIGAR string.
    """
    cigar = ""
    for pos in range(0, len(seq)):
        if seq[pos] == mature[pos]:
            cigar += "M"
        elif seq[pos] != mature[pos] and seq[pos] != "-" and mature[pos] != "-":
            cigar += seq[pos]
        elif seq[pos] == "-":
            cigar += "D"
        elif mature[pos] == "-":
            cigar += "I"

    cache = ""
    counter = 1
    short = ""
    for c in cigar:
        if c != cache and cache != "" and cache == "M":
            short += _add_cigar_char(counter, cache)
            counter = 1
            cache = c
        if c != "M":
            short += c
        if c == cache and c == "M":
            counter += 1
        cache = c

    if cache == "M":
        short += _add_cigar_char(counter, cache)
    return short


def cigar_correction(cigarLine, query, target):
    """
    Read from CIGAR in BAM file to define mismatches.

    Args:
        *cirgarLine(str)*: CIGAR string from BAM file.

        *query(str)*: read sequence.

        *target(str)*: target sequence.

    Returns:
        *(list)*: [query_nts, target_nts]
    """
    query_pos = 0
    target_pos = 0
    query_fixed = []
    target_fixed = []
    for (cigarType, cigarLength) in cigarLine:
        if cigarType == 0:  # match
            query_fixed.append(query[query_pos:query_pos+cigarLength])
            target_fixed.append(target[target_pos:target_pos+cigarLength])
            query_pos = query_pos + cigarLength
            target_pos = target_pos + cigarLength
        elif cigarType == 1:  # insertions
            query_fixed.append(query[query_pos:query_pos+cigarLength])
            target_fixed.append("".join(["-"] * cigarLength))
            query_pos = query_pos + cigarLength
        elif cigarType == 2:  # deletion
            target_fixed.append(target[target_pos:target_pos+cigarLength])
            query_fixed.append("".join(["-"] * cigarLength))
            target_pos = target_pos + cigarLength
    return ["".join(query_fixed), "".join(target_fixed)]


def expand_cigar(cigar):
    """
    From short CIGAR version to long CIGAR version
    where each character is each nts in the sequence.

    Args:
        *cigar(str)*: CIGAR string.

        >>> 10MA3M

    Returns:
        *cigar_long(str)*: CIGAR long.

        >>> MMMMMMMMMMAMMM
    """
    cigar_long = ""
    n = 0
    for nt in cigar:
        if nt in ["D", "M", "I", "A", "T", "C", "G"]:
            if n > 0:
                cigar_long += nt * int(n)
            else:
                cigar_long += nt
            n = 0
        else:
            if n > 0:
                n = int("%s%s" % (n, nt))
            else:
                n = int(nt)
    return cigar_long


def cigar2snp(cigar, reference):
    """
    From a CIGAR string and reference sequence
    detect mistmatches positions and reference and
    target nucleotides.

    Args:
        *cigar(str)*: CIGAR string.

        *reference(str)*: reference sequence.

    Returns:
        *snp(list)*: position of mismatches (indels included) as:

        >>> [pos, seq_nt, ref_nt]
    """
    snp = []
    pos_seq = 0
    pos_ref = 0
    for nt in expand_cigar(cigar):
        if nt != "M":
            if nt == "I":
                snp.append([pos_seq, nt, "-"])
                pos_seq += 1
            elif nt == "D":
                snp.append([pos_seq, "-", reference[pos_ref]])
                pos_ref += 1
            else:
                snp.append([pos_seq, nt, reference[pos_ref]])
                pos_ref += 1
                pos_seq += 1
        else:
            pos_ref += 1
            pos_seq += 1
    return snp


def reverse_complement(seq):
    """Get reverse complement of a sequences

    Args:
        *seq(str)*: sequence.

        >>> GCAT

    Returns:
        *(str)*: reverse complemente sequence:

        >>> ATGC
    """
    return str(Seq(seq).reverse_complement())


def get_mature_sequence(precursor, mature, exact=False, nt = 5):
    """
    From precursor and mature positions
       get mature sequence with +/- 4 flanking nts.

    Args:
       *precursor(str)*: long sequence.

       *mature(list)*: [start, end].

       *exact(boolean)*: not add 4+/- flanking nts.

       *nt(int)*: number of nts to get.

    Returns:
        *(str)*: mature sequence.
    """
    p = "%s%s" % ("".join(["N"]*nt), precursor)
    s = mature[0] + nt
    e = mature[1] + nt
    if exact:
        return p[s:e + 1]
    return p[s - (nt - 1) :e + nt]


def align_from_variants(sequence, mature, variants):
    """
    Giving the sequence read,
       the mature from get_mature_sequence,
       and the variant GFF annotation:
       get a list of substitutions

    Args:
        *sequence(str)*: read sequence.

        *mature(str)*: mature sequence from
            *mirtop.mirna.realing.get_mature_sequence()*.

        *variants(str)*: string from Variant attribute in GFF file.

    Returns:
        *snp(list)*: [[pos, target, reference]]
    """
    init_log = "iso:%s -> %s\nref:%s" % (sequence, variants, mature)
    snps = []
    k = [v.split(":")[0] for v in variants.split(",") if v.find(":") > -1]
    v = [int(v.split(":")[1]) for v in variants.split(",") if v.find(":") > -1]
    var_dict = dict(zip(k, v))
    logger.debug("realign::align_from_variants::sequence %s" % sequence)
    logger.debug("realign::align_from_variants::mature %s" % mature)
    logger.debug("realign::align_from_variants::variants %s" % variants)
    snp = ["iso_snv" for v in variants.split(",") if v.find("snv") > -1]
    fix_5p = 7
    if "iso_5p" in k:
        fix_5p = 7 + var_dict["iso_5p"]
    mature = mature[fix_5p:]
    if "iso_add3p" in k:
        sequence = sequence[:-1 * var_dict["iso_add3p"]]
    if "iso_3p" in k:
        mature = mature[:-(7 + (-1 * var_dict["iso_3p"]))]
    else:
        mature = mature[:-7]
    logger.debug("realign::align_from_variants::snp %s" % snp)
    logger.debug("realign::align_from_variants::sequence %s" % sequence)
    logger.debug("realign::align_from_variants::mature %s" % mature)

    if len(sequence) != len(mature):  # in case of indels, align again
        a = align(sequence, mature)
        sequence = a[0]
        mature = a[1]

    if len(sequence) > len(mature):
        logger.warning("Invalid isomiR definition:\n%s\niso:%s\nref:%s" % (init_log, sequence, mature))
        return "Invalid"

    for p in range(0, len(sequence)):
        if sequence[p] != mature[p]:
            if mature[p] == "N":
                continue
            value = ""
            pos = p + 1
            value = "iso_snv"
            logger.debug("realign::align_from_variants::value %s at %s" % (value, pos))
            if value in snp:
                snps.append([pos, sequence[p], mature[p]])
    logger.debug("realign::align_from_variants::snps %s" % snps)
    return snps


def variant_to_5p(hairpin, pos, variant):
    """
    From a sequence and a start position get the nts
       +/- indicated by iso_5p. Pos option is 0-base-index

    Args:
       *hairpin(str)*: long sequence:
            >>> AAATTTT

       *position(int)*: >>> 3

       *variant(int)*: number of nts involved in the variant:
            >>> -1
    Returns:
       *(str)*: nucleotide involved in the variant:
            >>> T
    """
    pos = pos[0]
    iso_t5 = [v for v in variant.split(",") if v.startswith("iso_5p")]
    if iso_t5:
        t5 = int(iso_t5[0].split(":")[-1][-1])
        direction_t5 = int(iso_t5[0].split(":")[-1]) * -1
        if direction_t5 > 0:
            return hairpin[pos - t5:pos]
        elif direction_t5 < 0:
            return hairpin[pos:pos + t5].lower()
    return "0"


def variant_to_3p(hairpin, pos, variant):
    """
    From a sequence and a start position get the nts
       +/- indicated by iso_3p. Pos option is 0-base-index

    Args:
       *hairpin(str)*: long sequence:
            >>> AAATTTT

       *position(int)*: >>> 3

       *variant(int)*: number of nts involved in the variant:
            >>> -1
    Returns:
       *(str)*: nucleotide involved in the variant:
            >>> A
    """
    pos = pos[1]
    iso_t3 = [v for v in variant.split(",") if v.startswith("iso_3p")]
    if iso_t3:
        t3 = int(iso_t3[0].split(":")[-1][-1])
        direction_t3 = int(iso_t3[0].split(":")[-1])
        if direction_t3 > 0:
            return hairpin[pos + 1:pos + t3 + 1]
        elif direction_t3 < 0:
            return hairpin[pos - t3 + 1:pos + 1].lower()
    return "0"


def variant_to_add(read, variant):
    """
    From a sequence and a start position get the nts
       +/- indicated by iso_3p. Pos option is 0-base-index

    Args:
       *hairpin(str)*: long sequence:
            >>> AAATTTT

       *position(int)*: >>> 3

       *variant(int)*: number of nts involved in the variant:
            >>> 2
    Returns:
       *(str)*: nucleotide involved in the variant:
            >>> TT
    """
    add = [v for v in variant.split(",") if v.startswith("iso_add")]
    if add:
        add = int(add[0].split(":")[-1][-1]) * -1
        return read[add:]
    return "0"
