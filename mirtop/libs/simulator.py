"""simulate cluster over the genome"""
from __future__ import print_function
import random
from mirtop.libs.read import get_fasta


def simulate(args):
    """Main function that manage simulation of small RNAs"""
    reads = dict()
    if args.fasta:
        name = None
        seq = ""
        with open(args.fasta) as in_handle:
            for line in in_handle:
                if line.startswith(">"):
                    if name:
                        reads.update(_generate_reads(seq, name))
                    seq = ""
                    name = line[1:-1]
                else:
                    seq += line.strip()

        reads.update(_generate_reads(seq, name))
    if reads:
         _write_reads(reads, args.out)


def _generate_reads(seq, name):
    """Main function that create reads from precursors"""
    reads = dict()
    if len(seq) < 130 and len(seq) > 70:
        reads.update(_mature(seq[:40], 0, name))
        reads.update(_mature(seq[-40:], len(seq) - 40, name))
        reads.update(_noise(seq, name))
        reads.update(_noise(seq, name, 25))
    return reads


def _mature(subseq, absolute, c,  size=33, total=5000):
    """Create mature sequences around start/end"""
    reads = dict()
    probs = [0.1, 0.2, 0.4, 0.2, 0.1]
    end = 5 + size
    error = [-2, -1, 0, 1, 2]
    for error5 in error:
        for error3 in error:
            s = 5 - error5
            e = end - error3
            seen = subseq[s:e]
            counts = int(probs[error5 + 2] * probs[error3 + 2] * total) + 1
            name = "seq_%s_%s_%s_x%s" % (c, s + absolute, e + absolute, counts)
            reads[name] = (seen, counts)
    return reads


def _noise(seq, c, size=33, total=1000):
    """Create mature sequences around start/end"""
    reads = dict()
    seen = 0
    while seen < total:
        s = random.randint(0, len(seq) - size)
        e = s + size + random.randint(-5,5)
        p = random.uniform(0, 0.1)
        counts = int(p * total) + 1
        seen += counts
        name = "seq_%s_%s_%s_x%s" % (c, s, e, counts)
        reads[name] = (seq[s:e], counts)
    return reads


def _write_reads(reads, prefix):
    """
    Write fasta file, ma file and real position
    """
    out_ma = prefix + ".ma"
    out_fasta = prefix + ".fasta"
    out_real = prefix + ".txt"
    with open(out_ma, 'w') as ma_handle:
        print("id\tseq\tsample", file=ma_handle, end="")
        with open(out_fasta, 'w') as fa_handle:
            with open(out_real, 'w') as read_handle:
                for idx, r in enumerate(reads):
                    info = r.split("_")
                    print("seq_%s\t%s\t%s" % (idx, reads[r][0], reads[r][1]), file=ma_handle, end="")
                    print(">seq_%s\n%s" % (idx, reads[r][0]), file=fa_handle, end="")
                    print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (idx, r, reads[r][0], reads[r][1], info[1], info[2], info[3]), file=read_handle, end="")


def _get_precursor(bed_file, reference, out_fa):
    """
    get sequence precursor from position
    """
    get_fasta(bed_file, reference, out_fa)
    return 0


def _get_spot(precursor):
    """
    get spot that will be enriched
    """
    return 0


def _get_type(pob):
    """
    randomly decide if is small rna or degradation
    """
    return 0


def _random_sequences(precursor, start= None, end=None):
    """
    randomly get sequences around some nucleotides.
    It could be enriched in some positions
    """
    return 0
