from __future__ import print_function

import argparse
import os
import random
from collections import defaultdict

import pysam

import mirtop.libs.logger as mylog
import mirtop.libs.do as runner

parser = argparse.ArgumentParser()
parser.add_argument("--fa",
                    help="File with mature sequences.", required=True)
parser.add_argument("-o", "--out", default="spikeins.fa",
                    help="Name used for output files.")
parser.add_argument("--seed", help="set up seed for reproducibility.",
                    default=42)
parser.add_argument("--max_size", help="maximum size allowed in the final output.",
                    default=25)

args = parser.parse_args()
random.seed(args.seed)

def _sam_to_bam(bam_fn):
    bam_out = "%s.bam" % os.path.splitext(bam_fn)[0]
    cmd = "samtools view -Sbh {bam_fn} -o {bam_out}"
    runner.run(cmd.format(**locals()))
    return bam_fn


def _bam_sort(bam_fn):
    bam_sort_by_n = os.path.splitext(bam_fn)[0] + "_sort.bam"
    runner.run(("samtools sort -n -o {bam_sort_by_n} {bam_fn}").format(
               **locals()))
    return bam_sort_by_n

def _read_fasta(fa):
    source = dict()
    with open(fa) as inh:
        for line in inh:
            if line.startswith(">"):
                name = line.strip().split()[0].replace(">", "")
            else:
                source.update({name: line.strip()})
    return source


def _write_fasta(sequences, filename, max=25):
    with open(filename, 'w') as outh:
        for name in sequences:
            if sequences[name]:
                if len(sequences[name]) < max:
                    print(">%s\n%s" % (name, sequences[name]), file=outh)
    return filename


def _parse_hits(sam, source):
    uniques = defaultdict(list)
 #   bam_fn = _sam_to_bam(sam)
 #   bam_fn = _bam_sort(bam_fn)
 # read sequences and score hits (ignore same sequence)
    handle = pysam.Samfile(sam, "rb")
    for line in handle:
        reference = handle.getrname(line.reference_id)
        name = line.query_name
        # sequence = line.query_sequence if not line.is_reverse else reverse_complement(line.query_sequence)
        if reference == name:
            continue
        # print([reference, name, line.get_tag("NM")])
        distance = line.get_tag("NM")
        uniques[name].append(distance)
        uniques[reference].append(distance)
    # read parsed data and keep the ones with score > 10 edit distance
    for name in uniques:
        if min(uniques[name]) < 5:
            if name in source:
                source[name] = None
    return source

# Map all vs all with razers3
source = _read_fasta(args.fa)
sam = os.path.join(os.path.dirname(args.out), "modified.bam")
runner.run(("razers3 -dr 5 -i 75 -rr 80 -f -so 1 -o {output} {target} {query}").format(output=sam, target=args.fa, query=args.fa))
uniques = _parse_hits(sam, source)

# Write uniques to fasta
_write_fasta(uniques, args.out, args.max_size)
