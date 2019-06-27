from __future__ import print_function

import argparse
import os
import random
from collections import defaultdict

import pysam

import mirtop.libs.logger as mylog
import mirtop.libs.do as runner


def _read_fasta(fa, size):
    source = dict()
    with open(fa) as inh:
        for line in inh:
            if line.startswith(">"):
                name = line.strip().split()[0].replace(">", "")
            else:
                if len(line.strip()) >= size:
                    source.update({name: line.strip()[0:size]})
    return source


def _update_ends(source):
    nts = ["A", "T", "C", "G"]
    start_idx = 0
    end_idx = 0
    for name in source:
        source[name] = nts[start_idx] + source[name] + nts[end_idx]
        if end_idx == 3 and start_idx == 3:
            end_idx = -1
            start_idx = 0
        if end_idx == 3:
            start_idx += 1
            end_idx = 0
        end_idx += 1
    return source


def _write_fasta(sequences, filename):
    with open(filename, 'w') as outh:
        for name in sequences:
            if sequences[name]:
                print(">%s\n%s" % (name, sequences[name]), file=outh)
    return filename


def _parse_hits(sam, source):
    uniques = defaultdict(list)
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
        if min(uniques[name]) < 4:
            if name in source:
                source[name] = None
    return source


parser = argparse.ArgumentParser()
parser.add_argument("--fa",
                    help="File with mature sequences.", required=True)
parser.add_argument("-s", "--size", default=22,
                    help="Size of spike-ins to generate.")
parser.add_argument("-n", "--number", default=16,
                    help="Number of spike-ins to generate.")
parser.add_argument("-o", "--out", default="spikeins.fa",
                    help="Name used for output files.")
parser.add_argument("--seed", help="set up seed for reproducibility.",
                    default=42)
parser.add_argument("--universe", help="Set up universe sequences to avoid duplication.",
                    default=None)

args = parser.parse_args()
random.seed(args.seed)

mylog.initialize_logger(os.path.dirname(os.path.abspath(args.out)))
logger = mylog.getLogger(__name__)

# Read file to get all sequences longer than size - 2
size = args.size - 2
source = _read_fasta(args.fa, size)
logger.info("%s was read: %s sequences were loaded" % (args.fa, len(source)))
source = _update_ends(source)
logger.info("source updated with extended nts: %s" % source)

# Map all vs all with razers3
modified = _write_fasta(source, os.path.join(os.path.dirname(args.out), "modified.fa"))
sam = os.path.join(os.path.dirname(args.out), "modified.bam")
runner.run(("razers3 -i 75 -rr 80 -f -so 1 -o {output} {target} {query}").format(output=sam, target=modified, query=modified))
uniques = _parse_hits(sam, source)
print(uniques)
if args.universe:
    sam = os.path.join(os.path.dirname(args.out), "modified_vs_universe.sam")
    runner.run(("razers3 -i 75 -rr 80 -f -o {output} {target} {query}").format(output=sam, target=args.universe, query=modified))
    uniques = _parse_hits(sam, uniques)
print(uniques)

# Write uniques to fasta
_write_fasta(uniques, args.out)
