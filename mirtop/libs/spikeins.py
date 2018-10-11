"""This code helps to work with spikeins"""
import os

from collections import defaultdict
from mirtop.libs import fastq
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def convert(args):
    out_dir = args.out if args.out else os.path.dirname(args.file)
    gff_file = os.path.join(out_dir, "spikeins_pre.gff")
    fasta_file = os.path.join(out_dir, "spikeins_pre.fasta")
    spikes = read_spikeins(args.file)
    write_gff(spikes, gff_file)
    logger.info("GFF file is at %s" % gff_file)
    write_precursors(spikes, fasta_file)
    logger.info("FASTA file is at %s" % fasta_file)


def read_spikeins(in_file):
    """Read FASTA file containing small spikeins

        Args:
            *in_file(str)*: file name with sequences.

        Returns:
            *spikeins(dict)*: dictionary with keys being names
             and values being: given spike-in, synthetic precursor and
             position on the synthetic precursor.
    """
    spikeins = defaultdict(dict)
    with fastq.open_fastq(in_file) as inh:
        for line in inh:
            if line.startswith(">"):
                spike = line.strip()[1:].split()[0]
            else:
                precursor = "GGGGG" + line.strip() + "GGGGG"
                pos = [5, 4 + len(line.strip())]
                logger.debug("SPIKE::NAME::%s -> %s" % (spike, line.strip()))
                logger.debug("SPIKE::GENERATED::precursor::%s" % precursor)
                logger.debug("SPIKE::GENERATED::position::%s" % pos)
                spikeins[spike] = {'mature': line.strip(),
                                   'precursor': precursor,
                                   'position': pos}
    return spikeins


def write_precursors(spikeins, out_file):
    """
    Write FASTA file for precursors to be used to annotate isomiRs

        Args:
            *spikeins(dict)*: output from read_spikeins
            *out_file*: file name to write the FASTA information
    """
    with open(out_file, 'w') as outh:
        for spike in spikeins:
            outh.write(">pre-{0}\n{1}\n".format(spike,
                                                 spikeins[spike]['precursor']))


def write_gff(spikeins, out_file):
    """
    Write GTF file for precursors to be used to annotate isomiRs

        Args:
            *spikeins(dict)*: output from read_spikeins
            *out_file*: file name to write the GFF information
    """
    with open(out_file, 'w') as outh:
        outh.write("#microRNAs:               spikeins_v1\n")
        for spike in spikeins:
            outh.write("\t".join(["chr%s" % spike, ".",
                                  "miRNA_primary_transcript",
                                  "1",
                                  str(len(spikeins[spike]['precursor'])),
                                  ".", "+", ".",
                                  "Name=pre-%s;ID=pre-%s" % (spike, spike),
                                  "\n"]))
            outh.write("\t".join(["chr%s" % spike, ".", "miRNA",
                                  str(spikeins[spike]['position'][0] + 1),
                                  str(spikeins[spike]['position'][1] + 1),
                                  ".", "+", ".",
                                  "Derives_from=pre-%s;"
                                  "Name=%s;ID=%s;" % (spike, spike, spike),
                                  "\n"]))
