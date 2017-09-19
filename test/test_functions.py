"""This directory is setup with configurations to run the main functional test.

Inspired in bcbio-nextgen code
"""
import os
import subprocess
import unittest
import shutil
import contextlib
import collections
import functools

from nose import SkipTest
from nose.plugins.attrib import attr
import yaml


class FunctionsTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    @attr(read=True)
    def test_read(self):
        from mirtop.mirna import map, fasta
        from mirtop.libs import logger
        logger.initialize_logger("test_read_files", True, True)
        map_mir = map.read_gtf_to_precursor("data/examples/annotate/hsa.gff3")
        if map_mir["hsa-let-7a-1"]["hsa-let-7a-5p"][0] != 4:
            raise ValueError("GFF is not loaded correctly.")
        fasta_precursor = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        # read data/aligments/let7-perfect.bam
        return True

    @attr(cigar=True)
    def test_cigar(self):
        """testing cigar correction function"""
        cigar = [[0, 14], [1, 1], [0, 5]]
        from mirtop.mirna.realign import cigar_correction
        fixed = cigar_correction(cigar, "AAAAGCTGGGTTGAGGAGGA", "AAAAGCTGGGTTGAGAGGA")
        print "\n testing cigar correction"
        print fixed[0]
        print fixed[1]

    @attr(alignment=True)
    def test_alignment(self):
        """testing alignments function"""
        from mirtop.libs import logger
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        from mirtop.mirna import fasta
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        matures = {}
        # matures = mirtop.mirna.read_mature("data/examples/annotate/mirnas.gff", "hsa")
        def annotate(fn, precursors, matures):
            from mirtop.bam import bam
            reads = bam.read_bam(fn, precursors)
            # ann = mirtop.bam.bam.annotate(reads, matures, precursors)
            return True
        print "\nlast1D\n"
        annotate("data/aligments/let7-last1D.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   AGAGGTAGTAGGTTGTA
        print "\n1D\n"
        annotate("data/aligments/let7-1D.sam", precursors, matures)
        #mirna TGAGGTAG-TAGGTTGTATAGTT
        #seq   TGAGGTAGGTAGGTTGTATAGTTA
        print "\nlast7M1I\n"
        annotate("data/aligments/let7-last7M1I.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTA-AGT
        print "\nmiddle1D\n"
        annotate("data/aligments/let7-middle1D.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTATAGTT
        print "\nperfect\n"
        annotate("data/aligments/let7-perfect.sam", precursors, matures)
        #mirna TGAGGTAGTAGGTTGTATAGTT
        #seq   TGAGGTAGTAGGTTGTATAG (3tt 3TT)
        print "\ntriming\n"
        annotate ("data/aligments/let7-triming.sam", precursors, matures)

    @attr(gff=True)
    def test_gff(self):
        """testing GFF function"""
        from mirtop.libs import logger
        from mirtop.mirna import map, fasta
        logger.initialize_logger("test", True, True)
        logger = logger.getLogger(__name__)
        precursors = fasta.read_precursor("data/examples/annotate/hairpin.fa", "hsa")
        # depend on https://github.com/miRTop/mirtop/issues/6
        matures = map.read_gtf_to_precursor("data/examples/annotate/hsa.gff3")
        # matures = mirtop.mirna.read_mature("data/examples/annotate/mirnas.gff", "hsa")
        from mirtop.bam import bam
        reads = bam.read_bam("data/aligments/let7-perfect.sam", precursors)
        ann = bam.annotate(reads, matures, precursors)
        # gff = mirtop.gff.body.create(ann)
        return True


