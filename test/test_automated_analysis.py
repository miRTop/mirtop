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


@contextlib.contextmanager
def make_workdir():
    remove_old_dir = True
    dirname = os.path.join(os.path.dirname(__file__), "test_automated_output")
    if remove_old_dir:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.makedirs(dirname)
    orig_dir = os.getcwd()
    try:
        os.chdir(dirname)
        yield dirname
    finally:
        os.chdir(orig_dir)

def expected_failure(test):
    """Small decorator to mark tests as expected failure.
    Useful for tests that are work-in-progress.
    """
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise SkipTest
        else:
            raise AssertionError('Failure expected')
    return inner

class AutomatedAnalysisTest(unittest.TestCase):
    """Setup a full automated analysis and run the pipeline.
    """
    def setUp(self):
        self.data_dir = os.path.join(os.path.dirname(__file__), "data", "automated")

    def _install_test_files(self, data_dir):
        """Download required sequence and reference files.
        """
        #       self._download_to_dir(url, dirname)

    def _download_to_dir(self, url, dirname):
        print dirname
        cl = ["wget", url]
        subprocess.check_call(cl)
        cl = ["tar", "-xzvpf", os.path.basename(url)]
        subprocess.check_call(cl)
        shutil.move(os.path.basename(dirname), dirname)
        os.remove(os.path.basename(url))

    @attr(simulate=True)
    def test_simulate(self):
        """Check simulated data"""
        mirna = "TGAGGTAGTAGGTTGTATAGTT"
        correct = 0
        n = 0
        with open("data/examples/simulation/res/reads.mirna") as inh:
            header = inh.readline()
            for line in inh:
                cols = line.strip().split()
                mut, add, t5, t3 = cols[6:10]
                seq = cols[0]
                if mut!="0":
                    pos = int(mut[:-2])
                    nt1 = mut[-2]
                    nt2 = mut[-1]
                    seql = list(seq)
                    seql[pos] = nt2
                    seq = "".join(seql)
                if t5!="0" and t5.islower():
                    seq = "%s%s" % (t5.upper(), seq)
                elif t5!="0" and t5.isupper():
                    seq = seq[len(t5):]
                if add!="0":
                    seq = seq[:-len(add)]
                if t3!="0" and t3.islower():
                    seq = "%s%s" % (seq, t3.upper())
                elif t3!="0" and t3.isupper():
                    seq = seq[:-len(t3)]
                if seq == mirna:
                    correct += 1
                else:
                    print "\nerror:\n%s\n%s" % (seq, mirna)
                n += 1
        print "rate %s/%s" % (correct, n)

    @attr(complete=True)
    @attr(annotate=True)
    def test_srnaseq_annotation_bam(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs_fasta",
                      "../../data/examples/annotate/sim_isomir.sam"]
            print ""
            print " ".join(clcode)
            subprocess.check_call(clcode)

    @attr(complete=True)
    @attr(annotate=True)
    def test_srnaseq_annotation_seqbuster(self):
        """Run miraligner analysis
        """
        with make_workdir():
            clcode = ["mirtop",
                      "gff",
                      "--format", "seqbuster",
                      "--sps", "hsa",
                      "--hairpin", "../../data/examples/annotate/hairpin.fa",
                      "--gtf", "../../data/examples/annotate/hsa.gff3",
                      "-o", "test_out_mirs_fasta",
                      "../../data/examples/seqbuster/reads.mirna"]
            print ""
            print " ".join(clcode)
            subprocess.check_call(clcode)
