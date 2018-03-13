mirtop
---------

[![](https://travis-ci.org/miRTop/mirtop.svg?branch=master)](https://travis-ci.org/miRTop/mirtop#)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

Command line tool to annotate with a standard naming miRNAs e isomiRs.

This tool adapt the miRNA GFF3 format agreed on here: https://github.com/miRTop/incubator/blob/master/format/definition.md

Chat
----

[Ask question, ideas](https://gitter.im/mirtop/Lobby#)
[Contributors to code](https://gitter.im/mirtop/devel)

Cite
---------

http://mirtop.github.io

About
-----

Join the team: https://orgmanager.miguelpiedrafita.com/join/15463928

Read more: http://mirtop.github.io

Installation
------------

### Bioconda

`conda install mirtop -c bioconda`

### PIP

`pip install mirtop`

### develop version

Thes best solution is to install conda to get an independent enviroment.

```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

bash Miniconda-latest-Linux-x86_64.sh -b -p ~/mirtop_env

export PATH=$PATH:~/mirtop_env

conda install -c bioconda pysam pybedtools pandas biopython samtools

git clone http://github.com/miRTop/mirtop

cd mirtop

python setup.py develop
```

Quick start
-----------

Read complete commands at: https://mirtop.readthedocs.org

```
git clone mirtop
cd mirtop/data/example/annotate
mirtop gff --sps hsa --hairpin hairpin.fa --gtf hsa.gff3 -o test_out sim_isomir.bam
```

Output
------

The `mirtop gff` generates the GFF3 adapter format to capture miRNA variations. The output is explained [here](https://github.com/miRTop/incubator/blob/master/format/definition.md).

Contributors
------------

* [Lorena Pantano](https://github.com/lpantano) (Bioinformatic Core, Harvard Chan School, Boston, USA)
* [Shruthi Bhat Bandyadka](https://github.com/sbb25) (Partners Personalized Medicine, Cambridge MA, USA)
