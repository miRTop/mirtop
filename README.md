mirtop
---------

[![](https://travis-ci.org/miRTop/mirtop.svg?branch=master)](https://travis-ci.org/miRTop/mirtop#)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

Command line tool to annotate with a standard naming miRNAs e isomiRs.

See more at: [isomiRs naming discussion](https://github.com/miRTop/incubator/blob/master/isomirs/isomir_naming.md)

Cite
---------

http://mirtop.github.io

About
-----

Join the team: https://orgmanager.miguelpiedrafita.com/join/15463928

Read more: http://mirtop.github.io

Installation
------------
Thes best solution is to install conda to get an independent enviroment.

`wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh`

`bash Miniconda-latest-Linux-x86_64.sh -b -p ~/mirtop_env`

`export PATH=$PATH:~/mirtop_env`

`git clone http://github.com/miRTop/mirtop`

`cd mirtop`

`python setup.py develop`

Quick start
-----------

`cd data`
`mirtop gff --sps hsa --hairpin examples/annotate/hairpin.fa --mirna examples/annotate/hsa.gff3 -o test_out_mirs_fasta examples/annotate/sim_isomir.sam

Output
------

There will be a `*.mirna` file that is a tabular file with this format:

Naming of isomiRs follow these rules: https://github.com/miRTop/incubator/blob/master/format/definition.md

Contributors
------------

* [Lorena Pantano](https://github.com/lpantano) (Bioinformatic Core, Harvard Chan School, Boston, USA)

