# Installation

## bioconda

`conda install mirtop -c bioconda`

## pypi

`pip install mirtop`

## update to develop version from pip

```
pip install --upgrade --no-deps git+https://github.com/miRTop/mirtop.git#egg=mirtop
```

## install develop version

Thes best solution is to install conda to get an independent enviroment.

```
wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh

bash Miniconda-latest-Linux-x86_64.sh -b -p ~/mirtop_env

export PATH=$PATH:~/mirtop_env

conda install -c bioconda bioconda bedtools samtools pip nose pysam pandas dateutil pyyaml pybedtools biopython setuptools

git clone http://github.com/miRTop/mirtop
cd mirtop
git fetch origin dev
git checkout dev

python setup.py develop
```