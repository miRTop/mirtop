import os
import argparse
import sqlite3
import re

conn = sqlite3.connect('mirtop.db')
parser = argparse.ArgumentParser(description='Create and query the contents of GFF3 through sqlite')
subparsers = parser.add_subparsers(help='commands')

# Options for creating database
create_parser = subparsers.add_parser('create', help="Create new database")
create_parser.add_argument('--gff', metavar='', action='store', required=True, help="/path/to/GFF/file/file.gff")
create_parser.add_argument('--db', metavar='', action='store', help='Database name to create. (default: mirtop.db)')

# Query from database
query_parser = subparsers.add_parser('query', help='Query from database')
query_parser.add_argument('--db', metavar='', action='store', required=True, help='Database name to query from ...')
query_parser.add_argument('--ex', metavar='', action='store', help='Expression of SQL query')

args = parser.parse_args()
user_options = vars(args)

'''
header:

    (R) small RNA GFF version ## VERSION: 1.2 
    (R) database: ##source-ontology using FAIRSharing.org: miRBase: (FAIRsharing) doi:10.25504/fairsharing.hmgte8 
    mirGeneDB: http://mirgenedb.org mirCarta: https://mircarta.cs.uni-saarland.de/ Custom database: please,
     provide a link to an archive release if this is the case 
    (R) tools used starting with the label ## TOOLS: and followed by tools used to call isomiRs separated by 
    comma (,). 
    (O) commands used to generate the file. At least information about adapter removal, filtering, 
    aligner, mirna tool. All of them starting like: ## CMD: . Can be multiple lines starting with this tag. 
    (O) genome/database version used (maybe try to get from BAM file if GFF3 generated from it): ## REFERENCE: (R) sample 
    names used in attribute:Expression: ## COLDATA: separated by comma: ,. 
    (O) Filter tags meaning: See Filter attribute below. Different filter tags should be separated by , character. 
    Example: ## FILTER: and example would be ## FILTER: PASS(is ok), REJECT(false positive), REJECT 
    lowcount(rejected due to low count in data). 

'''


def insert_sql():
    with open(args.gff, 'r') as f:
        for text in f:
            # HEADER INFORMATION
            if re.search("## VERSION", text):
                version = (text.strip().split(' ')[-1])
            if re.search("## source-ontology", text):
                source = (text.strip().split(' ')[-1])
            if re.search("## COLDATA", text):
                data_sets = (text.strip().split(' ')[-1])  # Might contain more than one data set, need to edit in furute
            if re.search("## TOOLS", text):
                tools = (text.strip().split(' ')[-1])
            if re.search("## CMD", text):
                commands_exec = (text.strip().split(' ')[-1])
            if re.search("## FILTER", text):
                filter_tags = (text.strip().split(' ')[-1])
            if re.search("## REFERENCE", text):
                references = (text.strip().split(' ')[-1])
            if not text.startswith('#'):
                text = f.readline().strip().split('\t')
                print(text, end='')
                break
        pass
    print()
    print("I am in create db")
    print(args.gff)


def query_sql():
    print("I am in query")
    print(args.db)
    pass


def main():
    if 'gff' in user_options.keys():
        insert_sql()
    elif 'db' in user_options.keys():
        query_sql()
    else:
        print("Version: v1.0.0")
        print()
        print("usage: gff2sql [-h] {create,query} ... ")


if __name__ == "__main__":
    main()
