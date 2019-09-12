# import os
import argparse
import sqlite3
import re


def insert_sql(args):
    with open(args.gff, 'r') as f:
        for text in f:
            # HEADER INFORMATION
            if re.search("## VERSION", text):
                version = (text.strip().split(' ')[-1])
            if re.search("## source-ontology", text):
                source = (text.strip().split(' ')[-1])
            if re.search("## COLDATA", text):
                data_sets = (
                text.strip().split(' ')[-1])  # Might contain more than one data set, need to edit in furute
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


def query_sql(args):
    print("I am in query")
    print(args.db)
    pass


def gff2sql():
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
    if 'gff' in user_options.keys():
        insert_sql(args)
    elif 'db' in user_options.keys():
        query_sql(args)
    else:
        print("Version: v1.0.0")
        print()
        print("usage: gff2sql [-h] {create,query} ... ")


gff2sql()
