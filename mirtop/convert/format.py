# import os
import argparse
import sqlite3
import re
from datetime import datetime
import time
import os.path as op

now = datetime.now()
# dd/mm/YY H:M:S
d2 = now.strftime("%B %d, %Y %H:%M:%S")


# print("date and time =", d2)

def insert_sql(args):
    if args.db:
        out_file = op.join(args.out, (args.db))
        conn = sqlite3.connect(out_file)
        c = conn.cursor()
    else:
        out_file = op.join(args.out, "mirtop.db")
        conn = sqlite3.connect(out_file)
        c = conn.cursor()

    with open(args.gff, 'r') as f:
        version = source = data_sets = tools = commands_exec = filter_tags = citation = num_records = ""
        cnt = 0
        c.execute('''CREATE TABLE IF NOT EXISTS data_sets(seqID text, source_file text, type text, start real, 
        end real, score text, strand text, phase text, UID text, Read text, Name text, Parent text, Variant text, 
        iso_5p real, iso_3p real, iso_add3p real, iso_snp	real, iso_5p_nt	real, iso_3p_nt	real, iso_add3p_nt real, 
        iso_snp_nt real, source text, cigar text, hits real, alias text, genomic_pos text, expression text, filter, 
        seed_fam text)''')
        conn.commit()
        for text in f:
            # HEADER INFORMATION
            if re.search("^## VERSION", text):  # (R)
                version = (text.strip().split(' ')[-1])
            elif re.search("^## source-ontology", text):  # (R)
                source = (text.strip().split(' ')[-1])
            elif re.search("^## COLDATA", text):  # (R)
                data_sets = (
                    text.strip().split(' ')[-1])  # Might contain more than one data set, need to edit in future
            elif re.search("^## TOOLS", text):  # (R)
                tools = (text.strip().split(' ')[-1])
            elif re.search("^## CMD", text):  # (O)
                commands_exec = (text.strip().split(' ')[-1])
            elif re.search("^## FILTER", text):  # (O)
                filter_tags = (text.strip().split(' ')[-1])
            elif re.search("^## REFERENCE", text):  # (O)
                citation = (text.strip().split(' ')[-1])
            # BODY - INFORMATION
            elif not re.search("^#", text):
                cnt += 1
                # print(str(cnt) + "    " + text)
                lines = text.strip().split('\t')
                if '=' in lines[-1]:
                    lines[-1].replace('=', ' ')

                info = lines[-1].split('; ')
                info_dict = variant_dict = dict()
                for elements in info:
                    (k, v) = elements.split(' ')
                    info_dict.update([(k, v)])
                    # iso_snp, iso_add: +4, iso_5p: -1;
                    if 'Variant' in k and ":" in v:
                        value_list = v.split(',')
                        for iso_vars in value_list:
                            if ":" in iso_vars:
                                (sub_k, sub_v) = iso_vars.split(':')
                                # print(sub_k + " ------> Here I am " + sub_v)
                                # time.sleep(3)
                                variant_dict.update([(str(sub_k), str(sub_v))])

                c.execute(
                    "INSERT INTO data_sets(seqID, source_file, type, start, end, score, strand, phase, UID, Read, "
                    "Name, Parent, Variant, iso_5p, iso_3p, iso_add3p, iso_snp, iso_5p_nt, iso_3p_nt, iso_add3p_nt, "
                    "iso_snp_nt, source, cigar, hits, alias, genomic_pos, expression, filter, seed_fam) VALUES (?, "
                    "?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (lines[0], lines[1], lines[2], lines[3], lines[4], lines[5], lines[6], lines[7],
                     info_dict.get('UID'), info_dict.get('Read'), info_dict.get('Name'), info_dict.get('Parent'),
                     info_dict.get('Variant'), str(info_dict.setdefault('iso_5p', None)),
                     str(info_dict.setdefault('iso_3p', None)), str(info_dict.setdefault('iso_add3p', None)),
                     str(info_dict.setdefault('iso_snp', None)), str(info_dict.setdefault('iso_5p_nt', None)),
                     str(info_dict.setdefault('iso_3p_nt', None)), str(info_dict.setdefault('iso_add3p_nt', None)),
                     str(info_dict.setdefault('iso_snp_nt', None)), source, info_dict.setdefault('Cigar', None),
                     info_dict.setdefault('Hits', None), info_dict.setdefault('Alias', None),
                     info_dict.setdefault('Genomic', None), info_dict.setdefault('Expression', 0),
                     info_dict.setdefault('Filter', None), info_dict.setdefault('Seed_fam', None)))
                conn.commit()
                # break

        c.execute('''CREATE TABLE IF NOT EXISTS summary(version text, source text, data_sets text, tools text,
         commands_exec text, filter_tags text, citation text, records real, date_stamp text)''')
        c.execute("INSERT INTO summary(version, source, data_sets, tools, commands_exec, filter_tags, citation, "
                  "records, date_stamp) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (version, source, data_sets, tools, commands_exec, filter_tags, citation, cnt, d2))

        # info_dict.setdefault('Sex', None)

        conn.commit()

    print()
    print("I am in create db")
    print(args.gff)


def query_sql(args):
    print("I am in query")
    print(args.db)
    pass


def sql(args):
    #insert_sql(args)
    user_options = vars(args)
    if args.gff: 
    #if 'gff' in user_options.keys():
        insert_sql(args)
    elif args.db:
        query_sql(args)
    else:
        print("Usage: mirtop convert --format sql --gff <input.gff> --db <input_dabasename>")
        print()
    print(user_options)
    print("Passed argument was SQL and requirement fulfilled")
    pass


def format_option(args):
    print(args)
    if args.format == 'sql':
        sql(args)
    else:
        pass

#    parser = argparse.ArgumentParser(description='Create and query the contents of GFF3 through sqlite')
#    subparsers = parser.add_subparsers(help='commands')
#
#    # Options for creating database
#    create_parser = subparsers.add_parser('create', help="Create new database")
#    create_parser.add_argument('--gff', metavar='', action='store', required=True, help="/path/to/GFF/file/file.gff")
#    create_parser.add_argument('--db', metavar='', action='store', help='Database name to create. (default: mirtop.db)')
#
#    # Query from database
#    query_parser = subparsers.add_parser('query', help='Query from database')
#    query_parser.add_argument('--db', metavar='', action='store', required=True, help='Database name to query from ...')
#    query_parser.add_argument('--ex', metavar='', action='store', help='Expression of SQL query')
#    query_parser.add_argument('--overview', metavar='', action='store', help='Expression of SQL query')
#
#    args = parser.parse_args()
#    user_options = vars(args)
#    else:
#        print("Version: v1.0.0")
#        print()
#        print("usage: gff2sql [-h] {create,query} ... ")
#

