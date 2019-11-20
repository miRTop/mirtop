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


def create_table(conn, sample_names):
    c = conn.cursor()
    data_columns = ['seqID text', 'source_file text', 'type text', 'start real',
    'end real', 'score text', 'strand text', 'phase text', 'UID text', 'Read text', 'Name text', 'Parent text', 'Variant text',
    'iso_5p real', 'iso_3p real', 'iso_add3p real', 'iso_snp real', 'iso_5p_nt	real', 'iso_3p_nt real', 'iso_add3p_nt real',
    'iso_snp_nt real', 'source text', 'cigar text', 'hits real', 'alias text', 'genomic_pos text', 'filter text',
    'seed_fam text']
    complete_headers = data_columns + sample_names
    q = "CREATE TABLE IF NOT EXISTS data_sets(%s)" % ", ".join(complete_headers)
    c.execute(q)
    conn.commit()


def gff_insert_values(conn, complete_list):
    try:
        conn.execute('INSERT INTO data_sets VALUES(' + ','.join("?" * len(complete_list)) + ')', complete_list)
        conn.commit()
    except sqlite3.OperationalError as e:
        print()
        print("ERROR:")
        print("sqlite3.OperationalError: {0}".format(e))
        print("Help: Make sure to delete any existing database with tables of different schema")
        exit()


# print("date and time =", d2)

def insert_sql(args):
    if args.db:
        out_file = op.join(args.out, args.db)
        conn = sqlite3.connect(out_file)
        c = conn.cursor()
    else:
        out_file = op.join(args.out, "mirtop.db")
        conn = sqlite3.connect(out_file)
        c = conn.cursor()

    with open(args.gff, 'r') as f:
        version = source = data_sets = tools = commands_exec = filter_tags = citation = num_records = ""
        cnt = 0
        for text in f:
            # HEADER INFORMATION
            if re.search("^## .* VERSION", text):  # (R)
                version = (text.strip().split(' ')[-1])
            elif re.search("^## source-ontology", text):  # (R)
                source = (text.strip().split(' ')[-1])
            elif re.search("^## COLDATA", text):  # (R)
                data_sets = (
                    text.strip().split(' ')[-1])  # Might contain more than one data set
                sample_names = data_sets.split(',')
                sample_names = [w.replace('-', '_') for w in sample_names]
                string_text = "text"
                output_sample_names = ["{} {}".format(i, string_text) for i in sample_names]
                create_table(conn, output_sample_names)
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
                lines = text.strip().split('\t')
                if '=' in lines[-1]:
                    lines_info_array = lines[-1].replace("=", " ")
                else:
                    lines_info_array = lines[-1]

                info = lines_info_array.split('; ')
                info_dict = variant_dict = dict()
                for elements in info:
                    (k, v) = elements.split(' ')
                    info_dict.update([(k, v)])
                    if 'Variant' in k and ":" in v:
                        value_list = v.split(',')
                        for iso_vars in value_list:
                            if ":" in iso_vars:
                                (sub_k, sub_v) = iso_vars.split(':')
                                variant_dict.update([(str(sub_k), str(sub_v))])

                prefix_list = [lines[0], lines[1], lines[2], lines[3], lines[4], lines[5], lines[6], lines[7],
                               info_dict.get('UID'), info_dict.get('Read'), info_dict.get('Name'),
                               info_dict.get('Parent'),
                               info_dict.get('Variant'), str(info_dict.setdefault('iso_5p', None)),
                               str(info_dict.setdefault('iso_3p', None)), str(info_dict.setdefault('iso_add3p', None)),
                               str(info_dict.setdefault('iso_snp', None)), str(info_dict.setdefault('iso_5p_nt', None)),
                               str(info_dict.setdefault('iso_3p_nt', None)),
                               str(info_dict.setdefault('iso_add3p_nt', None)),
                               str(info_dict.setdefault('iso_snp_nt', None)), source,
                               info_dict.setdefault('Cigar', None),
                               info_dict.setdefault('Hits', None), info_dict.setdefault('Alias', None),
                               info_dict.setdefault('Genomic', None),
                               info_dict.setdefault('Filter', None), info_dict.setdefault('Seed_fam', None)]
                expression_list = info_dict.get('Expression').split(',')
                complete_list = prefix_list + expression_list
                gff_insert_values(conn, complete_list)

        c.execute('''CREATE TABLE IF NOT EXISTS summary(version text, source text, data_sets text, tools text,
         commands_exec text, filter_tags text, citation text, records real, date_stamp text)''')
        c.execute("INSERT INTO summary(version, source, data_sets, tools, commands_exec, filter_tags, citation, "
                  "records, date_stamp) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (version, source, data_sets, tools, commands_exec, filter_tags, citation, cnt, d2))

        # info_dict.setdefault('Sex', None)

        conn.commit()


def query_sql(args):
    print("Function query is being implemented, will be updated soon!!!")
    print(args)
    if args.db:
        out_file = op.join(args.db)
        conn = sqlite3.connect(out_file)
        c = conn.cursor()
    else:
        out_file = op.join("mirtop.db")
        conn = sqlite3.connect(out_file)
        c = conn.cursor()

    c.execute("SELECT name FROM sqlite_master WHERE type = 'table';")
    record = c.fetchall()
    print(args.db)
    if args.expr == "show-tables":
        show_tables(conn)
    if args.expr == "show-schema":
        if args.table:
            show_schema(conn, args.table)
        else:
            print("Error: Require table name")
            print("Usage: mirtop sql --query --db <input_database> -e show-schema -t <table_name>")
    pass


def show_tables(connection):
    print()
    print(" +" + 25 * "-" + " +")
    print(' | Tables                   |')
    print(" +" + 25 * "-" + " +")
    for (tableName,) in connection.execute(
            """
        select NAME from SQLITE_MASTER where TYPE='table' order by NAME;
        """
    ):
        tn_name = len(tableName)
        req_format_space = 25 - tn_name
        print(" | {tn}{format_space}|".format(
            tn=tableName,
            format_space=req_format_space * " "
        ))  # Table name (for each table)
    print(" +" + 25 * "-" + " +")


def show_schema(connection, table_name):
    for (tableName,) in connection.execute(
            """
        select NAME from SQLITE_MASTER where TYPE='table' order by NAME;
        """
    ):
        #  print("{}:".format(tableName))  # Table name (for each table)
        if tableName == table_name:
            print(" +" + 57 * "-" + "+")
            print(' | Sl | Field                         | Type | NULL | Key  |')
            print(" +" + 57 * "-" + "+")
            for (
                    columnID, columnName, columnType,
                    columnNotNull, columnDefault, columnPK,
            ) in connection.execute("pragma table_info('{}');".format(tableName)):
                name_size = len(columnName)
                columnID_size = len(str(columnID))
                required_len = 30 - name_size
                req_col_size = 2 - columnID_size

                print(" | {colSpace}{id} | {name}{space}| {type} | {null} | {pk} |".format(
                    colSpace=req_col_size * " ",
                    id=columnID,
                    name=columnName,
                    space=required_len * " ",
                    type=columnType if columnType else "NULL",
                    null=" not null" if columnNotNull else " NO ",
                    default=" [{}]".format(columnDefault) if columnDefault else "NULL",
                    pk=" *{}".format(columnPK) if columnPK else "    ",
                ))
            print(" +" + 57 * "-" + "+")


def sql_options(args):
    user_options = vars(args)
    if args.create:
        if args.gff:
            insert_sql(args)
        else:
            print("Usage: mirtop sql --create --gff <input.gff> --db <new_db_name> \(Default: mirtop.db\)")
    elif args.query:
        if args.expr:
            #  print("Usage: mirtop sql --query --db <input_database> -e <user_query>")
            query_sql(args)
        else:
            print("Usage: mirtop sql --query --db <input_database> -e <user_query>")
    else:
        print("Usage: mirtop sql -h")
