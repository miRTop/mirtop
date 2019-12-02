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
    'iso_5p real', 'iso_3p real', 'iso_add3p real', 'iso_add5p real', 'iso_snv real', 'iso_snv_seed real', 'iso_snv_central real', 'iso_snv_central_offset real',
    'iso_snv_central_supp real', 'source text', 'cigar text', 'hits real', 'alias text', 'genomic_pos text', 'filter text',
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
                info_dict = dict()
                for elements in info:
                    (k, v) = elements.split(' ')
                    info_dict.update([(k, v)])
                    if 'Variant' in k and ":" in v:
                        value_list = v.split(',')
                        for iso_vars in value_list:
                            if ":" in iso_vars:
                                (sub_k, sub_v) = iso_vars.split(':')
                                info_dict.update([(str(sub_k), str(sub_v))])
                            else:
                                ### Exception for miRge format START 
                                if iso_vars == "iso_snp":
                                    iso_vars = "iso_snv"
                                elif iso_vars == "iso_add": 
                                    iso_vars = "iso_add3p"
                                ### Exception for miRge format END
                                info_dict['iso_snv'] = "1" if iso_vars == 'iso_snv' else 0
                                info_dict['iso_snv_seed'] = "1" if iso_vars == 'iso_snv_seed' else 0
                                info_dict['iso_snv_central'] = "1" if iso_vars == 'iso_snv_central' else 0
                                info_dict['iso_snv_central_offset'] = "1" if iso_vars == 'iso_snv_central_offset' else 0
                                info_dict['iso_snv_central_supp'] = "1" if iso_vars == 'iso_snv_central_supp' else 0
                                 
                prefix_list = [lines[0], lines[1], lines[2], lines[3], lines[4], lines[5], lines[6], lines[7],
                               info_dict.get('UID'), info_dict.get('Read'), info_dict.get('Name'),
                               info_dict.get('Parent'),
                               info_dict.get('Variant'), 
                               str(info_dict.setdefault('iso_5p', None)),
                               str(info_dict.setdefault('iso_3p', None)), 
                               str(info_dict.setdefault('iso_add3p', None)),
                               str(info_dict.setdefault('iso_add5p', None)),
                               str(info_dict.setdefault('iso_snv', "0")), 
                               str(info_dict.setdefault('iso_snv_seed', "0")), 
                               str(info_dict.setdefault('iso_snv_central', "0")), 
                               str(info_dict.setdefault('iso_snv_central_offset', "0")), 
                               str(info_dict.setdefault('iso_snv_central_supp', "0")), 
                               source,
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
    #print("Function query is being implemented, will be updated soon!!!")
    #print(args)
    if args.db:
        out_file = op.join(args.db)
        conn = sqlite3.connect(out_file)
        c = conn.cursor()
    else:
        out_file = op.join("mirtop.db")
        conn = sqlite3.connect(out_file)
        c = conn.cursor()

    #c.execute("SELECT name FROM sqlite_master WHERE type = 'table';")
    #record = c.fetchall()
    #print(args.db)
    if args.expr == "show-tables":
        show_tables(conn)
    if args.expr == "show-schema":
        if args.table:
            show_schema(conn, args.table)
        else:
            print("Error: Require table name")
            print("Usage: mirtop sql --query --db <input_database> -e show-schema -t <table_name>")
    if args.expr == "show-columns":
        show_columns(conn, args)
    if args.expr == "describe-gff":
        describe_gff_info(conn, args)
    if args.expr == "isomirs-per-mirna":
        if args.miRNA:
            stats_isomiR_per_miRNA(conn, args.miRNA, args)
        else:
            print("Error: Require miRNA name")
            print("Usage: mirtop sql --query --db <input_database> -e isomirs-per-mirna -miR <miRNA>")
    if args.expr == "select":
        select_query(conn, args)
    pass


def show_tables(connection):
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


def show_columns(connection, args):
    cur = connection.cursor()
    query="SELECT * FROM data_sets LIMIT 2"
    cur.execute(query)
    rows = cur.fetchall()
    col_name_list = [tuple[0] for tuple in cur.description]
    sl_no=1
    print("\nSerial\tColumn names")
    for col in col_name_list:
        print("  "+str(sl_no)+"\t"+str(col))
        sl_no+=1
    print()


def describe_gff_info(connection, args):
    cur = connection.cursor()
    query="SELECT * FROM summary"
    cur.execute(query)
    rows = cur.fetchone()
    col_name_list = [tuple[0] for tuple in cur.description]
    sl_no=1
    print("\nSerial\tColumn names\tDescription")
    for i, col in enumerate(col_name_list):
        desc_g = rows[i]
        if (desc_g ==""):
            desc_g = "--"
        print("  "+str(sl_no)+"\t"+str(col)+"\t\""+str(desc_g)+"\"")
        sl_no+=1
    print()


def stats_isomiR_per_miRNA(connection, miRNA_name, args):
    cur = connection.cursor()
    miR_array = add_mirnas(args)
    #cur.execute('SELECT * FROM data_sets WHERE seqID=?', (miRNA_name,))
    query="SELECT COUNT(*) FROM data_sets WHERE seqID=? AND type='isomiR' "
    query = add_filter(query, args)
    print()
    stat_counts=0
    if args.txtout:
        print("The results are being fetched and formated to be written to "+ args.txtout)
        #format_results()
        #with open(args.txtout, 'w') as w_stat:
        w_stat = open(args.txtout, 'w')
        w_stat.write("Serial number\tmiRNA\tisomiR Count\n")
    else:
        print("OUTPUT:")
    
    for miRs in miR_array:
        t=(miRs, )
        cur.execute(query, t)
        #cur.execute("SELECT COUNT(*) FROM data_sets WHERE seqID=? AND type='isomiR'", t)
        rows = cur.fetchall()
        for row in rows:
            stat_counts+=1
            row = row[0]
            if args.txtout:
                w_stat.write(str(stat_counts) +"\t"+miRs+"\t"+str(row)+"\n")
            else:
                print(str(stat_counts) +". " +"isomiRs for miRNA "+ miRs + ": "+ str(row))
    print()
    if args.txtout:
        w_stat.close()
    pass


def WHERE_CLAUSE(query, args):
    if "WHERE" in query: 
        query = query + " AND "
        return(query)
    else:
        query = query + " WHERE " 
        return(query)


# ALWAYS EXECUTE THIS LIMIT FUNCTION AT THE END 
def add_limit(query, args):
    if args.limit:
        query = query + " LIMIT "+ args.limit 
        return query
    else:
        return query


def add_filter(query, args):
    if args.filter:
        query = WHERE_CLAUSE(query, args)
        query = query + " filter='" + args.filter +"' "
        return query
    else:
        return query


def add_variants(query, args):
    my_var_dict = {'iso_5p': 'iso_5p != "None"', 'iso_3p': 'iso_3p != "None"', 'iso_add3p': 'iso_add3p != "None"', 'iso_add5p': 'iso_add5p != "None"', 
            'iso_snv_seed':'iso_snv_seed != 0', 'iso_snv_central_offset':'iso_snv_central_offset != 0', 'iso_snv_central':'iso_snv_central != 0', 
            'iso_snv_central_supp':'iso_snv_central_supp != 0', 'iso_snv':'iso_snv != 0'}
    user_req_var = args.variant.split(',')
    values_req_var =[]
    for eachVar in user_req_var:
        try:
            values_req_var.append(my_var_dict[eachVar])
        except KeyError:
            print("\nError: \"" + eachVar + "\" does not exist in the choices supported by (-var , --variant)\n")
            print("use: mirtop sql -qh for more options")
            exit()
    #print(values_req_var)
    insert_betwn = " AND "
    query_suffix = (insert_betwn.join( values_req_var ))
    query = WHERE_CLAUSE(query, args)
    query = query + query_suffix
    return query
    #if args.filter:
        #query = query + " AND " + query_suffix 
        #return query
    #else:
        #query = query + "WHERE " + query_suffix 
        #return query

def add_mirnas(args):
    if args.miRNA.endswith('.txt'):
        #print("I am called and I am safe here to read from a file")
        with open(args.miRNA, 'r') as miList:
            miR_array = miList.read().splitlines()
        return(miR_array)
    else:
        miR_array=args.miRNA.split(',')
        return(miR_array)


def perform_execution(conn, query, args):
    cur = conn.cursor()
    #print("QUERY: \n"+ query + "\n")
    cur.execute(query)
    rows = cur.fetchall()
    col_name_list = [tuple[0] for tuple in cur.description]
    if args.miRNA:
        return(col_name_list, rows)
    else:
        format_results(col_name_list, rows, args)

def format_results(header, output, args):
    header = '\t'.join(str(col) for col in header)
    if args.txtout:
        outList = open(args.txtout, 'w')
        print("\nWriting data to file: "+ args.txtout + "\n")
        outList.write(header+"\n")
        write_to_file(output, args, outList)
    else:
        print(header)
        if args.count:
            output = list(output[0])
            if args.columns:
                print("Unique counts for "+ str(args.columns) + " is: " + str(output[0]))
            else:
                print("Unique counts for all rows is: " + str(output[0]))
        else:
            for eachrow in output:
                row_tab = '\t'.join(str(items) for items in eachrow)
                print(row_tab)


def write_to_file(output, args, fileHandler):
        if args.miRNA:
            fileHandler.write(output + "\n")
        elif args.count:
            output = list(output[0])
            fileHandler.write("Unique counts for "+ args.columns + " is:\t" + str(output))
        else:
            for eachrow in output:
                row_tab = '\t'.join(str(items) for items in eachrow)
                fileHandler.write(row_tab+"\n")


def select_query(connection, args):
    if args.columns:
        if args.count: 
            if args.count == "T":
                query = "SELECT COUNT(" + args.columns + ") FROM data_sets "
            else:
                print("\nERROR: -n is incorrect!. \nPlease use -n T and optionally specify any one column in -col.\nFor more options see mirtop sql -h")
                exit()
        else:
            query = "SELECT " + args.columns + " FROM data_sets "
    elif args.count:
        query = "SELECT COUNT(*) FROM data_sets "
    else:
        query = "SELECT * FROM data_sets "
    query = add_filter(query, args)
    if args.variant:
        query = add_variants(query, args)
    if args.miRNA:
        miR_array = add_mirnas(args)
        query = WHERE_CLAUSE(query, args)
        query_series = query +  "seqID= "
        header_var= ""
        j=0
        if args.txtout:
            outList = open(args.txtout, 'w')
            print("\nWriting data to file: "+ args.txtout + "\n")

        for miRs in miR_array:
            query = query_series + "\"" + miRs + "\" "
            query = add_limit(query, args)
            (header, rows) = perform_execution(connection, query, args)
            j += 1
            for i, row in enumerate(rows):
                if (i == 0):
                    if header_var == "":
                        header_var = '\t'.join(str(col) for col in header)
                        if args.txtout:
                            outList.write(header_var+"\n")
                        else:
                            print(header_var)
                row_tab = '\t'.join(str(items) for items in row)
                if args.count:
                    newOut = str(j) + ". "+ miRs + ":\t" + row_tab
                    if args.txtout:
                        write_to_file(newOut, args, outList)
                    else:
                        print(newOut)
                else:
                    if args.txtout:
                        write_to_file(row_tab, args, outList)
                    else:
                        print(row_tab)
    else:
        query = add_limit(query, args)
        perform_execution(connection, query, args)


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
