import sqlite3
from sqlite3 import Error


def create_connection():
    """ create a database connection to the SQLite database
        specified by the db_file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(":memory:")
        return conn
    except Error as e:
        print(e)
    return None


def create_reads_table(conn, key="sequence"):
    """ create a table on the SQLite database
    :param conn: connection of database
    :return: Connection object or None
    """
    c = conn.cursor()
    c.execute("CREATE TABLE reads"
              " (name text, sequence text,"
              " chrom text, start int,"
              " PRIMARY KEY(%s, chrom, start))" % key)
    conn.commit()


def insert_row_in_reads_table(cur, fields):
    """ create a table on the SQLite database
    :param cur: connection of database
    :param fields: list with columns to fill table
    :return: Connection object or None
    """
    # c = conn.cursor()
    cur.execute("INSERT INTO reads VALUES"
                " (\"%s\", \"%s\", \"%s\", %s)" % (fields[0],
                                                   fields[1],
                                                   fields[2],
                                                   fields[3]))


def select_all_reads(conn):
    """
    Query all rows in the reads table
    :param conn: the Connection object
    :return:
    """
    cur = conn.cursor()
    cur.execute("SELECT * FROM reads")

    rows = cur.fetchall()

    return rows
