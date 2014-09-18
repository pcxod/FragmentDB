#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3 as lite
import sys
import os
print(sys.version)

con = None

cars = (
    (1, 'Audi', 52642),
    (2, 'Mercedes', 57127),
    (3, 'Skoda', 9000),
    (4, 'Volvo', 29000),
    (5, 'Bentley', 350000),
    (6, 'Hummer', 41400),
    (7, 'Volkswagen', 21600)
)

con = lite.connect('test.db')

def fill_table(con, cars):
    with con:
        cur = con.cursor()    
    
        cur.execute("DROP TABLE IF EXISTS Cars")
        cur.execute("CREATE TABLE Cars(Id INT, Name TEXT, Price INT)")
        cur.executemany("INSERT INTO Cars VALUES(?, ?, ?)", cars)


def print_table(con):
    with con:

        con.row_factory = lite.Row

        cur = con.cursor()    
        cur.execute("SELECT * FROM Cars")
    
        rows = cur.fetchall()
    
        for row in rows:
            print "%s %s %s" % (str(row["Id"]), str(row["Name"]), str(row["Price"]))
    
print_table(con)
