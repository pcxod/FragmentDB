#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3 as lite
import sys
import os
print(sys.version)

con = None

cars = (
    (1, 'Audi', 52643),
    (2, 'Mercedes', 57127),
    (3, 'Skoda', 9000),
    (4, 'Volvo', 29000),
    (5, 'Bentley', 350000),
    (6, 'Hummer', 41400),
    (7, 'Volkswagen', 21600),
    (8, 'Foo', 123)
)

con = lite.connect('test.db')

def fill_table(con, cars):
    with con:
        cur = con.cursor()    
    
        cur.execute("DROP TABLE IF EXISTS Cars")
        cur.execute('''CREATE TABLE Cars
                    (Id int, Name TEXT, Price INT)''')
        cur.executemany("INSERT INTO Cars VALUES(?, ?, ?)", cars)
        cur.execute("INSERT INTO Cars VALUES ('9','Blub', '10035.14')")


def print_table(con):
    with con:

        con.row_factory = lite.Row #erzeugt ein tuple

        cur = con.cursor()    
        #cur.execute('''SELECT * FROM Cars''')
        cur.execute('SELECT * FROM Cars WHERE Name like "Vo%"')
        rows = cur.fetchall()
        print(rows)
        cur.execute('''SELECT * FROM Cars''')
        rows = cur.fetchall()
        for row in rows:
        #    print(row)
            idd = row["Id"]
            name = row["Name"]
            pri = row["Price"]
            print "{} {} {}".format(idd, name, pri)

    
fill_table(con, cars)
print_table(con)
