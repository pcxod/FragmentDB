import olex
import olx
#import sys
from fragmentdb import FragmentTable
from olexFunctions import OlexFunctions
OV = OlexFunctions()


'''
To run this example script, type spy.example("Hello") in Olex2
'''

def example(text="No Text Provided"):
  print "Example Function is now printing your text: %s" %text

OV.registerFunction(example)

def dktest():
  datadir = olex.f("DataDir()")
  basedir = olex.f("BaseDir()")
  print('hello world')
  print(datadir)
  print(basedir)

OV.registerFunction(dktest)


def match_tol(frag='toluene'):
  dbfile = 'F:\Programme\Olex2-1.2-dev\etc\scripts\dk-database.sqlite'
  db = FragmentTable(dbfile)
  #print(db.get_all_fragment_names())
  for i in db.get_fragment(49)[1]:
    print(i)
    print('{}  {}  {}  {}  {}'.format(*i))

OV.registerFunction(match_tol)

