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


def match_dbfrag(fragId=17):
  dbfile = 'C:\Program Files\Olex2-1.2\etc\scripts\dk-database.sqlite'
  db = FragmentTable(dbfile)
  #print(db.get_all_fragment_names())
  for i in db.get_fragment(fragId)[1]:
    #print('adding {}'.format(i[0]))
    print('adding {}  {}  {}  {}  {}'.format(*i))
    label = str(i[0])
    x, y, z = float(i[2]), float(i[3]), float(i[4])
    x, y, z = olx.xf.au.Fractionalise(x,y,z).split(',')
    id = olx.xf.au.NewAtom(label, x, y, z)
    # set occupancy using olx.xf.au.SetAtomOccupancy(id, value)
    #olx.xf.au.SetAtomOccupancy(id, 1)
    # or U/Uiso using olx.xf.au.SetAtomU(id, *[U values])
    olx.xf.au.SetAtomU(id, 0.04)
  olx.xf.EndUpdate()

OV.registerFunction(match_dbfrag)

def find_frag(name):
  dbfile = 'C:\Program Files\Olex2-1.2\etc\scripts\dk-database.sqlite'
  db = FragmentTable(dbfile)
  for num, name in db.find_fragment_by_name(name):
    print(num, name)


OV.registerFunction(find_frag)