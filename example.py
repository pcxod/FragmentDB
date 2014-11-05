import olex
import olx
#import sys
from fragmentdb import FragmentTable
from olexFunctions import OlexFunctions
import sqlite3
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
  print(olx.xf.GetFormula('list'))


OV.registerFunction(dktest)


def olex_functions():
  #print(olx.xf.latt.IsGrown())
  print(help(olx.xf.au.NewAtom))



OV.registerFunction(olex_functions)

dbfile = 'F:\Programme\Olex2-1.2-dev\etc\scripts\dk-database.sqlite'
#dbfile = 'C:\Program Files\Olex2-1.2-dev\Olex2-1.2-alpha\etc\scripts\dk-database.sqlite'
#dbfile = r'C:\Program Files\Olex2-1.2-dev\etc\scripts\dk-database.sqlite'
db = FragmentTable(dbfile)

def all_frags():
  for i in db.get_all_fragment_names():
    print(i)

OV.registerFunction(all_frags)

def dblookup():
  con = sqlite3.connect(dbfile)
  con.execute("PRAGMA foreign_keys = ON")
  #self.con.text_factory = str
  request = '''SELECT * from fragment'''# WHERE fragment.id = 2'''
  with con:
    print('')
    # set the database cursor
    cur = con.cursor()
    cur.execute(request)
    rows = cur.fetchall()
  print(rows)
  
OV.registerFunction(dblookup)

def match_dbfrag(fragId=17):
  atoms = []
  labeldict = {}
  for i in db[fragId]:
    label = str(i[0])
    x, y, z = olx.xf.au.Fractionalise(i[2],i[3],i[4]).split(',')
    id = olx.xf.au.NewAtom(label, x, y, z, False)
    labeldict[label.upper()] = id 
    olx.xf.au.SetAtomPart(id, -1)
    olx.xf.au.SetAtomOccu(id, 1)
    olx.xf.au.SetAtomU(id, 0.04)
    name = olx.xf.au.GetAtomName(id)
    print('adding {}, name: {}, Id: {}, coords: {} {} {}'.format(i[0], name, id, x, y, z))
    atoms.append(id)
  olx.xf.EndUpdate()
  OV.cmd("sel #c{}".format(' #c'.join(atoms)))
  OV.cmd("RESI tst1 9")
  print(labeldict)
  for i in db.get_restraints(fragId):
    if '>' in i[1]:
      continue
    line = []
    for at in i[1].split():
      if at[0].isalpha():
        line.append('#c'+labeldict[at])
      else:
        line.append(at)
    OV.cmd("{} {}".format(i[0], ' '.join(line)))
  OV.cmd("sel #c{}".format(' #c'.join(atoms)))
  OV.cmd("mode fit")
  #print('Now you can fit the fragment with "mode fit"')


OV.registerFunction(match_dbfrag)

def find_frag(name):
  for num, name in db.find_fragment_by_name(name):
    print(num, name)

OV.registerFunction(find_frag)

