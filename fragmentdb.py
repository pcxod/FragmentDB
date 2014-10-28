'''
Created on 09.10.2014

@author: Daniel Kratzert

Usage:
dbfile = 'database_file.sqlite'
db = FragmentTable(dbfile)

db[2]   # get database fragment number 2
----------
del db[3]  # deletes fragment 3.
----------
if 'Benzene' in db: # test if name 'Benzene' is in the database
  pass
----------
for i in db:
  print(i)  # get all fragment names and their id number.
----------
for i in db[2]:
  print(i)  # get all atoms of fragment number 2
----------
for i in db.get_restraints(15):
  print(i)   # get restraints of fragment number 15
----------  
len(db)  # number of fragments in the database
----------
# store a fragment into the database.
id = db.store_fragment(fragment_name, atoms, restraints, tag)
if id:
  print('fragment is stored successfully')
----------  
db.find_fragment_by_name('super', selection=3) # find a fragment with 
                                                 foult tolerance search

def match_dbfrag(fragId):
  for i in db[fragId]:
    label = str(i[0])
    x, y, z = olx.xf.au.Fractionalise(i[2],i[3],i[4]).split(',')
    id = olx.xf.au.NewAtom(label, x, y, z, False)
    print('adding {}, Id: {}, coords: {} {} {}'.format(i[0], id, x, y, z))
    olx.xf.au.SetAtomPart(id, -1)
    olx.xf.au.SetAtomOccu(id, 1)
    olx.xf.au.SetAtomU(id, 0.04)
  olx.xf.EndUpdate()
  print('Now you can fit the fragment with "mode fit"')                                                 
'''

__metaclass__ = type  # use new-style classes
import sqlite3
import sys
from sqlite3 import OperationalError
#print(sys.version)
#print


__all__ = ['DatabaseRequest', 'FragmentTable', 'Restraints' ]

def dice_coefficient(a, b):
  '''
  dice coefficient 2nt/na + nb
  Compares the similarity of a and b
  :param a: string
  :param b: string
  '''
  a = a.lower()
  b = b.lower()
  if not len(a) or not len(b): return 0.0
  if len(a) == 1:  a=a+'.'
  if len(b) == 1:  b=b+'.'
  a_bigram_list=[]
  for i in range(len(a)-1):
    a_bigram_list.append(a[i:i+2])
  b_bigram_list=[]
  for i in range(len(b)-1):
    b_bigram_list.append(b[i:i+2])
  a_bigrams = set(a_bigram_list)
  b_bigrams = set(b_bigram_list)
  overlap = len(a_bigrams & b_bigrams)
  dice_coeff = overlap * 2.0/(len(a_bigrams) + len(b_bigrams))
  dice_coeff = 1-dice_coeff # invert the result
  if dice_coeff < 0.5:  # make a cutoff for the best matches
    return 0.0
  return round(dice_coeff, 6)

SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST',
             'L.S.', 'CGLS', 'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF',
             'SIMU', 'RIGU', 'WGHT', 'FVAR', 'DELU', 'SAME', 'DISP', 'LAUE',
             'REM',  'MORE', 'TIME', 'END',  'HKLF', 'OMIT', 'SHEL', 'BASF',
             'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV',
             'CONN', 'BIND', 'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT',
             'DEFS', 'ISOR', 'NCSY', 'SUMP', 'BLOC', 'DAMP', 'STIR', 'MPLA',
             'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE', 'XNPD', 'REST',
             'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST',
             'PART', 'DANG')

class DatabaseRequest():
  def __init__(self, dbfile):
    '''
    creates a connection and the cursor to the SQLite3 database file "dbfile".
    :param dbfile: database file
    :type dbfile: str
    '''
    # open the database
    self.con = sqlite3.connect(dbfile)
    self.con.execute("PRAGMA foreign_keys = ON")
    #self.con.text_factory = str
    with self.con:
      # set the database cursor
      self.cur = self.con.cursor()


  def db_request(self, request, *args):
    '''
    Performs a SQLite3 database request with "request" and optional arguments
    to insert parameters via "?" into the database request.
    A push request will return the last row-Id.
    A pull request will return the requested rows
    :param request: sqlite database request like: """SELECT fragment.id FROM fragment"""
    :type request: str
    '''
    try:
      if isinstance(args[0], (list, tuple)):
        args = args[0]
    except IndexError:
      pass
#      print request
    try:
      self.cur.execute(request, args)
      last_rowid = self.cur.lastrowid
    except OperationalError as e:
      print(e)
      return False
    rows = self.cur.fetchall()
    if not rows:
      self.con.commit()
      return last_rowid
    else:
      return rows


class FragmentTable():
  def __init__(self, dbfile):
    '''
    Class to modify the database tables of the fragment database in "dbfile"
    :param dbfile: database file path
    :type dbfile: str
    '''
    self.database = DatabaseRequest(dbfile)

  def __contains__(self, name):
    '''
    Returns a database fragment if its name contains "name".
    if 'benzene' in db:
      print('yes')
    :param name: (partial) name of a database fragment.
    :type name: str
    '''
    all_fragments = self.get_all_fragment_names()
    if not all_fragments:
      print('No database items found!')
      return False
    all_names = (i[1].lower() for i in all_fragments)
    found = False
    if isinstance(name, (basestring)):
      name = name.lower()
      found = any(name in s for s in all_names)
      return found
    else:
      raise TypeError('Name must be string.')

  def __len__(self):
    '''
    Called to implement the built-in function len().
    Should return the number of database entrys.
    '''
    req = '''SELECT fragment.id FROM fragment'''
    rows = self.database.db_request(req)
    if rows:
      return len(rows)
    else:
      return False

  def __getitem__(self, fragment_id):
    '''
    Called to implement evaluation of self[fragment_id].
    print FragmentTable[fragment_id]
    :param fragment_id: Id number of fragment to return.
    :type fragment_id: int
    '''
    if self.fragment_id_is_valid(fragment_id):
      found = self._get_fragment(fragment_id)
      if found:
        return found
      else:
        raise IndexError('Database fragment not found.')
    else:
      return False

  def __delitem__(self, fragment_id):
    '''
    Called to implement deletion of self[fragment_id].
    del FragmentTable[3]
    :param fragment_id: Id number of fragment to delete.
    :type fragment_id: int
    '''
    deleted = False
    req = '''DELETE FROM fragment WHERE fragment.id = ?'''
    if self.fragment_id_is_valid(fragment_id):
      deleted = self.database.db_request(req, fragment_id)
    return deleted

  def __iter__(self):
    '''
    This method is called when an iterator is required for FragmentTable.
    Returns the Id and the Name as tuple. 
    e.g.:
    (21, u'Octane, C8H18')
    (33, u'PM4, C37H32B2N8O9')
    for i in db:
      print(i)
    '''
    all_fragments = self.get_all_fragment_names()
    if all_fragments:
      #return (i[1] for i in all_fragments)
      return (i for i in all_fragments)
    else:
      return False

  def fragment_id_is_valid(self, fragment_id):
    '''
    Returns True if the supplied fragment id is valid for a database request.
    :param fragment_id: a primary database key
    :type fragment_id: int (str if convertable to int)
    '''
    valid = False
    if isinstance(fragment_id, int):
      valid = True
    if isinstance(fragment_id, float):
      valid = False
    if isinstance(fragment_id, (basestring)):
      try:
        int(fragment_id)
      except(SyntaxError, KeyError, ValueError):
        valid = False
        return valid
      valid = True
    return valid

  def get_all_fragment_names(self):
    '''
    returns all fragment names in the database, sorted by name
    '''
    req = '''SELECT fragment.id, fragment.name FROM fragment ORDER BY name'''
    rows = self.database.db_request(req)
    return rows

  def _get_fragment(self, fragment_id):
    '''
    returns a full fragment with all atoms, atom types as a tuple.
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    '''
    req_atoms = '''SELECT Atoms.name, Atoms.element, Atoms.x, Atoms.y, Atoms.z
      FROM Fragment, Atoms on Fragment.Id=Atoms.FragmentId WHERE
      Fragment.Id = {}'''.format(fragment_id)
    atomrows = self.database.db_request(req_atoms)
    return (atomrows)

  def get_restraints(self, fragment_id):
    '''
    returns the restraints for Fragment(Id) from the database.
    :param fragment_Id: id of the fragment in the database
    :type fragment_Id: int
    '''
    req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
            FROM Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return (restraintrows)

  def find_fragment_by_name(self, name, selection=5):
    '''
    find a fragment by its name in the database. This method will output a
    selection of (default=5) best hits.
    :param name: (part of) the name of a fragment to find
    :type name: str
    :return fragment_id: list of id numbers of the found fragments e.g. [1, 5]
    :type fragment_id: int
    '''
    req = '''SELECT fragment.id, fragment.name FROM fragment'''
    frags = self.database.db_request(req)
    fragment_ids = self._search_name(name, frags, selection)
    return fragment_ids

  def _search_name(self, search_string, frags, selection=5):
    '''
    searches the names in the database for a given name
    :param search_string: search for this string
    :type search_string: str
    :param frags: fragments in the database
    :type frags: list
    :param selection: return this amount of results
    :type selection: int
    '''
    search_results = {}
    for i in frags:
      db_entry = i[1]
      coefficient = dice_coefficient(search_string, db_entry)
      search_results[coefficient] = i
    # select the best [selection] results:
    selected_results = [search_results[i] for i in sorted(search_results)[0:selection]]
    return selected_results

  def store_fragment(self, fragment_name, atoms, restraints=None, tag=None,
                     reference=None, comment=None):
    '''
    Store a complete new fragment into the database. Minimal requirement is a
    fragment name (Full chemical name) and a list of atoms. Restraints, short
    name tag, reference and comments are optional.

    :param fragment_name: full chemical name of the fragment
    :type fragment_name: string
    :param atoms: atoms of the fragment ['name', 'atomic number', 'x', 'y', 'z']
    :type atoms: list
    :param restraints: [['DFIX', '1.564', 'C1 C2'], ['SADI', 'C2 C3 C4 C5']]
    :type restraints: list of list
    :param tag: short name tag (not mandatory)
    :type tag: string
    :param reference: optional short description where it came from
    :type reference: string
    :param comment: optional any comment about the fragment
    :type comment: string
    '''
    # first stores the meta-information in the Fragment table:
    FragmentId = self._fill_fragment_table(fragment_name, tag, reference, comment)
    if not FragmentId:
      raise Exception('No Id obtained during fragment storage.')
    # then stores atoms with the previously obtained FragmentId
    self._fill_atom_table(FragmentId, atoms)
    # in case of supplied restraints store them also:
    if restraints:
      self._fill_restraint_table(FragmentId, restraints)
    return FragmentId


  def _fill_fragment_table(self, fragment_name, tag=None, reference=None, comment=None):
    '''
    Fills a fragment into the database.
    :param fragment_name: Nam eof the Fragment
    :type fragment_name: str
    :param tag: short name tag (not mandatory)
    :type tag: str
    :param comment: any comment about the fragment
    :type comment: str
    '''
    table = (tag, fragment_name, reference, comment)
    req = '''INSERT INTO Fragment (tag, name, reference, comment) VALUES(?, ?, ?, ?)'''
    return self.database.db_request(req, table)


  def _fill_atom_table(self, FragmentId, atom_table, tag=None):
    '''
    Fills atoms into the Atoms table.
    [('C1', '6', 1.2, -0.023, 3.615), ('C2', '6', 1.203, -0.012, 2.106), ...]
    or
    [('C1 6 1.2 -0.023 3.615'), ('C2 6 1.203 -0.012 2.106'), ...]
    :param FragmentId: Id of the respective Fragment
    :type FragmentId: int or str
    :param atom_table: list of lits or list of strings
    :type atom_table: list
    '''
    # test wether atom_table is a list or list of list, because we want no string
    # in a list here.
    if not isinstance(atom_table[0], (list, tuple)):
      if isinstance(atom_table[0], str):
        atom_table = [i.split() for i in atom_table]
      else:
        raise Exception('wrong data type "{}" for atom list.'.format(type(atom_table[0])))
    for line in atom_table:
        Name = line[0]
        element = line[1]
        x = line[2]
        y = line[3]
        z = line[4]
        req = '''INSERT INTO atoms (FragmentId, tag, Name, element,
                x, y, z) VALUES(?, ?, ?, ?, ?, ?, ?)'''
        self.database.db_request(req, (FragmentId, tag, Name, element, x, y, z))


  def _fill_restraint_table(self, FragmentId, restraints_list):
    '''
    Fills the restraints table with restraints. restraints_list must be a list 
    or tuple of string lists like:
    [['SADI C1 F1 C2 F2'],
    ['SADI 0.04 F2 C3 F1 C6 F2 C1 F1 C2'],
    ['SAME C2 > C6 C1']]
    :param FragmentId: Fragment Id where the restraints belong to.
    :type FragmentId: int or str
    :param restraints_list: list with restraints
    :type restraints_list:
    '''
    # test if restraint_list is a list of strings. we dont want list of list here.
    try:
      restraints_list[0] # is there even one restraint?
    except KeyError:
      print('No restraints found in input.')
      return False
    if isinstance(restraints_list[0], (list, tuple)):
      # convert to list of stings
      restraints_list = [' '.join(['{}'.format(a) for a in i]) for i in restraints_list]
    if isinstance(restraints_list[0], str):
      pass
    else:
      raise Exception('wrong data type "{}" for restraint list.'.format(
                                                  type(restraints_list[0])))
    for line in restraints_list:
      restr_table = []
      if line[:4] in SHX_CARDS:
        restr_table.append(str(FragmentId))
        restr_table.append(line[:4])
        restr_table.append(line[5:])
        req = '''INSERT INTO Restraints (FragmentId, ShelxName, atoms)
                  VALUES(?, ?, ?)'''
        self.database.db_request(req, restr_table)


class Restraints():
  def __init__(self, dbfile):
    self.database = DatabaseRequest(dbfile)

  def get_restraints_from_fragmentId(self, fragment_id):
    req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
      FROM Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return restraintrows



if __name__ == '__main__':
  #dbfile = 'F:\GitHub\DSR-db\dk-database_2.sqlite'
  #dbfile = 'C:\Users\daniel\Documents\GitHub\DSR-db\dk-database_2.sqlite'
  dbfile = 'dk-database_2.sqlite'
  db = FragmentTable(dbfile)

  #print(db[5])
  #def match_dbfrag(fragId=17):
  #  for i in db[fragId]:
  #    print(i)

  #match_dbfrag(57)

  #print(db.get_restraints(57))
  #atoms = [[u'C1', u'6', 1.2, -0.023, 3.615], (u'C2', u'6', 1.203, -0.012, 2.106), (u'C3', u'6', 0.015, -0.011, 1.39), (u'C4', u'6', 0.015, -0.001, 0.005), (u'C5', u'6', 1.208, 0.008, -0.688), (u'C6', u'6', 2.398, 0.006, 0.009), (u'C7', u'6', 2.394, -0.004, 1.394)]
  atoms = ['C1 6 1.2 -0.023 3.615', 'C2 6 1.203 -0.012 2.106', 'C3 6 0.015 -0.011 1.39', 'C4 6 0.015 -0.001 0.005', 'C5 6 1.208 0.008 -0.688', 'C6 6 2.398 0.006 0.009', 'C7 6 2.394 -0.004 1.394']

  #print([' '.join(['{}'.format(a) for a in i]) for i in atoms])

  table = ['sBenzene', 'super Benzene',
           'Name: Trisxdfhdcxh, [(CH3)3Si]3Si, Sudfhdle\nSrc: CCDC CEYMID']
  restraints2 = (('SADI C1 F1 C2 F2'),
                ('SADI 0.04 F2 C3 F1 C6 F2 C1 F1 C2'),
                ('SAME C2 > C6 C1'),
                ('FLAT C1 > F2'),
                ('SIMU C1 > F2'),
                ('RIGU C1 > F2'))
  restraints = (('SADI', 'C1', 'F1', 'C2', 'F2'),
                ('SADI', 0.04, 'F2', 'C3', 'F1', 'C6', 'F2', 'C1', 'F1', 'C2'),
                ('SAME', 'C2 > C6', 'C1'),
                ('FLAT', 'C1 > F2'),
                ('SIMU', 'C1 > F2'),
                ('RIGU' 'C1 > F2'))
  reference = 'sdfg ayrfgrawg adrgaegh ef'
  fragment_name=table[1]
  tag= 'benz'
  comment = 'asfgagr'
  #id = db.store_fragment(fragment_name, atoms, restraints2, tag, reference, comment)
  #if id:
  #  print('stored', id)
  
  #del db[id]
  id = '2'
  print(db[id])

  #for i in db:
  #  print(i)

  # print('len', len(db))
  # if 'benzene' in db:
  #   print('yes')
  
  print(len(db))
  #dbr = DatabaseRequest(dbfile)
  #dbr.db_request('''SELECT * FROM asert''')

  #import cProfile
  #cProfile.run("allnames()", "foo.profile")
  #   res = Restraints(dbfile)

  #for r in db.get_restraints(15):
  #  pass
  #  print(r)

  #   print(hasattr(db, '__iter__'))
 # if 'OC(CF3)3' in db:
 #   print('yes')
 # else:
 #   print('no benz')

  #for i in db:
  #  print(i)

  #print('##################')
  #found = db.find_fragment_by_name('super', selection=3)
  #print(found)
  #for i in db._get_fragment(fragment_id=2):
  #  print(i)
