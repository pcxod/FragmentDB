'''
Created on 09.10.2014

@author: Daniel Kratzert

'''
import sys

__metaclass__ = type  # use new-style classes
import sqlite3
from sqlite3 import OperationalError

__all__ = ['DatabaseRequest', 'FragmentTable', 'Restraints', 'restraint_check']


SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST',
             'L.S.', 'CGLS', 'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF',
             'SIMU', 'RIGU', 'WGHT', 'FVAR', 'DELU', 'SAME', 'DISP', 'LAUE',
             'REM', 'MORE', 'TIME', 'END', 'HKLF', 'OMIT', 'SHEL', 'BASF',
             'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV',
             'CONN', 'BIND', 'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT',
             'DEFS', 'ISOR', 'NCSY', 'SUMP', 'BLOC', 'DAMP', 'STIR', 'MPLA',
             'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE', 'XNPD', 'REST',
             'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST',
             'PART', 'DANG')

RESTRAINT_CARDS = ('SIMU', 'RIGU', 'DELU', 'SAME', 'FREE', 'DFIX', 'BUMP',
                'HFIX', 'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'DANG')


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
  if len(a) == 1:
      a = a+'.'
  if len(b) == 1:
      b = b + '.'
  a_bigram_list = []
  for i in range(len(a)-1):
    a_bigram_list.append(a[i:i+2])
  b_bigram_list = []
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


def restraint_check(restraints, atoms, fragment_name):
  '''
  - Checks if the Atomnames in the restraints are also in
    the atom list of the respective dbentry.
  - Checks wether restraints cards are valid.
  not used atm and not ready
  '''
  atoms = [i[0].upper() for i in atoms]
  restraint_atoms_list = set([])
  for n, line in enumerate(restraints):
    if not line:
        continue
    if line[:4].upper() not in SHX_CARDS:  # only the first 4 characters, because SADI_TOL would be bad
      print('Bad line in header of database entry "{}" found!'.format(n, fragment_name))
      print(' '.join(line))
      sys.exit(False)
    if line[:4].upper() in RESTRAINT_CARDS:
      line = line[5:].split()
      for i in line:
        if i in ('>', '<'):
          continue
        try:
          float(i)
        except(ValueError):
          restraint_atoms_list.add(i)
  for atom in restraint_atoms_list:
    if not atom.upper() in atoms:
      print('Bad atom "{}" in restraints of "{}"'.format(atom, fragment_name))

def call_profile(dbfile):
  import cProfile
  import pstats
  cp = cProfile.Profile()
  db=FragmentTable(dbfile)
  cp.runcall(db._get_fragment, 17)
  pstats.Stats(cp).strip_dirs().sort_stats('time').print_stats(20)
      
      
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
    self.con.text_factory = str
    #self.con.text_factory = bytes
    with self.con:
      # set the database cursor
      self.cur = self.con.cursor()

  def db_request(self, request, *args):
    '''
    Performs a SQLite3 database request with "request" and optional arguments
    to insert parameters via "?" into the database request.
    A push request will return the last row-Id.
    A pull request will return the requested rows
    :param request: sqlite database request like:
                """SELECT fragment.id FROM fragment"""
    :type request: str
    '''
    try:
      if isinstance(args[0], (list, tuple)):
        args = args[0]
    except IndexError:
      pass
    try:
      self.cur.execute(request, args)
      last_rowid = self.cur.lastrowid
    except OperationalError as e:
      #print(e)
      return False
    rows = self.cur.fetchall()
    if not rows:
      return last_rowid
    else:
      return rows

  def __del__(self):
    # commit is very slow:
    self.con.commit()
    self.con.close()
  
class FragmentTable():
  '''
  >>> dbfile = 'fragment-database.sqlite'
  >>> db = FragmentTable(dbfile)
  >>> print db[2]
  [(u'N1', u'7', 20.6124, 12.15, 22.3497), (u'N2', u'7', 20.9206, 10.5397, 20.3473), (u'C1', u'6', 20.425, 12.9255, 23.3244), (u'C2', u'6', 19.2147, 13.3242, 23.9359), (u'C3', u'6', 17.9873, 12.5825, 23.4736), (u'C4', u'6', 18.2353, 11.6971, 22.5007), (u'C5', u'6', 19.5035, 11.446, 21.9744), (u'C6', u'6', 19.6805, 10.5356, 20.7944), (u'C7', u'6', 18.6945, 9.6218, 20.3676), (u'C8', u'6', 18.9823, 8.7613, 19.3477), (u'C9', u'6', 20.2843, 8.7959, 18.8213), (u'C10', u'6', 21.23, 9.6533, 19.2507)]

  >>> for num, i in enumerate(db):
  ...   print(i)
  ...   if num >= 5: break
  (63, u'(1-methyl-1H-imidazol-2-yl)methanol, C5H8N2O')
  (58, u'1,2-Dichlorobenzene, C6H4Cl2')
  (54, u'1,2-Difluorobenzene, C6H4F2')
  (55, u'1,2-Dimethoxyethane, coordinated to Na, C4H10O2, DME')
  (21, u'1,2-Dimethoxyethane, not coordinated, C4H10O2, DME')
  (25, u'1,4-Diazabicyclo[2.2.2]octane, DABCO')

  '''
  def __init__(self, dbfile):
    '''
    Class to modify the database tables of the fragment database in "dbfile"
    :param dbfile: database file path
    :type dbfile: str
    '''
    self.database = DatabaseRequest(dbfile)

  def __contains__(self, fragment_id):
    '''
    Returns a database fragment if its name contains element of type int.
    E.g. db[2]

    >>> dbfile = 'fragment-database.sqlite'
    >>> 2 in FragmentTable(dbfile)
    True

    >>> db = FragmentTable(dbfile)
    >>> if 2 in db:
    ...   print('yes')
    yes

    >>> 'benzene' in FragmentTable(dbfile)
    Wrong type. Expected integer.
    False

    >>> '2' in FragmentTable(dbfile)
    True
    
    >>> 2 in FragmentTable(dbfile)
    True

    :param name: (partial) name of a database fragment.
    :type name: str
    '''
    try:
      fragment_id = int(fragment_id)
    except(ValueError, TypeError):
      print('Wrong type. Expected integer.')
    if self.has_index(fragment_id):
      return True
    else:
      return False

  def __len__(self):
    '''
    Called to implement the built-in function len().
    Should return the number of database entrys.
    
    # number of fragments in the database:
    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> len(db)  
    69
    
    :rtype: int
    '''
    req = '''SELECT Fragment.id FROM Fragment'''
    rows = self.database.db_request(req)
    if rows:
      return len(rows)
    else:
      raise IndexError('Could not determine database size')

  def __getitem__(self, fragment_id):
    '''
    Called to implement evaluation of self[fragment_id].
    print FragmentTable[fragment_id]

    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> print(db[0])
    Traceback (most recent call last):
      ...
    IndexError: Database fragment not found.
    
    >>> for i in db[2]:
    ...   print(i)
    (u'N1', u'7', 20.6124, 12.15, 22.3497)
    (u'N2', u'7', 20.9206, 10.5397, 20.3473)
    (u'C1', u'6', 20.425, 12.9255, 23.3244)
    (u'C2', u'6', 19.2147, 13.3242, 23.9359)
    (u'C3', u'6', 17.9873, 12.5825, 23.4736)
    (u'C4', u'6', 18.2353, 11.6971, 22.5007)
    (u'C5', u'6', 19.5035, 11.446, 21.9744)
    (u'C6', u'6', 19.6805, 10.5356, 20.7944)
    (u'C7', u'6', 18.6945, 9.6218, 20.3676)
    (u'C8', u'6', 18.9823, 8.7613, 19.3477)
    (u'C9', u'6', 20.2843, 8.7959, 18.8213)
    (u'C10', u'6', 21.23, 9.6533, 19.2507)

    >>> dblen = len(db)
    >>> print(db[dblen-1])
    [(u'O1', u'8', -2.0934, 1.1341, -0.085), (u'C1', u'6', -0.7437, 1.2282, -0.5749), (u'C2', u'6', -0.0311, 0.005, 0.0066), (u'C3', u'6', -0.7498, 1.181, -2.1104), (u'C4', u'6', -0.1025, 2.529, -0.0678)]

    >>> print(db[-1])
    [(u'O1', u'8', -2.0934, 1.1341, -0.085), (u'C1', u'6', -0.7437, 1.2282, -0.5749), (u'C2', u'6', -0.0311, 0.005, 0.0066), (u'C3', u'6', -0.7498, 1.181, -2.1104), (u'C4', u'6', -0.1025, 2.529, -0.0678)]

    :param fragment_id: Id number of fragment to return.
    :type fragment_id: int
    '''
    try:
      fragment_id = int(fragment_id)
    except(ValueError, KeyError):
      print('Wrong type. Integer expected')
      sys.exit()
    if fragment_id < 0:
      fragment_id = len(self)-abs(fragment_id)
    found = self._get_fragment(fragment_id)
    if found:
      return found
    else:
      raise IndexError('Database fragment not found.')


  def __delitem__(self, fragment_id):
    '''
    Called to implement deletion of self[fragment_id].
    
    # have to create a copy of the db before I delete the entry:
    >>> import shutil
    >>> shutil.copyfile('fragment-database.sqlite', 'tst.sqlite')
    >>> dbfile = 'tst.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> del db[2]
    >>> db[2]
    Traceback (most recent call last):
      ...
    IndexError: Database fragment not found.
    
    >>> print 'before:', len(db)
    before: 68
    
    # del db[0] deletes nothing:
    
    >>> del db[0]
    >>> print 'after:', len(db)
    after: 68

    >>> del db[3]
    >>> print 'after:', len(db)
    after: 67
    
    >>> print db[-3][1]
    (u'O1', u'8', -1.9563, -0.0969, 0.0001)
    
    >>> del db[-3]
    >>> print db[-3][1]
    (u'C2', u'6', 0.0, 1.4155, 0.0)
    
    :param fragment_id: Id number of fragment to delete.
    :type fragment_id: int
    '''
    # in case of negative id, get a list of the ids and access the id through
    # that list: 
    if fragment_id < 0:
      fragment_id = self.get_all_rowids()[fragment_id]
    req = '''DELETE FROM Fragment WHERE rowid = ?'''
    try:
      fragment_id = int(fragment_id)
    except(ValueError, TypeError):
      print('Wrong type. Expected integer.')
    # actually delete the item:
    deleted = self.database.db_request(req, fragment_id)
    return deleted

  def __iter__(self):
    '''
    This method is called when an iterator is required for FragmentTable.
    Returns the Id and the Name as tuple.

    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> for num, i in enumerate(db):
    ...   print(i)
    ...   if num > 1:
    ...     break
    (63, u'(1-methyl-1H-imidazol-2-yl)methanol, C5H8N2O')
    (58, u'1,2-Dichlorobenzene, C6H4Cl2')
    (54, u'1,2-Difluorobenzene, C6H4F2')
    '''
    all_fragments = self.get_all_fragment_names()
    return iter(all_fragments)

  def has_name(self, name):
    '''
    Returns True if a partial name is found in the DB.
    '''
    req = '''SELECT Name FROM Fragment WHERE Fragment.Name like "%{}%" '''.format(name)
    if self.database.db_request(req):
      return True

  def has_index(self, Id):
    '''
    Returns True if db has index Id
    :param Id: Id of the respective fragment
    :type Id: int
    :rtype: bool
    '''
    req = '''SELECT Id FROM Fragment WHERE Fragment.Id = {}'''.format(Id)
    if self.database.db_request(req):
      return True
  
  def get_all_rowids(self):
    '''
    returns all Ids in the database as list.
    :rtype: list 
    '''
    ids = []
    req = '''SELECT Id FROM Fragment ORDER BY Id'''
    rows = self.database.db_request(req)
    for i in rows:
      ids.append(i[0])
    return ids
  
  def get_all_fragment_names(self):
    '''
    returns all fragment names in the database, sorted by name
    '''
    req = '''SELECT Fragment.Id, Fragment.name FROM Fragment ORDER BY Name'''
    rows = self.database.db_request(req)
    return rows

  def _get_fragment(self, fragment_id):
    '''
    returns a full fragment with all atoms, atom types as a tuple.
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    :rtype: tuple
    '''
    req_atoms = '''SELECT Atoms.name, Atoms.element, Atoms.x, Atoms.y, Atoms.z
      FROM Fragment, Atoms on Fragment.Id=Atoms.FragmentId WHERE
      Fragment.Id = {}'''.format(fragment_id)
    atomrows = self.database.db_request(req_atoms)
    return atomrows

  def get_fragment_name(self, fragment_id):
    '''
    returns the "Name" column entry of fragment with id "fragment_id"

    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)    
    >>> db.get_fragment_name(2)
    [(u"2,2'-Bipyridine, C10H8N2, bipy",)]
    
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    '''
    req_name = '''SELECT Fragment.Name FROM Fragment WHERE Fragment.Id = {}
               '''.format(fragment_id)
    name = self.database.db_request(req_name)
    return name
  
  def get_picture(self, fragment_id):
    '''
    returns a picture of the fragment if one exist in the database. Otherwise 
    it returns False.
    
    dbfile = 'fragment-database.sqlite'
    db = FragmentTable(dbfile)    
    pic = db.get_picture(2)
    '''
    req_picture = '''SELECT Fragment.picture FROM Fragment WHERE Fragment.Id = {}
               '''.format(fragment_id)
    picture = self.database.db_request(req_picture)[0][0]
    return picture
  
  def get_residue_class(self, fragment_id):
    '''
    returns the "class" column entry of fragment with id "fragment_id"

    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)    
    >>> db.get_residue_class(2)
    'BIP'
    
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    '''
    req_class = '''SELECT Fragment.class FROM Fragment WHERE Fragment.Id = {}
               '''.format(fragment_id)
    classname = self.database.db_request(req_class)
    try:
      classname = classname[0][0]
    except(IndexError):
      print('Could not find residue class.')
      return
    return str(classname) 

  def get_restraints(self, fragment_id):
    '''
    returns the restraints for Fragment(Id) from the database.
    
    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)    
    >>> db.get_restraints(5)
    [(u'DFIX', u'O1 C3 1.2623'), (u'DFIX', u'O2 C3 1.2623'), (u'DFIX', u'C3 C4 1.5633'), (u'DFIX', u'O1 O2 2.2755'), (u'DFIX', u'O1 C4 2.3970'), (u'DFIX', u'O2 C4 2.3970'), (u'FLAT', u'O1 > C4'), (u'SIMU', u'O1 > C4'), (u'RIGU', u'O1 > C4')]
    
    :param fragment_Id: id of the fragment in the database
    :type fragment_Id: int
    '''
    req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
            FROM Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return restraintrows

  def find_fragment_by_name(self, name, selection=5):
    '''
    find a fragment by its name in the database. This method will output a
    selection of (default=5) best hits.
    
    >>> dbfile = 'fragment-database.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> db.find_fragment_by_name('cf3', selection=3)
    [(3, u'Trifluoroethanol, OCH2CF3-'), (36, u'Nonafluoro-tert-butoxy, [(CF3)3CO]-'), (49, u'Trifluoromethanesulfonate, CF3SO3-, Triflate')]
    
    :param name: (part of) the name of a fragment to find
    :type name: str
    :return fragment_id: list of id numbers of the found fragments e.g. [1, 5]
    :type fragment_id: int
    '''
    return self._search_name(name, selection)

  def _search_name(self, search_string, selection=5):
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
    for i in self:
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
    :rtype int: FragmentId -> last_rowid
    '''
    # first stores the meta-information in the Fragment table:
    # The FragmentId is the last_rowid from sqlite
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
    :rtype list: last_rowid
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
  import doctest
  failed, attempted = doctest.testmod()
  if failed == 0:
    print('passed all tests!')

  # import cProfile
  dbfile = 'fragment-database.sqlite'
#  call_profile(dbfile)
  db = FragmentTable(dbfile)
  db.get_picture(2)

  atoms = [[u'C1', u'6', 1.2, -0.023, 3.615], (u'C2', u'6', 1.203, -0.012, 2.106), (u'C3', u'6', 0.015, -0.011, 1.39), (u'C4', u'6', 0.015, -0.001, 0.005), (u'C5', u'6', 1.208, 0.008, -0.688), (u'C6', u'6', 2.398, 0.006, 0.009), (u'C7', u'6', 2.394, -0.004, 1.394)]
  #atoms = ['C1 6 1.2 -0.023 3.615', 'C2 6 1.203 -0.012 2.106', 'C3 6 0.015 -0.011 1.39', 'C4 6 0.015 -0.001 0.005', 'C5 6 1.208 0.008 -0.688', 'C6 6 2.398 0.006 0.009', 'C7 6 2.394 -0.004 1.394']
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

 
