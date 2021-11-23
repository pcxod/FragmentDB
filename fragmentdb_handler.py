"""
Created on 09.10.2014

@author: Daniel Kratzert

"""
from __future__ import print_function

import sys

from FragmentDB.helper_functions import dice_coefficient2, SHX_CARDS, make_sortkey

__metaclass__ = type  # use new-style classes
import sqlite3
from sqlite3 import OperationalError

__all__ = ['DatabaseRequest', 'FragmentTable', 'Restraints']


class DatabaseRequest():
  def __init__(self, dbfile, userdb_path=False):
    """
    creates a connection and the cursor to the SQLite3 database file "dbfile".
    :param dbfile: database file
    :type dbfile: str
    """
    # open the database
    self.con = sqlite3.connect(dbfile, check_same_thread=False)
    if userdb_path:
      self.con.execute('ATTACH "{}" AS userdb'.format(userdb_path))
    self.con.execute("PRAGMA foreign_keys = ON")
    # self.con.text_factory = str
    # self.con.text_factory = sqlite3.OptimizedUnicode
    with self.con:
      # set the database cursor
      self.cur = self.con.cursor()

  def db_request(self, request, *args):
    """
    Performs a SQLite3 database request with "request" and optional arguments
    to insert parameters via "?" into the database request.
    A push request will return the last row-Id.
    A pull request will return the requested rows
    :param request: sqlite database request like:
    """
    try:
      if isinstance(args[0], (list, tuple)):
        args = args[0]
    except IndexError as e:
      # print(e, request, args)
      pass
    try:
      self.cur.execute(request, args)
      last_rowid = self.cur.lastrowid
    except OperationalError as e:
      print(e)
      print(request, args)
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
  """
  >>> dbfile = 'tests/tst.sqlite'
  >>> db = FragmentTable(dbfile, './tests/tst-usr.sqlite')
  >>> print(db[3])
  [(u'O1', u'8', -0.7562, 1.6521, -0.3348), (u'C1', u'6', -2.1051, 1.9121, -0.4223), (u'C2', u'6', -2.5884, 1.9919, -1.8717), (u'F1', u'9', -1.9571, 2.9653, -2.5781), (u'F2', u'9', -3.9264, 2.2817, -1.9122), (u'F3', u'9', -2.413, 0.8316, -2.546)]

  >>> for num, i in enumerate(db):
  ...   print(i)
  ...   if num >= 5: break
  [4, u'Acetate anion, C2H3O2-']
  [1000004, u'Acetate anion, C2H3O2-  *user*']
  [47, u'Acetone, C3H6O']
  [28, u'Acetonitrile, C2H3N, NMe']
  [63, u'Adamantane, C10H16']
  [37, u'Adamantane-N, C10H16N']
  """

  def __init__(self, dbfile, userdb_path=''):
    """
    Class to modify the database tables of the fragment database in "dbfile"
    :param dbfile: database file path
    :type dbfile: str
    """
    self.userdb = userdb_path
    self.database = DatabaseRequest(dbfile, userdb_path)

  def __contains__(self, fragment_id):
    """
    Returns a database fragment if its name contains element of type int.
    E.g. db[2]

    >>> dbfile = 'tests/tst.sqlite'
    >>> 3 in FragmentTable(dbfile)
    True

    >>> db = FragmentTable(dbfile)
    >>> if 3 in db:
    ...   print('yes')
    yes

    >>> 'benzene' in FragmentTable(dbfile)
    Traceback (most recent call last):
    ...
    ValueError: invalid literal for int() with base 10: 'benzene'

    >>> '3' in FragmentTable(dbfile)
    True

    >>> 3 in FragmentTable(dbfile)
    True

    >>> 999 in FragmentTable(dbfile)
    False

    >>> "999" in FragmentTable(dbfile)
    False

    :param name: (partial) name of a database fragment.
    :type name: str
    """
    try:
      fragment_id = int(fragment_id)
    except(ValueError, TypeError):
      print('Wrong type. Expected integer.')
    if self.has_index(fragment_id):
      return True
    else:
      return False

  def __len__(self):
    """
    Called to implement the built-in function len().
    Should return the number of database entrys.

    # number of fragments in the database:
    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> len(db)
    71

    :rtype: int
    """
    rows = 0
    num = 0
    for table_name in ['Fragment', 'userdb.Fragment']:
      req = 'SELECT COUNT(*) FROM {}'.format(table_name)
      try:
        rows = self.database.db_request(req)[0][0]
      except TypeError:
        pass
      num = num + rows
    if num:
      return num
    else:
      raise IndexError('Could not determine database size')

  def __getitem__(self, fragment_id):
    """
    Called to implement evaluation of self[fragment_id].
    print FragmentTable[fragment_id]

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> print(db[0])
    False

    >>> print(db[5])
    [(u'C1', u'6', 0.5337, 0.6968, 0.0), (u'C2', u'6', 0.5337, -0.6968, 0.0), (u'C3', u'6', -0.6618, -1.4038, 0.0), (u'C4', u'6', -1.8695, -0.6988, 0.0), (u'C5', u'6', -1.8695, 0.6988, 0.0), (u'C6', u'6', -0.6618, 1.4037, 0.0), (u'F1', u'9', 1.7137, 1.3554, 0.0), (u'F2', u'9', 1.7137, -1.3554, 0.0)]

    >>> for i in db[8]:
    ...   print(i)
    (u'P1', u'15', 0.0, 0.0, -0.0001)
    (u'F1', u'9', -1.5721, 0.4748, 0.1103)
    (u'F2', u'9', -0.1279, -0.042, -1.6402)
    (u'F3', u'9', 1.5722, -0.4749, -0.1103)
    (u'F4', u'9', 0.1278, 0.042, 1.6403)
    (u'F5', u'9', 0.4704, 1.5755, -0.077)
    (u'F6', u'9', -0.4705, -1.5754, 0.077)

    >>> print(db[-1])
    Only positive index numbers allowed!
    False

    :param fragment_id: Id number of fragment to return.
    :type fragment_id: int
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 0:
      print('Only positive index numbers allowed!')
      return False
    elif fragment_id == 0:
      return False
      # fragment_id = len(self)-abs(fragment_id)
    found = self._get_fragment(fragment_id)
    if found:
      return found
    else:
      raise IndexError('Database fragment not found.')

  def __delitem__(self, fragment_id):
    """
    Called to implement deletion of self[fragment_id].

    # have to create a copy of the db before I delete the entry:
    >>> import shutil
    >>> shutil.copyfile('tests/tst.sqlite', 'tests/tst1.sqlite')
    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> del db[2]
    >>> db[2]
    [(u'N1', u'7', 5.2916, 0.1111, 10.84), (u'N2', u'7', 4.1672, 1.8088, 9.1734), (u'C1', u'6', 5.8369, -0.7123, 11.7478), (u'C2', u'6', 5.9368, -2.0826, 11.5626), (u'C3', u'6', 5.4664, -2.6325, 10.3782), (u'C4', u'6', 4.9264, -1.7885, 9.4165), (u'C5', u'6', 4.8535, -0.4237, 9.6802), (u'C6', u'6', 4.2675, 0.5478, 8.7222), (u'C7', u'6', 3.8265, 0.1802, 7.4491), (u'C8', u'6', 3.2516, 1.1452, 6.6304), (u'C9', u'6', 3.1267, 2.4421, 7.1103), (u'C10', u'6', 3.5991, 2.7331, 8.3809)]
    >>> print('before: ' + str(len(db)))
    before: 71

    # del db[0] deletes nothing:
    >>> del db[0]
    >>> print('after: ' + str(len(db)))
    after: 71
    >>> del db[3]
    >>> print('after: ' + str(len(db)))
    after: 71
    >>> print(db[-3][1])
    Traceback (most recent call last):
    ...
    TypeError: 'bool' object has no attribute '__getitem__'

    :param fragment_id: Id number of fragment to delete.
    :type fragment_id: int
    """
    # in case of negative id, get a list of the ids and access the id through
    # that list:
    if fragment_id < 0:
      fragment_id = self.get_all_rowids()[fragment_id]
    # req = '''DELETE FROM Fragment WHERE rowid = ?'''
    req_usr = '''DELETE FROM userdb.Fragment WHERE rowid = ?'''
    try:
      fragment_id = int(fragment_id)
    except(ValueError, TypeError):
      print('Wrong type. Expected integer.')
    # actually delete the item:
    if fragment_id < 1000000:
      # print('can not delete fragment of main database.')
      return False
      # deleted = self.database.db_request(req, fragment_id)
    else:
      fragment_id = fragment_id - 1000000
      deleted = self.database.db_request(req_usr, fragment_id)
    return deleted

  def __iter__(self):
    """
    This method is called when an iterator is required for FragmentTable.
    Returns the Id and the Name as tuple.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> for num, i in enumerate(db):
    ...   print(i)
    ...   if num > 1:
    ...     break
    [4, u'Acetate anion, C2H3O2-']
    [1000004, u'Acetate anion, C2H3O2-  *user*']
    [47, u'Acetone, C3H6O']
    """
    all_fragments = self.get_all_fragment_names()
    if all_fragments:
      return iter(all_fragments)
    else:
      return False

  def has_name(self, name):
    """
    Returns True if a partial name is found in the DB.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.has_name('Benzene')
    'userdb'
    >>> db.has_name('Benzil')
    False
    """
    req = '''SELECT Name FROM Fragment WHERE Name like "%{}%" '''.format(name)
    req_usr = '''SELECT Name FROM userdb.Fragment WHERE Name like "%{}%" '''.format(name)
    # user modifications have highest priority:
    if self.database.db_request(req_usr):
      return 'userdb'
    elif self.database.db_request(req):
      return True
    else:
      return False

  def has_exact_name(self, name):
    """
    Returns True if an exact name is found in the DB. Returns 'userdb if
    the name is found in the userdb.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.has_exact_name('1,2-Dichlorobenzene, C6H4Cl2')
    True

    >>> db.has_exact_name('Benzene')
    False
    """
    req = '''SELECT Name FROM Fragment WHERE Fragment.Name = "{}" '''.format(name)
    req_usr = '''SELECT Name FROM userdb.Fragment WHERE Fragment.Name = "{}" '''.format(name)
    # user modifications have highest priority:
    if self.database.db_request(req_usr):
      return 'userdb'
    elif self.database.db_request(req):
      return True
    else:
      return False

  def has_exact_resi_class(self, resi_class):
    """
    Returns True if an exact residue class is found in the DB.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.has_exact_resi_class('BENZ')
    True
    >>> db.has_exact_resi_class('Benzene')
    False
    """
    req = '''SELECT class FROM Fragment WHERE Fragment.class = "{}" '''.format(resi_class)
    req_usr = '''SELECT class FROM userdb.Fragment WHERE Fragment.class = "{}" '''.format(resi_class)
    # user modifications have highest priority:
    if self.database.db_request(req_usr):
      return True
    elif self.database.db_request(req):
      return True
    else:
      return False

  def has_index(self, fragment_id):
    """
    Returns True if db has index Id
    :param Id: Id of the respective fragment
    :type Id: int
    :rtype: bool

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> db.has_index('5')
    True
    >>> db.has_index('999')
    False
    """
    try:
      fragment_id = int(fragment_id)
    except TypeError:
      print('Wrong type for database index. Only numbers allowed.')
      return False
    req = '''SELECT Id FROM Fragment WHERE Fragment.Id = ?'''
    req_usr = '''SELECT Id FROM userdb.Fragment WHERE Id = ?'''
    if fragment_id < 1000000:
      rows = self.database.db_request(req, fragment_id)
    else:
      fragment_id = fragment_id - 1000000
      rows = self.database.db_request(req_usr, fragment_id)
    if rows:
      return True
    else:
      return False

  def get_all_rowids(self):
    """
    returns all Ids in the database as list.
    :rtype: list

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.get_all_rowids()
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 1000001, 1000002, 1000003, 1000004, 1000005, 1000006]
    """
    req = '''SELECT Id FROM Fragment ORDER BY Id'''
    req_usr = '''SELECT Id FROM userdb.Fragment ORDER BY Id'''
    rows = [i[0] for i in self.database.db_request(req)]
    if not rows:
      return False
    if self.userdb:
      rows_usr = self.database.db_request(req_usr)
      if rows_usr:
        rows_usr = [i[0] + 1000000 for i in rows_usr]
        rows = rows + rows_usr
    return rows

  def get_all_fragment_names(self):
    """
    returns all fragment names in the database, sorted by name
    """
    req = '''SELECT Id, name FROM Fragment'''
    req_usr = '''SELECT Id, name FROM userdb.Fragment'''
    rows = [list(i) for i in self.database.db_request(req)]
    if self.userdb:
      rows_usr = self.database.db_request(req_usr)
      if rows_usr:
        rows_usr = [[i[0] + 1000000, i[1] + '  *user*'] for i in rows_usr]
        rows.extend(rows_usr)
    for num, i in enumerate(rows):
      # searchkey also adds sum formula etc to sortkey:
      key = make_sortkey(i[1], searchkey=False)
      rows[num].append(key)
    if rows:
      rows.sort(key=lambda x: (x[2][0], x[2][1]))
      for i in rows:
        del i[-1]
      return rows
    else:
      return False

  def fragid_toint(self, fragment_id):
    try:
      return int(fragment_id)
    except ValueError as e:
      return 0

  def _get_fragment(self, fragment_id):
    """
    returns a full fragment with all atoms, atom types as a tuple.
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    :rtype: tuple
    """
    fragment_id = self.fragid_toint(fragment_id)
    req_atoms = ("""SELECT Atoms.name, Atoms.element, Atoms.x, Atoms.y, Atoms.z
                       FROM Fragment, Atoms ON Fragment.Id=Atoms.FragmentId WHERE
                       Fragment.Id = ?""")
    req_atoms_usr = '''SELECT userdb.Atoms.name, userdb.Atoms.element, userdb.Atoms.x, 
            userdb.Atoms.y, userdb.Atoms.z FROM userdb.Fragment, userdb.Atoms ON 
            userdb.Fragment.Id=userdb.Atoms.FragmentId WHERE Fragment.Id = ?'''
    if fragment_id < 1000000:
      atomrows = self.database.db_request(req_atoms, fragment_id)
    else:
      fragment_id = fragment_id - 1000000
      atomrows = self.database.db_request(req_atoms_usr, fragment_id)
    return atomrows

  def get_fragment_name(self, fragment_id):
    """
    returns the "Name" column entry of fragment with id "fragment_id"
    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> db.get_fragment_name(8)
    u'Hexafluorophosphate, PF6'

    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_name = '''SELECT Fragment.Name FROM Fragment WHERE Fragment.Id = ? '''
    else:
      fragment_id = fragment_id - 1000000
      req_name = '''SELECT userdb.Fragment.Name FROM userdb.Fragment WHERE Fragment.Id = ? '''
    name = self.database.db_request(req_name, fragment_id)[0]
    return name[0]

  def get_picture(self, fragment_id):
    """
    returns a picture of the fragment if one exist in the database. Otherwise
    it returns False.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> pic = db.get_picture(2)
    >>> print(str(pic[1:4]))
    PNG
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_picture = '''SELECT picture FROM Fragment WHERE Id = ? '''
    else:
      fragment_id = fragment_id - 1000000
      req_picture = '''SELECT picture FROM userdb.Fragment WHERE Id = ? '''
    try:
      picture = self.database.db_request(req_picture, fragment_id)[0][0]
    except TypeError:
      return None
    return picture

  def get_residue_class(self, fragment_id):
    """
    returns the "class" column entry of fragment with id "fragment_id"

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> db.get_residue_class(47)
    'ACE'

    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_class = '''SELECT class FROM Fragment WHERE Id = ?'''
    else:
      fragment_id = fragment_id - 1000000
      req_class = '''SELECT class FROM userdb.Fragment WHERE Id = ?'''
    classname = self.database.db_request(req_class, fragment_id)
    try:
      classname = classname[0][0]
    except(IndexError, TypeError):
      print('Could not find residue class.')
      return ''
    return str(classname)

  def get_restraints(self, fragment_id):
    """
    returns the restraints for Fragment(Id) from the database.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile)
    >>> db.get_restraints(5)
    [(u'SADI', u'0.02 C1 F1 C2 F2'), (u'SADI', u'0.02 C1 C6 C5 C6 C1 C2 C4 C5 C3 C4 C2 C3'), (u'SADI', u'0.04 C1 F2 C3 F2 C2 F1 C6 F1'), (u'SADI', u'0.04 C4 C6 C3 C5 C2 C4 C1 C3 C2 C6 C1 C5'), (u'FLAT', u'C1 > F2'), (u'SIMU', u'C1 > F2'), (u'RIGU', u'C1 > F2')]

    :param fragment_Id: id of the fragment in the database
    :type fragment_Id: int
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms FROM Restraints WHERE FragmentId = ?'''
    else:
      fragment_id = fragment_id - 1000000
      req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms FROM userdb.Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return restraintrows

  def get_reference(self, fragment_id):
    """
    returns the reference for Fragment(Id) from the database.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.get_reference(4)
    u'CCDC LIBXUR'
    >>> db.get_reference(999)
    'no reference found'

    :param fragment_Id: id of the fragment in the database
    :type fragment_Id: int
    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_ref = '''SELECT Reference FROM Fragment WHERE Fragment.Id = ? '''
    else:
      fragment_id = fragment_id - 1000000
      req_ref = '''SELECT Reference FROM userdb.Fragment WHERE Fragment.Id = ? '''
    rows = self.database.db_request(req_ref, fragment_id)
    try:
      ref = rows[0][0]
      return ref
    except(TypeError):
      return 'no reference found'

  def find_fragment_by_name(self, name, selection=5):
    """
    :type selection: int
    :param selection: return only this number of hits

    find a fragment by its name in the database. This method will output a
    selection of (default=5) best hits.

    >>> dbfile = 'tests/tst1.sqlite'
    >>> db = FragmentTable(dbfile, 'tests/tst-usr.sqlite')
    >>> db.find_fragment_by_name('nona', selection=3)
    [[50, u'n-Nonane, C8H18', [0.571429, '']], [34, u'Nitrate anion, [NO3]-', [0.8, '']], [58, u'Nonafluoro-tert-butoxy, [(CF3)3CO]-', [0.8, '']]]

    :param name: (part of) the name of a fragment to find
    :type name: str
    :return fragment_id: list of id numbers of the found fragments e.g. [1, 5]
    :type fragment_id: int
    """
    return self._search_name(name, selection)

  def _search_name(self, search_string, selection=5):
    """
    searches the names in the database for a given name
    :param search_string: search for this string
    :type search_string: str
    :param selection: return this amount of results
    :type selection: int
    """
    search_results = []
    for i in self:
      key = make_sortkey(i[1], searchkey=True)
      coefficient = dice_coefficient2(search_string, key[0] + key[1])
      i.append([coefficient, key[1]])
      search_results.append(i)
    # select the best n results:
    selected_results = sorted(search_results, key=lambda coeff: (coeff[-1][0], coeff[-1][1]), reverse=False)[:selection]
    return selected_results

  def store_fragment(self, fragment_name=None, atoms=None, resiclass=None, restraints=None,
                     reference=None, comment=None, picture=None):
    """
    Store a complete new fragment into the database. Minimal requirement is a
    fragment name (Full chemical name) and a list of atoms. Restraints,
    reference and comments are optional.

    :param fragment_name: full chemical name of the fragment
    :type fragment_name: string
    :param atoms: atoms of the fragment ['name', 'atomic number', 'x', 'y', 'z']
    :type atoms: list
    :param restraints: [['DFIX', '1.564', 'C1 C2'], ['SADI', 'C2 C3 C4 C5']]
    :type restraints: list of list
    :param reference: optional short description where it came from
    :type reference: string
    :param comment: optional any comment about the fragment
    :type comment: string
    :param picture: a picture of the molecule
    :type picture: binary
    :rtype int: FragmentId -> last_rowid
    """
    # first stores the meta-information in the Fragment table:
    # The FragmentId is the last_rowid from sqlite
    if not fragment_name or fragment_name == '':
      return None
    fragmentid = self._fill_fragment_table(fragment_name, resiclass,
                                           reference, comment, picture)
    if not fragmentid:
      raise Exception('No Id obtained during fragment storage.')
    # then stores atoms with the previously obtained FragmentId
    self._fill_atom_table(fragmentid, atoms)
    # in case of supplied restraints store them also:
    if restraints:
      self._fill_restraint_table(fragmentid, restraints)
    return fragmentid

  def _fill_fragment_table(self, fragment_name, resiclass=None,
                           reference=None, comment=None, picture=None):
    """
    Fills a fragment into the database.
    :param fragment_name: Nam eof the Fragment
    :type fragment_name: str
    :param comment: any comment about the fragment
    :type comment: str
    :rtype list: last_rowid
    """
    if picture:
      picture = sqlite3.Binary(picture)
    table = (fragment_name, resiclass, reference, comment, picture)
    req = '''INSERT INTO userdb.Fragment (name, class, reference, comment, picture)  VALUES(?, ?, ?, ?, ?)'''
    fragid = self.database.db_request(req, table)
    return fragid + 1000000

  def _fill_atom_table(self, fragment_id, atom_table):
    """
    Fills atoms into the Atoms table.
    [('C1', 1.2, -0.023, 3.615), ('C2', 1.203, -0.012, 2.106), ...]
    or
    [('C1 1.2 -0.023 3.615'), ('C2 1.203 -0.012 2.106'), ...]
    :param fragment_id: Id of the respective Fragment
    :type fragment_id: int or str
    :param atom_table: list of lits or list of strings
    :type atom_table: list
    """
    # test wether atom_table is a list or list of list, because we want no string
    # in a list here.
    fragment_id = int(fragment_id)
    if fragment_id > 1000000:
      fragment_id = fragment_id - 1000000
    if not atom_table or not fragment_id:
      print('No atoms supplied! Doing nothing')
      return
    if not isinstance(atom_table[0], (list, tuple)):
      if isinstance(atom_table[0], str):
        atom_table = [i.split() for i in atom_table]
      else:
        raise Exception('wrong data type "{}" for atom list.'.format(type(atom_table[0])))
    for line in atom_table:
      Name = line[0]
      element = 999  # line[1]
      x = line[1]
      y = line[2]
      z = line[3]
      req = '''INSERT INTO userdb.atoms (FragmentId, Name, element, x, y, z) VALUES(?, ?, ?, ?, ?, ?)'''
      self.database.db_request(req, (fragment_id, Name, element, x, y, z))

  def _fill_restraint_table(self, fragment_id, restraints_list):
    """
    Fills the restraints table with restraints. restraints_list must be a list
    or tuple of string lists like:
    [['SADI C1 F1 C2 F2'],
    ['SADI 0.04 F2 C3 F1 C6 F2 C1 F1 C2'],
    ['SAME C2 > C6 C1']]
    :param FragmentId: Fragment Id where the restraints belong to.
    :type FragmentId: int or str
    :param restraints_list: list with restraints
    :type restraints_list:
    """
    # test if restraint_list is a list of strings. we dont want list of list here.
    fragment_id = int(fragment_id)
    if fragment_id > 1000000:
      fragment_id = fragment_id - 1000000
    try:
      restraints_list[0]  # is there even one restraint?
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
        restr_table.append(str(fragment_id))
        restr_table.append(line[:4])
        restr_table.append(line[5:])
        req = '''INSERT INTO userdb.Restraints (FragmentId, ShelxName, atoms) VALUES(?, ?, ?)'''
        self.database.db_request(req, restr_table)


class Restraints():
  def __init__(self, dbfile, userdb):
    self.database = DatabaseRequest(dbfile, userdb)

  def fragid_toint(self, fragment_id):
    try:
      int(fragment_id)
    except ValueError as e:
      print(e)
      return
    return int(fragment_id)

  def get_restraints_from_fragmentId(self, fragment_id):
    """
    returns the restraints from a database entry

    >>> dbfile = 'tests/tst1.sqlite'
    >>> res = Restraints(dbfile, 'tests/tst-usr.sqlite')
    >>> res.get_restraints_from_fragmentId(7)
    [(u'DFIX', u'1.783 C1 CL1 C1 CL2 C1 CL3'), (u'DANG', u'2.946 CL1 CL2 CL2 CL3 CL3 CL1'), (u'RIGU', u'C1 > CL3'), (u'SIMU', u'C1 > CL3')]

    >>> res.get_restraints_from_fragmentId(1000001)
    [(u'SADI', u'0.02 C1 C2 C2 C3 C3 C4 C4 C5 C5 C6 C6 C1'), (u'SADI', u'0.04 C1 C5 C1 C5 C4 C2 C4 C6 C6 C2 C5 C3'), (u'FLAT', u'C1 > C6'), (u'SIMU', u'C1 > C6'), (u'RIGU', u'C1 > C6')]

    """
    fragment_id = self.fragid_toint(fragment_id)
    if fragment_id < 1000000:
      req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms FROM Restraints WHERE FragmentId = ?'''
    else:
      fragment_id = fragment_id - 1000000
      req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms FROM userdb.Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return restraintrows


if __name__ == '__main__':
  import os
  import doctest

  failed, attempted = doctest.testmod()  # verbose=True)
  if failed == 0:
    print('passed all {} tests!'.format(attempted))

  # import cProfile
  dbname = 'tst.sqlite'
  dbfile = os.path.abspath('tests/{}').format(dbname)
  print(dbfile)

  # dbfile = 'tst.sqlite'
  #  call_profile(dbfile)
  # db = FragmentTable(dbfile, userdb='tests/tst-usr.sqlite')
  # picture = db.get_picture(2)
  # print(db.get_residue_class(1000005))
  # names = db.get_all_fragment_names()
  # print(db[1000005])
  # for i in names:
  #  print(i)

  sys.exit()
  atoms = [[u'C1', u'6', 1.2, -0.023, 3.615], (u'C2', u'6', 1.203, -0.012, 2.106), (u'C3', u'6', 0.015, -0.011, 1.39),
           (u'C4', u'6', 0.015, -0.001, 0.005), (u'C5', u'6', 1.208, 0.008, -0.688), (u'C6', u'6', 2.398, 0.006, 0.009),
           (u'C7', u'6', 2.394, -0.004, 1.394)]
  # atoms = ['C1 6 1.2 -0.023 3.615', 'C2 6 1.203 -0.012 2.106', 'C3 6 0.015 -0.011 1.39', 'C4 6 0.015 -0.001 0.005', 'C5 6 1.208 0.008 -0.688', 'C6 6 2.398 0.006 0.009', 'C7 6 2.394 -0.004 1.394']
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
  fragment_name = table[1]
  comment = 'asfgagr'
  resiclass = 'bnzr'

  # database = DatabaseRequest(dbfile)
  # database.table_info('userdb.')
  # database.table_info('')
  # database.values_in_col('Fragment')
  # idf = db.has_index('1')
  # print(idf)
  # print(db[1])
  # print(db.has_name('Benzene'))

  # fid = db.store_fragment(fragment_name, atoms, resiclass, restraints, reference, comment, picture=False)
  # if fid:
  #   print('stored', fid)
  # print(db[5])

  dbfile = "fragment-database.sqlite"


  def call_profile(dbfile):
    import cProfile
    import pstats
    cp = cProfile.Profile()
    db = FragmentTable(dbfile)
    cp.runcall(db._get_fragment, 17)
    pstats.Stats(cp).strip_dirs().sort_stats('time').print_stats(20)


  call_profile(dbfile)
