'''
Created on 09.10.2014

@author: Daniel Kratzert
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


class DatabaseRequest():
  def __init__(self, dbfile):
    '''
    creates a connection to the SQLite3 database file "dbfile".
    :param dbfile: database file
    :type dbfile: str
    '''
    # open the database
    self.con = sqlite3.connect(dbfile)
    self.con.text_factory = str
    with self.con:
      # set the database cursor
      self.cur = self.con.cursor()

  def fetch_all(self):
    '''fetch all rows'''
    rows = self.cur.fetchall()
    return rows
  
  def db_request(self, request, no_fetch=False, *args):
    '''
    Performs a SQLite3 database request with "request" and optional arguments
    to insert parameters via "?" into the database request.
    :param request:
    :type request:
    '''
    # make the request
    try:
      self.cur.execute(request, args)
    except OperationalError as e:
      print(e)
      return False
    if not no_fetch:
      return self.fetch_all()

class FragmentTable():
  def __init__(self, dbfile):
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
      self.no_db()
    all_names = (i[1].lower() for i in all_fragments)
    found = False
    if isinstance(name, (str, unicode)):
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
    return len(rows)
  
  def __getitem__(self, fragment_id):
    '''
    Called to implement evaluation of self[fragment_id]. 
    print FragmentTable[fragment_id]
    :param fragment_id: Id number of fragment to return.
    :type fragment_id: int
    '''
    if isinstance(fragment_id, int):
      found = self.get_fragment(fragment_id)
      if found:
        return found 
      else:
        raise IndexError
    if isinstance(fragment_id, (str, unicode)):
      found = self.get_fragment(fragment_id)
      if found:
        return found 
      else:
        raise IndexError
  
  def __delitem__(self, fragment_id):
    '''
    Called to implement deletion of self[fragment_id].
    del FragmentTable[3]
    :param fragment_id: Id number of fragment to delete.
    :type fragment_id: int
    '''
    req = '''DELETE FROM fragment WHERE fragment.id = ?'''
    deleted = self.database.db_request(req, fragment_id)
    return deleted
  
  def __iter__(self):
    '''
    This method is called when an iterator is required for FragmentTable. 
    for i in db:
      print(i)
    '''
    all_fragments = self.get_all_fragment_names()
    return (i[1] for i in all_fragments)

  def no_db(self):
    print('No database items found!')
    sys.exit()
  
  def get_all_fragment_names(self):
    '''
    returns all fragment names in the database, sorted by name
    '''
    req = '''SELECT fragment.id, fragment.name FROM fragment ORDER BY name'''
    rows = self.database.db_request(req)
    return rows

  def get_fragment(self, fragment_id):
    '''
    returns a full fragment with all atoms, atom types as a dictionary
    :param fragment_id: id of the fragment in the database
    :type fragment_id: int
    '''
    req_atoms = '''SELECT atoms.name, atoms.element, atoms.x, atoms.y, atoms.z
      FROM fragment, atoms on fragment.id=atoms.fragmentid WHERE
      fragment.id = {}'''.format(fragment_id)
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
    restraintrows = self.db_request(req_restr, fragment_id)
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
    frags = self.db_request(req)
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
    Store a new fragment into the database.
    
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
    fragment_id = 999
    return fragment_id


class Restraints(DatabaseRequest):
  def __init__(self, dbfile):
    self.database = DatabaseRequest(dbfile)

  def get_restraints_from_fragmentId(self, fragment_id):
    req_restr = '''SELECT Restraints.ShelxName, Restraints.Atoms
      FROM Restraints WHERE FragmentId = ?'''
    restraintrows = self.database.db_request(req_restr, fragment_id)
    return(restraintrows)



if __name__ == '__main__':
  dbfile = 'F:\GitHub\DSR-db\dk-database_2.sqlite'
  db = FragmentTable(dbfile)
  
  print(db[5])
  def match_dbfrag(fragId=17):
    for i in db[fragId]:
      print(i)
  
  match_dbfrag(6)
  
  #tst = {'name': 'hallofrag', 'atoms': [['C1', '6', '0.123', '1.324', '0.345'], 
  #                                  ['C1', '6', '0.6123', '1.3624', '0.3645']]}
  #print(tst['name'])
  #db[99, tst]
  
  #for i in db:
  #  print(i)
  
 # print('len', len(db))
 # if 'benzene' in db:
 #   print('yes')

  
  
  #import cProfile
  #cProfile.run("allnames()", "foo.profile")
  #   res = Restraints(dbfile)
  #  for r in res.get_restraints_from_fragmentId(15):
  #    pass
    #print(r)

  #   print(hasattr(db, '__iter__'))
  #   if 'OC(CF3)3' in db:
  #     print('yes')
  #   else:
  #     print('no benz')

  # for i in db:
  #   print(i)

  #print('##################')
  #found = db.find_fragment_by_name('CF3', selection=3)
  #print(found)
  #for i in db.get_fragment(fragment_id=2):
  #  print(i)
