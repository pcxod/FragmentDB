'''
Created on 03.05.2015

@author: daniel

here are only functions that are completely independent of olex
'''
from collections import Counter

SHX_CARDS = ('TITL', 'CELL', 'ZERR', 'LATT', 'SYMM', 'SFAC', 'UNIT', 'LIST',
             'L.S.', 'CGLS', 'BOND', 'FMAP', 'PLAN', 'TEMP', 'ACTA', 'CONF',
             'SIMU', 'RIGU', 'WGHT', 'FVAR', 'DELU', 'SAME', 'DISP', 'LAUE',
             'REM' , 'MORE', 'TIME',  'END', 'HKLF', 'OMIT', 'SHEL', 'BASF',
             'TWIN', 'EXTI', 'SWAT', 'HOPE', 'MERG', 'SPEC', 'RESI', 'MOVE',
             'ANIS', 'AFIX', 'HFIX', 'FRAG', 'FEND', 'EXYZ', 'EADP', 'EQIV',
             'CONN', 'BIND', 'FREE', 'DFIX', 'BUMP', 'SADI', 'CHIV', 'FLAT',
             'DEFS', 'ISOR', 'NCSY', 'SUMP', 'BLOC', 'DAMP', 'STIR', 'MPLA',
             'RTAB', 'HTAB', 'SIZE', 'WPDB', 'GRID', 'MOLE', 'XNPD', 'REST',
             'CHAN', 'FLAP', 'RNUM', 'SOCC', 'PRIG', 'WIGL', 'RANG', 'TANG',
             'ADDA', 'STAG', 'NEUT', 'ABIN', 'ANSC', 'ANSR', 'NOTR', 'TWST',
             'PART', 'DANG')

RESTRAINT_CARDS = ('SIMU', 'RIGU', 'DELU', 'SAME', 'FREE', 'DFIX', 'BUMP', 'HFIX', 
                   'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'DANG')

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

def check_restraints_consistency(restraints, atoms, fragment_name):
  '''
  - Checks if the Atomnames in the restraints of the dbhead are also in
    the list of the atoms of the respective dbentry.
  - Checks wether restraints cards are vaid.
  - checks for duplicated atoms in the atoms list
  :param restraints: list of strings like ['SADI C1 C2', '', ...]
  :type restraints: list of strings
  :param atoms: list of atoms in the fragment
  :type atoms: list
  :param fragment_name: fragment name
  :type fragment_name: string
  '''
  if not restraints:
    print('No restraints found!')
    return True
  status = True
  atoms = [i[0].upper() for i in atoms]
  # check for duplicates:
  if len(set(atoms)) != len(atoms):
            c1 = Counter(atoms)
            c2 = Counter(set(atoms))
            diff = c1 - c2
            duplicates = list(diff.elements())
            for i in duplicates:
                print('\nDuplicate atom "{}" found!\n'.format(duplicates.pop()))
                status = False
  # check if restraint cards are valid
  restraint_atoms_list = set([])
  for n, line in enumerate(restraints):
    if not line:
      continue
    line = line.upper()
    line2 = line.split()
    # only the first 4 characters, because SADI_TOL would be bad:
    if line2[0] not in SHX_CARDS:  
      status = False
      print('Invalid line in header of database entry "{}" found!'.format(n, fragment_name))
      print(line)
    if line[:4] in RESTRAINT_CARDS:
      line = line[5:].split()
      for i in line:
        if i in ('>', '<'):
          continue
        try:
          float(i)
        except(ValueError):
          restraint_atoms_list.add(i)
  # check if restrained atoms are in th eatom list:
  for atom in restraint_atoms_list:
    atom = atom.upper()
    if not atom in atoms:
      status = False
      print('Unknown atom "{}" in restraints of "{}".'.format(atom, fragment_name))
  if not status:
    print('Check database entry.\n')
  return status




def initialize_user_db(user_dbpath):
  '''
  initializes an empty user database in DataDir()/db/
  :param user_dbpath: path with filename where to create the database.
  :type user_dbpath: string
  '''
  import sqlite3 as lite
  con = lite.connect(user_dbpath)
  con.text_factory = str
  cur = con.cursor()
  print('initializing FragmentDB user databse.')
  con.execute("PRAGMA foreign_keys = ON")
  cur.execute("DROP TABLE IF EXISTS fragment")
  cur.execute("DROP TABLE IF EXISTS atoms")
  cur.execute("DROP TABLE IF EXISTS atom")
  cur.execute("DROP TABLE IF EXISTS FragmentRestraints")
  cur.execute("DROP TABLE IF EXISTS Restraints")
  try:
    cur.execute("DROP INDEX Atoms_FK")
  except:
    pass
  try:
    cur.execute("DROP INDEX Restraint_FK")
  except:
    pass
  try:
    cur.execute("DROP INDEX Fragment_Name")
  except:
    pass
  try:
    cur.execute("DROP INDEX AtomId")
  except:
    pass
  cur.execute('''
              CREATE TABLE Fragment (
                  Id    INTEGER NOT NULL,
                  class  VARCHAR(4),
                  version TEXT,
                  Name    TEXT,
                  Reference    TEXT,
                  comment    TEXT,
                  picture    BLOB,
                  PRIMARY KEY(Id));
              ''')
  cur.execute('''
              CREATE TABLE Atoms (
                  Id    INTEGER NOT NULL,
                  FragmentId    INTEGER NOT NULL,
                  version TEXT,
                  Name    VARCHAR(255),
                  element    VARCHAR(2),
                  x    FLOAT,
                  y    FLOAT,
                  z    FLOAT,
              PRIMARY KEY(Id),
                FOREIGN KEY(FragmentId)
                  REFERENCES Fragment(Id)
                    ON DELETE CASCADE
                    ON UPDATE NO ACTION);
              ''')
  cur.execute('''
                 CREATE TABLE Restraints (
                    Id INTEGER  NOT NULL,
                    FragmentId  INTEGER NOT NULL,
                    version TEXT,
                    ShelxName CHAR(4),
                    Atoms TEXT,
                  PRIMARY KEY(Id),
                    FOREIGN KEY(FragmentId)
                      REFERENCES Fragment(Id)
                      ON DELETE CASCADE
                      ON UPDATE NO ACTION);
                  ''')
  cur.execute('''
              CREATE INDEX Atoms_FK ON Atoms(FragmentId);
              ''')
  cur.execute('''
              CREATE INDEX Restraint_FK ON Restraints(FragmentId);
              ''')
  cur.execute('''
              CREATE INDEX Fragment_Name ON Atoms(Name);
              ''')
  cur.execute('''
              CREATE INDEX AtomId ON Atoms(Id);
              ''')
  con.execute("PRAGMA foreign_keys = ON")
  con.commit()



if __name__ == '__main__':
  pass