'''
Created on 03.05.2015

@author: daniel

here are only functions that are completely independent of olex
'''
from __future__ import print_function
from collections import Counter
from math import radians, cos, sqrt, sin
import string
from copy import deepcopy
#import ast

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
             'PART', 'DANG', 'TRIA')

RESTRAINT_CARDS = ('SIMU', 'RIGU', 'DELU', 'SAME', 'FREE', 'DFIX', 'BUMP', 'HFIX', 
                   'SADI', 'CHIV', 'FLAT', 'DEFS', 'ISOR', 'NCSY', 'DANG')

# these are supported by -i option in olex2:
IMPL_RESTRAINT_CARDS = ('SAME', 'SADI', 'DFIX', 'BUMP', 'DANG', 'FLAT', 'TRIA',
                        'CHIV', 'DELU', 'SIMU', 'ISOR') 

ABS_RESTR_CARDS = ('DFIX', 'DANG', 'BUMP', 'TRIA', 'CHIV')
REL_RESTR_CARDS = ('SAME', 'SADI', 'SIMU', 'RIGU', 'ISOR', 'NCSY', 'FLAT')

def atomic_distance(p1, p2, cell):
    '''
    p1 and p2 are x, y , z coordinates as list ['x', 'y', 'z']
    cell are the cell parameters as list: ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    returns the distance between the two points.
    
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
    >>> coord1 = (-0.186843,   0.282708,   0.526803) 
    >>> coord2 = (-0.155278,   0.264593,   0.600644) 
    >>> atomic_distance(coord1, coord2, cell)
    1.5729229943265979
    '''
    cell = [float(y) for y in cell]
    a , b, c =  cell[:3]
    al = radians(cell[3])
    be = radians(cell[4])
    ga = radians(cell[5])
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    dx = (x1-x2)
    dy = (y1-y2)
    dz = (z1-z2)
    dsq = (a*dx)**2+(b*dy)**2+\
          (c*dz)**2+2*b*c*cos(al)*dy*dz+\
          2*dx*dz*a*c*cos(be)+2*dx*dy*a*b*cos(ga)
    return(sqrt(dsq))

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
      status = 'False'
  # check if restraint cards are valid
  restraint_atoms_list = set([])
  for line in restraints:
    if not line:
      continue
    line = [i.upper() for i in line]
    # only the first 4 characters, because SADI_TOL would be bad:
    if line[0][:4] not in SHX_CARDS:  
      status = False
      print('\nInvalid line in restraints list of "{}" found!'.format(fragment_name))
      print(' '.join(line))
    if line[0][:4].upper() in RESTRAINT_CARDS:
      for i in line[1:]:
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
      print('\nUnknown atom "{}" in restraints of "{}".'.format(atom, fragment_name))
      print('Please remove non-existent atoms from the restraint list!')
  if not status:
    print('Check database entry.\n')
  return status

def pairwise(iterable):
    '''
     s -> (s0,s1), (s2,s3), (s4, s5), ...
     
     >>> liste = ['C1', 'C2', 'C2', 'C3', 'C4', 'C5', 'C5', 'C6']
     >>> pairwise(liste)
     [('C1', 'C2'), ('C2', 'C3'), ('C4', 'C5'), ('C5', 'C6')]
     '''
    a = iter(iterable)
    return zip(a, a)

def mean(values):
    '''
    returns mean value of a list of numbers
    
    >>> mean([1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> round(mean([1, 2, 3, 4, 1, 2, 3, 4.1, 1000000]), 4)
    111113.3444
    '''
    mean = sum(values) / float(len(values)) 
    return mean

def median(nums):
    """
    calculates the median of a list of numbers
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4])
    2.5
    >>> median([1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4.1, 1000000])
    3
    >>> median([])
    Traceback (most recent call last):
    ...
    ValueError: Need a non-empty iterable
    """
    ls = sorted(nums)
    n = len(ls)
    if n == 0:
        raise ValueError("Need a non-empty iterable")
    # for uneven list length:
    elif n % 2 == 1:
        # // is floordiv:
        return ls[n // 2]
    else:
        return sum(ls[(n / 2 - 1):(n / 2 + 1)]) / 2.0

def std_dev(data):
    '''
    returns standard deviation of values rounded to pl decimal places
    S = sqrt( (sum(x-xm)^2) / n-1 )
    xm = sum(x)/n
    :param values: list with integer or float values
    :type values: list  
    :param pl: round to n places
    :type pl: integer
    
    >>> round(std_dev([1.234, 1.222, 1.345, 1.451, 1.000, 1.234, 1.321, 1.222]), 8)
    0.1303522
    '''
    if len(data) == 0:
        return 0
    K = data[0]
    n = 0
    Sum = 0
    Sum_sqr = 0
    for x in data:
        n = n + 1
        Sum += x - K
        Sum_sqr += (x - K) * (x - K)
    variance = (Sum_sqr - (Sum * Sum) / n) / (n - 1)
    # use n instead of (n-1) if want to compute the exact variance of the given data
    # use (n-1) if data are samples of a larger population
    return sqrt(variance)

def check_sadi_consistence(atoms, restr, cell, fragment):
  '''
  check if same distance restraints make sense. Each length of an atom
  pair is tested agains the deviation from the mean of each restraint.
  The deviation must ly in factor times the rmsd.        
  :param atoms: atoms list of the fragment ['C1', x, y, z]
  :param restraints: restraints list
  :param fragment: frag name
  :param factor: factor for confidence interval
  '''
  restraints = deepcopy(restr)
  atnames = [i[0].upper() for i in atoms]
  for num, line in enumerate(restraints):
    if not line:
      continue
    if line[0].upper() == 'SADI':
      del line[0]
      try:
        if not str(line[0][0]).isalpha():
          del line[0] # delete standard deviation
      except(IndexError):
        return
      if len(line)%2 == 1: # test for uneven atoms count
        print('Inconsistent SADI restraint line {} of "{}". Not all atoms form a pair.'.format(num, fragment))   
      pairs = pairwise(line)
      distances = []
      pairlist = []
      for i in pairs:  
        pairlist.append(i)
        try:
          a = atoms[atnames.index(i[0])][1:4]
          b = atoms[atnames.index(i[1])][1:4]
        except(ValueError):
          return
        dist = atomic_distance(a, b, cell)
        distances.append(dist)
      # factor time standard deviation of the SADI distances
      outliers = nalimov_test(distances)
      if outliers:
        for x in outliers:
          print("\n{}:".format(fragment))
          pair = ' '.join(pairlist[x])
          print('Suspicious deviation in atom pair "{}" ({:4.4f} A) of SADI line {}.'.format(pair, distances[x], num+1))
          print(' '.join(restraints[num])[:40], '...')

def nalimov_test(data):
  '''
  returns a index list of outliers base on the Nalimov test for data.
  Modified implementation of:
  "R. Kaiser, G. Gottschalk, Elementare Tests zur Beurteilung von Messdaten
  Bibliographisches Institut, Mannheim 1972." It is tuned to detect only 
  the really bad outliers.
  
  >>> data = [1.120, 1.234, 1.224, 1.469, 1.145, 1.222, 1.123, 1.223, 1.2654, 1.221, 1.215]
  >>> nalimov_test(data)
  [3]
  '''
  # q-values for degrees of freedom:
  f = {1:1.409, 2:1.645, 3:1.757, 4:1.814, 5:1.848, 6:1.870, 7:1.885, 8:1.895,
       9:1.903, 10:1.910, 11:1.916, 12:1.920, 13:1.923, 14:1.926, 15:1.928, 
       16:1.931, 17:1.933, 18:1.935, 19:1.936, 20:1.937, 30:1.945}
  fact = sqrt(float(len(data))/(len(data)-0.4))
  fval = len(data)-2
  if fval < 2:
    return []
  outliers = []
  if fval in f:
    # less strict than the original:
    q_crit = f[fval]+(f[fval]*0.4)
  else:
    q_crit = 1.95+(1.95*0.4)
  for num, i in enumerate(data):
    q = abs(((i-median(data))/std_dev(data))*fact)
    if q > q_crit:
      outliers.append(num)
  return outliers

def invert_atomlist_coordinates(atomst):
    '''
    Inverts atom coordinates
    :param atoms: list of atom list
    >>> c1 = [['c1', '1', '1', '-2', '3'], ['c1', 2, 1, -2, 3]]
    >>> invert_atomlist_coordinates(c1)
    [['c1', '1', -1.0, 2.0, -3.0], ['c1', 2, -1.0, 2.0, -3.0]]
    '''
    atoms = []
    for line in atomst:
      line = list(line)
      try:
        inv_coord = [ -float(i) for i in line[2:] ]
      except:
        print('Unable to invert fragment coordinates.')
        return False
      line[2:] = inv_coord
      atoms.append(line)
    return atoms

def remove_partsymbol(atom):
    '''
    strips the part symbol like C1_4b from an atom name
    :param atom: 'C1_4b'
    :type atom: string
    
    >>> remove_partsymbol('C2_4b')
    'C2_4'
    >>> remove_partsymbol('C22_b')
    'C22'
    >>> remove_partsymbol('C_5')
    'C_5'
    >>> remove_partsymbol('SAME/SADI')
    'SAME/SADI'
    '''
    if '_' in atom:
        presuff = atom.split('_')
        prefix, suffix = presuff[0], presuff[-1].strip(string.ascii_letters)
        if not suffix:
            atom = prefix
        else:
            if suffix == '0':
                atom = prefix
            else:
                atom = prefix+'_'+suffix
    return atom


def flatten(nested):
  '''
  flattens a nested list
  >>> flatten([['wer', 234, 'brdt5'], ['dfg'], [[21, 34,5], ['fhg', 4]]])
  ['wer', 234, 'brdt5', 'dfg', 21, 34, 5, 'fhg', 4]
  '''
  result = []
  try:
    # dont iterate over string-like objects:
    try: nested + ''
    except(TypeError): pass
    else: raise TypeError
    for sublist in nested:
      for element in flatten(sublist):
        result.append(element)
  except(TypeError):
    result.append(nested)
  return result

def frac_to_cart(frac_coord, cell):
  '''
  Converts fractional coordinates to cartesian coodinates
  :param frac_coord: [float, float, float]
  :param cell:       [float, float, float, float, float, float]
    >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
  >>> coord1 = (-0.186843,   0.282708,   0.526803)
  >>> print(frac_to_cart(coord1, cell))
  (-2.741505423999065, 5.909586678000002, 10.775200700893734)
  '''
  a, b, c, alpha, beta, gamma = cell
  x, y, z = frac_coord
  alpha = radians(alpha)
  beta  = radians(beta)
  gamma = radians(gamma)
  cosastar = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma))
  sinastar = sqrt(1-cosastar**2)
  Xc = a*x + (b*cos(gamma))*y + (c*cos(beta))*z
  Yc = 0   + (b*sin(gamma))*y + (-c*sin(beta)*cosastar)*z
  Zc = 0   +  0               + (c*sin(beta)*sinastar)*z
  return(Xc, Yc, Zc)
  
############################################################################
# Experimental:
######################################################################

"""
def get_overlapped_chunks(ring, size):
  '''
  returns a list of chunks of size 'size' which overlap with one field.
  If the last chunk is smaller than size, the last 'size' chunks are returned as last chunk.
  >>> get_overlapped_chunks(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], 4)
  [['a', 'b', 'c', 'd'], ['d', 'e', 'f', 'g'], ['e', 'f', 'g', 'h']]
  '''
  chunks = []
  for i in range(0, len(ring)-size+3, 3):
    chunk = ring[i:i+size]
    if len(chunk) < 4:
      chunk = ring[-size:]
    chunks.append(chunk)
  return chunks

def is_flat(chunk):
  '''
  check if four atoms are flat
  '''
  tetrahedron_atoms = []
  for atom in chunk:
    ### how do I get the coordinates of these atoms?
    single_atom_coordinate = coords_dict[atom]
    tetrahedron_atoms.append(single_atom_coordinate)
  a, b, c, d = tetrahedron_atoms
  volume = (vol_tetrahedron(a, b, c, d))
  if volume < 0.085:
    return True
  else:
    #print('volume of', chunk, 'too big:', volume)
    return False

def get_formated_flats():
  '''
  formats the FLAT restraints and removes the part symbol
  '''
  flats = make_flat_restraints()
  if not flats:
    return ['']
  flat_format = []
  for i in flats:
    #i = [misc.remove_partsymbol(x) for x in i]
    flat_format.append('FLAT {}\n'.format(' '.join(i)))
  return flat_format

def vol_tetrahedron(a, b, c, d, cell=None):
  '''
  returns the volume of a terahedron spanned by four points:
  e.g. A = (3, 2, 1), B = (1, 2, 4), C = (4, 0, 3), D = (1, 1, 7)
          |u1 u2 u3|
  v = 1/6*|v1 v2 v3|
          |w1 w2 w3|
  AB = (1-3, 2-2, 4-1) = (-2, 0, 3)
  AC = ...
  AD = ...
  V = 1/6[u,v,w]
            |-2,  0, 3|
  [u,v,w] = | 1, -2, 2| = 24-3-12 = 5
            |-2, -1, 6|
  V = 1/6*5
  >>> cell = (10.5086, 20.9035, 20.5072, 90, 94.13, 90)
  >>> a = (0.838817,   0.484526,   0.190081) # a ist um 0.01 ausgelenkt
  >>> b = (0.875251,   0.478410,   0.256955)
  >>> c = (0.789290,   0.456520,   0.301616)
  >>> d = (0.674054,   0.430194,   0.280727)
  >>> print('volume of Benzene ring atoms:')
  volume of Benzene ring atoms:
  >>> print(round(vol_tetrahedron(a, b, c, d, cell), 5))
  0.06335
  '''
  A = [float(i) for i in a]
  B = [float(i) for i in b]
  C = [float(i) for i in c]
  D = [float(i) for i in d]
  if cell:
    A = frac_to_cart(a, cell)
    B = frac_to_cart(b, cell)
    C = frac_to_cart(c, cell)
    D = frac_to_cart(d, cell)
  AB = subtract_vect(A, B)
  AC = subtract_vect(A, C)
  AD = subtract_vect(A, D)
  D = determinante([AB, AC, AD])
  volume = abs((D/6))
  return volume

def subtract_vect(a, b):
  '''
  subtract vector b from vector a
  Deprecated, use mpmath instead!!!
  :param a: [float, float, float]
  :param b: [float, float, float]
  
  >>> subtract_vect([1, 2, 3], [3, 2, 2])
  (-2, 0, 1)
  '''
  return (a[0] - b[0],
          a[1] - b[1],
          a[2] - b[2])

def determinante(a):
  '''
  return determinant of 3x3 matrix
  Deprecated, use mpmath instead!!!
  
  >>> m1 = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
  >>> determinante(m1)
  8
  '''
  return (a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
         -a[1][0] * (a[0][1] * a[2][2] - a[2][1] * a[0][2])
         +a[2][0] * (a[0][1] * a[1][2] - a[1][1] * a[0][2]))
"""

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
  import doctest
  failed, attempted = doctest.testmod()#verbose=True)
  if failed == 0:
    print('passed all {} tests!'.format(attempted))