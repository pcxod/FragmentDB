from olexFunctions import OlexFunctions
from collections import OrderedDict
OV = OlexFunctions()

import os
import htmlTools
import olex
import olx
from FragmentDB_handler import FragmentTable


instance_path = OV.DataDir()

p_path = os.path.dirname(os.path.abspath(__file__))
OV.SetVar('FragmentDB_plugin_path', p_path)
p_name = "Fragment_DB"
p_scope = "fragment_DB"
p_htm = "Fragment_DB"
p_img = [("FragmentDB",'h1')]

from PluginTools import PluginTools as PT

class FragmentDB(PT):

  def __init__(self):
    super(FragmentDB, self).__init__()
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    self.setup_gui()
    self.dbfile  = os.sep.join([self.p_path, "fragment-database.sqlite"])
    #self.db = FragmentTable(self.dbfile) # why is it so slow to make the db instance here?
    OV.registerFunction(self.list_fragments,True,"FragmentDB")
    #OV.registerFunction(self.range_resolver,True,"FragmentDB")
    OV.registerFunction(self.make_residue,True,"FragmentDB")
    OV.registerFunction(self.make_restraints,True,"FragmentDB")
    OV.registerFunction(self.fit_db_fragment,True,"FragmentDB")

  def list_fragments(self):
    db = FragmentTable(self.dbfile)
    items = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in db])
    return items

  def run(self):
    db = FragmentTable(self.dbfile)
    db[2]   # get database fragment number 2
    for i in db:
      print(i)  # get all fragment names and their id number.
    for i in db[2]:
      print(i)  # get all atoms of fragment number 2
    for i in db.get_restraints(15):
      print(i)   # get restraints of fragment number 15
    len(db)  # number of fragments in the database
    db.find_fragment_by_name('super', selection=3) # find a fragment with default tolerance search


  def fit_db_fragment(self, fragId=None):
    db = FragmentTable(self.dbfile)
    if not fragId:
      try:
        fragId = olx.GetVar('fragment_ID')
      except(RuntimeError):
        # no fragment chosen-> do nothing
        return
    resinum = OV.GetParam('fragment_DB.fragment.resinum')
    resiclass = OV.GetParam('fragment_DB.fragment.resi_class')
    partnum = OV.GetParam('fragment_DB.fragment.frag_part')
    occupancy = OV.GetParam('fragment_DB.fragment.frag_occ')
    freevar = OV.GetParam('fragment_DB.fragment.frag_fvar')
    print('#'*40)
    print('resinum, resiclass, partnum, freevar, occupancy:', resinum, resiclass, partnum, freevar, occupancy)
    print('#'*40)
    atoms = []
    atom_names = []
    labeldict = OrderedDict()
    # adding atoms to structure:
    for i in db[fragId]:
      label = str(i[0])
      x, y, z = olx.xf.au.Fractionalise(i[2],i[3],i[4]).split(',')
      id = olx.xf.au.NewAtom(label, x, y, z, False)
      labeldict[label.upper()] = id
      olx.xf.au.SetAtomPart(id, partnum)
      olx.xf.au.SetAtomU(id, 0.045)
      name = olx.xf.au.GetAtomName(id)
      atom_names.append(name)
      print('adding {}, name: {}, Id: {}, coords: {} {} {}'.format(i[0], name, id, x, y, z))
      atoms.append(id)
    olx.xf.EndUpdate()
    # now residues and otgher stuff:
    if resiclass and resinum:
      self.make_residue(atoms, resiclass, resinum)
    # Placing restraints:
    self.make_restraints(atoms, db, labeldict, fragId, atom_names)
    # select all atoms to do the fit:
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("fvar {} {}".format(freevar, occupancy))
    # select again, because fvar deselects the fragment
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("mode fit")

  def make_residue(self, atoms, resiclass, resinum):
    '''
    selects the atoms and applies "RESI class number" to them
    '''
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("RESI {} {}".format(resiclass, resinum))

  def make_restraints(self, atoms, db, labeldict, fragId, atom_names):
    '''
    applies restraints to atoms
    TODO:
    - resolve ranges like 'C1 > C5' or 'C5 < C2'
    '''
    for num, i in enumerate(db.get_restraints(fragId)):
      # i[0] is restraint like SADI or DFIX
      # i[1] is a string of atoms like 'C1 C2'
      restraintat = i[1]
      if '>' in restraintat or '<' in restraintat:
        restraintat = self.range_resolver(restraintat, atom_names)
      line = []
      for at in restraintat.split():
        if at[0].isalpha():
          try:
            line.append('#c'+labeldict[at])
          except(KeyError):
            # in this case, an atom name in the restraint does not 
            # exist in the fragments atom list!
            print('\nBad restraint found in line {}.\n'.format(num))
        else:
          line.append(at)
      # applies the restraint to atoms in line
      OV.cmd("{} {}".format(i[0], ' '.join(line)))

  def range_resolver(self, restraintat, atom_names):
    '''
    resolves the atom names of ranges like "C1 > C5"
    works for each restraint line separately.
    TODO:
    - does not work for SAME, need to resolve SADI!
    '''
    restraintat = restraintat.split()
    # dict with lists of indexes of the > or < sign:
    rightleft = {'>':[], '<': []}
    for rl in rightleft:
      for num, i in enumerate(restraintat):
        if rl == i:
          # fill the dictionary:
          rightleft[rl].append(num)
    for rl in rightleft:
      # for each sign:
      for i in rightleft[rl]:
        # for each position of < or >:
        if rl == '>':
          # forward range
          left = atom_names.index(restraintat[i-1])+1
          right = atom_names.index(restraintat[i+1])
          restraintat[i:i+1] = atom_names[left:right]
        else:
          # backward range
          left = atom_names.index(restraintat[i-1])
          right = atom_names.index(restraintat[i+1])+1
          names = atom_names[right:left]
          names.reverse()
          restraintat[i:i+1] = names
    return ' '.join(restraintat)
    

FragmentDB_instance = FragmentDB()

'''
    import olex_core
    for r in olex_core.GetRefinementModel(true)['aunit']['residues']:
      print(r)      
'''


