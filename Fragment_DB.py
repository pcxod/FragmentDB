from olexFunctions import OlexFunctions
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
    self.db = FragmentTable(self.dbfile)
    OV.registerFunction(self.list_fragments,True,"FragmentDB")
    OV.registerFunction(self.run,True,"FragmentDB")
    OV.registerFunction(self.fit_db_fragment,True,"FragmentDB")

    
  def list_fragments(self):
    db = FragmentTable(self.dbfile)
    items = ""
    for fragment in db:
      _ = fragment[1]#.replace(',', ' ')
      ID = fragment[0]
      items += "%s<-%s;" %(_, ID)
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

    
  def fit_db_fragment(self):
    fragId = olx.GetVar('fragment_ID')
    db = FragmentTable(self.dbfile)
    resinum = olx.GetVar('resinum')
    resiclass = olx.GetVar('resi_class')
    partnum = olx.GetVar('frag_part')
    occupancy = olx.GetVar('frag_occ')
    atoms = []
    labeldict = {}
    for i in db[fragId]:
      label = str(i[0])
      x, y, z = olx.xf.au.Fractionalise(i[2],i[3],i[4]).split(',')
      id = olx.xf.au.NewAtom(label, x, y, z, False)
      labeldict[label.upper()] = id 
      olx.xf.au.SetAtomPart(id, partnum)
      olx.xf.au.SetAtomOccu(id, occupancy)
      olx.xf.au.SetAtomU(id, 0.04)
      name = olx.xf.au.GetAtomName(id)
      print('adding {}, name: {}, Id: {}, coords: {} {} {}'.format(i[0], name, id, x, y, z))
      atoms.append(id)
    olx.xf.EndUpdate()
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    if resinum:
      OV.cmd("RESI {} {}".format(resiclass, resinum))
    for num, i in enumerate(db.get_restraints(fragId)):
      if '>' in i[1]: # needs a range resolving method
        continue      # ignore ranges for now
      line = []
      for at in i[1].split():
        if at[0].isalpha():
          try:
            line.append('#c'+labeldict[at])
          except(KeyError):
            print('\nBad restraint found in line {}.\n'.format(num))
        else:
          line.append(at)
      OV.cmd("{} {}".format(i[0], ' '.join(line)))
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("mode fit")

  def make_restraints():
    pass
    
FragmentDB_instance = FragmentDB()
