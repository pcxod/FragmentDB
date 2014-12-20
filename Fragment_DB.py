print('halllo DBrunnner')
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
    OV.registerFunction(self.run,True,"FragmentDB")
    OV.registerFunction(self.match_dbfrag,True,"FragmentDB")
    OV.registerFunction(self.list_fragments,True,"FragmentDB")
    OV.registerFunction(self.fit_frag,True,"FragmentDB")
    OV.registerFunction(self.set_part,True,"FragmentDB")
    OV.registerFunction(self.set_resi,True,"FragmentDB")
    OV.registerFunction(self.set_occ,True,"FragmentDB")
    OV.registerFunction(self.set_fvar,True,"FragmentDB")
  
   
  def set_part(self):
    frag_part = olx.GetVar('frag_part')
    print(frag_part)
  
  def set_resi(self):
    frag_resi = olx.GetVar('frag_resi')
    print(frag_resi)
  
  def set_occ(self):
    frag_occ = olx.GetVar('frag_occ')
    print(frag_occ)
  
  def set_fvar(self):
    frag_fvar = olx.GetVar('frag_fvar')
    print(frag_fvar)

    
  def list_fragments(self):
    dbfile  = os.sep.join([self.p_path, "dk-database.sqlite"])
    db = FragmentTable(dbfile)
    items = ""
    for fragment in db:
      _ = fragment[1].replace(',', ' ')
      ID = fragment[0]
      items += "%s<-%s;" %(_, ID)
    #olx.html.SetItems("LIST_FRAGMENTS", items)
    return items

  def fit_frag(self):
    fragId = olx.GetVar('fragment_ID')
    self.match_dbfrag(fragId)
    
  def run(self):
    dbfile  = os.sep.join([self.p_path, "dk-database.sqlite"])
    db = FragmentTable(dbfile)
    db[2]   # get database fragment number 2
    for i in db:
      print(i)  # get all fragment names and their id number.
    for i in db[2]:
      print(i)  # get all atoms of fragment number 2
    for i in db.get_restraints(15):
      print(i)   # get restraints of fragment number 15
    len(db)  # number of fragments in the database

    ## store a fragment into the database.
    #id = db.store_fragment(fragment_name, atoms, restraints, tag)
    #if id:
      #print('fragment is stored successfully')

    db.find_fragment_by_name('super', selection=3) # find a fragment with default tolerance search

  def match_dbfrag(self, fragId=None):
    dbfile  = os.sep.join([self.p_path, "dk-database.sqlite"])
    db = FragmentTable(dbfile)
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
    for num, i in enumerate(db.get_restraints(fragId)):
      if '>' in i[1]: # needs a renge resolving method
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


    
FragmentDB_instance = FragmentDB()
