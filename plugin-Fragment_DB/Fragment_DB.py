from olexFunctions import OlexFunctions
from collections import OrderedDict
OV = OlexFunctions()

'''

Fragen und Ideen:

- Make GUI more compact.
- how can I color already used residue numbers red in the spinner?
- Problem solved by using olx.xf.au.SetAtomOccu(id, occupancy), but now Olex2 tells me an 
  occupancy of 2 directly after the fit, but the value in the file and after refinement is the correct 
  value of 1.
- I would like to see the residue numbers of the atoms on changing the number 
  like with part numbers
- fit fragment to or near selected atoms/Q-peaks (cctbx model_matches())
- If atom is near other atom (< 1/2*wavelength) and has same name (if in same resi class) or same atom type:
  make them eadp. Maybe an extra button? Or is it possible to start something after "mode fit"?
- observe disagreeable restraints and warn user
- can the state of the plugin be updated after fit to initialize e.g. the residue number again?
- How should I handle hydrogen atoms from water? They should get constraints for vibrations!
- having to put symmetric fragments into negative part is a problem!
- build routines to insert new fragments to db for users
- a working SAME_resiclass atomlist would be great for repeating residues!
- If I fit a fragment (e.g. tert-butyl-n) to an already existing nitrogen, the new nitrogen is not fitted
  (which is ok) but the restraints like SIMU N1 C1 C2 C3 C4 miss the nitrogen.
- in mode fit: clicking two times the same atom pair makes strange things.
- If placing a fragment (e.g. toluene) into a negative part the restraints should kept integral for this fragment 
  and not expand to symmetry equivalent atoms like 
  EQIV $4 -1+X,1+Y,+Z
  DFIX 1.509 0.011 C1 C1_$1 C1_$2 C1_$3 C1_$4 C2 C2_$1 C2_$2 C2_$3 C2_$4
- Is it possible to use javascript in Olex2? For JSmol for example?
- How can I set a default value for snippets/input-combo out of the list of fragments?

- das bild in ein <a href="function"> BILD</a> einbinden

- maybe first apply relative restraints and then analyze the residuals. If they are bad try automatic 
  generated direct restraints. 
  
'''


import os
import htmlTools
import olex
import olx
import OlexVFS
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
    self.params = OV.GuiParams()
    self.dbfile  = os.sep.join([self.p_path, "fragment-database.sqlite"])
    #self.db = FragmentTable(self.dbfile) # why is it so slow to make the db instance here?
    OV.registerFunction(self.list_fragments,True,"FragmentDB")
    OV.registerFunction(self.fit_db_fragment,True,"FragmentDB")
    OV.registerFunction(self.resi_class,True,"FragmentDB")
    OV.registerFunction(self.find_free_residue_num,True,"FragmentDB")
    OV.registerFunction(self.set_occu,True,"FragmentDB")
    OV.registerFunction(self.set_resiclass,True,"FragmentDB")
    OV.registerFunction(self.set_fragment_picture,True,"FragmentDB")
    #OV.registerFunction(self.print_func,True,"FragmentDB")
    #self.print_func()

  def print_func(self):
    import olex_core
    for i in olex_core.ExportFunctionList():
      for y in i:
        try:
          print(y)
        except:
          pass
          
          
  def set_occu(self, occ):
    '''
    sets the occupancy, even if you enter a comma valueinstead of point as 
    dcimal separator.
    '''
    if ',' in occ:
      occ = occ.replace(',', '.')
    else:
      try:
        float(occ)
      except(SyntaxError, NameError, ValueError):
        print('bad value for occupancy provided')
        return
    OV.SetParam('fragment_DB.fragment.frag_occ', occ)

  def set_resiclass(self, resiclass):
    '''
    sets the residue class and ensures that is of len 4 
    and .isalpha is the first char.
    '''
    if not resiclass[0].isalpha():
      # resiclass does not startt with a char:
      OV.SetParam('fragment_DB.fragment.resi_class', '')
    # force 4 characters:
    elif len(resiclass) > 4:
      OV.SetParam('fragment_DB.fragment.resi_class', resiclass[:4])
    else:
      OV.SetParam('fragment_DB.fragment.resi_class', resiclass)
    
  def list_fragments(self):
    '''
    returns the available fragments in the database
    the list of names is separated by semicolon
    '''
    db = FragmentTable(self.dbfile)
    items = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in db])
    return items

  def fit_db_fragment(self, fragId=None):
    '''
    fit a molecular fragment from the database into olex2
    '''
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
      # if label is H... then SetAtomU == -1.3
      olx.xf.au.SetAtomU(id, 0.045)
      olx.xf.au.SetAtomOccu(id, occupancy)
      name = olx.xf.au.GetAtomName(id)
      atom_names.append(name.upper())
      print('adding {}, Id: {}, coords: {} {} {}'.format(i[0], id, x, y, z))
      atoms.append(id)
    olx.xf.EndUpdate()
    # now residues and otgher stuff:
    if resiclass and resinum:
      self.make_residue(atoms, resiclass, resinum)
    # Placing restraints:
    self.make_restraints(atoms, db, labeldict, fragId, atom_names)
    # select all atoms to do the fit:
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("fvar {}".format(freevar))
    # select again, because fvar deselects the fragment
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("mode fit")
    return atoms

  def make_residue(self, atoms, resiclass, resinum):
    '''
    selects the atoms and applies "RESI class number" to them
    '''
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("RESI {} {}".format(resiclass, resinum))

  def make_restraints(self, atoms, db, labeldict, fragId, frag_atom_names):
    '''
    applies restraints to atoms
    '''
    for num, i in enumerate(db.get_restraints(fragId)):
      # i[0] is restraint like SADI or DFIX
      # i[1] is a string of atoms like 'C1 C2'
      restraint_atoms = i[1]
      if '>' in restraint_atoms or '<' in restraint_atoms:
        restraint_atoms = self.range_resolver(restraint_atoms, frag_atom_names)
      line = []
      for at in restraint_atoms.split():
        # is it a potential atom:
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

  
  def set_fragment_picture(self, max_size=150):
    '''
    displays a picture of the fragment from the database in Olex2
    '''
    max_size = int(max_size)
    from PIL import Image, ImageFile
    import StringIO
    import OlexVFS
    fragId = olx.GetVar('fragment_ID')
    db = FragmentTable(self.dbfile)
    pic = db.get_picture(fragId)
    im = Image.open(StringIO.StringIO(pic))
    im = im.convert(mode="RGBA")
    img_w, img_h = im.size
    ratio = float(max_size)/float(max(im.size))
    # just an empirical value:
    if float(max_size)/float(max(im.size)) > 0.6:
      ratio = 0.6
    # resize equally to fit in max_size 
    im = im.resize((int(img_w*ratio), int(img_h*ratio)), Image.ANTIALIAS)
    # empty image of max_size
    IM = Image.new('RGBA', (max_size,max_size), self.params.html.table_bg_colour.rgb)
    bg_w, bg_h = IM.size
    img_w, img_h = im.size
    # offset for image placement
    offset = ((bg_w - img_w) / 2, (bg_h - img_h) / 2)
    # place image in center of background image:
    IM.paste(im, offset)
    # save it
    OlexVFS.save_image_to_olex(IM, 'pic.png', 0)
    # display it
    olx.html.SetImage('MOLEPIC', 'pic.png')
  
  def range_resolver(self, restraintat, atom_names):
    '''
    resolves the atom names of ranges like "C1 > C5"
    works for each restraint line separately.
    TODO:
    - does not work for SAME, need to resolve SADI!
    '''
    restraintat = restraintat.split()
    # dict with lists of positions of the > or < sign:
    rightleft = {'>':[], '<': []}
    for rl in rightleft:
      for num, i in enumerate(restraintat):
        i = i.upper()
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
          names.reverse() # counting backwards
          restraintat[i:i+1] = names
    return ' '.join(restraintat)

  def find_free_residue_num(self):
    '''
    Determines which residue number is unused.
    Either finds the first unused residue in the list of residue numbers
    or returns the last used + 1.
    '''
    residues = self.get_residue_numbers()
    # find unused numbers in the list:
    resi = False
    try:
      for i, j in zip(residues, range(1, residues[-1]+1)):
        # gap in list found:
        if i != j:
          resi = j
          break
    except(IndexError):
      return 1
    # no gap, thus use next number:
    if not resi:
      try:
        resi = residues[-1]+1
      except(IndexError):
        # no residue at all in the structure
        return 1
    return resi
    
  def get_residue_numbers(self):
    '''
    returns a list of residue numbers in the structure
    '''
    import olex_core
    residues = []
    # get the properties of the atoms:
    for r in olex_core.GetRefinementModel(True)['aunit']['residues']:
      try:
        # atoms in residue 0 have no 'number'
        residues.append(r['number'])
      except:
        pass
    residues.sort()
    return residues
  
  def resi_class(self):
    '''
    sets the residue class from the respective database fragment.
    '''
    db = FragmentTable(self.dbfile)
    try:
      fragId = olx.GetVar('fragment_ID')
    except(RuntimeError):
      #print('could not get fragments ID!')
      return
    resiclass = db.get_residue_class(fragId)
    if OV.IsControl('RESIDUE_CLASS'):
      olx.html.SetValue('RESIDUE_CLASS', resiclass)
    OV.SetParam('fragment_DB.fragment.resi_class', resiclass)




FragmentDB_instance = FragmentDB()
