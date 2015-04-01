from olexFunctions import OlexFunctions
from collections import OrderedDict
OV = OlexFunctions()

'''

Fragen und Ideen:

- add search field which decreases the options in the input combobox.
- possibility to add sump to the free variable
- Checkbox for "use DFIX"
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
- On MAC systems the drop down selector can not be toggeled with the arrow keys.
- das bild in ein <a href="function"> BILD</a> einbinden

- maybe first apply relative restraints and then analyze the residuals. If they are bad try automatic 
  generated direct restraints. 
  
- which format should the atoms in the "insert atoms" windows have?
  Name   element  x  y  z
  C1     6        1  2  3
- what should we do about duplicated atom labels? should we do anything? seem to work fine.
- how can I set a fragment as default upon startup and activate it's picture? 
- Olex2 should not delete any atom while fragment fit/AddAtom()
- I might insert new atoms just by selected atoms? with a button "from selection"?
- I would like to replace atoms in 1 A around the fitting fragment. I should implement a replace mode.
- OV.registerFunction() should warn about initialisation of non-existent functions
- How can I reload the db after deletion of a fragment?
'''


import os
import olex
import gui
import olx
import OlexVFS
from FragmentDB_handler import FragmentTable, SHX_CARDS


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
    # for edited fragments:
    self.frag_cell = []
    #self.db = FragmentTable(self.dbfile)
    #OV.registerFunction(self.print_func,True,"FragmentDB")
    #self.print_func()



  def print_func(self):
    import olex_core
    #l = olex_core.ExportFunctionList()
    #for i in l:
    #  print(i)
    for i in olex_core.ExportFunctionList():
      for y in i:
        try:
          print(y)
        except:
          pass
          
          
  def set_occu(self, occ):
    '''
    sets the occupancy, even if you enter a comma value instead of point as 
    decimal separator.
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
    sets the residue class and ensures that it is of len(4) 
    and .isalpha is the first character.
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
    labeldict = OrderedDict()
    #OV.GetCurrentSelection()
    # adding atoms to structure:
    for i in db[fragId]:
      label = str(i[0])
      x, y, z = olx.xf.au.Fractionalise(i[2],i[3],i[4]).split(',')
      id = olx.xf.au.NewAtom(label, x, y, z, False)
      olx.xf.au.SetAtomPart(id, partnum)
      # if label is H... then SetAtomU == -1.3
      olx.xf.au.SetAtomU(id, 0.045)
      olx.xf.au.SetAtomOccu(id, occupancy)
      name = olx.xf.au.GetAtomName(id)
      labeldict[name.upper()] = id
      print('adding {}, Id: {}, coords: {} {} {}'.format(i[0], id, x, y, z))
      atoms.append(id)
    olx.xf.EndUpdate()
    # now residues and otgher stuff:
    if resiclass and resinum:
      self.make_residue(atoms, resiclass, resinum)
    # Placing restraints:
    self.make_restraints(db, labeldict, fragId)
    # select all atoms to do the fit:
    if freevar != 1:
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

  def make_restraints(self, db, labeldict, fragId):
    '''
    applies restraints to atoms
    '''
    for num, i in enumerate(db.get_restraints(fragId)):
      # i[0] is restraint like SADI or DFIX
      # i[1] is a string of atoms like 'C1 C2'
      restraint_atoms = i[1]
      if '>' in restraint_atoms or '<' in restraint_atoms:
        restraint_atoms = self.range_resolver(restraint_atoms, labeldict.keys())
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
            # I must exit here!
            return
        else:
          line.append(at)
      # applies the restraint to atoms in line
      OV.cmd("{} {}".format(i[0], ' '.join(line)))

  
  def set_fragment_picture(self, name, max_size=150):
    '''
    displays a picture of the fragment from the database in Olex2
    '''
    max_size = int(max_size)
    from PIL import Image, ImageFile
    import StringIO
    import OlexVFS
    import random
    fragId = olx.GetVar('fragment_ID')
    db = FragmentTable(self.dbfile)
    pic = db.get_picture(fragId)
    if not pic:
      print('No fragment picture found.')
      return
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
    #randnum = random.randint(0, 999) 
    OlexVFS.save_image_to_olex(IM, 'pic_{0}.png'.format(fragId), 0)
    # display it.
    olx.html.SetImage(name, 'pic_{0}.png'.format(fragId))
  
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
  
  def get_resi_class(self):
    '''
    sets the residue class from the respective database fragment.
    '''
    db = FragmentTable(self.dbfile)
    try:
      fragId = olx.GetVar('fragment_ID')
    except(RuntimeError):
      return
    try:
      int(fragId)
    except ValueError:
      return
    resiclass = db.get_residue_class(fragId)
    # set the class in the text field of the gui:
    if OV.IsControl('RESIDUE_CLASS'):
      olx.html.SetValue('RESIDUE_CLASS', resiclass)
    OV.SetParam('fragment_DB.fragment.resi_class', resiclass)

###############################################################################

  def frac_to_cart(self, frac_coord, cell):
    '''
    Converts fractional coordinates to cartesian coodinates
    :param frac_coord: [float, float, float]
    :param cell:       [float, float, float, float, float, float]
    '''
    from math import cos, sin, sqrt, radians
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
    return (Xc, Yc, Zc)
  
  def open_edit_fragment_window(self):
    '''
    opens a new window to input/update a database fragment
    '''
    pop_name = "Inputfrag"
    screen_height = int(olx.GetWindowSize('gl').split(',')[3])
    screen_width = int(olx.GetWindowSize('gl').split(',')[2])
    box_x = int(screen_width*0.1)
    box_y = int(screen_height*0.1)
    width, height = 530, 680
    path = "%s/inputfrag.htm" % (self.p_path)
    olx.Popup(pop_name, path,  b="tcrp", t="Create/Edit Fragments", w=width, h=height,
              x=box_x, y=box_y)
    #if res == '1':
    #  affiliation.name = olx.html.GetValue('Affiliation.AFFILIATION_NAME')

  
  def set_frag_name(self, enable_check=True):
    '''
    handles the name of a new/edited fragment
    '''
    fragname = ""
    fragname = OV.GetParam('fragment_DB.new_fragment.frag_name')
    if self.check_name(fragname) and enable_check:
      print('\n{} is already in the database. \nPlease choose a different name.\n'.format(fragname))
      return False
    else:
      OV.SetParam('fragment_DB.new_fragment.frag_name', fragname)
      return True
    
  def check_name(self, name):
    '''
    check if name is already present in the db
    Acetone, C3H6O
    '''
    db = FragmentTable(self.dbfile)
    if db.has_exact_name(name):
      return True
    return False
  
  def set_frag_cell(self):
    '''
    set the unit cell of a new fragment to convert its coordinates to cartesian
    '''
    self.frag_cell = ''
    frag_cell = OV.GetParam('fragment_DB.new_fragment.frag_cell')
    if frag_cell:
      try:
        cell = [float(i) for i in frag_cell.split()]
      except ValueError, TypeError:
        print('Bad unit cell given!')
    self.frag_cell = cell


  def set_frag_atoms(self):
    '''
    handles the fragment atoms of a new/edited fragment
    '''
    atlines = []
    atoms = OV.GetParam('fragment_DB.new_fragment.frag_atoms')
    try:
      atoms = atoms.split()
    except AttributeError:
      atoms = None
      return
    try:
      for i in range(0, len(atoms), 5):
          atomline = atoms[i:i+5]
          if len(atomline) < 5:
            continue
          atlines.append(' '.join(atomline))
    except TypeError:
      pass
    return atlines

  

  def set_frag_restraints(self):
    '''
    handles the restraints of a new/edited fragment
    '''
    restraints = []
    restr = OV.GetParam('fragment_DB.new_fragment.frag_restraints')
    if restr:
      line = restr.split()
    else:
      return
    for n, i in enumerate(line):
      if i[:4] in SHX_CARDS:
        line[n] = '\n'+line[n]
    restraints = ' '.join(line).strip().split('\n')
    return restraints


  def set_frag_resiclass(self):
    '''
    set the residue class of a new fragment
    '''
    resi_class = OV.GetParam('fragment_DB.new_fragment.frag_resiclass')
    if resi_class:
      OV.SetParam('fragment_DB.new_fragment.frag_resiclass', resi_class)
      return resi_class
    
  """
  def check_resiclass_name(self):
          if not self.resiclass[0].isalpha():
        # resiclass does not startt with a char:
        OV.SetParam('fragment_DB.new_fragment.frag_resiclass', '')
        if OV.IsControl('residue'):
          olx.html.SetValue('residue', '')
  """
  
  
  """
  def check_residue(self, residue):
    '''
    check if residue class is already in the db
    '''
    pass
  """
 
  def prepare_atoms_list(self, db, fragId):
    '''
    prepare the atom list to display in a multiline edit field
    '''
    try:
      atoms_list = db[fragId]
    except IndexError:
      print('Database Fragment {} not found.'.format(fragId))
    if not atoms_list:
      return
    atoms_list = [[str(i) for i in y] for y in atoms_list]
    at = ' \n'.join([' '.join(i) for i in atoms_list])
    return at

  def prepare_fragname(self, db, fragId):
    '''
    prepare the fragment name to display in a multiline edit field
    '''
    try:
      name = str(db.get_fragment_name(fragId)[0])
    except:
      name = 'Could not get a name from the database.'
    return name
      
  
  def prepare_restraints(self, db, fragId):
    '''
    prepare the fragment restraints to display in a multiline edit field
    ''' 
    restr_list = db.get_restraints(fragId)
    restr_list = [[str(i) for i in y] for y in restr_list]
    restr = ' \n'.join(['  '.join(i) for i in restr_list])
    return restr
  
  
  def add_new_frag(self):
    '''
    execute this to add a new fragment
    '''
    # add _new_frag must check the restraints vor validity!
    state = self.set_frag_name()
    if not state:
      return
    self.set_frag_cell()
    atoms = self.set_frag_atoms()
    restraints = self.set_frag_restraints()
    self.set_frag_resiclass()
    self.store_new_fragment(atoms, restraints)

  
  def update_fragment(self):
    '''
    execute this to update a fragment

    updates the database information of a fragment
    '''    
    state = self.set_frag_name(enable_check=False)
    if not state:
      return
    self.set_frag_cell()
    atoms = self.set_frag_atoms()
    restraints = self.set_frag_restraints()
    self.set_frag_resiclass()
    self.store_new_fragment(atoms, restraints)
  
  def get_frag_for_gui(self):
    '''
    get the fragment to display in the multiline edit field
    '''
    fragId = olx.GetVar('fragment_ID')
    db = FragmentTable(self.dbfile)
    at = self.prepare_atoms_list(db, fragId)
    if not at:
      return
    name = self.prepare_fragname(db, fragId)
    restr = self.prepare_restraints(db, fragId)
    residue = db.get_residue_class(fragId)
    cell = '1  1  1  90  90  90'
    olx.html.SetValue('Inputfrag.SET_ATOM', at)
    olx.html.SetValue('Inputfrag.set_cell', cell)
    olx.html.SetValue('Inputfrag.set_name', name)
    olx.html.SetValue('Inputfrag.restraints', restr)
    olx.html.SetValue('Inputfrag.residue', residue)
    OV.SetParam('fragment_DB.new_fragment.frag_name', name)
    OV.SetParam('fragment_DB.new_fragment.frag_atoms', at)
    OV.SetParam('fragment_DB.new_fragment.frag_cell', cell)
    OV.SetParam('fragment_DB.new_fragment.frag_restraints', restr)
    OV.SetParam('fragment_DB.new_fragment.frag_resiclass', residue)



  def display_image(self, zimg):
    '''
    display the zimg
    '''
    olx.html.SetImage(self, zimg_name, image_file)
    
      
  def store_new_fragment(self, atlines, restraints):
    '''
    add a new fragment to the database
    '''
    # check if the given name already exist in the database
    # store fragment with a new number
    db = FragmentTable(self.dbfile)
    fragname = OV.GetParam('fragment_DB.new_fragment.frag_name')
    resiclass = OV.GetParam('fragment_DB.new_fragment.frag_resiclass')
    frag_cell = OV.GetParam('fragment_DB.new_fragment.frag_cell')
    coords = []
    for line in atlines:
      line = line.split()
      frac_coord = [ float(i) for i in line[2:5] ]
      coord = self.frac_to_cart(frac_coord, self.frag_cell)
      line[2:5] = coord
      coords.append(line)
    print('Adding fragment "{0}" to the database.'.format(fragname))
    #print('Atoms:', coords)
    #print('Restraints:', restraints)
    #print('Residue:', resiclass)
    id = db.store_fragment(fragname, coords, resiclass, restraints)
    if id:
      olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
      olx.html.SetItems('Inputfrag.LIST_INPFRAGMENTS', self.list_fragments())
    # now get the fragment back from the db to display the new cell:
    olx.SetVar('fragment_ID', id)
    #self.get_frag_for_gui()
    
  def delete_fragment(self):
    '''
    deletes a fragment from the database
    '''
    db = FragmentTable(self.dbfile)
    fragId = olx.GetVar('fragment_ID')
    del db[fragId]
    olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    olx.html.SetItems('Inputfrag.LIST_INPFRAGMENTS', self.list_fragments())
    #reload(FragmentDB)
                   #'Inputfrag.LIST_INPFRAGMENTS'
    
# olx.html.SetEnabled(self.ctrl("UpdatePerson"), True)    

fdb = FragmentDB()
OV.registerFunction(fdb.open_edit_fragment_window,False,"FragmentDB")
OV.registerFunction(fdb.list_fragments,False,"FragmentDB")
OV.registerFunction(fdb.fit_db_fragment,False,"FragmentDB")
OV.registerFunction(fdb.get_resi_class,False,"FragmentDB")
OV.registerFunction(fdb.find_free_residue_num,False,"FragmentDB")
OV.registerFunction(fdb.get_frag_for_gui,False,"FragmentDB")
OV.registerFunction(fdb.set_occu,False,"FragmentDB")
OV.registerFunction(fdb.set_resiclass,False,"FragmentDB")
OV.registerFunction(fdb.store_new_fragment,False,"FragmentDB")
OV.registerFunction(fdb.set_fragment_picture,False,"FragmentDB")

OV.registerFunction(fdb.check_name,False,"FragmentDB")
OV.registerFunction(fdb.set_frag_name,False,"FragmentDB")
OV.registerFunction(fdb.set_frag_cell,False,"FragmentDB")

OV.registerFunction(fdb.add_new_frag,False,"FragmentDB")
OV.registerFunction(fdb.update_fragment,False,"FragmentDB")

OV.registerFunction(fdb.set_frag_atoms,False,"FragmentDB")
OV.registerFunction(fdb.set_frag_restraints,False,"FragmentDB")
OV.registerFunction(fdb.set_frag_resiclass,False,"FragmentDB")

OV.registerFunction(fdb.delete_fragment,False,"FragmentDB")
OV.registerFunction(fdb.update_fragment,False,"FragmentDB")


