from olexFunctions import OlexFunctions
from collections import OrderedDict
from ImageTools import ImageTools
import StringIO
from PIL import Image, ImageFile, ImageDraw
from helper_functions import check_restraints_consistency, initialize_user_db
import time

OV = OlexFunctions()
IT = ImageTools()

'''
Fragen und Ideen:

- If I fit a fragment (e.g. tert-butyl-n) to an already existing nitrogen, the new nitrogen
  gets deleted and the restraints like SIMU N1 C1 C2 C3 C4 break.
  I need to replace target positions, or no atom!

- If placing a fragment (e.g. toluene) into a negative part the restraints should kept integral for
  this fragment and not expand to symmetry equivalent atoms like
  EQIV $4 -1+X,1+Y,+Z
  DFIX 1.509 0.011 C1 C1_$1 C1_$2 C1_$3 C1_$4 C2 C2_$1 C2_$2 C2_$3 C2_$4

- I would like to replace atoms in 1.3 A around the fitting fragment.

- can the state of the plugin be updated after fit to initialize e.g. the
  residue number again?
  -Y Yes with "mode -e fit"

- check_same_thread=False ?

- Can I get the restraints that are involved with selected atoms?

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
    self.userdb = 'user-fragment-database.sqlite'
    self.userdbfile = os.sep.join([instance_path, os.sep,'db', os.sep, self.userdb])
    if not os.path.exists(self.userdbfile) or os.path.getsize(self.userdbfile) < 100:
      initialize_user_db(self.userdbfile)
    # for edited fragments:
    self.frag_cell = []
    self.db = FragmentTable(self.dbfile, self.userdbfile)

  def init_plugin(self):
    '''
    initialize the plugins main form
    '''
    try:
      int(olx.GetVar('fragment_ID'))
    except(RuntimeError, ValueError):
      return
    self.get_resi_class()
    self.set_fragment_picture(100)
    self.display_image('FDBMOLEPIC', 'displayimg.png')
    self.show_reference()
    resinum = self.find_free_residue_num()
    olx.html.SetValue('RESIDUE', resinum)
    OV.SetParam('fragment_DB.fragment.resinum', resinum)

  def blink_field(self, name, duration=0.1, colour='#aa0000', repeat=5):
    '''
    lets an input field of Olex2 blink
    :param name: name of the html item to blink
    :param duration: time for each blink event
    :param colour: blink colour
    '''
    startcolour = OV.GetParam(name+'.bgcolor')
    colorfield = name#+'.bgcolor'
    # this does not work:
    for n, i in enumerate(range(repeat)):  # @UnusedVariable
      olx.html.SetBG(colorfield,'#aa0000')
      time.sleep(duration)
      olx.html.SetBG(colorfield,'#0000aa')
      time.sleep(duration)
    olx.html.SetBG(colorfield, startcolour)
    
  def set_occu(self, occ):
    '''
    sets the occupancy, even if you enter a comma value instead of point as
    decimal separator.
    '''
    if ',' in occ:
      occ = occ.replace(',', '.')
    try:
      float(occ)
    except(SyntaxError, NameError, ValueError):
      print('Invalid value for occupancy provided')
      return
    varname = 'fragment_DB.fragment.frag_occ'
    OV.SetParam(varname, occ)
    # I want to do this, but it immediately forces the cursor to the left
    # position:
    #olx.html.SetValue('FRAG_OCCUPANCY', OV.GetParam(varname))

  def set_resiclass(self, resiclass, name):
    '''
    sets the residue class and ensures that it is of len(4)
    and .isalpha is the first character.
    '''
    varname = ''
    if name.upper() == 'RESIDUE_CLASS'.upper():
      varname = 'fragment_DB.fragment.resi_class'
    if name.upper() == 'Inputfrag.residue'.upper():
      varname = 'fragment_DB.new_fragment.frag_resiclass'
    if not resiclass:
      return
    if not resiclass[0].isalpha():
      # resiclass does not start with a char:
      OV.SetParam('fragment_DB.fragment.resi_class', resiclass[1:].upper())
      olx.html.SetValue(name, OV.GetParam(varname))
    # force 4 characters:
    elif len(resiclass) > 4:
      OV.SetParam(varname, resiclass[:4].upper())
      olx.html.SetValue(name, OV.GetParam(varname))
    else:
      OV.SetParam(varname, resiclass.upper())
      #olx.html.SetValue(name, OV.GetParam(varname))

  def list_fragments(self):
    '''
    returns the available fragments in the database
    the list of names is separated by semicolon
    i[0] => number
    i[1] => name
    '''
    items = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in self.db])
    return items

  def search_fragments(self, search_string):
    '''
    performs a search for an unsharp name in a list
    '''
    selected_list = ''
    if not search_string:
      selected_results = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in self.db])
    else:
      selected_results = self.db.find_fragment_by_name(search_string, 7)
      selected_list = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in selected_results])
    # propagate the smaller list to the combo-box:
    olx.html.SetItems('LIST_FRAGMENTS', selected_list)
    # show the first result in combo box and intialize the fragment:
    olx.html.SetValue('LIST_FRAGMENTS', '{}'.format(selected_results[0][1]))
    frag_id = selected_results[0][0]
    olx.SetVar('fragment_ID', frag_id)
    self.init_plugin()

  def format_atoms_for_importfrag(self, atoms):
    '''
    format the input atoms to use them with importfarg
    '''
    newlist = []
    finallist = []
    for i in atoms:
      # atoms without sfac:
      newlist.append('{:4.4s} {:>7.4f}  {:>7.4f}  {:>7.4f}'.format(i[0], i[2], i[3], i[4]))
    atoms = '\n'.join(newlist)
    finallist.append('FRAG')
    finallist.append('\n'+atoms)
    finallist.append('\nFEND')
    text = ' '.join(finallist)
    return text

  def onInport(self, atoms):
    '''
    Function is called when ImportFrag exits  
    :param atoms:
    :type atoms:
    '''
    print('Imported atom ids: {}'.format(atoms))
    resinum = OV.GetParam('fragment_DB.fragment.resinum')
    resiclass = OV.GetParam('fragment_DB.fragment.resi_class')
    freevar = OV.GetParam('fragment_DB.fragment.frag_fvar')
    atoms = atoms.split()
    # define the other properties:
    self.define_atom_properties(atoms, resinum, resiclass, freevar)
    OV.unregisterCallback('onFragmentImport', self.onInport)

  def insert_frag_with_ImportFrag(self, fragId, part=1, occ=1,  dfix=False):
    '''
    Input a fragment with ImportFrag
    fvar and resi things are currently not possible in ImportFrag!
    :param fragId: FragmentId
    :param part: SHELX part
    :param occ: occupancy
    :param dfix: generate dfix restraints or not after fit.
    '''
    # this callback runs in the moment when ImportFrag is finished. onInport
    # then defines the further properties of the fragment:
    OV.registerCallback('onFragmentImport', self.onInport)
    fragpath = os.sep.join(['.olex', 'fragment.txt'])
    atoms = self.format_atoms_for_importfrag([ i for i in self.db[fragId]])
    with open(fragpath, 'w') as f:
      f.write(atoms)
    if dfix:
      OV.cmd(r'ImportFrag -p={0} -o={1} -d {2}'.format(part, occ, fragpath))
    else:
      OV.cmd(r'ImportFrag -p={0} -o={1} {2}'.format(part, occ, fragpath))
    return
    
  def fit_db_fragment(self, fragId=None):
    '''
    fit a molecular fragment from the database into olex2
    '''
    if not fragId:
      try:
        fragId = olx.GetVar('fragment_ID')
      except(RuntimeError):
        # no fragment chosen-> do nothing
        return
    partnum = OV.GetParam('fragment_DB.fragment.frag_part')
    occupancy = OV.GetParam('fragment_DB.fragment.frag_occ')
    if OV.GetParam('fragment_DB.fragment.use_dfix'):
      # adding atomids with ImportFrag and DFIX to structure:
      atomids = self.insert_frag_with_ImportFrag(fragId, part=partnum, occ=occupancy, dfix=True)
      return atomids
    else:
      atomids = self.insert_frag_with_ImportFrag(fragId, part=partnum, occ=occupancy, dfix=False)
      return atomids
  
      
  def define_atom_properties(self, atomids, resinum, resiclass, freevar, fragId=None):
    '''
    Defines the atoms properties of the fitted fragment after ImportFrag
    :param atomids:
    :type atomids:
    :param resinum:
    :type resinum:
    :param resiclass:
    :type resiclass:
    :param freevar:
    :type freevar:
    '''
    print('Applying fragment properties:')
    if not fragId:
      try:
        fragId = olx.GetVar('fragment_ID')
      except(RuntimeError):
        # no fragment chosen-> do nothing
        return
    labeldict = OrderedDict()
    atomids = [str(i) for i in atomids]
    dbatom_names = [ i[0] for i in self.db[fragId]]
    for at_id, name in zip(atomids, dbatom_names):
      labeldict[name.upper()] = at_id
    # now residues and otgher stuff:
    if resiclass and resinum:
      self.make_residue(atomids, resiclass, resinum)
    # Placing restraints:
    self.make_restraints(labeldict, fragId)
    # select all atomids to do the fit:
    if freevar != 1:
      OV.cmd("sel #c{}".format(' #c'.join(atomids)))
      OV.cmd("fvar {}".format(freevar))
    # select again, because fvar deselects the fragment
    OV.cmd("sel #c{}".format(' #c'.join(atomids)))
    resinum = self.find_free_residue_num()
    olx.html.SetValue('RESIDUE', resinum)
    OV.SetParam('fragment_DB.fragment.resinum', resinum)
    OV.cmd("sel #c{}".format(' #c'.join(atomids)))
    return atomids


  def is_near_atoms(self):
    '''
    for atom in olex_core.GetRefinementModel(True)['atoms']:
      coord = atom['crd'][0]
    '''
    pass

  def make_residue(self, atoms, resiclass, resinum):
    '''
    selects the atoms and applies "RESI class number" to them
    '''
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("RESI {} {}".format(resiclass, resinum))

  def make_restraints(self, labeldict, fragId):
    '''
    Applies restraints to atoms.
    :param labeldict: a dictionary where the keys are atom names and the
                      values are Olex2 atom Ids
    :type labeldict: dictionary
    :param fragId: the database Id of the fragment
    :type fragId: integer
    '''
    #dfix = OV.GetParam('fragment_DB.fragment.use_dfix')
    #if dfix:
    #  print('calculating DFIX restraints is not implemented in FragmentDB.')
    #  return
    restraints = self.db.get_restraints(fragId)
    if not restraints:
      return
    for num, i in enumerate(restraints):
      # i[0] is restraint like SADI or DFIX
      # i[1] is a string of atoms like 'C1 C2'
      restraint_atoms = i[1]
      if '>' in restraint_atoms or '<' in restraint_atoms:
        restraint_atoms = self.range_resolver(restraint_atoms.split(), labeldict.keys())
      line = []
      for at in restraint_atoms.split():
        # is it a potential atom (starts with alphabetic character):
        if at[0].isalpha():
          try:
            line.append('#c'+labeldict[at])
          except(KeyError):
            # in this case, an atom name in the restraint does not
            # exist in the fragments atom list!
            print('\nUnknown restraint found in line {}.\n'.format(num+1))
            # I must exit here!
            return
        else:
          line.append(at)
      # applies the restraint to atoms in line
      OV.cmd("{} {}".format(i[0], ' '.join(line)))
      #olx.xf.rm.NewRestraint(i[0], ' '.join(line))


  def prepare_picture(self, im, max_size=100):
    '''
    resizes and colorizes the picture to diplay it in olex2
    needs a PIL Image instance
    :param im: PIL Image instance
    :type im: PIL Image instance
    :param max_size: maximum size of the picture in pixels
    :type max_size: int
    :param control: html control to get the background color right
    :type control: string
    '''
    for i in im.size:
      if i < 50:
        print('This picture is too small.')
    im = im.convert(mode="RGBA")
    img_w, img_h = im.size
    ratio = float(max_size) / float(max(im.size))
    # just an empirical value:
    if float(max_size) / float(max(im.size)) > 0.6:
      ratio = 0.6
      # resize equally to fit in max_size
    im = im.resize((int(img_w * ratio), int(img_h * ratio)), Image.ANTIALIAS)
    # empty image of max_size
    bgcolor = self.params.html.table_bg_colour.rgb
    #bgcolor = self.params.html.bg_colour.rgb
    IM = Image.new('RGBA', (max_size, max_size), bgcolor)
    bg_w, bg_h = IM.size
    img_w, img_h = im.size
    # offset for image placement
    offset = (bg_w - img_w) / 2, (bg_h - img_h) / 2
    # place image in center of background image:
    IM.paste(im, offset)
    return IM

  def set_fragment_picture(self, max_size=100):
    '''
    displays a picture of the fragment from the database in Olex2
    :param name: name of the zimg html name
    :type name: string
    :param max_size: maximum size of the picture in pixels
    :type max_size: int
    :param control: name of the htmnl control
    :type control: string
    '''
    max_size = int(max_size)
    fragId = olx.GetVar('fragment_ID')
    pic = self.db.get_picture(fragId)
    if not pic:
      print('No fragment picture found.')
      return False
    im = Image.open(StringIO.StringIO(pic))
    # save it as raw and small pic:
    OlexVFS.save_image_to_olex(im, 'storepic.png', 0)
    im = self.prepare_picture(im, max_size=100)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)

  def display_image(self, zimgname, image_file):
    '''
    display an image in zimgname
    '''
    try:
      olx.html.SetImage(zimgname, image_file)
    except RuntimeError:
      pass

  def store_picture(self, title, ffilter, location, default_name=''):
    '''
    opens a file dialog and stores the selected picture in the olex-vfs
    :param title: Title of the file dialog
    :type title: string
    :param filter: file endings filter like '*.png; *.tiff; *.tif; *.jpg; *.jpeg'
    :type filter: sting
    :param location: location of the path of file dialog
    :type location: string
    :param default_name: default preselected file name for the dialog
    :type default_name: string
    '''
    picfile = olx.FileOpen(title, ffilter, location, default_name)
    if not picfile:
      return
    fsize = os.stat(picfile).st_size
    # do not allow picture of more than 1MB (1000000 bytes):
    if fsize > 1000000:
      print('This file is too large!')
      return
    if fsize <=1:
      print('No valid picture found!')
      return
    im = Image.open(picfile)
    #im = self.prepare_picture(im, max_size=100)
    OlexVFS.save_image_to_olex(im, 'storepic.png', 0)
    # display it.
    im = self.prepare_picture(im, max_size=100)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)
    olx.html.SetImage('Inputfrag.MOLEPIC2', 'displayimg.png')

  def range_resolver(self, restraintat, atom_names):
    '''
    resolves the atom names of ranges like "C1 > C5"
    works for each restraint line separately.
    :param restraintat: atoms with a range definition
    :type restraintat: list
    :param atom_names: names of atoms in the fragment
    :type atom_names: list
    '''
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
    import olex_core  # @UnresolvedImport
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
    try:
      fragId = olx.GetVar('fragment_ID')
    except(RuntimeError):
      return
    try:
      int(fragId)
    except ValueError:
      return
    resiclass = self.db.get_residue_class(int(fragId))
    OV.SetParam('fragment_DB.fragment.resi_class', resiclass)
    OV.SetParam('fragment_DB.new_fragment.frag_resiclass', resiclass)
    # set the class in the text field of the gui:
    olx.html.SetValue('RESIDUE_CLASS', resiclass.upper())

  def show_reference(self, edit=False):
    '''
    show the reference of a fragment in the GUI
    '''
    try:
      fragId = olx.GetVar('fragment_ID')
    except(RuntimeError):
      return
    ref = self.db.get_reference(fragId)
    if not edit:
      olx.html.SetValue('REFERENCE', ref)
    else:
      # reference of the edit window
      olx.html.SetValue('REFERENCE_edit', ref)

  def get_selected_atoms(self):
    '''
    returns the currently selected atoms for the atoms field.
    Hydrogen atoms are discarded.
    '''
    atlist = []
    atoms = []
    atoms_all = olex.f("sel(a)").split()
    if not atoms_all:
      print('No atoms selected!')
      return
    # make a list without H atoms:
    for at in atoms_all:
      if at[0] == 'H' and at[:2] not in ('He', 'Hf', 'Ho', 'Hg'):
        continue
      atoms.append(at)
    atoms.sort()
    crd = [ olx.Crd(x) for x in atoms ]
    # now I want to remove the residue number:
    #atoms = [ y.split('_')[0] for y in atoms]
    for atom, coord in zip(atoms, crd):
      atlist.append('{:4.4s} {:>8.6s} {:>8.6s} {:>8.6s}'.format(atom, *coord.split()))
    at = ' \n'.join(atlist)
    olx.html.SetValue('Inputfrag.SET_ATOM', at)
    self.cell = '1 1 1 90 90 90'
    olx.html.SetValue('Inputfrag.set_cell', self.cell)
    OV.SetParam('fragment_DB.new_fragment.frag_cell', self.cell)
    OV.SetParam('fragment_DB.new_fragment.frag_atoms', at)
    return atlist


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
    return(round(Xc, 6), round(Yc, 6), round(Zc, 6))

  def open_edit_fragment_window(self):
    '''
    opens a new window to input/update a database fragment
    '''
    blank = False
    try:
      olx.GetVar('fragment_ID')
    except RuntimeError:
      blank = True
    pop_name = "Inputfrag"
    screen_height = int(olx.GetWindowSize('gl').split(',')[3])
    screen_width = int(olx.GetWindowSize('gl').split(',')[2])
    box_x = int(screen_width*0.1)
    box_y = int(screen_height*0.1)
    width, height = 540, 660
    path = "{}/inputfrag.htm".format(self.p_path)
    olx.Popup(pop_name, path,  b="tcrp", t="Create/Edit Fragments", w=width,
              h=height, x=box_x, y=box_y)
    if blank:
      self.blank_state()
      return
    else:
      frag = self.get_frag_for_gui()
      if frag:
        self.display_image('Inputfrag.MOLEPIC2', 'displayimg.png')

  def set_frag_name(self, enable_check=True):
    '''
    handles the name of a new/edited fragment
    '''
    fragname = OV.GetParam('fragment_DB.new_fragment.frag_name')
    if self.db.has_exact_name(fragname) and enable_check:
      print('\n{} is already in the database. \nPlease choose a different name.\n'.format(fragname))
      return False
    else:
      OV.SetParam('fragment_DB.new_fragment.frag_name', fragname)
      return True

  def set_frag_cell(self):
    '''
    set the unit cell of a new fragment to convert its coordinates to cartesian
    the resulting cell is retrieved from phil and stored in self.frag_cell
    '''
    self.frag_cell = ''
    try:
      frag_cell = OV.GetParam('fragment_DB.new_fragment.frag_cell').split()
    except(AttributeError):
      print('No unit cell defined. You nedd to define unit cell parameters!')
      return False
    if frag_cell:
      if len(frag_cell) < 6:
        print('\nUnknown unit cell. Only {} values instead of 6.'.format(len(frag_cell)))
        return False
      if len(frag_cell) > 6:
        print('\nUnknown unit cell defined. {} values instead of 6.'.format(len(frag_cell)))
        return False
      try:
        cell = [float(i) for i in frag_cell]
      except(ValueError, TypeError):
        print('Invalid unit cell given!')
        return False
    else:
      print('\nNo unit cell found!')
      print('You have to supply a unit cell!')
      return False
    self.frag_cell = cell
    return True


  def set_frag_atoms(self):
    '''
    -handles the fragment atoms of a new/edited fragment
    -returns a list of lists:
     [['C4', '0.282212', '0.368636', '0.575493'], ...]
    '''
    atoms = OV.GetParam('fragment_DB.new_fragment.frag_atoms')
    try:
      atoms = atoms.split('\n')
      # remove empty lines:
      atoms = [i for i in atoms if i ]
      atoms = [i.split() for i in atoms]
    except AttributeError:
      atoms = None
      return
    finalatoms = self.atoms_parser(atoms)
    return finalatoms

  def atoms_parser(self, atoms):
    '''
    formats the atoms from shelx as long list to a list of list with
    exactly four items: Atom  x y z
    :param atoms: line of atoms
    :type atoms: list
    '''
    # go through all atoms and cut their line to 4 list elements At  x y z:
    for num, line in enumerate(atoms):
      # line with sfac
      if len(line) > 4 and len(line[1]) < 3:
        atoms[num] = [line[0]]+line[2:5]
      # line without sfac, but long
      if len(line) > 4 and len(line[1]) > 3:
        atoms[num] = line[:4]
      if len(line) == 4:
        atoms[num] = line[:4]
      if len(line) < 4:
        # too short, parameters missing
        print('Invalid atom line found!! Parameter(s) missing.')
        del atoms[num]
        continue
      for x in line[1:4]:
        # check if each is a rea number except for the atom:
        try:
          float(x)
        except:
          # if something is wrong, determine the bad guy:
          for i in x:
            if not i.isdigit() and i != '.':
              print('Invalid charachter {} in line.'.format(i))
              del atoms[num]
              continue
          print('Invalid atom line found!')
          del atoms[num]
          continue
    # now rename all the atoms that have no number:
    atomnames = [i[0] for i in atoms]
    for num, line in enumerate(atoms):
      if len(line[0]) == 1:
        # make sure no name is doubled
        while line[0] in atomnames:
          num = num+1
          line[0] = line[0]+str(num)
    return atoms


  def set_frag_restraints(self):
    '''
    handles the restraints of a new/edited fragment
    '''
    restraints = []  # @UnusedVariable
    restr = OV.GetParam('fragment_DB.new_fragment.frag_restraints')
    try:
      restr = restr.split('\n')
      # remove empty lines:
      restr = [i.upper() for i in restr if i ]
      restr = [i.split() for i in restr]
    except AttributeError, e:  # @UnusedVariable
      #print(e)
      restr = None
      return
    return restr


  def prepare_residue_class(self):
    '''
    set the residue class of a new fragment
    '''
    resi_class = OV.GetParam('fragment_DB.new_fragment.frag_resiclass')
    if resi_class:
      return resi_class

  def prepare_reference(self):
    '''
    prepare the reference for a new/edited fragment
    '''
    reference = OV.GetParam('fragment_DB.new_fragment.frag_reference')
    if reference:
      return reference

  def prepare_atoms_list(self, fragId):
    '''
    prepare the atom list to display in a multiline edit field
    '''
    atlist = []
    atoms_list = []
    try:
      atoms_list = self.db[fragId]
    except IndexError:
      print('Database Fragment {} not found.'.format(fragId))
    if not atoms_list:
      return
    atoms_list = [[i for i in y] for y in atoms_list]
    for i in atoms_list:
      atlist.append('{:4.4s} {:>8.4f} {:>8.4f} {:>8.4f}'.format(i[0], i[2], i[3], i[4]))
    at = ' \n'.join(atlist)
    return at

  def prepare_fragname(self, fragId):
    '''
    prepare the fragment name to display in a multiline edit field
    '''
    try:
      name = str(self.db.get_fragment_name(fragId)[0])
    except:
      name = 'Could not get a name from the database.'
    return name


  def prepare_restraints(self, fragId):
    '''
    prepare the fragment restraints to display in a multiline edit field
    '''
    restr_list = self.db.get_restraints(fragId)
    if not restr_list:
      return
    restr_list = [[str(i) for i in y] for y in restr_list]
    restr = '\n'.join(['  '.join(i) for i in restr_list])
    return restr


  def add_new_frag(self):
    '''
    execute this to add a new fragment
    '''
    # add _new_frag must check the restraints vor validity!
    state = self.set_frag_name()
    if not state:
      return
    if not self.set_frag_cell():
      return False
    atoms = self.set_frag_atoms()
    restraints = self.set_frag_restraints()
    resiclass = self.prepare_residue_class()
    reference = self.prepare_reference()
    self.store_new_fragment(atoms, restraints, resiclass, reference)


  def update_fragment(self):
    '''
    execute this to update a fragment
    updates the database information of a fragment
    '''
    fragname = OV.GetParam('fragment_DB.new_fragment.frag_name')
    state = self.set_frag_name(enable_check=False)
    if not state:
      return
    if not self.set_frag_cell():
      return False
    atlines = self.set_frag_atoms()
    restraints = self.set_frag_restraints()
    resiclass = self.prepare_residue_class()
    reference = self.prepare_reference()
    # frag_cell = OV.GetParam('fragment_DB.new_fragment.frag_cell')
    try:
      pic_data = OlexVFS.read_from_olex('storepic.png')
    except TypeError:
      print('No picture found')
      pic_data = ''
    coords = self.prepare_coords_for_storage(atlines)
    if not check_restraints_consistency(restraints, atlines, fragname):
      #self.blink_field('Inputfrag.restraints')
      #OV.Alert('Invalid restraint', 'One of the restraints is invalid. \nNo changes to the database were performed.', 'IROK')
      print('\nFragment was not added to the database!')
      return
    self.delete_fragment(reset=False)
    frag_id = self.db.store_fragment(fragname, coords, resiclass, restraints,
                                      reference, picture=pic_data)
    print('Updated fragment "{0}".'.format(fragname))
    if frag_id:
      olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
      olx.SetVar('fragment_ID', frag_id)
    else:
      print('Something is wrong with fragment storage.')
    self.get_frag_for_gui()
    self.set_fragment_picture(100)
    olx.html.SetValue('RESIDUE_CLASS', '')
    self.show_reference()
    olx.html.SetImage('FDBMOLEPIC', 'blank.png')


  def prepare_coords_for_storage(self, atlines):
    '''
    transform the string formated atoms into a list of atoms lists
    :param atlines: texte lines with atoms from the input field
    :type atlines: string
    '''
    coords = []
    for line in atlines:
      try:
        frac_coord = [ float(i) for i in line[1:4] ]
      except(ValueError):
        print('Invalid coordinate defined in line "{}"'.format(' '.join(line)))
        return
      if len(frac_coord) < 3:
        print('Coordinate value missing in "{}".'.format(' '.join(line)))
        continue
      coord = self.frac_to_cart(frac_coord, self.frag_cell)
      line[1:4] = coord
      coords.append(line)
    return coords

  def get_frag_for_gui(self):
    '''
    get the fragment to display in the multiline edit field
    '''
    fragId = olx.GetVar('fragment_ID')
    at = self.prepare_atoms_list(fragId)
    if not at:
      return False
    name = self.prepare_fragname(fragId)
    restr = self.prepare_restraints(fragId)
    residue = self.prepare_residue_class()
    reference = self.db.get_reference(fragId)
    cell = '1  1  1  90  90  90'
    olx.html.SetValue('Inputfrag.SET_ATOM', at)
    olx.html.SetValue('Inputfrag.set_cell', cell)
    olx.html.SetValue('Inputfrag.set_name', name)
    olx.html.SetValue('Inputfrag.restraints', restr)
    olx.html.SetValue('Inputfrag.residue', residue)
    olx.html.SetValue('Inputfrag.reference_ed', reference)
    OV.SetParam('fragment_DB.new_fragment.frag_name', name)
    OV.SetParam('fragment_DB.new_fragment.frag_atoms', at)
    OV.SetParam('fragment_DB.new_fragment.frag_cell', cell)
    OV.SetParam('fragment_DB.new_fragment.frag_restraints', restr)
    OV.SetParam('fragment_DB.new_fragment.frag_resiclass', residue)
    OV.SetParam('fragment_DB.new_fragment.frag_reference', reference)
    return True


  def store_new_fragment(self, atlines, restraints, resiclass, reference):
    '''
    add a new fragment to the database
    '''
    # check if the given name already exist in the database
    # store fragment with a new number
    fragname = OV.GetParam('fragment_DB.new_fragment.frag_name')
    self.set_frag_cell()
    #self.frag_cell = OV.GetParam('fragment_DB.new_fragment.frag_cell')
    try:
      pic_data = OlexVFS.read_from_olex('storepic.png')
    except TypeError:
      pic_data = ''
    coords = self.prepare_coords_for_storage(atlines)
    print('Adding fragment "{0}" to the database.'.format(fragname))
    if not check_restraints_consistency(restraints, atlines, fragname):
      print('Fragment was not added to the database!')
      #self.blink_field('Inputfrag.restraints')
      #OV.Alert('Invalid restraint', 'One of the restraints is invalid. \nNo changes to the database were performed.', 'OK')
      return
    at_id = self.db.store_fragment(fragname, coords, resiclass, restraints,
                                reference, picture=pic_data)
    if at_id:
      olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    # now get the fragment back from the db to display the new cell:
    olx.SetVar('fragment_ID', at_id)
    self.init_plugin()
    #olx.html.SetValue('RESIDUE_CLASS', '')
    #self.get_resi_class()
    #olx.html.SetImage('FDBMOLEPIC', 'blank.png')
    #self.get_frag_for_gui()


  def blank_state(self):
    olx.html.SetValue('Inputfrag.SET_ATOM', '')
    olx.html.SetValue('Inputfrag.set_cell', '')
    olx.html.SetValue('Inputfrag.set_name', '')
    olx.html.SetValue('Inputfrag.restraints', '')
    olx.html.SetValue('Inputfrag.residue', '')
    olx.html.SetValue('Inputfrag.reference_ed', '')
    olx.html.SetValue('RESIDUE_CLASS', '')
    olx.html.SetImage('Inputfrag.MOLEPIC2', 'blank.png')
    olx.html.SetImage('FDBMOLEPIC', 'blank.png')
    if olx.fs.Exists('displayimg.png'):
      OV.CopyVFSFile('blank.png', 'displayimg.png')

  def delete_fragment(self, reset=True):
    '''
    deletes a fragment from the database
    # Todo: reset all fields (oic, atoms, ...) after deltion
    '''
    fragId = olx.GetVar('fragment_ID')
    del self.db[fragId]
    olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    # Now delete the fields:
    if reset:
      self.blank_state()

  def get_chemdrawstyle(self):
    '''
    opens a file dialog window to save the standard chemdraw style
    '''
    from shutil import copyfile
    title = "Save Chemdrawstyle file"
    ffilter = ''
    location = '.'
    stylename = None  # @UnusedVariable
    default_name = 'chemdraw_style.cds'
    stylename = olx.FileSave(title, ffilter, location, default_name)
    if not stylename:
      return
    spath = "%s/drawstyle.cds" % (self.p_path)
    copyfile(spath, stylename)

  def make_selctions_picture(self):
    '''
    creates a picture from the currently selected fragment in Olex2 and
    stores it in 'storepic.png' as well as 'displaypic.png'
    '''
    # "select with mouse"
    if not olex.f("sel()").split():
      return
    picfile = "fdb_tmp.png"
    OV.cmd("sel atom bonds")
    OV.cmd("showh a False")
    OV.cmd("sel -i")
    OV.cmd("hide")
    OV.cmd("center")
    OV.cmd("sel -i")
    OV.cmd("mpln -n")
    OV.cmd("sel -i")
    OV.cmd("sel atom bonds -i")
    OV.cmd('label -a')
    OV.cmd('pers')
    OV.cmd("pict fdb_tmp.png -nbg")
    OV.cmd('kill labels')
    OV.cmd("showP -m")
    OV.cmd("showh a True")
    OV.cmd('telp 50')
    im = Image.open(picfile)
    im = IT.trim_image(im)
    OlexVFS.save_image_to_olex(im, 'storepic.png', 0)
    # display it.
    im = self.prepare_picture(im, max_size=100)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)
    olx.html.SetImage('Inputfrag.MOLEPIC2', 'displayimg.png')
    try:
      os.remove(picfile)
    except:
      pass

  def write_text_on_image(self, text):
    '''
    write a text on an image to display it in an olex2 zimg
    :param text: text to draw
    :type text: string
    :param zimg_name:
    :type zimg_name:
    '''
    imsize = (600, 600)
    bg_color = self.params.html.table_bg_colour.rgb
    grey_IM = Image.new("RGB", imsize, bg_color) # Create new white image
    draw = ImageDraw.Draw(grey_IM)
    IT = ImageTools()
    IT.write_text_to_draw(draw, text, font_size=60)
    OlexVFS.save_image_to_olex(grey_IM, 'fdb_intro.png', 0)
    #self.display_image(zimg_name, 'fdb_intro.png')



fdb = FragmentDB()

OV.registerFunction(fdb.init_plugin, False, "FragmentDB")
OV.registerFunction(fdb.search_fragments, False, "FragmentDB")
OV.registerFunction(fdb.show_reference,False,"FragmentDB")
OV.registerFunction(fdb.make_selctions_picture,False,"FragmentDB")
OV.registerFunction(fdb.get_selected_atoms,False,"FragmentDB")
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
OV.registerFunction(fdb.get_chemdrawstyle,False,"FragmentDB")
OV.registerFunction(fdb.add_new_frag,False,"FragmentDB")
OV.registerFunction(fdb.update_fragment,False,"FragmentDB")
OV.registerFunction(fdb.delete_fragment,False,"FragmentDB")
OV.registerFunction(fdb.store_picture,False,"FragmentDB")
OV.registerFunction(fdb.display_image,False,"FragmentDB")
