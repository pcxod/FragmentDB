from __future__ import print_function

import StringIO
import os
import pprint
from collections import OrderedDict

import gui.maps
import olex
import olex_core
import olx
from PIL import Image, ImageChops
from olexFunctions import OlexFunctions

import OlexVFS
import helper_functions
from ImageTools import ImageTools
from fragmentdb_handler import FragmentTable
from helper_functions import check_restraints_consistency, initialize_user_db, \
  invert_atomlist_coordinates, frac_to_cart, atomic_distance
from refine_model_tasks import Refmod

OV = OlexFunctions()
IT = ImageTools()
# FragmentDB version number:
FDB_VERSION = 13

r'''
Ideas:
- onreturn="html.Call(~name~.onchange)"
- list of rings for FLAT restraints
- http://interactivepython.org/courselib/static/pythonds/Graphs/Implementation.html
- spy.styles.set_rim_colour(#hex)
- labels -p does not work anymore in dev.
- sel resi 1 does not work anymore in dev.
- what should we do with SADI C1 C2
- residues as default?
- live invert during key press
- better preview needed!
'''

instance_path = OV.DataDir()

p_path = os.path.dirname(os.path.abspath(__file__))
OV.SetVar('FragmentDB_plugin_path', p_path)
p_name = "FragmentDB"
p_scope = "FragmentDB"
p_htm = "fragmentdb"
p_img = [("FragmentDB", 'h1')]

from PluginTools import PluginTools as PT


class FragmentDB(PT):

  def __init__(self):
    super(FragmentDB, self).__init__()
    self.cite_str = "D. Kratzert, I. Krossing, (2018), J. Appl. Cryst. 51, 928-934. "
    self.cite_str2 = "Kratzert, D., Holstein, J.J. & Krossing, I. (2015). J. Appl. Cryst. 48, 933-938."
    self.p_name = p_name
    self.p_path = p_path
    self.p_scope = p_scope
    self.p_htm = p_htm
    self.p_img = p_img
    self.deal_with_phil(operation='read')
    self.print_version_date()
    self.setup_gui()
    self.params = OV.GuiParams()
    self.dbfile = os.sep.join([self.p_path, "fragment-database.sqlite"])
    self.userdb = 'user-fragment-database.sqlite'
    self.userdbfile = os.sep.join([instance_path, os.sep, 'db', os.sep, self.userdb])
    if not os.path.exists(self.userdbfile) or os.path.getsize(self.userdbfile) < 100:
      initialize_user_db(self.userdbfile)
    # for edited fragments:
    self.cell = []
    print(' FragmentDB version:', FDB_VERSION)

  def get_cell(self):
    """
    returns the cell from the refinement model
    """
    precell = olex_core.GetRefinementModel(False)['aunit']['cell']
    cell = [precell['a'][0], precell['b'][0], precell['c'][0], precell['alpha'][0],
            precell['beta'][0], precell['gamma'][0]]
    return cell

  def init_plugin(self):
    """
    initialize the plugins main form
    """
    try:
      fragid = int(OV.GetParam('FragmentDB.fragment.fragId'))
    except(RuntimeError, ValueError):
      return
    if fragid == 0:
      return
    self.get_resi_class()
    self.toggle_resinum(OV.GetParam('FragmentDB.fragment.resitoggle'))
    self.set_fragment_picture()
    self.display_image('FDBMOLEPIC', 'displayimg.png')
    self.show_reference()
    resinum = self.find_free_residue_num()
    # olx.html.SetValue('RESIDUE', resinum)
    OV.SetParam('FragmentDB.fragment.resinum', resinum)
    # olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    # self.guess_values()

  def set_id(self, fragid=0):
    """
    Sets the fragment id in the phil for the search field
    """
    try:
      int(fragid)
    except(ValueError):
      return
    OV.SetParam("FragmentDB.fragment.fragId", fragid)

  def clear_mainvalues(self):
    """
    clears the state of the main interface
    """
    olx.html.SetImage('FDBMOLEPIC', 'blank.png')
    if olx.fs.Exists('largefdbimg.png') == 'true':
      im = Image.new('RGBA', (1, 1), self.params.html.table_bg_colour.rgb)
      OlexVFS.save_image_to_olex(im, 'largefdbimg.png', 0)
    olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    resinum = self.find_free_residue_num()
    olx.html.SetValue('RESIDUE', resinum)
    OV.SetParam('FragmentDB.fragment.fragId', 0)
    # olx.SetVar('fragment_ID', 0)
    OV.SetParam('FragmentDB.fragment.frag_part', 0)
    olx.html.SetValue('FRAGMENT_PART', 0)
    OV.SetParam('FragmentDB.fragment.frag_fvar', 1)
    olx.html.SetValue('FRAG_FVAR', 1)
    self.set_occu('1.0')
    self.get_fvar_occ()
    OV.SetParam('FragmentDB.fragment.resi_class', "")
    OV.SetParam('FragmentDB.new_fragment.frag_resiclass', "")
    olx.html.SetValue('RESIDUE_CLASS', "")

  def set_occu(self, occ):
    """
    sets the occupancy, even if you enter a comma value instead of point as
    decimal separator.
    """
    if ',' in occ:
      occ = occ.replace(',', '.')
    try:
      # do not throw the warning below in case of empty number:
      if occ == '':
        return
      float(occ)
    except(SyntaxError, NameError, ValueError):
      print('Invalid value for occupancy provided')
      return
    varname = 'FragmentDB.fragment.frag_occ'
    olx.html.SetValue('FRAG_OCCUPANCY', occ)
    OV.SetParam(varname, occ)

  def toggle_resinum(self, resinum):
    """
    :param resinum: decide wether on or off
    :type resinum: bool
    """
    OV.SetParam('FragmentDB.fragment.resitoggle', resinum)

  def set_resiclass(self, resiclass, name):
    """
    sets the residue class and ensures that it is of len(4)
    and .isalpha is the first character.
    """
    varname = ''
    if name.upper() == 'RESIDUE_CLASS'.upper():
      varname = 'FragmentDB.fragment.resi_class'
    if name.upper() == 'Inputfrag.residue'.upper():
      varname = 'FragmentDB.new_fragment.frag_resiclass'
    if resiclass is None:
      resiclass = ''
    if not resiclass:
      return
    if not resiclass[0].isalpha():
      if not resiclass[1:]:
        print('The residue class has to start with a letter.')
        OV.SetParam(varname, '')
        olx.html.SetValue(name, '')
        return
      # resiclass does not start with a char:
      print('The residue class has to start with a letter.')
      OV.SetParam('FragmentDB.fragment.resi_class', resiclass[1:].upper())
      olx.html.SetValue(name, OV.GetParam(varname))
    # force 4 characters:
    elif len(resiclass) > 4:
      OV.SetParam(varname, resiclass[:4].upper())
      olx.html.SetValue(name, OV.GetParam(varname))
    else:
      OV.SetParam(varname, resiclass.upper())
      # olx.html.SetValue(name, OV.GetParam(varname))

  def list_fragments(self):
    """
    returns the available fragments in the database
    the list of names is separated by semicolon
    i[0] => number
    i[1] => name
    """
    items = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in FragmentTable(self.dbfile, self.userdbfile)])
    return items

  def search_fragments(self, search_string):
    """
    performs a search for an unsharp name in a list
    """
    selected_list = ''
    db = FragmentTable(self.dbfile, self.userdbfile)
    if not search_string:
      print("Empty search string.")
      selected_list = self.list_fragments()
      olx.html.SetItems('LIST_FRAGMENTS', selected_list)
      return
    else:
      selected_results = db.find_fragment_by_name(search_string, 8)
      # do nothing in case we searched with the original name (exact one hit):
      if len(selected_results) == 1:
        print('Nothing found...')
        return
      selected_list = ';'.join(['{}<-{}'.format(i[1], i[0]) for i in selected_results])
    # propagate the smaller list to the combo-box:
    olx.html.SetItems('LIST_FRAGMENTS', selected_list)
    # show the first result in combo box and intialize the fragment:
    olx.html.SetValue('LIST_FRAGMENTS', '{}'.format(selected_results[0][1]))
    frag_id = int(selected_results[0][0])
    OV.SetParam('FragmentDB.fragment.fragId', frag_id)
    # olx.SetVar('fragment_ID', frag_id)
    self.init_plugin()

  def format_atoms_for_importfrag(self, atoms):
    """
    format the input atoms to use them with importfarg
    """
    newlist = []
    finallist = []
    if OV.GetParam('FragmentDB.fragment.invert'):
      atoms = invert_atomlist_coordinates(atoms)
    for i in atoms:
      # atoms without sfac:
      newlist.append('{:4.4s} {:>7.4f}  {:>7.4f}  {:>7.4f}'.format(i[0], i[2], i[3], i[4]))
    atoms = '\n'.join(newlist)
    finallist.append('FRAG')
    finallist.append('\n' + atoms)
    finallist.append('\nFEND')
    text = ' '.join(finallist)
    return text

  def onImport(self, atom_ids):
    """
    Function is called when ImportFrag terminates after import. It gets the
    atomIds from ImportFrag via the onFragmentImport callback.

    Essentially, this runs during fragment fit!

    :param atom_ids: list of atoms as olex2 atom ids
    :type atom_ids: str
    """
    print('Imported atom ids: {}'.format(atom_ids))
    atom_ids = atom_ids.split()
    replatoms = None
    if OV.GetParam('FragmentDB.fragment.replace'):
      replatoms = self.find_atoms_to_replace(atom_ids)
    # now replace atoms in a certain distance in part 0:
    if replatoms:
      OV.cmd("sel u")
      OV.cmd("sel #c{}".format(' #c'.join(replatoms)))
      OV.cmd('KILL')
      print('Deleting: {}'.format(' '.join(replatoms)))
    # define the other properties:
    self.define_atom_properties(atom_ids)
    OV.unregisterCallback('onFragmentImport', self.onImport)

  def find_atoms_to_replace(self, frag_atoms, remdist=1.22):
    """
    this method looks around every atom of the fitted fragment and removes
    atoms that are near a certain distance.
    go through all fragment atoms and look for each atom for all atoms in
    part 0 if they are nearer as remdist to them
    :param frag_atoms: atom ids of the fitting fragment ['21', '22', '23', '24']
    :type frag_atoms: list
    :param remdist: distance below atoms shoud be deleted.
    :type remdist: float
    """
    cell = self.get_cell()
    atoms_to_delete = []
    all_atoms_dict = self.get_atoms_list(part=0, notype='')
    frag_crd_dict = {}
    for i in frag_atoms:
      # create fragment dict:
      frag_crd_dict[i] = all_atoms_dict[int(i)]
      # remove fragment atoms from structure:
      all_atoms_dict.pop(int(i), None)
    for aa_id in all_atoms_dict:
      if all_atoms_dict[aa_id][2] == 0:
        for f_id in frag_crd_dict:
          at1 = all_atoms_dict[aa_id][1]  # coordinates
          at2 = frag_crd_dict[f_id][1]  # coordinates
          d = atomic_distance(at1, at2, cell)
          # now the atoms inside the remdist go into deltion list
          if d < remdist:
            atoms_to_delete.append(str(aa_id))
    return list(set(atoms_to_delete))

  def insert_frag_with_ImportFrag(self, fragId, part=-1, occ=1, afix=False):
    """
    Input a fragment with ImportFrag
    :param afix: rigid afix group or not
    :type afix: bool
    :param fragId: FragmentId
    :type fragId: int
    :param part: SHELX part
    :type part: int
    :param occ: occupancy
    :param dfix: generate dfix restraints or not after fit.
    """
    # this callback runs in the moment when ImportFrag is finished. onImport
    # then defines the further properties of the fragment:
    OV.registerCallback('onFragmentImport', self.onImport)
    fragpath = os.sep.join([OV.StrDir(), 'fragment.txt'])
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      atoms = self.format_atoms_for_importfrag([i for i in db[fragId]])
    except(IndexError):
      return
    try:
      with open(fragpath, 'w') as f:
        f.write(atoms)
    except:
      print('Unable to write fragment to temporary file. Not fit possible!')
      return
    OV.cmd("file")
    if afix:
      afix = "6"
    else:
      afix = '0'
    if OV.GetParam('FragmentDB.fragment.use_dfix'):
      print('Applying DFIX restraints ...')
      olx.ImportFrag(fragpath, a=afix, p=part, o=occ, d=True)
      print('Finished.')
      # onImport() runs after ImportFrag
    else:
      olx.ImportFrag(fragpath, a=afix, p=part, o=occ, d=False)
    return

  def fit_db_fragment(self):
    """
    fit a molecular fragment from the database into olex2
    """
    try:
      fragId = int(OV.GetParam('FragmentDB.fragment.fragId'))
    except(RuntimeError, ValueError):
      # no fragment chosen-> do nothing
      return False
    OV.cmd("labels false")
    occupancy = OV.GetParam('FragmentDB.fragment.frag_occ')
    afix = False
    if OV.GetParam('FragmentDB.fragment.rigid'):
      afix = True
    # alway use "part -1" to prevent atom deletion:
    if fragId != 0 and fragId is not None:
      self.insert_frag_with_ImportFrag(fragId, part=-1, occ=occupancy, afix=afix)
    else:
      print('Please select a fragment first, or type text and hit Enter key to search.')
      return False
    gui.report.add_to_citation_list(self.cite_str)
    gui.report.add_to_citation_list(self.cite_str2)
    return True

  def atomrenamer(self, labeldict):
    """
    Renames atoms according to the database names after ImportFrag.
    :param labeldict: dictionary with names and ids of the atoms
    :type labeldict: dict
    """
    for at in labeldict:
      olx.Name('#c' + labeldict[at.upper()], at)

  def define_atom_properties(self, atomids, fragId=None):
    """
    Defines the atoms properties of the fitted fragment after ImportFrag
    :param atomids: atomic olex2 ids of the atoms
    :type atomids: list
    :param fragId: fragment id
    :type fragId: int
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    resiclass = OV.GetParam('FragmentDB.fragment.resi_class')
    freevar = int(OV.GetParam('FragmentDB.fragment.frag_fvar'))
    partnum = OV.GetParam('FragmentDB.fragment.frag_part')
    if not OV.GetParam('FragmentDB.fragment.use_dfix'):
      print('Applying fragment properties:')
    if not fragId:
      try:
        fragId = int(OV.GetParam('FragmentDB.fragment.fragId'))  # olx.GetVar('fragment_ID')
      except(RuntimeError, ValueError):
        # no fragment chosen-> do nothing
        return
    labeldict = OrderedDict()
    atomids = [str(x) for x in atomids]
    dbatom_names = [i[0].upper() for i in db[fragId]]
    for at_id, name in zip(atomids, dbatom_names):
      labeldict[name.upper()] = at_id
    # select all atomids to do the fit:
    if freevar != 1:
      OV.cmd("sel #c{}".format(' #c'.join(atomids)))
      OV.cmd("fvar {}".format(freevar))
    resi_on = bool(OV.GetParam('FragmentDB.fragment.resitoggle'))
    resinum = 0
    if resi_on:
      resinum = self.find_free_residue_num()
      if resiclass and resinum or not resiclass and resinum:
        self.make_residue(atomids, resiclass, resinum)
        self.atomrenamer(labeldict)
    # Placing restraints:
    if not OV.GetParam('FragmentDB.fragment.use_dfix') and not OV.GetParam('FragmentDB.fragment.roff'):
      self.make_restraints(labeldict, fragId, resinum, resiclass)
    self.make_part(atomids, partnum)
    return atomids

  def make_part(self, atoms, partnum):
    """
    Assign part number to a fragment
    :param atoms: list of atoms
    :type atoms: list
    :param partnum: SHELX part number
    :type partnum: int
    """
    OV.cmd("sel #c{}".format(' #c'.join(atoms)))
    OV.cmd("PART {}".format(partnum))

  def make_residue(self, atoms, resiclass, resinum):
    """
    selects the atoms and applies "RESI class number" to them
    """
    if resiclass:
      OV.cmd("sel #c{}".format(' #c'.join(atoms)))
      OV.cmd("RESI {} {}".format(resiclass, resinum))
      # Applies a soft SIMU to all neughbouring atoms (not the bonded), because SHELXL behaves different
      # for residues.
      # In SHELXL, SIMU_NAME At1 > At5 applies only to each residue itself, not near neighbours of
      # other residues!
      OV.cmd('SIMU 0.04 0.08 1')
    else:
      OV.cmd("sel #c{}".format(' #c'.join(atoms)))
      OV.cmd("RESI {}".format(resinum))

  def make_restraints(self, labeldict, fragId, resinum=0, resiclass=''):
    """
    Applies restraints to atoms.
    :param labeldict: a dictionary where the keys are atom names and the
                      values are Olex2 atom Ids
    :type labeldict: dictionary
    :param fragId: the database Id of the fragment
    :type fragId: integer
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    restraints = db.get_restraints(fragId)
    if not restraints:
      return
    for num, i in enumerate(restraints):
      # i[0] is restraint like SADI or DFIX
      # i[1] is a string of atoms like 'C1 C2'
      restraint_atoms = i[1].upper()
      if '>' in restraint_atoms or '<' in restraint_atoms:
        restraint_atoms = self.range_resolver(restraint_atoms.split(), labeldict.keys())
      line = []
      for at in restraint_atoms.split():
        # is it a potential atom (starts with alphabetic character):
        if at[0].isalpha():
          try:
            line.append('#c' + labeldict[at.upper()])
          except KeyError:
            # in this case, an atom name in the restraint does not
            # exist in the fragments atom list!
            print('\nUnknown restraint or atom found in restraints line {}.\n'.format(num + 1))
            # I must exit here!
            return
        else:
          line.append(at)
      # Applies the implicit restraint to atoms in line:
      # I disabled implicit restraints because they cause trouble during atom
      # renaming. Update: It seems to work now, enabling again.
      if i[0] in helper_functions.IMPL_RESTRAINT_CARDS and resinum != 0 and resiclass:
        OV.cmd("{} -i {}".format(i[0], ' '.join(line)))
      else:
        # applies direct restraints:
        OV.cmd("{} {}".format(i[0], ' '.join(line)))

  def prepare_picture(self, im, max_size=120, ratiolim=0.6):
    """
    resizes and colorizes the picture to diplay it in olex2
    needs a PIL Image instance
    :param im: PIL Image instance
    :type im: PIL Image instance
    :param max_size: maximum size of the picture in pixels
    :type max_size: int
    :param control: html control to get the background color right
    :type control: string
    :param ratiolim: limits the resize for small fragments
    :type ratiolim: float
    """
    for i in im.size:
      if i < 50:
        print('This picture is too small.')
    im = im.convert(mode="RGBA")
    img_w, img_h = im.size
    ratio = abs(float(max_size) / float(max(im.size)))
    # just an empirical value:
    if float(max_size) / float(max(im.size)) > ratiolim:
      ratio = ratiolim
      # resize equally to fit in max_size
    im = im.resize((int(img_w * ratio), int(img_h * ratio)), Image.ANTIALIAS)
    # empty image of max_size
    bgcolor = self.params.html.table_bg_colour.rgb
    # bgcolor = self.params.html.bg_colour.rgb
    IM = Image.new('RGBA', (max_size, max_size), bgcolor)
    bg_w, bg_h = IM.size
    img_w, img_h = im.size
    # offset for image placement
    offset = (bg_w - img_w) / 2, (bg_h - img_h) / 2
    # place image in center of background image:
    IM.paste(im, offset)
    return IM

  def set_fragment_picture(self, max_size=120):
    """
    displays a picture of the fragment from the database in Olex2
    :param max_size: maximum size of the picture in pixels
    :type max_size: int
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    max_size = int(max_size)
    fragId = OV.GetParam('FragmentDB.fragment.fragId')  # olx.GetVar('fragment_ID')
    pic = db.get_picture(fragId)
    if not pic:
      print('No fragment picture found.')
      return False
    imo = Image.open(StringIO.StringIO(pic))
    # save it as raw and small pic:
    OlexVFS.save_image_to_olex(imo, 'storepic.png', 0)
    im = self.prepare_picture(imo, max_size)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)
    iml = self.prepare_picture(imo, max_size=450, ratiolim=1.0)
    OlexVFS.save_image_to_olex(iml, 'largefdbimg.png', 0)

  def display_image(self, zimgname, image_file):
    """
    display an image in zimgname
    """
    try:
      olx.html.SetImage(zimgname, image_file)
    except RuntimeError:
      pass

  def store_picture(self, title, ffilter, location, default_name=''):
    """
    opens a file dialog and stores the selected picture in the olex-vfs
    :param title: Title of the file dialog
    :type title: string
    :param filter: file endings filter like '*.png; *.tiff; *.tif; *.jpg; *.jpeg'
    :type filter: sting
    :param location: location of the path of file dialog
    :type location: string
    :param default_name: default preselected file name for the dialog
    :type default_name: string
    """
    picfile = olx.FileOpen(title, ffilter, location, default_name)
    if not picfile:
      return
    fsize = os.stat(picfile).st_size
    # do not allow picture of more than 1MB (1000000 bytes):
    if fsize > 1000000:
      print('This file is too large!')
      return
    if fsize <= 1:
      print('No valid picture found!')
      return
    im = Image.open(picfile)
    # TODO: better use set_fragment_picture() here
    OlexVFS.save_image_to_olex(im, 'storepic.png', 0)
    iml = self.prepare_picture(im, max_size=450, ratiolim=1.0)
    OlexVFS.save_image_to_olex(iml, 'largefdbimg.png', 0)
    # display it.
    im = self.prepare_picture(im)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)
    olx.html.SetImage('Inputfrag.MOLEPIC2', 'displayimg.png')

  def range_resolver(self, restraintat, atom_names):
    """
    resolves the atom names of ranges like "C1 > C5"
    works for each restraint line separately.
    :param restraintat: atoms with a range definition
    :type restraintat: list
    :param atom_names: names of atoms in the fragment
    :type atom_names: list
    """
    # dict with lists of positions of the > or < sign:
    rightleft = {'>': [], '<': []}
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
          left = atom_names.index(restraintat[i - 1]) + 1
          right = atom_names.index(restraintat[i + 1])
          restraintat[i:i + 1] = atom_names[left:right]
        else:
          # backward range
          left = atom_names.index(restraintat[i - 1])
          right = atom_names.index(restraintat[i + 1]) + 1
          names = atom_names[right:left]
          names.reverse()  # counting backwards
          restraintat[i:i + 1] = names
    return ' '.join(restraintat)

  def find_free_residue_num(self):
    """
    Determines which residue number is unused.
    Either finds the first unused residue in the list of residue numbers
    or returns the last used + 1.
    :rtype: int
    """
    residues = self.get_residue_numbers()
    # find unused numbers in the list:
    resi = False
    try:
      for i, j in zip(residues, range(1, residues[-1] + 1)):
        # gap in list found:
        if i != j:
          resi = j
          break
    except IndexError:
      return 1
    # no gap, thus use next number:
    if not resi:
      try:
        resi = residues[-1] + 1
      except IndexError:
        # no residue at all in the structure
        return 1
    return resi

  def get_residue_numbers(self):
    """
    returns a list of residue numbers in the structure
    :return: list of residue number
    :rtype: list
    """
    residues = []
    # get the properties of the atoms:
    for resi in olex_core.GetRefinementModel(False)['aunit']['residues']:
      try:
        # atoms in residue 0 have no 'number'
        residues.append(resi['number'])
      except:
        pass
    residues.sort()
    return residues

  def get_resi_class(self):
    """
    sets the residue class from the respective database fragment.
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      fragId = int(OV.GetParam('FragmentDB.fragment.fragId'))  # olx.GetVar('fragment_ID')
    except(RuntimeError, ValueError) as e:
      print('Unable to get fragment Id', e)
      return
    resiclass = db.get_residue_class(int(fragId))
    OV.SetParam('FragmentDB.fragment.resi_class', resiclass)
    OV.SetParam('FragmentDB.new_fragment.frag_resiclass', resiclass)
    # set the class in the text field of the gui:
    olx.html.SetValue('RESIDUE_CLASS', resiclass.upper())

  def show_reference(self, edit=False):
    """
    show the reference of a fragment in the GUI
    :type edit: bool
    """
    try:
      fragId = OV.GetParam('FragmentDB.fragment.fragId')  # olx.GetVar('fragment_ID')
    except(RuntimeError, ValueError):
      return
    db = FragmentTable(self.dbfile, self.userdbfile)
    ref = db.get_reference(fragId)
    if not edit:
      pass
      # disabled, because replace checkbox is there
      # olx.html.SetValue('REFERENCE', ref)
    else:
      # reference of the edit window
      olx.html.SetValue('REFERENCE_edit', ref)

  def get_selected_atom_names(self):
    """
    returns the names of the currently selected atoms.
    Hydrogen atoms are discarded.
    """
    atoms_all = olex.f("sel(a)").split()
    if not atoms_all:
      print('No atoms selected!')
      return []
    return atoms_all

  def prepare_selected_atoms(self):
    """
    prepares atoms for FragmentDB.new_fragment.frag_atoms
    atom field in inputfrag.
    :rtype: list
    """
    atoms_all = self.get_selected_atom_names()
    if len(atoms_all) < 3:
      print('Please select at least three atoms.')
      return
    atlist = []
    atoms = []
    # make a list without H atoms:
    for at in atoms_all:
      if at[0] == 'H' and at[:2] not in ('He', 'Hf', 'Ho', 'Hg'):
        continue
      atoms.append(at)
    atoms.sort()
    # natsort(atoms)
    crd = [olx.Crd(x) for x in atoms]
    # now I want to remove the residue number:
    atoms = [y.split('_')[0] for y in atoms]
    for atom, coord in zip(atoms, crd):
      atlist.append('{:4.4s} {:>8.6s} {:>8.6s} {:>8.6s}'.format(atom, *coord.split()))
    at = ' \n'.join(atlist)
    olx.html.SetValue('Inputfrag.SET_ATOM', at)
    self.cell = '1  1  1  90  90  90'
    olx.html.SetValue('Inputfrag.set_cell', self.cell)
    OV.SetParam('FragmentDB.new_fragment.frag_cell', self.cell)
    OV.SetParam('FragmentDB.new_fragment.frag_atoms', at)
    return atlist

  def open_edit_fragment_window(self):
    """
    opens a new window to input/update a database fragment
    """
    blank = False
    try:
      fragId = OV.GetParam('FragmentDB.fragment.fragId')
    except RuntimeError:
      blank = True
    pop_name = "Inputfrag"
    screen_height = int(olx.GetWindowSize('gl').split(',')[3])
    screen_width = int(olx.GetWindowSize('gl').split(',')[2])
    box_x = int(screen_width * 0.2)
    box_y = int(screen_height * 0.1)
    width, height = 550, 710
    path = "{}/inputfrag.htm".format(self.p_path)
    olx.Popup(pop_name, path, b="tcrp", t="Create/Edit Fragments", w=width,
              h=height, x=box_x, y=box_y)
    frag = self.get_frag_for_gui()
    if frag:
      self.display_image('Inputfrag.MOLEPIC2', 'displayimg.png')
    if blank:
      self.display_image('Inputfrag.MOLEPIC2', 'blank.png')

  def display_large_image(self):
    """
    display the current fragment image in large to be able to read the labels
    """
    if olx.fs.Exists("largefdbimg.png") == 'false':
      print('No larger picture available.')
      return
    pop_name = "View Fragment"
    screen_height = int(olx.GetWindowSize('gl').split(',')[3])
    screen_width = int(olx.GetWindowSize('gl').split(',')[2])
    box_x = int(screen_width * 0.1)
    box_y = int(screen_height * 0.1)
    width, height = 500, 520
    html = """
    <a target="" href="spy.FragmentDB.save_picture()">
    <zimg name="LMOLEPIC" 
        border="0" 
        src="largefdbimg.png" 
        height=450 
        width=450 
        align="center">
    """
    OV.write_to_olex('large_fdb_image.htm', html)
    olx.Popup(pop_name, "large_fdb_image.htm", b="tcrp", t="View Fragment", w=width,
              h=height, x=box_x, y=box_y)

  def save_picture(self):
    """
    save the enlarged picture to a file
    """
    title = 'Save Pictue'
    ffilter = '*.png;*.PNG'
    location = ''
    default_name = 'molecule.png'
    img_name = olx.FileSave(title, ffilter, location, default_name)
    if not img_name:
      return
    olx.fs.Dump("largefdbimg.png", img_name)

  def set_frag_name(self, enable_check=True):
    """
    handles the name of a new/edited fragment
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    fragname = OV.GetParam('FragmentDB.new_fragment.frag_name')
    if fragname == '':
      return False
    if db.has_exact_name(fragname) and enable_check:
      print('\n{} is already in the database. \nPlease choose a different name.\n'.format(fragname))
      return False
    else:
      OV.SetParam('FragmentDB.new_fragment.frag_name', fragname)
      return True

  def set_frag_cell(self):
    """
    set the unit cell of a new fragment to convert its coordinates to cartesian
    the resulting cell is retrieved from phil and stored in self.frag_cell
    :rtype: bool
    """
    self.frag_cell = ''
    try:
      frag_cell = olx.html.GetValue('Inputfrag.set_cell').split()
    except(AttributeError):
      print('No unit cell defined. You nedd to define unit cell parameters!')
      return
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
    """
    -handles the fragment atoms of a new/edited fragment
    -returns a list of lists:
      [['C4', '0.282212', '0.368636', '0.575493'], ...]
    :rtype: list
    """
    atoms = OV.GetParam('FragmentDB.new_fragment.frag_atoms')
    try:
      atoms = atoms.split('\n')
      # remove empty lines:
      atoms = [i for i in atoms if i]
      atoms = [i.split() for i in atoms]
    except AttributeError:
      atoms = None
      return
    finalatoms = self.atoms_parser(atoms)
    return finalatoms

  def atoms_parser(self, atoms):
    """
    formats the atoms from shelx as long list to a list of list with
    exactly four items: Atom  x y z
    :param atoms: line of atoms
    :type atoms: list
    :rtype list
    """
    # go through all atoms and cut their line to 4 list elements At  x y z:
    for num, line in enumerate(atoms):
      # line with sfac
      if len(line) > 4 and len(line[1]) < 3:
        atoms[num] = [line[0]] + line[2:5]
      # line without sfac, but long
      if len(line) > 4 and len(line[1]) > 3:
        atoms[num] = line[:4]
      if len(line) == 4:  # name and x,y,z
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
          num = num + 1
          line[0] = line[0] + str(num)
    return atoms

  def set_frag_restraints(self):
    """
    handles the restraints of a new/edited fragment
    :rtype: list
    """
    restraints = []  # @UnusedVariable
    restr = OV.GetParam('FragmentDB.new_fragment.frag_restraints')
    try:
      restr = restr.split('\n')
      # remove empty lines:
      restr = [i.upper() for i in restr if i]
      restr = [i.split() for i in restr]
    except AttributeError as e:  # @UnusedVariable
      # print(e)
      restr = None
      return
    return restr

  def prepare_residue_class(self):
    """
    set the residue class of a new fragment
    :rtype: str
    """
    resi_class = OV.GetParam('FragmentDB.new_fragment.frag_resiclass')
    if resi_class:
      return resi_class

  def prepare_reference(self):
    """
    prepare the reference for a new/edited fragment
    :rtype: str
    """
    reference = OV.GetParam('FragmentDB.new_fragment.frag_reference')
    if reference:
      return reference

  def prepare_atoms_list(self, fragId, element=False, as_html=False):
    """
    prepare the atom list to display in a multiline edit field
    :type fragId: int
    :type element: bool
    :type as_html: bool
    """
    atlist = []
    atoms_list = []
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      atoms_list = db[fragId]
    except IndexError:
      return False
      # print('Database Fragment {} not found.'.format(fragId))
    if not atoms_list:
      return False
    atoms_list = [[i for i in y] for y in atoms_list]
    for num, i in enumerate(atoms_list):
      try:
        if element:
          # Have to remove this, because a table is bad to copy&paste:
          # if as_html:
          #  atlist.append('<tr width=60%><td align=left  width=20%>{:4.4s}&nbsp;</td>  '
          #                    '<td align=right>1 &nbsp;&nbsp;</td> '
          #                    '<td align=right>{:>10.4f}&nbsp;&nbsp;</td> '
          #                    '<td align=right>{:>10.4f}&nbsp;&nbsp;</td> '
          #                    '<td align=right>{:>10.4f}&nbsp;&nbsp;</td></tr>'.format(i[0], i[2], i[3], i[4]))
          # else:
          atlist.append('{:4.4s}  1  {:>8.4f} {:>8.4f} {:>8.4f}'.format(i[0], i[2], i[3], i[4]))
        else:
          atlist.append('{:4.4s} {:>8.4f} {:>8.4f} {:>8.4f}'.format(i[0], i[2], i[3], i[4]))
      except UnicodeEncodeError:
        print('Invalid atomline found. Non-ASCII character in line {}.'.format(num + 1))
    if not as_html:
      at = ' \n'.join(atlist)
    else:
      at = '<br>'.join(atlist)
    return at

  def prepare_fragname(self, fragId):
    """
    prepare the fragment name to display in a multiline edit field
    :type fragId: int
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      name = str(db.get_fragment_name(fragId)[0])
    except:
      print('Could not get a name from the database.')
      return False
    return name

  def prepare_restraints(self, fragId):
    """
    prepare the fragment restraints to display in a multiline edit field
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    restr_list = db.get_restraints(fragId)
    if not restr_list:
      return False
    restr_list = [[str(i) for i in y] for y in restr_list]
    restr = '\n'.join(['  '.join(i) for i in restr_list])
    return restr

  def add_new_frag(self):
    """
    execute this to add a new fragment
    """
    # add_new_frag must check the restraints vor validity!
    state = self.set_frag_name()
    if not state:
      print('Please define a fragment name. Fragment is not stored!')
      return
    if not self.set_frag_cell():
      return False
    atoms = self.set_frag_atoms()
    if not atoms:
      print('Please add atoms!!!')
      return
    restraints = self.set_frag_restraints()
    resiclass = self.prepare_residue_class()
    reference = self.prepare_reference()
    self.store_new_fragment(atoms, restraints, resiclass, reference)

  def update_fragment(self):
    """
    execute this to update a fragment
    updates the database information of a fragment
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    fragname = OV.GetParam('FragmentDB.new_fragment.frag_name')
    state = self.set_frag_name(enable_check=False)
    if not state:
      return
    if not self.set_frag_cell():
      return False
    atlines = self.set_frag_atoms()
    restraints = self.set_frag_restraints()
    resiclass = self.prepare_residue_class()
    reference = self.prepare_reference()
    try:
      pic_data = OlexVFS.read_from_olex('storepic.png')
    except TypeError:
      print('No picture found')
      pic_data = ''
    coords = self.prepare_coords_for_storage(atlines)
    if not check_restraints_consistency(restraints, atlines, fragname):
      print('\nFragment was not added to the database!')
      return
    if restraints:
      helper_functions.check_sadi_consistence(atlines, restraints, self.frag_cell,
                                              fragname)
    self.delete_fragment(reset=False)
    frag_id = db.store_fragment(fragname, coords, resiclass, restraints,
                                reference, picture=pic_data)
    print('Updated fragment "{0}".'.format(fragname))
    if frag_id:
      olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
      olx.SetVar('fragment_ID', frag_id)
    else:
      print('Something went wrong during fragment storage.')
    self.get_frag_for_gui()
    self.set_fragment_picture()
    olx.html.SetValue('RESIDUE_CLASS', '')
    self.show_reference()
    olx.html.SetImage('FDBMOLEPIC', 'blank.png')

  def prepare_coords_for_storage(self, atlines):
    """
    transform the string formated atoms into a list of atoms lists

    :param atlines: atoms from the input field
    :type atlines: list
    :rtype: list
    """
    coords = []
    for line in atlines:
      try:
        frac_coord = [float(i) for i in line[1:4]]
      except ValueError:
        print('Invalid coordinate defined in line "{}"'.format(' '.join(line)))
        return
      if len(frac_coord) < 3:
        print('Coordinate value missing in "{}".'.format(' '.join(line)))
        continue
      coord = frac_to_cart(frac_coord, self.frag_cell)
      line[1:4] = coord
      coords.append(line)
    return coords

  def get_frag_for_gui(self):
    """
    get the fragment to display in the multiline edit field
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      fragId = OV.GetParam('FragmentDB.fragment.fragId')
    except(RuntimeError):
      return
    at = self.prepare_atoms_list(fragId)
    if not at:
      return False
    name = self.prepare_fragname(fragId)
    if name == False:
      name = "Could not get a fragment name from the database"
    restr = self.prepare_restraints(fragId)
    residue = self.prepare_residue_class()
    reference = db.get_reference(fragId)
    cell = '1  1  1  90  90  90'
    olx.html.SetValue('Inputfrag.SET_ATOM', at)
    olx.html.SetValue('Inputfrag.set_cell', cell)
    olx.html.SetValue('Inputfrag.set_name', name)
    olx.html.SetValue('Inputfrag.restraints', restr)
    olx.html.SetValue('Inputfrag.residue', residue)
    olx.html.SetValue('Inputfrag.reference_ed', reference)
    OV.SetParam('FragmentDB.new_fragment.frag_name', name)
    OV.SetParam('FragmentDB.new_fragment.frag_atoms', at)
    OV.SetParam('FragmentDB.new_fragment.frag_cell', cell)
    OV.SetParam('FragmentDB.new_fragment.frag_restraints', restr)
    OV.SetParam('FragmentDB.new_fragment.frag_resiclass', residue)
    OV.SetParam('FragmentDB.new_fragment.frag_reference', reference)
    return True

  def store_new_fragment(self, atlines, restraints, resiclass, reference):
    """
    add a new fragment to the database
    :type atlines: list
    :type restraints: list
    :type resiclass: str
    :type reference: str
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    # check if the given name already exist in the database
    # store fragment with a new number
    fragname = OV.GetParam('FragmentDB.new_fragment.frag_name')
    self.set_frag_cell()
    try:
      pic_data = OlexVFS.read_from_olex('storepic.png')
    except TypeError:
      pic_data = ''
    coords = self.prepare_coords_for_storage(atlines)
    print('Adding fragment "{0}" to the database.'.format(fragname))
    if not check_restraints_consistency(restraints, atlines, fragname):
      print('Fragment was not added to the database!')
      # self.blink_field('Inputfrag.restraints')
      # OV.Alert('Invalid restraint', 'One of the restraints is invalid. \nNo changes to the database were performed.', 'OK')
      return
    if restraints:
      helper_functions.check_sadi_consistence(atlines, restraints, self.frag_cell,
                                              fragname)
    frag_id = db.store_fragment(fragname, coords, resiclass, restraints,
                                reference, picture=pic_data)
    if frag_id:
      olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    else:
      print('Something went wrong during fragment storage.')
    # now get the fragment back from the db to display the new cell:
    olx.SetVar('fragment_ID', frag_id)
    self.init_plugin()
    self.clear_mainvalues()

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
    if olx.fs.Exists('displayimg.png') == 'true':
      OV.CopyVFSFile('blank.png', 'displayimg.png')
    if olx.fs.Exists('largefdbimg.png') == 'true':
      im = Image.new('RGBA', (1, 1), self.params.html.table_bg_colour.rgb)
      OlexVFS.save_image_to_olex(im, 'largefdbimg.png', 0)

  def delete_fragment(self, reset=True):
    """
    deletes a fragment from the database
    # Todo: reset all fields (oic, atoms, ...) after deltion
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    try:
      fragId = OV.GetParam('FragmentDB.fragment.fragId')
    except(RuntimeError):
      print('Could not delete a fragment. No fragment is selected.')
      return
    del db[fragId]
    olx.html.SetItems('LIST_FRAGMENTS', self.list_fragments())
    olx.SetVar('fragment_ID', '')
    # Now delete the fields:
    if reset:
      self.blank_state()

  def get_chemdrawstyle(self):
    """
    opens a file dialog window to save the standard chemdraw style
    """
    from shutil import copyfile
    title = "Save Chemdrawstyle file"
    ffilter = ''
    location = '.'
    stylename = None  # @UnusedVariable
    default_name = 'chemdraw_style.cds'
    stylename = olx.FileSave(title, ffilter, location, default_name)
    if not stylename:
      return
    spath = "{}/drawstyle.cds".format(self.p_path)
    copyfile(spath, stylename)

  def trim(self, im):
    """
    Trims a given PIL picture to remove whitespace
    :type im: Image
    :rtype: Image
    """
    bg = Image.new(im.mode, im.size, im.getpixel((0, 0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
      return im.crop(bbox)
    else:
      print("Trim of image failed.")
      return im

  def make_selctions_picture(self):
    """
    creates a picture from the currently selected fragment in Olex2 and
    stores it in 'storepic.png' as well as 'displaypic.png'
    """
    if not olex.f("sel()").split():
      return
    if len(olex.f("sel()").split()) < 3:
      # print('Please select at least three atoms.')
      return
    picfile = "fdb_tmp.png"
    # OV.cmd('save model "fragdb"') # does not work!
    gui.maps.mu.MapView('off')
    OV.cmd("Legend False")
    OV.cmd("sel atom bonds")
    OV.cmd("ShowQ a false")
    OV.cmd("labels false")
    OV.cmd("showh a False")
    OV.cmd("sel -i")
    OV.cmd("hide")
    OV.cmd("center")
    OV.cmd("sel -i")
    OV.cmd("mpln -n")
    OV.cmd("sel -a")
    OV.cmd("sel bonds -u")
    OV.cmd('label -a')
    OV.cmd('pers')
    OV.cmd("pict fdb_tmp.png -nbg")
    OV.cmd('kill labels')
    OV.cmd("showP -m")
    OV.cmd("showh a True")
    OV.cmd('telp 50')
    OV.cmd("legend true")
    im = Image.open(picfile)
    im = self.trim(im)
    # im = IT.trim_image(im)  # works not so well as trim above
    OlexVFS.save_image_to_olex(im, 'storepic.png', 0)
    iml = self.prepare_picture(im, max_size=450, ratiolim=1.0)
    OlexVFS.save_image_to_olex(iml, 'largefdbimg.png', 0)
    # display it.
    im = self.prepare_picture(im)
    OlexVFS.save_image_to_olex(im, 'displayimg.png', 0)
    olx.html.SetImage('Inputfrag.MOLEPIC2', 'displayimg.png')
    try:
      os.remove(picfile)
    except:
      pass

  def get_fvar_occ(self):
    var = OV.GetParam('FragmentDB.fragment.frag_fvar')
    occ = OV.GetParam('FragmentDB.fragment.frag_occ')
    if float(var) < 0:
      fvar = -(float(abs(var)) * 10 + float(occ))
    else:
      fvar = float(var) * 10 + float(occ)
    olx.html.SetValue('FVAROCC', fvar)
    OV.SetParam('FragmentDB.fragment.fvarocc', fvar)
    return fvar

  def get_atoms_list(self, part=None, notype=''):
    """
    returns all atoms of the crystal
    if part is defined, only collect atoms of defined part
    {0: [u'F1', (0.16726, 0.42638, -0.23772), 0, 0, 'F'],
     id: [u'label', (x, y, z), part, resinum, type], ...}
    """
    model = olex_core.GetRefinementModel(False)
    asym_unit = model['aunit']
    atoms = {}
    for residue in asym_unit['residues']:
      for atom in residue['atoms']:
        if atom['type'] == notype:
          continue
        if part:
          if atom['part'] != part:
            continue
        try:
          resnum = residue['number']
        except KeyError:
          resnum = 0
        atoms[atom['aunit_id']] = [atom['label'], atom['crd'][0], atom['part'], resnum, atom['type']]
    return atoms

  def exportfrag(self):
    """
    print the fragment details on screen for DSR

    Type "spy.FragmentDB.exportfrag" to export the current fragment.
    """
    db = FragmentTable(self.dbfile, self.userdbfile)
    fragtext = []
    htm = ' '
    try:
      fragId = OV.GetParam('FragmentDB.fragment.fragId')
    except RuntimeError:
      return False
    if fragId == 0:
      print("No fragment selected!")
      return False
    pop_name = "Inputfrag"
    screen_height = int(olx.GetWindowSize('gl').split(',')[3])
    screen_width = int(olx.GetWindowSize('gl').split(',')[2])
    box_x = int(screen_width * 0.2)
    box_y = int(screen_height * 0.1)
    width, height = 650, 710
    name = self.prepare_fragname(fragId)
    if name == False:
      return False
    restr = self.prepare_restraints(fragId)
    residue = self.prepare_residue_class()
    reference = db.get_reference(fragId)
    cell = 'FRAG 17 1 1 1 90 90 90'
    fragtext.append('<font face="courier" size="11pt"> ')
    fragtext.append('REM Name: {}'.format(name))
    fragtext.append('REM Src: {}'.format(reference))
    fragtext.append(' <br>'.join(restr.split('\n')))
    fragtext.append('RESI {}'.format(residue))
    fragtext.append(cell)
    at = self.prepare_atoms_list(fragId, element=True, as_html=True)
    if not at:
      print('Not atoms found in fragment!')
      return
    fragtext.append("<table cellspacing=2 cellpadding=1>" + at + "</table>")
    fragtext.append('</font>')
    fragtext.append('<br><b> Please send new fragments to Daniel Kratzert (dkratzert@gmx.de) '
                    '<br> if you whish them to be in the FragmentDB database.</b>')
    htm = ' <br>\n'.join(fragtext)
    OV.write_to_olex('exportfragPop.htm', htm)
    olx.Popup(pop_name, 'exportfragPop.htm', b="tcrp", t="{}".format(name), w=width, h=height, x=box_x, y=box_y)

  def imagedisp(self, name, height=120):
    """
    todo:
    make an object for each picture place with a clear state for every
    state of the plugin.
    """
    if olx.fs.Exists('displayimg.png') == 'false':
      imgname = 'blank.png'
    else:
      imgname = 'displayimg.png'
    html = '''
      <zimg name="{}" 
        border="0" 
        src="{}"
        height={} 
        width=120 
        align="center">
        '''.format(name, imgname, height)
    return html

  def revert_last(self):
    """
    Revertes last structure model.
    """
    pass

  def det_refmodel(self):
    """
    For the graph:
    - each node the atom id
    - the weight is the distance
    - additional property are the coordinates and atom type
    - keep it open for more properties
    {adp: ((0.040012128291554504,
            0.028449999999999993,
            0.018489999999999993,
            0.0038800000000000006,
            -0.0008756170278833539,
            -0.003693364226525346),
           (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)),
   'aunit_id': 11,
   'charge': 0,
   'crd': ((0.268021, 0.476066, 0.394047), (0.0, 0.0, 0.0)),
   'label': u'C13,
   'neighbours': (6, 12, 13, 14),
   'occu': (1.0, 0.0),
   'part': 0,
   'tag': 11,
   'type': u'C}
    """
    model = olex_core.GetRefinementModel(True)
    asym_unit = model['aunit']
    atoms = {}
    for residue in asym_unit['residues']:
      for atom in residue['atoms']:
        pprint.pprint(atom)
        # atoms[atom['aunit_id']] = [ atom['label'], atom['crd'][0], atom['part'], resnum, atom['type'] ]
    return atoms


fdb = FragmentDB()
ref = Refmod()
OV.registerFunction(fdb.det_refmodel, False, "FragmentDB")

OV.registerFunction(fdb.set_id, False, "FragmentDB")
OV.registerFunction(fdb.imagedisp, False, "FragmentDB")
OV.registerFunction(fdb.prepare_selected_atoms, False, "FragmentDB")
OV.registerFunction(fdb.exportfrag, False, "FragmentDB")
OV.registerFunction(fdb.init_plugin, False, "FragmentDB")
OV.registerFunction(fdb.get_fvar_occ, False, "FragmentDB")
OV.registerFunction(fdb.search_fragments, False, "FragmentDB")
OV.registerFunction(fdb.show_reference, False, "FragmentDB")
OV.registerFunction(fdb.make_selctions_picture, False, "FragmentDB")
OV.registerFunction(fdb.set_frag_atoms, False, "FragmentDB")
OV.registerFunction(fdb.open_edit_fragment_window, False, "FragmentDB")
OV.registerFunction(fdb.list_fragments, False, "FragmentDB")
OV.registerFunction(fdb.fit_db_fragment, False, "FragmentDB")
OV.registerFunction(fdb.get_resi_class, False, "FragmentDB")
OV.registerFunction(fdb.find_free_residue_num, False, "FragmentDB")
OV.registerFunction(fdb.get_frag_for_gui, False, "FragmentDB")
OV.registerFunction(fdb.set_occu, False, "FragmentDB")
OV.registerFunction(fdb.set_resiclass, False, "FragmentDB")
OV.registerFunction(fdb.toggle_resinum, False, "FragmentDB")
OV.registerFunction(fdb.store_new_fragment, False, "FragmentDB")
OV.registerFunction(fdb.set_fragment_picture, False, "FragmentDB")
OV.registerFunction(fdb.get_chemdrawstyle, False, "FragmentDB")
OV.registerFunction(fdb.add_new_frag, False, "FragmentDB")
OV.registerFunction(fdb.update_fragment, False, "FragmentDB")
OV.registerFunction(fdb.delete_fragment, False, "FragmentDB")
OV.registerFunction(fdb.display_large_image, False, "FragmentDB")
OV.registerFunction(fdb.save_picture, False, "FragmentDB")
OV.registerFunction(fdb.store_picture, False, "FragmentDB")
OV.registerFunction(fdb.display_image, False, "FragmentDB")
# OV.registerFunction(fdb.revert_last, False, "FragmentDB")

OV.registerFunction(ref.results, False, "FragmentDB")
OV.registerFunction(fdb.clear_mainvalues, False, "FragmentDB")
