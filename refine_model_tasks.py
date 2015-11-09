'''
Created on 23.10.2015

@author: daniel
'''
#import pprint
from helper_functions import REL_RESTR_CARDS, SHX_CARDS, remove_partsymbol
import os
try:
  import olx  # @UnresolvedImport
  from olexFunctions import OlexFunctions
  from ImageTools import ImageTools
except:
  pass
try:
  OV = OlexFunctions()
  IT = ImageTools()
except:
  pass
#import cProfile
#import pstats
#cp = cProfile.Profile()
#cp.enable(subcalls=True, builtins=True)

class Refmod(object):
    '''
    handles the refinement model of olex2
    '''
    def __init__(self):
      '''
      Constructor
      '''
      self.htm = html_Table()
      
    def fileparser(self, lstfile):
      '''
      gathers the residuals of the lst file
      It searches for the final cycle summary and then for " Disagreeable restraints".
      End is reached with ' Summary of restraints'
      '''  
      if not lstfile:
        return
      disag = False
      final = False
      disargeelist = []
      num = 0
      with open(lstfile, 'r') as f:
        for line in f:
          splitline = line.split()
          if not splitline:
            continue
          if splitline[0].startswith('Observed'):
            continue
          if line.startswith(" Final Structure Factor"):
            final = True
          if final and line.startswith(' Disagreeable restraints'):
            disag = True
            continue
          if line.startswith(' Summary of restraints'):
            final = False
            disag = False
          if disag: 
            # in this case, the desired line is found:
            fline = self.lineformatter(splitline)
            # fline is a list of list
            disargeelist.append(fline)
            num = num+1
            if num > 1000:
              print('Cutting restraints list. Too many restraints...')
              return disargeelist
        return disargeelist
    
    def lineformatter(self, line):
      '''
      takes care of some extra things with different restraints. For example
      RIGU xy should be in one column
      :type line: list
      '''
      for num, i in enumerate(line):
        i = i.replace('/', ' ')
        if i[0].isalpha():
          pos = num
          break
      # remove the part symbol from e.g. F1_2a:
      for num, i in enumerate(line):
        line[num] = remove_partsymbol(i)
      # joining columns without numbers:
      line[pos:] = [' '.join(line[pos:])]
      tline = ' '.join(line)
      for n in REL_RESTR_CARDS:
        if n in tline:
          # adding placeholders for empty fields:
          line = ['-', '-']+line
          break
      return line
    
    def get_listfile(self):
      '''
      returns the path of the current SHELXL list file
      '''
      try:
        lstfile = os.path.abspath(OV.FilePath()+os.path.sep+OV.FileName()+'.lst')
        if not os.path.isfile(lstfile):
          #print('No list file found.')
          return ''
      except:
        print('somethin is wrong with the lst file path.')
        return ''
      return lstfile
    
    def results(self):  
      '''
      prepare the results for the plugin 
      '''
      filedata = self.fileparser(self.get_listfile()) # raw table
      if not filedata:
        filedata = []
      html = self.htm.table_maker(filedata)
      return html

    
class html_Table(object):
  '''
  html table generator
  '''
  def __init__(self):
    try:
      # more than two colors here are too crystmas treelike:
      grade_2_colour = OV.GetParam('gui.skin.diagnostics.colour_grade2')
      self.grade_2_colour = self.rgb2hex(IT.adjust_colour(grade_2_colour, luminosity=1.8)) 
      grade_4_colour = OV.GetParam('gui.skin.diagnostics.colour_grade4')
      self.grade_4_colour = self.rgb2hex(IT.adjust_colour(grade_4_colour, luminosity=1.8))
    except(ImportError, NameError):
      self.grade_2_colour = '#FFD100'
      self.grade_4_colour = '#FF1030'

  def rgb2hex(self, rgb):
    """
    return the hexadecimal string representation of an rgb colour
    """
    return '#%02x%02x%02x' % rgb
  
  def table_maker(self, tabledata=[]):
    '''
    builds a html table out of a datalist from the final 
    cycle summary of a shelxl list file.
    '''
    table=[]
    for line in tabledata:
      table.append(self.row(line))
    if not table:
      return ''
    header = r"""
        <table> 
        <tr> 
          <td align='left'> 
          <b>List of most disagreeable restraints:</b>
              &nbsp;
          </td>
          <td align='right'>
            $spy.MakeHoverButton('small-Short@bitmap', delins more>>addins more -1>>refine 4)
                &nbsp; 
            $spy.MakeHoverButton('small-Full@bitmap', delins more>>addins more -4>>refine 4)            
          </td>
        </tr>
        </table>"""
    footer = ""
    html = r"""
      {0}
    <table border="0" cellpadding="0" cellspacing="6" width="100%" > 
      <tr>
         <td align='center'> Observed </td>
         <td align='center'> Target   </td>
         <td align='center'> Error    </td>
         <td align='center'> Sigma    </td>
         <td align='left'> Restraint  </td> 
      </tr>

      {1}
    </table>
      {2}
      """.format(header, '\n'.join(table), footer)
    return html  
    

  def row(self, rowdata):
    '''
    creates a table row for the restraints list.
    :type rowdata: list
    '''
    td = []
    bgcolor = ''
    try:
      if abs(float(rowdata[2])) > 2.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(self.grade_2_colour)
      if abs(float(rowdata[2])) > 3.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(self.grade_4_colour)        
    except():
      pass
    for num, item in enumerate(rowdata): 
      try:
        # align right for numbers:
        float(item)
        if num < 2: 
          # do not colorize the first two columns:
          td.append(r"""<td align='right'> {} </td>""".format(item))
        else:
          td.append(r"""<td align='right' {0}> {1} </td>""".format(bgcolor, item))
      except:
        if item.startswith('-'):  
          # only a minus sign
          td.append(r"""<td align='center'> {} </td>""".format(item))
        else:
          if num < 4:
            td.append(r"""<td align='right'> {} </td>""".format(item))
            continue
          # align left for words:
          td.append(r"""
            <td align='left'>  
              {}
              <a href="spy.Refmod.edit_restraints({})">  edit </a> 
            </td>""".format(item, item))
    if not td:
      row = "<tr> No (disagreeable) restraints found in .lst file. </tr>"
    else:
      row = "<tr> {} </tr>".format(''.join(td))
    return row

  def edit_restraints(self, restr):
    '''
    this method gets the atom list of an disagreeable restraint in olex2.
    The text string is then formated so that "editatom" can handle it.
    editatom opens an editor window with the respective atoms. 
    :type restr: string like "SAME/SADI Co N11 Co N22"
    '''
    # separate SAME/SADI:
    restr = restr.replace('/', ' ')
    restrlist = restr.split()
    atoms = []
    for i in restrlist:
      if i in SHX_CARDS:
        continue
      if i in ['xz', 'yz', 'xy', 'etc.']:
        continue
      else:
        atoms.append(remove_partsymbol(i))
    OV.cmd('editatom {}'.format(' '.join(atoms)))
    
  
  
  
htm = html_Table()
try:
  OV.registerFunction(htm.edit_restraints, False, "Refmod")
except:
  pass


if __name__ == '__main__':
  import cProfile
  import pstats
  cp = cProfile.Profile()
  cp.enable(subcalls=True, builtins=True)
  ref = Refmod()
  try:
    lst = ref.fileparser(r'D:\Programme\DSR\example\p21c.lst')
    lst = ref.fileparser(r'D:\tmp\big_strukt\p-1.lst')
  except:
    lst = ref.fileparser('/Users/daniel/Documents/DSR/example/p21c.lst')
  tab = htm.table_maker(lst)
  print(tab)
  cp.disable()
  
  pstats.Stats(cp).strip_dirs().sort_stats('cumtime').print_stats(30)
  