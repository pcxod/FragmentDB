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
  #from ImageTools import ImageTools
except:
  pass
try:
  OV = OlexFunctions()
  #IT = ImageTools()
except:
  pass



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
      final = False
      disag = False
      disargeelist = []
      with open(lstfile, 'r') as f:
        for line in f:
          splitline = line.split()
          if not splitline:
            continue
          if splitline[0].startswith('Observed'):
            continue
          if line.startswith(" Final Structure Factor"):
            final = True
#          if line.startswith(' Total number of'):
#            parameters = line.split()[6]
#          if final and line.startswith(' wR2 ='):
#            data = line.split()[7]
            #disargeelist.append(['Results', 'data:', data, 'parameters:', parameters])
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
      line2 = []
      # remove the part symbol from e.g. F1_2a:
      for i in line:
        line2.append(remove_partsymbol(i))
      # joining columns without numbers:
      line[pos:] = [' '.join(line2[pos:])]
      tline = ' '.join(line)
      for n in REL_RESTR_CARDS:
        if n in tline:
          # adding placeholders for empty fields:
          line = ['-', '-']+line
          break
      return line
        
    def open_listfile(self, title = "Select a .lst file", 
            ffilter = '*.lst; *.LST', location = '', default_name = '' ):
      lstfile = olx.FileOpen(title, ffilter, location, default_name)
      return lstfile 
    
    def get_listfile(self):
      '''
      returns the path of the current SHELXL list file
      '''
      try:
        lstfile = os.path.abspath(OV.FilePath()+os.path.sep+OV.FileName()+'.lst')
        if os.path.isfile(lstfile):
          pass
        else:
          print('No list file found.')
          return ''
      except:
        print('somethin is wrong with the lst file path.')
        return ''
      return lstfile
    
    def results(self):  
      '''
      prepare the results for the plugin 
      '''
      filedata = self.fileparser(self.get_listfile())
      if not filedata:
        filedata =[]
      html = self.htm.table_maker(filedata)
      return html

    
class html_Table(object):
  '''
  html table generator
  '''
  def __init__(self):
    pass

  def table_maker(self, tabledata=[]):
    '''
    builds a html table out of a datalist from the final 
    cycle summary of a shelxl list file.
    '''
    table=[]
    for line in tabledata:
      table.append(self.row(line))
    if not table:
      return 'The restraints seem to be well appied.'
    header = r"""
        <table> 
        <tr> 
          <td> 
          <b>List of most disagreeable restraints:</b>
              &nbsp;
          </td>
          <td>
            $spy.MakeHoverButton('small-Short_List@bitmap', delins more>>addins more -1>>refine)
                &nbsp; 
            $spy.MakeHoverButton('small-Full_List@bitmap', delins more>>addins more -4>>refine)            
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
    

  def row(self, rowdata, yellow='#FFD100', red='#FF1030'):
    '''
    creates a table row
    :type rowdata: list
    '''
    td = []
    bgcolor = ''
    try:
      if abs(float(rowdata[2])) > 2.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(yellow)
      if abs(float(rowdata[2])) > 3.5*float(rowdata[3]):
        bgcolor = r"""bgcolor='{}'""".format(red)
    except:
      pass
    for num, item in enumerate(rowdata): 
      try:
        # align right for numbers:
        float(item)
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
          td.append(r"""<td align='left'> 
              <a href="spy.Refmod.edit_restraints({})" style="text-decoration:none"> 
              {} </a> </td>""".format(item, item))
    if not td:
      row = "<tr> No (disagreeable) restraints found in .lst file. </tr>"
    else:
      row = "<tr> {} </tr>".format(''.join(td))
    return row

  def edit_restraints(self, restr):
    '''
    this method gets the atom list of an disagreeable restraint and
    should be able to edit the restraint. 
    '''
    restr = restr.replace('/', ' ')
    restrlist = restr.split()
    atoms = []
    for i in restrlist:
      if i in SHX_CARDS:
        continue
      if i in ['xz', 'yz', 'xy', 'etc.']:
        continue
      else:
        #print(remove_partsymbol(i))
        atoms.append(remove_partsymbol(i))
    OV.cmd('editatom {}'.format(' '.join(atoms)))
    
  
  
  
htm = html_Table()
try:
  OV.registerFunction(htm.edit_restraints, False, "Refmod")
except:
  pass


if __name__ == '__main__':
  ref = Refmod()
  #lst = ref.fileparser('d:\Programme\DSR\example\p21c-test.lst')
  lst = ref.fileparser('/Users/daniel/Documents/DSR/example/p21c.lst')
  print(htm.table_maker(lst))