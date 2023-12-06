from __future__ import print_function

"""
     CCP4SelectionTree.py: 
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2011 University of York
     Copyright (C) 2012 STFC

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
     March 2012 Liz Potterton Copied from CCP4mg model_selection.py
"""

import re,string,types
from core.CCP4ErrorHandling import *
import ccp4mg
import mmdb2 as mmdb
import mmut

#-------------------------------------------------------------------
#-------------------------------------------------------------------
class CMolData:
#-------------------------------------------------------------------
#-------------------------------------------------------------------

  ERROR_CODES = { 101 : { 'description' : 'Atom selection failed. Failed to import MMDB' },
                  102 : { 'description' : 'Atom selection failed. Failed creating mmdb2.Manager object' },
                  103 : { 'description' : 'Atom selection failed. Faied reading coordinate file.' },
                  104 : { 'description' : 'Atom selection failed. Failed parsing command' },
                  105 : { 'description' : 'Atom selection failed. Error creating PPCAtom' },
                  106 : { 'description' : 'Atom selection failed. Error in GetSelIndex' },
                  107 : { 'description' : 'Atom selection failed. Error loading selection tree' },
                  108 : { 'description' : 'Atom selection failed. Error applying selection tree' }
                }

  def __init__(self,fileName):
      try:
        self.molHnd = mmdb.Manager()
        self.molHnd.SetFlag(mmdb.MMDBF_AutoSerials)
        self.molHnd.SetFlag(mmdb.MMDBF_IgnoreRemarks)
        self.molHnd.SetFlag(mmdb.MMDBF_IgnoreBlankLines)
        self.molHnd.SetFlag(mmdb.MMDBF_IgnoreHash)
      
      except:
        raise CException(self.__class__,102)
      else:
        try:
          RC = self.molHnd.ReadCoorFile(fileName)
        except:
          raise CException(self.__class__,103,str(fileName))
      

  def interpretSelection(self,command,fileName=None):
    try:
      toks = SelectionParser().tokenise(command)
    except CException as e:
      raise e
    except:
      raise CException(self.__class__,104,str(command))
    print('interpretSelection',toks)
    try:
      tree = Cselect_tree('TOP',mol=self)
      tree.import_command_string(toks[0],toks[1],toks[2])
    except:
      raise CException(self.__class__,107,str(command))
    try:
      status,selHnd = tree.apply()
    except:
      raise CException(self.__class__,108,str(command))
    print('interpetSelection',status,selHnd)
    try:
      selAtoms = mmdb.newPPCAtom()
    except:
      raise CException(self.__class__,105,str(command))
    try:
      try:
        selindexp = mmut.intp()
        selAtoms = mmut.GetAtomSelIndex(self.molHnd,selHnd,selindexp)
        nselatoms = selindexp.value()
      except:
        nselatoms = self.molHnd.GetSelIndex(selHnd,selAtoms)
    except:
      raise CException(self.__class__,106,str(command))
    #print 'interpetSelection nselatoms',nselatoms

    if fileName is not None:
      self.writeSelection(selHnd,fileName)

  def writeSelection(self,selHnd, fileName, format='PDB' ):
    # Does an atom-by-atom copy so potentially slow
    mmdb2 = mmdb.Manager()
    mmdb2.Copy (self.molHnd ,mmdb.MMDBFCM_Title | mmdb.MMDBFCM_Cryst  )
  
    if self.molHnd.GetNumberOfModels()>1 and False:
      #This method does not exist!
      self.molHnd.CopySelection(selHnd,mmdb2)
    else:
      try:
        selindexp = mmut.intp()
        selAtoms = mmut.GetAtomSelIndex(self.molHnd,selHnd,selindexp)
        nSelAtoms = selindexp.value()
        for i in range(0,nSelAtoms):
          pcat = mmdb.getPCAtom(selAtoms,i)
          if self.molHnd.GetNumberOfModels() == 1 or pcat.GetModelNum() == 1:
              mmdb2.PutAtom ( i+1,pcat ) 
      except:
        selAtoms = mmdb.newPPCAtom()
        nSelAtoms = self.molHnd.GetSelIndex(selHnd,selAtoms)
        for i in range(0,nSelAtoms):
          pcat = mmdb.CAtomPtr(mmdb.getPCAtom(selAtoms,i))
          mmdb2.PutAtom ( i+1,pcat ) 

    if format == "PDB":
      RC = mmdb2.WritePDBASCII( fileName )
    else:
      RC = mmdb2.WriteCIFASCII( fileName )
    print(format,'file, containing',nSelAtoms,'atoms, written to',fileName)
    del mmdb2
    return RC

    
#-------------------------------------------------------------------
#-------------------------------------------------------------------
class SelectionParser:
#-------------------------------------------------------------------
#-------------------------------------------------------------------

  ERROR_CODES = { 101: { 'description' : 'Error in selection command' }
                }

#-------------------------------------------------------------------
  def command_to_list( self,command ):
#-------------------------------------------------------------------
    '''
    Convert an input command string to a python list selection
    '''
    toks = self.tokenise(command)
    #toks=self.expand(toks[1],toks[2],toks[3])
    t = Cselect_tree('TOP','',self)
    rv = t.import_command_string(toks[0],toks[1],toks[2])
    l = t.get()
    t.cleanup()
    return l
    
#-------------------------------------------------------------------
  def list_to_command(self,list):
#-------------------------------------------------------------------
    t = Cselect_tree('TOP','',self)
    t.set(list)
    rv = t.export_command_string()
    t.cleanup()
    return rv

#-------------------------------------------------------------------------
  def tokenise(self,inp_str,report=1):
#-------------------------------------------------------------------------
    stdops = ['and','or','xor','excl','not']
    elements = []
    ops = []
    levels = []
    wk_str = re.sub('{',' { ',inp_str)
    wk_str = re.sub('}',' } ',wk_str)
    wk_str = re.sub('not',' not ',wk_str)
    # This is a kludge to hide operaters within quotes
    quoted = re.findall('".*?"',wk_str)
    #print "quoted",quoted
    if len(quoted)>0:wk_str = re.sub('".*?"','QuOTE',wk_str)
      
    rx = re.compile('\{| and | or | xor | excl | not |\}')
    nq = -1
    for ele in rx.split(wk_str):
      qu = re.findall('QuOTE',ele)
      
      for ii in range(0,len(qu)):
        nq = nq+1
        ele = re.sub('QuOTE',quoted[nq],ele,1)
      elements.append(ele.strip())
    for op in rx.findall(wk_str): 
      ops.append(op.strip())
    levels = self.getlevels(ops)
    

    #print "elements",elements
    #print "ops",ops
    #print "levels",levels

    if len(ops) > 0:

# Test the first operator
      if ops[0] == '{':
        if elements[0] != '':
          #return [3,'Missing operator before open brace at ' \
          #      + "'" + elements[0] + ' ' + ops[0] + "'"]
          raise CException(self.__class__,101,
            'Missing operator before open brace at '+ "'" + elements[0] + ' ' + ops[0] + "'")
      else:
        if elements[0] == '' and ops[0] != 'not':
          #return [4,'Missing element before ' + ops[0]]
          raise CException(self.__class__,101,'Missing element before ' + ops[0])

      lop = stdops.count(ops[0]) 
      for i in range(1,len(ops)):
        cop = stdops.count(ops[i])
        dlev = levels[i] - levels[i-1]
        if lop and cop:
# test if two consecutive operators
          if elements[i] == '' and ops[i] != 'not':
            #return [4,'Missing element between ' + ops[i-1] + ' ' + ops[i]]
            raise CException(self.__class__,101,'Missing element between ' + ops[i-1] + ' ' + ops[i])
        if lop:
# test if there is operator before/after braced group
          if dlev > 0:
            if elements[i] != '':
              #return [3,'Missing operator before open braces at ' \
              #    + "'" + ops[i-1] + ' ' + elements[i] + ' ' + ops[i] + "'"]
              raise CException(self.__class__,101,
                 'Missing operator before open braces at '+ "'" + ops[i-1] + ' ' + elements[i] + ' ' + ops[i] + "'")
          elif ops[i] != 'not':
            if elements[i] == '':
              #return [4,'Missing element between ' + ops[i-1] + ' ' + ops[i]]
              raise CException(self.__class__,101,
                        'Missing element between ' + ops[i-1] + ' ' + ops[i] )
        if cop:
          if dlev < 0:
            if elements[i] != '':
              #return [3,'Missing operator after close braces at ' \
              #  + "'" + ops[i-1] + ' ' + elements[i] + ' ' + ops[i] + "'"]
              raise CException(self.__class__,101,
                 'Missing operator after close braces at ' + "'" + ops[i-1] + ' ' + elements[i] + ' ' + ops[i] + "'" )
          elif ops[i] != 'not':
            if elements[i] == '':
              #return [4,'Missing element between ' + ops[i-1] + ' ' + ops[i]]
              raise CException(self.__class__,101,
                 'Missing element between ' + ops[i-1] + ' ' + ops[i] )
        else:
# test for missing operator between braced groups
          if ops[i-1] == '}' and ops[i] == '{':
            #return [5,'Missing operator between braces } { ']
            raise CException(self.__class__,101, 'Missing operator between braces } { ' )
        lop = cop
      if lop:
        if elements[-1] == '':
          #return [4,'Missing element after final ' + ops[-1]]
          raise CException(self.__class__,101, 'Missing element after final ' + ops[-1] )
      else:
        if elements[-1] != '':
          #return [3,'Missing operator after braces at ' \
          #      + "'" + ops[-1] + ' ' + elements[-1] + "'"]
          raise CException(self.__class__,101,
               'Missing operator after braces at '+ "'" + ops[-1] + ' ' + elements[-1] + "'")
    #print 'tokenise',[elements,ops,levels]
    return [elements,ops,levels]

#-------------------------------------------------------------------------
  def expand(self,elements,ops,levels):
#-------------------------------------------------------------------------
    # Expand a command string to replace any aliases with their
    # full command
    initl = len(elements)
    i = 0
    while i < len(elements):
      toks = self.expand0(elements[i])
      if len(toks)==0:
        pass
      else:
        elements.pop(i)
# There is an expansion which needs to be inserted into the lists
        if len(toks[0]) == 1:
          elements.insert(i,toks[0][0])
        else:
          elements.insert(i,'')
          elements0 =  toks[0]
          elements0.reverse()
          for e in elements0: elements.insert(i,e)
          elements.insert(i,'')
          ops0 = toks[1]
          ops0.reverse()
          ops.insert(i,'}')
          for o in ops0: ops.insert(i,o)
          ops.insert(i,'{')
#      print 'expanded elements',elements,'ops',ops
   
    
    if len(elements) != initl:
      levels = self.getlevels(ops)
    return [elements,ops,levels]

#-------------------------------------------------------------------------
  def expand0(self,element):
#-------------------------------------------------------------------------
# Test if first word of element corresponds to a defined protocol
# or a user selection alias
    if element == '': return []
    words = splitWords(element)
    fword = words[0]
    expanded = self.expand_selection_alias(fword)
    if expanded:
      return self.tokenise(expanded)
    else:
      return []

#-------------------------------------------------------------------------
  def getlevels(self,ops):
#-------------------------------------------------------------------------
    level = 0
    levels =[]
    for op in ops:
      if op == '{':
        level = level +1
        levels.append(level)
      elif op == '}':
        levels.append(level)
        level = level - 1
        if level < 0:
          #return (2,'Unmatched closing brace',len(levels))
          raise CException(self.__class__,101,'Unmatched closing brace',len(levels))
      else:
        levels.append(level)
    if level > 0:
      #return (1,'Unmatched opening brace',0)
      raise CException(self.__class__,101,'Unmatched opening brace')
    else:
      return levels


    
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
class Cselect_tree:
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

  aminoResidues = ["GLY","ALA","VAL","PRO","SER","THR","LEU","ILE","CYS","ASP","GLU","ASN","GLN","ARG","LYS","MET","MSE","HIS","PHE","TYR","TRP","HCS","ALO","PDD","UNK"]
  metalElements = ["LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "UN", "UU", "UB", "UQ", "UH", "UO"];
  saccharideResidues = ["BMA","MAN","NAG","GLC","BGC","GCS","GAL","NGA","MGC","NG1","NG6","A2G","6MN","GLP","GP4","BEM","KDN","XLS","CXY","RBY","TDX","XYL","XYS","XYP","FCA","FCB","FUC","GTR","ADA","DGU","SI3","NCC","NGF","NGE","NGC","GC4","GCD","GCU","GCV","GCW","IDS","REL","IDR","SIA"];
  nucleicResidues = [ "DC","DT","U","T","C","PSU","DA","DG","A","G","I" ]
  nucleotideResidues = [ "Ad","Cd","Gd","Td","ADE","CYT","GUA","INO","THY","URA","AMP","ADP","ATP","CMP","CDP","CTP","GMP","GDP","GTP","TMP","TDP","TTP" ]
  waterResidues = [ "HOH","H2O","WAT","SOL","DOD","D2O" ]
  metalResidues = [ "LI","BE","NA","MG","AL","K","CA","SC","TI","V","MN","FE","CO","NI","CU","ZN","GA","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR","RF","DB","SG","BH","HS","MT","UN","UU","UB","UQ","UH","UO" ]
  soluteResidues = [ "SUL","SO4","NO3","MOH","EOH","GOL","ACT","CL","BR","PG4","PG5","PG6","1PE","2PE","7PE","PE3","PE4","PE5","PE6","PE7","PGE","XPE","C10","CE1","CE9","CXE","EDO","N8E","P33","P4C","12P","15P","DMS","IOD","MRD","BE7","MHA","BU3","PGO","BU2","PDO","BU1","1BO","TFP","DHD","PEU","TRS","TAU","SBT","SAL","MPD","IOH","IPA","BU2","PIG","B3P","BTB","B3P","NHE","C8E","OTE","PE8","2OS","1PS","CPS","DMX","MPO","DXG","CM5","ACA","ACN","CCN","DR6","NH4","AZI","LAK","BCN","BRO","CAC","CBX","FMT","ACY","CBM","CLO","FCL","CIT","3CO","NCO","CU1","CYN","CYN","MA4","BTC","TAR","MTL","DPR","SOR","SYL","DDQ","DMF","DIO","DOX","SDS","EEE","EGL","FLO","TRT","FCY","FRU","GBL","GPX","HTO","HTG","B7G","16D","HEZ","IDO","ICI","ICT","TLA","LDA","MRY","BEQ","C15","MG8","POL","JEF","DIA","IPH","PIN","CRY","PGR","PGQ","SPD","SPK","SPM","TBU","TMA","TEP","SCN","ETF","144","UMQ","URE","CN","EPE","HEPES","HEPE","MES","PEG","ETA","HED","IMD","BO3","IUM","PO4" ]
#-------------------------------------------------------------------------
  def __init__(self,name,parent=None,mol=None,dispobj=''):
#-------------------------------------------------------------------------
    self.name = name
    self.parent = parent
    self.mol = mol
    self.dispobj = dispobj
    self.element = None
    self.op = None
    self.children = []
    self.selHnd = None
    self.com_string = ""
  

#-------------------------------------------------------------------------
  def set(self,defn):
#-------------------------------------------------------------------------
    #print "set defn",defn
    if len(defn) < 1:
      return 1
    self.op = defn[0]
    if not self.op :
      self.element = defn[1]
    elif self.op == 'not':
      cl = Cselect_tree(self.name + '_l',self.name,self.mol,dispobj=self.dispobj)
      cl.set(defn[1])
      self.children.append(cl)      
    elif len(defn)>2:
      cl = Cselect_tree(self.name + '_l',self.name,self.mol,dispobj=self.dispobj)
      cl.set(defn[1])
      self.children.append(cl)
      cr = Cselect_tree(self.name + '_r',self.name,self.mol,dispobj=self.dispobj)
      cr.set( defn[2])
      self.children.append(cr)
    else:
      return 1
    #print "RETURN 0"
    return 0
#-------------------------------------------------------------------------
  def import_command_string(self,elements,ops,levels):
#-------------------------------------------------------------------------

    #print 'elements',elements,'ops',ops,'levels',levels

# check if down to one component
    if len(elements) == 1:
      self.element = elements[0]
#### USEFUL DIAGNOSTIC ON TREE STRUCTURE
#      print self.name, self.element
      return 0
# if enclosed in braces then strip them off
    while levels.count(0) == 0:
      l = len(levels)
      levs = []
      for lev in levels[1:l-1]: levs.append(lev -1)
      levels = levs
      ops = ops[1:l-1]
      l = len(elements)
      elements = elements[1:l-1]
# check again if we are down to 1 component
      if len(ops) <= 0:
        self.element = elements[0]
#        print self.name, self.element
        return 0
#      print 'stripped levels',levels,'elements',elements,'ops',ops

# find the *last* unnested logical operator
    iop = len(levels)-1
    while (levels[iop] != 0 or ops[iop] == 'not') and iop > 0: 
      iop = iop -1
    self.op = ops[iop]
#    print self.name, self.op
    if self.op != 'not':
      cl = Cselect_tree(self.name + '_l',self.name,self.mol,dispobj=self.dispobj)
      cl.import_command_string(elements[0:(iop+1)],ops[0:iop],levels[0:iop])
      self.children.append(cl)
    cr = Cselect_tree(self.name + '_r',self.name,self.mol,dispobj=self.dispobj)
    cr.import_command_string( elements[iop+1:], ops[iop+1:], levels[iop+1:])
    self.children.append(cr)
    return 0
  
#-------------------------------------------------------------------------
  def cleanup(self,save_selHnd=0):
#-------------------------------------------------------------------------
# Close the selHnd for all but the top object
    if self.selHnd and not save_selHnd: self.delete_selhnd (self.selHnd)
    for child in self.children:
      if not save_selHnd or child.selHnd != self.selHnd:
        child.cleanup()

#-------------------------------------------------------------------------
  def __del__(self):
#-------------------------------------------------------------------------
    for child in self.children:
      del child

#-------------------------------------------------------------------------
  def skey(self):
#-------------------------------------------------------------------------
    if self.op == 'and':
      return mmdb.SKEY_AND
    elif self.op == 'or':
      return mmdb.SKEY_OR
    elif self.op == 'xor':
      return mmdb.SKEY_XOR
    elif self.op == 'excl':
      return mmdb.SKEY_CLR
    

#-------------------------------------------------------------------------
  def apply(self):
#-------------------------------------------------------------------------

# propagate the 'apply' down through all nodes in the tree

    for child in self.children:
      if not child.selHnd and child.op:
        child.apply()

# For an 'ops' object apply the selection criteria of its 
# children.  Create a 'selHnd' for this object
#    print "in apply", self.name

    if len(self.children) == 0:
# There is one simple selection command
      self.selHnd = self.new_selHnd()
      rv = self.interpret(mmdb.SKEY_NEW,self.element)
      return rv
 
    elif self.op == 'not':
      if not self.children[0].selHnd:
        self.selHnd = self.new_selHnd()
        rv = self.interpret(mmdb.SKEY_NEW,self.children[0].element)
        if rv[0] > 0: return rv
      else:
        self.selHnd = self.children[0].selHnd
      rv = self.interpret(mmdb.SKEY_XOR,'/0')
      return rv

    if not self.children[0].selHnd and \
       not self.children[1].selHnd:
# Neither child selection is currently evaluated
      self.selHnd = self.new_selHnd()
      rv = self.interpret(mmdb.SKEY_NEW,self.children[0].element)
      if rv[0] > 0: return rv
      rv = self.interpret(self.skey(),self.children[1].element)
      if rv[0] > 0: return rv

    elif self.children[0].selHnd and self.children[1].selHnd:
# both child selections are evaluated
      rv = self.merge_selection(self.skey(), \
        self.children[0].selHnd,self.children[1].selHnd)
      self.delete_selhnd(self.children[1].selHnd)
      if rv[0] > 0: return rv
      self.selHnd = self.children[0].selHnd

    elif self.children[0].selHnd:
# the left side selection is already evaluated
      self.selHnd = self.children[0].selHnd
      rv = self.interpret(self.skey(),self.children[1].element)
      if rv[0] > 0: return rv
    else:
# the right side selection is already evaluated
      if self.op != 'excl':
        self.selHnd = self.children[1].selHnd
        rv = self.interpret(self.skey(),self.children[0].element)
        if rv[0] > 0: return rv
      else:
# excl operation not cummutative
        self.selHnd = self.new_selHnd()
        rv = self.interpret(mmdb.SKEY_NEW,self.children[0].element)
        if rv[0] > 0: return rv
        self.merge_selection(self.skey(),self.selHnd,self.children[1].selHnd)
        self.delete_selhnd(self.children[1].selHnd)
        if rv[0] > 0: return rv
    return [0,self.selHnd]


#-------------------------------------------------------------------------
  def new_selHnd(self):
#-------------------------------------------------------------------------
    selHnd = self.mol.molHnd.NewSelection()
    return selHnd

  '''
#-------------------------------------------------------------------------
  def interpret(self,skey,element):
#-------------------------------------------------------------------------
#   

    #print "interpret element",element
    words = splitWords0(element)
    #print "interpret words",words
    if len(words) <= 0:  return [1,'No selection']
    # Test if first word is an alias and
    # parse the expansion for the alias
    expansion = self.mol.expand_selection_alias(words[0])
    #print "expansion",words[0],expansion
    if expansion:
      if skey == SKEY_NEW: self.delete_selhnd(self.selHnd)
      rv = self.mol.parse_selection(expansion)
      if rv[0]:
        return [1,'Error interpreting expansion of '+words[0]]
      elif skey == SKEY_NEW:
        self.selHnd = rv[1]
      else:
        self.merge_selection(skey,self.selHnd,rv[1])
        self.delete_selhnd(rv[1])
      return (0,self.selHnd)
        
    kw = matchKeyword(string.lower(words[0]),self.mol.selection_commands)
    #print "word,kw",words[0],kw
    if kw[0] >= 1: 
      cmd = 'rv  = self.mol.' + self.mol.selection_methods[kw[1][0]] + '('
      for w in words[1:]:
        ww = re.split('=', w)
        if len(ww) > 1:
          if ( re.match("\'", ww[1]) and re.search("'$",ww[1]) ) or \
            ( re.match('\"', ww[1]) and re.search('"$',ww[1]) ):
            cmd = cmd + ww[0] + "=" + ww[1] + ","
          else:
            cmd = cmd + ww[0] + "='" + ww[1] + "',"
        else:
          cmd = cmd + "'" + w + "',"
      cmd = cmd + 'selhnd=' +  str(self.selHnd) + ',' + 'skey=' + str(skey)  + ')'
      #print "interpret cmd",cmd,self.nSelAtoms()
      if DEVELOPER():
        exec cmd
        if rv[0] == 0:
          self.selHnd = rv[1]
          return rv
        else:
          return [rv[0],'Error interpreting command ' + self.mol.selection_commands[kw[1][0]] + ': ' + rv[1]]
      else:
        try:
          exec cmd
          if rv[0] == 0:
            self.selHnd = rv[1]
            return rv
          else:
            return [rv[0],'Error interpreting command ' + self.mol.selection_commands[kw[1][0]] + ': ' + rv[1]]
        except:
           #import sys
           #print 'ModelAnalysis.interpret',sys.exc_info()
           return [1,'Error interpreting command ' + self.mol.selection_commands[kw[1][0]]]

    # Does the element contain operators - if so attempt to
    # interpret
    elif re.search(r'[=><]', element):
      #print "matched operator"
      sel_ops = []
      ele = element
      while ele:
        t = re.search(r'([^=><]*)(==|!=|<=|>=|<|>)([^=><]*)(.*)',ele)
        #if not t:
        #  t = re.search(r'([^=><]*)(=|>|<)([^=><]*)(.*)',ele)
        if not t:
          # Have an operator symbol but it dont make sense
          return[1,'Can not interpret '+element]
        
        tt = t.groups()
        #print "tt",tt
        # tt[1] should be operator and tt[0] and tt[2] are parameter
        # and limit value (either way round)
        tt0 = tt[0].strip()
        tt2 = tt[2].strip()
        param = self.interpret_operator_param(tt0)
        if param:
          op = ['==','!=','<=','<','>=','>','='].index(tt[1])
          if op == 6: op = 0
          value = tt2
        else:
          param = self.interpret_operator_param(tt2)
          if param:
            op = ['==','!=','>=','>','<=','<','='].index(tt[1])
            if op == 6: op = 0
            value = tt0
          else:
            return[1,'No parameter recognised in '+element]
        try:
          if param[1]=='float':
            value=float(value)
          elif param[1]=='int':
            # Weird this - something like int('6.0') breaks
            value = int(float(value))
          else:
            pass
        except:
           return[1,'Inappropriate parameter value '+value]
        sel_ops.append([param[0],param[1],op,value])
        #print "sel_ops",sel_ops
        if tt[3]:
          ele = param[0]+tt[3]
        else:
          ele = ''

      # Success! put return fail for now
      rv_op = self.apply_operators(sel_ops)
      #print "from apply_operators rv",rv_op
      if rv_op[0]==0:
        # merge the selection from the operator and
        # delete the selHnd
        rv = self.merge_selection(skey,self.selHnd,rv_op[1])
        #print "merge_selection rv",rv,self.selHnd,self.nSelAtoms(rv[1])
        self.delete_selhnd(rv_op[1])
        if (rv[0] == 0 ):
          return [0,self.selHnd]
        else:
          return [1,'Error merging '+element]
      else:        
        return rv_op
        
    else:
      #element0 = self.interpretSymop(element)
      rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,element,skey)
      #print 'interpret selection',element,self.nSelAtoms()
      if rv == 0:
        return [rv,self.selHnd]
      else:
        return [rv,'Error applying selection CID ' + element]
#------------------------------------------------------------------------
  def interpretSymop(self,element):
#------------------------------------------------------------------------
    m = re.match(r"(.*)(\d+)_(\d+)(.*)",element)
    if not m: return element
    gps = m.groups()
    #print 'interpretSymop',element,gps
    op = int(gps[1])
    tran = int(gps[2])
    
    if not ['','/'].count(gps[0]) or \
       op<1 or op>self.mol.model_symmetry.nofCrystalSymops() or \
       tran<111 or tran>999 :
      return element
    
    imodel = self.mol.model_symmetry.getPendingSymmetryMate(gps[1]+'_'+gps[2])
    if imodel<0:
      return element
    else:
      return '/'+str(imodel)+gps[3]
    
#------------------------------------------------------------------------
  def interpret_operator_param( self, tt):
#------------------------------------------------------------------------
    #print "interpret_operator_param",tt
    tt = string.upper(tt)
    for par in ['B','OCC','CHARGE','X','Y','Z','ATOM_SAS','RES_SAS', \
                'ATOM_BURIED','RES_BURIED']:
      if re.match(par,tt): return [par,'float']    
    for par in ['SERIAL']:
      if re.match(par,tt): return [par,'int']
    for par in ['SEC']:
      if re.match(par,tt): return [par,'string']
    return []
  '''
      
#-------------------------------------------------------------------------
  def interpret(self,skey,element):
#-------------------------------------------------------------------------
      print('Cselect_tree.interpret',self.selHnd,mmdb.STYPE_ATOM,element,skey)

      if element == "peptide":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.aminoResidues))+")/*:*",skey)
      elif element == "nucleic":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.nucleicResidues))+")/*:*",skey)
      elif element == "nucleotide":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.nucleotideResidues))+")/*:*",skey)
      elif element == "ligands":
        tempSelHnd = self.mol.molHnd.NewSelection()
        ligandSelHnd = self.mol.molHnd.NewSelection()
        rv = self.mol.molHnd.Select(ligandSelHnd,mmdb.STYPE_ATOM,"/*/*/*/*:*",mmdb.SKEY_NEW)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.aminoResidues))+")/*:*",mmdb.SKEY_NEW)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.nucleicResidues))+")/*:*",mmdb.SKEY_OR)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.nucleotideResidues))+")/*:*",mmdb.SKEY_OR)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.waterResidues))+")/*:*",mmdb.SKEY_OR)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.soluteResidues))+")/*:*",mmdb.SKEY_OR)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.saccharideResidues))+")/*:*",mmdb.SKEY_OR)
        rv = self.mol.molHnd.Select(tempSelHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.metalResidues))+")/*["+(",".join(self.metalElements))+"]:*",mmdb.SKEY_OR)
        self.mol.molHnd.Select(ligandSelHnd,mmdb.STYPE_ATOM,tempSelHnd,mmdb.SKEY_XOR)
        self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,ligandSelHnd,skey)
        self.mol.molHnd.DeleteSelection(tempSelHnd)
        self.mol.molHnd.DeleteSelection(ligandSelHnd)
      elif element == "water":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.waterResidues))+")/*:*",skey)
      elif element == "solute":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.soluteResidues))+")/*:*",skey)
      elif element == "saccharide":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.saccharideResidues))+")/*:*",skey)
      elif element == "metal":
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,"/*/*/("+(",".join(self.metalResidues))+")/*["+(",".join(self.metalElements))+"]:*",skey)
      else:
        rv = self.mol.molHnd.Select(self.selHnd,mmdb.STYPE_ATOM,element,skey)
      #print 'interpret selection',element,self.nSelAtoms()
      if rv == 0:
        return [rv,self.selHnd]
      else:
        return [rv,'Error applying selection CID ' + element]


#------------------------------------------------------------------------
  def apply_operators( self, sel_ops ):
#------------------------------------------------------------------------
    import model
    import re
    import string
    import utils
    selHnd=self.mol.molHnd.NewSelection()
    skey = mmdb.SKEY_NEW
    for op in sel_ops:
      #print "apply_operators op",op
      
      if not SelHandle.property_alias.count(op[0]): break
      # Make sure SAS calculated
      
      if ['ATOM_SAS','RES_SAS'].count(op[0]):
        self.mol.update_sasarea(reapply=0)
      if re.search('BURIED',op[0]):
        #n = string.split(op[0],'-')[1]
        #print "contact n",n
        if re.match('ATOM',op[0]):
          udd =  self.mol.molHnd.GetUDDHandle( mmdb.UDR_ATOM,'atom_contact'+self.dispobj)
        else:
          udd =  self.mol.molHnd.GetUDDHandle( mmdb.UDR_RESIDUE,'residue_contact'+self.dispobj)
        #print "apply_operators contact UDD",udd
      else:
        if op[0] == 'SERIAL':
          udd = self.mol.molHnd.GetUDDHandle( mmdb.UDR_ATOM,"atomSerial")
          if udd<=0:
            serial_molHnd = mmdb.Manager()
            serial_molHnd.SetFlag(MMDBF_IgnoreRemarks)
            serial_molHnd.SetFlag(MMDBF_IgnoreBlankLines)
            serial_molHnd.SetFlag(MMDBF_IgnoreHash)
            #serial_molHnd.SetFlag(MMDBF_EnforceUniqueChainID)
            #serial_molHnd.ReadCoorFile(self.mol.filename[2])
            utils.ReadCoorFile(serial_molHnd,self.mol.filename[2])
            udd = self.mol.molHnd.LoadSerial(serial_molHnd)
            del serial_molHnd
        else:
          #Get the enum code for the type of data
          prop_enum = SelHandle.property_enum[ \
             SelHandle.property_alias.index(op[0])]
          #print "in apply_operators prop_enum",prop_enum
          udd = self.mol.molHnd.LoadUDDData(prop_enum )
          #print "udd",udd
      if udd < 0:
        self.mol.molHnd.DeleteSelection(selHnd)
        return [1,"Error loading data "+op[0]]
      
      if op[0] == 'SEC':
        #print "op[3]",op[3]
        # Beware pre 1.111 used single letter code secstr_code
        # but now use the secstr_alias
        sec_code = model.MolData.secstr_code
        sec_alias = model.MolData.secstr_alias
        ss = -1
        if sec_code.count(op[3]):
          ss = sec_code.index(op[3])
        elif sec_alias.count(op[3]):
          ss = sec_alias.index(op[3])
        if ss>=0:
          min = int(ss)
          max = int(ss)
        else:
          min=0
          max=6
      elif op[1] == 'int':
        if op[2] == 0:
          max = int(op[3])
          min =  int(op[3])
        elif op[2] == 1:
          max = int(op[3])-1
          min = -99999999
        elif op[2] == 2:
          min = -99999999
          max =  int(op[3])
        elif op[2] == 3:
          min = -99999999
          max =  int(op[3])-1
        elif op[2] == 4:
          max = 99999999
          min =  int(op[3])
        elif op[2] == 5:
          max = 99999999
          min =  int(op[3])+1
      elif op[1] == 'float':
        if op[2] == 0:
          if op[3]>0.0:
            max = float(op[3])*1.0001
            min =  float(op[3])*0.9999
          else:
            max = float(op[3])*0.9999
            min =  float(op[3])*1.0001
        elif op[2] == 1:
          if op[3]>0.0:
            max = float(op[3])*0.9999
          else:
            max = float(op[3])*1.0001
          min = -9999999.9
        elif op[2] == 2:
          min = -9999999999.9
          max =  float(op[3])
        elif op[2] == 3:
          min = -9999999999.9
          max =  float(op[3])
        elif op[2] == 4:
          max = 9999999999.9
          min =  float(op[3])
        elif op[2] == 5:
          max = 9999999999.9
          min =  float(op[3])
          
      #print "min,max",min,max,op[0],udd
      if re.match('RES',op[0]) or op[0]=='SEC':
        resSelHnd=self.mol.molHnd.NewSelection()
        self.mol.molHnd.SelectUDD (resSelHnd,mmdb.STYPE_RESIDUE,udd, \
            min,max,mmdb.SKEY_NEW)
        self.mol.molHnd.Select(selHnd,mmdb.STYPE_ATOM,resSelHnd,skey)
        self.mol.molHnd.DeleteSelection(resSelHnd)
        resSelHnd = -1
      else:
        self.mol.molHnd.SelectUDD (selHnd,mmdb.STYPE_ATOM,udd, \
            min,max,skey)
      skey = mmdb.SKEY_AND
      if op[2] == 1:
      # Need a second selection to get the atoms above the value
        if op[1] == 'int':
          min = int(op[3])+1
          max = 99999999
        elif op[1] == 'float':
          if op[3]>0:
            min = float(op[3])*1.0001
          else:
            min = float(op[3])*0.9999
          max = 9999999999.9
        self.mol.molHnd.SelectUDD (selHnd,mmdb.STYPE_ATOM,udd, \
            min,max,mmdb.SKEY_OR)
      #print "apply_operator nSelAtoms",self.nSelAtoms(selHnd)
    return [0,selHnd]
  
#-------------------------------------------------------------------------
  def merge_selection(self,skey,selHnd1,selHnd2):
#-------------------------------------------------------------------------
    #print 'CSelectTree.merge_selection',skey,selHnd1,selHnd2
    # Beware Select_propagate does not return a status flag
    #print "SelHandle.merge_selection skey",skey
    self.mol.molHnd.Select(selHnd1,mmdb.STYPE_ATOM,selHnd2,skey)
    return [0,selHnd1]

#-------------------------------------------------------------------------
  def delete_selhnd(self,selHnd):
#-------------------------------------------------------------------------
    if selHnd >= 0: self.mol.molHnd.DeleteSelection(selHnd)
    selHnd = -1

#-------------------------------------------------------------------------
  def nSelAtoms(self,selhnd=None):
#-------------------------------------------------------------------------
    if not selhnd: selhnd = self.selHnd
    #if self.selAtoms: delPPCAtom(self.selAtoms)
    try:
      selindexp = mmut.intp()
      self.selAtoms = mmut.GetAtomSelIndex(self.mol.molHnd,selhnd,selindexp)
      nselatoms = selindexp.value()
    except:
      self.selAtoms = newPPCAtom()
      nselatoms =self.mol.molHnd.GetSelIndex(selhnd,self.selAtoms)
#    print "Number of selected atoms",nselatoms
    #delPPCAtom(self.selAtoms)
    return nselatoms

#-------------------------------------------------------------------------
  def export_command_string(self):
#-------------------------------------------------------------------------
    for child in self.children:
      if not child.com_string:
        child.export_command_string()

    if len(self.children) == 0:
      self.com_string = self.element

    elif self.op == 'not':
      #print "export not",self.children[0].element
      if len(self.children[0].children) > 1:
        self.com_string = 'not {'+self.children[0].com_string+'}'
      else:
        self.com_string = 'not '+ self.children[0].com_string

    else:
      if not self.children[0].op or self.children[0].op == self.op:
        self.com_string =  self.children[0].com_string + ' ' + self.op + ' '
      else:
        #print self.children, self.op
        self.com_string = ' {'+self.children[0].com_string+'} '+self.op+' '

      if not self.children[1].op or self.children[1].op == self.op:
        #print "export_command_string",self.com_string,'*',self.children[1].com_string
        self.com_string = self.com_string+self.children[1].com_string
      else:
        self.com_string = self.com_string+'{'+self.children[1].com_string+'}'

    #print "com_string",self.com_string
    return [0,self.com_string]

#------------------------------------------------------------------------
  def search(self,keyword):
#------------------------------------------------------------------------
    hits = self.search0(keyword)
    # search0 returns a list of Cselect_tree objects
    # convert this to a list of selection definitions
    hitdefns = []
    for hit in hits:
      hitdefns.append(hit[1].get())

    #print "search",keyword,hitdefns
    return hitdefns


#-------------------------------------------------------------------------
  def search0(self,keyword):
#-------------------------------------------------------------------------
    # At a leaf - does the com_string match the keyword
    #print "search0",self.element,self.op
    if len(self.children) == 0:
      if self.element == keyword:
        return [[1,self]]
      else:
        return []

    # Loop over all children
    hits = []
    for child in self.children:
      hits.extend(child.search0(keyword))

    newhits=[]
    if len(hits)>0:
      if self.op == 'and':
        for hit in hits:
          if hit[0] == 2:
            newhits.append(hit)
          elif hit[0] == 1:
            newhits.append([1,self])
      else:
        for hit in hits: newhits.append([2,hit[1]])

    #print "search0",self.element,newhits
    return newhits

#----------------------------------------------------------------------
  def get(self):
#----------------------------------------------------------------------
    if not self.op:
      return [None,self.element]
    elif self.op == 'not':
      return ['not',self.children[0].get()]
    else:
      return [self.op,self.children[0].get(),self.children[1].get()]


#-------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------
#-------------------------------------------------------------------
def list_to_command(list):
#-------------------------------------------------------------------
  t = Cselect_tree('TOP','')
  t.set(list)
  rv = t.export_command_string()
  t.cleanup()
  return rv

#-------------------------------------------------------------------
def command_to_list( command ):
#-------------------------------------------------------------------
  '''
  Convert an input command string to a python list selection
  '''
  parser = SelectionParser()
  toks = parser.tokenise(command)
  if toks[0] > 0:
    return [1,'Error parsing command']
  else:
    #toks=self.expand(toks[1],toks[2],toks[3])
    t = Cselect_tree('TOP','')
    rv = t.import_command_string(toks[1],toks[2],toks[3])
    l = t.get()
    t.cleanup()
    return [0,l]



def interpretSelection(fileName,command,fileOut=None):
  mol = MolData(fileName)
  toks = SelectionParser().tokenise(command)
  print('interpretSelection',toks)
  if toks[0] > 0:
    return [1,'Error parsing command']
  tree = Cselect_tree('TOP',mol=mol)
  tree.import_command_string(toks[1],toks[2],toks[3])
  status,selHnd = tree.apply()
  print('interpetSelection',status,selHnd)  
  try:
    selindexp = mmut.intp()
    selAtoms = mmut.GetAtomSelIndex(mol.molHnd,selHnd,selindexp)
    nselatoms = selindexp.value()
  except:
    selAtoms = newPPCAtom()
    nselatoms = mol.molHnd.GetSelIndex(selHnd,selAtoms)
  print('interpetSelection nselatoms',nselatoms)

