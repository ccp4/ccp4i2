"""
Copyright (C) 2001-2008 University of York, CCLRC
Copyright (C) 2009-2011 University of York
Copyright (C) 2012 STFC
March 2012 Liz Potterton Copied from CCP4mg model_selection.py
"""

import re

from .CCP4ErrorHandling import CException
from .CCP4MgImports import mmdb2 as mmdb


class SelectionParser:

  ERROR_CODES = { 101: { 'description' : 'Error in selection command' }}

  def tokenise(self,inp_str,report=1):
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
      
    rx = re.compile(r'\{| and | or | xor | excl | not |\}')
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

  def getlevels(self,ops):
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


class Cselect_tree:

  aminoResidues = ["GLY","ALA","VAL","PRO","SER","THR","LEU","ILE","CYS","ASP","GLU","ASN","GLN","ARG","LYS","MET","MSE","HIS","PHE","TYR","TRP","HCS","ALO","PDD","UNK"]
  metalElements = ["LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "UN", "UU", "UB", "UQ", "UH", "UO"];
  saccharideResidues = ["BMA","MAN","NAG","GLC","BGC","GCS","GAL","NGA","MGC","NG1","NG6","A2G","6MN","GLP","GP4","BEM","KDN","XLS","CXY","RBY","TDX","XYL","XYS","XYP","FCA","FCB","FUC","GTR","ADA","DGU","SI3","NCC","NGF","NGE","NGC","GC4","GCD","GCU","GCV","GCW","IDS","REL","IDR","SIA"];
  nucleicResidues = [ "DC","DT","U","T","C","PSU","DA","DG","A","G","I" ]
  nucleotideResidues = [ "Ad","Cd","Gd","Td","ADE","CYT","GUA","INO","THY","URA","AMP","ADP","ATP","CMP","CDP","CTP","GMP","GDP","GTP","TMP","TDP","TTP" ]
  waterResidues = [ "HOH","H2O","WAT","SOL","DOD","D2O" ]
  metalResidues = [ "LI","BE","NA","MG","AL","K","CA","SC","TI","V","MN","FE","CO","NI","CU","ZN","GA","RB","SR","Y","ZR","NB","MO","TC","RU","RH","PD","AG","CD","IN","SN","SB","CS","BA","LA","CE","PR","ND","PM","SM","EU","GD","TB","DY","HO","ER","TM","YB","LU","HF","TA","W","RE","OS","IR","PT","AU","HG","TL","PB","BI","PO","FR","RA","AC","TH","PA","U","NP","PU","AM","CM","BK","CF","ES","FM","MD","NO","LR","RF","DB","SG","BH","HS","MT","UN","UU","UB","UQ","UH","UO" ]
  soluteResidues = [ "SUL","SO4","NO3","MOH","EOH","GOL","ACT","CL","BR","PG4","PG5","PG6","1PE","2PE","7PE","PE3","PE4","PE5","PE6","PE7","PGE","XPE","C10","CE1","CE9","CXE","EDO","N8E","P33","P4C","12P","15P","DMS","IOD","MRD","BE7","MHA","BU3","PGO","BU2","PDO","BU1","1BO","TFP","DHD","PEU","TRS","TAU","SBT","SAL","MPD","IOH","IPA","BU2","PIG","B3P","BTB","B3P","NHE","C8E","OTE","PE8","2OS","1PS","CPS","DMX","MPO","DXG","CM5","ACA","ACN","CCN","DR6","NH4","AZI","LAK","BCN","BRO","CAC","CBX","FMT","ACY","CBM","CLO","FCL","CIT","3CO","NCO","CU1","CYN","CYN","MA4","BTC","TAR","MTL","DPR","SOR","SYL","DDQ","DMF","DIO","DOX","SDS","EEE","EGL","FLO","TRT","FCY","FRU","GBL","GPX","HTO","HTG","B7G","16D","HEZ","IDO","ICI","ICT","TLA","LDA","MRY","BEQ","C15","MG8","POL","JEF","DIA","IPH","PIN","CRY","PGR","PGQ","SPD","SPK","SPM","TBU","TMA","TEP","SCN","ETF","144","UMQ","URE","CN","EPE","HEPES","HEPE","MES","PEG","ETA","HED","IMD","BO3","IUM","PO4" ]

  def __init__(self,name,parent=None,mol=None,dispobj=''):
    self.name = name
    self.parent = parent
    self.mol = mol
    self.dispobj = dispobj
    self.element = None
    self.op = None
    self.children = []
    self.selHnd = None
    self.com_string = ""

  def import_command_string(self,elements,ops,levels):

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

  def __del__(self):
    for child in self.children:
      del child

  def skey(self):
    if self.op == 'and':
      return mmdb.SKEY_AND
    elif self.op == 'or':
      return mmdb.SKEY_OR
    elif self.op == 'xor':
      return mmdb.SKEY_XOR
    elif self.op == 'excl':
      return mmdb.SKEY_CLR

  def apply(self):
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

  def new_selHnd(self):
    selHnd = self.mol.molHnd.NewSelection()
    return selHnd

  def interpret(self,skey,element):
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
      if rv == 0:
        return [rv,self.selHnd]
      else:
        return [rv,'Error applying selection CID ' + element]

  def merge_selection(self,skey,selHnd1,selHnd2):
    #print 'CSelectTree.merge_selection',skey,selHnd1,selHnd2
    # Beware Select_propagate does not return a status flag
    #print "SelHandle.merge_selection skey",skey
    self.mol.molHnd.Select(selHnd1,mmdb.STYPE_ATOM,selHnd2,skey)
    return [0,selHnd1]

  def delete_selhnd(self,selHnd):
    if selHnd >= 0: self.mol.molHnd.DeleteSelection(selHnd)
    selHnd = -1
