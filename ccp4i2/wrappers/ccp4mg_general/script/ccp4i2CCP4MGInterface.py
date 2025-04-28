from inspect import getsourcefile
import functools
import glob
import os
import shutil

from PySide2 import QtCore, QtGui, QtWidgets

from ....core.CCP4MgImports import displayTableObjects
from ....core.CCP4MgImports import get_dispobj
from ....core.CCP4MgImports import global_definitions
from ....core.CCP4MgImports import MGApplication
from ....core.CCP4MgImports import mmdb2 as mmdb
from ....core.CCP4MgImports import mmut
from ....core.CCP4MgImports import MolLabel
from ....core.CCP4MgImports import point_funcs
from ....core.CCP4MgImports import pygl_coord
from ....core.CCP4MgImports import sequence_util
from ....core.CCP4MgImports import SequenceViewer


def InstallSaveToi2MenuItem(workDirectory):
  def SequenceView__init__(self,parent=None):
      SequenceViewer.SequenceView.old__init__(self,parent)
      children = self.findChildren(QtWidgets.QAction)
      for child in children:
        if str(child.text()) == "Load alignment from file":
           child.setText("Load sequence/alignment from file")

  def setDelegateACVs(self):
      atomColourVectors = []
      for sd in self.sequence_displays:
            acv = []
            if hasattr(sd.atomColourVector,"GetNumberOfColours"):
              for i in range(sd.atomColourVector.GetNumberOfColours()):
                p = sd.atomColourVector.GetRGB(i)
                a = point_funcs.doublea(0)
                a = a.frompointer(p)
                acv.append(QtGui.QColor(int(a[0]*255), int(a[1]*255), int(a[2]*255)))
            elif hasattr(sd.atomColourVector,"__str__") and hasattr(sd.atomColourVector,"seq"):
              for i in range(len(sd.atomColourVector.seq)):
                a = sd.atomColourVector.GetRGB(i)
                acv.append(QtGui.QColor(int(a[0]*255), int(a[1]*255), int(a[2]*255)))
            else:
              acv = sd.atomColourVector
              
            atomColourVectors.append(acv)
      self.sequenceView.itemDelegate().setAtomColourVectorsSize(len(atomColourVectors))
      for i in range(len(atomColourVectors)):
        self.sequenceView.itemDelegate().setAtomColourVector(atomColourVectors[i],i)
      self.sequenceView.itemDelegate().setOffsets(self.sequenceView.rowOffset,self.sequenceView.columnOffset)


  def GetColourByNucleotideAtomTable():
    colourByNucleotideAtomTable = {}
    colourByNucleotideAtomTable["A"] = QtGui.QColor(255,0,0);     # red
    colourByNucleotideAtomTable["T"] = QtGui.QColor(255,255,0);   # yellow
    colourByNucleotideAtomTable["G"] = QtGui.QColor(0,255,0);     # green
    colourByNucleotideAtomTable["C"] = QtGui.QColor(0,0,255);     # blue
    colourByNucleotideAtomTable["U"] = QtGui.QColor(255,0,255);   # magenta
    colourByNucleotideAtomTable["DA"] = QtGui.QColor(255,0,0);    # red
    colourByNucleotideAtomTable["DT"] = QtGui.QColor(255,255,0);  # yellow
    colourByNucleotideAtomTable["DG"] = QtGui.QColor(0,255,0);    # green
    colourByNucleotideAtomTable["DC"] = QtGui.QColor(0,0,255);    # blue
    colourByNucleotideAtomTable["ADE"] = QtGui.QColor(255,0,0);   # red
    return colourByNucleotideAtomTable

  def GetColourByAtomTable():
    colourByAtomTable = {}
    colourByAtomTable["A"] = QtGui.QColor(255,127,80);  # coral
    colourByAtomTable["R"] = QtGui.QColor(0,0,255);     # blue
    colourByAtomTable["N"] = QtGui.QColor(0,255,255);   # cyan
    colourByAtomTable["D"] = QtGui.QColor(255,0,0);     # red
    colourByAtomTable["C"] = QtGui.QColor(255,255,0);   # yellow
    colourByAtomTable["Q"] = QtGui.QColor(0,255,255);   # cyan
    colourByAtomTable["E"] = QtGui.QColor(255,0,0);     # red
    colourByAtomTable["G"] = QtGui.QColor(255,255,255); # white
    colourByAtomTable["H"] = QtGui.QColor(65,154,225);  # light blue
    colourByAtomTable["I"] = QtGui.QColor(255,127,80);  # coral
    colourByAtomTable["L"] = QtGui.QColor(255,127,80);  # coral
    colourByAtomTable["K"] = QtGui.QColor(0,0,255);     # blue
    colourByAtomTable["M"] = QtGui.QColor(255,255,0);   # yellow
    colourByAtomTable["F"] = QtGui.QColor(255,0,255);   # magenta
    colourByAtomTable["P"] = QtGui.QColor(255,127,80);  # coral
    colourByAtomTable["S"] = QtGui.QColor(0,255,255);   # cyan
    colourByAtomTable["T"] = QtGui.QColor(0,255,255);   # cyan
    colourByAtomTable["W"] = QtGui.QColor(255,0,255);   # magenta
    colourByAtomTable["Y"] = QtGui.QColor(255,0,255);   # magenta
    colourByAtomTable["V"] = QtGui.QColor(255,127,80);  # coral
    colourByAtomTable["X"] = QtGui.QColor(128,128,128); # grey
    colourByAtomTable["O"] = QtGui.QColor(128,128,128); # grey
    colourByAtomTable["U"] = QtGui.QColor(128,128,128); # grey
    colourByAtomTable["Z"] = QtGui.QColor(128,128,128); # grey
    return colourByAtomTable

  def alignmentToSequenceDisplay(self,new_sequences,checkStates={},align=True):
                  mappings = []
                  mappedNew = {}
                  acv = None
                  for new_seq in new_sequences:
                        if new_seq.seqtype != 'ALIGNMENT':
                              mappedNew[new_seq.name] = False
                              """
                              if not align:
                                  print "Might be nice to make AtomColourVector for",new_seq,"here!!"
                                  print new_seq.seqtype 
                                  acv = acvFromSequence(new_seq.seq,new_seq.seqtype)
                              """
                        else:
                              acv = alignmentAtomColourVector(new_seq.seq)
                  
                  contConservScores = self.computeConservationScore(new_sequences)

                  pm = global_definitions.PM("sequence_prefs")
                  continuous = False
                  if pm is not None:
                    continuous = pm.getparams()["colouring_mode"]
                  print("Got continuous",continuous)
                    
                  for seq in self.sequence_displays:
                        seq.isAligned = False
                        if seq.name+"_"+seq.chain in self.mappings:
                              mapping = self.mappings[seq.name+"_"+seq.chain]
                        else:
                              mapping = []
                        if hasattr(seq,"atomColourVector_orig"):
                              seq.atomColourVector = seq.atomColourVector_orig
                        for new_seq in new_sequences:
                              if new_seq.seqtype == 'ALIGNMENT': continue
                              key = new_seq.name
                              val = new_seq.seq
                              #print "val",val
                              if key == seq.name+"_"+seq.chain:
                                    #print key
                                    if 1 and acv:
                                       if continuous:
                                          seq.atomColourVector_orig = seq.atomColourVector
                                          seq.atomColourVectorAlign = acv
                                          seq.contConservScores = contConservScores
                                          seq.atomColourVectorAlign = sequence_util.GetSequenceConservationACVFromScores(contConservScores)
                                          seq.atomColourVector = seq.atomColourVectorAlign
                                       else:
                                          seq.atomColourVector_orig = seq.atomColourVector
                                          seq.atomColourVectorAlign = acv
                                          seq.atomColourVector = acv
                                    mapping = []
                                    mappedNew[key] = True
                                    seq.sequence = val
                                    seq.isAligned = True
                                    i = 0
                                    #print "acv",acv
                                    if not acv:
                                          for c in val:
                                                if c == '-':
                                                      mapping.append(-1)
                                                else:
                                                      mapping.append(i)
                                                      i = i + 1
                              self.mappings[seq.name+"_"+seq.chain] = mapping
                        #print mappings
                        mappings.append(mapping)

                  for new_seq in new_sequences:
                        key = new_seq.name
                        val = new_seq.seq.upper()
                        if new_seq.seqtype == 'ALIGNMENT':
                              continue
                        if not mappedNew[key]:
                              mapping = []
                              if acv is None:
                                print("acv is None")
                                nposs_nuc = val.count("A") + val.count("G") + val.count("C") + val.count("T") + val.count("U") + val.count("N");
                                if float(nposs_nuc)/len(val)<0.9:
                                  colourTable = GetColourByAtomTable()
                                  isNucleotide = False
                                else:
                                  colourTable = GetColourByNucleotideAtomTable()
                                  isNucleotide = True
 
                                acv = []
                                for seqi in val:
                                  if seqi in colourTable:
                                    acv.append(colourTable[seqi])
                                  else:
                                    acv.append(QtGui.QColor(128,128,128));
                              print("acv is ...",acv)
                              sd = SequenceViewer.SequenceDisplay(sequence=val,name=key,chain='',atomColourVector=acv,isMolData=False)
                              self.sequence_displays.append(sd)
                              self.mappings[val] = mapping
                              mappings.append(mapping)
                              if 1 and acv:
                                    sd.atomColourVector = acv

                  self.sequenceView.itemDelegate().setColumnMappingsSize(len(mappings))
                  #print "Set",len(mappings),"mappings"
                  for i in range(len(mappings)):
                    self.sequenceView.itemDelegate().setColumnMapping(mappings[i],i)
                    #print "Set",i,len(mappings[i])
                  self.resetModelandView(checkStates)

                  self.updateUDDs()
                  self.applySelectionsFromDispobjs()
                  global_definitions.HISTORY().SavePluginStatus('SequenceViewer',self.getParams())

  oldBoring_getGuiDef = displayTableObjects().GMolDisp.getGuiDef
  def getGuiDef(self,name='row',pickevent=None):
    target = get_dispobj(name=self.objectName())
    if not target:
      print('Can not find GMolDisp target',self)
      return []
  
    if name == "icon":
        labels = [self.tr('Atom labels')]
        labels.extend(MolLabel.LABEL_SELECTION_ALIAS)
        labels.extend(['sep','labels_details'])
        menu = ['visible','centre_on','flash','transparency',labels,['Symmetry','apply_symmetry','draw_central','apply_biomolecule'],'custom_drawing_style']
        menu.extend(['list','savelist'])
        menu.extend(['savetoi2'])
        if ['atom_sas','res_sas','atom_buried','res_buried','secstr','interface'].count(target.model_colour.colour_mode): menu.append('list_colour_data')
        menu.extend(['clone','delete'])
        #menu.extend(['clone','delete'])
        return menu
    else:
      return oldBoring_getGuiDef(self,name,pickevent)
  
  oldBoring_getActionDef = displayTableObjects().GMolDisp.getActionDef
  def getActionDef(self,name,**info):
   if name =='savetoi2':
     return dict (text = self.tr('Save selected atoms to ccp4i2'), slot = self.saveToi2, enabled = 1 )
   else:
     return oldBoring_getActionDef(self,name,**info)
  
  def saveToi2(self):

    target = get_dispobj(name=self.objectName())

    if not target:
      print("No target!!!!")
      return

    rv = target.parent.list_data_file(target.name)
    if rv[0] == 0:
     dropDir = os.path.join(workDirectory,"CCP4MG_FILE_DROP")
     outList = glob.glob(os.path.join(dropDir,'output*.pdb'))
     maxIndx = 0
     for f in outList:
       fpath,fname = os.path.split(f)
       maxIndx =  max(maxIndx,int(fname[6:-4]))

     dlfname = rv[1]
     fname = os.path.join(dropDir,'output'+str(maxIndx+1)+'.pdb')
     shutil.copy2(dlfname,fname)
    else:
      print("Some failure!!!!")
  
  displayTableObjects().GMolDisp._getGuiDef = displayTableObjects().GMolDisp.getGuiDef
  displayTableObjects().GMolDisp._getActionDef = displayTableObjects().GMolDisp.getActionDef
  
  displayTableObjects().GMolDisp.getGuiDef = getGuiDef
  displayTableObjects().GMolDisp.getActionDef = getActionDef
  displayTableObjects().GMolDisp.saveToi2 = saveToi2

  #SequenceViewer.SequenceView.alignmentToSequenceDisplay = alignmentToSequenceDisplay
  #SequenceViewer.SequenceView.setDelegateACVs = setDelegateACVs
  SequenceViewer.SequenceView.old__init__ = SequenceViewer.SequenceView.__init__
  SequenceViewer.SequenceView.__init__ = SequenceView__init__

SequenceViewer_initialized = 0
def SetupSequenceLoadingFromI2():
  SEQUENCE_SUFFIXES = [".pir",".fasta",".pfam",".gde",".rsf",".gcg",".cd",".amps",".gb",".msf",".clw",".afa",".seq"]
  
  def openSequence(args):
    global SequenceViewer_initialized
    if not SequenceViewer_initialized:
      initSequenceViewer()
      SequenceViewer_initialized = 1
  
    print("Treat",args,"as sequence")
    global_definitions.MAINWINDOW().sequence_dialog.loadExternalAlignment(args)
  
  def initSequenceViewer():
    if not hasattr(global_definitions.MAINWINDOW(),"sequence_dialog"):
      SequenceViewer.initializePlugin()
      SequenceViewer.handleSequenceDialog()
    global_definitions.MAINWINDOW().sequence_dialog.ClearAlign()
    global_definitions.MAINWINDOW().sequence_dialog.sequence_displays = []
    global_definitions.MAINWINDOW().sequence_dialog.resetModelandView()
    global_definitions.MAINWINDOW().sequence_dialog.updateUDDs()
    global_definitions.MAINWINDOW().sequence_dialog.applySelectionsFromDispobjs()
    global_definitions.MAINWINDOW().UnClose()
  
  def loadSequencesFromCommandLine():
  
    for f in global_definitions.HISTORY().command_line_other_files:
      suffix = os.path.splitext(f)[1]
      if suffix.lower() in SEQUENCE_SUFFIXES:
        openSequence(f)
  
  loadSequencesFromCommandLine()

@QtCore.Slot()
def saveEnsembleToI2(workDirectory):
      newManager = mmdb.Manager()
      model = mmdb.Model()
      model.thisown = 0
      newManager.AddModel(model)
      theDataobjs = global_definitions.get_dataobj(object_type='MolData')
      for obj in theDataobjs:
        if hasattr(obj,"visible") and obj.visible:
           nChains = pygl_coord.intp()
           chainTable = mmut.GetChainTable(obj.molHnd,obj.first_nmr_model,nChains)
           selHnd = obj.molHnd.NewSelection()
           obj.molHnd.SelectAtoms(selHnd,obj.first_nmr_model,"XXX",mmdb.ANY_RES,"X",mmdb.ANY_RES,"X","X","X","X","X",mmdb.SKEY_NEW )
           for dispobj in obj.get_dispobj('MolDisp'):
              if hasattr(dispobj,"visible") and dispobj.visible:
                do_selHnd = dispobj.SelHandle.getSelHnd()
                obj.molHnd.Select(selHnd,mmdb.STYPE_ATOM,do_selHnd,mmdb.SKEY_OR)
                obj.molHnd.SelectAtoms(selHnd,obj.first_nmr_model,"*",mmdb.ANY_RES,"*",mmdb.ANY_RES,"*","*","*","*","*",mmdb.SKEY_AND )
           obj.molHnd.MakeSelIndex(selHnd)
           #selindexp = pygl_coord.intp()
           #SelAtoms = mmut.GetAtomSelIndex(obj.molHnd,selHnd,selindexp)
           #print "Selected",selindexp.value(),"atoms in total",obj.name
           print("There are",nChains.value(),"chains")
           newchid = 65
           for i in range(nChains.value()):
              selHndCh = obj.molHnd.NewSelection()
              chid = mmut.getPCChain(chainTable,i).GetChainID()
              obj.molHnd.SelectAtoms(selHndCh,obj.first_nmr_model,chid,mmdb.ANY_RES,"*",mmdb.ANY_RES,"*","*","*","*","*",mmdb.SKEY_NEW )
              obj.molHnd.Select(selHndCh,mmdb.STYPE_ATOM,selHnd,mmdb.SKEY_AND )
              selindexp = pygl_coord.intp()
              SelAtoms = mmut.GetAtomSelIndex(obj.molHnd,selHndCh,selindexp)
              print("Selected",selindexp.value(),"atoms from",obj.name,chid)
              if selindexp.value()>0:
                print("Creating chain",chr(newchid))
                newChain = model.GetChainCreate(chr(newchid),True);
                print("Copying selection")
                ncopied = obj.molHnd.CopySelectedAtomsToChain(selHndCh,newChain)
                model = mmdb.Model()
                model.thisown = 0
                newManager.AddModel(model)
                print("Copied",ncopied,"to new chain")
                newchid += 1
              obj.molHnd.DeleteSelection(selHndCh)
           obj.molHnd.DeleteSelection(selHnd)
      newManager.FinishStructEdit()

      dropDir = os.path.join(workDirectory,"CCP4MG_FILE_DROP")
      outList = glob.glob(os.path.join(dropDir,'output*.pdb'))

      maxIndx = 0
      for f in outList:
        fpath,fname = os.path.split(f)
        maxIndx =  max(maxIndx,int(fname[6:-4]))

      fname = os.path.join(dropDir,'output'+str(maxIndx+1)+'.pdb')

      newManager.WriteCIFASCII(fname)


def InstallSaveEnsembleToi2MenuItem(workDir):
    mainwin = MGApplication.GetMainWindow()

    menu_defn = {}
    menu_defn['text'] = "Save all visible to ccp4i2 database "
    window_item = mainwin.addMenuDefinition(menu_defn['text'],menu_defn,'_file_save',draw_menu=True)
    window_item.triggered.connect(functools.partial(saveEnsembleToI2,workDir))

if __name__ == "__main__" or __name__ == "builtins":
  workDir = os.path.dirname(getsourcefile(lambda:0))
  InstallSaveToi2MenuItem(workDir)
  InstallSaveEnsembleToi2MenuItem(workDir)
  SetupSequenceLoadingFromI2()
