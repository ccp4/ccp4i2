import sys
import os
import shutil
import glob
import functools
from inspect import getsourcefile

from ccp4i2.baselayer import QtCore, QtWidgets
from ccp4i2.core.mgimports import displayTableObjects
from ccp4i2.core.mgimports import get_dispobj
from ccp4i2.core.mgimports import global_definitions
from ccp4i2.core.mgimports import MGApplication
from ccp4i2.core.mgimports import mmdb2 as mmdb
from ccp4i2.core.mgimports import mmut
from ccp4i2.core.mgimports import MolLabel
from ccp4i2.core.mgimports import pygl_coord
from ccp4i2.core.mgimports import SequenceViewer


def InstallSaveToi2MenuItem(workDirectory):

  def SequenceView__init__(self,parent=None):
      SequenceViewer.SequenceView.old__init__(self,parent)
      children = self.findChildren(QtWidgets.QAction)
      for child in children:
        if str(child.text()) == "Load alignment from file":
           child.setText("Load sequence/alignment from file")

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
        return menu
    else:
      return displayTableObjects().GMolDisp._getGuiDef(self,name,pickevent)
  
  def getActionDef(self,name,**info):
   if name =='savetoi2':
     return dict (text = self.tr('Save selected atoms to ccp4i2'), slot = self.saveToi2, enabled = 1 )
   else:
     return displayTableObjects().GMolDisp._getActionDef(self,name,**info)
  
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
      import SequenceViewer
      SequenceViewer.initializePlugin()
      SequenceViewer.handleSequenceDialog()
    global_definitions.MAINWINDOW().sequence_dialog.ClearAlign()
    global_definitions.MAINWINDOW().sequence_dialog.sequence_displays = []
    global_definitions.MAINWINDOW().sequence_dialog.resetModelandView()
    global_definitions.MAINWINDOW().sequence_dialog.updateUDDs()
    global_definitions.MAINWINDOW().sequence_dialog.applySelectionsFromDispobjs()
  
  def loadSequencesFromCommandLine():
  
    for f in global_definitions.HISTORY().command_line_other_files:
      suffix = os.path.splitext(f)[1]
      if suffix.lower() in SEQUENCE_SUFFIXES:
        openSequence(f)
  
  loadSequencesFromCommandLine()

@QtCore.Slot(str)
def saveEnsembleToI2(workDirectory):
      newManager = mmdb.Manager()
      iModel = 0
      preamble1 = ""
      preamble2 = ""
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
                model = mmdb.Model()
                model.thisown = 0
                newManager.AddModel(model)
                iModel += 1
                newChain = model.GetChainCreate(chr(newchid),True);
                ncopied = obj.molHnd.CopySelectedAtomsToChain(selHndCh,newChain)
                print("Copied",ncopied,"to new chain")
                udd_rmsd_model=obj.molHnd.GetUDDHandle(mmdb.UDR_HIERARCHY,"mrbump_gesamt_multi_rmsd")
                if udd_rmsd_model>0:
                   udd_rmsd = mmdb.doublep()
                   obj.molHnd.GetUDData(udd_rmsd_model,udd_rmsd)
                   rmsdval = udd_rmsd.value()
                   ftype, shortName, longName = obj.filename
                   preamble1 += "REMARK PHASER ENSEMBLE MODEL "+str(iModel)+" ID "+str(rmsdval) + "\n"
                   preamble2 += "REMARK MODEL "+str(iModel)+": "+os.path.basename(shortName)+", MODEL '', CHAIN "+chid + "\n"
                else:
                   print("Clearing preambles")
                   preamble1 = None
                   preamble2 = None
                   
                sys.stdout.flush()
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

      if (preamble1 is not None) and (preamble2 is not None):
          import tempfile
          tfile = tempfile.NamedTemporaryFile(suffix=".pdb",prefix="ccp4mg"+str(os.getpid()),delete=False)
          tname = tfile.name
          tfile.close()
          os.unlink(tname)
          newManager.WritePDBASCII(tname)
          rfile = open(tname)
          fcontents = rfile.read()
          rfile.close()
          ofile = open(fname,"wb+")
          ofile.write(preamble1)
          ofile.write(preamble2)
          ofile.write(fcontents)
          ofile.close()
      else:
          newManager.WritePDBASCII(fname)


def InstallSaveEnsembleToi2MenuItem(workDir):
    mainwin = MGApplication.GetMainWindow()

    menu_defn = {}
    menu_defn['text'] = "Save all visible to ccp4i2 database "
    window_item = mainwin.addMenuDefinition(menu_defn['text'],menu_defn,'_file_save',draw_menu=True)
    window_item.triggered.connect(functools.partial(saveEnsembleToI2,workDir))

if __name__ == "__main__" or __name__ == "builtins":
  workDir = os.path.dirname(getsourcefile(lambda:0))
  print(workDir)
  InstallSaveToi2MenuItem(workDir)
  InstallSaveEnsembleToi2MenuItem(workDir)
  SetupSequenceLoadingFromI2()
