from __future__ import print_function

from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils

from ccp4i2.baselayer import QtCore
import os,re,time,sys

class coot_rebuild(CPluginScript):
#class coot_rebuild(CInternalPlugin):

    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Rebuild model with coot'     # A short title for gui menu
    TASKNAME = 'coot_rebuild'                                  # Task name - should be same as class name
    TASKCOMMAND = 'coot'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    ERROR_CODES = {  200 : { 'description' : 'Coot exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        self.dropDir = os.path.join(self.workDirectory,'COOT_FILE_DROP')
        if not os.path.exists(self.dropDir):
          try:
            os.mkdir(self.dropDir)
          except:
            self.dropDir = self.workDirectory
            print('Could not make dropDir reset to',self.dropDir)

        '''
        # Look for a state file with name like the input coordinate file -
        # only works if the coord file is from coot
        if not self.container.inputData.COOTSTATEFILE.isSet() and self.container.inputData.XYZIN.fullPath.count('COOT_FILE_DROP')>0:
          statePath,fileName = os.path.split(self.container.inputData.XYZIN.fullPath)
          statePath = os.path.join(statePath,re.sub('output','state',fileName))
          print 'coot_rebuild derived state file',statePath, os.path.exists(statePath)
          if os.path.exists(statePath):
              self.container.inputData.COOTSTATEFILE.setFullPath(statePath)
        '''
        # Copy a startup state file to the drop directory
        if self.container.inputData.COOTSTATEFILE.isSet(): self.copyStateFile()

        '''
          import shutil
          try:
            shutil.copyfile(self.container.inputData.COOTSTATEFILE.__str__(),os.path.join(self.dropDir,'0-coot-history.py'))
          except:
            pass
        '''

        # Make a script file with additional menu options to save to i2
        self.cootScriptPath = os.path.normpath(os.path.join(self.workDirectory,'script.py'))
        if sys.platform == 'win32':
          self.cootScriptPath = re.sub(r'\\\\',r'\\',self.cootScriptPath)
        # Declare script text then re.sub in the variables
        if sys.platform == "win32":
          i2dir = CCP4Utils.getCCP4I2Dir().replace('\\','/')
        else:
          i2dir = CCP4Utils.getCCP4I2Dir()
        script = """
from __future__ import division
import os, sys
import __builtin__
pythonDefsFile = os.path.normpath('"""+i2dir+"""')
sys.path.append(pythonDefsFile)
#Here provide that module with "global function calls it will need
for cootFunction in [coot_menubar_menu, add_simple_coot_menu_menuitem, read_pdb, set_molecule_name, make_and_draw_map, auto_read_make_and_draw_maps, read_cif_dictionary, save_coordinates, save_state_file_py, molecule_chooser_gui, residue_centre, set_rotation_centre, handle_read_draw_molecule_with_recentre, refine_zone, model_molecule_number_list, interesting_things_gui, file_to_preferences, parse_wwpdb_validation_xml, sort_subgroups, validation_to_gui, add_status_bar_text, molecule_name]:
    setattr(__builtin__, cootFunction.__name__, cootFunction)

from wrappers.coot_rebuild.script import ccp4i2CootInterface
print('Managed to load ccp4i2CootInterface')
ccp4i2Interface = ccp4i2CootInterface.ccp4i2CootInterface(dropDir=r'"""+self.dropDir+"""')
print(ccp4i2CootInterface, ccp4i2Interface)
try:
    ccp4i2Interface.installMenus()
except Exception as e:
    print(e)
"""

        if self.container.inputData.USEKEYBINDINGS.isSet() and self.container.inputData.USEKEYBINDINGS:
            script+="""
file_to_preferences('template_key_bindings.py')
"""

        if self.container.inputData.XYZIN_LIST.isSet():
            try:
                iFile = 1
                for XYZIN in self.container.inputData.XYZIN_LIST:
                    if os.path.isfile(XYZIN.__str__()):
                        if XYZIN.__str__().lower().endswith("cif"):
                            script += "try:\n"
                            script += ("  filePath = r'"+XYZIN.__str__()+"'\n")
                            script += ("  MolHandle_"+str(iFile)+"=read_pdb(filePath)\n")
                            script += ("  ccp4i2Interface.patchMoleculeName(MolHandle_"+str(iFile)+", filePath, True)\n")
                            script += "except Exception as err:\n  print('Error {} loading coordinates {}'.format(err, filePath))\n  pass\n"
                        else:
                            script += "try:\n"
                            script += ("  filePath = r'"+XYZIN.__str__()+"'\n")
                            script += ("  MolHandle_"+str(iFile)+"=read_pdb(filePath)\n")
                            script += ("  ccp4i2Interface.patchMoleculeName(MolHandle_"+str(iFile)+", filePath)\n")
                            script += "except Exception as err:\n  print('Error {} loading coordinates {}'.format(err, filePath))\n  pass\n"
                    else:
                        print('coot_rebuild.makeCommandAndScript XYZIN does not exist:',XYZIN.__str__())
                    iFile += 1
            except:
                #an issue with the existence of files
                pass
        if self.container.inputData.FPHIIN_LIST.isSet():
            try:
                iFile = 1
                for FPHIIN in self.container.inputData.FPHIIN_LIST:
                    print(' reading file number ' + str ( iFile ))
                    if os.path.isfile(FPHIIN.__str__()):
                        script += "try:\n"
                        script += ("  filePath = r'"+FPHIIN.__str__()+"'\n")
                        script += ("  MapHandle_"+str(iFile)+"=make_and_draw_map(filePath, 'F', 'PHI', 'PHI', 0, 0)\n")
                        script += ("  ccp4i2Interface.patchMoleculeName(MapHandle_"+str(iFile)+", filePath)\n")
                        script += "except Exception as err:\n  print('Error {} loading map {}'.format(err, filePath))\n  pass\n"
                    else:
                        print('coot_rebuild.makeCommandAndScript FPHIIN does not exist:',FPHIIN.__str__())
                    iFile += 1
            except:
                print(' Exception ')
                #an issue with the existence of files
                pass
        if self.container.inputData.DELFPHIIN_LIST.isSet():
            try:
                iFile = 1
                for DELFPHIIN in self.container.inputData.DELFPHIIN_LIST:
                    print(' reading diff file number ' + str ( iFile ))
                    if os.path.isfile(DELFPHIIN.__str__()):
                        script += "try:\n"
                        script += ("  filePath = r'"+DELFPHIIN.__str__()+"'\n")
                        script += ("  DifMapHandle_"+str(iFile)+"=make_and_draw_map(filePath, 'F', 'PHI', 'PHI', 0, 1)\n")
                        script += ("  ccp4i2Interface.patchMoleculeName(DifMapHandle_"+str(iFile)+",filePath)\n")
#Make anomolous difference maps white.
                        if DELFPHIIN.subType == 3:
                            script += ("  set_map_colour(DifMapHandle_"+str(iFile)+",0.75,0.9,0.75)\n")
                        script += "except Exception as err:\n  print('Error {} loading difmap {}'.format(err, filePath))\n  pass\n"
                    else:
                        print('coot_rebuild.makeCommandAndScript FPHIIN does not exist:',DELFPHIIN.__str__())
                    iFile += 1
            except:
                print(' Exception ')
                #an issue with the existence of files
                pass

        if self.container.inputData.DELFPHIINANOM_LIST.isSet():
            try:
                iFile = 1
                for DELFPHIIN in self.container.inputData.DELFPHIINANOM_LIST:
                    print(' reading anomalous diff file number ' + str ( iFile ))
                    if os.path.isfile(DELFPHIIN.__str__()):
                        script += "try:\n"
                        script += ("  filePath = r'"+DELFPHIIN.__str__()+"'\n")
                        script += ("  DifMapHandle_"+str(iFile)+"=make_and_draw_map(filePath, 'F', 'PHI', 'PHI', 0, 1)\n")
                        script += ("  ccp4i2Interface.patchMoleculeName(DifMapHandle_"+str(iFile)+",filePath)\n")
#Make anomolous difference maps white.
                        if DELFPHIIN.subType == 3:
                            script += ("  set_map_colour(DifMapHandle_"+str(iFile)+",0.75,0.9,0.75)\n")
                        script += "except Exception as err:\n  print('Error {} loading difmap {}'.format(err, filePath))\n  pass\n"
                    else:
                        print('coot_rebuild.makeCommandAndScript FPHIIN does not exist:',DELFPHIIN.__str__())
                    iFile += 1
            except:
                print(' Exception ')
                #an issue with the existence of files
                pass
                
        if self.container.inputData.COOTSCRIPTFILE.isSet():
            scriptLines = open(self.container.inputData.COOTSCRIPTFILE.fullPath.__str__()).readlines()
            if len(scriptLines)>0:
              script += 'try:\n'
              for line in scriptLines:
                  #MN from __future__ imports have to be at top of a module
                  if not "from __future__ import" in line:
                    #SJM - horrible kludge to deal with hopefully temporary problem with the Python in CCP4 Coot.
                    #SJM - horrible kludge to deal with some filenames not being Windows compatible
                    if "six.moves" in line:
                        script += ('    '+line.replace("from six.moves import","from rdkit.six.moves import") + '\n')
                    elif "with open(" in line and "win" in sys.platform:
                        script += ('    '+line.replace("with open(","with open(r") + '\n')
                    elif "with gzip.open(" in line and "win" in sys.platform:
                        script += ('    '+line.replace("with gzip.open(","with gzip.open(r") + '\n')
                    elif "handle_read_draw_probe_dots_unformatted(" in line and "win" in sys.platform:
                        script += ('    '+line.replace("handle_read_draw_probe_dots_unformatted(","handle_read_draw_probe_dots_unformatted(r") + '\n')
                    elif "os.remove (" in line and "win" in sys.platform:
                        script += ('    '+line.replace("os.remove (","os.remove (r") + '\n')
                    else:
                        script += ('    '+line + '\n')
              script += 'except:\n    pass\n'

        CCP4Utils.saveFile(self.cootScriptPath,script)

        clArgs = ['--no-state-script','--python']

        #clArgs = ['--python','--pdb',self.container.inputData.XYZIN.fullPath.__str__()]



        dict_is_meaningful = True
        if self.container.inputData.DICT.isSet():

            try:
                from Bio.PDB.MMCIF2Dict import MMCIF2Dict

                mmcif_dict = MMCIF2Dict ( str ( self.container.inputData.DICT.fullPath.__str__() ) )
                lib_name    = mmcif_dict['_lib_name']
                lib_version = mmcif_dict['_lib_version']
                lib_update  = mmcif_dict['_lib_update']

                print(lib_name[0])
                print(lib_version[0])
                print(lib_update[0])
                #MN This test does not work on dicts made by ACEDRG, which show up as ??? for these properties
                #if lib_name[0] == '?' and lib_version[0] == '?' and lib_update[0] == '?' :
                #    dict_is_meaningful = False
            except:
                print('Bio python probably not available in this build of ccp4')

        if self.container.inputData.DICT.isSet() and dict_is_meaningful :
            clArgs += ['--dictionary']
            clArgs += [self.container.inputData.DICT.fullPath.__str__()]


        #MN Please talk to me before changing the below.  COOTSTATEFILE has almost no place in how
        #i2 is used, but is incorporated into script.py above.
        clArgs += ['--script',self.cootScriptPath ]

        '''if  self.container.inputData.COOTSTATEFILE.exists():
          from core import CCP4Utils
          contents = CCP4Utils.readFile(self.container.inputData.COOTSTATEFILE.__str__())
          if len(contents)>5:
            clArgs += ['--script',self.container.inputData.COOTSTATEFILE.__str__()]
          else:
            clArgs += ['--script',self.cootScriptPath ]
        else:
          clArgs += ['--script',self.cootScriptPath ]
        '''

        self.appendCommandLine(clArgs)
        # Use Qt class to watch the drop directory
        self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=self)
        self.fileSystemWatcher.addPath(self.dropDir)
        self.fileSystemWatcher.directoryChanged.connect(self.handleFileDrop)

        return CPluginScript.SUCCEEDED


    def numberOfOutputFiles(self):
        import glob
        outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.pdb')))
        outList += glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.cif')))
        #print 'numberOfOutputFiles outList',os.path.join(self.dropDir,'output*.pdb'),outList
        #print 'numberOfOutputFiles xmlList',glob.glob(os.path.normpath(os.path.join(self.workDirectory,'*.xml')))
        maxIndx = 0
        for f in outList:
           fpath,fname = os.path.split(f)
           #print 'numberOfOutputFiles  fpath,fname', fpath,fname
           maxIndx =  max(maxIndx,int(fname[6:-4]))
        return maxIndx

    @QtCore.Slot(str)
    def handleFileDrop(self,directory):
        import time,glob
        print('coot_rebuild',time.time())
        print('coot_rebuild',glob.glob(os.path.join(self.workDirectory,'*.*')))
        #print 'handleFileDrop',directory
        #Note that I don't copy the file to the appropriate xyzout filename here, since the file may not yet
        #be closed and/or flushed


    def processOutputFiles(self):
        try:
            # First up import PDB files that have been output

            import os, glob, shutil
            outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.pdb')))
            outList += glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.cif')))

            xyzoutList = self.container.outputData.XYZOUT
            for outputPDB in outList:
                fname = os.path.split(outputPDB)[1]
                iFile = int(fname[6:-4])
                newPath = str(xyzoutList[iFile].fullPath)
                if fname.endswith(".cif") and newPath.endswith(".pdb"):
                    newPath = newPath[:-4] + ".cif"
                    xyzoutList[iFile].setFullPath(newPath)
                os.rename(outputPDB, newPath)
                xyzoutList[iFile].annotation = "Coot output file number"+str(iFile)
                xyzoutList[iFile].subType = 1
            #Here truncate the xyzoutList back to the numberof files that were actually found
            xyzoutList = xyzoutList[0:len(outList)]

            #Ligand builder places output cifs in the coot-cif directory as prorg-out.cif
            #'prodrgify this residue' places output cifs in the coot-cif directory as prodrg-???.cif
            #pyrogen create "TLC"_pyrogen.cif
            cifOutList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'coot-ccp4', 'prodrg-*.cif')))
            cifOutList += glob.glob(os.path.normpath(os.path.join(self.workDirectory,'coot-ccp4', 'prodrg-*.cif')))
            cifOutList += glob.glob(os.path.normpath(os.path.join(self.workDirectory,'*pyrogen.cif')))
            cifOutList += glob.glob(os.path.normpath(os.path.join(self.workDirectory,'acedrg-*.cif')))

            dictoutList = self.container.outputData.DICTOUT
            for iFile, outputCIF in enumerate(cifOutList):
                fpath,fname = os.path.split(outputCIF)
                os.rename(outputCIF, dictoutList[iFile].fullPath.__str__())
                if 'acedrg' in fname: annotation='Coot/Acedrg created geometry for ligand'
                elif 'pyrogen' in fname: annotation='Coot/Pyrogen created geometry for ligand'
                elif 'prodrg' in fname: annotation='Coot/Prodrg created geometry for ligand'
                else: annotation='Coot/Prodrg created geometry for ligand'
                dictoutList[iFile].annotation = annotation
            #Here truncate the dictoutList back to the numberof files that were actually found
            dictoutList = dictoutList[0:len(cifOutList)]

            # Create a trivial xml output file
            from lxml import etree
            self.xmlroot = etree.Element('coot_rebuild')
            e = etree.Element('number_output_files')
            e.text = str(self.numberOfOutputFiles())
            e = etree.Element('number_output_dicts')
            e.text = str(len(dictoutList))
            self.xmlroot.append(e)

            #Separate out here activity to attempt merge into project dictionary....this seems flakey,
            #but is needed for ongoing work, so I ammaking it give an report a warning in case of failure, rather than
            #offer the sad face ofdoom
            try:
                for dictFile in dictoutList:
                    try:
                        self.mergeDictToProjectLib(fileName=dictFile.__str__())
                    except:
                        self.addReportWarning('mergeDictToProjectLib raised exception: Does not compromise output Dict')

                    ligNodes = self.xmlroot.xpath('//LIGANDS')
                    if len(ligNodes) == 0: ligNode = etree.SubElement(self.xmlroot,'LIGANDS')
                    else: ligNode = ligNodes[0]
                    try:
                        annotation='Coot/Prodrg created geometry for'
                        for item in dictFile.fileContent.monomerList:
                            lig =  etree.SubElement(ligNode,'ligand')
                            lig.text = str(item.three_letter_code)
                            annotation += (' ' + str(item.three_letter_code))
                        dictFile.annotation = annotation
                    except:
                        self.addReportWarning('fileContent.monomerList raised exception: Does not compromise output Dict')
            except:
                self.addReportWarning('failed elsewhere in merging/analysing dicts: Does not compromise output Dict')
        except:
            self.appendErrorReport(202,'Data harvesting failed')

        CCP4Utils.saveEtreeToFile(self.xmlroot,self.makeFileName('PROGRAMXML'))
        if ( len(outList) + len(cifOutList) ) > 0:
          return CPluginScript.SUCCEEDED
        else:
          return CPluginScript.MARK_TO_DELETE

    def clearCootWorkingDir(self):
        # Remove the working directory state and history files
        import glob
        zeroFileList = glob.glob(os.path.normpath(os.path.join(self.projectDirectory(),'CCP4_COOT','0-coot*')))
        for filn in zeroFileList:
            os.remove(filn)

    def copyStateFile(self):
        newText = ''
        text = CCP4Utils.readFile(self.container.inputData.COOTSTATEFILE.fullPath.__str__())
        for line in text.split('\n'):
          if line.count('handle-read-draw-molecule'):
            newText = newText + '(handle-read-draw-molecule "'+self.container.inputData.XYZIN_LIST[0].__str__()+'" 1)\n'
          else:
            newText = newText + line +'\n'
        CCP4Utils.saveFile(os.path.normpath(os.path.join(self.dropDir,'0-coot-history.scm')),text)

    def addReportWarning(self, text):
        from lxml import etree
        warningsNode = None
        warningsNodes = self.xmlroot.xpath('//Warnings')
        if len(warningsNodes) == 0: warningsNode = etree.SubElement(self.xmlroot, 'Warnings')
        else: warningsNode = warningsNodes[0]
        warningNode = etree.SubElement(warningsNode,'Warning')
        warningNode.text = text
