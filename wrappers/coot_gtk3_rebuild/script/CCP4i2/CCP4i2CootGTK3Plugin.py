
import coot_gui_api
import coot_gui
import coot
from gi.repository import Gtk
import gi
import os
import sys
import glob
from xml.etree import ElementTree as ET

import django
from pathlib import Path
from django.conf import settings

BASE_DIR = Path(__file__).resolve().parent.parent
sys.path.append(BASE_DIR.__str__())
BASE_DIR = Path(__file__).resolve().parent.parent.parent.parent.parent
sys.path.append(BASE_DIR.__str__())
print(sys.path)

settings.configure(
    INSTALLED_APPS=(
        'CCP4i2',
    ),
    DATABASES={
        'default': {
            'ENGINE': 'django.db.backends.sqlite3',
            'NAME': Path(os.environ['HOME']) / '.CCP4I2' / 'db' / 'database.sqlite',
            'CONN_MAX_AGE': 0
        }
    },
)
django.setup()

from dbapi.CCP4i2 import models

gi.require_version('Gtk', '3.0')

# example for main_hbox and main_toolbar:
# label = Gtk.Label(label="Test Label")
# test_button = Gtk.Button(label="Test Button")
# coot_gui_api.main_hbox().append(label) # works
# coot_gui_api.main_toolbar().append(test_button)

def fileValuesToFileInfo(fileValues):
    filePath = Path(fileValues['jobid__projectid__projectdirectory'])
    if fileValues['pathflag'] is None or fileValues['pathflag'] == 1:
        filePath = filePath / 'CCP4_JOBS' / os.path.sep.join(
            [f"job_{part}" for part in fileValues["jobid__jobnumber"].split('.')]) / fileValues['filename']
    elif fileValues['pathflag'] == 2:
        filePath = filePath / 'CCP4_IMPORTED_FILES' / fileValues['filename']
    
    return dict(path=filePath.__str__(), 
                jobNumber=fileValues['jobid__jobnumber'], 
                annotation=fileValues['annotation'],
                mimetype=fileValues['filetypeid__filetypename'],
                subType=fileValues['filesubtype'],
                content=fileValues['filecontent'])

def paramFileElementToFileInfo(paramFile):
    # Does the file have a dbFileId ?
    dbFileIdElement = paramFile.find('./dbFileId')
    
    if dbFileIdElement is not None:
        fileValues = models.Files.objects\
            .filter(fileid=dbFileIdElement.text)\
            .values(*CCP4i2Menu.fileFields)[0]
        return fileValuesToFileInfo(fileValues)

    else:
        projectElement = paramFile.find('./project')
        baseNameElement = paramFile.find('./baseName')
        relPathElement = paramFile.find('./relPath')
        subTypeElement = paramFile.find('./subType')
        contentFlagElement = paramFile.find('./contentFlag')
        annotationElement = paramFile.find('./annotation')

        fileProject = models.Projects.objects.get(
            projectid=projectElement.text)
        filePath = os.path.join(
            fileProject.projectdirectory, relPathElement.text, baseNameElement.text)
        
        subType = 0
        if subTypeElement is not None and len(subTypeElement.text.strip()) != 0: 
            subType =  int(subTypeElement.text)

        contentFlag = 0
        if contentFlagElement is not None and len(contentFlagElement.text.strip()) != 0: 
            contentFlag =  int(contentFlagElement.text)

        annotation = ""
        if annotationElement and len(annotationElement.text.strip()) > 0:
            annotation = annotationElement.text

        if paramFile.tag == "CPdbDataFile":
            mimetype="chemical/x-pdb"

        return dict(path=filePath.__str__(), 
                jobNumber="", 
                annotation=annotation,
                mimetype=mimetype,
                subType=subType,
                content=contentFlag)

class CCP4i2Menu:

    fileFields = ['filename', 'filesubtype', 'filecontent', 'fileid', 'annotation', 'pathflag', 'jobparamname', 'filetypeid__filetypename',
                  'jobid__projectid__projectname', 'jobid__projectid__projectdirectory', 'jobid__jobnumber']

    menuFileEntries = [
        {'label': "Import coordinates from project",
         'predicate': {'filetypeid__filetypename': "chemical/x-pdb"}},
        {'label': "Import 2Fo-Fc maps from project",
         'predicate': {
             'filetypeid__filetypename': "application/CCP4-mtz-map",
             'filesubtype': 1
         }},
        {'label': "Import Fo-Fc maps from project",
         'predicate': {
             'filetypeid__filetypename': "application/CCP4-mtz-map",
             'filesubtype': 2
         }},
        {'label': "Import anomalous maps from project",
         'predicate': {
             'filetypeid__filetypename': "application/CCP4-mtz-map",
             'filesubtype': 3
         }},
        {'label': "Import anomalous double difference maps from project",
         'predicate': {
             'filetypeid__filetypename': "application/CCP4-mtz-map",
             'filesubtype': 4
         }},
        {'label': "Import cif dictionary",
         'predicate': {'filetypeid__filetypename': "application/refmac-dictionary"}},
    ]

    def setProjectName(self, newName):
        print('in set projectName', newName)
        self.currentProject = models.Projects.objects.get(projectname=newName)

    def readFile(self, fileInfo):
        print('Read file ', fileInfo)
        if fileInfo['mimetype'] == "chemical/x-pdb":
            newMol = coot.read_pdb(fileInfo['path'])
            coot.set_molecule_name(newMol, f"{fileInfo['jobNumber']}: {fileInfo['annotation']}")
        elif fileInfo['mimetype'] == "application/CCP4-mtz-map":
            isDiff = 0
            if fileInfo['subType'] in [2,3,4]:
                isDiff = 1
            newMap=coot.make_and_draw_map(fileInfo['path'], 'F', 'PHI', 'PHI', 0, isDiff)
            coot.set_molecule_name(newMap, f"{fileInfo['jobNumber']}: {fileInfo['annotation']}")
        elif fileInfo['mimetype'] == "application/refmac-dictionary":
            coot.read_cif_dictionary(fileInfo['path'])

    def readValuesFile(self, fileValues):
        print('Read values file', fileValues)
        fileInfo = fileValuesToFileInfo(fileValues)
        self.readFile(fileInfo)

    def readParamFile(self, paramFile):
        print('Read param file', ET.dump(paramFile))
        fileInfo = paramFileElementToFileInfo(paramFile)
        self.readFile(fileInfo)

    def saveToDatabase(self, molNo):
        outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.pdb')))
        #print 'saveToDatabase outList:'
        maxIndx = 0
        for f in outList:
            fpath,fname = os.path.split(f)
            maxIndx =  max(maxIndx,int(fname[6:-4]))
        #print 'saveToDatabase',dropDir,maxIndx
        coordPath = os.path.normpath(os.path.join (self.dropDir,'output'+str(maxIndx+1)+'.pdb') )
        coot.save_coordinates(molNo, coordPath)
        print('Saved file to ', coordPath)
        return

    def selectAndSave(self, arg):
        coot.molecule_chooser_gui ("Molecule to capture: ", lambda imol : self.saveToDatabase(imol))
        return

    def __init__(self, jobId=None):
        print('Hello!!!!!!!!!!!!!')
        self.dropDir = os.environ['PWD']
        if jobId is not None:
            self.job = models.Jobs.objects.get(jobid=jobId)
            print(
                f"Starting coot to run job {jobId} in directory {self.job.jobDirectory}")
            self.dropDir = os.path.join(self.job.jobDirectory, "COOT_FILE_DROP")
            self.currentProject = self.job.projectid
            paramsXmlPath = os.path.join(
                self.job.jobDirectory, "input_params.xml")
            paramsXmlRoot = ET.parse(paramsXmlPath)
            # Iterate over input coordinate files
            for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/XYZIN_LIST/CPdbDataFile'):
                self.readParamFile(paramFile)
            for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/FPHIIN_LIST/CMapCoeffsDataFile'):
                self.readParamFile(paramFile)
            for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/DELFPHIIN_LIST/CMapCoeffsDataFile'):
                self.readParamFile(paramFile)
            for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/DICT'):
                self.readParamFile(paramFile)

        def populateMenuItem(w, paramDict):
            # for item in coordMenu.children():
            #    coordMenu.remove(item)
            print('In populate', paramDict)
            parentMenu = paramDict['parentMenu']
            predicate = paramDict['predicate']
            for dbFile in models.Files.objects\
                    .filter(**predicate, jobid__projectid__projectname=self.currentProject.projectname)\
                    .values(*CCP4i2Menu.fileFields):
                coordMenuItem = Gtk.MenuItem(
                    label=f"{dbFile['jobid__jobnumber']}: {dbFile['annotation']}")
                coordMenuItem.connect(
                    "activate", lambda widget, arg: self.readValuesFile(arg), dbFile)
                parentMenu.append(coordMenuItem)
                coordMenuItem.show()
            parentMenu.show()

        def populateSaveMoleculeList(widget, parentMenu):
            print('In populate MoleculeList', parentMenu)
            for iMol in range(coot.graphics_n_molecules()):
                print(iMol,coot.is_valid_model_molecule(iMol) )
                if coot.is_valid_model_molecule(iMol): 
                    menuItem = Gtk.MenuItem(
                        label=f"{iMol}: {coot.molecule_name(iMol)}")
                    menuItem.connect(
                        "activate", lambda widget, arg: self.saveToDatabase(arg), iMol)
                    parentMenu.append(menuItem)
                    menuItem.show()
            parentMenu.show()

        if coot_gui_api.main_menubar():
            main_menubar = coot_gui_api.main_menubar()
            menu = coot_gui.coot_menubar_menu("CCP4i2")

            sub_menuitem = Gtk.MenuItem(label="Import")
            menu.append(sub_menuitem)
            sub_menuitem.show()

            for menuFileEntry in CCP4i2Menu.menuFileEntries:
                menu_sub = Gtk.Menu()
                menuitem_sub = Gtk.MenuItem(label=menuFileEntry["label"])
                menuitem_sub.connect("activate", populateMenuItem,
                                     {'parentMenu': menu_sub,
                                      'predicate': menuFileEntry["predicate"]})
                menuitem_sub.set_submenu(menu_sub)
                menu.append(menuitem_sub)
                menuitem_sub.show()

            h_sep = Gtk.SeparatorMenuItem()
            menu.append(h_sep)
            h_sep.show()

            sub_menuitem = Gtk.MenuItem(label="Other")
            menu.append(sub_menuitem)
            sub_menuitem.show()

            menu_sub = Gtk.Menu()
            for project in models.Projects.objects.all():
                print('Known project', project.projectname)
                projectItem = Gtk.MenuItem(label=project.projectname)
                projectItem.connect(
                    "activate", lambda widget, arg: self.setProjectName(arg), project.projectname)
                menu_sub.append(projectItem)
                projectItem.show()
            menuitem_sub = Gtk.MenuItem(label="Other Projects")
            menuitem_sub.set_submenu(menu_sub)
            menu.append(menuitem_sub)
            menuitem_sub.show()

            h_sep = Gtk.SeparatorMenuItem()
            menu.append(h_sep)
            h_sep.show()

            sub_menuitem = Gtk.MenuItem(label="Save")
            menu.append(sub_menuitem)
            sub_menuitem.show()

            menuitem_sub = Gtk.MenuItem(label="Save to CCP4i2")
            menuitem_sub.connect(
                "activate", lambda widget, arg: self.saveToDatabase(arg), 0)
            menu.append(menuitem_sub)
            menuitem_sub.show()

            menuitem_sub = Gtk.MenuItem(label="Save mol to CCP4i2")
            menu_sub = Gtk.Menu()
            menuitem_sub.connect(
                "activate", lambda widget, arg: populateSaveMoleculeList(widget, arg), menu_sub)
            menuitem_sub.set_submenu(menu_sub)
            menu.append(menuitem_sub)
            menuitem_sub.show()

            h_sep = Gtk.SeparatorMenuItem()
            menu.append(h_sep)
            h_sep.show()

if __name__ == "__main__":
    a = CCP4i2Menu()
