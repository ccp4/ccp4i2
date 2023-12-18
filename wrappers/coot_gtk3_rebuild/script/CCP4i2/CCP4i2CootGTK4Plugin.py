
import urllib.request
import urllib.parse
from xml.etree import ElementTree as ET
from pathlib import Path
import os
import sys
import json
import argparse

import coot
import gi
gi.require_version('Gtk', '4.0')

from gi.repository import Gtk, Gio, GLib, GObject
import pathlib
sys.path.append(pathlib.Path(__file__).absolute().parents[4].__str__())
print (sys.path)
from dbapi import CCP4DjangoApiClient

fileFields = ['filename', 'filesubtype', 'filecontent', 'fileid', 'annotation', 'pathflag', 'jobparamname', 'filetypeid__filetypename',
              'jobid__projectid__projectname', 'jobid__projectid__projectdirectory', 'jobid__jobnumber']

menuFileEntries = [
    {'label': "Import coordinates from project",
        'predicate': {'filetypeid__filetypename': "chemical/x-pdb",
                      "jobid__parentjobid__isnull": True}},
    {'label': "Import 2Fo-Fc maps from project",
        'predicate': {
            'filetypeid__filetypename': "application/CCP4-mtz-map",
            'filesubtype': 1,
            "jobid__parentjobid__isnull": True
        }},
    {'label': "Import Fo-Fc maps from project",
        'predicate': {
            'filetypeid__filetypename': "application/CCP4-mtz-map",
            'filesubtype': 2,
            "jobid__parentjobid__isnull": True
        }},
    {'label': "Import anomalous maps from project",
        'predicate': {
            'filetypeid__filetypename': "application/CCP4-mtz-map",
            'filesubtype': 3,
            "jobid__parentjobid__isnull": True
        }},
    {'label': "Import anomalous double difference maps from project",
        'predicate': {
            'filetypeid__filetypename': "application/CCP4-mtz-map",
            'filesubtype': 4,
            "jobid__parentjobid__isnull": True
        }},
    {'label': "Import cif dictionary",
        'predicate': {'filetypeid__filetypename': "application/refmac-dictionary"}},
]

def labelOfFile(fileValues):
    return f"{fileValues['jobid__jobnumber']}:{fileValues['annotation']}"


class CCP4i2Menu:

    def __init__(self, baseUrl=f"http://127.0.0.1:43434/database",
                 jobId=None, menuBar=None, app=None, projectName=None, jobNumber=None):
        self.baseUrl = baseUrl
        self.djangoClient = CCP4DjangoApiClient.CCP4DjangoApiClient(baseUrl)
        if jobId is not None:
            self.jobId = jobId
        elif projectName is not None and jobNumber is not None:
            candidateJobs = self.djangoClient.modelValues("Jobs", projectid__projectname=projectName, jobnumber=jobNumber)["results"]
            print("candidateJobs", candidateJobs)
            self.jobId = candidateJobs[0]["jobid"]
        self.app = app

        #Create the importFile action
        action = Gio.SimpleAction.new("importFile", GLib.VariantType.new('s'))
        action.connect("activate", lambda a, b: print(a, b))
        self.app.add_action(action)

        self.drawMenu(menuBar)
        print(f"jobId is {jobId}")
        if jobId is not None:
            self.loadCootJob(jobId)
        self.makeBrowserWindow()

    def makeBrowserWindow(self):
        self.dbBrowserWindow = DbBrowserWindow()
        self.dbBrowserWindow.set_hide_on_close(True)

    def fileValuesToFileInfo(self, fileValues):
        filePath = Path(fileValues['jobid__projectid__projectdirectory'])
        if fileValues['pathflag'] is None or fileValues['pathflag'] == 1:
            filePath = filePath / 'CCP4_JOBS' / os.path.sep.join(
                [f"job_{part}" for part in fileValues["jobid__jobnumber"].split('.')]) / fileValues['filename']
        elif fileValues['pathflag'] == 2:
            filePath = filePath / 'CCP4_IMPORTED_FILES' / \
                fileValues['filename']

        return dict(path=filePath.__str__(),
                    jobNumber=fileValues['jobid__jobnumber'],
                    annotation=fileValues['annotation'],
                    mimetype=fileValues['filetypeid__filetypename'],
                    subType=fileValues['filesubtype'],
                    content=fileValues['filecontent'])

    def paramFileElementToFileInfo(self, paramFile):
        # Does the file have a dbFileId ?
        dbFileIdElement = paramFile.find('./dbFileId')
        if dbFileIdElement is not None:
            fileValues = self.djangoClient.modelValues(
                "Files",
                values=fileFields,
                fileid=dbFileIdElement.text
            )["results"]
            return self.fileValuesToFileInfo(fileValues[0])

        else:
            projectElement = paramFile.find('./project')
            baseNameElement = paramFile.find('./baseName')
            relPathElement = paramFile.find('./relPath')
            subTypeElement = paramFile.find('./subType')
            contentFlagElement = paramFile.find('./contentFlag')
            annotationElement = paramFile.find('./annotation')

            fileProject = self.djangoClient.modelValues(
                "Projects",
                projectid=projectElement.text)["results"][0]

            filePath = os.path.join(
                fileProject['projectdirectory'], relPathElement.text, baseNameElement.text)

            subType = 0
            if subTypeElement is not None and len(subTypeElement.text.strip()) != 0:
                subType = int(subTypeElement.text)

            contentFlag = 0
            if contentFlagElement is not None and len(contentFlagElement.text.strip()) != 0:
                contentFlag = int(contentFlagElement.text)

            annotation = ""
            if annotationElement and len(annotationElement.text.strip()) > 0:
                annotation = annotationElement.text

            if paramFile.tag == "CPdbDataFile":
                mimetype = "chemical/x-pdb"

            return dict(path=filePath.__str__(),
                        jobNumber="",
                        annotation=annotation,
                        mimetype=mimetype,
                        subType=subType,
                        content=contentFlag)

    def readFile(self, fileInfo):
        if fileInfo['mimetype'] == "chemical/x-pdb":
            newMol = coot.read_pdb(fileInfo['path'])
            # coot.set_molecule_name(newMol, f"{fileInfo['jobNumber']}: {fileInfo['annotation']}")
        elif fileInfo['mimetype'] == "application/CCP4-mtz-map":
            isDiff = 0
            if fileInfo['subType'] in [2, 3, 4]:
                isDiff = 1
            newMap = coot.make_and_draw_map(
                fileInfo['path'], 'F', 'PHI', 'PHI', 0, isDiff)
            coot.set_molecule_name(
                newMap, f"{fileInfo['jobNumber']}: {fileInfo['annotation']}")
        elif fileInfo['mimetype'] == "application/refmac-dictionary":
            coot.read_cif_dictionary(fileInfo['path'])

    def readParamFile(self, paramFile):
        fileInfo = self.paramFileElementToFileInfo(paramFile)
        self.readFile(fileInfo)

    def loadCootJob(self, jobId):
        jobXml = self.getJobFile(jobId, "input_params.xml")
        paramsXmlRoot = ET.fromstring(jobXml)
        # Iterate over input coordinate files
        for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/XYZIN_LIST/CPdbDataFile'):
            self.readParamFile(paramFile)
        for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/FPHIIN_LIST/CMapCoeffsDataFile'):
            self.readParamFile(paramFile)
        for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/DELFPHIIN_LIST/CMapCoeffsDataFile'):
            self.readParamFile(paramFile)
        for paramFile in paramsXmlRoot.findall('./ccp4i2_body/inputData/DICT'):
            self.readParamFile(paramFile)

    def populateMenuItem(self, menu, menuFileEntry):
        files = self.djangoClient.modelValues(
            "Files", values=fileFields, **menuFileEntry["predicate"])["results"]
        for aFile in files:
            menuItem = Gio.MenuItem.new(labelOfFile(aFile))
            menuItem.set_action_and_target_value(
                "app.importFile", GLib.Variant.new_string(aFile["fileid"]))
            menu.append_item(menuItem)

    def setFilePredicate(self, clickedButton):
        jsonString = clickedButton.get_action_target_value().get_string()
        predicate = json.loads(jsonString)
        print("predicate", predicate)
        self.dbBrowserWindow.setFilePredicate(self.djangoClient, predicate)

    def makePopoverPanel(self):
        popoverPanel = Gtk.Popover()
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        popoverPanel.set_child(vbox)
        for menuFileEntry in menuFileEntries:
            browseButton = Gtk.Button(label=menuFileEntry["label"])
            browseButton.set_action_target_value(GLib.Variant.new_string(json.dumps(menuFileEntry["predicate"])))
            browseButton.connect("clicked", self.setFilePredicate)
            vbox.append(browseButton)
        return popoverPanel

    def drawMenu(self, menuBar):
        menu = Gio.Menu.new()

        section_import = Gio.Menu.new()
        section_other = Gio.Menu.new()
        section_save = Gio.Menu.new()

        section_import.connect("notify", lambda a: print(a))
        for menuFileEntry in menuFileEntries:
            menuForThisModel = Gio.Menu.new()
            self.populateMenuItem(menuForThisModel, menuFileEntry)
            section_import.append_submenu(
                menuFileEntry["label"], menuForThisModel)

        menu.append_section("Import", section_import)
        menu.append_section("Other",  section_other)
        menu.append_section("Save",   section_save)

        popover = Gtk.PopoverMenu()
        popover.set_menu_model(menu)

        popoverPanel = self.makePopoverPanel()

        ccp4i2_menu_button = Gtk.MenuButton(label="CCP4i2")
        ccp4i2_menu_button.connect(
            "activate", lambda a: print(a))
        ccp4i2_menu_button.set_popover(popoverPanel)

        menuBar.append(ccp4i2_menu_button)

    def getJobFile(self, jobId, fileName="input_params.xml"):
        encodedQuery = urllib.parse.urlencode(
            dict(jobId=jobId, fileName=fileName))
        req = urllib.request.Request(
            f"{self.baseUrl}/getJobFile?{encodedQuery}")
        print(f"Req is {req.full_url}")
        with urllib.request.urlopen(req) as response:
            fileContent = response.read()
            return fileContent


class DbItem(GObject.Object):
    label = GObject.Property(type=str)
    jsonRep = GObject.Property(type=str)
    def __init__(self, name, jsonRep):
        super().__init__()
        self.name = name
        self.jsonRep = jsonRep

class DbBrowserWindow(Gtk.Window):
    def __init__(self, *args, **kwargs):
        self.listStore = Gio.ListStore() # Label and json encoding of data
        super(DbBrowserWindow, self).__init__(*args, **kwargs)
        singleSelection = Gtk.SingleSelection()
        singleSelection.set_model(self.listStore)
        treeView = Gtk.ListView(model=singleSelection)

        factory = Gtk.SignalListItemFactory()

        def f_setup(fact, item):
            label = Gtk.Label(halign=Gtk.Align.START)
            label.set_selectable(True)
            item.set_child(label)

        factory.connect("setup", f_setup)

        def f_bind(fact, item):
            item.get_child().set_label(item.get_item().name)

        factory.connect("bind", f_bind)    

        treeView.set_factory(factory)
        sw = Gtk.ScrolledWindow()
        sw.set_min_content_height(500)
        sw.set_min_content_width(500)
        sw.set_child(treeView)

        self.set_child(sw)
        self.connect("close-request",self.close_request)

    def close_request(self, userData):
        print("In destroy")    
        self.hide()
        return True
    
    def setFilePredicate(self, djangoClient, predicate):
        self.listStore.remove_all()
        dbEntities = djangoClient.modelValues("Files", values=fileFields, **predicate)["results"]
        self.set_title(f"Files of type {predicate['filetypeid__filetypename']}")
        for dbEntity in dbEntities:
            self.listStore.append(DbItem(labelOfFile(dbEntity), json.dumps(dbEntity)))
        self.present()

class MainWindow(Gtk.ApplicationWindow):

    def __init__(self, *args, jobId=None, projectName=None, jobNumber=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.header = Gtk.HeaderBar()
        self.set_titlebar(self.header)
        self.hBox = Gtk.Box()
        self.header.pack_start(self.hBox)
        # Add menu button to the header bar
        self.ccp4i2Menu = CCP4i2Menu(
            menuBar=self.hBox, app=self.get_application(), jobId=jobId, projectName=projectName, jobNumber=jobNumber)        

class MyApp(Gtk.Application):
    def __init__(self, jobId=None, projectName=None, jobNumber=None, **kwargs):
        self.jobId = jobId
        self.jobNumber = jobNumber
        self.projectName = projectName
        super().__init__(**kwargs)
        self.connect('activate', self.on_activate)

    def on_activate(self, app):
        self.win = MainWindow(application=app, jobId=self.jobId, projectName=self.projectName, jobNumber=self.jobNumber)
        self.win.present()


if __name__ == "__main__":
    print("__file__",pathlib.Path(__file__).absolute())
    sys.path.append(pathlib.Path(__file__).absolute().parents[4])
    parser = argparse.ArgumentParser(
        description='CCP4i2 plugin for GTK4Coot')
    parser.add_argument('--baseUrl',
                        default="http://127.0.0.1:43434/database",
                        help='baseUrl of the HTTP server')
    parser.add_argument('--jobId',
                        help='id of a coot or moorhen job')
    parser.add_argument('--projectName',
                        help='id of a coot or moorhen job')
    parser.add_argument('--jobNumber',
                        help='id of a coot or moorhen job')
    args = parser.parse_args()

    app = MyApp(projectName=args.projectName, jobNumber=args.jobNumber, jobId=args.jobId)
    for arg in ["--projectName", "--jobId", "--jobNumber"]:
        if arg in sys.argv:
            iArg = sys.argv.index(arg)
            sys.argv.pop(iArg)
            sys.argv.pop(iArg)

    app.run(sys.argv)

