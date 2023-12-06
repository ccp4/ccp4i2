from __future__ import print_function
import gtk
from core import CCP4Modules
from core import CCP4Utils
from wrappers.SyncToDjango.script import CCP4i2DjangoSession
import os
from xml.etree import ElementTree as ET

class DjangoInteractionWidget:

    # close the window and quit
    def delete_event(self, widget, event, data=None):
         gtk.main_quit()
         return False
    
    def __init__(self):
        # Create a new window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.vbox = gtk.VBox()
        self.window.add(self.vbox)
        self.vbox.pack_start(gtk.Label("Username:"), expand=False)
        self.uname=gtk.Entry()
        self.vbox.pack_start(self.uname, expand=False)
        self.vbox.pack_start(gtk.Label("Password:"), expand=False)
        self.pword=gtk.Entry()
        self.pword.set_visibility(gtk.FALSE)
        self.vbox.pack_start(self.pword, expand=False)
        updateButton = gtk.Button("Update")
        updateButton.connect("clicked",self.loadFromDjango,"Woohoo")
        self.vbox.pack_start(updateButton, expand=False)
        self.vbox.pack_start(gtk.Label("Filter:"), expand=False)
        self.filterEntry = gtk.Entry()
        self.filterEntry.connect("changed", self.filterChanged)
        self.vbox.pack_start(self.filterEntry, expand=False)
        
        bulkLoadLigandButton = gtk.Button("Download Ligand Pipelines")
        bulkLoadLigandButton.connect("clicked",self.bulkLoadFromDjango,"LigandPipeline")
        self.vbox.pack_start(bulkLoadLigandButton, expand=False)

        bulkLoadRefmacButton = gtk.Button("Download Refmac runs")
        bulkLoadRefmacButton.connect("clicked",self.bulkLoadFromDjango,"REFMAC")
        self.vbox.pack_start(bulkLoadRefmacButton, expand=False)

        self.scrolledWindow = gtk.ScrolledWindow()
        self.vbox.add(self.scrolledWindow)
        self.window.set_title("Basic TreeView Example")

        self.window.set_size_request(300, 600)

        self.window.connect("delete_event", self.delete_event)
        self.window.show_all()

    def filterChanged(self, inEntry):
        self.filteredTreestore.refilter()

    def bulkLoadFromDjango(self, whichButton, userData):
        print("Bulk load Ligand")
        print("Self", self)
        print("whichButton",whichButton)
        print("userData",userData)
        print("Done")
        selection = self.treeview.get_selection()
        print("selection", selection)
        model, selectedPaths = selection.get_selected_rows()
        #print "selectedPaths", selectedPaths
        for selectedPath in selectedPaths:
            #print "model, selectedPath", model, selectedPath
            treeiter = model.get_iter(selectedPath)
            #print "treeiter", treeiter
            if treeiter != None:
                projectId = model.get_value(treeiter, 0)
                #print "projectId", projectId
                projectName = model.get_value(treeiter, 1)
                #print "projectName", projectName
                self.addProjectToCoot(projectId, projectName, userData)
                #print "You selected", projectName, projectId
    
    def loadFromDjango(self, whichButton, userData):
        username = self.uname.get_text()
        password = self.pword.get_text()
        
        filename = str(os.path.join(CCP4Utils.getDotDirectory(),'configs','serverSetup.params.xml'))
        if os.path.exists(filename):
            serversEtree = ET.parse(filename)
        else:
            filename = str(os.path.join(CCP4Utils.getCCP4I2Dir(),'local_setup','serverSetup.params.xml'))
            serversEtree = ET.parse(filename)
        
        allServerNames = [node.text for node in serversEtree.findall(".//CHostname")]
        archiveServerNames = [name for name in allServerNames if "ManageCCP4i2Archive" in name]
        self.serverURL = "http://"+archiveServerNames[0]
        self.djangoSession = CCP4i2DjangoSession.CCP4i2DjangoSession(self.serverURL, username, password)
        response = self.djangoSession.getURLWithValues(self.serverURL+"?listProjects",{})
        responseText = response.read()
        import json
        self.projectList = json.loads(responseText)
        
        # create a TreeStore with one string column to use as the model
        self.treestore = gtk.TreeStore(str, str, str, str)
        self.filteredTreestore = self.treestore.filter_new()
        self.filteredTreestore.set_visible_func(self.visible_func, "Erm")

        # we'll add some data now - 4 rows with 3 child rows each
        idToIterMapping = {}

        for project in self.projectList:
            piter = self.treestore.append(None, project)
        # create the TreeView using treestore
        self.treeview = gtk.TreeView(self.filteredTreestore)
        self.treeview.get_selection().set_mode(gtk.SELECTION_MULTIPLE)

        # create the TreeViewColumn to display the data
        #self.idColumn = gtk.TreeViewColumn('Id')
        self.nameColumn = gtk.TreeViewColumn('Name')
        #self.directoryColumn = gtk.TreeViewColumn('Directory')
        #self.parentIdColumn = gtk.TreeViewColumn('Parent Id')

        # add tvcolumn to treeview
        self.treeview.append_column(self.nameColumn)
        #self.treeview.append_column(self.idColumn)
        #self.treeview.append_column(self.directoryColumn)
        #self.treeview.append_column(self.parentIdColumn)

        # create a CellRendererText to render the data
        #self.idCell = gtk.CellRendererText()
        self.nameCell = gtk.CellRendererText()
        #self.directoryCell = gtk.CellRendererText()
        #self.parentIdCell = gtk.CellRendererText()

        # add the cell to the tvcolumn and allow it to expand
        #self.idColumn.pack_start(self.idCell, True)
        self.nameColumn.pack_start(self.nameCell, True)
        #self.directoryColumn.pack_start(self.directoryCell, True)
        #self.parentIdColumn.pack_start(self.parentIdCell, True)

        # set the cell "text" attribute to column 0 - retrieve text
        # from that column in treestore
        #self.idColumn.add_attribute(self.idCell, 'text', 0)
        self.nameColumn.add_attribute(self.nameCell, 'text', 1)
        #self.directoryColumn.add_attribute(self.directoryCell, 'text', 2)
        #self.parentIdColumn.add_attribute(self.parentIdCell, 'text', 3)

        # make it searchable
        self.treeview.set_search_column(1)
        self.treeview.set_enable_search(True)
        self.treeview.connect("row-activated",self.rowActivated, "OfTree")
        
        # Allow sorting on the column
        self.nameColumn.set_sort_column_id(0)

        # Allow drag and drop reordering of rows
        self.treeview.set_reorderable(True)
        
        self.scrolledWindow.add(self.treeview)
        self.window.show_all()

    def rowActivated(self, treeview, path, view_column, user_param1):
        iter = self.filteredTreestore.get_iter(path)
        projectId = self.filteredTreestore.get_value(iter, 0)
        projectName = self.filteredTreestore.get_value(iter, 1)
        self.addProjectToCoot(projectId, projectName)
    
    def addProjectToCoot(self, projectId, projectName, refmacOrLigandPipeline="LigandPipeline"):
        if refmacOrLigandPipeline=="LigandPipeline":
            archiveName = self.djangoSession.lastLigandPipelineExportForProjectId(projectId)
        elif refmacOrLigandPipeline=="REFMAC":
            archiveName = self.djangoSession.lastRefmacExportForProjectId(projectId)
        currentDirectory = os.path.realpath(os.getcwd())
        import zipfile
        try:
            zip_ref = zipfile.ZipFile(archiveName, 'r')
            zip_ref.extractall(currentDirectory)
            filePaths = [filePath for filePath in zip_ref.namelist()]
            print('filePaths',filePaths)
            pdbFiles = [filePath for filePath in zip_ref.namelist() if filePath.endswith(".pdb")]
            hkloutFiles = [filePath for filePath in zip_ref.namelist() if filePath.endswith("hklout.mtz")]
            finalMTZFiles = [filePath for filePath in zip_ref.namelist() if filePath.endswith("final.mtz")]
            filenames = [os.path.basename(filePath) for filePath in zip_ref.namelist()]
            dirnames = [os.path.dirname(filePath) for filePath in zip_ref.namelist()]
            zip_ref.close()
            for pdbFile in pdbFiles:
                newMol = handle_read_draw_molecule_with_recentre(pdbFile, 0)
                set_molecule_name(newMol, projectName)
            for hkloutFile in hkloutFiles:
                print(hkloutFile)
                iMap1, iMap2 = auto_read_make_and_draw_maps(hkloutFile)
                print(iMap1, iMap2)
                set_molecule_name(iMap1, "2FoFc_"+projectName)
                set_molecule_name(iMap2, "FoFc_"+projectName)
            for finalMTZFile in finalMTZFiles:
                iMap1, iMap2 = auto_read_make_and_draw_maps(finalMTZFile)
                set_molecule_name(iMap1, "2FoFc_"+projectName)
                set_molecule_name(iMap2, "FoFc_"+projectName)
        except zipfile.BadZipfile:
            message = gtk.MessageDialog(parent=None,
                                        flags=0,
                                        type=gtk.MESSAGE_INFO,
                                        buttons=gtk.BUTTONS_NONE,
                                        message_format=None)
            message.set_markup("Bad zip file for project "+ projectName)
            message.show()
        except:
            message = gtk.MessageDialog(parent=None,
                                        flags=0,
                                        type=gtk.MESSAGE_INFO,
                                        buttons=gtk.BUTTONS_NONE,
                                        message_format=None)
            message.set_markup("Failed for some other reason for project "+ projectName)
            message.show()

    def visible_func(self, model, iter, user_data):
        textToFilterWith = self.filterEntry.get_text()
        return (textToFilterWith.strip()=="") or (textToFilterWith.upper() in self.treestore.get_value(iter,1).upper())

