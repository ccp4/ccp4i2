from __future__ import print_function


import sys, os
if sys.version_info >= (3,0):
    import urllib.request, urllib.error, urllib.parse
else:
    import urllib
    import urllib2
import traceback

class ccp4i2CootInterface():
    FileTypes = [('chemical/x-pdb',None,'Import coordinates from project'),
     ('application/CCP4-mtz-map', 1, 'Import 2Fo-Fc maps from project'),
     ('application/CCP4-mtz-map', 2, 'Import Fo-Fc maps from project'),
     ('application/CCP4-mtz-map', 3, 'Import Anomolous difference maps from project'),
     ('application/refmac-dictionary', None, 'Import CIF dictionary'),
     ]

    def __init__(self, dropDir=None):
        self.dropDir = dropDir
        print('** Drop directory for output files is ', self.dropDir)
        #Cache relevant environment variables
        self.ccp4i2HTTPPort = os.environ.get('CCP4I2_HTTP_PORT',None)
        print('** CCP4I2_HTTP_PORT is ', self.ccp4i2HTTPPort)
        if self.ccp4i2HTTPPort is not None:
          self.serverRoot = "http://127.0.0.1:%s/database"%self.ccp4i2HTTPPort
        else:
          self.serverRoot = None

        self.ccp4i2ProjectName = os.environ.get('CCP4I2_PROJECT_NAME',None)
        print('** CCP4I2_PROJECT_NAME is ')
        print(self.ccp4i2ProjectName)

        self.ccp4i2ProjectID = os.environ.get('CCP4I2_PROJECT_ID',None)
        print('** CCP4I2_PROJECT_ID is ')
        print(self.ccp4i2ProjectID)

    def addSeparator(self, menu):
        import gtk
        sep = gtk.MenuItem()
        menu.add(sep)
        sep.show()

    def installMenus(self):
        import gtk
        menu = coot_menubar_menu("CCP4i2 extensions")
        self.ccp4i2ExtensionsMenuItem = menu.get_attach_widget()
        self.ccp4i2ExtensionsMenuItem.connect('activate', self.populateMenus)

        if self.serverRoot is not None:
          #Add menu button to import coordinates
          #First coordinates

          self.projectDirectory = self.performDatabaseLookup("getProjectDirectory",projectId=self.ccp4i2ProjectID)
          for fileType, subType, menuText in self.FileTypes:
              menuitem = gtk.MenuItem(menuText)
              menu.append(menuitem)
              menuitem.show()
        self.addSeparator(menu)
        # Add broader project interrogating item
        self.otherProjectsMenuItem = gtk.MenuItem('Other projects')
        menu.append(self.otherProjectsMenuItem)
        self.otherProjectsMenuItem.show()
        self.addSeparator(menu)
        # Add simple save options
        add_simple_coot_menu_menuitem(menu, "Save mol to CCP4i2...", lambda func : self.selectAndSave() )
        add_simple_coot_menu_menuitem(menu, "Save to CCP4i2", lambda func : self.saveToDatabase(0) )
        #add_simple_coot_menu_menuitem(menu, "Launch Django interaction", lambda func : self.launchDjangoWidget() )

        fileMenu = coot_menubar_menu("File")
        self.addSeparator(fileMenu)
        add_simple_coot_menu_menuitem(fileMenu, "Save mol to CCP4i2...", lambda func : self.selectAndSave() )
        add_simple_coot_menu_menuitem(fileMenu, "Save to CCP4i2", lambda func : self.saveToDatabase(0) )

    # def launchDjangoWidget(self):
    #     from ccp4i2.wrappers.SyncToDjango.script import testTree
    #     self.djangoWidget = testTree.DjangoInteractionWidget()

    def populateMenus(self, menuitem):
        import gtk
        self.files = []
        self.populateSubmenus(menuitem, [self.ccp4i2ProjectID])

        projectList  = self.performDatabaseLookup('listProjects')
        #print projectList
        #sys.stdout.flush()
        projectSubMenu = gtk.Menu()
        projectSubMenu.show()
        for projectInfo in projectList:
            projectMenuItem = gtk.MenuItem(projectInfo[1])
            projectMenuItem.connect('activate', self.populateSubmenus, projectInfo)
            projectMenuItem.show()
            fileTypesMenu = gtk.Menu()
            fileTypesMenu.show()
            projectMenuItem.set_submenu(fileTypesMenu)
            for fileType, subType, menuText in self.FileTypes:
                fileTypeMenuItem = gtk.MenuItem(menuText)
                fileTypeMenuItem.show()
                fileTypesMenu.append(fileTypeMenuItem)
            projectSubMenu.append(projectMenuItem)
            projectMenuItem.show()
        self.otherProjectsMenuItem.set_submenu(projectSubMenu)
        self.otherProjectsMenuItem.show()

    def addFilesOfProject_ofType_subType_toMenu(self, projectId, fileType, subType, menu):
        import gtk
        self.projectDirectory = self.performDatabaseLookup("getProjectDirectory",projectId=projectId)
        fileTypeAndSubType = fileType
        if subType is not None: fileTypeAndSubType+= ('_'+str(subType))
        getFilesArgs  = {'projectId':projectId, 'fileType':fileType, 'topLevelOnly':True}
        if subType is not None: getFilesArgs['subType'] = subType
        newFiles2 = self.getProjectFiles(**getFilesArgs)
        newFiles1 = sorted(newFiles2, key=lambda fileJob: int(str(fileJob[6]).split('.')[0]), reverse=True)

        newFiles = []
        for newFile1 in newFiles1:
            file = {'fileid':newFile1[1], 'jobnumber':newFile1[6], 'filename':newFile1[4],'annotation':newFile1[5]}
            file['fullpath'] = os.path.join(self.projectDirectory,'CCP4_JOBS',
                                            'job_'+str(file['jobnumber']),file['filename'])
            textToUse = 'Default label'
            try: textToUse = file['jobnumber']
            except: pass
            try: textToUse += file['annotation']
            except: pass
            menuItem = gtk.MenuItem(textToUse)
            menuItem.connect('activate',self.loadFileClicked, file['fileid'], fileTypeAndSubType)
            menuItem.show()
            menu.append(menuItem)
            newFiles.append(file)
        self.files += newFiles

    def loadFileClicked(self, menuitem, fileId, type):
        self.retrieveAndLoadFileByID(fileId=fileId, type=type)

    def populateSubmenus(self, menuitem, projectInfo):
        import gtk
        sys.stdout.flush()
        for child in menuitem.get_submenu().get_children():
            for fileType, subType, menuText in self.FileTypes:
                if  child.get_label().__str__() == menuText:
                    subMenu = gtk.Menu()
                    subMenu.show()
                    self.addFilesOfProject_ofType_subType_toMenu(projectInfo[0], fileType, subType, subMenu)
                    child.set_submenu(subMenu)
                    continue

        sys.stdout.flush()

    def retrieveAndLoadFileByID(self, fileId=None, type='PDB'):
        if fileId is None: return
        outputFilename = 'ProjectFile_'+fileId
        if type == 'chemical/x-pdb': outputFilename += '.pdb'
        elif type == 'application/CCP4-mtz-map_1': outputFilename += '_2FoFc.mtz'
        elif type == 'application/CCP4-mtz-map_2': outputFilename += '_FoFc.mtz'
        destination = os.path.normpath(os.path.join(self.dropDir, outputFilename))
        downloadedFilePath = self.retrieveFileByID(fileId=fileId, destination=outputFilename)
        if downloadedFilePath is not None:
            fileInfo = self.getFileInfo(fileId=fileId,mode='all')
            if fileInfo['annotation'] is not None:
                annotation = fileInfo['annotation'].replace('/','_')
            else:
                annotation = ''
            jobnumber = self.getFileInfo(fileId=fileId,mode='jobnumber')['jobnumber']
            projectname = self.getFileInfo(fileId=fileId,mode='projectname')['projectname']
            annotationRoot = ""
            if projectname != self.ccp4i2ProjectName: annotationRoot = projectname+":"
            if type == 'chemical/x-pdb':
                newMol = read_pdb(downloadedFilePath)
                set_molecule_name(newMol,annotationRoot + 'Job '+str(jobnumber)+': '+annotation)
            elif type == 'application/CCP4-mtz-map_1':
                newMap=make_and_draw_map(downloadedFilePath, 'F', 'PHI', 'PHI', 0, 0)
                set_molecule_name(newMap,annotationRoot + 'Job '+str(jobnumber)+': '+annotation)
            elif type == 'application/CCP4-mtz-map_2':
                newMap=make_and_draw_map(downloadedFilePath, 'F', 'PHI', 'PHI', 0, 1)
                set_molecule_name(newMap,annotationRoot + 'Job '+str(jobnumber)+': '+annotation)
            elif type == 'application/refmac-dictionary':
                read_cif_dictionary(downloadedFilePath)

    def retrieveFileByID(self, fileId=None, destination=None):
        if fileId is None: return
        if destination is None: return
        url = self.serverRoot + "?File?fileId="+fileId
        try:
            with open(destination,"w") as f:
                #proxy = urllib2.ProxyHandler({})
                #opener = urllib2.build_opener(proxy)
                #urllib2.install_opener(opener)
                if sys.version_info >= (3,0):
                    f.write(urllib.request.urlopen(url,proxies={}).read())
                else:
                    f.write(urllib.urlopen(url,proxies={}).read())
        except IOError:
            print('Unable to retrieve '+url)
            return None
        print('Ostensibly recovered url '+url+' to '+destination)
        if os.path.isfile(destination):
            return destination
        else:
            print('But not found it at '+destination+':-(')
            return None

    def getJobsWithOutputFiles(self, **kwargs):
        #Premissible arguments: projectId=None, fileTypeId=None,projectName=None,fileType=None,fileTypeClass=None,
        #subType=None,contentFlag=None,topLevelOnly=True,importFiles=False):
        return self.performDatabaseLookup("getJobsWithOutputFiles",**kwargs)

    def getProjectFiles(self, **kwargs):
        #Premissible arguments: projectId=None, fileTypeId=None,projectName=None,fileType=None,fileTypeClass=None,
        #subType=None,contentFlag=None,topLevelOnly=True,importFiles=False):
        return self.performDatabaseLookup("getProjectFiles",**kwargs)

    def getFullPath(self, **kwargs):
        return self.performDatabaseLookup("getFullPath",**kwargs)

    def getFileInfo(self, **kwargs):
        print(type(kwargs))
        return self.performDatabaseLookup("getFileInfo",**kwargs)

    def performDatabaseLookup(self, command, **kwargs):
        url = self.serverRoot + "?" + command
        for key in kwargs: url+="?"+str(key)+"="+str(kwargs[key])
        # fetch the url
        if sys.version_info >= (3,0):
            json = urllib.request.urlopen(url,proxies={}).read()
        else:
            json = urllib.urlopen(url,proxies={}).read()
        # convert to a native python object
        (true,false,null) = (True,False,None)
        fileInfo = eval(json)
        return fileInfo

    def saveToDatabase(self, molNo):
        import os,glob
        outList = glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.pdb')))
        outList += glob.glob(os.path.normpath(os.path.join(self.dropDir,'output*.cif')))
        #print 'saveToDatabase outList:'
        maxIndx = 0
        for f in outList:
            fpath,fname = os.path.split(f)
            maxIndx =  max(maxIndx,int(fname[6:-4]))
        #print 'saveToDatabase',dropDir,maxIndx
        molName = molecule_name(molNo)
        if molName.endswith(" (CIF)"):
            save_coordinates(molNo,os.path.normpath(os.path.join (self.dropDir,'output'+str(maxIndx+1)+'.cif') ))
        else:
            save_coordinates(molNo,os.path.normpath(os.path.join (self.dropDir,'output'+str(maxIndx+1)+'.pdb') ))
        save_state_file_py(os.path.normpath(os.path.join (self.dropDir,'state'+str(maxIndx+1)+'.py') ) )
        return

    def selectAndSave(self):
        molecule_chooser_gui ("Molecule to capture: ", lambda imol : self.saveToDatabase(imol))
        return

    def addInterestingBitsMenu(self, title='Interesting bits', interestingBits = []):
        import gtk
        menu = coot_menubar_menu("CCP4i2 extensions")
        self.addSeparator(menu)
        import gtk
        #Create submenu for the "interesting bts"
        menuitem = gtk.MenuItem(title)
        menu.append(menuitem)
        submenu = gtk.Menu()
        menuitem.set_submenu(submenu)
        menuitem.show()
        for interestingBit in interestingBits:
            textToUse = 'Chain %s %s-%s'%( interestingBit['chain'], interestingBit['firstResidue'], interestingBit['lastResidue'] )
            menuItem = gtk.MenuItem ( textToUse )
            submenu.append ( menuItem )
            menuItem.connect ( 'activate', self.interestingBitItemClicked, interestingBit['chain'], interestingBit['firstResidue'] )
            menuItem.show ( )

    def interestingBitItemClicked(self, menuitem, chain, residue):
        self.centreOn(chain=chain, residue=residue, iMol = 0)

    def centreOn(self, chain='A', residue='1', iMol=0):
        import re
        literals = re.findall ( r'\D' , str(residue) )
        number   = int ( re.findall ( r'\d+', str(residue) )[0] )

        if len ( literals ) > 0 :
            literal = literals[0]
        else :
            literal = ""

        centreList = residue_centre_py(int(iMol),str(chain),int(number),str(literal))
        set_rotation_centre(centreList[0],centreList[1],centreList[2])

    def patchMoleculeName(self, molNo, fileName, isCif=False):
        if not hasattr(self,'files'): self.populateMenus(self.ccp4i2ExtensionsMenuItem)
        label = 'Molecule '+str(molNo)
        for file in self.files:
            if 'fullpath' in file and file['fullpath'] == fileName:
                label = ('Job '+str(file['jobnumber'])+': '+file['annotation'].replace('/','_')).replace("(","").replace(")","")
                if isCif:
                    label += " (CIF)"
        set_molecule_name(molNo, label)

    def convert ( self, item ):
        residue_names={ 'UNK', 'ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU' }
        ls = model_molecule_number_list()
        imol = ls[-1]
        ins_code = '' # do we need to fill this really?
        alt_conf = '' # or this???
        chain_id = item['chain']
        res_no   = int ( item['id'] )
        name_aa  =  item['name']
        target = ''

        if name_aa in residue_names :
            target = ' CA '
        elif name_aa == "HOH" :
            target = ' O '
        else :
            target = ' C1 '

        atom=[ imol, chain_id, res_no, ins_code, target, alt_conf ]
        label = name_aa + item['id'] + " in chain " + chain_id
        r = atom
        r.insert(0, label)
        diagnostic = item['diagnostic']
        potential_fix = item['potential_fix']

        fix_func = [ refine_zone, imol, chain_id, res_no, res_no, ""]

        r.append(fix_func)
        r.append("Refine!")
        r.append("Refine!")

        return r


    def add_bits_menu ( self, title, interesting_bits ):
        nls = [ self.convert(item) for item in interesting_bits ]
        interesting_things_gui ( title, nls )


    def add_consolidated_menu ( self, title, interesting_bits ):
        menu = coot_menubar_menu ( "EDSTATS" )
        add_simple_coot_menu_menuitem ( menu, title, lambda func: self.add_bits_menu(title, interesting_bits) )
