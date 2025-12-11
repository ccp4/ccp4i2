import os
import glob
from ccp4i2.core import CCP4Modules
from ccp4i2.core import CCP4File
from ccp4i2.core.CCP4QtObject import CObject
from ccp4i2.core.CCP4ErrorHandling import *


class CCustomManager(CObject):
    ERROR_CODES = {101 : { 'description' : 'Error attempting to create directory to save customisation'},
                   104 : { 'description' : 'Error attempting to overwrite directory to save customisation'},
                   105 : { 'description' : 'Error attempting to open compressed file for write'},
                   106 : { 'description' : 'Error attempting to save customisation to compressed file'},
                   107 : { 'description' : 'Error opening compressed file'},
                   108 : { 'description' : 'Error extracting from compressed file'},
                   109 : { 'description' : 'Compressed file does not have expected content'},
                   110 : { 'description' : 'Customisation in compressed file has same name as existing customisation'},
                   111 : { 'description' : 'Error attemting to remove a customisation directory'},
                   112 : { 'description' : 'Error attemting to clone a customisation directory'},
                   113 : { 'description' : 'Error opening definition file in newly cloned customisation directory'},
                   114 : { 'description' : 'Error updating definition file in newly cloned customisation directory'}}

    insts = None

    def __init__(self, parent=None, mode=None):
        if parent is None:
            parent = CCP4Modules.QTAPPLICATION()
        CObject.__init__(self,parent)
        if mode is not None:
            self.mode = str(mode)
        else:
            self.mode = None

    def getList(self):
        from ccp4i2.core import CCP4Utils
        dirList = glob.glob(os.path.join(CCP4Utils.getDotDirectory(), 'custom', self.mode + 's','*'))
        #print 'CCustomManager.getList',self.mode,os.path.join(CCP4Utils.getDotDirectory(),'custom',self.mode+'s','*'),dirList
        titleList = []
        for dr in dirList:
            titleList.append(os.path.split(dr)[1])
        return titleList

    def getTitleList(self):
        # Return a list of (name,title)
        nameList = self.getList()
        titleList = []
        for name in nameList:
            titleList.append((name,self.getTitle(name)))
        return titleList

    def getTitle(self, name):
        fileName = self.getCustomFile(name)
        #print 'CCustomManager.getTitle',name,fileName
        if fileName is None:
            return name
        if not os.path.exists(fileName):
            return name
        header = CCP4File.CI2XmlHeader(parent=self)
        header.loadFromXml(fileName)
        if header.pluginTitle.isSet():
            return header.pluginTitle.__str__()
        else:
            return name

    def getDirectory(self, name=None):
        from ccp4i2.core import CCP4Utils
        if name is None:
            return os.path.join(CCP4Utils.getDotDirectory(), 'custom', self.mode + 's')
        else:
            return os.path.join(CCP4Utils.getDotDirectory(), 'custom', self.mode + 's',str(name))

    def createDirectory(self, name, overwrite=False):
        newDir = self.getDirectory(name=name)
        if overwrite and os.path.exists(newDir):
            try:
                import shutil
                shutil.rmtree(newDir)
            except:
                raise CException(self.__class__, 104, newDir)
        try:
            os.mkdir(newDir)
        except:
            raise CException(self.__class__, 101, newDir)
        return newDir

    def getCustomFile(self, name, mustExist=True):
        fileName = os.path.join(self.getDirectory(name=name), self.mode + '.xml')
        if (not mustExist) or os.path.exists(fileName):
            return fileName
        else:
            return None

    
    def openManagerDialog(self):
        self.openCreateDialog()

    def delete(self,name):
        dr = self.getDirectory(name)
        import shutil
        try:
            shutil.rmtree(dr)
        except:
            return 1
        else:
            return 0

    def export(self, name, fileName):
        try:
            import tarfile
            tf = tarfile.open(fileName, mode='w:gz')
        except:
            err = CErrorReport(self.__class__, 105, fileName)
            return err
        flowDir = self.getDirectory(name)
        try:
            tf.add(flowDir, arcname=name)
        except:
            err = CErrorReport(self.__class__, 106, name)
            return err
        tf.close()
        return CErrorReport()

    def clone(self, original, new, title=None):
        import shutil
        diry = self.getDirectory(original)
        #print 'CCustomisationGui.clone',original,new,title,diry
        try:
            newDir = os.path.join(os.path.split(diry)[0],new)
            shutil.copytree(diry, newDir)
        except:
            raise CException(self.__class__, 112, newDir)
        fileName = self.getCustomFile(new)
        if not os.path.exists(fileName):
            raise CException(self.__class__, 113, fileName)
        try:
            fileObj = CCP4File.CI2XmlDataFile(fileName)
            fileObj.loadFile()
            fileObj.header.pluginName.set(new)
            if title is not None:
                fileObj.header.pluginTitle.set(title)
            fileObj.saveFile(fileObj.getBodyEtree())
        except:
            raise CException(self.__class__, 114, fileName)

    def testImport(self, fileName):
        import tempfile
        tmpDir = tempfile.mkdtemp()
        self.uncompress(fileName, tmpDir)
        globList = glob.glob(os.path.join(tmpDir,'*'))
        if len(globList) < 1:
            raise CException(self.__class__, 109, 'Extracting from ' + fileName + ' to ' + tmpDir)
        customXml = os.path.join(globList[0], self.mode + '.xml')
        if not os.path.exists(customXml):
            raise CException(self.__class__, 109, 'Extracting from ' + fileName + ' to ' + tmpDir)
        name = os.path.split(globList[0])[1]
        customList = self.getList()
        if name in customList:
            raise CException(self.__class__, 110, name)
        return name

    def import_(self, fileName, name, overwrite=False, rename=None):
        if overwrite:
            import shutil
            try:
                shutil.rmtree(self.getDirectory(name))
            except:
                raise CException(self.__class__, 111)

    def uncompress(self, fileName, targetDir, rename=None):
        import tarfile
        try:
            tf = tarfile.open(fileName, mode='r:gz')
        except:
            raise CException(self.__class__, 107, fileName)
        if rename is not None:
            for tarinfo in tf.members:
                #print 'uncompress tarinfo',tarinfo.name
                tarinfo.name = rename
        try:
            tf.extractall(path=targetDir)
        except:
            raise CException(self.__class__, 108,'Extracting from ' + fileName + ' to ' + targetDir)
        tf.close()
        #print 'uncompress done', os.path.samefile(targetDir,self.getDirectory())
