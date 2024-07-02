from __future__ import print_function


import os
import sys

from PySide6 import QtCore
from core.CCP4ErrorHandling import *
from core.CCP4Modules import PROJECTSMANAGER

COMPRESSED_SUFFIX = 'ccp4_project.zip'
COMPRESSION_MODE = 'zip'
ALLOWZIP64 = True

#  See discussion on handling busted zip files..
#  http://stackoverflow.com/questions/20890950/python-extracting-files-from-a-large-6gb-zip-file

class ExportProjectThread(QtCore.QThread):

    startSavingJobData = QtCore.Signal(int)
    savingJobData = QtCore.Signal(tuple)

    ERROR_CODES = {170 : {'description' : 'Export unfinished or failed'},
                   171 : {'description' : 'Error searching database for project jobs'},
                   172 : {'description' : 'Error saving to project tar compressed file'},
                   173 : {'description' : 'Error deleting temporary database xml file'},
                   174 : {'description' : 'Error creating temporary database xml file'},
                   175 : {'description' : 'Error creating project compressed tar file'},
                   180 : {'description' : 'Error creating project compressed zip file'},
                   181 : {'description' : 'Error saving job to project zip compressed file'},
                   182 : {'description' : 'Error saving input file to project zip compressed file'},
                   183 : {'description' : 'Error saving directory to project zip compressed file'},
                   184 : {'description' : 'Error saving database to project zip compressed file'},
                   185 : {'description' : 'Error closing project zip compressed file'}}

    def __init__(self,parent=None,projectDir=None,dbxml=None,target=None,jobList=[],inputFilesList=[],directoriesList=['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES','CCP4_TEST_SYSTEM'],extraJobList=[]):
        QtCore.QThread.__init__(self,parent)    
        self.projectDir = projectDir
        self.dbxml=dbxml
        self.target=target
        self.jobList = jobList
        self.extraJobList = extraJobList
        self.inputFilesList = inputFilesList
        self.directoriesList = directoriesList
        self.errorReport = CErrorReport(self.__class__,170)

    def run(self):    
        self.errorReport = self.compressProject()
        if len(self.errorReport) > 0:
            print('ExportProjectThread errors',self.errorReport.report())
        return

    def compressProject(self):
        if COMPRESSION_MODE == 'zip':
            return self.zipCompressProject()
        else:
            return self.tarCompressProject()

    def tarCompressProject(self):
        try:
            import tarfile
            tf = tarfile.open(self.target,mode='w:gz')
        except:
            err = CErrorReport(self.__class__, 175, fileName)  # KJS : fileName ? Looks like an issue here.
            return err
        self.startSavingJobData.emit(len(self.jobList))
        savedInputRelPaths0 = []
        savedInputRelPaths1 = []
        try:
            done=0
            saveDir = ''
            for job in self.jobList:
                savedInputRelPaths0.append(os.path.join('CCP4_JOBS','job_'+job))
                saveDir = os.path.join(self.projectDir,'CCP4_JOBS','job_'+job)
                print('Saving job number:',job,saveDir);sys.stdout.flush()
                PROJECTSMANAGER().cleanupJob(jobDirectory=saveDir)
                tf.add(saveDir,arcname='CCP4_JOBS/job_'+job)
                done+=1      
                self.savingJobData.emit((job,done))
            for relPath,baseName in self.inputFilesList:
                if relPath not in savedInputRelPaths0:
                    # Its not been saved from the jobList
                    if not relPath in savedInputRelPaths1:
                        # Make a directory in the archive - necessary to jog the import to
                        # know it has got that job in the archive
                        saveDir = os.path.join(self.projectDir,relPath)
                        tf.add(saveDir,arcname=relPath,recursive=False)
                        savedInputRelPaths1.append(relPath)
                    saveDir = os.path.join(self.projectDir,relPath,baseName)
                    #print 'ExportProjectThread.compressProject',relPath,baseName,saveDir
                    tf.add(saveDir,arcname=os.path.join(relPath,baseName))
            for dirName in self.directoriesList:
                saveDir = os.path.join(self.projectDir,dirName)
                if os.path.exists(saveDir):
                    self.savingJobData.emit(('IMPORT',0))
                    tf.add(saveDir,arcname=dirName)
                    self.savingJobData.emit(('IMPORT',1))
            saveDir = self.dbxml
            self.savingJobData.emit(('DATABASE',0))
            tf.add(saveDir,arcname='DATABASE.db.xml')
            self.savingJobData.emit(('DATABASE',1))
            tf.close()
        except:
            #print 'ERROR saving to tarfile',saveDir
            err = CErrorReport(self.__class__, 172, saveDir)
            return err
        try:
            os.remove(self.dbxml)
        except:
            err = CErrorReport(self.__class__, 173)
            return err
        return CErrorReport()
  
    def zipCompressProject(self):
        from core import CCP4Utils
        # zipfile does not recurse over directory so use CCP4Utils.zipDirectory()
        try:
            import zipfile
            zip = zipfile.ZipFile(self.target, mode='w', allowZip64=ALLOWZIP64)
        except:
            return CErrorReport(self.__class__, 180, self.target)
        self.startSavingJobData.emit(len(self.jobList))
        savedInputRelPaths0 = []
        savedInputRelPaths1 = []
        done=0
        saveDir = ''
        err = CErrorReport()
        for job in self.jobList:
            try:
                savedInputRelPaths0.append(os.path.join('CCP4_JOBS', 'job_' + job))
                saveDir = os.path.join(self.projectDir, 'CCP4_JOBS', 'job_' + job)
                print('Saving job number:', job, saveDir); sys.stdout.flush()
                PROJECTSMANAGER().cleanupJob(jobDirectory=saveDir)
                CCP4Utils.zipDirectory(zip, saveDir, rootRelPath=self.projectDir)
                done += 1
                self.savingJobData.emit((job, done))
            except Exception as e:
                print(e)
                err.append(self.__class__, 181, 'Job number: ' + str(job) + '\n'+str(e))
                if len(err) > 3:
                    return err
        for relPath, baseName in self.inputFilesList:
            try:
                if relPath not in savedInputRelPaths0:
                    # Its not been saved from the jobList
                    if not relPath in savedInputRelPaths1:
                        # Make a directory in the archive - necessary to jog the import to
                        # know it has got that job in the archive
                        saveDir = os.path.join(self.projectDir,relPath)
                        zip.write(saveDir,arcname=relPath)
                        savedInputRelPaths1.append(relPath)
                    saveDir = os.path.join(self.projectDir,relPath,baseName)
                    #print 'ExportProjectThread.compressProject',relPath,baseName,saveDir
                    CCP4Utils.zipDirectory(zip,saveDir,rootRelPath=self.projectDir)
            except:
                err.append(self.__class__,182,'File: '+str(relPath)+' '+str(baseName))
                if len(err) > 3:
                    return err
        # deal with problem case of imported file being in the IMPORTED_FILES directory
        # but the project import mechanism expecting a job directory
        for job in self.extraJobList:
            relPath = os.path.join('CCP4_JOBS','job_'+job)
            if not relPath in savedInputRelPaths0 and not relPath in savedInputRelPaths1:
                savedInputRelPaths1.append(relPath)
                saveDir = os.path.join(self.projectDir,relPath)
                zip.write(saveDir,arcname=relPath)
        for dirName in self.directoriesList:
            try:
                saveDir = os.path.join(self.projectDir,dirName)
                if os.path.exists(saveDir):
                    self.savingJobData.emit(('IMPORT',0))
                    CCP4Utils.zipDirectory(zip,saveDir,rootRelPath=self.projectDir)
                    self.savingJobData.emit(('IMPORT',1))
            except:
                err.append(self.__class__,183,'Directory: '+str(dirName))
                if len(err) > 3:
                    return err
        saveDir = self.dbxml
        self.savingJobData.emit(('DATABASE',0))
        try:
            zip.write(saveDir, arcname='DATABASE.db.xml')
            self.savingJobData.emit(('DATABASE',1))
        except:
            err.append(self.__class__, 184)
            return err
        try:
            zip.close()
        except:
            return CErrorReport(self.__class__, 185, saveDir)
        try:
            os.remove(self.dbxml)
        except:
            err.append(self.__class__, 173, str(self.dbxml))
        return err


#class ImportProjectThread(QtCore.QThread):
class ImportProjectThread(QtCore.QObject):

    extractingJobData = QtCore.Signal(tuple)

    ERROR_CODES = {170 : {'description' : 'Importing files unfinished or failed'},
                   176 : {'description' : 'Error opening project tar compressed file'},
                   177 : {'description' : 'Error reading job files from project tar compressed file'},
                   178 : {'severity' : SEVERITY_WARNING, 'description' : 'Error reading from project tar compressed file'},
                   190 : {'description' : 'Error opening project zip compressed file'},
                   191 : {'description' : 'Error reading job files from project zip compressed file'},
                   192 : {'severity' : SEVERITY_WARNING, 'description' : 'Error reading from project zip compressed file'},
                   191 : {'description' : 'Error reading compressed file - file extension unrecognised (expects tar.gz or zip)'}}

    def __init__(self,parent=None, projectDir=None, compressedFile=None, dbImport=None, diagnostic=False):
        #QtCore.QThread.__init__(self,parent)
        QtCore.QObject.__init__(self, parent)
        self.projectDir = projectDir
        self.compressedFile=compressedFile
        self.dbImport = dbImport
        self.diagnostic = diagnostic
        self.errReport = CErrorReport()

    def run(self):
        if self.diagnostic:
            print('Starting to extract files from: ', self.compressedFile)
            print('to', self.projectDir)
        self.extractJobs(self.compressedFile, self.projectDir, dbImport=self.dbImport)
        for mode in ['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES','CCP4_TEST_FILES']:
            self.extractingJobData.emit((mode,0))
            self.extractProjectDir(self.compressedFile, self.projectDir, mode=mode)
            self.extractingJobData.emit((mode,1))
        # Ensure all other directories are created
        PROJECTSMANAGER().makeProjectDirectory(directory=self.projectDir)
        return self.errReport
  
    def extractJobs(self,compressedFile,targetDir,dbImport=None):
        if compressedFile.endswith('tar.gz'):
            self.extractJobsFromTar(compressedFile, targetDir=targetDir, dbImport=dbImport)
        elif compressedFile.endswith('zip'):
            self.extractJobsFromZip(compressedFile, targetDir=targetDir, dbImport=dbImport)
        else:
            self.errReport.append(self.__class__, 193, str(compressedFile))
    
    def extractJobsFromTar(self, compressedFile, targetDir, dbImport=None):
        #print 'extractAllJobs nextJobNumber',nextJobNumber,targetDir, dbImport
        # dbImport (a CDbXml) is passed if we are not importing the entire project
        # - it provides dbImport.importThisFile() function which returns ifImport flag and newJobNumber
        # Note that dbImport.importThisFile() accesses the (probably sqlite) database which
        # has the limitation that it will only run in one thread
        # Only if we have a dbImport use the tarfile.extractall() argument members to call the
        # copyThisFile() function that accesses  dbImport.importThisFile()
        import tarfile
        try:
            tf = tarfile.open(compressedFile, mode='r:gz')
        except:
            raise CException(self.__class__, 176, compressedFile)
        #tf.list()
        done = 0

    def copyThisFile(self,members, dbImport=None): # KJS : Hmmm.. think it's missing a self in there.
        lastJobNo = None                        # Also this function is absolutely riddled with errors
        done = -1
        for tarinfo in members:
            #print 'ImportProjectThread.copyThisFile tarinfo.name',tarinfo.name
            if tarinfo.name[0:9] == "CCP4_JOBS":
                # Interpret job number and fileName from tarinfo.name
                newName = None
                dirSplit = tarinfo.name.split('/')
                jobNo = dirSplit[1].split('_')[1]
                level = 2
                fileName = None
                while level<len(dirSplit):
                    dirName = dirSplit[level]
                    if dirName[0:4] == 'job_':
                        jobNo = jobNo + '.' + dirName.split('_')[1]
                        level += 1
                    else:
                        fileName = dirName
                        break
                ifImport, newJobNumber = dbImport.importThisFile( jobNumber = jobNo, fileName = fileName )
                #print 'Import.importThisFile',tarinfo.name,jobNo,ifImport,newJobNumber,level
                if ifImport:
                    done +=1
                    if newJobNumber is not None:
                        newJobSplit = newJobNumber.split('.')
                        newName = "CCP4_JOBS/job_"+newJobSplit[0]
                        for jn in newJobSplit[1:]:
                            newName = newName + '/job_'+jn
                        for item in dirSplit[level:]:
                            newName = newName + '/' + item
                        self.extractingJobData.emit((newJobNumber,done))
                    else:
                        self.extractingJobData.emit((jobNo,done))
                    if newName is not None:
                        #print 'Renaming file from',tarinfo.name,'to',newName
                        tarinfo.name = newName
                    yield tarinfo
        try:
            if dbImport is not None:
                tf.extractall(path=targetDir,members=copyThisFile(tf,dbImport))
            else:
                tf.extractall(path=targetDir)
        except:
            tf.close()
            self.errReport.append(self.__class__,177,'Extracting from '+compressedFile+' to '+targetDir)
        if dbImport is not None:
            dbImport.db.commit()
        tf.close()


    def copyThisFile(self, zinfo, dbImport=None): # KJS : Aww this is just marvellous, a borked duplicate.
        lastJobNo = None
        done = -1
        #print 'ImportProjectThread.copyThisFile zinfo.filename',zinfo.filename
        if zinfo.filename[0:9] == "CCP4_JOBS":
            # Interpret job number and fileName from zinfo.filename
            newName = None
            dirSplit = zinfo.filename.strip('/').split('/')
            jobNo = dirSplit[1].split('_')[1]
            level = 2
            fileName = None
            while level<len(dirSplit):
                dirName = dirSplit[level]
                if dirName[0:4] == 'job_':
                    jobNo = jobNo + '.' + dirName.split('_')[1]
                    level += 1
                else:
                    fileName = dirName
                    break
            ifImport,newJobNumber = dbImport.importThisFile( jobNumber = jobNo, fileName = fileName )
            #print 'Import.importThisFile',zinfo.filename,jobNo,ifImport,newJobNumber,level
            if ifImport:
                done +=1
                if newJobNumber is not None:
                    newJobSplit = newJobNumber.split('.')
                    newName = "CCP4_JOBS/job_"+newJobSplit[0]
                    for jn in newJobSplit[1:]: newName = newName + '/job_'+jn
                    for item in dirSplit[level:]: newName = newName + '/' + item
                    self.extractingJobData.emit((newJobNumber,done))
                else:
                    self.extractingJobData.emit((jobNo,done))
                      
                return  zinfo.filename,newName
            else:
                return None,None
        else:
            return zinfo.filename,None

    def extractJobsFromZip(self, compressedFile, targetDir, dbImport=None):
        #print 'extractAllJobs nextJobNumber',nextJobNumber,targetDir, dbImport
        # dbImport (a CDbXml) is passed if we are not importing the entire project
        # - it provides dbImport.importThisFile() function which returns ifImport flag and newJobNumber
        # Note that dbImport.importThisFile() accesses the (probably sqlite) database which
        # has the limitation that it will only run in one thread
        # Only if we have a dbImport use the zipfile.extractall() argeument members to call the
        # copyThisFile() function that accesses  dbImport.importThisFile()
        import zipfile
        try:
            zip = zipfile.ZipFile(compressedFile,mode='r',allowZip64=ALLOWZIP64)
        except:
            if self.diagnostic:
                print('Error reading zip file',compressedFile)
            raise CException(self.__class__,190,compressedFile)
        #print 'extractJobsFromZip'; zip.list()
        if dbImport is not None:
            tmpDir = None
            for zinfo in zip.infolist():
                #There is special code in copyThisFile if filename starts with CCP4_JOBS, which assumes plain directory will not exist in zip file. ?!
                if zinfo.filename == "CCP4_JOBS" or zinfo.filename == "CCP4_JOBS/":
                    continue
                try:
                    filename,newName = self.copyThisFile(zinfo,dbImport)
                    #print 'copyThisFile',filename,newName
                    if not newName:
                        zip.extract(zinfo,targetDir)
                    else:
                        if tmpDir is None:
                            import tempfile,shutil
                            tmpDir = tempfile.mkdtemp()
                            print('Extracting files to temp directory:',tmpDir)
                        zip.extract(zinfo,tmpDir)
                        if filename.count('CCP4_JOBS/job_223/job_1'): print('Extracting',filename,'to',newName)
                        try:
                            targetObj = os.path.join(targetDir,newName)
                            shutil.move(os.path.join(tmpDir,filename),targetObj)
                            print('Creating',targetObj,os.path.isdir(targetObj))
                        except:
                            pass
                except:
                    print('Error extracting',zinfo.filename)
        else:
            zip.extractall(path=targetDir)
        #except:
        #  zip.close()
        #  self.errReport.append(self.__class__,191,'Extracting from '+compressedFile+' to '+targetDir)
        if dbImport is not None:
            dbImport.db.commit()
        zip.close()

    def extractProjectDir(self,compressedFile,targetDir,mode='CCP4_IMPORTED_FILES',selectedFiles=None):
        if compressedFile.endswith('tar.gz'):
            self.extractProjectDirFromTar(compressedFile,targetDir,mode=mode,selectedFiles=selectedFiles)
        elif compressedFile.endswith('zip'):
            self.extractProjectDirFromZip(compressedFile,targetDir,mode=mode,selectedFiles=selectedFiles)
        else:
            self.errReport.append(self.__class__,193,str(compressedFile))

    def extractProjectDirFromTar(self,compressedFile,targetDir,mode='CCP4_IMPORTED_FILES',selectedFiles=None):
        import tarfile
        try:
            tf = tarfile.open(compressedFile,mode='r:gz')
        except:
            raise CException(self.__class__,176,compressedFile)
        #---------------------------------------
        def isSelected(members):
            for tarinfo in members:
                if tarinfo.name[0:len(mode)] == mode:
                    if selectedFiles is None:
                        yield tarinfo
                    else:
                        #try:
                        fileName = tarinfo.name.split('/')[1]
                        #except:
                        #   yield tarinfo
                        #else:
                        if fileName in selectedFiles:
                            yield tarinfo
        #---------------------------------------
        try:
            tf.extractall(path=targetDir,members=isSelected(tf))
        except:
            self.errorReport.append(self.__class__,178,'Extracting '+str(mode)+' from '+compressedFile+' to '+targetDir)
        tf.close()
    
    def extractProjectDirFromZip(self,compressedFile,targetDir,mode='CCP4_IMPORTED_FILES',selectedFiles=None):
        import zipfile
        try:
            zip = zipfile.ZipFile(compressedFile,mode='r',allowZip64=ALLOWZIP64)
        except:
            raise CException(self.__class__,190,compressedFile)
        #print 'extractProjectDirFromZip selectedFiles',selectedFiles
        #-----------------------------------------------------------------
        def isSelected(nameList):
            #print 'extractProjectDirFromZip.isSelected nameList',nameList
            for name in nameList:
                if name[0:len(mode)] == mode:
                    if selectedFiles is None:
                        yield name
                    else:
                        #try:
                        fileName = name.split('/')[1]
                        #except:
                        #   yield tarinfo
                        #else:
                        if fileName in selectedFiles:
                            yield name
        #----------------------------------------------------------------
        try:
            zip.extractall(path=targetDir,members=isSelected(zip.namelist()))
        except:
            self.errReport.append(self.__class__,192,'Extracting '+str(mode)+' from '+compressedFile+' to '+targetDir)   
        zip.close()
