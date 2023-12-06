from __future__ import print_function


"""
     CCP4Utils.py: CCP4 GUI Project
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009-2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Jan 2010 - Copied from ccp4mg python/ui/utils.py and converted to Qt
"""

import os
import sys
import re
import types
import pickle
import time
import glob
import copy
import shutil
import tarfile
from lxml import etree
from core.CCP4Config import DEVELOPER
from core.CCP4ErrorHandling import *

def writeXML(f,t):
    if sys.version_info > (3,0):
        f.write(t.decode("utf-8"))
    else:
        try:
            f.write(t)
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')
            import traceback
            traceback.print_tb(exc_tb)
            sys.stdout.flush()
            

class CUtils:

    ERROR_CODES = {101 : {'description' :'Input is not etree element for output to' },
                   102 : {'description' :'Failed creating XML text' },
                   103 : {'description' :'Failed writing file'},
                   104 : {'description' :'Failed parsing file'},
                   105 : {'description' :'Error finding root in'},
                   106 : {'description' :'Failed opening file to write'},
                   107 : {'description' :'Error writing file'},
                   108 : {'description' :'Failed opening file to read'},
                   109 : {'description' :'Error reading file'},
                   110 : {'description' :'Error opening file'},
                   111 : {'description' :'Error reading limited number of lines from file'}}
    pass

##Convert input string to int
#@param value Input string
#@param default Optional default value

def safeInt(value, default=None):
    if isinstance(value,str):
        value = value.strip()
    if not value:
        return default
    try:
        i = int(value)
    except:
        i = default
    return i

##Convert input string to float
#@param value Input string
#@param default Optional default value

def safeFloat(value, default=None):
    if isinstance(value,str):
        value = value.strip()
    if not value:
        return default
    try:
        i = float(value)
    except:
        i = default
    return i

##Convert input string to Boolean
#@param value Input string
#@param default Optional default value

def safeBoolean(value, default=None):
    if isinstance(value,str):
        value = value.strip()
    if not value:
        return default
    try:
        i = bool(int(value))
    except:
        i = default
    return i

def safeOneWord(value):
    # Replace spaces and any 'odd' characters with underscore
    # Resulting word should be safe in a filename
    one_word_name = re.sub('[^a-zA-Z0-9_-]', '_', value)
    return one_word_name

def newLineFormat(fileName):
    f = openFile(fileName)
    line = f.readline()
    f.close()
    if line.endswith('\r\n'):
        return 'windows'
    else:
        return 'unix'

def dos2unix(fileName, fileOut=None):
    txt = readFile(fileName)
    txtout = re.sub(r'\r\n', '\n', txt)
    if fileOut is None:
        fSplit = os.path.splitext(fileName)
        if sys.platform.startswith('darwin'):
            fileOut = fSplit[0] + '_macosx' + fSplit[1]
        else:
            fileOut = fSplit[0] + '_linux' + fSplit[1]
    saveFile(fileOut, txtout)
    return fileOut

def unix2dos(fileName, fileOut=None):
    txt = readFile(fileName)
    txtout = re.sub(r'\n', '\r\n', txt)
    if fileOut is None:
        fSplit = os.path.splitext(fileName)
        fileOut = fSplit[0] + '_windows' + fSplit[1]
    saveFile(fileOut, txtout)
    return fileOut

#-----------------------------------------------------------------
def readFile(fileName=None, limit=None):
#-----------------------------------------------------------------
    fileName = os.path.normpath(fileName)
    try:
        f = open(fileName,'r')
    except:
        raise CException(CUtils, 108, str(fileName))
    if limit is None:
        try:
            text = f.read()
            f.close()
        except:
            raise CException(CUtils, 109, str(fileName))
    else:
        #try:
        text = f.read(limit)
        f.close()
        #except:
        #  raise CException(CUtils,111,str(fileName))
    return text


#-----------------------------------------------------------------
def openFile(fileName=None, mode='r', overwrite=1):
#-----------------------------------------------------------------
    fileName = os.path.normpath(fileName)
    if mode == 'w' and os.path.exists(fileName):
        if overwrite == 1:
            pass
        elif overwrite == -1:
            return ''
        elif overwrite == 0 :
            #NEED TO DELETE FILE
            pass
    if ['win32'].count(sys.platform):
        if re.search('b',mode) is None:
            mode = mode + 'b'
    try:
        f = open(fileName,mode)
    except:
        raise CException(CUtils, 110, str(fileName) + " mode:" + str(mode) + " overwrite:" + str(overwrite) + " exists:" + str(os.path.exists(fileName)))
    return f

#-----------------------------------------------------------------
def saveFile(fileName=None, text=None, text_list=[], overwrite=0):
#-----------------------------------------------------------------
    # If fileName is a CFilePath or CDataFile this will fix it
    fileName = os.path.normpath(str(fileName))
    if os.path.exists(fileName):
        mode = 'w+'
    else:
        mode = 'w'
    f = openFile(fileName,mode,overwrite=overwrite)
    # Try writing to file
    try:
        if len(text_list) > 0:
            if "win32" in sys.platform:
                for l in text_list:
                     if hasattr(l,"encode"):
                         f.write(l.encode("utf-8"))
                     else:
                         f.write(l)
            else:
                f.writelines(text_list)
        elif text is not None:
            if not "win32" in sys.platform:
                if hasattr(text,"decode"):
                    f.write(text.decode())
                else:
                    f.write(text)
            else:
                if hasattr(text,"encode"):
                    f.write(text.encode("utf-8"))
                else:
                    f.write(text)
        f.close()
    except:
        raise CException(CUtils, 107, fileName)

def saveEtreeToFile(tree=None, fileName=None):
    if tree is None:
        raise CException(CUtils, 101, fileName)
    if DEVELOPER():  ##  KJS ** This is the cause of the coupling to CCP4Config.
        # No error trapping - we just want it to crash so we have to fix it!
        text = etree.tostring(tree, pretty_print=True, xml_declaration=True)
        saveFile(fileName=fileName, text=text, overwrite=1)
    else:
        try:
            text = etree.tostring(tree, pretty_print=True, xml_declaration=True)
        except:
            raise CException(CUtils, 102, fileName)
        try:
            saveFile(fileName=fileName, text=text, overwrite=1)
        except:
            raise CException(CUtils, 103, fileName)

utf8_parser = etree.XMLParser(encoding='utf-8')

def parse_from_unicode(unicode_str,useLXML=True):
    if useLXML:
        s = unicode_str.encode('utf-8')
        return etree.fromstring(s, parser=utf8_parser)
    else:
        import xml.etree.ElementTree as etree_xml
        return etree_xml.fromstring(unicode_str)

def openFileToEtree(fileName=None, printout=False,useLXML=True):
    # Use this as etree.parse() seg faults on some Linux
    try:
        f = open(os.path.normpath(fileName))
        s = f.read()
        f.close()
        tree = parse_from_unicode(s,useLXML=useLXML)
    except:
        raise CException(CUtils, 104, fileName)
    else:
        if printout:
            print(etree.tostring(tree, pretty_print=True))
        return tree
    '''
    if DEVELOPER():
        tree = etree.parse(fileName)
    else:
        try:
          tree = etree.parse(fileName)
        except:
          raise CException(CUtils,104,filename)
      try:
        root = tree.getroot()
    except:
        raise CException(CUtils,105,filename)
    if printout: print etree.tostring(root,pretty_print=True)

    return root
    '''

def getHostName():
    # From http://stackoverflow.com/questions/4271740/how-can-i-use-python-to-get-the-system-name
    import socket
    if socket.gethostname().find('.') >= 0:
        name = socket.gethostname()
    else:
        try:
            name = socket.gethostbyaddr(socket.gethostname())[0]
        except:
            try:
                name = socket.gethostbyaddr('localhost')[0]
            except:
                name = socket.gethostbyaddr('127.0.0.1')[0]
    return name

def getUserId():
    name = os.environ.get('LOGNAME', None)
    if name is not None:
        return name
    name = os.environ.get('USERNAME', None)
    if name is not None:
        return name
    import getpass
    name = getpass.getuser()
    if name is not None:
        return name
    try:
        return os.getlogin()
    except:
        return None

#--------------------------------------------------------------------
def getHOME():
#--------------------------------------------------------------------
    '''
    print "HOME", os.environ.get('HOME')
    print "HOMEDRIVE", os.environ.get('HOMEDRIVE')
    print "HOMEPATH",  os.environ.get('HOMEPATH')
    print "USERPROFILE", os.environ.get('USERPROFILE')

    # Expect HOME defined in LINUX and use it if exists in windows
    #homedir = os.environ.get('HOME')
    #if homedir and homedir != '' and os.path.exists(homedir):
    #  return homedir
    if sys.platform == 'win32':
        # Try using HOMEDRIVE & HOMEPATH
        if os.environ.get('HOMEDRIVE'):
          if  os.environ.get('HOMEPATH'):
            homedir = os.path.join(os.environ.get('HOMEDRIVE'), os.environ.get('HOMEPATH'))
          else:
            homedir = os.environ.get('HOMEDRIVE')
        else:
          # Last ditch effort - use USERPROFILE
          homedir = os.environ.get('USERPROFILE')
        if homedir and os.path.exists(homedir):
          print "get_HOME",homedir
          return homedir
    
    # Arghh - no HOME
    return ""
    '''
    homedir = os.environ.get('CCP4_LOCAL_HOME', None)
    if homedir is not None:
        return homedir
    try:
        if sys.platform == 'win32':
            homedir = os.path.normpath(os.environ.get('USERPROFILE'))
        else:
            homedir = os.environ.get('HOME')
    except:
        homedir = ''
    if homedir and os.path.exists(homedir):
        return homedir
    else:
        return ''

def getTMP():
    # Return a temp directory
    #return os.path.join(getDotDirectory(),'tmp')
    import tempfile
    return tempfile.gettempdir()

def getTestTmpDir():
    tmp = os.environ.get('CCP4I2_TEST', None)
    if tmp is not None and not os.path.exists(tmp):
        try:
            os.mkdir(tmp)
        except:
            pass
    if tmp is None or not os.path.exists(tmp):
        tmp = os.path.join(getHOME(),'CCP4I2_TEST')
    if not os.path.exists(tmp):
        os.mkdir(tmp)
    for subDir in ['CCP4_JOBS','CCP4_PROJECT_FILES']:
        if not os.path.exists(os.path.join(tmp, subDir)):
            os.mkdir(os.path.join(tmp, subDir))
    return tmp

def makeTmpFile(name='tmp', extension='tmp', mode=None, cdir=False):
    dotExtension = ''
    if extension:
        dotExtension = '.' + extension
    t = int(time.time())
    n = 0
    fileName = None
    while fileName is None or os.path.exists(fileName):
        n = n + 1
        fileName = os.path.join(getTMP(), name + '_' + str(t) + '_' + str(n) + dotExtension)
    # Make the file to be sure if is reserved
    #print 'CCP4Utils.makeTmpFile',fileName
    if cdir:
        os.mkdir(fileName)
        return fileName
    elif mode is not None:
        fd = open(fileName, mode)
        return (fd, fileName)
    else:
        return fileName

def backupFile(fileName=None, delete=False):
    fileName = os.path.normpath(fileName)
    label = 'previous'
    labelLen = len(label)
    if not os.path.exists(fileName):
        return None
    backup = 0
    ext1 = ''
    ext2 = ''
    base,ext0 = os.path.splitext(fileName)
    base,ext1 = os.path.splitext(base)
    if ext1[0:] == label:
        try:
            backup = int(ext1[labelLen+1:])
        except:
            pass
    elif len(ext1) > 0:
        ext2 = os.path.splitext(base)
        if  ext2[0:labelLen] == label:
            try:
                backup = int(ext2[labelLen+1:])
            except:
                pass
    #print 'CCP4Utils.backupFile',base,ext0,ext1,ext2,backup
    backup = backup + 1
    newName = base + '.' + label + '_' + str(backup)+ext1+ext0
    while os.path.exists(newName):
        backup = backup + 1
        newName = base + '.' + label + '_' + str(backup)+ext1+ext0
    #print 'CCP4Utils.backupFile',newName
    if sys.platform == "win32" or not delete:
        shutil.copyfile(fileName, newName)
    else:
        shutil.move(fileName, newName)
    return newName

def splitPath(path):
    pathList= []
    while len(path) > 0:
        path,tail = os.path.split(path)
        pathList.append(tail)
    pathList.reverse()
    return 

def pythonExecutable():
    print('into CCP4Utils.pythonExecutable')
    # There is also os.path.join(os.environ["CCP4"],"libexec","ccp4i2")
    if os.path.exists(os.path.join(os.environ["CCP4"], "bin", "ccp4-python")):
        return os.path.join(os.environ["CCP4"], "bin", "ccp4-python")
    elif 'PYTHONHOME' in os.environ:
        return os.path.join(os.environ['PYTHONHOME'], 'bin', 'python')
    else:
        return "python"
    '''
    if sys.platform == 'darwin':
        return os.path.join(getCCP4I2Dir(),'bin','Python')
    else:
        return os.path.join(os.environ['CCP4MG'],'pythondist','bin','python')
    '''

'''  
def getCCP4I2Dir(up=1):
  target = os.environ.get('CCP4I2_TOP',None)
  #print 'CCP4Utils.getCCP4I2Dir environ',target
  if target is not None:
    return os.path.abspath(target)
  else:
    target = os.path.join(os.path.abspath(sys.argv[0]),"..")
    abstarget = os.path.abspath(target)
    splittarget = abstarget.split('/')
    if splittarget.count('ccp4i2'):
      splittarget.reverse()
      up = splittarget.index('ccp4i2')
    while up>0:
      abstarget = os.path.dirname(abstarget)
      up = up -1
    #print 'CCP4Utils.getCCP4I2Dir Environment variable CCP4I2_TOP not defined - trying ',abstarget
    return abstarget
''' 

def getCCP4I2Dir(**kw):
    f = os.path.normpath(__import__('core.CCP4Utils').__file__)
    return os.path.split(os.path.split(f)[0])[0]

def getOSDir():
    if sys.platform == 'darwin':
        return os.path.join(getCCP4I2Dir(), 'MacOSX')
    elif sys.platform[0:5] == 'linux':
        return os.path.join(getCCP4I2Dir(), 'Linux')
    elif sys.platform[0:3] == 'win':
        return os.path.normpath(os.path.join(getCCP4I2Dir(), 'Windows'))
    else:
        return ''

def getTestSysDir():
    return os.path.join(getCCP4I2Dir(), 'testsysdefs')

def getCCP4Dir():
    try:
        target = os.path.normpath(os.environ.get('CCP4', None))
    except:
        target = ''
    return target

def interpretPath(path='', currentDir=None):
    path = path.strip()
    # print 'interpretPath',path
    m = re.search(r'(.*?)\$(.*?)(\/|$)(.*)', path)
    while m is not None:
        #print 'interpretPath m',m.groups()
        reGroups = m.groups()
        envPath = os.environ.get(reGroups[1], '')
        path = reGroups[0] + os.environ.get(reGroups[1], '') + reGroups[2] + reGroups[3]
        m = re.search(r'(.*?)\$(.*?)(\/|$)(.*)', path)
    if currentDir is not None and path[0] == '.':
        if currentDir[-1] == '/' and path[0] == '/':
            path = currentDir[0:-1] + path
        elif currentDir[-1] != '/' and path[0] != '/':
            path = currentDir + '/' + path
        else:
            path = currentDir + path
    return os.path.abspath(path)

def getDotDirectory():
    home = os.environ.get('CCP4_LOCAL_DOTDIR', None)
    if home is None:
        home = getHOME()
    ccp4i2 = os.path.join (home,'.CCP4I2')
    if not os.path.exists(ccp4i2) and sys.platform == 'win32':
        ccp4i2 = os.path.normpath(os.path.join(home, 'CCP4I2'))
    if not os.path.exists(ccp4i2):
        try:
            os.mkdir(ccp4i2)
        except:
            print("ERROR creating .CCP4I2 directory")
            return ""
    # Make sure that we also have the subdirectories
    for subd in ['status', 'logs', 'configs', 'db', 'tmp', 'custom', 'demo_data', 'i1supplement']:
        path = os.path.join(ccp4i2, subd)
        if  not os.path.exists(path):
            try:
                os.mkdir(path)
            except:
                print("ERROR creating subdirectories in ", ccp4i2)
    for subd in ['workflows', 'comfilepatchs', 'importedjobs', 'tasks']:
        path = os.path.join(ccp4i2, 'custom', subd)
        if  not os.path.exists(path):
            try:
                os.mkdir(path)
            except:
                print("ERROR creating subdirectories in ",ccp4i2)
    #target =  os.path.join(ccp4i2 ,'configs', 'ccp4i2_config.params.xml')
    #if not os.path.exists(target):
    #  source = os.path.join(getCCP4I2Dir(),'data','config_templates','ccp4i2_config.params.xml')
    #  shutil.copyfile(source,target)
    return ccp4i2

def getProjectDirectory():
    ccp4i2 =  os.path.join(getHOME(),'CCP4I2_PROJECTS')
    if not os.path.exists(ccp4i2):
        try:
            os.mkdir(ccp4i2)
        except:
            print("ERROR creating CCP4I2_PROJECTS directory")
            return None
    return ccp4i2

def globSearchPath(searchPath=[], cfile='*'):
    fileList = []
    for path in searchPath:
        fileList.extend(glob.glob(os.path.join(path, cfile)))
    return fileList

def importFileModule(pyFile, report=False):
    try:
        import sys, os
        extraPath = os.path.split(pyFile)[0]
        sys.path.insert(0, extraPath)
        moduleName = os.path.splitext(os.path.basename(pyFile))[0]
        #print 'importFileModule moduleName', moduleName
        module = __import__(moduleName)
        sys.path.remove(extraPath)
        return module, None
    except Exception as e:
        if report:
            print('ERROR attempting to import module: ', moduleName, e)
        module = None
        return module, str(e)

def writeTarGzip(directory=None, tarFile=None):
    if tarFile is None:
        tarFile = directory + '.tar.gz'
    try:
        tarObj = tarfile.open(tarFile,'w:gz')
        #tarObj.add(directory,arcname=os.path.split(directory)[1],exclude=excludeFromTarGzip)
        tarObj.add(directory,arcname=os.path.split(directory)[1])
        tarObj.close()
    except:
        return None
    return tarFile

def excludeFromTarGzip(fileName):
    #print 'excludeFromTarGzip',fileName
    return False

def readTarGzip(fileName, destination=None):
    if destination is None:
        destination = os.path.split(fileName)[0]
    tarObj = tarfile.open(fileName)
    next = tarObj.next().name
    print('CCP4Utils.  readTarGzip next',next)
    tarObj.extractall(path=destination)
    tarObj.close()
    return next

def getProgramVersion(programName, mode='version'):
    # Run CCP4 program with -i to get either explicit version info
    # or version as part of the program banner header which appears before program fails
    from core import CCP4Modules
    bin = copy.deepcopy(programName)
    programName = programName.lower()
    logFile = makeTmpFile(extension='log')
    if programName == 'ccp4i2':
        from core import CCP4File
        versionHeader = CCP4File.CI2XmlHeader()
        versionHeader.loadFromXml(os.path.join(getCCP4I2Dir(), 'core', 'version.params.xml'))
        if mode == 'version':
            return versionHeader.ccp4iVersion.__str__()
        elif mode == 'date':
            return versionHeader.creationTime.date()
        elif mode == 'revision':
            text = readFile(os.path.join(getCCP4I2Dir(), 'core', 'version.params.xml'))
            tree = etree.fromstring(text)
            if os.path.exists(os.path.join(getCCP4I2Dir(),".bzr","branch","last-revision")):
                try:
                    bzrrevf = open(os.path.join(getCCP4I2Dir(),".bzr","branch","last-revision"))
                    bzrrevt = bzrrevf.read()
                    bzrrevf.close()
                    return bzrrevt.split(" ")[0]
                except:
                    pass #Read from version.params.xml
            #Aaargh, pluginVersion is a CVersion in CI2XmlHeader, but does not adhere to CVersion rules...
            try:
                plugVer = tree.xpath('/ccp4i2/ccp4i2_header/pluginVersion')[0].text
                return plugVer
            except:
                return ""
    elif programName == 'ccp4':
        CCP4Modules.PROCESSMANAGER().startProcess('fft', ['-i'], logFile=logFile)
        text = readFile(logFile)
        m1 = re.search('(.*)patch level(.*)', text)
        if m1 is None:
            return None
        text = m1.groups()[1].strip()
        return text
    elif programName == 'python':
        return sys.version.split()[0]
    elif programName == 'qt':
        from PySide2 import QtCore
        return QtCore.qVersion()
    elif programName == 'arp_warp':
        # Need to get the it exit
        return None
    elif programName.startswith('phenix'):
        testArgList = [['--version'],['-version']]
    else:
        testArgList = [['-i'],['--version'],['-version']]
    for argList in testArgList:
        CCP4Modules.PROCESSMANAGER().startProcess(bin, argList, logFile=logFile)
        #print 'CCP4Utils.getProgramVersion programName',programName,logFile
        text = readFile(logFile)
        if argList == ['-i']:
            ret = searchVersion(text, programName)
        else:
            ret = searchVersion(text)
        if ret is not None:
            return ret

def searchVersion(text, programName=None):
    if programName is not None:
        #Split text at program name - for CCP4 progs expect the version on same line
        m1 = re.search('(.*)' + programName + '(.*)', text)
        if m1 is not None:
            text = m1.groups()[1]
        print('CCP4Utils.getProgramVersion m1', m1.groups())
    m2 = re.search('(.*)(V|v)ersion([ :\.]*)([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[3]) > 0:
        return m2.groups()[3]
    m2 = re.search('(.*)(V|v)er([ :\.]*)([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[3]) > 0:
        return m2.groups()[3]
    m2 = re.search('(.*)([V|v])([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[2]) > 0:
        return m2.groups()[2]
    return None

def versionLogHeader():
    text = 'CCP4i2 version: ' + getProgramVersion('ccp4i2') + '\n' + \
            'Running CCP4 version: ' + getProgramVersion('ccp4') + '\n' + \
            'Using Python version: ' + sys.version + '\n' + \
            'Using Qt version: ' + getProgramVersion('qt') + '\n'
    return text

def listReMatch(lst,reExp):
    # search list of strings with reg expression and return hit string and index
    # More genius code brought to you by
    # http://stackoverflow.com/questions/23229675/search-list-of-string-with-wildcard-and-return-matches-indexes-python
    regex=re.compile(reExp)
    locs, matches = list(zip(*[(idx, string) for idx, string in enumerate(lst) if re.match(regex, string)]))
    #print 'listSearch matches', matches
    #print 'listSearch locs', locs
    return locs, matches

def isAlive(qobj):
    import shiboken2
    return shiboken2.isValid(qobj)

# Slightly modified from http://timgolden.me.uk/python/win32_how_do_i/see_if_two_files_are_the_same_file.html
def get_read_handle(filename):
    import win32file
    if os.path.isdir(filename):
        dwFlagsAndAttributes = win32file.FILE_FLAG_BACKUP_SEMANTICS
    else:
        dwFlagsAndAttributes = 0
    return win32file.CreateFile (
        filename,
        win32file.GENERIC_READ,
        win32file.FILE_SHARE_READ,
        None,
        win32file.OPEN_EXISTING,
        dwFlagsAndAttributes,
        None)

def get_unique_id(hFile):
    import win32file
    (attributes, created_at, accessed_at, written_at, volume, file_hi, file_lo, n_links, index_hi, index_lo) = win32file.GetFileInformationByHandle (hFile)
    return volume, index_hi, index_lo

def files_are_equal (filename1, filename2, default=False):
    try:
        hFile1 = get_read_handle (filename1)
    except:
        return default
    try:
        hFile2 = get_read_handle (filename2)
    except:
        return default
    try:
        are_equal = (get_unique_id (hFile1) == get_unique_id (hFile2))
        hFile2.Close ()
        hFile1.Close ()
    except:
        return default
    return are_equal

def samefile(f1, f2, default=False):
    if sys.platform != 'win32':
        return os.path.samefile(f1,f2)
    else:
        return files_are_equal(f1,f2,default=default)

def zipDirectory(czip, sourceDirectory, rootRelPath=None):
    if rootRelPath is None:
        rootRelPath=sourceDirectory
    # Beware the 'sourceDirectory' could be a file
    if os.path.isfile(sourceDirectory):
        root, cfile = os.path.split(sourceDirectory)
        arcname = os.path.join(os.path.relpath(root, rootRelPath ), cfile)
        czip.write(sourceDirectory, arcname)
    else:
        for root, dirs, files in os.walk(sourceDirectory):
            # add directory (needed for empty dirs)
            # Beware relpath can barf on Windows if comparing paths on differnt drives
            # - this should not be the case here
            czip.write(root, os.path.relpath(root, rootRelPath))
            for cfile in files:
                filename = os.path.join(root, cfile)
                if os.path.isfile(filename): # regular files only
                    arcname = os.path.join(os.path.relpath(root,rootRelPath ), cfile)
                    czip.write(filename, arcname)

def findCootPath():
    # This is the default install path for ccp4 6.4.0
    # Except windows release does not include coot and this is default install for coot
    if path is not None and os.path.exists(path):   # KJS : Well this isn't right. Not used anywhere either.
        return path
    if sys.platform == 'win32':
        path = os.path.normpath(os.path.join('C:','WinCoot','runwincoot'))
    elif sys.platform == 'darwin':
        path = os.path.normpath(os.path.join(os.environ['CCP4'],'coot'))
    else:
        path = os.path.normpath(os.path.join(os.environ['CCP4'], 'bin', 'coot'))
    if os.path.exists(path):
        return path
    return None

def getCCP4Exe(program):
    if sys.platform == 'win32':
        exe = program + '.exe'
    else:
        exe = program
    ccp4bin = os.path.join( getOSDir(), 'bin', exe)
    if not os.path.exists(bin):
        ccp4bin = os.path.join(getCCP4Dir(), 'bin', exe)
    return ccp4bin

        
def which(program, mode=os.F_OK | os.X_OK, path=None):
    '''
    From shutil.which in python 3.3 back ported for local use here in 2.7
    via Tom Burnley and the CCP4EM GUI.
    Given a command, mode, and a PATH string, return the path which
    conforms to the given mode on the PATH, or None if there is no such
    file.
    `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
    of os.environ.get("PATH"), or can be overridden with a custom search
    path.
    '''
    # Check ccpem settings for specified bins.  E.g. if bin is set via alias
    # it will not be located via which function.
    cmd = program
    # Check that a given file can be accessed with the correct mode.
    # Additionally check that `file` is not a directory, as on Windows
    # directories pass the os.access check.
    def _access_check(fn, mode):
        return (os.path.exists(fn) and os.access(fn, mode)
                and not os.path.isdir(fn))
    # If we're given a path with a directory part, look it up directly rather
    # than referring to PATH directories. This includes checking relative to the
    # current directory, e.g. ./script
    if os.path.dirname(cmd):
        if _access_check(cmd, mode):
            return cmd
        return None
    if path is None:
        path = os.environ.get("PATH", os.defpath)
    if not path:
        return None
    path = path.split(os.pathsep)
    if sys.platform == "win32":
        # The current directory takes precedence on Windows.
        if not os.curdir in path:
            path.insert(0, os.curdir)
        # PATHEXT is necessary to check on Windows.
        pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
        # See if the given file matches any of the expected path extensions.
        # This will allow us to short circuit when given "python.exe".
        # If it does match, only test that one, otherwise we have to try
        # others.
        if any([cmd.lower().endswith(ext.lower()) for ext in pathext]):
            files = [cmd]
        else:
            files = [cmd + ext for ext in pathext]
    else:
        # On other platforms you don't have things like PATHEXT to tell you
        # what file suffixes are executable, so just pass on cmd as-is.
        files = [cmd]
    seen = set()
    for cdir in path:
        normdir = os.path.normcase(cdir)
        if not normdir in seen:
            seen.add(normdir)
            for thefile in files:
                name = os.path.join(cdir, thefile)
                if _access_check(name, mode):
                    return name
    return None

def nonWhiteDifferences(file1, file2):
    retDiffs = []
    retCount = []
    if sys.version_info > (3,0):
        from googlecode import diff_match_patch_py3
        dmp =  diff_match_patch_py3.diff_match_patch()
    else:
        from googlecode import diff_match_patch
        dmp =  diff_match_patch.diff_match_patch()
    
    text1= readFile(file1)
    text2 = readFile(file2)
    text1 = re.sub(r' +\n','\n',text1)
    text2 = re.sub(r' +\n','\n',text2)
    diffs = dmp.diff_main(text1,text2)
    #dmp.diff_cleanupSemantic(diffs)
    for n in range(len(diffs)):
        if diffs[n][0] != 0 and len(diffs[n][1].strip())>0 :
            d = [ diffs[n][0],diffs[n][1] ]
            if retDiffs.count(d):
                idx = retDiffs.index(d)
                retCount[idx] += 1
            else:
                retDiffs.append(d)
                retCount.append(1)
    for n in range(len(retDiffs)):
        retDiffs[n].append(retCount[n])
    #print 'nonWhiteDifferences',retDiffs
    return retDiffs

#===================================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testEtreeTools)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testFileOpen))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testBackup))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testEtreeTools(unittest.TestCase):

    def test1(self):
        fileName = os.path.join(getCCP4I2Dir(),'test','data','test_job_104.def.xml')
        root = openFileToEtree(fileName=fileName)
        self.assertEqual(root.tag,'ccp4i2','Failed to read first item in eTree')
 
    def test2(self):
        fileName = os.path.join(getTestTmpDir(),'CCP4Utils_test.def.xml')
        from core import CCP4Data
        x = CCP4Data.CFloat(42,default=0)
        ele = x.getEtree()
        saveEtreeToFile(tree=ele,fileName=fileName)
        if not os.path.exists(fileName):
            self.fail('Failed to write eTree file')

class testFileOpen(unittest.TestCase):

    def test1(self):
        fileName = os.path.join(getCCP4I2Dir(),'test','data','test_job_104.def.xml')
        text = readFile(fileName)
        textLines = text.split('\n')
        self.assertEqual(textLines[0],"<?xml version='1.0'?>",'Error reading file with readFile')

    def test2(self):
        text = 'The quick brown fox jumped over the lazy dog'
        fileName = os.path.join(getTMP(),'fox.dog')
        saveFile(fileName,text)
        self.assertTrue(os.path.exists(fileName),'Failed to save file with saveFile')

    def test3(self):
        text = 'The quick brown fox jumped over the lazy dog'
        fileName = '/foo/bar/fox.dog'
        try:
            saveFile(fileName,text)
        except CException as e:
            #print 'testFileOpen.test3',e.report()
            self.assertEqual(e[0]['code'],110,'saveFile gives unexpected CException code')
        except:
            self.fail('saveFile gives unexpected Python exception')

class testBackup(unittest.TestCase):
    def test1(self):
        tempdir = getTMP()
        fileName = os.path.join(tempdir,'fox.PARAMS.xml')
        try:
            os.remove(os.path.join(tempdir,'fox.backup_1.PARAMS.xml'))
        except:
            pass
        saveFile(fileName,'whatever')
        newFile = backupFile(fileName=fileName,delete=True)
        self.assertEqual(newFile,os.path.join(tempdir,'fox.backup_1.PARAMS.xml'),'backupFile gives wrong file name')   
        if not os.path.exists(newFile):
            self.fail('backupFile does not create new file')
        if os.path.exists(fileName):
            self.fail('backupFile old file still exists')
