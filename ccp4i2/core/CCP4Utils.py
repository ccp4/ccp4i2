"""
Liz Potterton Jan 2010 - Copied from ccp4mg python/ui/utils.py and converted to Qt
"""

import importlib
import getpass
import glob
import os
import re
import shutil
import socket
import sys
import tarfile
import tempfile
import time
import xml.etree.ElementTree as ET

import shiboken2

from .. import I2_TOP
from ..googlecode import diff_match_patch_py3
from .CCP4ErrorHandling import CException


def printXml(element, pretty_print=True):
    if isinstance(element, ET.ElementTree):
        element = element.getroot()
    if pretty_print:
        ET.indent(element)
    print(ET.tostring(element))


def writeXml(tree, file_or_filename, pretty_print=True, encoding=None, xml_declaration=None):
    if isinstance(tree, ET.Element):
        tree = ET.ElementTree(tree)
    if pretty_print:
        ET.indent(tree)
    tree.write(file_or_filename, encoding=encoding, xml_declaration=xml_declaration)


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


def safeOneWord(value):
    # Replace spaces and any 'odd' characters with underscore
    # Resulting word should be safe in a filename
    one_word_name = re.sub('[^a-zA-Z0-9_-]', '_', value)
    return one_word_name


def readFile(fileName=None, limit=None):
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


def openFile(fileName=None, mode='r', overwrite=1):
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


def saveFile(fileName=None, text=None, text_list=[], overwrite=0):
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
    try:
        ET.indent(tree)
        text = ET.tostring(tree, xml_declaration=True)
    except:
        raise CException(CUtils, 102, fileName)
    try:
        saveFile(fileName=fileName, text=text, overwrite=1)
    except:
        raise CException(CUtils, 103, fileName)


def getHostName():
    # From http://stackoverflow.com/questions/4271740/how-can-i-use-python-to-get-the-system-name
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
    name = getpass.getuser()
    if name is not None:
        return name
    try:
        return os.getlogin()
    except:
        return None


def getHOME():
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


def pythonExecutable():
    print('into CCP4Utils.pythonExecutable')
    # There is also os.path.join(os.environ["CCP4"],"libexec","ccp4i2")
    if os.path.exists(os.path.join(os.environ["CCP4"], "bin", "ccp4-python")):
        return os.path.join(os.environ["CCP4"], "bin", "ccp4-python")
    if 'PYTHONHOME' in os.environ:
        return os.path.join(os.environ['PYTHONHOME'], 'bin', 'python')
    return "python"


def getCCP4I2Dir():
    return str(I2_TOP)


def getCCP4Dir():
    try:
        return os.path.normpath(os.environ.get('CCP4'))
    except:
        return ''


def getDotDirectory():
    home = os.environ.get('CCP4_LOCAL_DOTDIR', default=getHOME())
    ccp4i2 = os.path.join(home, '.CCP4I2')
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
        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except:
                print("ERROR creating subdirectories in", ccp4i2)
    for subd in ['workflows', 'comfilepatchs', 'importedjobs', 'tasks']:
        path = os.path.join(ccp4i2, 'custom', subd)
        if  not os.path.exists(path):
            try:
                os.mkdir(path)
            except:
                print("ERROR creating subdirectories in ", ccp4i2)
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


def importModule(name):
    return importlib.import_module(f"ccp4i2.{name}")


def importFileModule(pyFile, report=False):
    try:
        extraPath = os.path.split(pyFile)[0]
        sys.path.insert(0, extraPath)
        moduleName = os.path.splitext(os.path.basename(pyFile))[0]
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
        tarObj.add(directory,arcname=os.path.split(directory)[1])
        tarObj.close()
    except:
        return None
    return tarFile


def searchVersion(text, programName=None):
    if programName is not None:
        #Split text at program name - for CCP4 progs expect the version on same line
        m1 = re.search('(.*)' + programName + '(.*)', text)
        if m1 is not None:
            text = m1.groups()[1]
        print('CCP4Utils.getProgramVersion m1', m1.groups())
    m2 = re.search(r'(.*)(V|v)ersion([ :\.]*)([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[3]) > 0:
        return m2.groups()[3]
    m2 = re.search(r'(.*)(V|v)er([ :\.]*)([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[3]) > 0:
        return m2.groups()[3]
    m2 = re.search('(.*)([V|v])([0123456789.]*)',text)
    if m2 is not None and len(m2.groups()[2]) > 0:
        return m2.groups()[2]
    return None


def isAlive(qobj):
    return shiboken2.isValid(qobj)


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


def nonWhiteDifferences(file1, file2):
    retDiffs = []
    retCount = []
    dmp =  diff_match_patch_py3.diff_match_patch()

    text1= readFile(file1)
    text2 = readFile(file2)
    text1 = re.sub(r' +\n','\n', text1)
    text2 = re.sub(r' +\n','\n', text2)
    diffs = dmp.diff_main(text1, text2)
    for diff in diffs:
        if diff[0] != 0 and len(diff[1].strip()) > 0:
            d = [diff[0], diff[1]]
            if retDiffs.count(d):
                idx = retDiffs.index(d)
                retCount[idx] += 1
            else:
                retDiffs.append(d)
                retCount.append(1)
    for n, retDiff in enumerate(retDiffs):
        retDiff.append(retCount[n])
    return retDiffs
