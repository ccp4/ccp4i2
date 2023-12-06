from __future__ import print_function

"""
     CCP4CustomMimeTypes.py: CCP4 GUI Project
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
     Liz Potterton Jan 2010 - Create CCP4MimeTypes. Copy MGAbstractViewer and subclasses from MG project
"""
## @package CCP4CustomMimeTypes (QtWebKit) Manager for file MIME types - access via CCP4Modules.MIMETYPESHANDLER()
import os
from core.CCP4Config import GRAPHICAL
from core.CCP4ErrorHandling import *
from PySide2 import QtCore
from core import CCP4Utils

#-------------------------------------------------------------------------------------

# Extension of Qts MimeType class - beware they may be developing this
# further.  We really just want to know whether to use our viewer class
# or the Qt DesktpServices class to handle display

if GRAPHICAL():
    from PySide2 import QtGui

if False:
    from PyQt4 import QtWebKit, QtGui

    class CMimeType(QtWebKit.QWebPluginFactory.MimeType):
        def __init__(self):
            QtWebKit.QWebPluginFactory.MimeType.__init__(self)
            self.viewer = ''                   # our viewer class - subclass of CAbstractViewer
            self.useDesktopServices = 0        # use QtDesktopServices to display
            self.useWebBrowser = 0             # use CWebBrowser to display
            self.icon = ''                     # icon from ccp4i2/qticons
            self.fixedWidthFont = False        # display with fixed width font
            self.fileValidity = None           # method to test file validity - returns CErrorReport
            self.contentLabel = None           # terminal part of basename indicating file content
            self.className = None              # name (without initial C) of the CDataFile subclass

else:
    class CMimeType:
        def __init__(self):
            self.name = ''
            self.description = ''
            self.fileExtensions = []
            self.viewers = []
            self.viewer = ''                   # our viewer class - subclass of CAbstractViewer
            self.useDesktopServices = 0        # use QtDesktopServices to display
            self.useWebBrowser = 0             # use CWebBrowser to display
            self.icon = ''                     # icon from ccp4i2/qticons
            self.fixedWidthFont = False        # display with fixed width font
            self.fileValidity = None           # method to test file validity - returns CErrorReport
            self.contentLabel = None           # terminal part of basename indicating file content
            self.className = None              # name (without initial C) of the CDataFile subclass

#--------------------------------------------------------------------------------------

# Container for CMimeType classes
class CCustomMimeTypes(QtCore.QObject):
    insts = None
    ICONS = {}
    ERROR_CODES = {
          101 : {'description' : 'File format not recognised'},
          102 : {'severity' : SEVERITY_WARNING, 'description' : 'No validity test method for format'},
          103 : {'description' : 'File does not exist'},
          104 : {'severity' : SEVERITY_WARNING, 'description' : 'UNKNOWN ERROR running validity test'},
          105 : {'severity' : SEVERITY_WARNING, 'description' : 'File has inappropriate extension'},
          106 : {'description' : 'File is not this format'},
          107 : {'severity' : SEVERITY_WARNING, 'description' : 'File not strictly valid but may be usable'}}

    def __init__(self):
        QtCore.QObject.__init__(self)
        if not CCustomMimeTypes.insts:
            CCustomMimeTypes.insts = self
        self.mimeTypes = {}
        self.setupStandards()
        self.miniMtzMimeTypes = ["application/CCP4-mtz-observed","application/CCP4-mtz-phases","application/CCP4-mtz-map","application/CCP4-mtz-freerflag"]
        # Cache mimeType of files that need opening to find the file type
        self.knownCifFiles = {}
        self.knownEntFiles = {}

    def setupStandards(self):
        if GRAPHICAL():
            from qtgui import CCP4TextViewer
            from qtgui import CCP4FileSystemView
            from qtgui import CCP4ImageViewer

        mimeType = CMimeType()
        mimeType.name = "text/html"
        mimeType.description = "Hypertext markup language"
        mimeType.fileExtensions = ['html','htm']
        mimeType.viewers = []
        mimeType.useWebBrowser = 1
        self.mimeTypes["text/html"] = mimeType

        mimeType = CMimeType()
        if GRAPHICAL():
            mimeType.viewers = [CCP4FileSystemView.CFileSystemView]
        mimeType.name = "dir"
        mimeType.description = "Directory"
        mimeType.fileExtensions = []
        mimeType.useWebBrowser = 0
        self.mimeTypes["dir"] = mimeType

        mimeType = CMimeType()
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.name = "text/plain"
        mimeType.description = "Standard plain text"
        mimeType.fileExtensions = ['txt','log','csv','py','scm','def','sh']
        mimeType.fixedWidthFont = True
        self.mimeTypes["text/plain"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/postscript"
        mimeType.description = "Postscript file"
        mimeType.fileExtensions = ['ps']
        mimeType.icon = 'document_blank'
        mimeType.useDesktopServices = 1
        self.mimeTypes["application/postscript"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/x-pdf"
        mimeType.description = "PDF file"
        mimeType.fileExtensions = ['pdf']
        mimeType.icon = 'PDFDataFile'
        mimeType.useDesktopServices = 1
        self.mimeTypes[mimeType.name] = mimeType

        mimeType = CMimeType()
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.name = "application/coot-script"
        mimeType.icon = 'CootHistoryDataFile'
        mimeType.description = "Coot script"
        mimeType.fileExtensions = ['py','scm']
        mimeType.fixedWidthFont = False
        self.mimeTypes["application/coot-script"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "chemical/x-pdb"
        mimeType.description = "Model coordinates"
        mimeType.fileExtensions = ['pdb','cif','ent']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CCoordsViewer]
        mimeType.icon = 'PdbDataFile'
        mimeType.fixedWidthFont = True
        mimeType.className = 'PdbDataFile'
        self.mimeTypes["chemical/x-pdb"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "chemical/x-pdb-ensemble"
        mimeType.description = "Ensemble of model coordinates"
        mimeType.fileExtensions = ['pdb','cif','ent']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CCoordsViewer]
        mimeType.icon = 'EnsemblePdbDataFile'
        mimeType.fixedWidthFont = True
        mimeType.className = 'EnsemblePdbDataFile'
        self.mimeTypes["chemical/x-pdb-ensemble"] = mimeType

        if GRAPHICAL():
            image_formats_c = QtGui.QImageReader.supportedImageFormats()
            image_formats = []
            for im in image_formats_c:
                if not image_formats.count(str(im).lower()):
                    image_formats.append(str(im))

            mimeType = CMimeType()
            mimeType.name = "image"
            mimeType.description = "2D image"
            mimeType.fileExtensions = image_formats
            mimeType.viewers = [CCP4ImageViewer.CImageViewer]
            self.mimeTypes["image"] = mimeType

            image_formats_c = QtGui.QMovie.supportedFormats()
            image_formats = []
            for im in image_formats_c:
                if not image_formats.count(str(im).lower()):
                    image_formats.append(str(im))
            #print 'CCustomMimeTypes.setupStandards image_formats',image_formats_c,image_formats
            # The QMovie widget supports animated gif and png
            mimeType = CMimeType()
            mimeType.name = "animation"
            mimeType.description = "Animation"
            mimeType.fileExtensions = image_formats
            mimeType.viewers = [CCP4ImageViewer.CAnimationViewer]
            self.mimeTypes["image"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/pdf"
        mimeType.description = "Portable Document Format document"
        mimeType.fileExtensions = ['pdf']
        mimeType.viewers = []
        mimeType.useDesktopServices = 1
        self.mimeTypes["application/pdf"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "video/mp4"
        mimeType.description = "MP4 Video"
        mimeType.fileExtensions = ['mp4']
        mimeType.viewers = []
        mimeType.useDesktopServices = 1
        self.mimeTypes["video/mp4"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "video/mpg"
        mimeType.description = "MPG Video"
        mimeType.fileExtensions = ['mpg']
        mimeType.viewers = []
        mimeType.useDesktopServices = 1
        self.mimeTypes["video/mp4"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz"
        mimeType.description = "MTZ experimental data"
        mimeType.fileExtensions = ['mtz','cif','ent']
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'MTZDataFile'
        mimeType.className = 'MtzDataFile'
        self.mimeTypes["application/CCP4-mtz"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz-mini"
        mimeType.description = "Experimental data object"
        mimeType.fileExtensions = ['mtz','cif','ent']
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'MTZDataFile'
        mimeType.className = 'MtzDataFile'
        self.mimeTypes["application/CCP4-mtz-mini"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz-observed"
        mimeType.description = "Reflections"
        mimeType.fileExtensions = ['mtz','cif']
        mimeType.content = 'observed_data'
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'ObsDataFile'
        mimeType.className = 'ObsDataFile'
        self.mimeTypes["application/CCP4-mtz-observed"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz-phases"
        mimeType.description = "Phases"
        mimeType.contentLabel = 'phases'
        mimeType.fileExtensions = ['mtz','cif']
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'PhsDataFile'
        mimeType.className = 'PhsDataFile'
        self.mimeTypes["application/CCP4-mtz-phases"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz-map"
        mimeType.description = "Map coefficients"
        mimeType.fileExtensions = ['mtz','cif']
        mimeType.contentLabel = 'map_coefficients'
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'MapDataFile'
        mimeType.className = 'MapDataFile'
        self.mimeTypes["application/CCP4-mtz-map"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-mtz-freerflag"
        mimeType.description = "Free R set"
        mimeType.fileExtensions = ['mtz','cif']
        mimeType.contentLabel = 'freeRflag'
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'FreeRDataFile'
        mimeType.className = 'FreeRDataFile'
        self.mimeTypes["application/CCP4-mtz-freerflag"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-shelx-FA"
        mimeType.description = "Shelx FA"
        mimeType.fileExtensions = ['hkl']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = ''
        mimeType.className = 'ShelxFADataFile'
        self.mimeTypes["application/CCP4-shelx-FA"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-unmerged-mtz"
        mimeType.description = "MTZ unmerged experimental data"
        mimeType.fileExtensions = ['mtz']
        mimeType.viewers =['viewhkl']
        mimeType.icon = 'MTZDataFile'
        self.mimeTypes["application/CCP4-unmerged-mtz"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-unmerged-experimental"
        mimeType.description = 'Unmerged experimental data'
        mimeType.fileExtensions = ['mtz','hkl','HKL','sca','SCA']
        mimeType.viewers =['viewhkl',CCP4TextViewer.CTextViewer]
        mimeType.fixedWidthFont = True
        mimeType.icon = 'UnmergedDataFile'
        mimeType.className = 'UnmergedDataFile'
        self.mimeTypes["application/CCP4-unmerged-experimental"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-generic-reflections"
        mimeType.description = 'Merged experimental data'
        mimeType.fileExtensions = ['mtz','hkl','HKL','sca','SCA','mmcif','cif','ent']
        mimeType.viewers =['viewhkl',CCP4TextViewer.CTextViewer]
        mimeType.fixedWidthFont = True
        mimeType.icon = 'GenericReflDataFile'
        mimeType.className = 'GenericReflDataFile'
        self.mimeTypes["application/CCP4-generic-reflections"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "chemical/x-cif"
        mimeType.description = 'mmCif reflection data'
        mimeType.fileExtensions = ['cif','mmcif']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.fixedWidthFont = True
        self.mimeTypes["chemical/x-cif"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-map"
        mimeType.description = "Electron density map"
        mimeType.fileExtensions = ['map']
        mimeType.viewers = ['coot','ccp4mg']
        mimeType.icon = 'MapDataFile'
        mimeType.className = 'MapDataFile'
        self.mimeTypes["application/CCP4-map"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-seq"
        mimeType.description = "Sequence file"
        mimeType.fileExtensions = ['pir','seq','fas','fsa','fa','fasta']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'SeqDataFile'
        mimeType.className = 'SeqDataFile'
        mimeType.fixedWidthFont = True
        self.mimeTypes["application/CCP4-seq"] = mimeType

        mimeType = CMimeType()
        mimeType.name = 'application/CCP4-asu-content'
        mimeType.description = "Asu content file"
        mimeType.fileExtensions = ['asu.xml']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'AsuDataFile'
        mimeType.className = 'AsuDataFile'
        mimeType.fixedWidthFont = True
        self.mimeTypes["application/CCP4-asu-content"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-seqalign"
        mimeType.description = "Sequence alignment file"
        mimeType.fileExtensions = ['fas','fasta','pir','aln','msf','phy','bla']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'SeqAlignDataFile'
        mimeType.className = 'SeqAlignDataFile'
        mimeType.fixedWidthFont = True
        self.mimeTypes["application/CCP4-seqalign"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/refmac-dictionary"
        mimeType.description = "Dictionary file"
        mimeType.fileExtensions = ['cif','dict']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'DictDataFile'
        mimeType.className = 'DictDataFile'
        self.mimeTypes["application/refmac-dictionary"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/refmac-external-restraints"
        mimeType.description = "External restraints file"
        mimeType.fileExtensions = ['txt']
        mimeType.fixedWidthFont = True
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'RestraintsFile'
        self.mimeTypes["application/refmac-external-restraints"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/refmac-keywords"
        mimeType.description = "Refmac5 keyword file"
        mimeType.fileExtensions = ['txt']
        mimeType.fixedWidthFont = True
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'RefmacKeywordFile'
        mimeType.className = 'CRefmacKeywordFile'
        self.mimeTypes["application/refmac-keywords"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "text/xml"
        mimeType.description = "XML parameters file"
        mimeType.fileExtensions = ['xml']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'task_file'
        self.mimeTypes["text/xml"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-com"
        mimeType.description = "Program command file"
        mimeType.fileExtensions = ['com']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.icon = 'task_file'
        self.mimeTypes["application/CCP4-com"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-compressed-db"
        mimeType.description = "Compressed database file"
        mimeType.fileExtensions = ['ccp4db.zip','ccp4_project.zip']
        mimeType.icon = 'db_file'
        self.mimeTypes["application/CCP4-compressed-db"] = mimeType

        '''
        mimeType = CMimeType()
        mimeType.name = "application/CCP4-task-data"
        mimeType.description = "Task parameters file"
        mimeType.fileExtensions = ['params.xml','input_params.xml']
        mimeType.viewers = []
        mimeType.icon = 'task_file'
        self.mimeTypes["application/CCP4-task-data"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-task-def"
        mimeType.description = "Task definition file"
        mimeType.fileExtensions = ['def.xml','def']
        mimeType.viewers = []
        mimeType.icon = 'task_file'
        self.mimeTypes["application/CCP4-task-def"] = mimeType
        '''

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-I2-command"
        mimeType.description = "Kludge to get i2 command done"
        mimeType.fileExtensions = ['i2com']
        mimeType.viewers = []
        mimeType.icon = 'task_file'
        self.mimeTypes["application/CCP4-I2-command"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-scene"
        mimeType.description = "CCP4mg scene file"
        mimeType.fileExtensions = ['scene.xml']
        mimeType.viewers = ['ccp4mg']
        mimeType.icon = 'SceneDataFile'
        self.mimeTypes["application/CCP4-scene"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/phaser-sol"
        mimeType.description = "Phaser solution file"
        mimeType.fileExtensions = ['phaser_sol.pkl']
        mimeType.viewers = []
        mimeType.icon = 'PhaserSolDataFile'
        mimeType.className = 'PhaserSolDataFile'
        self.mimeTypes["application/phaser-sol"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "chemical/x-mdl-molfile"
        mimeType.description = "MDL Molfile"
        mimeType.fileExtensions = ['mol']
        mimeType.viewers = ['lidia']
        mimeType.className = 'MDLMolDataFile'
        self.mimeTypes["chemical/x-mdl-molfile"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/iMosflm-xml"
        mimeType.description = "iMosflm data"
        mimeType.fileExtensions = ['imosflm.xml']
        mimeType.viewers = []
        mimeType.className = 'ImosflmXmlDataFile'
        self.mimeTypes["application/iMosflm-xml"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/CCP4-image"
        mimeType.description = "Image file"
        mimeType.fileExtensions = ['img']
        mimeType.viewers = []
        mimeType.className = 'ImageFile'
        self.mimeTypes["application/CCP4-image"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/HHPred-alignments"
        mimeType.description = "HHPred sequence search results"
        mimeType.fileExtensions = ['hhr']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.className = 'HhpredDataFile'
        self.mimeTypes["application/Hhpred-alignments"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/Blast-alignments"
        mimeType.description = "Blast sequence search results"
        mimeType.fileExtensions = ['bla']
        if GRAPHICAL():
            mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.className = 'BlastDataFile'
        self.mimeTypes["application/Blast-alignments"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/dials-jfile"
        mimeType.description = "Dials json data file"
        mimeType.fileExtensions = ['json','expt','jsn']
        mimeType.viewers = [CCP4TextViewer.CTextViewer]
        mimeType.className = 'DialsJsonFile'
        self.mimeTypes["application/dials-jfile"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/dials-pfile"
        mimeType.description = "Dials pickle data file"
        mimeType.fileExtensions = ['pickle','refl']
        mimeType.viewers = []
        mimeType.className = 'DialsPickleFile'
        self.mimeTypes["application/dials-pfile"] = mimeType

        mimeType = CMimeType()
        mimeType.name = "application/phaser-rfile"
        mimeType.description = "Phaser solution file"
        mimeType.fileExtensions = ['phaser_rfile.pkl']
        mimeType.viewers = []
        mimeType.icon = 'PhaserRFileDataFile'
        mimeType.className = 'PhaserRFileDataFile'
        self.mimeTypes["application/phaser-rfile"] = mimeType

    def addMimeType(self,name='',mimeType=None):
        #print 'CCustomMimeTypes.addMimeType',name,mimeType
        if name and mimeType and isinstance(mimeType,CMimeType):
            self.mimeTypes[name] = mimeType
            return 0
        else:
            return 1

    def removeMimeType(self,name=''):
        if name in self.mimeTypes:
            del self.mimeTypes[name]
            return 0
        else:
            return 1

    def getMimeTypes(self):
        return list(self.mimeTypes.keys())

    def getMimeTypeInfo(self,name='',info='',singleExtension=False):
        if name not in self.mimeTypes:
            return ''
        if info == 'filter':
            result = str(getattr(self.mimeTypes[name],'description'))+' ('
            exts = getattr(self.mimeTypes[name],'fileExtensions')
            if singleExtension:
                for e in exts: result = result+'*.'+e.split('.')[1]
            else:
                for e in exts: result = result+'*.'+str(e)+' '
            result = result[0:-1] + ')'
            return result
        else:
            result = getattr(self.mimeTypes[name],info,'')
            if info == 'fileExtensions':
                result0 = []
                for item in result: result0.append(str(item))
                return result0
            else:
                return result

    def getCustomMimeExtensions(self):
        ext_list = []
        for key, mime_type in list(self.mimeTypes.items()):
            ext_list.extend(mime_type.fileExtensions)
        return ext_list

    def isSupportedFormat(self, cformat):
        return cformat in self.mimeTypes

    def classFromFileName(self,fileName=''):
        mimeType = self.formatFromFileExt(fileName=fileName)
        if mimeType is None: return None
        return self.mimeTypes[mimeType].className

    def classListFromFileName(self,fileName=''):
        mimeTypeList = self.formatListFromFileExt(fileName=fileName)
        classList = []
        for mimeType in mimeTypeList:
            if self.mimeTypes[mimeType].className is not None:
                classList.append(self.mimeTypes[mimeType].className)
        return classList

    def getFileExt(self,fileName='',ext='',contentLabel=''):
        ext2 = None
        root = None
        if fileName and not ext:
            root, ext = os.path.splitext(str(fileName))
            if ext:
                ext=ext[1:]
            if ext == 'xml':
                ext2 = os.path.splitext(os.path.splitext(str(fileName))[0])[1]
                if ext2: ext2 = ext2[1:]+'.'+ext
            elif ext == 'mtz':
                from core import CCP4File
                base = os.path.split(root)[1]
                if base.count(CCP4File.CDataFile.SEPARATOR):
                    contentLabel = base.split(CCP4File.CDataFile.SEPARATOR)[-1]
        else:
            ext = ext.strip('.')
        return root,contentLabel,ext,ext2

    def formatFromFileExt(self,fileName='',ext='',contentLabel=''):
        root,contentLabel,ext,ext2 = self.getFileExt(fileName=fileName,ext=ext,contentLabel=contentLabel)
        #print 'CCustomMimeTypes.formatFromFileExt',fileName,root,contentLabel,ext,ext2
        if ext2 is not None:
            for name,mimeType in list(self.mimeTypes.items()):
                if mimeType.fileExtensions.count(ext2):
                    return name
        if root is not None and  ext in ['ent','cif']:
            if ext == 'ent':
                name = self.disambiguateEnt(fileName)
                if name is not None:
                    return name[0]
            elif ext == 'cif':
                name = self.disambiguateCif(fileName)
                if name is not None:
                    return name[0]
        if ext == 'mtz' and len(contentLabel)>0:
            for name in self.miniMtzMimeTypes:
                mimeType = self.mimeTypes[name]
                if mimeType.contentLabel == contentLabel:
                    return name
            return 'application/CCP4-mtz'
        else:
            for name,mimeType in list(self.mimeTypes.items()):
                if mimeType.fileExtensions.count(ext):
                    #print 'CCustomMimeTypes.formatFromFileExt',name
                    return name
        return None

    def formatListFromFileExt(self,fileName='',ext='',contentLabel=''):
        root,contentLabel,ext,ext2 = self.getFileExt(fileName=fileName,ext=ext,contentLabel=contentLabel)
        formatList = []
        #print 'CCustomMimeTypes.formatListFromFileExt',fileName,root,contentLabel,ext,ext2
        if root is not None and  ext in ['ent','cif']:
            path, base = os.path.split(root)
            if ext == 'ent':
                if base[-2:] == 'sf':
                    return ["application/CCP4-mtz-observed"]
                elif base[0:3] == 'pdb':
                    return ["chemical/x-pdb"]
                else:
                    return ["chemical/x-pdb"]
            elif ext == 'cif':
                return self.disambiguateCif(fileName)
        if ext2 is not None:
            for name,mimeType in list(self.mimeTypes.items()):
                if mimeType.fileExtensions.count(ext2):
                    formatList.append(name)
        if ext == 'mtz' and len(contentLabel)>0:
            formatList.append('application/CCP4-mtz')
            for name in self.miniMtzMimeTypes:
                mimeType = self.mimeTypes[name]
                if mimeType.contentLabel == contentLabel:
                    formatList[-1] = name
                    break
        else:
            for name,mimeType in list(self.mimeTypes.items()):
                if mimeType.fileExtensions.count(ext):
                    #print 'CCustomMimeTypes.formatListFromFileExt',name
                    formatList.append(name)
        return formatList

    def disambiguateCif(self,fileName):

        import gemmi
        try:
            tryPdb = gemmi.read_structure(fileName)
            if len(tryPdb) > 0:
                return ["chemical/x-pdb"]
        except RuntimeError:
            try:
                tryDict = gemmi.read_monomer_cif(fileName)
                if len(tryDict.monomers) > 0 or len(tryDict.links) > 0 or len(tryDict.modifications) > 0:
                    return ["chemical/x-cif","application/refmac-dictionary"]
            except:
                return ["text/plain"]
        except:
            return ["text/plain"]
        return ["text/plain"]

    def disambiguateEnt(self, fileName, base=None):
        if base is None:
            base = os.path.splitext(os.path.split(fileName)[1])[0]
        if base[-2:] == 'sf':
            return ["application/CCP4-mtz-observed"]
        elif base[0:3] == 'pdb':
            return ["chemical/x-pdb"]
        else:
            return ["chemical/x-pdb"]

    def getViewers(self, cformat):
        if cformat in self.mimeTypes:
            return getattr(self.mimeTypes[cformat],'viewers',[])
        else:
            return []

    def getMimeTypeForViewer(self, viewer):
        for key,mimetype in list(self.mimeTypes.items()):
            if mimetype.viewers.count(viewer)>0:
                return key
        return ''

    def useDesktopServices(self, cformat):
        print('useDesktopServices',cformat,type(self.mimeTypes))
        if cformat in self.mimeTypes:
            return self.mimeTypes[cformat].useDesktopServices
        else:
            return 0

    def useWebBrowser(self, cformat):
        if cformat in self.mimeTypes:
            return self.mimeTypes[cformat].useWebBrowser
        else:
            return 0

    def getIconsForFileFilters(self):
        # Used by file browser to set icons
        filters = {}
        from core import CCP4Utils
        path = os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons')
        for key, mimetype in list(self.mimeTypes.items()):
            if mimetype.icon is not None:
                iconFile = os.path.join(path, mimetype.icon+'.png')
                if not os.path.exists(iconFile):
                    iconFile = os.path.join(path, mimetype.icon+'.svg')
                if os.path.exists(iconFile):
                    extList = []
                    for item in mimetype.fileExtensions:
                        extList.append('.'+item)
                    filters[iconFile] = extList
        #print 'getIconsForFileFilters',filters
        return filters

    def icon(self,mimeType,modifier=None):
        MODIFIERS = {'import' : ['import_arrow.png',10,12]}
        if mimeType not in self.mimeTypes:
            iconName = 'DataFile'
        else:
            iconName = self.mimeTypes[mimeType].icon
        # It loads with mimetype & key is by mimetype. Changed the if to match properly.
        if not mimeType in CCustomMimeTypes.ICONS:
            fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), 'qticons', iconName+'.png')
            if os.path.exists(fileName):
                CCustomMimeTypes.ICONS[mimeType] = QtGui.QIcon(QtGui.QPixmap(fileName))
        if not mimeType in CCustomMimeTypes.ICONS:
            return QtGui.QIcon()
        if modifier is None or modifier not in MODIFIERS:
            return CCustomMimeTypes.ICONS[mimeType]
        else:
            # Add an overlay to the icon image
            overlayPix = QtGui.QPixmap(os.path.join(CCP4Utils.getCCP4I2Dir(),'qticons',MODIFIERS[modifier][0]))
            #print 'CCustomMimeTypes.icon overlayPix',overlayPix.height(),overlayPix.width()
            #return QtGui.QIcon(overlayPix)
            ICONBUTTONSIZE = 24
            basePix = CCustomMimeTypes.ICONS[mimeType].pixmap(ICONBUTTONSIZE-6,ICONBUTTONSIZE-6)
            pix = QtGui.QPixmap( ICONBUTTONSIZE, ICONBUTTONSIZE)
            painter = QtGui.QPainter( pix )
            # White background before overlayingimages
            painter.fillRect( QtCore.QRect( 0 , 0, ICONBUTTONSIZE, ICONBUTTONSIZE ), QtGui.QColor('white'))
            painter.drawPixmap( basePix.rect(), basePix )
            painter.drawPixmap( MODIFIERS[modifier][1],MODIFIERS[modifier][2],overlayPix.rect().width(),overlayPix.rect().height(), overlayPix )
            ico = QtGui.QIcon(pix)
            # NB if dont del the painter program crashes!
            del painter
            #print 'done modified icon',mimeType,modifier
            return ico

    def fileValidity(self, cformat=None, fileName=None):
        # Bunch of traps for failure to test validity
        if not os.path.exists(str(fileName)):
            return CCP4ErrorHandling.CErrorReport(self.__class__, 103, str(fileName))
        if cformat not in self.mimeTypes:
            return CCP4ErrorHandling.CErrorReport(self.__class__, 101, str(cformat))
        elif self.mimeTypes[cformat].fileValidity is None:
            return CCP4ErrorHandling.CErrorReport(self.__class__, 102, str(cformat))
        try:
            rv = self.mimeTypes[cformat].fileValidity(fileName)
        except:
            return CCP4ErrorHandling.CErrorReport(self.__class__,104)
        # Potentially add warning that file has inapt extension
        if os.path.splitext[1] not in self.mimeTypes[cformat].fileExtensions:
            rv.append(self.__class__,105)
        # Validity function (mimeType.fileValidity) should take fileName as argument and return a CErrorReport
        # The CErrorReport is a list of errors which may be emtpy or contain error code 106, 107
        # (see ERROR_CODES at top of file) or others specific to the file type
        # If creating new error codes consider carefully if they should just have severity set to warning
        # which implies that it would be worth trying to use the file
        # This is all implemented by analogy to the CData.validity() method.
        # Is there a case for a fix() method to deal with broken files?
        # The fileValidity() function could be appended to this file or implemented elsewhere in
        # CCP4i2 core or qtcore directories
        return rv
