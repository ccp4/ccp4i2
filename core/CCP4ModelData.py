from __future__ import print_function

"""
     CCP4ModelData.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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
   Liz Potterton Nov 2010 -  CCP4Data sub-classes relating to xtallographic model
"""

import os
import re
import sys
import types
import shutil
import functools
import gemmi

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from core.CCP4ErrorHandling import *
from core.CCP4Config import QT,XMLPARSER,GRAPHICAL
if QT():
    from core.CCP4QtObject import CObject
else:
    from core.CCP4Object import CObject
if GRAPHICAL():
    from PySide6 import QtCore,QtGui, QtWidgets
else:
    from PySide6 import QtCore
if XMLPARSER() == 'lxml':
    from lxml import etree

from core import CCP4Data
from core import CCP4File
from core import CCP4MathsData


try:
    import Bio.SeqIO, Bio.AlignIO, Bio.Align
    BIOPYTHON = True
except:
    print('FAILED CCP4ModelData imported Bio.SeqIO')
    BIOPYTHON = False

EXTLIST = {'txt' : 'fasta', 'fasta':'fasta', 'fa' :'fasta',
           'fsa' : 'fasta', 'faa' : 'fasta', 'seq' : 'fasta',
           'pir':'pir', 'xml':'seqxml', 'embl':'embl',
           'aln':'clustal', 'clustal' : 'clustal', 'sth':'stockholm',
           'pfam':'stockholm', 'phy':'phylip', 'bla':'blast',
           'hhr' : 'hhpred'}
SEQFORMATLIST = ['pir', 'fasta', 'seqxml', 'embl', 'ncbi']
ALIGNFORMATLIST = ['clustal', 'pir', 'fasta', 'stockholm', 'phylip']
BLASTFORMATLIST = ['blast']
HHPREDFORMATLIST = ['hhpred']

#  PIR format description from http://www.bioinformatics.nl/tools/crab_pir.html

PIR_DESCRIPTION = """A  PIR (sometimes called NBRF) file contains one or more sequences in the form:
        Line 1:  ">P1;" which includes a two-letter code defining the sequence type (P1, F1, DL, DC, RL, RC, or XX)
           followed by the database ID code.
        Line 2: text description of the sequence.
        Lines 3+:  the sequence, which can include white space and '-' (that will be ignored) ending with "*" character.
        Optionally these can be followed by more lines describing the sequence.

    Example:

>P1; RNASE
Chain A for Rnase
DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDGVVFQNRESVLPTQSYGYYHEYTVITPGARTRGTRRIICGE
ATQEDYYTGDHYATFSLIDQTC*
    """


def MODELINFO(*keys):
    # Extract data from selection_protocols.py - this is a kludge
    if not CModelInfo.insts:
        CModelInfo()
    return CModelInfo.insts.info(keys)

class CModelInfo:
    insts = None

    def __init__(self):
        from core import CCP4Utils
        CModelInfo.insts = self
        globalVar = {}
        localVar = {}
        fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), 'data', 'model_description', 'selection_protocols.py')
        exec(compile(open(fileName).read(), fileName, 'exec'),globalVar,localVar)
        self.testDict = localVar

    def info(self, keys=[]):
        if keys[0] in self.testDict:
            subD = self.testDict[keys[0]]
        else:
            return None
        for key in keys[1:]:
            if key in subD:
                subD = subD[key]
            else:
                subD = None
            break
        return subD

class CBioPythonSeqInterface:
    '''Interface to BioPython'''
    ERROR_CODES = {401 : {'description' : 'Attempting to load from non-existent file'},
                   402 : {'description' : 'Error reading from file'},
                   403 : {'description' : 'Unknown sequence file format'},
                   405 : {'description' : 'Error reading identifiers from multi-record file'},
                   406 : {'description' : 'Error opening file'},
                   407 : {'description' : "The 'PIR' file did not have the correct format"},
                   408 : {'severity' : SEVERITY_WARNING, 'description' : "The 'PIR' file format was corrected"},
                   409 : {'description' : "Error opening file to write"},
                   410 : {'description' : "Error attempting to write out sequence file"},
                   411 : {'description' : "Error attempting to create a temporary sequence file"},
                   412 : {'description' : "Sequence file is empty"},
                   413 : {'description' : "Unable to read BLAST format file"},
                   414 : {'description' : "Unable to read hhpred format file"}}

    def __init__(self):
        pass

    def simpleFormatTest(self, filename):
        from core import CCP4Utils
        formt = 'unknown'
        text = CCP4Utils.readFile(str(filename))
        segments = text.split('>')
        if len(segments[0]) != 0:
            return formt
        lines = ('>'+segments[1]).split('\n')
        nonAlphaList = []
        nonAlphaTot = 0
        il = 0
        for l in lines:
#If only non-alpha character is a "-" and it is not the first line we assume this is a gap signifier and OK.
            if il >0  and ((len(list(set(re.findall('[^(a-z,A-Z, )]', l)))) == 1 and list(set(re.findall('[^(a-z,A-Z, )]', l)))[0] == "-") or l.endswith("*")):
                nonAlphaList.append(0)
            else:
                nonAlphaList.append(len(re.findall('[^(a-z,A-Z, )]', l)))
            nonAlphaTot = nonAlphaTot + nonAlphaList[-1]
            il += 1
        #print 'simpleFormatTest nonAlphaList',nonAlphaList,nonAlphaTot
        if nonAlphaTot == 0:
            formt = 'noformat'
        elif len(lines) > 0 and len(lines[0]) > 0 and '>' ==  lines[0][0]:
            if ';' in lines[0] and (nonAlphaTot - (nonAlphaList[0] + nonAlphaList[1] + nonAlphaList[-1]))==0:
                formt = 'pir'
            elif (nonAlphaTot - nonAlphaList[0]) == 0:
                formt = 'fasta'
        return formt

    def loadExternalFile(self,filename=None, format=None, record=0, diagnostic=False):
        #print 'CBioPythonSeqInterface.loadExternalFile filename', filename, format, record
        from core import CCP4Utils
        filename = str(filename)
        if not os.path.exists(filename):
            return
        #import traceback
        #traceback.print_stack(limit=10)
        # Unset now so if file reading fails at least we dont have misleading data
        for item in self.CONTENTS_ORDER:
            self.__dict__['_value'][item].unSet()
        self.__dict__['lastLoadedFormat'] = None
        if 'loadWarning' in self.__dict__:
            del  self.__dict__['loadWarning']
        # Beware CDataObj input
        self.blockSignals(True)
        if self.__dict__.get('lastLoadedFile',None) is not None and not os.path.exists(self.__dict__['lastLoadedFile']):
            self.unSet()
            raise CException(CBioPythonSeqInterface, 401, str(filename), name=self.objectPath())
        else:
            if BIOPYTHON:
                # guess file format from extension
                # Formats listed at http://biopython.org/wiki/SeqIO (and AlignIO)
                if format is  None or format == 'unknown':
                    ext = os.path.splitext(str(filename))[1][1:]
                    format = EXTLIST.get(ext, None)
                # Put expected formated first in testOrder
                testOrder = SEQFORMATLIST + ALIGNFORMATLIST
                if format is not None:
                    testOrder.remove(format)
                    testOrder.insert(0, format)
                #print 'CBioPythonSeqInterface.loadExternalFile', format, testOrder
                n = 0
                while n < len(testOrder):
                    if diagnostic:
                        print('CBioPythonSeqInterface.loadExternalFile trying', testOrder[n], 'record',record)
                    err,rv = self.bioLoadSeqFile(filename, testOrder[n], record=record)
                    if diagnostic:
                        print('CBioPythonSeqInterface.loadExternalFile', testOrder[n], err, rv)
                    if err.maxSeverity() <= SEVERITY_WARNING:
                        break
                    if n == 0 and format == 'pir':
                        fixedPirFile,err0,rv0 = self.fixPirFile(filename)
                        #print 'from fixPirFile', fixedPirFile, err0, rv0
                        if fixedPirFile is None:
                            if err0 is not None:
                                self.__dict__['loadWarning'] = err0
                            else:
                                self.__dict__['loadWarning'] = CErrorReport(self.__class__, 407, filename, stack=False)
                        else:
                            self.__dict__['loadWarning'] = CErrorReport(self.__class__, 408, fixedPirFile, stack=False)
                            self.__dict__['validatedFile'] = fixedPirFile
                            err = err0
                            rv  = rv0
                            break
                    n = n + 1
                if n < len(testOrder):
                    if diagnostic:
                        print('CBioPythonSeqInterface.loadExternalFile Sequence file read as format:', testOrder[n], rv)
                    self.__dict__['lastLoadedFormat'] = testOrder[n]
                    self.__dict__['lastLoadedFile'] = str(filename)
                    if testOrder[n] == 'fasta':
                        try:
                            referenceDb, reference, identifier = rv['identifier'].split('|')
                        except:
                            pass
                        else:
                            # .. and beware biopython sets identifier to  first word of first line but we want to allow spaces
                            # in the identifier .. so reread
                            try:
                                text = CCP4Utils.readFile(fileName=filename)
                                lines = text.split('\n')
                                db, ref, identifier = lines[0][1:].split('|')
                            except:
                                print('Failed rereading identifier for',filename)
                            if referenceDb in self.referenceDb.qualifiers('enumerators'):
                                rv['referenceDb'] = referenceDb
                            elif referenceDb in self.referenceDb.qualifiers('menuText'):
                                rv['referenceDb'] = self.referenceDb.qualifiers('enumerators')[self.referenceDb.qualifiers('menuText').index(referenceDb)]
                            elif referenceDb.lower().count('swi'):
                                rv['referenceDb'] = 'sp'
                            elif referenceDb.lower().count('uni'):
                                rv['referenceDb'] = 'tr'
                            reference = reference.strip()
                            if len(reference) > 0:
                                rv['reference'] = reference
                            # Dont change the identifier - its becoming useless too often
                            #rv['identifier'] = identifier
                    if diagnostic:
                        print('CBioPythonSeqInterface.loadExternalFile rv', filename, rv)
                    for item in ['identifier', 'sequence', 'reference', 'referenceDb', 'name', 'description']:
                        if item in rv and item in self.__dict__['_value']:
                            self.__dict__['_value'][item].set(rv[item])
                    self.blockSignals(False)
                    self.dataChanged.emit()
                return
            # Try to do it without biopython
            #try:
            if 1:
                text = CCP4Utils.readFile(str(filename))
                lines = text.split('\n')
                if '>P1;' in lines[0]:
                    format = 'pir'
                    try:
                        for idx in range(2, len(lines)):
                            if lines[idx][-1] == '*':
                                lines[idx] = lines[idx][0:-1]
                    except:
                        pass
                    l0 = 2
                elif '>' in lines[0]:
                    format = 'fasta'
                    l0 = 1
                else:
                    format = 'adhoc'
                    l0 = 0
                seq = lines[l0]
                for line in lines[l0 + 1:]:
                    seq = seq + '\n' + line
                if 'sequence' in self.__dict__['_value']:
                    self.__dict__['_value']['sequence'].set(seq)
                print('Loaded sequence without BioPython assuming format:', format)
            #except:
            #  raise CException(self.__class__, 402, str(filename), name=self.objectPath())
            self.__dict__['lastLoadedFile'] = str(filename)
            self.__dict__['lastLoadedFormat'] = format
            self.__dict__['_value']['identifier'].set(os.path.splitext(os.path.basename(filename))[0])
            self.blockSignals(False)
            self.dataChanged.emit()

    def fixPirFile(self, fileName, importedFile=None):
        import tempfile
        from core import CCP4Utils
        try:
            text = CCP4Utils.readFile(str(fileName))
        except:
            return None, None, None
        if len(text.strip()) == 0:
            return None, CErrorReport(self.__class__, 412, stack=False), None
        fragments = text.split('\n>')
        #print 'fixPirFile fragments',fragments
        if len(fragments) > 0 and len(fragments[0]) > 0 and fragments[0][0] == '>':
            fragments[0] = fragments[0][1:]
        output = ''
        for text in fragments:
            text = text.strip()
            if len(text) > 2 and text[2] != ';':
                text = 'P1;' + text
            if not text.count('*'):
                text = text + '*'
            #print 'fixPirFile fragmant',text
            output = output + '>' + text + '\n'
        # Write to temporary file
        f1 = tempfile.mkstemp()
        if sys.version_info > (3,0):
            os.write(f1[0], output.encode("utf-8"))
        else:
            os.write(f1[0], output)
        os.close(f1[0])
        err, data = self.bioLoadSeqFile(f1[1], 'pir')
        #print 'fixPirFile',f1[1],err,data
        if err.maxSeverity() <= SEVERITY_WARNING:
            if importedFile is not None:
                shutil.move(f1[1], importedFile)
                return importedFile, err, data
            else:
                return f1[1], err, data
        else:
            return None, None, None

    def fixFastaFile(self, fileName, importedFile=None):
        # Attempting to fix a 'seq' file containing only sequence & no formatting
        import tempfile
        from core import CCP4Utils
        try:
            text = CCP4Utils.readFile(str(fileName))
        except:
            return None, None, None
        if len(text.strip()) == 0:
            return None, CErrorReport(self.__class__, 412, stack=False), None
        lines = text.split('\n')
        nonAlphaList = []
        for l in lines:
            nonAlphaList.extend(re.findall('[^(a-z,A-Z, )]', l))
        if len(nonAlphaList) > 0:
            return None, None, None
        output = '>' + os.path.splitext(os.path.basename(fileName))[0] + '\n' + text
        # Write to temporary file
        f1 = tempfile.mkstemp()
        if sys.version_info > (3,0):
            os.write(f1[0], output.encode("utf-8"))
        else:
            os.write(f1[0], output)
        os.close(f1[0])
        err,data = self.bioLoadSeqFile(f1[1], 'fasta')
        #print 'fixFastaFile',f1[1],err,data
        if err.maxSeverity() <= SEVERITY_WARNING:
            if importedFile is not None:
                shutil.move(f1[1], importedFile)
                return importedFile, err, data
            else:
                return f1[1], err, data
        else:
            return None, None, None

    def bioLoadSeqFile(self, filename, format, record=0, saveFormat=None, saveFilename=None, diagnostic=False):
        err = CErrorReport()
        idList = []
        try:
            f = open(str(filename), 'r')
        except:
            err.append(self.__class__, 406, filename, stack=False)
            return err, {}
        # Set up to save (probably either chosing one record or changing format)
        fout = None
        if saveFormat is not None and saveFormat in SEQFORMATLIST:
            if saveFilename is None:
                try:
                    import tempfile
                    fout, saveFilename = tempfile.mkstemp()
                except:
                    err.append(CBioPythonSeqInterface, 411, stack=False)
            else:
                try:
                    fout = open(str(saveFilename), 'w')
                except:
                    err.append(CBioPythonSeqInterface, 409, saveFilename, stack=False)
        try:
            if format in SEQFORMATLIST:
                seq_records = list(Bio.SeqIO.parse(f, format))
            elif format in ALIGNFORMATLIST:
                ali_records = list(Bio.AlignIO.parse(f, format))
                seq_records = ali_records[0]
        except Exception as e:
            err.append(CBioPythonSeqInterface, 402, saveFilename, exc_info=sys.exc_info())
            return err, {}
        else:
            for rec in seq_records:
                idList.append(rec.id)
        if format in BLASTFORMATLIST:
            err.append(CBioPythonSeqInterface, self.__class__, 413, stack=False)
            return err, {}
        if format in HHPREDFORMATLIST:
            err.append(CBioPythonSeqInterface, self.__class__, 414, stack=False)
            return err, {}
        '''
        if format in BLASTFORMATLIST:
          from Bio.Blast import NCBIStandalone
          titleList = []
          queryList = []
          sbjctList = []
          blast_parser = NCBIStandalone.BlastParser()
          blast_record = blast_parser.parse(f)
          for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
              titleList.append(alignment.__str__().split('\n')[0])
              queryList.append(hsp.query)
              sbjctList.append(hsp.sbjct)
        '''
        """
        if len(seq_records)>0:
          print 'sequence record dir',dir(seq_records[0])
          for item in ['annotations', 'dbxrefs', 'description', 'features', 'format', 'id', 'letter_annotations', 'lower', 'name', 'reverse_complement', 'seq', 'upper']:
            print item,seq_records[0].__getattribute__(item)
        """
        if fout is not None:
            try:
                nRecOut = Bio.SeqIO.write(seq_records[record], fout, saveFormat)
            except:
                err.append(CBioPythonSeqInterface, 410, filename, exc_info=sys.exc_info())
        try:
            #print 'bioLoadSeqFile',seq_records[record].seq
            return err, {'format' : format, 'identifier': seq_records[record].id, 'sequence': seq_records[record].seq,
                         'idList' : idList , 'name' : seq_records[record].name, 'description': seq_records[record].description}
        except:
            err.append(CBioPythonSeqInterface, 402, filename, stack=False)
            return err, {}


class CSequenceMeta(CCP4Data.CData):

    CONTENTS = {'uniprotId' : {'class' : CCP4Data.CString},
                'organism' : {'class' : CCP4Data.CString},
                'expressionSystem' :{'class' : CCP4Data.CString}}
    ERROR_CODES = {401 : {'description' : 'No uniprot id available'},
                   402 : {'description' : 'No uniprot xml file available to read'},
                   403 : {'description' : 'No project id provided to determine uniprot xml filename'},
                   404 : {'description' : 'Reading uniprot xml file failed'}}

    def getUniprotUrl(self):
        if not self.uniprotId.isSet():
            return None
        else:
            code = str(self.uniprotId)
            return "http://www.uniprot.org/uniprot/" + str(self.uniprotId)

    def getUniprotXml(self, projectId=None):
        from core import CCP4Modules
        if projectId is None:
            raise CErrorReport(self.__class__, 403)
        if not self.uniprotId.isSet():
            raise CErrorReport(self.__class__, 401)
        tmpDir = os.path.join(CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectId=projectId),'CCP4_DOWNLOADED_FILES')
        if not os.path.exists(tmpDir):
            try:
                os.mkdir(tmpDir)
            except:
                from core import CCP4Utils
                tmpDir = CCP4Utils.getTMP()
        targetFile = os.path.join(tmpDir, 'uniprot_' + str(self.uniprotId) + '.xml')
        return targetFile

    def downloadUniprotXml(self, projectId):
        #can test download with something like:
        #from core.CCP4ModelData import *; s = CSequenceMeta(uniprotId='P12345'); s.loadFromUniprotXml(projectId='efaebd47a70a11e493bf9cf387d93af8')
        from qtgui import CCP4FileBrowser
        import ccp4mg
        import UtilityThread
        if not self.uniprotId.isSet():
            raise CErrorReport(self.__class__, 401)
        mode = 'uniprotXml'
        code = str(self.uniprotId)
        urlname = "http://www.uniprot.org/uniprot/" + code + ".xml"
        targetFile = self.getUniprotXml(projectId=projectId)
        if 'downloader' not in self.__dict__:
            self.__dict__['downloader'] =  CCP4FileBrowser.CDownloader()
        self.__dict__['downloader'].Finished.connect(functools.partial(self.handleDownloadFinished,code,mode,targetFile))
        self.__dict__['downloader'].Error.connect(functools.partial(self.handleDownloadError,code))
        self.__dict__['downloadThread'] = UtilityThread.UtilityThread(functools.partial(self.__dict__['downloader'].download, urlname))
        self.__dict__['downloadThread'].start()

    @QtCore.Slot(str,str,str,str)
    def handleDownloadFinished(self, code, mode, targetFile, tmpFile):
        #print 'handleDownloadFinished',code,mode,targetFile,tmpFile
        if os.path.exists(tmpFile):
            shutil.copyfile(tmpFile, targetFile)

    @QtCore.Slot(str)
    def handleDownloadError(self, code):
        print('handleDownloadError', code)

    def loadFromUniprotXml(self, fileName=None, projectId=None):
        from lxml import etree
        from core import CCP4Utils
        if fileName is None:
            fileName = self.getUniprotXml(projectId=projectId)
        if fileName is None or not os.path.exists(fileName):
            raise CErrorReport(self.__class__, 402, str(fileName))
        # Every element in these files has the namespace prefix
        nsmap = {'u' : 'http://uniprot.org/uniprot'}
        ret = {}
        root = CCP4Utils.openFileToEtree(fileName)
        try:
            eleList = root.xpath('./u:entry', namespaces=nsmap)
            if len(eleList) > 0:
                ret['db'] = eleList[0].get('dataset')
            eleList = root.xpath('./u:entry/u:accession', namespaces=nsmap)
            if len(eleList) > 0:
                ret['accessionList'] = []
                for ele in eleList:
                    ret['accessionList'].append(str(ele.text))
        except:
            raise CErrorReport(self.__class__, 404, str(fileName))
        eleList = root.xpath('./u:entry/u:protein', namespaces=nsmap)
        if len(eleList) > 0:
            e = eleList[0].xpath('./u:recommendedName/u:fullName', namespaces=nsmap)
            print('e', e)
            if len(e) > 0: 
                ret['proteinFullName'] = str(e[0].text)
        eleList = root.xpath('./u:entry/u:sequence', namespaces=nsmap)
        if len(eleList) > 0:
            ret['sequence'] = str(eleList[0].text)
        return ret

class CSequence(CCP4Data.CData, CBioPythonSeqInterface):
    '''
    A string of sequence one-letter codes
    Need to be able to parse common seq file formats
    Do we need to support alternative residues
    What about nucleic/polysach?
    '''
    ERROR_CODES = {201 : {'description' : 'Sequence undefined', 'severity' : SEVERITY_UNDEFINED},
                   202 : {'description' : 'error reading from file'},
                   203 : {'description' : 'Comparing sequences: Sequence item different'},
                   204 : {'description' : 'Comparing sequences: One item set - the other is unset'}}
    ERROR_CODES.update(CBioPythonSeqInterface.ERROR_CODES)
    CONTENTS = {'identifier' : {'class' : CCP4Data.CString,
                                'qualifiers' : {'toolTip' : 'Description of sequence' , 'minlength' : 4}} ,
                'referenceDb' : {'class' : CCP4Data.CString ,
                                 'qualifiers' : {'onlyEnumerators' : False,
                                                 'default' : 'unk',
                                                 'enumerators' : ['unk', 'sp' , 'tr' , 'pdb'],
                                                 'menuText' : ['Unknown' , 'UniProt/Swiss-Prot', 'UniProt/TrEMBL', 'ProteinDatabank']}},
                'reference' : {'class' : CCP4Data.CString,
                               'qualifiers' : {'toolTip' : 'Optional reference for sequence'}},
                'name' : {'class' : CCP4Data.CString,
                          'qualifiers' : {'toolTip' : 'User friendly name of sequence'}},
                'description' : {'class' : CCP4Data.CString,
                                 'qualifiers' : {'toolTip' : 'User friendly description of sequence'}},
                'sequence' : {'class' : CCP4Data.CString,
                              'qualifiers' : {'toolTip' : 'Single letter sequence (white space and dash ignored)'}},
                'moleculeType' : {'class' : CCP4Data.CString,
                                  'qualifiers' : {'onlyEnumerators' : True, 'enumerators' : ['PROTEIN','NUCLEIC'],
                                                  'menuText' : ['protein','nucleic acid'], 'default' : 'PROTEIN',
                                                  'toolTip' : 'Molecule type'}}}
    CONTENTS_ORDER = ['identifier', 'name', 'description', 'referenceDb', 'reference', 'moleculeType', 'sequence']

    def saveFile(self, fileName=None):
        from core import CCP4Utils
        '''
        if self.reference.isSet():
          ref = self.reference.__str__().strip()
        else:
          ref = ' '
        text = '>'+str(self.referenceDb)+ '|' +ref + '|' + self.identifier.__str__().strip() + \
              '\n'+self.sequence.__str__()
        '''
        text = '>'+self.identifier.__str__().strip() + '\n' + self.sequence.__str__()
        CCP4Utils.saveFile(fileName=fileName, text=text)

    def loadFile(self, fileName, format='unknown'):
        #print 'CSequence.loadFile',fileName,format
        #import traceback
        #traceback.print_stack(limit=5)
        if format == 'internal':
            self.loadInternalFile(fileName)
        elif format == 'uniprot':
            #Try treating it as a uniprot file
            m = CSequenceMeta()
            try:
                ret = m.loadFromUniprotXml(fileName)
            except:
                pass
            else:
                self.sequence = ret['sequence']
                if ret['db'] == 'Swiss-Prot':
                    self.referenceDb = 'sp'
                if len(ret['accessionList'][0]) > 0:
                    self.reference = ret['accessionList'][0]
                self.identifier = ret['proteinFullName']
                self.__dict__['lastLoadedFile'] = str(fileName)
                self.__dict__['lastLoadedFormat'] = 'uniprot'
                return
        else:
            self.unSet()
            CBioPythonSeqInterface.loadExternalFile(self, fileName, format=format)

    def loadInternalFile(self, fileName):
        self.unSet()
        from core import CCP4Utils
        try:
            text = CCP4Utils.readFile(fileName=fileName)
        except Exception as e:
            print('ERROR loading sequence file', e)
            return
        lines = text.split('\n')
        try:
            splitList = lines[0][1:].split('|')
            if len(splitList) == 3:
                self.referenceDb.set(splitList[0])
                if len(splitList[1].strip()) > 0:
                    self.reference.set(splitList[1])
            self.identifier.set(lines[0][1:])
            seq = lines[1]
            for l in lines[2:]:
                seq = seq + l
            self.sequence.set(seq)
            #print 'CSequence.loadInternalFile', fileName, seq
        except Exception as e:
            print('ERROR loading sequence file', e)

    def getAnalysis(self, mode='molecularWeight'):
        if mode == 'molecularWeight':
            if not self.__dict__['_value']['sequence'].isSet():
                return 0
            from Bio.SeqUtils.ProtParam import ProteinAnalysis
            # Beware BioPython not robust to bad sequences
            seq = re.sub('[^GALMFWKQESPVICYHRNDT]', '', str(self.__dict__['_value']['sequence']))
            pa = ProteinAnalysis(seq)
            print('CSequence.getAnalysis', str(self.__dict__['_value']['sequence']))
            print('CSequence.getAnalysis', pa.molecular_weight())
            return pa.molecular_weight()

    def assertSame(self, other):
        err = CErrorReport()
        for item in ['identifier', 'reference', 'referenceDb', 'sequence']:
            if self.__dict__['_value'][item].isSet() != other.__dict__['_value'][item].isSet() and item != 'reference':
                #print 'CSequence.assertSame',item,self.__dict__['_value'][item].isSet(),other.__dict__['_value'][item].isSet()
                err.append(self.__class__, 204, 'For ' + item + ': *' + str(self.__dict__['_value'][item]) + '*  *' + str(other.__dict__['_value'][item]) + '*', stack=False)
            elif item == 'sequence':
                mySeq = re.sub('\n', '', self.__dict__['_value'][item].__str__())
                otherSeq = re.sub('\n', '', other.__dict__['_value'][item].__str__())
                if mySeq != otherSeq:
                    err.append(self.__class__, 203, 'For ' + item + ': *' + str(self.__dict__['_value'][item]) + '*  *' + str(other.__dict__['_value'][item]) + '*' , stack=False)
            elif len(self.__dict__['_value'][item]) > 0 and len( other.__dict__['_value'][item]) > 0 and \
                     self.__dict__['_value'][item] != other.__dict__['_value'][item]:
                err.append(self.__class__, 203, 'For ' + item + ': *' + str(self.__dict__['_value'][item]) + '*  *' + str(other.__dict__['_value'][item]) + '*' , stack=False)
        #print 'CSequence.assertSame',err.report()
        return err

    def guiLabel(self):
        if self.identifier.isSet():
            return str(self.identifier)
        #elif self.reference.isSet():
        #  return str(self.reference)
        elif self.sequence.isSet():
            return str(self.sequence)[0:20]
        else:
            return str(self.objectName())

    def getTextItem(self):
        return self.guiLabel()

    def validity(self, args):
        err = CErrorReport()
        #print 'CSequence.validity args', args
        if ('sequence' not in args) or args['sequence'] is None:
            err.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif  isinstance(args['sequence'],CCP4Data.CString):
            if not args['sequence'].isSet():
                err.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        elif len(args['sequence']) == 0:
            err.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err

    def __eq__(self, arg):
        # Considered equal if name same identifier
        return self.__dict__['_value']['identifier'].__eq__(arg.get('identifier'))


class CChemComp(CCP4Data.CData):
    '''Component of CDictDataFile contents'''
    CONTENTS = {'id' : {'class' : CCP4Data.COneWord}, 'three_letter_code' : {'class' : CCP4Data.COneWord},
                'name' : {'class' : CCP4Data.CString}, 'group' : {'class' : CCP4Data.CString},
                'number_atoms_all' : {'class' : CCP4Data.CInt}, 'number_atoms_nh' : {'class' : CCP4Data.CInt},
                'desc_level' : {'class' : CCP4Data.CInt},}
    ERROR_CODES = {201 : {'description' : 'Error reading monomer id and name'},
                   202 : {'description' : 'Error writing monomer id and name'}}

    def load(self, loop, loopIndex):
        # Beware GetReal/GetString needed swig hack to work and have
        # returned parameters order reversed
        errCount = 0
        for key in ['id', 'three_letter_code', 'name', 'group']:
            try:
                ret = loop.GetString(key, loopIndex)
                if isinstance(ret, list):
                    value, rv = ret
                else:
                    rv = ret
                #print 'CChemComp.load',loopIndex,key,value,rv
                self.__dict__['_value'][key] = value.strip()
            except:
                #print 'CChemComp.load error reading',key
                errCount += 1
        for key in ['number_atoms_all', 'number_atoms_nh']:
            try:
                value, rc = loop.GetInteger(key, loopIndex)
                #print 'CChemComp.load', loopIndex, key, value, rc
                if rc == 0:
                    self.__dict__['_value'][key] = value
            except:
                pass
        if self.__dict__['_value']['id'] is None or self.__dict__['_value']['id'] == '.':
            self.__dict__['_value']['id'] =  self.__dict__['_value']['three_letter_code']
        if errCount > 0:
            return CException(self.__class__, 201, name=self.objectPath())
        else:
            return CException()

    def writeToLoop(self, loop):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        errCount = 0
        indx = loop.GetLoopLength()
        for key in ['id', 'three_letter_code', 'name', 'group']:
            try:
                if self.__dict__['_value'][key].isSet():
                    loop.PutString(self.__dict__['_value'][key].__str__(), key, indx)
                else:
                    loop.PutString('.', key, indx)
            except:
                errCount += 1
        for key in ['number_atoms_all', 'number_atoms_nh']:
            try:
                if self.__dict__['_value'][key].isSet():
                    loop.PutInteger(self.__dict__['_value'][key].__int__(), key, indx)
            except:
                errCount += 1
        try:
            loop.PutNoDataType(mmdb.CIF_NODATA_DOT, 'desc_level', indx)
        except:
            pass
        if errCount > 0:
            return CException(self.__class__, 202, name=self.objectPath())
        else:
            return CException()


class CDictData(CCP4Data.CData):
    # content is a list of the chem_comp
    CONTENTS = {'monomerList' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CChemComp}}}

    ERROR_CODES = {101 : {'description' : 'Error opening MMCIF format file'},
                   102 : {'description' : 'Error merging data - monomer already in geometry file'},
                   103 : {'severity' : SEVERITY_WARNING, 'description' : 'Warning merging data - overwriting geometry for monomer with same id'},
                   104 : {'description' : 'Error reading geometry cif file - does not contain expected data'},
                   105 : {'description' : 'Unknown error reading geometry file'},
                   106 : {'description' : '_chem_comp section not found in geometry file'},
                   110 : {'description' : 'Attemting to delete unrecognised chem_comp.id'}}

    def openCifFile(self, fileName=None):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        print('#CDictData.openCifFile fileName', fileName)
        cifFile = mmdb.File()
        print('#CDictData.openCifFile cifFile', cifFile)
        rc = cifFile.ReadMMCIFFile(fileName)
        print('#CDictData.openCifFile rc', rc)
        if rc !=0:
            raise CException(self.__class__, 101, 'Code: ' + str(rc) + ' file: ' + str(fileName), name=self.objectPath())
        return cifFile

    def getCompData(self, id):
        cifFile = self.openCifFile(self.parent().__str__())
        dataBlock = cifFile.GetCIFData('comp_' + id)
        return dataBlock

    def getCompDataText(self, id=None):
        return None

    def monomerIdList(self):
        idList = []
        for obj in self.__dict__['_value']['monomerList']:
            idList.append(obj.id.__str__())
        return idList

    def loadFile(self, fileName=None):
        # Beware fileName is CFilePath
        #print 'CDictData.loadFile fileName', fileName
        if fileName is None:
            self.__dict__['_value']['monomerList'].unSet()
            return
        err = CException()
        cifFile=self.openCifFile(str(fileName))
        #print 'CDictData.loadFile cifFile', cifFile
        dataBlock = cifFile.GetCIFData('comp_list')
        #print 'CDictData.loadFile dataBlock', dataBlock
        if dataBlock is None:
            raise CException(self.__class__, 104, stack=False)
        self.__dict__['_value']['monomerList'].unSet()
        print('CDictData.loadFile monomerList', self.__dict__['_value'])
        print(dir(dataBlock))
        loop = dataBlock.GetLoop('_chem_comp')
        #print 'CDictData.loadFile loop', loop
        if loop is not None:
            loopLength = loop.GetLoopLength()
            #print 'CDictData.loadFile loopLength', loopLength
            for n in range(loopLength):
                if loop.GetString('id', n) != 0:
                    chemComp = self.__dict__['_value']['monomerList'].addItem()
                    err.extend(chemComp.load(loop=loop, loopIndex=n))
            if err.maxSeverity() > SEVERITY_WARNING:
                raise err
        #print 'CDictData.loadFile', self.monomerList
        self.dataChanged.emit()

    def getChemComp(self, id=None):
        indx = 0
        for obj in self.__dict__['_value']['monomerList']:
            if obj.id == id:
                return obj, indx
            indx += 1
        return None,-1

    def mergeFile(self, fileName, overwrite=False):
        #print 'CDictData.mergeFile',fileName
        err = CException()
        dictObj = CDictDataFile(fullPath=fileName)
        try:
            dictObj.loadFile()
        except CException as e:
            err.extend(e)
        except Exception as e:
            err.append(self.__class__, 105, exc_info=sys.exc_info())
        #print 'CDictData.mergeFile', err.report()
        if err.maxSeverity() > SEVERITY_WARNING:
            return err
        err.extend(self.merge(dictData=dictObj.fileContent, overwrite=overwrite))
        return err

    def merge(self, dictData=None, idList=None, overwrite=False):
        #print 'CDictData.merge', dictData, idList
        err = CException()
        otherCifFile = dictData.openCifFile(dictData.parent().__str__())
        cifFile = self.openCifFile(self.parent().__str__())
        if idList is None:
            idList = dictData.monomerIdList()
        compLoop = cifFile.GetCIFData('comp_list').GetLoop('_chem_comp')
        #print 'CDictData.merge', overwrite, compLoop
        for idd in idList:
            doMerge = True
            if self.getChemComp(idd)[0] is not None:
                if not overwrite:
                    doMerge = False
                    err.append(self.__class__, 102, idd, name=self.objectPath())
                else:
                    err.append(self.__class__, 103, idd, name=self.objectPath())
                    self.delete(id)
            #print 'merge',idd,doMerge
            if doMerge:
                chemCompObj, chemCompIndex = dictData.getChemComp(idd)
                if chemCompObj is not None:
                    newObj = self.__dict__['_value']['monomerList'].addItem()
                    newObj.set(chemCompObj)
                    newObj.writeToLoop(compLoop)
                otherDataBlock = otherCifFile.GetCIFData('comp_' + idd)
                #print 'merge otherDataBlock',otherDataBlock
                if otherDataBlock is not None:
                    rc = cifFile.AddCIFData('comp_' + idd)
                    myDataBlock = cifFile.GetCIFData('comp_' + idd)
                    #print 'copy',rc,myDataBlock
                    myDataBlock.Copy(otherDataBlock)
            err.extend(self.save(cifFile=cifFile))
        return err

    def delete(self, id):
        cifFile = self.openCifFile(self.parent().__str__())
        # get the monomer row - and delete from monomerList
        chemCompObj, chemCompIndex = self.getChemComp(id)
        if chemCompObj is None:
            return CException(self.__class__, 110, id, name=self.objectPath())
        self.__dict__['_value']['monomerList'].__delitem__(chemCompIndex)
        # Del the monomer from list
        compLoop = cifFile.GetCIFData('comp_list').GetLoop('_chem_comp')
        if compLoop.GetLoopLength() == 1:
            # Don't delete the last row - does not write out the loop and screws up subsequent reading 
            for key in ['id', 'three_letter_code', 'name','group']:
                compLoop.PutString('.', key, 1)
        else:
            compLoop.DeleteRow(chemCompIndex)
        compLoop.Optimize()
        # Del the monomer block
        dataBlock = cifFile.GetCIFData('comp_' + id)
        if dataBlock is not None:
            for lname in ['_chem_comp_atom', '_chem_comp_tree', '_chem_comp_bond', '_chem_comp_angle', '_chem_comp_tor', '_chem_comp_chir', '_chem_comp_plane_atom']:
                dataBlock.DeleteLoop(lname)
                dataBlock.Optimize()
        return self.save(cifFile=cifFile)

    def save(self, cifFile=None, fileName=None):
        from core import CCP4Utils
        if fileName is None:
            fileName = self.parent().__str__()
        try:
            CCP4Utils.backupFile(fileName=fileName, delete=False)
            cifFile.WriteMMCIFFile(fileName)
        except:
            self.loadFile(fileName=fileName)
            return CException(self.__class__, 104, fileName, name=self.objectPath())
        else:
            #print 'CDictData.save emitSignal dataChanged'
            self.dataChanged.emit()
            return CException()

class CDictDataFile(CCP4File.CDataFile):
    '''A refmac dictionary file'''

    QUALIFIERS = {'fileLabel' : 'dictionary', 'mimeTypeName' : 'application/refmac-dictionary',
                  'mimeTypeDescription' : 'Geometry file', 'guiLabel' : 'Geometry dictionary',
                  'toolTip' : 'Idealised geometry of ligands for refinement',
                  'fileExtensions' : ['cif'], 'fileContentClassName' : 'CDictData',
                  'helpFile': 'model_data#ligand_geometry'}
    ERROR_CODES = {201 : {'description' : 'Error attempting to merge geometry files - no libcheck script'},
                   202 : {'description' : 'Error attempting to merge geometry files - failed creating working directory'},
                   203 : {'description' : 'Error attempting to merge geometry files - setting libcheck parameters'},
                   204 : {'description' : 'Error attempting to merge geometry files - running libcheck'},
                   205 : {'description' : 'Error attempting to merge geometry files - failed to run libcheck'}}
  
    def defaultProjectDict(self, projectId=None, projectName=None, create=True):
        from core import CCP4Modules
        projectDir = CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectId=projectId, projectName=projectName)
        if projectDir is None:
            return None
        projectFilesDir = os.path.join(projectDir, 'CCP4_PROJECT_FILES')
        if not os.path.exists(projectFilesDir):
            try:
                os.mkdir(projectFilesDir)
            except:
                return None
        projectDict = os.path.join(projectDir, 'CCP4_PROJECT_FILES', 'refmac_dictionary.cif')
        if not os.path.exists(projectDict):
            from core import CCP4Utils
            shutil.copyfile(os.path.join(CCP4Utils.getCCP4I2Dir(), 'data', 'refmac_dictionary.cif'), projectDict)
        return projectDict

    def saveToDb(self):
        # Dont save (ie import if it is the project monomer library file)
        if self.__str__().count('CCP4_PROJECT_FILES'):
            return [], None,{}
        else:
            return CCP4File.CDataFile.saveToDb(self)

    def mergeInDictFiles(self, dictFileList=[], parentWorkDirectory=None):
        #print 'CDictDataFile.mergeInDictFiles', dictFileList, parentWorkDirectory
        err = CErrorReport()
        try:
            import libcheck
            wrapper = libcheck.libcheck(self)
        except:
            err.append(self.__class__, 201, name=self.objectName())
            return None, err
        try:
            indx = 1
            myDir = os.path.join(parentWorkDirectory, 'libcheck_' + str(indx))
            while os.path.exists(myDir):
                indx += 1
                myDir = os.path.join(parentWorkDirectory, 'libcheck_' + str(indx))
            os.mkdir(myDir)
            wrapper.workDirectory = myDir
        except:
            err.append(self.__class__, 202, name=self.objectPath())
            return None, err
        try:
            wrapper.container.controlParameters.RUN_MODE = 'MERGE'
            wrapper.container.inputData.DICTLIB.setFullPath(self.__str__())
            for dictFile in dictFileList:
                if isinstance(dictFile, CCP4File.CDataFile):
                    dictFile=dictFile.fullPath.__str__()
                wrapper.container.inputData.MERGELIST.append(dictFile)
        except:
            err.append(self.__class__, 203, name=self.objectPath())
            return None, err
        try:
            from core import CCP4PluginScript
            status = wrapper.process()
            if status != CCP4PluginScript.CPluginScript.SUCCEEDED:
                err.append(self.__class__, 204, name=self.objectPath())
                return None, err
            else:
                return self.__str__(), err
        except CException as e:
            err.append(e)
            return None,err
        except Exception as e:
            print(e)
            err.append(self.__class__, 205, name=self.objectPath())
            return self.__str__(),err

class CTLSDataFile(CCP4File.CDataFile):
    '''A refmac TLS file'''

    QUALIFIERS = {'fileLabel' : 'tls', 'mimeTypeName' : 'application/refmac-TLS',
                  'mimeTypeDescription' : 'Refmac TLS file', 'guiLabel' : 'TLS coefficients',
                  'toolTip' : 'Definition of model domains for TLS refinement',
                  'fileExtensions' : ['tls'], 'fileContentClassName' : None,
                  'helpFile' : 'model_data#tls_file' }


class CMDLMolDataFile(CCP4File.CDataFile):
    '''A molecule definition file (MDL)'''

    QUALIFIERS = {'fileLabel' : 'mol', 'mimeTypeName' : 'chemical/x-mdl-molfile',
                  'mimeTypeDescription' : 'MDL Molfile', 'guiLabel' : 'Mol file',
                  'toolTip' : 'Structure geometry of ligands for refinement in MDL mol format',
                  'fileExtensions' : ['mol'], 'fileContentClassName' : None,
                  'helpFile' : 'model_data#mol_file' }


class CSeqDataFile(CCP4File.CDataFile):
    '''A sequence file'''

    QUALIFIERS = {'fileLabel' : 'sequence', 'mimeTypeName' : 'application/CCP4-seq',
                  'mimeTypeDescription' : 'Sequence file', 'guiLabel' : 'Sequence',
                  'tooltip' : 'Sequence in any of the common formats (pir,fasta..)',
                  'fileExtensions' : ['seq','pir','fasta'], 'fileContentClassName' : 'CSequence',
                  'downloadModes' : ['uniprotFasta'], 'helpFile' : 'model_data#sequences'}

    ERROR_CODES = {201 : {'description' : 'Error reading sequence file'},
                   202 : {'description' : 'Error in BioPython attempting to identify file type'}}

    def __init__(self, **kw):
        CCP4File.CDataFile.__init__(self, **kw)
        self.__dict__['format'] = None
        self.__dict__['identifiers'] = []

    def qualifiers(self, name=None, default=True, custom=True, contentQualifiers=True):
        if name is not None and name == 'fileExtensions':
            keys = list(EXTLIST.keys())
            keys.remove('seq')
            keys.insert(0, 'seq')
            return keys
        else:
            return CCP4File.CDataFile.qualifiers(self, name=name, default=default, custom=custom, contentQualifiers=contentQualifiers)

    def get(self, name=None):
        rv = CCP4Data.CData.get(self, name)
        if name is None:
            rv['fileContent'] = self.__dict__['_fileContent']
        return rv

    def validity(self, arg):
        v = CCP4File.CDataFile.validity(self, arg)
        '''
        if v.maxSeverity() in [SEVERITY_UNDEFINED,SEVERITY_UNDEFINED_ERROR]:
          if arg.has_key('fileContent') and arg['fileContent'] is not None:
            v = arg['fileContent'].validity(arg['fileContent'].get())
            #print 'CSeqDataFile.validity content',content.get(),v
        print 'CSeqDataFile.validity',len(v),v.report()
        '''
        return v

    def updateData(self):
        self.identifyFile()
        if self.__dict__.get('_fileContent', None) is not None:
            self.loadFile()
        self.dataChanged.emit()

    def loadFile(self):
        if self.__dict__.get('_fileContent',None) is None:
            self.__dict__['_fileContent'] = CSequence(parent=self, name='fileContent')
        if self.__dict__.get('format', 'unset') == 'unset':
            self.identifyFile()
        self.__dict__['_fileContent'].loadFile(self.__str__(), format= self.__dict__['format'])

    def saveSequence(self, jobId=None):
        annotation = self.annotation.__str__()
        if annotation[-8:] != '(edited)':
            annotation = annotation +' (edited)'
        self.importFile(jobId=jobId, annotation=annotation, edited=True)
        
    def importFile(self, fileName=None, jobId=None, annotation=None, validatedFile=None, edited=False, jobNumber=None):
        #print 'CSeqDataFile.importFile',fileName,self.fileContent.identifier,'validatedFile',validatedFile,'fileContent.sequence',self.fileContent.sequence.isSet()
        #print 'CSeqDataFile.importFile','dict validatedFile',self.__dict__.get('validatedFile',None)
        #import traceback
        #traceback.print_stack(limit=5)
        self.blockSignals(True)
        #if edited:
        #  self.__dict__['sourceFileName'] = 'COPYDONE'
        #else:
        self.__dict__['sourceFileName'] = self.__str__()
        if fileName is None:
            fileName = self.importFileName(jobId=jobId,ext='.fasta')
        if 'validatedFile' in self.__dict__:
            validatedFile = self.__dict__.get('validatedFile',None)
            del self.__dict__['validatedFile']
        if validatedFile is not None:
            shutil.copyfile(validatedFile,fileName)
        elif self.fileContent.sequence.isSet():
            self.fileContent.saveFile(fileName=fileName)
        else:
            shutil.copyfile(self.__str__(), fileName)
        self.setFullPath(fileName)
        # Set source information up for the PROJECTSMANAGER.importFiles()
        if annotation is not None:
            print('importFile setting annotation', annotation)
            self.annotation = annotation
        elif not self.annotation.isSet():
            self.annotation = re.sub('[^a-zA-Z0-9_-]', '_', str(self.fileContent.identifier))
        self.__dict__['sourceFileAnnotation'] = self.annotation.__str__()
        self.__dict__['format'] = 'internal'
        self.dbFileId.unSet()
        self.blockSignals(False)
        self.dataChanged.emit()
    
    def identifyFile(self, filename=None):
        self.__dict__['format'] = 'unknown'
        self.__dict__['identifiers'] = []
        if self.dbFileId.isSet():
            self.__dict__['format'] = 'internal'
            return
        if filename is None: filename = self.__str__()
        # Try if its a uniprot xml
        nsmap = {'u' : 'http://uniprot.org/uniprot'}
        try:
            #print 'identifyFile trying uniprot',filename
            from core import CCP4Utils
            root = CCP4Utils.openFileToEtree(filename)
            eleList = root.xpath('./u:entry/u:accession', namespaces=nsmap)
            #print 'identifyFile eleList', eleList
        except:
            pass
        else:
            if len(eleList) > 0:
                self.__dict__['format'] = 'uniprot'
                self.__dict__['identifiers'] = [str(eleList[0].text)]
                return
        # Return the file format and list of seq identifiers in the file
        # This requires BioPython tools
        if not BIOPYTHON:
            text = CCP4Utils.readFile(str(filename))
            lines = text.split('\n')
            if '>P1;' in lines[0]:
                formt = 'pir'
            elif '>' in lines[0]:
                formt = 'fasta'
            else:
                formt = 'unknown'
            self.__dict__['format'] = formt
            self.__dict__['identifiers'] = []
        else:
            if not os.path.exists(filename):
                return
            ext = os.path.splitext(str(filename))[1][1:]
            firstFormat = EXTLIST.get(ext, None)
            testOrder = SEQFORMATLIST + ALIGNFORMATLIST
            if firstFormat is not None:
                testOrder.remove(firstFormat)
                testOrder.insert(0, firstFormat)
            #print 'identifyFile firstFormat', firstFormat, testOrder
            n = 0
            formt = 'unknown'
            while formt == 'unknown' and n < len(testOrder):
                try:
                    err, rv = self.bioGetSeqIdentifiers(filename, testOrder[n])
                except:
                    err = CErrorReport(self.__class__, 'Error attempting to read as: ' + testOrder[n])
                print('identifyFile', err, 'rv=', rv)
                if err.maxSeverity() <= SEVERITY_WARNING:
                    if testOrder[n] == 'fasta':
                        #Beware biopython allows a 'noformat' file as fasta
                        format0 = self.fileContent.simpleFormatTest(filename)
                        #print 'simpleFormatTest',format0
                        if format0 == 'noformat':
                            tmpFile, err0, rv0 = self.fileContent.fixFastaFile(filename)
                            rv = [os.path.splitext(os.path.basename(filename))[0]]
                            if tmpFile is not None:
                                formt = testOrder[n]
                                self.__dict__['validatedFile'] = tmpFile
                        elif format0 == 'fasta':
                            formt = testOrder[n]
                    else:
                        formt = testOrder[n]
                    n = n + 1
                elif testOrder[n] == 'pir':
                    tmpFile, err0, rv0 = self.fileContent.fixPirFile(filename)
                    if tmpFile is not None:
                        try:
                            err,rv = self.bioGetSeqIdentifiers(tmpFile, testOrder[n])
                        except:
                            err = CErrorReport(self.__class__, 'Error attempting to read as: ' + testOrder[n])
                        if err.maxSeverity() <= SEVERITY_WARNING:
                            formt = testOrder[n]
                        else:
                            n = n + 1
                    else: 
                        n = n + 1
                else:
                    n = n + 1
            self.__dict__['format'] = formt
            self.__dict__['identifiers'] = rv
        #print 'CSeqDataFile.identifyFile', filename, n, self.__dict__['format'], self.__dict__['identifiers']
        return

    def bioGetSeqIdentifiers(self, filename, format):
        #print 'bioGetSeqIdentifiers',filename,formatss
        err = CErrorReport()
        try:
            if format in SEQFORMATLIST:
                seq_records = list(Bio.SeqIO.parse(filename, format))
                #print 'bioGetSeqIdentifiers', seq_records
                ali_records = []
            else:
                ali_records = list(Bio.AlignIO.parse(filename, format))
                seq_records = []
                for ali in ali_records:
                    seq_records.extend(ali)
        except Exception as e:
            print('Error extracting sequence identifiers from', filename, 'reading as', format)
            print(e)
            err.append(self.__class__, 201, 'Reading ' + filename+' as ' + str(format) + '\n' + str(e), stack=False)
            return err, []
        #print 'bioGetSeqIdentifiers',filename,format,seq_records,ali_records
        idList = []
        for r in seq_records:
            try:
                idList.append(r.id)
            except:
                err.append(self.__class__, 105, filename, stack=False)
        #print 'bioGetSeqIdentifiers',idList
        return err, idList


class CSeqDataFileList(CCP4Data.CList):
    SUBITEM = {'class' : CSeqDataFile , 'qualifiers' : {'allowUndefined' : False}}
    ERROR_CODES = {150 : {'description' : 'No file content information'},
                   151 : {'description' : 'Two sequences have the same identifier'},
                   152 : {'description' : 'Failed in merging sequence files to read sequence file'},
                   153 : {'description' : 'Failed in merging sequence files to write merged file'}}

    def validity(self, arg):
        err = CCP4Data.CList.validity(self, arg)
        for ii in range(1, len(arg)):
            #print 'CSeqDataFileList.validity', ii, arg
            if isinstance(arg[ii], CSeqDataFile):
                iContent = arg[ii].fileContent
            else:
                iContent = arg[ii].get('fileContent', None)
            if iContent is None:
                err.append(self.__class__, 150, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                for jj in range(0, ii):
                    if isinstance(arg[jj], CSeqDataFile):
                        jContent = arg[jj].fileContent
                    else:
                        jContent = arg[jj].get('fileContent', None)
                    if jContent is not None:
                        #print 'CSeqDataFileList.validity', ii, iContent.identifier, jj, jContent.identifier
                        if iContent.identifier.__str__() == jContent.identifier.__str__():
                            err.append(self.__class__, 151, jContent.identifier.__str__(), name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err

    def mergeFiles(self, outputFile):
        from core import CCP4Utils
        err = CErrorReport()
        text = ''
        for seqFileObj in self.__dict__['_value']:
            try:
                text = text + CCP4Utils.readFile(str(seqFileObj)) + '\n'
            except:
                err.append(self.__class__, 152, str(seqFileObj))
        #print 'CSeqDataFileList.mergeFiles', self.__dict__['_value'], text
        try:
            CCP4Utils.saveFile(text=text, fileName=outputFile)
        except:
            err.append(self.__class__, 153, str(outputFile))
        return err


class CSequenceAlignment(CCP4Data.CData, CBioPythonSeqInterface):
    '''
    An alignment of two or more sequences.
    Each sequence is obviously related to class CSequence, but
    will also contain gaps relevant to the alignment. We could
    implement the contents as a list of CSequence objects?
    The alignment is typically formatted in a file as consecutive 
    or interleaved sequences.
    '''

    CONTENTS = {'identifier' : {'class' : CCP4Data.CString,
                                'qualifiers' : { 'toolTip' : 'Optional convenient name for sequence alignment'}},
                'moleculeType' : {'class' : CCP4Data.CString,
                                  'qualifiers' : {'onlyEnumerators' : True, 'enumerators' : ['PROTEIN','NUCLEIC'],
                                                  'menuText' : ['protein','nucleic acid'], 'default' : 'PROTEIN',
                                                  'toolTip' : 'Molecule type'}}}
    CONTENTS_ORDER = ['identifier', 'moleculeType']

    def loadFile(self, fileName, format='unknown'):
        '''
        import os
        if format == 'internal':
          self.loadInternalFile(fileName)
        elif format == 'uniprot':
          #Try treating it as a uniprot file
          m = CSequenceMeta()
          try:
            ret = m.loadFromUniprotXml(fileName)
          except:
            pass
          else:
            self.sequence = ret['sequence']
            if ret['db']=='Swiss-Prot':
              self.referenceDb = 'sp'
            if len(ret['accessionList'][0])>0:
              self.reference = ret['accessionList'][0]
            self.identifier = ret['proteinFullName']
            self.__dict__['lastLoadedFile'] = str(fileName)
            self.__dict__['lastLoadedFormat'] = 'uniprot'
            return
        else:
        '''
        self.unSet()
        CBioPythonSeqInterface.loadExternalFile(self, fileName, format=format)


class CSeqAlignDataFile(CCP4File.CDataFile):
    '''A (multiple) sequence alignment file'''
    QUALIFIERS = { }
    QUALIFIERS.update(CCP4File.CDataFile.QUALIFIERS)
    QUALIFIERS.update({'mimeTypeName' : 'application/CCP4-seqalign',
                       'mimeTypeDescription' : 'Sequence alignment file',
                       'fileExtensions' : ['aln', 'pir', 'fasta', 'msf', 'phy'],
                       'fileContentClassName' : 'CSequenceAlignment',
                       'guiLabel' : 'Aligned sequences',
                       'toolTip' : 'Multiple sequence alignment in any of the common formats (pir,fasta..)',
                       'helpFile' : 'model_data#alignments',
                       'requiredSequences' : NotImplemented})
    QUALIFIERS_ORDER = ['requiredSequences']
    QUALIFIERS_DEFINITION = {'requiredSequences' : {'type' : list, 'listItemType' : int,
                             'description' : 'A list of allowed numbers of sequences in file (usually [2])'}}
    ERROR_CODES = {202 : {'description' : 'Error reading from file'},
                   203 : {'description' : 'Unknown alignment file format'},
                   204 : {'description' : 'Can not read Blast or HHPred file format'},
                   205 : {'description' : 'Error reading identifiers from multi-record file'},
                   206 : {'description' : 'Error attempting to identify file format'},
                   250 : {'description' : 'Alignment file format not recognised - can not convert'},
                   251 : {'description' : 'Alignment file conversion failed to overwrite existing file'},
                   252 : {'description' : 'Alignment file conversion failed writing file'},
                   260 : {'description' : 'Alignment file does not contain required number of sequences'}}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, keywords={}, **kw):
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, fullPath=fullPath, keywords=keywords, **kw)
        self.__dict__['format'] = None
        self.__dict__['identifiers'] = []

    def updateData(self, **kw):
        self.identifyFile()
        self.loadFile()
        #print 'CSeqAlignDataFile.updateData', self.__dict__['format'], self.__dict__['identifiers']
        CCP4File.CDataFile.updateData(self, **kw)
    
    def qualifiers(self, name=None, default=True, custom=True, contentQualifiers=True):
        if name is not None and name == 'fileExtensions':
            keys = list(EXTLIST.keys())
            keys.remove('aln')
            keys.insert(0, 'aln')
            return keys
        else:
            return CCP4File.CDataFile.qualifiers(self, name=name, default=default, custom=custom, contentQualifiers=contentQualifiers)

    def identifyFile(self, filename=None):
        self.__dict__['format'] = 'unknown'
        self.__dict__['identifiers'] = []
        if filename is None:
            filename = self.__str__()
        '''
        # Try if its a uniprot xml
        nsmap = { 'u' : 'http://uniprot.org/uniprot' }   
        try:
          #print 'identifyFile trying uniprot',filename
          from core import CCP4Utils
          root = CCP4Utils.openFileToEtree(filename)
          eleList = root.xpath('./u:entry/u:accession',namespaces=nsmap)
          #print 'identifyFile eleList',eleList
        except:      
          pass
        else:
          if len(eleList)>0:
            self.__dict__['format'] = 'uniprot'
            self.__dict__['identifiers'] = [ str(eleList[0].text) ]
            return
        # Return the file format and list of seq identifiers in the file
        # This requires BioPython tools
        if not BIOPYTHON:
          text = CCP4Utils.readFile(str(filename))
          lines = text.split('\n')
          if '>P1;' in lines[0]:
            formt = 'pir'
          elif '>' in lines[0]:
            formt = 'fasta'
          else:
            formt = 'unknown'
          self.__dict__['format'] = formt
          self.__dict__['identifiers'] = []
        else:
        '''
        allErrors = CErrorReport()
        if 1:
            if not os.path.exists(filename):
                return
            ext = os.path.splitext(str(filename))[1][1:]
            firstFormat = EXTLIST.get(ext, None)
            if firstFormat == 'blast':
                self.__dict__['format'] = firstFormat
                self.__dict__['identifiers'] = []
                return self.__dict__['format'], self.__dict__['identifiers']
            testOrder = ALIGNFORMATLIST
            if firstFormat is not None and firstFormat in testOrder:
                testOrder.remove(firstFormat)
                testOrder.insert(0, firstFormat)
            #print 'testOrder', testOrder, 'firstFormat', firstFormat
            n = 0
            formt = 'unknown'
            rv = []
            while formt == 'unknown' and n < len(testOrder):
                try:
                    err, rv = self.bioGetSeqIdentifiers(filename, testOrder[n])
                except Exception as e:
                    err = CErrorReport(self.__class__, 206, details='Attempting to read as format:' + testOrder[n] + '\n' + str(e), stack=False)
                    #print 'bioGetSeqIdentifiers', testOrder[n], str(err)
                if err.maxSeverity() <= SEVERITY_WARNING:
                    formt = testOrder[n]
                else:
                    allErrors.extend(err)
                    n = n + 1
            self.__dict__['format'] = formt
            self.__dict__['identifiers'] = rv
        #print 'CSeqAlignDataFile.identifyFile errors', allErrors.report()
        #print 'CSeqAlignDataFile.identifyFile', self.__dict__['format'], self.__dict__['identifiers']
        return self.__dict__['format'],self.__dict__['identifiers']

    def bioGetSeqIdentifiers(self, fileName=None, format=None):
        if fileName is None:
            fileName = self.__str__()
        if format is None:
            if self.__dict__['format'] is None:
                self.identifyFile()
            format = self.__dict__['format']
        err = CErrorReport()
        #print 'into bioGetSeqIdentifiers', fileName, format
        seq_records = []
        #try:
        if 1:
            if format in BLASTFORMATLIST or format in HHPREDFORMATLIST:
                err.append(self.__class__, 204, str(fileName))
            elif format in ALIGNFORMATLIST:
                try:
                    ali_records = list(Bio.AlignIO.parse(fileName, format))
                    for ali in ali_records:
                        seq_records.extend(ali)
                except Exception as e:
                    print('Error extracting sequence identifiers from', fileName, 'reading as', format)
                    print(e)
                    err.append(self.__class__, 203, str(fileName))
            else:
                try:
                    seq_records = list(Bio.SeqIO.parse(fileName, format))
                    #print 'Extracting sequence identifiers:', seq_records
                except Exception as e:
                    print('Error extracting sequence identifiers from', fileName, 'reading as', format)
                    print(e)
                    err.append(self.__class__, 203, str(fileName))
            '''
            except Exception as e:
            print 'CSeqAlignData.bioGetSeqIdentifiers',format,e
            err.append(self.__class__,202,str(fileName)+' as '+str(format)+': '+str(e),stack=False)
            return err,[]
            else:
            '''
            #print 'seq_records',seq_records
            #if len(seq_records)>0: print 'seq_records[0]',seq_records[0]
            idList = []
            for r in seq_records:
                try:
                    idList.append(r.id)
                except:
                    err.append(self.__class__, 205, fileName, stack=False)
        return err, idList

    def convertFormat(self, toFormat, fileName, reorder=None):
        if self.__dict__['format'] is None or self.__dict__['format'] == 'unknown':
            self.identifyFile()
        if self.__dict__['format'] == 'unknown':
            return CErrorReport(self.__class__, 250, self.__str__() + ' to ' + str(fileName), stack=False)
        #print 'CSeqAlignDataFile.convertFormat', self.__str__(), self.__dict__['format'], reorder
        # Try reading the input file
        try:
            input_handle = open(self.__str__(), "rU")
            alignments = Bio.AlignIO.read(self.__str__(), self.__dict__['format'])
        except:
            return CErrorReport(self.__class__, 202, self.__str__(), stack=False)
        if reorder is None:
            pass
        if reorder == 'reverse':
            aliout = Bio.Align.MultipleSeqAlignment([])
            for ii in range(len(alignments) - 1, -1, -1):
                aliout.append(alignments[ii])
            alignments = aliout
        elif isinstance(reorder, list):
            pass
        try:
            if os.path.exists(fileName):
                os.remove(fileName)
        except:
            return CErrorReport(self.__class__, 251, fileName)
        out = open(fileName, "w")
        try:
            #print 'toFormat',toFormat
            #print 'alignments',alignments
            #print 'CSeqAlignDataFile.convertFormat output',fileName
            Bio.AlignIO.write(alignments, out ,toFormat )
        except Exception as e:
            err = CErrorReport(self.__class__, 252, 'format: ' + toFormat, exc_info=sys.exc_info(), stack=False)
        else:
            err = CErrorReport()
        try:
            input_handle.close()
        except:
            pass
        try:
            out.close()
        except:
            pass
        return err

    def validity(self, arg):
        v = CCP4File.CDataFile.validity(self, arg)
        if v.maxSeverity() > SEVERITY_WARNING:
            return v
        reqNSeq = self.qualifiers('requiredSequences')
        #print 'CSeqAlignDataFile.validity',reqNSeq,type(reqNSeq)
        if reqNSeq is NotImplemented or len(reqNSeq) == 0:
            return v
        err, idList = self.bioGetSeqIdentifiers()
        #print 'CSeqAlignDataFile.validity',err.maxSeverity(),idList
        if err.maxSeverity() > SEVERITY_WARNING:
            v.extend(err)
            return v
        if not len(idList) in reqNSeq:
            v.append(self.__class__, 260, self.stripedName() + ' contains ' + str(len(idList)) + ' sequences, requires ' + str(reqNSeq))
        return v


class CHhpredDataFile(CCP4File.CDataFile):

    QUALIFIERS = {'fileLabel' : 'HHPred sequence search', 'mimeTypeName' : 'application/HHPred-alignments',
                  'mimeTypeDescription' : 'HHPred sequence search results',
                  'guiLabel' : 'HHPred results', 'tooltip' : 'Output from HHPred search',
                  'fileExtensions' : ['hhr'], 'fileContentClassName' : 'CHhpredData',
                  'helpFile' : 'model_data#ali'}

class CHhpredItem(CCP4Data.CData):

    CONTENTS = {'annotation' : {'class' : CCP4Data.CString},
                'identifier' : {'class' : CCP4Data.CString},
                'chain' : {'class' : CCP4Data.CString}}

    # Support Qt QAbstractItemModel
    def data(self, role):
        if role == QtCore.Qt.DisplayRole:
            return self.__dict__['_value']['annotation'].__str__()
        elif role == QtCore.Qt.ToolTipRole:
            #print 'CHhpredItem parent',type(self.parent()),type(self.parent().parent())
            return "<pre>" + self.parent().parent().getAlignmentText(self.row()) + '</pre>'
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def child(self, row):
        return None

    def rowCount(self, mI):
        return 0

    def columnCount(self):
        return 1

    def abstractModelParent(self):
        #return self.parent().__dict__['_abstractModelParent']
        return QtCore.QModelIndex()

    def row(self):
        return self.parent().__dict__['_value'].index(self)

    def index(self):
        return self.parent().__dict__['_value'].index(self)


class CHhpredData(CCP4File.CDataFileContent):

    ERROR_CODES = {201 : {'description' : 'Failed to read HHPred file'},
                   202 : {'description' : 'Failed to load iotbx software to read HHPred file'}}

    CONTENTS = {'alignmentList' : {'class' :CCP4Data.CList, 'subItem' : {'class' :CHhpredItem}}}

    def hhpredParser(self, fileName):
        from core import CCP4Utils
        try:
            text = CCP4Utils.readFile(str(fileName))
        except Exception as e:
            raise CException(self.__class__, 201, fileName)
        #print 'loadFile',text
        try:
            import iotbx
            from iotbx import bioinformatics
            hp = iotbx.bioinformatics.hhsearch_parser(text)
        except Exception as e:
            print("FAILED TO PARSE HHPRED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            raise CException(self.__class__, 202)
        return hp

    def loadFile(self, fileName=None):
        print('CHhpredData.loadFile', fileName)
        hp = self.hhpredParser(fileName)
        self.__dict__['_value']['alignmentList'].unSet()
        for hit in hp.hits():
            self.__dict__['_value']['alignmentList'].addItem()
            self.__dict__['_value']['alignmentList'][-1].annotation.set(str(hit.annotation))
            self.__dict__['_value']['alignmentList'][-1].chain.set(str(hit.chain))
            self.__dict__['_value']['alignmentList'][-1].identifier.set(str(hit.identifier))

    def getAlignmentText(self, indx):
        indx = int(indx)
        hp = self.hhpredParser(self.parent().__str__())
        ii = -1
        for hit in hp.hits():
            ii += 1
            if ii == indx:
                return 'CLUSTAL from ' + str(hit.alignment)
        return ' '


class CBlastDataFile(CCP4File.CDataFile):

    QUALIFIERS = {'fileLabel' : 'Blast sequence search', 'mimeTypeName' : 'application/Blast-alignments',
                  'mimeTypeDescription' : 'Blast sequence search results',
                  'guiLabel' : 'Blast results', 'tooltip' : 'Output from Blast search',
                  'fileExtensions' : ['bla', 'blast', 'xml'], 'fileContentClassName' : 'CBlastData',
                  'helpFile' : 'model_data#ali'}


class CBlastItem(CCP4Data.CData):

    CONTENTS = {'hitId' : {'class' : CCP4Data.CString},
                'querySequence' : {'class' : CCP4Data.CString},
                'hitSequence' : {'class' : CCP4Data.CString}}

    # Support Qt QAbstractItemModel used for selection in ProvideAlignment task
    def data(self, role):
        #print 'CBlastItem.data',self.__dict__['_value']['hitId']
        if role == QtCore.Qt.DisplayRole:
            return self.__dict__['_value']['hitId'].__str__()
        elif role ==  QtCore.Qt.ToolTipRole:
            return "<pre>" + self.parent().parent().getAlignmentText(self.row()) + '</pre>'
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def child(self, row):
        return None

    def rowCount(self, mI):
        return 0

    def columnCount(self):
        return 1

    def abstractModelParent(self):
        #return self.parent().__dict__['_abstractModelParent']
        return QtCore.QModelIndex()

    def row(self):
        return self.parent().__dict__['_value'].index(self)

    def index(self):
        return self.parent().__dict__['_value'].index(self)


class CBlastData(CCP4File.CDataFileContent):

    CONTENTS = {'queryId' : {'class' : CCP4Data.CString },
                'alignmentList' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CBlastItem}}}
    ERROR_CODES = {201 : {'description' : 'Failed reading blast file'},
                   202 : {'description' : 'Blast file contains results of more than one query - only the first is read', 'severity' : SEVERITY_WARNING},
                   203 : {'description' : 'Failed parsing Blast file'}}

    def loadFile(self, fileName=None):
        # Bio.SearchIo formats: "hmmer3-tab', 'blast-xml', 'phmmer3-domtab', 'exonerate-text', 'hmmer3-text',
        # 'blast-text', 'hmmer2-text', 'exonerate-vulgar', 'exonerate-cigar', 'blast-tab', 'fasta-m10', 'hmmsearch3-domtab', 'blat-psl', 'hmmscan3-domtab"
        import Bio.SearchIO
        mode = None
        fileName = str(fileName)
        try:
            f = open(fileName, 'r')
        except:
            raise CException(self.__class__, 201, fileName)
        if mode is None:
            if fileName.endswith('xml'):
                mode = 'blast-xml'
            else:
                mode = 'blast-text'
        try:
            qresults = Bio.SearchIO.parse(fileName, mode)
        except:
            raise CException(self.__class__, 203, fileName)
        first = True
        for qresult in qresults:
            if first:
                first = False
                self.queryId.set(str(qresult.id))
                for hit in qresult.hits:
                    self.alignmentList.addItem()
                    self.alignmentList[-1].hitId = hit.fragments[0].hit_id
                    self.alignmentList[-1].querySequence.set(hit.fragments[0].query.seq)
                    self.alignmentList[-1].hitSequence.set(hit.fragments[0].hit.seq)
                    print(self.alignmentList[-1])
            else:
                raise  CException(self.__class__, 202, fileName)

    def getBioPyVersion(self):
        from Bio import __version__ as bioversion
        try:
            version=float(bioversion)
            return version
        except:
            return None

    def getAlignmentText(self, indx):

        bioversion=self.getBioPyVersion()
        if bioversion is not None:
            if bioversion < 1.79:
                from Bio.Alphabet import generic_protein
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Align import MultipleSeqAlignment
        try:
            if bioversion is not None:
                if bioversion < 1.79:
                    a = SeqRecord(Seq(self.alignmentList[indx].querySequence.__str__(), generic_protein), id=self.queryId.__str__())
                    b = SeqRecord(Seq(self.alignmentList[indx].hitSequence.__str__(), generic_protein), id=self.alignmentList[indx].hitId.__str__())
                else:
                    a = SeqRecord(Seq(self.alignmentList[indx].querySequence.__str__(), id="prot1", annotations={"molecule_type": "protein"}), id=self.queryId.__str__())
                    b = SeqRecord(Seq(self.alignmentList[indx].hitSequence.__str__(), id="prot2", annotations={"molecule_type": "protein"}), id=self.alignmentList[indx].hitId.__str__())
            else:
                a = SeqRecord(Seq(self.alignmentList[indx].querySequence.__str__(), id="prot1", annotations={"molecule_type": "protein"}), id=self.queryId.__str__())
                b = SeqRecord(Seq(self.alignmentList[indx].hitSequence.__str__(), id="prot2", annotations={"molecule_type": "protein"}), id=self.alignmentList[indx].hitId.__str__())
            #a = SeqRecord(Seq(self.alignmentList[indx].querySequence.__str__(), generic_protein), id=self.queryId.__str__())
            #b = SeqRecord(Seq(self.alignmentList[indx].hitSequence.__str__(), generic_protein), id=self.alignmentList[indx].hitId.__str__())
            align = MultipleSeqAlignment([a, b], annotations={"tool": "CCP4i2"})
            f = StringIO()
            count = Bio.AlignIO.write(align, f, "clustal")
        except:
            return 'Error creating alignment text'
        else:  
            return f.getvalue()


class CElement(CCP4Data.COneWord):
    ''' Chemical element '''

    QUALIFIERS = {'onlyEnumerators' : True,
                  'enumerators' : ['H','He',
                   'Li','Be','B','C','N','O','F','Ne',
                   'Na','Mg','Al','Si','P','S','Cl','Ar',
                   'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                   'Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',
                   'Cs','Ba',
                        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
                        'Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',
                   'Fr','Ra',
                        'Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
                        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn']}


class CMonomer(CCP4Data.CData):
    '''A monomer compound. ?smiles'''
    CONTENTS = {'identifier' : {'class' : CCP4Data.CString,
                                'qualifiers' : {'toolTip' : 'The name you use for the monomer'}},
                'formula' : {'class' : CCP4Data.CString,
                             'qualifiers' : { 'toolTip' : 'The formula for the monomer'}},
                'dictionaryName' : {'class' : CCP4Data.CString,
                                    'qualifiers' : {'toolTip' : 'The REFMAC dictionary name if not the same as the name'}},
                'smiles' : {'class' : CCP4Data.CString,
                            'qualifiers' : {'toolTip' : 'The smiles string for the monomer'}},}
    CONTENTS_ORDER = ['identifier', 'formula', 'dictionaryName', 'smiles']

    def guiLabel(self):
        #print 'CMonomer.getGuiLabel'
        if self.identifier.isSet():
            return str(self.identifier)
        elif self.formula.isSet():
            return str(self.formula)
        else:
            return str(self.objectName())
        self.identifier.dataChanged.connect(self.bleep)

    @QtCore.Slot()
    def bleep(self):
        print('CMonomer.identifier dataChanged')

    def getTextItem(self):
        return self.guiLabel()

    def __eq__(self, arg):
        # Considered equal if name same identifier
        #print 'CMonomer.__eq__',arg.get('identifier')
        return self.__dict__['_value']['identifier'].__eq__(arg.get('identifier'))
  
class CPdbData(CCP4File.CDataFileContent):
    '''Contents of a PDB file - a subset with functionality for GUI'''

    ERROR_CODES = {101 : {'description' : 'Unable to load mmdb - ensure LD_LIBRARY_PATH is set'},
                   102 : {'description' : 'Error reading PDB file into MMDB object'},
                   103 : {'description' : 'Residue range selection does not specify chain'},
                   104 : {'description' : 'Residue range selection specifies non-existant chain id'},
                   105 : {'description' : 'Residue range selection - no residues selected'},
                   106 : {'description' : 'Residue range selection - residue number is not an integer'},
                   112 : {'description' : 'Atom selection failed. Failed creating CMMDBManager object'},
                   113 : {'description' : 'Atom selection failed. Faied reading coordinate file.'},
                   114 : {'description' : 'Atom selection failed. Failed parsing command'},
                   115 : {'description' : 'Atom selection failed. Error creating PPCAtom'},
                   116 : {'description' : 'Atom selection failed. Error in GetSelIndex'},
                   117 : {'description' : 'Atom selection failed. Error loading selection tree'},
                   118 : {'description' : 'Atom selection failed. Error applying selection tree'},
                   119 : {'description' : 'Creating new PDB file failed on writing file'},
                   120 : {'description' : 'Creating new PDB file failed converting from fractional coordinates'}}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, **kw):
        #print 'CPdbData.__init__'
        qualis = {}
        if len(qualifiers) > 0:
            qualis.update(qualifiers)
        if len(kw) > 0:
            qualis.update(kw)
        CCP4Data.CData.__init__(self, qualifiers=qualis, parent=parent, name=name)
        self.__dict__['_molHnd'] = None
        self.__dict__['_composition'] = None
        self.__dict__['_sequences'] = None
        self.__dict__['lastLoadedFile'] = None

    def loadFile(self, fileName=None):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        #print 'CPdbData.loadFile',fileName
        from core import CCP4Utils
        if fileName is None:
            raise CException(self.__class__, 101, str(fileName), name=self.objectPath())
        # Convert CDataFile or CFilePath to string
        fileName = str(fileName)
        if not os.path.exists(fileName):
            raise CException(self.__class__, 101, str(fileName), name=self.objectPath())
        try:
            import ccp4mg
        except:
            raise CException(self.__class__, 101, name=self.objectPath())
        if self.__dict__['_molHnd'] is not None:
            del self.__dict__['_molHnd']
            self.__dict__['_molHnd'] = None
        try:
            self.__dict__['_molHnd'] = mmdb.Manager()
            #print fileName, type(fileName) 
            #print  self.__dict__['_molHnd']
            self.__dict__['_molHnd'].ReadCoorFile(str(fileName))
        except:
            raise CException(self.__class__, 102, str(fileName), name=self.objectPath())
        #print 'CPdbData.loadFile nof atoms',self.__dict__['_molHnd'].GetNumberOfAtoms()
        if self.__dict__['_molHnd'] is not None:
            self.__dict__['_composition'] = CPdbDataComposition(self.__dict__['_molHnd'])
        self.dataChanged.emit()

    def __getattr__(self, name):
        if name in ['mmdbManager', 'molHnd']:
            return self.__dict__['_molHnd']
        elif name == 'composition':
            if self.__dict__['_composition'] is None:
                if self.__dict__['_molHnd'] is not None:
                    self.__dict__['_composition'] = CPdbDataComposition(self.__dict__['_molHnd'])
            return self.__dict__['_composition']
        elif name == 'sequences':
            if self.__dict__['_sequences'] is None and self.__dict__['_molHnd'] is not None:
                self.loadSequences(self.__dict__['_molHnd'])
            return self.__dict__['_sequences']
        else:
            return CCP4Data.CData.__getattr__(self, name)

    def interpretResidueInput(self, resid=None):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        insCode = '*'
        if resid is None:
            resNum = mmdb.ANY_RES
        else:
            resid = resid.strip().strip('/')
            if len(resid) == 0:
                resNum = mmdb.ANY_RES
            else:
                rr = resid.split('.')
                try:
                    resNum = int(rr[0])
                except:
                    raise CException(self.__class__, 106, resid, name=self.objectPath())
                if len(rr) > 1:
                    insCode = rr[1]
        return resNum, insCode

    def validRangeSelection(self, chainId=None, firstRes=None, lastRes=None, selModel=1):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        import mmut
        if chainId is None:
            if len(self.__dict__['_composition'].chains) > 0:
                raise CException (self.__class__, 103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                chainId = self.__dict__['_composition'].chains[0]
        else:
            if not self.__dict__['_composition'].chains.count(chainId):
                raise CException (self.__class__, 104, chainId, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        '''
        if firstRes is None and lastRes is None:
          selCom = '/'+str(selModel)+'/'+chainId
        elif lastRes is None:
          selCom = chainId + '/' + str(firstRes)
        else:
          selCom = chainId + '/' + str(firstRes) + '-' + str(lastRes)
        print 'validRangeSelection selCom',selCom
          
        resHnd = self.__dict__['_molHnd'].NewSelection()
        resSel = mmdb.newPPCResidue()
        self.__dict__['_molHnd'].Select(resHnd,mmdb.STYPE_RESIDUE,selCom,mmdb.SKEY_NEW)
        nSel = self.__dict__['_molHnd'].GetSelIndex(resHnd,resSel)
        self.__dict__['_molHnd'].DeleteSelection(resHnd)
        print 'validRangeSelection nSel',nSel
        '''
        if chainId is None or len(chainId.strip()) == 0:
            chainId ='*'
        else:
            chainId = chainId.strip('/')
        firstNum, firstAlt = self.interpretResidueInput(firstRes)
        lastNum, lastAlt = self.interpretResidueInput(lastRes)
        resHnd = self.__dict__['_molHnd'].NewSelection()
        resSel = mmdb.newPPCResidue()
        self.__dict__['_molHnd'].Select (resHnd, mmdb.STYPE_RESIDUE, selModel, chainId, firstNum,
                                         firstAlt, lastNum, lastAlt, '*', '*', '*', '*', mmdb.SKEY_NEW)
        try:
            selindexp = mmut.intp()
            resSel = mmut.GetAtomSelIndex(self.__dict__['_molHnd'], resHnd, selindexp)
            nSel = selindexp.value()
        except:
            nSel = self.__dict__['_molHnd'].GetSelIndex(resHnd, resSel)
        self.__dict__['_molHnd'].DeleteSelection(resHnd)
        #print 'validRangeSelection nSel', nSel
        if nSel == 0:
            raise CException(self.__class__, 105, selCom, name=self.objectPath())   # KJS : Need to fix (no selCom)
        return nSel

    def splitAtomId(self, atomId='', splitRes=False, resRange=False):
        model = ''
        chain = ''
        res = ''
        res2 = ''
        atom = ''
        aid = ''
        atomId = atomId.strip()
        if atomId:
            aid = atomId.split('/')
            if len(aid) <= 1:
                # there are no separators
                atom = atomId
            else:
                # If atomId begins with a slash and model number then expect
                # first item (before first slash) to be the blank
                # If this is not so then insert a model number 
                if aid[0] and len(aid) < 5:
                    aid.insert(0, '')
                    aid.insert(1, '1')
                if len(aid) >= 2:
                    model = aid[1]
                if len(aid) >= 3:
                    chain = aid[2]
                if len(aid) >= 5:
                    atom = aid[4].split('[')[0]
                    # mmdb GetAtomID returns the altLoc after the element type
                    if len(aid[4].split(']'))>1:
                        atom = atom + aid[4].split(']')[1]
                if len(aid) >= 4:
                    if aid[3].count('-'):
                        res, res2 = aid[3].split('-')[0:2]
                    else:
                        if splitRes:
                            res = aid[3].split('(')[0].split('.')
                            if len(res)==1:
                                res.append('')
                            return [model, chain, res[0], res[1], atom]
                        else:
                            res = aid[3]
        #print 'splitAtomID',atomId,aid,[model,chain,res,res2,atom]
        if resRange:
            return [model, chain, res, res2, atom]
        else:
            return [model, chain, res, atom]

    def interpretSelection(self, command, fileName=None):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        from model import CCP4SelectionTree
        import mmut
        try:
            toks = CCP4SelectionTree.SelectionParser().tokenise(command)
        except CException as e:
            raise e
        except:
            raise CException(self.__class__, 114, str(command), name=self.objectPath())
        print('CPdbData.interpretSelection',toks)
        try:
            tree = CCP4SelectionTree.Cselect_tree('TOP' ,mol=self)
            tree.import_command_string(toks[0], toks[1], toks[2])
        except:
            raise CException(self.__class__, 117, str(command), name=self.objectPath())
        try:
            status, selHnd = tree.apply()
        except:
            raise
            raise CException(self.__class__, 118, str(command), name=self.objectPath())
        #print 'CPdbData.interpetSelection',status,selHnd
        try:
            selAtoms = mmdb.newPPCAtom()
        except:
            raise CException(self.__class__, 115, str(command), name=self.objectPath())
        try:
            try:
                selindexp = mmut.intp()
                selAtoms = mmut.GetAtomSelIndex(self.molHnd, selHnd, selindexp)
                nselatoms = selindexp.value()
            except:
                nselatoms = self.molHnd.GetSelIndex(selHnd, selAtoms)
        except:
            raise CException(self.__class__, 116, str(command), name=self.objectPath())
        #print 'CPdbData.interpetSelection nselatoms',nselatoms
        if fileName is not None:
            self.writeSelection(selHnd, fileName)
        # Beware the calling function must clean up the selHnd
        return nselatoms, selHnd

    def writeSelection(self, selHnd, fileName, format='PDB' ):
        # Does an atom-by-atom copy so potentially slow
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        import mmut
        thisMolHnd = mmdb.Manager()
        thisMolHnd.Copy(self.molHnd, mmdb.MMDBFCM_Title | mmdb.MMDBFCM_Cryst)
        if self.molHnd.GetNumberOfModels() > 1 and False:
            #This method does not exist!
            self.molHnd.CopySelection(selHnd, thisMolHnd)
        else:
            try:
                selindexp = mmut.intp()
                selAtoms = mmut.GetAtomSelIndex(self.molHnd, selHnd, selindexp)
                nSelAtoms = selindexp.value()
                for i in range(0,nSelAtoms):
                    pcat = mmdb.getPCAtom(selAtoms, i)
                    if self.molHnd.GetNumberOfModels() == 1 or pcat.GetModelNum() == 1:
                        thisMolHnd.PutAtom(i+1, pcat)
            except:
                selAtoms = mmdb.newPPCAtom()
                nSelAtoms = self.molHnd.GetSelIndex(selHnd, selAtoms)
                for i in range(0, nSelAtoms):
                    pcat = mmdb.CAtomPtr(mmdb.getPCAtom(selAtoms, i))
                    thisMolHnd.PutAtom(i+1, pcat)
        if format == "PDB":
            RC = thisMolHnd.WritePDBASCII(fileName)
        else:
            RC = thisMolHnd.WriteCIFASCII(fileName)
        print(format,'file, containing', nSelAtoms, 'atoms, written to', fileName)
        del thisMolHnd
        return RC

    def makeOneResPDB(self, resName='UNK', atomDefList=[], cell=None, spaceGroup='P 1', fileName=None):
        import ccp4mg
        import mmdb2 as mmdb
        #Trivial test data
        #atomDefList = [ { 'name' : ' SE ', 'element' : 'SE' , 'x' : 45.0, 'y': 45.0, 'z' : 78.0, 'occupancy' : 1.0, 'tempFactor' : 10.0 } ]
        #atomDefList = [ { 'name' : ' SE ', 'element' : 'SE' , 'xFrac' : 0.45, 'yFrac': 0.45, 'zFrac' : 0.78, 'occupancy' : 1.0, 'tempFactor' : 10.0 } ]
        m = mmdb.Manager()
        if cell is None:
            pass
        elif isinstance(cell, CCP4Data.CData):
            m.SetCell(str(cell.x), str(cell.y), str(cell.z), str(cell.alpha), str(cell.beta), str(cell.gamma))
        elif isinstance(cell, list):
            if len(cell) == 3:
                m.SetCell(cell[0], cell[1], cell[2], 90, 90, 90)
            elif len(cell) == 6:
                m.SetCell(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5])
            else:
                m.SetCell(100, 100, 100, 90, 90, 90)
        m.SetSpaceGroup(str(spaceGroup))
        mod = mmdb.Model()
        mod.thisown = 0
        chain = mod.CreateChain("A")
        chain.thisown = 0
        res = mmdb.Residue()
        res.thisown = 0
        res.SetResName(resName)
        res.SetResID(resName, 1, "")
        res.SetChain(chain)
        chain.AddResidue(res)
        for atomDef in atomDefList:
            atom = mmdb.Atom()
            atom.thisown = 0
            res.AddAtom(atom)
            atom.SetAtomName(atomDef.get('atomName'))
            atom.SetElementName(atomDef.get('element', 'C'))
            atom.SetCharge(atomDef.get('charge', 0.0))
            if 'x' in atomDef:
                atom.SetCoordinates(atomDef['x'], atomDef['y'], atomDef['z'], atomDef.get('occupancy', 1.0), atomDef.get('tempFactor', 99))
            elif 'xFrac' in atomDef:
                x, y, z = self.frac2Orth(m, atomDef['xFrac'], atomDef['yFrac'], atomDef['zFrac'])
                if x is None:
                    return CErrorReport(self.__class__, 120, fileName)
                atom.SetCoordinates(x, y, z, atomDef.get('occupancy', 1.0), atomDef.get('tempFactor', 99))
        chain.SetModel(mod)
        m.AddModel(mod)
        m.FinishStructEdit()
        m.PDBCleanup(mmdb.PDBCLEAN_SERIAL | mmdb.PDBCLEAN_INDEX)
        print('CPdbData.makeOneResPDB number of atoms written:', m.GetNumberOfAtoms())
        if fileName is not None:
            rv = m.WritePDBASCII(fileName)
            print('WritePDBASCII', rv)
            if rv != 0:
                return CErrorReport(self.__class__, 119, fileName)
        del m
        return CErrorReport()

    def frac2Orth(self, m, xFrac, yFrac, zFrac):
        try:
            status, x, y, z = m.Frac2Orth(xFrac, yFrac, zFrac)
            return x, y, z
        except:
            try:
                import pygl_coord
                xptr = pygl_coord.doublep()
                yptr = pygl_coord.doublep()
                zptr = pygl_coord.doublep()
                rv = m.Frac2Orth(xFrac, yFrac, zFrac, xptr, yptr, zptr)
                if not rv:
                    return None, None, None
                xp = xptr.value()
                yp = yptr.value()
                zp = zptr.value()
                return xp, yp, zp
            except:
                return None, None, None

    def loadSequences(self, molHnd):
        import collections
        self.__dict__['_sequences'] = collections.OrderedDict()
        amino_acid_code = MODELINFO('amino_acid_code')
        solventList = MODELINFO('solvent', 'ORDER')
        soluteList = MODELINFO('solute', 'ORDER')
        nucleic_acid_code = MODELINFO('nucleic_acid_code')
        for chainId in self.composition.peptides:
            self.__dict__['_sequences'][chainId] = ''
            peptideChain = molHnd.GetChain(1, chainId)
            for i in range(peptideChain.GetNumberOfResidues()):
                residueName = peptideChain.GetResidue(i).GetResName()
                if residueName in amino_acid_code:
                    self.__dict__['_sequences'][chainId]+= amino_acid_code[residueName]
                elif (residueName in solventList) or (residueName in soluteList):
                    pass
                else:
                    print('Residue name not recognised', residueName)
            if len(self.__dict__['_sequences'][chainId]) == 0:
                del self.__dict__['_sequences'][chainId]
        for chainId in self.composition.nucleics:
            if chainId in self.__dict__['_sequences']: continue
            self.__dict__['_sequences'][chainId] = ''
            nucleicChain = molHnd.GetChain(1, chainId)
            for i in range(nucleicChain.GetNumberOfResidues()):
                residueName = nucleicChain.GetResidue(i).GetResName()
                if residueName in nucleic_acid_code:
                    self.__dict__['_sequences'][chainId]+= nucleic_acid_code[residueName]
                elif (residueName in solventList) or (residueName in soluteList):
                    pass
                else:
                    print('Residue name not recognised', residueName)
            if len(self.__dict__['_sequences'][chainId]) == 0:
                del self.__dict__['_sequences'][chainId]

    def writeFasta(self, fileName):
        import copy
        from core import CCP4Utils
        text = ''
        for chainId,sequence in self.sequences:
            text += '>' + chainId + '\n'
            seq = copy.deepcopy(sequence)
            while len(seq) > 0:
                text += seq[0:60] + '\n'
                seq = seq[60:]
        CCP4Utils.saveFile(fileName, text)

class CPdbDataComposition:

    def __init__(self, molHnd):
        self.nModels=molHnd.GetNumberOfModels()
        self.chains = []
        self.peptides = []
        self.nucleics = []
        self.solventChains = []
        self.monomers = []
        self.nresSolvent = 0
        self.moleculeType = []
        self.containsHydrogen = False
        selModel = 1
        #try:
        self.analyseResTypes(molHnd, selModel)
        #except:
        #  print 'ERROR analysing coordinate file composition - maybe due to failure to import mmdb/mmut'
        try:
            self.elements = self.analyseElements(molHnd, selModel)
        except:
            print('ERROR analysing coordinate file composition - maybe due to failure to import mmdb/mmut')

        #print 'CPdbDataComposition.__init__',self.monomers,self.chains

    def analyseElements(self, molHnd, selModel=1):
        try:
            import ccp4mg
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        import mmut
        elements = []
        hydHnd = molHnd.NewSelection()
        hydSel = mmdb.newPPCAtom()
        # Remove common element types
        molHnd.Select(hydHnd, mmdb.STYPE_ATOM, selModel, '*', \
                      mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', '*', '*', mmdb.SKEY_NEW)
        for ele in ['C', 'N', 'O']:
            molHnd.Select(hydHnd,mmdb.STYPE_ATOM, selModel, '*', \
                          mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', ele, '*', mmdb.SKEY_CLR)
        # Are there any novel element types
        try:
            selindexp = mmut.intp()
            hydSel = mmut.GetAtomSelIndex(molHnd, hydHnd, selindexp)
            nother = selindexp.value()
            while nother > 0:
                pcat = mmdb.getPCAtom(hydSel, 0)
                ele = pcat.element.rstrip()
                #print "extra element ", ele
                elements.append(ele)
                if ['H', ' H', 'H '].count(ele):
                    self.contains_hydrogen = True
                molHnd.Select(hydHnd, mmdb.STYPE_ATOM, selModel, '*', mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', ele, '*', mmdb.SKEY_CLR)
                selindexp = mmut.intp()
                hydSel = mmut.GetAtomSelIndex(molHnd, hydHnd, selindexp)
                nother = selindexp.value()
        except:
            nother = molHnd.GetSelIndex(hydHnd, hydSel)
            while nother > 0:
                pcat = mmdb.CAtomPtr(mmdb.getPCAtom(hydSel, 0))
                ele = pcat.element.rstrip()
                #print "extra element ",ele
                elements.append(ele)
                if ['H', ' H', 'H '].count(ele):
                    self.contains_hydrogen = True
                molHnd.Select(hydHnd, mmdb.STYPE_ATOM, selModel, '*', mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', ele, '*', mmdb.SKEY_CLR)
                nother = molHnd.GetSelIndex(hydHnd, hydSel)
        mmdb.delPPCAtom(hydSel)
        molHnd.DeleteSelection(hydHnd)
        #print 'elements', elements
        return elements

    def resTypeSelCommand(self, resType):
        resList =  MODELINFO(resType, 'ORDER')
        com = resList[0]
        for item in resList[1:]:
            com = com + ',' + item
        return com

    def analyseResTypes(self, molHnd, selModel=1):
        try:
            import ccp4mg            # KJS : This function is too long.
            import mmdb2 as mmdb
        except:
            print('FAILED CCP4ModelData imported ccp4mg')
        import mmut
        new_swig_mmdb = True
        try:
            import version
            ccp4mg_version_tup = tuple([int(x) for x in version.ccp4mg_version.split('.')])
        except:
            ccp4mg_version_tup = (2, 7, 0)
        if ccp4mg_version_tup >= (2, 8, 0):
            nChainsp = mmut.intp()
            chainTable = mmut.GetChainTable(molHnd, selModel, nChainsp)
            self.nChains = nChainsp.value()
        else:
            chainTable = mmdb.newPPCChain()
            nChains = mmdb.intp()
            molHnd.GetChainTable(selModel, chainTable, nChains)
            self.nChains = nChains.value()
        hydHnd = molHnd.NewSelection()
        hydSel = mmdb.newPPCAtom()
        resHnd = molHnd.NewSelection()
        resSel = mmdb.newPPCResidue()
        self.chains = []
        self.chainInfo = []
        self.monomers = []
        self.nresSolvent = 0
        self.peptides = []
        self.nucleics = []
        self.solventChains = []
        #self.lipids = []
        self.saccharides = []
        self.nResidues = 0
        self.nAtoms = 0
        self.nAtoms = molHnd.GetNumberOfAtoms()
        for n in range(0, self.nChains):
            if ccp4mg_version_tup >= (2, 8, 0):
                #print 'analyseResTypes',
                pc = mmdb.getPCChain(chainTable, n)
            else:
                pc = mmdb.CChainPtr(mmdb.getPCChain(chainTable, n))
            chainID = pc.GetChainID()
            self.chains.append(chainID)
            nres = pc.GetNumberOfResidues()
            self.chainInfo.append([nres])
            #print 'CPdbData.analyseResTypes chain',chainID,nres
            molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                          mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', '*', '*', mmdb.SKEY_NEW)
            nresaa = 0
            nresna = 0
            nressac = 0
            try:
                selindexp = mmut.intp()
                resSel = mmut.GetResidueSelIndex(molHnd, resHnd, selindexp)
                nrestot = selindexp.value()
                if resSel is not None:
                    #  chainInfo - terminal residues
                    for ir in (0, nres - 1):
                        pcres = mmdb.getPCResidue(resSel, ir)
                        resid = str(pcres.GetResidueID()[0])
                        self.chainInfo[-1].append(resid)
                    com = self.resTypeSelCommand('amino_acid')
                    molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                                  mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                    resSel =  mmut.GetResidueSelIndex(molHnd, resHnd, selindexp)
                    nresaa = nrestot - selindexp.value()
                    com = self.resTypeSelCommand('nucleic_acid')
                    molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                                  mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                    resSel =  mmut.GetResidueSelIndex(molHnd, resHnd, selindexp)
                    nresna = (nrestot - nresaa) - selindexp.value()
                    #print "nresna", chainID,nresna
                    com = self.resTypeSelCommand('saccharide')
                    molHnd.Select(resHnd,mmdb.STYPE_RESIDUE, selModel, chainID, \
                                  mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                    resSel = mmut.GetResidueSelIndex(molHnd, resHnd, selindexp)
                    nressac = (nrestot - nresaa - nresna) - selindexp.value()
                    #print "nressac",chainID,nressac
            except:
                nrestot = molHnd.GetSelIndex(resHnd, resSel)
                com = self.resTypeSelCommand('amino_acid')
                molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                              mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                nresaa = nrestot - molHnd.GetSelIndex(resHnd, resSel)
                com = self.resTypeSelCommand('nucleic_acid')
                molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                              mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                nresna = (nrestot - nresaa) - molHnd.GetSelIndex(resHnd, resSel)
                #print "nresna", chainID, nresna
                com = self.resTypeSelCommand('saccharide')
                molHnd.Select(resHnd, mmdb.STYPE_RESIDUE, selModel, chainID, \
                              mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', com, '*', '*', '*', mmdb.SKEY_XOR)
                nressac = (nrestot - nresaa - nresna) - molHnd.GetSelIndex(resHnd, resSel)
                #print "nressac",chainID,nressac

            chainselHnd = molHnd.NewSelection()
            molHnd.Select(chainselHnd, mmdb.STYPE_ATOM, selModel, chainID, \
                          mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', '*', '*', mmdb.SKEY_NEW)

            solventCom = self.resTypeSelCommand('solvent')
            #soluteCom = self.resTypeSelCommand('solute')
            molHnd.Select(resHnd,mmdb.STYPE_RESIDUE, selModel, chainID, \
                          mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', solventCom, '*', '*', '*', mmdb.SKEY_XOR)
            new_swig_mmdb = False
            try:
                selindexp = mmut.intp()
                resSel = mmut.GetResidueSelIndex(molHnd, resHnd, selindexp)
                nreslig = selindexp.value()
                new_swig_mmdb = True
            except:
                nreslig = molHnd.GetSelIndex(resHnd, resSel)
            nressol = nrestot - (nresaa + nresna + nreslig)
            self.nresSolvent = self.nresSolvent + nressol
            if nresaa > 0:
                self.peptides.append(chainID)
            if nresna > 0:
                self.nucleics.append(chainID)
            if nressac > 0:
                self.saccharides.append(chainID)
            if nressol > 0:
                self.solventChains.append(chainID)
            self.nResidues += nrestot
            #if have_lipids> 0: self.lipids.append(chainID)  # This isn't used for anything yet though, I think.
            #Deal with the remaining residues - check if they are novel
            #amino acid/nucleic/solvent/solute and if so then add them      #to the selection definition for those types.
            #If not then add them to the list of monomers
            novel_restype = []
            if nreslig > 0:
                for ir in range(0, nreslig):
                    if new_swig_mmdb:
                        pcres = mmdb.getPCResidue(resSel, ir)
                    else:
                        pcres = mmdb.CResiduePtr(mmdb.getPCResidue(resSel, ir))
                    resname = pcres.name
                    resid = pcres.GetResidueID()[0]
                    #print "monomer", resname, resid
                    self.monomers.append(resid)

        #print "self.nucleics,peptides,chains", self.nucleics, self.peptides, self.chains
        if len(self.peptides) > 0:
            self.moleculeType.append('PROTEIN')
        if len(self.nucleics) > 0:
            self.moleculeType.append('NUCLEIC')
        if len(self.saccharides) > 0:
            self.moleculeType.append('SACCHARIDE')
        if len(self.monomers) > 0:
            self.moleculeType.append('MONOMER')
        #if len(self.lipids) > 0:
        #  self.molecule_type.append('LIPID')
        molHnd.DeleteSelection(resHnd)
        molHnd.DeleteSelection(hydHnd)

class CAtomSelection(CCP4Data.CData):
    CONTENTS = {'text' : {'class' : CCP4Data.CString}}
    QUALIFIERS = {'pdbFileKey' : ''}
    QUALIFIERS_DEFINITION = {'pdbFileKey' : {'type' :str,
                                             'description' : 'The key for a CPdbDataFile in the same CContainer'}}

    def __str__(self):
        return self.text.__str__()


class CPdbDataFile(CCP4File.CDataFile):

    selectionChanged = QtCore.Signal()
    '''A PDB coordinate file'''
    SUBTYPE_UNKNOWN = 0
    SUBTYPE_MODEL = 1
    SUBTYPE_HOMOLOG = 2
    SUBTYPE_FRAGMENT = 3
    SUBTYPE_HEAVY_ATOMS = 4
    CONTENT_FLAG_PDB = 1
    CONTENT_FLAG_MMCIF = 2
    CONTENTS = {}
    CONTENTS.update(CCP4File.CDataFile.CONTENTS)
    CONTENTS['selection'] = {'class' : CAtomSelection}
    CONTENTS['subType'] = {'class' : CCP4Data.CInt,
                           'qualifiers' : {'default' : 0, 'enumerators' : [0, 1, 2, 3, 4],
                                           'onlyEnumerators':True, 'menuText' : ['unknown', 'model', 'homolog', 'fragment', 'heavy atoms']}}
    QUALIFIERS = {}
    QUALIFIERS.update(CCP4File.CDataFile.QUALIFIERS)
    QUALIFIERS.update({'mimeTypeName' : 'chemical/x-pdb',
                       'mimeTypeDescription' : 'Model coordinates',
                       'fileExtensions' : ['pdb','cif','mmcif','ent'],
                       'fileContentClassName' : 'CPdbData',
                       'fileLabel' : 'coordinates',
                       'guiLabel': 'Atomic model',
                       'toolTip' : 'A model coordinate file in PDB or mmCIF format',
                       'ifInfo' : True,
                       'ifAtomSelection' : False,
                       'downloadModes' : ['ebiPdb','rcsbPdb','uniprotAFPdb'],
                       'helpFile' : 'model_data#coordinate_files'})
    QUALIFIERS_DEFINITION = {'ifAtomSelection' : {'type' :bool,
                                                  'description' : 'Atom selection option enabled'}}
    ERROR_CODES = {401 : {'description' : 'Failed running coord_format to fix coordinate file - is it a PDB file?'},
                   402 : {'severity' : SEVERITY_WARNING, 'description' : 'Badly formated PDB file fixed'},
                   403 : {'severity' : SEVERITY_WARNING,'description' : 'Fixed by removing text'},
                   404 : {'severity' : SEVERITY_WARNING,'description' : 'Fixed by adding text'},
                   405 : {'description' : 'There are no ATOM or HETATM lines in the PDB file'},
                   410 : {'description' : 'No file loaded - can not convert coordinate file format'},
                   411 : {'description' : 'Failed loading file - can not convert coordinate file format'},
                   412 : {'description' : 'Can not overwrite existing file - can not convert coordinate file format'},
                   413 : {'description' : 'Failed writing coordinate file'},
                   414 : {'description' : 'Failed to identify coordinate file format'}}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, **kw):
        #print 'CPdbDataFile.__init__'
        qualis = {}
        if len(qualifiers) > 0:
            qualis.update(qualifiers)
        if len(kw) > 0:
            qualis.update(kw)
        CCP4File.CDataFile.__init__(self, value=value, qualifiers=qualis, parent=parent, name=name, **kw)
        self.selection.dataChanged.connect(self.selectionChanged.emit)
        self.selection.dataChanged.connect(self.dataChanged.emit)

    def updateData(self):
        self.unsetFileContent()
        self.dataChanged.emit()

    def isSelectionSet(self):
        return self.__dict__['_value']['selection'].isSet()

    def requiredContent(self):
        # Return the allowed contentFlag values for CDataFileView.getJobsWithOutputFiles to call to db getJobsWithOutputFiles
        contentList = self.qualifiers('requiredContentFlag')
        if contentList is None or contentList is NotImplemented:
            return None
        else:
            return contentList

    def isMMCIF(self):
        try: 
            gemmi.cif.read(self.fullPath.__str__())
            return True
        except (ValueError, RuntimeError):
            try:
                gemmi.read_structure(self.fullPath.__str__())
                return False
            except (ValueError, RuntimeError):
                raise CException(self.__class__, 414, self.__str__())
    
    def isPDB(self):
        try: 
            gemmi.cif.read(self.fullPath.__str__())
            return False
        except (ValueError, RuntimeError):
            try:
                gemmi.read_structure(self.fullPath.__str__())
                return True
            except (ValueError, RuntimeError):
                raise CException(self.__class__, 414, self.__str__())

    def getSelectedAtomsPdbFile(self, fileName=None):
        if not self.isSelectionSet():
            shutil.copyfile(self.fullPath.__str__(), fileName)
        else:
            self.loadFile()
            nSelAtoms, selHnd = self.fileContent.interpretSelection(self.selection.__str__())
            RC = self.fileContent.writeSelection(selHnd, fileName)
            return RC

    def getTableTextItems(self):
        if self.annotation.isSet():
            text = self.annotation.__str__()
        else:
            text = os.path.split(self.__str__())[1]
        return [text, self.selection.__str__()]

    def assertSame(self, other, diagnostic=False, **kw):
        import tempfile
        from core import CCP4Utils
        report = CCP4File.CDataFile.assertSame(self, other, diagnostic=diagnostic, testChecksum=True)
        if not(len(report) == 1 and report[0]['code'] == 308):
            return report
        # If we've got a failed checksum test this may be due to trivia in the file
        # Try removing the first line 'HEADER' with date put in by Refmac
        # Other kludges may be necessary for files from other sources
        try:
            retest = True
            f1 = tempfile.mkstemp()
            f2 = tempfile.mkstemp()
            if diagnostic:
                print('CPdbDataFile.assertSame comparing temp files', f1[1], f2[1])
            for obj, fileObj in [[self, f1], [other, f2]]:
                text = CCP4Utils.readFile(obj.__str__())
                if text[0:6] == 'HEADER':
                    if sys.version_info > (3,0):
                        os.write(fileObj[0], text.split('\n', 1)[1].encode("utf-8"))
                    else:
                        os.write(fileObj[0], text.split('\n', 1)[1])
                else:
                    retest = False
                os.close(fileObj[0])
            if retest:
                # Beware need to set object name to get it printed out in report
                obj1 = CPdbDataFile(name=self.objectPath(False), fullPath=f1[1])
                obj2 = CPdbDataFile(fullPath=f2[1])
                otherSum = obj2.checksum()
                selfSum = obj1.checksum()
                if otherSum != selfSum:
                    report.append(self.__class__, 308, name=self.objectPath(False), details=str(self) + ' : ' + str(other))
            else:
                return report
        except:
            return report

    def fixFile(self, xyzout=None, jobId=None, overwrite=False): # unused. KJS
        self.runCoord_format(xyzout=xyzout, jobId=jobId, overwrite=overwrite)

    def runCoord_format(self, xyzout=None, outputFormat=None, jobId=None, overwrite=False):
        from core import CCP4Modules
        from core import CCP4Utils
        if xyzout is None:
            xyzout = self.importFileName(jobId=jobId)
        if jobId is not None:
            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
            logFile = os.path.normpath(os.path.join(jobDirectory, self.objectName() + '_coord_format.log'))
        else:
            logFile = None
        arglist = ['xyzin', self.__str__()]
        arglist.extend(['xyzout', xyzout])
        if outputFormat is not None:
            com = 'OUTPUT ' + outputFormat + '''\nFIXBLANK\nEND\n'''
        else:
            com = '''FIXBLANK\nEND\n'''
        pid = CCP4Modules.PROCESSMANAGER().startProcess('coord_format', arglist, inputText=com, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status != 0:
            return CErrorReport(self.__class__, 401, 'Exit status:' + str(status), name=self.objectPath(), stack=False)
        diffs = CCP4Utils.nonWhiteDifferences(self.__str__(), xyzout)
        if len(diffs) > 0:
            if overwrite:
                os.remove(xyzin)    # KJS : Something wrong here, can't see xyzin anywhere...
                os.rename(xyzout, xyzin)
            ret =  CErrorReport(self.__class__, 402, xyzout, name=self.objectPath(), stack=False)
            for diff in diffs:
                if diff[0] < 0:
                    ret.append(self.__class__, 403, diff[1] + '\nChange made ' + str(diff[2]) + ' times', stack=False)
                else:
                    ret.append(self.__class__, 404, diff[1], stack=False)
            return ret
        else:
            return CErrorReport()

    def importFile(self, jobId=None, sourceFileName=None, ext=None, annotation=None, jobNumber=None):
        from core import CCP4Utils
        if sourceFileName is None:
            sourceFileName = self.__str__()
        if annotation is None:
            annotation= self.qualifiers('guiLabel') + ' imported from ' + os.path.split(sourceFileName)[1]
            if jobNumber is not None:
                annotation += ' by job ' + str(jobNumber)
        if ext is None:
            ext=os.path.splitext(self.__str__())[1]
        if 'sourceFileName' in self.__dict__:
            del self.__dict__['sourceFileName']
        if os.path.splitext(sourceFileName)[1] in ['.pdb', '.ent']:
            text = CCP4Utils.readFile(sourceFileName)
            if text.count('\nATOM ') + text.count('\nHETATM ') == 0:
                return CErrorReport(self.__class__, 405, stack=False)
            filename = self.importFileName(jobId=jobId, ext=ext)
            #print 'CPdbDataFile.importFile copy', sourceFileName, filename
            err = self.runCoord_format(xyzout=filename)
            if len(err) == 1 and err[0]['code'] == 401:
                shutil.copyfile(sourceFileName, filename)
            self.blockSignals(True)
            self.setFullPath(filename)
            if annotation is not None:
                self.annotation.set(annotation)
            self.__dict__['sourceFileName'] = sourceFileName
            self.blockSignals(False)
            return err
        return CErrorReport()

    def convertFormat(self, toFormat, fileName):
        if not self.exists():
            raise CException(self.__class__, 410, self.__str__())
        if self.fileContent.__dict__.get('_molHnd', None) is None:
            self.fileContent.loadFile()
        if self.fileContent.__dict__.get('_molHnd', None) is None:
            raise CException(self.__class__, 411, self.__str__())
        if os.path.exists(fileName):
            try:
                os.remove(fileName)
            except:
                raise CException(self.__class__, 412, fileName)
        if toFormat.lower().count('cif'):
            self.fileContent.__dict__['_molHnd'].WriteCIFASCII(fileName)
        else:
            self.fileContent.__dict__['_molHnd'].WritePDBASCII(fileName)
        if not os.path.exists(fileName):
            raise CException(self.__class__, 413, fileName)

    def removeDummyAtoms(self, cleanedFileName=None, otherFileName=None, excludeResNames='dummy', excludeAtomNames='dummy'):
        ''' Remove 'DUM' or other problem residue names or atom names and optionally save to otheFileName '''
        if not self.exists():
            return 0
        import ccp4mg
        import mmdb2 as mmdb
        import mmut
        import copy
        molHnd = mmdb.Manager()
        molHnd.ReadCoorFile(str(self))
        print('Remove dummy atoms - file contains', molHnd.GetNumberOfAtoms(), 'atoms')
        selHnd = molHnd.NewSelection()
        selAtoms = mmdb.newPPCAtom()
        if excludeResNames == 'dummy':
            excludeResNames = self.fileContent.composition.resTypeSelCommand('dummy')
        print('Remove dummy atoms removing', excludeResNames, 'residues')
        molHnd.Select(selHnd, mmdb.STYPE_ATOM, 1, '*', \
                      mmdb.ANY_RES, '*', mmdb.ANY_RES, '*',excludeResNames, '*', '*', '*', mmdb.SKEY_NEW)
        if excludeAtomNames.count('dummy'):
            excludeAtomNames = self.fileContent.composition.resTypeSelCommand('dummy_atoms')
            print('Remove dummy atoms removing', excludeAtomNames, 'atoms')
        molHnd.Select(selHnd, mmdb.STYPE_ATOM, 1, '*',
                      mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', excludeAtomNames, '*', '*', mmdb.SKEY_OR)
        #nSelAtoms = molHnd.GetSelIndex(selHnd, selAtoms)
        selindexp = mmut.intp()
        selAtoms = mmut.GetAtomSelIndex(molHnd, selHnd, selindexp)
        nSelAtoms = selindexp.value()
        print('Found',nSelAtoms,'dummy atoms')
        exSelAtoms = copy.deepcopy(nSelAtoms)
        if nSelAtoms == 0:
            return 0
        #--------------------------------------------------------------------
        def saveFile(fileName):
            print('Saving to', fileName)
            thisMolHnd = mmdb.Manager()
            thisMolHnd.Copy(molHnd ,mmdb.MMDBFCM_Title | mmdb.MMDBFCM_Cryst)
            if molHnd.GetNumberOfModels() > 1:
                molHnd.CopySelection(selHnd, thisMolHnd)
            else:
                try:
                    for i in range(0, nSelAtoms):
                        pcat = mmdb.getPCAtom(selAtoms, i)
                        thisMolHnd.PutAtom(i+1, pcat) 
                except:
                    pass
                    '''
                    selAtoms = mmdb.newPPCAtom()
                    nSelAtoms = molHnd.GetSelIndex(selHnd,selAtoms)
                    for i in range(0,nSelAtoms):
                        pcat = mmdb.CAtomPtr(mmdb.getPCAtom(selAtoms,i))
                        thisMolHnd.PutAtom ( i+1,pcat )
                    '''
            if os.path.splitext(fileName)[1].count('cif'):
                RC = thisMolHnd.WriteCIFASCII(fileName)
            else:
                RC = thisMolHnd.WritePDBASCII(fileName)
        #--------------------------------------------------------------------
        if otherFileName is not None:
            saveFile(otherFileName)
        if cleanedFileName is not None:
            if cleanedFileName == 'overwrite':
                cleanedFileName = str(self)
            molHnd.Select(selHnd, mmdb.STYPE_ATOM, 1, '*',
                          mmdb.ANY_RES, '*', mmdb.ANY_RES, '*', '*', '*', '*', '*', mmdb.SKEY_XOR)
            #nSelAtoms = molHnd.GetSelIndex(selHnd, selAtoms)
            selindexp = mmut.intp()
            selAtoms = mmut.GetAtomSelIndex(molHnd, selHnd, selindexp)
            nSelAtoms = selindexp.value()
            saveFile(cleanedFileName)
        molHnd.DeleteSelection(selHnd)
        return exSelAtoms


class CPdbDataFileList(CCP4Data.CList):
    SUBITEM = {'class' : CPdbDataFile, 'qualifiers' : {'allowUndefined' : False, 'mustExist': True}}


class CPdbEnsembleItem(CCP4Data.CData):
    CONTENTS = {'structure' : {'class' : CPdbDataFile, 'qualifiers' : {'allowUndefined' : False, 'mustExist': True, 'fromPreviousJob' : True, 'ifAtomSelection' : True}},
                'identity_to_target' : {'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0, 'max': 1.0}},
                'rms_to_target' : {'class' : CCP4Data.CFloat, 'qualifiers' : {'min' : 0.0, 'max' : 100.0}}}
    CONTENTS_ORDER = ['structure', 'identity_to_target', 'rms_to_target']
    ERROR_CODES = {101 : {'description' : 'No sequence identity or structure RMS to target set'}}
    QUALIFIERS = {'guiLabel' : 'Structure in ensemble',
                  'toolTip' : 'Homologous model and its similarity to the target structure',
                  'allowUndefined' : False }

    def validity(self, arg):
        err = CErrorReport()
        allowUndefined = self.qualifiers('allowUndefined')
        #print 'CPdbEnsembleItem.validity', allowUndefined
        e = self.structure.validity(arg.get('structure', {}))
        if arg.get('structure', {}).get('baseName', None) is not None and arg.get('identity_to_target', None) is None and arg.get('rms_to_target', None) is None:
            err.append(self.__class__, 101, name=self.objectPath())
        else:
            err1 = self.identity_to_target.validity(arg.get('identity_to_target', None))
            err2 = self.rms_to_target.validity(arg.get('rms_to_target', None))
            if err1.maxSeverity() < err2.maxSeverity():
                err.extend(err1)
            else:
                err.extend(err2)
        #print 'CPdbEnsembleItem.validity', self.objectPath(), err.report(), arg
        return err

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        return (self.__dict__['_value']['structure'].isSet(allowUndefined=allowUndefined, allowDefault=allowDefault) and \
                (self.__dict__['_value']['identity_to_target'].isSet() or self.__dict__['_value']['rms_to_target'].isSet()))

    # Support Qt QAbstractItemModel
    def child(self, row):
        return None

    def childCount(self):
        return 0

    def columnCount(self):
        return 4

    def data(self, column, role):
        if role == QtCore.Qt.DisplayRole:
            if column == 0:
                if self.structure.isSet():
                    if self.structure.annotation.isSet():
                        return self.structure.annotation.__str__()
                    elif self.structure.baseName.isSet() and len(self.structure.baseName.__str__()) > 0:
                        return self.structure.baseName.__str__()
                return '--'
            elif column == 1:
                if self.__dict__['_value']['structure'].selection.isSet():
                    return self.__dict__['_value']['structure'].selection.__str__()
            elif column == 2:
                if self.__dict__['_value']['identity_to_target'].isSet():
                    return self.__dict__['_value']['identity_to_target'].__str__()
            elif column == 3:
                if self.__dict__['_value']['rms_to_target'].isSet():
                    return self.__dict__['_value']['rms_to_target'].__str__()
        elif role == QtCore.Qt.BackgroundRole:
            from qtgui import CCP4StyleSheet
            if self.__dict__.get('currentItem', False):
                return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.HIGHLIGHTCOLOUR))
        elif role ==  QtCore.Qt.UserRole:
            #print 'CPdbEnsembleItem.data UserRole',str(repr(self))
            return str(repr(self))
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def abstractModelParent(self):
        return self.parent().parent()

    def row(self):
        #print 'CPdbEnsembleItem.row',self.parent().__dict__['_value'],repr(self)
        return self.parent().__dict__['_value'].index(self)


class CEnsemblePdbDataFile(CPdbDataFile):
    '''A PDB coordinate file containing ensemble of structures as 'NMR' models'''
    QUALIFIERS = {}
    QUALIFIERS.update(CPdbDataFile.QUALIFIERS)
    QUALIFIERS.update({'mimeTypeName' : 'chemical/x-pdb', 'mimeTypeDescription' : 'Model coordinates',
                       'fileExtensions' : ['pdb','cif','mmcif','ent'], 'fileContentClassName' : 'CPdbData',
                       'fileLabel' : 'ensemble coordinates', 'guiLabel': 'Model ensemble',
                       'toolTip' : 'An ensemble of model coordinates in PDB or mmCIF format',
                       'downloadModes' : [], 'helpFile' : 'model_data#ensemble_coordinate_files' } )


class CEnsemble(CCP4Data.CData):
    '''An ensemble of models. Typically, this would be a set of related
     PDB files, but models could also be xtal or EM maps. This should
     be indicated by the types entry.
     A single ensemble is a CList of structures.'''
    CONTENTS = {'label' : {'class' : CCP4Data.COneWord},
                'number' : {'class' : CCP4Data.CInt, 'qualifiers' : {'min' : 0, 'default' : 1, 'enumerators' : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], 'menuText' : ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11','12']}},
                'use' : {'class' : CCP4Data.CBoolean, 'qualifiers' : {'default' : True}},
                'pdbItemList' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CPdbEnsembleItem} , 'qualifiers' : {'listMinLength' : 1}}}
    QUALIFIERS = {'guiLabel' : 'Ensemble', 'allowUndefined' : False}

    def greyOut(self):
        return not bool(self.use)

    # Support Qt QAbstractItemModel
    def data(self, column, role):
        if role == QtCore.Qt.DisplayRole:
            if column == 0:
                return self.__dict__['_value']['number'].__str__()+' x '+self.__dict__['_value']['label'].__str__()
        elif role ==  QtCore.Qt.BackgroundRole:
            if self.__dict__.get('currentItem', False):
                from qtgui import CCP4StyleSheet
                return QtGui.QBrush(QtGui.QColor(CCP4StyleSheet.HIGHLIGHTCOLOUR))
        elif role == QtCore.Qt.ForegroundRole:
            if self.greyOut():
                return QtGui.QBrush(QtCore.Qt.gray)
        elif role == QtCore.Qt.UserRole:
            #print 'CEnsemble.data UserRole',str(repr(self))
            return str(repr(self))
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def childListObject(self):
        return self.__dict__['_value']['pdbItemList']

    def child(self,row):
        return self.__dict__['_value']['pdbItemList'].__getitem__(row)

    def childCount(self):
        return self.__dict__['_value']['pdbItemList'].__len__()

    def columnCount(self):
        return 1

    def abstractModelParent(self):
        #return self.parent().__dict__['_abstractModelParent']
        return QtCore.QModelIndex()

    def row(self):
        return self.parent().__dict__['_value'].index(self)

    '''
    def isSet(self,allowUndefined=False,allowDefault=True,allSet=True):
        if allSet: return CData.isSet(self,allowUndefined=allowUndefined,allowDefault=allowDefault,allSet=True)
        for key in ['pdbItemList']:
          if self._value[key].isSet(allowUndefined=allowUndefined,allowDefault=allowDefault,allSet=allSet):
             return True
        return False
    '''


class CEnsembleList(CCP4Data.CList):

    QUALIFIERS = {'listMinLength' : 1}
    SUBITEM = {'class' : CEnsemble}

    def __init__(self, *args, **kws):
        super(CEnsembleList, self).__init__(*args, **kws)
        self[-1].label.set(self.firstFreeEnsembleName(-1))
        #print self[-1]

    def validity(self, arg):
        #print 'CEnsembleList.validity',len(self.__dict__['_value'][0].pdbItemList),self.__dict__['_value'][0].pdbItemList.isSet()
        if self.qualifiers('listMinLength') == 0 and len(self.__dict__['_value']) == 1 and not self.__dict__['_value'][0].pdbItemList.isSet():
            return CErrorReport()
        else:
            return CCP4Data.CList.validity(self, arg)

    def firstFreeEnsembleName(self, index=-1):
        if index < 0:
            index = self.__len__() - 1
        name = 'Ensemble_'+str(index + 1)
        while not self.isUniqueName(name):
            index += 1
            name = 'Ensemble_' + str(index + 1)
        return name

    def addItem(self, value={}, index=-1):
        obj = CCP4Data.CList.addItem(self, value=NotImplemented, index=index)
        obj.label.set(self.firstFreeEnsembleName(index))
        return obj

    def isUniqueName(self, name):
        for item in self.__dict__['_value']:
            if item.label == name:
                return False
        return True

    def saveToDb(self):
        saveList = []
        for obj0 in self.__dict__['_value']:
            for obj1 in obj0.pdbItemList:
                if obj1.structure.exists():
                    obj1.structure.subType = 2
                    saveList.append(obj1.structure)
        #print 'CEnsembleList.saveToDb',saveList
        return saveList, None, {}

class CAtomRefmacSelection(CCP4Data.CData):
    '''A residue range selection for rigid body groups'''

    CONTENTS = {'groupId' : {'class' : CCP4Data.CInt},
                'chainId' : {'class' : CCP4Data.COneWord},
                'firstRes' : {'class' : CCP4Data.CInt},
                'lastRes' : {'class' : CCP4Data.CInt}}
    CONTENTS_ORDER = ['groupId', 'chainId', 'firstRes', 'lastRes']

class CAtomRefmacSelectionOccupancy(CCP4Data.CData):
    '''A residue range selection for occupancy groups'''

    CONTENTS = {'groupId' : {'class' : CCP4Data.CInt},
                'chainIds' : {'class' : CCP4Data.CString},
                'firstRes' : {'class' : CCP4Data.CInt},
                'lastRes' : {'class' : CCP4Data.CInt},
                'atoms' : {'class' : CCP4Data.CString},
                'alt' : {'class' : CCP4Data.COneWord}}
    CONTENTS_ORDER = ['groupId', 'chainIds', 'firstRes', 'lastRes', 'atoms', 'alt']

class CAtomRefmacSelectionGroups(CCP4Data.CData):
    '''A group selection for occupancy groups'''
    CONTENTS = {'groupIds' : {'class' : CCP4Data.CString}}
    CONTENTS_ORDER = ['groupIds']

class CResidueRange(CCP4Data.CData):
    '''A residue range selection'''

    CONTENTS = {'chainId' : {'class' : CCP4Data.COneWord, 'qualifiers' : {'default' : ''}},
                'firstRes' : {'class' : CCP4Data.COneWord}, 'lastRes' : {'class' : CCP4Data.COneWord}}
    CONTENTS_ORDER = ['chainId', 'firstRes', 'lastRes']
    QUALIFIERS = {'pdbFileKey' : None}
    QUALIFIERS_DEFINITION = {'pdbFileKey' : {'type' : str,
                                             'description' : 'The key for a CPdbDataFile in the same CContainer'}}
    QUALIFIERS_ORDER = ['pdbFileKey']

    def getFileContent(self):
        pdbDataFile = self.getDataByKey('pdbFileKey')
        if pdbDataFile is None or pdbDataFile.fileContent is None:
            return None
        else:
            return pdbDataFile.fileContent

class CResidueRangeList(CCP4Data.CList):
    '''A list of residue range selections'''
    SUBITEM = {'class' : CResidueRange}

class CSequenceString(CCP4Data.CString):
    ERROR_CODES = {401 : {'description' : 'Non-alphabet character removed from sequence', 'severity' : SEVERITY_WARNING},
                   402 : {'description' : 'Invalid characters (BJOXZ) in sequence'},
                   403 : {'description' : 'Sequence undefined', 'severity' : SEVERITY_WARNING}}
    pass

    # This is used by CAsuContentSeqView.validate() to test validity of sequence
    # but we don't want its strict criteria blocking setting the model object value
    # in CAsuContentSeqView.updateModelFromView()
    def validity0(self, arg):
        err = CErrorReport()
        if arg is None or len(arg.strip()) == 0:
            err.append(self.__class__, 403)
            return err
        seq1 = re.sub(r'[BJOXZ]', '', arg)
        #print 'CSequenceString.validity',seq1
        #print 'CSequenceString.validity',len(seq1),len(arg)
        if len(seq1) < len(arg):
            err.append(self.__class__, 402)
        return err


class CAsuContentSeq(CCP4Data.CData):
    CONTENTS = {'sequence' : { 'class' : CSequenceString , 'qualifiers' : { 'allowUndefined' : False, 'minLength' : 1 } },
               'nCopies' : { 'class' : CCP4Data.CInt ,
                             'qualifiers' :  { 'enumerators' : [0,1,2,3,4,5,6,7,8,9,10,11,12], 'default' : 1 , 'min' : 0 }},
               'polymerType' :  { 'class' : CCP4Data.CString,
                        'qualifiers' : { 'onlyEnumerators' : True, 'enumerators' : ['PROTEIN','RNA','DNA'], 'default' : 'PROTEIN' }  },
               'name' : { 'class' : CCP4Data.CString,
                          'qualifiers' : { 'allowUndefined' : False, 'minLength' : 1, 'allowedCharsCode' : 1 } },
               'description' : { 'class' : CCP4Data.CString, 'qualifiers' : { 'allowUndefined' : True } },
               'source' : { 'class' : CCP4File.CDataFile }
               }

    def autoSetPolymerType(self):
        #print("CAsuContentSeq.autoSetPolymerType")
        #print(self.sequence)
        nposs_nuc = self.sequence.count('A') + self.sequence.count('G') + self.sequence.count('C') + self.sequence.count('T') + self.sequence.count('U') + self.sequence.count('N');
        if float(nposs_nuc)/len(self.sequence)>0.9:
            #print("This is possibly nucleic acid")
            if self.sequence.count('U') > 0:
                self.polymerType.set("RNA")
            else:
                self.polymerType.set("DNA")
        else:
            #print("This is probably protein")
            self.polymerType.set("PROTEIN")
        

    def getTableTextItems(self):
        line = []
        if self.name.isSet():
            line.append(str(self.name))
        else:
            line.append('')
        line.append(str(self.nCopies))
        if self.description.isSet():
            line.append(str(self.description))
        else:
            line.append('')
        if self.sequence.isSet():
            line.append(str(self.sequence))
        else:
            line.append('')
        if self.polymerType.isSet():
            line.append(str(self.polymerType))
        else:
            line.append('')
        return line

    def cleanupSequence(self,sequence=None):
        err = CErrorReport()
        seq1 = re.sub(r'[\r\n ]', '', str(sequence))
        seq2 = re.sub(r'[^A-Z]', '', seq1.upper())
        if len(seq2) < len(seq1):
            err.append(CSequenceString, 401)
        seq3 = re.sub(r'[BJOXZ]', '', seq2)
        if len(seq3) < len(seq2):
            err.append(CSequenceString, 402)
        #print 'cleanupSequence',sequence
        #print 'cleanupSequence3',seq3
        return seq3, err

    def formattedSequence(self):
        nGap = 10
        nLine = 60
        seqList = []
        import copy
        text = copy.deepcopy(str(self.sequence))
        while len(text) > 0:
            line = text[0:nLine]
            text = text[nLine:]
            while len(line) > 0:
                seqList.append(line[0:nGap] + ' ')
                line = line[nGap:]
            seqList[-1] = seqList[-1] + '\n'
        newText = ''
        newText = newText.join(seqList)
        return newText

    def molecularWeight(self,polymerType="PROTEIN"):
        import Bio.SeqUtils
        # Beware BioPython not robust to bad sequences
        if polymerType == "PROTEIN":
            seq = re.sub('[^GALMFWKQESPVICYHRNDT]', '', str(self.__dict__['_value']['sequence']))
            wt = Bio.SeqUtils.molecular_weight(seq,seq_type='protein')
        else:
            seq = re.sub('[^CAUGT]', '', str(self.__dict__['_value']['sequence']))
            wt = Bio.SeqUtils.molecular_weight(seq,seq_type=polymerType)
        return wt * float(self.__dict__['_value']['nCopies'])

    def numberOfResidues(self, countMulti=False):
        if not self.sequence.isSet():
            return 0
        nRes = len(self.sequence)
        if countMulti and self.nCopies.isSet():
            nRes = nRes * int(self.nCopies)
        return nRes


class CAtomRefmacSelectionList(CCP4Data.CList):
    SUBITEM = {'class' : CAtomRefmacSelection}
    def validity(self,arg):
        return CCP4Data.CList.validity(self, arg)

    def __init__(self, *args, **kws):
        super(CAtomRefmacSelectionList, self).__init__(*args, **kws)

class COccRefmacSelectionList(CCP4Data.CList):
    SUBITEM = {'class' : CAtomRefmacSelectionOccupancy}
    def validity(self,arg):
        return CCP4Data.CList.validity(self, arg)

    def __init__(self, *args, **kws):
        super(COccRefmacSelectionList, self).__init__(*args, **kws)

class COccRelationRefmacList(CCP4Data.CList):
    SUBITEM = {'class' : CAtomRefmacSelectionGroups}
    def validity(self,arg):
        return CCP4Data.CList.validity(self, arg)

    def __init__(self, *args, **kws):
        super(COccRelationRefmacList, self).__init__(*args, **kws)

class CAsuContentSeqList(CCP4Data.CList):
    SUBITEM = {'class' : CAsuContentSeq}
    ERROR_CODES = {401 : {'description' : 'Sequence the same as a sequence that is already loaded'},
                   402 : {'description' : 'Sequence names are not unique: '}}
    QUALIFIERS = {'listMinLength' : 0}

    def matchingSequence(self, compSeq, name=None):
        rv = CErrorReport()
        for idx in range(self.__len__()):
            cleanSeq, err = self[idx].cleanupSequence(compSeq)
            if self[idx].sequence == cleanSeq:
                rv.append(self.__class__, 401, str(name) + ' same as ' + str(self[idx].name))
                return rv
        return rv

    '''
    def loadSequenceFile(self,record=0,fileId=None):
        fileAnnotation = None
        if fileId is not None:
          from core import CCP4Modules
          fileName = CCP4Modules.PROJECTSMANAGER().db().getFullPath(fileId=fileId)
          fileAnnotation = CCP4Modules.PROJECTSMANAGER().db().getFileInfo(fileId=fileId,mode='annotation')
          #print 'CAsuContentSeqList.loadSequenceFile fileAnnotation',fileAnnotation
          seqFile = CCP4ModelData.CSeqDataFile(parent=self)
          seqFile.setFullPath(fileName)
          seqFile.__dict__['format'] = 'internal'
          seqFile.fileContent.loadInternalFile(str(self.seqFile))
        else:
          #print 'CAsuContentSeqList.loadSequenceFile to loadExternalFile',seqFile
          seqFile.fileContent.loadExternalFile(str(seqFile),seqFile.__dict__['format'],record=record)
          #self.model.source = self.seqFile.__str__()
        cleanSeq,err = self[0].cleanupSequence(seqFile.fileContent.sequence)
        if err.maxSeverity()>SEVERITY_WARNING:
          raise err
        self.model.sequence.set(cleanSeq)
        try:
          self.model.name.set(self.model.name.fix(self.seqFile.fileContent.name))
        except:
          pass
        if not self.model.name.isSet():
          if fileAnnotation is not None:
            self.model.name.set(self.model.name.fix(fileAnnotation))
          else:
            import os
            self.model.name.set(os.path.split(os.path.splitext(self.seqFile.__str__())[0])[1])
        self.model.description.set(self.seqFile.fileContent.description)
        self.extendSeqList(self.seqFile,recordList[1:])
    '''

    def extendSeqList(self, fileObject, recordList=[]):
        '''
        Reading sequence file into CAsuContentSeqView user seelcts multiple sequences
        so need to create and load more CAsuContentSeq to self.seqList
        '''
        ret = CErrorReport()
        if isinstance(fileObject, CSeqDataFile):
            for record in recordList:
                fileObject.fileContent.loadExternalFile(fileObject.__str__(), fileObject.__dict__['format'], record=record)
                rv =self.matchingSequence(str(fileObject.fileContent.sequence), str(fileObject.fileContent.name))
                if rv.maxSeverity() <= SEVERITY_WARNING: 
                    self.addItem()
                    self[-1].name.set(self[-1].name.fix(str(fileObject.fileContent.name)))
                    self[-1].description.set(fileObject.fileContent.description)
                    cleanSeq, err = self[-1].cleanupSequence(fileObject.fileContent.sequence)
                    self[-1].sequence.set(cleanSeq)
                    if hasattr(self[-1],"autoSetPolymerType"):
                        self[-1].autoSetPolymerType()
                else:
                    ret.append(rv)
        else:
            for record in recordList:
                chainId = list(fileObject.fileContent.sequences.keys())[record]
                rv =self.matchingSequence(fileObject.fileContent.sequences[chainId], chainId)
                if rv.maxSeverity() <= SEVERITY_WARNING:
                    self.addItem()
                    #print 'CAsuContentSeqList.extendSeqList chainId',chainId
                    self[-1].name.set(chainId)
                    cleanSeq, err = self[-1].cleanupSequence(fileObject.fileContent.sequences[chainId])
                    self[-1].sequence.set(cleanSeq)
                    if hasattr(self[-1],"autoSetPolymerType"):
                        self[-1].autoSetPolymerType()
                else:
                    ret.append(rv)
        self.dataChanged.emit()
        return rv

    def validity(self,value):
        rv = CCP4Data.CList.validity(self, value)
        for idx in range(1, len(value)):
            #print 'CAsuContentSeqList.validity', idx,value[idx], type(value[idx])
            for jdx in range(0, idx):
                if str(value[idx].get('name')) == str(value[jdx].get('name')):
                    rv.append(self.__class__, 402, value[idx].get('name'))
        return rv

    def importSourceFiles(self):
        pass

    def nameList(self):
        ret = []
        for seqObj in self.__dict__['_value']:
            ret.append(str(seqObj.name))
        return ret

    def molecularWeight(self):
        totWeight = 0.
        for seqObj in self:
            totWeight = totWeight + seqObj.molecularWeight(seqObj.polymerType)
        return totWeight

class CAsuContent(CCP4File.CDataFileContent):

    CONTENTS = {'seqList' : {'class' : CAsuContentSeqList}}
                #  'exptlMolWt' : {'class' : CCP4Data.CFloat, 'qualifiers' : {'min' : 0.0, 'allowUndefined' : True}}
    ERROR_CODES = {101 : {'description' : 'Failed reading file - is it correct file type?'},
                   102 : {'description' : 'Failed reading file - it is not AU contents file'}}

    def loadFile(self, fileName=None):
        self.unSet()
        fileName = str(fileName)
        if fileName is None or not os.path.exists(fileName):
            return
        xmlFileObject = CCP4File.CI2XmlDataFile(fileName)
        err = xmlFileObject.loadHeader()
        if err.maxSeverity() > SEVERITY_WARNING:
            raise CErrorReport(self.__class__, 101, str(fileName), name=self.objectPath())
        if xmlFileObject.header.function != 'ASUCONTENT':
            raise CErrorReport(self.__class__, 102, str(fileName), name=self.objectPath())
        root = xmlFileObject.getBodyEtree()
        self.setEtree(root)

    def saveFile(self, fileName=None, jobDetails={}):
        #print 'CAsuContent.saveFile', fileName, jobDetails
        self.seqList.importSourceFiles()
        xmlFileObject = CCP4File.CI2XmlDataFile(fileName)
        xmlFileObject.header.setCurrent()
        xmlFileObject.header.function.set('ASUCONTENT')
        for item in ['projectName', 'projectId', 'jobId', 'jobNumber', 'comment']:
            if item in jobDetails and jobDetails[item] is not None:
                xmlFileObject.header.get(item).set(jobDetails[item])
        xmlFileObject.saveFile(self.getEtree())

    def getChainNames(self):
        chainNamesList = []
        for item in self.seqList:
            chainNamesList.append(str(item.name))
        return chainNamesList


    def molecularWeight(self):
        return self.seqList.molecularWeight()

class CAsuDataFile(CCP4File.CI2XmlDataFile):
    CONTENTS = {}
    CONTENTS.update(CCP4File.CDataFile.CONTENTS)
    CONTENTS['selection'] = {'class' : CCP4Data.CDict, 'subItem' : {'class' : CCP4Data.CBoolean}}
    CONTENTS_ORDER = []
    CONTENTS_ORDER.extend(CCP4File.CDataFile.CONTENTS_ORDER)
    CONTENTS_ORDER.append('selection')
    QUALIFIERS = {'mimeTypeName' : 'application/CCP4-asu-content', 'mimeTypeDescription' : 'AU content',
                  'fileExtensions' : ['asu.xml'], 'fileContentClassName' : 'CAsuContent',
                  'fileLabel' : 'AU contents', 'guiLabel' : 'AU contents',
                  'toolTip' : 'A CCP4i2 file specifying AU contents',
                  'helpFile' : 'model_data#sequences', 'saveToDb' : True, 'selectionMode' : 0}
    QUALIFIERS_DEFINITION = {'selectionMode' : {'type' :int, 'description' : 'Chain selection options'}}
    #selectionMode = 0 - not chain selection
    #selectionMode = 1 - select one chain
    #selectionMode = 2 - select multi chains

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, fullPath=None, keywords={}, **kw):
        CCP4File.CI2XmlDataFile.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, fullPath=fullPath, keywords=keywords, **kw)

    def updateData(self):
        if self.qualifiers('selectionMode') > 0:
            #print 'CAsuDataFile.updateData unSet',self.selection,self.__dict__['_value']['baseName']
            if self.__dict__['_value']['baseName'].isSet():
                self.loadFile()
                self.buildSelection()
            else:
                self.selection.clear()
        CCP4File.CI2XmlDataFile.updateData(self)
  
    def buildSelection(self):
        ok = True
        #print 'buildSelection selection',repr(self.selection),self.selection.keys(),'seqList',self.fileContent.seqList.nameList()
        if len(self.fileContent.seqList) == len(self.selection):
            for idx in range(len(self.fileContent.seqList)):
                if self.fileContent.seqList[idx].name not in self.selection:
                    ok = False
        else:
            ok = False
        if ok:
            return
        selDict = {}
        status = True
        for idx in range(len(self.fileContent.seqList)):
            selDict[str(self.fileContent.seqList[idx].name)] = status
            if self.qualifiers('selectionMode') == 1:
                status = False
        #print 'buildSelection selDict',selDict
        self.selection.set(selDict)

    def saveFile(self,jobDetails={}):
        self.fileContent.saveFile(self.__str__(), jobDetails=jobDetails)

    def loadFile(self):
        # Circumvent the CI2XmlDataFile.loadFile() which throws complete wobbler of non-existant file
        CCP4File.CDataFile.loadFile(self)
        #print 'from CAsuDataFile.loadFile',self.selection

    def writeFasta(self, fileName, indx=-1, format='fasta', writeMulti=False, polymerTypes=["PROTEIN", "RNA", "DNA"]):
        # Write a fasta file or a bad pir file with extra blank line but not other
        # proper pir features
        from core import CCP4Utils
        self.loadFile()
        selectionMode = self.qualifiers('selectionMode')
        if indx < 0:
            text = ''
            for seqObj in self.fileContent.seqList:
                if seqObj.polymerType not in polymerTypes:
                    continue
                if writeMulti:
                    nCopies = int(seqObj.nCopies)
                else:
                    nCopies = min(1,int(seqObj.nCopies))
                for nC in range(nCopies):
                    name = seqObj.name.__str__()
                    if selectionMode == 0 or (not self.selection.isSet()) or self.selection[name]:
                        text += '>' + name + '\n'
                        if format == 'pir':
                            text += '\n'
                        seq = seqObj.sequence.__str__()
                        while len(seq) > 0:
                            text += seq[0:60] + '\n'
                            seq = seq[60:]
        else:
            seqObj = self.fileContent.seqList[indx]
            text = '>' + seqObj.name.__str__() + '\n'
            if format == 'pir':
                text += '\n'
            seq = seqObj.sequence.__str__()
            while len(seq) > 0:
                text += seq[0:60] + '\n'
                seq = seq[60:]
        CCP4Utils.saveFile(fileName, text)

    def writeArpPir(self, fileName, indx=-1, writeMulti=False):
        from core import CCP4Utils
        self.loadFile()
        selectionMode = self.qualifiers('selectionMode')
        if indx < 0:
            text = ''
            for seqObj in self.fileContent.seqList:
                if writeMulti:
                    nCopies = int(seqObj.nCopies)
                else:
                    nCopies = min(1, int(seqObj.nCopies))
                for nC in range(nCopies):
                    name = seqObj.name.__str__()
                    if selectionMode == 0 or (not self.selection.isSet()) or self.selection[name]:
                        if len(text) > 0:
                            # Use multiple A residues as chain separator    # KJS : Fix this
                            #text += 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
                            text += '> xyz\n'
                        else:
                            text += '>' + name + '\n\n'
                        seq = seqObj.sequence.__str__()
                        while len(seq) > 0:
                            text += seq[0:60] + '\n'
                            seq = seq[60:]
        else:
            seqObj = self.fileContent.seqList[indx]
            text = '>' + seqObj.name.__str__() + '\n'
            if format == 'arp':
                text += '\n'
            seq = seqObj.sequence.__str__()
            while len(seq) > 0:
                text += seq[0:60] + '\n'
                seq = seq[60:]
        CCP4Utils.saveFile(fileName, text)

    def molecularWeight(self, ignoreSelection=False):
        if not self.isSet():
            return 0.0
        selectionMode = self.qualifiers('selectionMode')
        totWeight = 0.0
        for seqObj in self.fileContent.seqList:
            if ignoreSelection or selectionMode == 0 or self.selection[seqObj.name]:
                totWeight = totWeight + seqObj.molecularWeight(seqObj.polymerType)
        return totWeight

    def numberOfResidues(self, countMulti=False, ignoreSelection=False):
        if not self.isSet():
            return 0
        selectionMode = self.qualifiers('selectionMode')
        nTot = 0
        for seqObj in self.fileContent.seqList:
            if ignoreSelection or selectionMode == 0 or self.selection[seqObj.name]:
                nTot += seqObj.numberOfResidues(countMulti=countMulti)
        return nTot

    def isSelected(self, seqObj=None):
        name = None
        if seqObj is not None:
            name = str(seqObj.name)
        if name is None:
            return False
        elif self.selection.get(name, None) is None:
            return False
        else:
            return bool(self.selection[name])


class CPairwiseAlignment:

    def __init__(self, sequence1=None, sequence2=None):
        if isinstance(sequence1, CAsuContentSeq):
            self.seq1 = str(sequence1.sequence)
        elif isinstance(sequence1, (CCP4Data.CString, str)):
            self.seq1 = str(sequence1)
        if isinstance(sequence2, CAsuContentSeq):
            self.seq2 = str(sequence2.sequence)
        elif isinstance(sequence2, (CCP4Data.CString, str)):
            self.seq2 = str(sequence2)

    def align(self):
        import Bio
        from Bio import pairwise2
        alignments = pairwise2.align.globalxx(self.seq1, self.seq2)
        if len(alignments) > 0:
            #print 'align', alignments[0]
            return alignments[0]
        else:
            return None

    def formattedAlignment(self, label1='', label2=''):
        import Bio
        from Bio import pairwise2
        ali = self.align()
        if ali is None:
            return 'No alignment found\n'
        text = pairwise2.format_alignment(*ali)
        labels = [(label1 + '              ')[0:15], '               ', (label2 + '              ')[0:15]]
        lines = text.split('\n')
        print('formattedAlignment', lines)
        ret = ''
        while len(lines[0]) > 0:
            for n in (0, 1, 2):
                ret += labels[n] + '    ' + lines[n][0:60] + '\n'
                lines[n] = lines[n][60:]
            ret += '\n'
        ret += '\n' + lines[3] + ' for ' + str(ali[4]) + ' residues\n'
        #print 'formattedAlignment', text
        return ret

class CChainMatch:

    def __init__(self, model1, model2):
        ''' Find best chain matches and sequence alignment for those chain matches where
        model1 and model2 can be any combination of CPdbDataFile or CAsuDataFile.
        The matching and output format is designed for model1 is CPdbDataFile and model2 is CAsuDataFile '''
        for chain, model in [['chains1', model1], ['chains2', model2]]:
            if isinstance(model, CPdbDataFile):
                setattr(self, chain, CAsuContentSeqList())
                addItem = True  # CAsuContentSeqList created with zero items in list
                for cid, seq in list(model.fileContent.sequences.items()):
                    #print 'CChainMatch',cid,seq
                    if addItem:
                        getattr(self, chain).addItem()
                    addItem = True
                    getattr(self, chain)[-1].name.set(cid)
                    getattr(self, chain)[-1].sequence.set(seq)
                    #print 'CChainMatch', getattr(self, chain)
            elif isinstance(model, CAsuDataFile):
                setattr(self, chain, model.fileContent.seqList)

    def score2(self):
        # Attempts to take into account the number of copies expected for sequence (model2)
        import copy
        # setup the score 'matrix'
        # http://stackoverflow.com/questions/6667201/how-to-define-two-dimensional-array-in-python
        nCh1 = nCh2 = 0
        nCh1 = len(self.chains1)
        for seqObj in self.chains2:
            nCh2 += max(1, int(seqObj.nCopies))
        # A matrix of scores - this has multiple elements if chain has multiple copies
        score = [[0]*nCh2 for i in range(nCh1)]
        # A matrix with same elements as score matrix providing indexing back to the chains
        chIndex = [[None]*nCh2 for i in range(nCh1)]
        print('CChainMatch nCh', nCh1, nCh2, score)
        # Get score for all chains in chains1 v. all chains in chains2
        # Where the chain has multiple copies make multiple row/column
        n1 = 0
        chIndx1 = 0
        for c1 in self.chains1:
            n2 = 0
            chIndx2 = 0
            for c2 in self.chains2:
                align = CPairwiseAlignment(str(c1.sequence), str(c2.sequence))
                ret = align.align()
                #print 'n1, n2, ret[2], ret[4]', n1, n2, ret[2], ret[4]
                if ret is not None:
                    print('CChainMatch ret', ret)
                    print('CChainMatch', n1, n2, score)
                    score[n1][n2] = float(ret[2])/float(ret[4])
                    chIndex[n1][n2] = [chIndx1,chIndx2]
                    print('score', n1, n2, score[n1][n2])
                for i2 in range(int(c2.nCopies) - 1):
                    #print 'Copying n2 score', n1, n2
                    n2 += 1
                    score[n1][n2] = score[n1][n2 - 1]
                    chIndex[n1][n2] = [chIndx1,chIndx2]
                n2 += 1
                chIndx2 += 1
            n1 += 1
            chIndx1 += 1
        #for i1 in range(nCh1): print 'CChainMatch.match', score[i1]
        #for i1 in range(nCh1): print 'CChainMatch.match', chIndex[i1]
        matchList = []
        maxMatch = min(nCh1, nCh2)
        while len(matchList) < maxMatch:
            bestValue = 0.0
            bestPair = [-1, -1]
            for n1 in range(nCh1):
                for n2 in range(nCh2):
                    if score[n1][n2] > bestValue:
                        bestValue = copy.deepcopy(score[n1][n2])
                        bestPair = [n1, n2]
            #print 'bestValue',bestValue,bestPair
            if bestValue > 0.0:
                # add best matched pair of chains to list and 'blank' them out of the score matrix
                matchList.append(chIndex[bestPair[0]][bestPair[1]])
                for i2 in range(nCh2):
                    score[bestPair[0]][i2] = -1.0
                for i1 in range(nCh1):
                    score[i1][bestPair[1]] = -1.0
            else:
                break
        print('matchList', matchList)
        return matchList

    def score(self):
        import copy
        # setup the score 'matrix'
        # http://stackoverflow.com/questions/6667201/how-to-define-two-dimensional-array-in-python
        nCh1 = len(self.chains1)
        nCh2 = len(self.chains2)
        score = [[0]*nCh2 for i in range(nCh1)]
        # Get score for all chains in chains1 v. all chains in chains2
        n1 = 0
        for c1 in self.chains1:
            n2 = 0
            for c2 in self.chains2:
                align = CPairwiseAlignment(str(c1.sequence), str(c2.sequence))
                ret = align.align()
                #print 'n1,n2,ret[2],ret[4]', n1, n2, ret[2], ret[4]
                if ret is not None:
                    print('CChainMatch ret', ret)
                    print('CChainMatch', n1, n2, score)
                    score[n1][n2] = float(ret[2])/float(ret[4])
                    print('score', n1, n2,score[n1][n2])
                n2 += 1
            n1 += 1
        #for i1 in range(nCh1): print 'CChainMatch.match', score[i1]
        #for i1 in range(nCh1): print 'CChainMatch.match', chIndex[i1]
        matchList = []
        for n1 in range(nCh1):
            bestValue = 0.0
            bestIndx = -1
            for n2 in range(nCh2):
                if score[n1][n2] > bestValue:
                    bestValue = copy.deepcopy(score[n1][n2])
                    bestIndx = n2
            #print 'bestValue',bestValue,bestPair
            if bestValue > 0.0:
                # add best matched pair of chains to list and 'blank' them out of the score matrix
                matchList.append((n1, bestIndx))
        print('matchList', matchList)
        return matchList

    def bestAlignments(self):
        matchList = self.score()
        for chIndx1, chIndx2 in matchList:
            #print 'chIndx1,chIndx2',chIndx1,chIndx2
            #print 'matching',self.chains1[chIndx1].name,self.chains2[chIndx2].name
            align = CPairwiseAlignment(str(self.chains1[chIndx1].sequence), str(self.chains2[chIndx2].sequence))
            ret = align.align()
            print(ret)

    def reportXmlAlignments(self):
        from lxml import etree
        matchList = self.score()
        root = etree.Element('ChainMatching')
        for chIndx1, chIndx2 in matchList:
            #print 'chIndx1,chIndx2', chIndx1, chIndx2
            ele = etree.SubElement(root, 'AlignChain')
            align = CPairwiseAlignment(str(self.chains1[chIndx1].sequence), str(self.chains2[chIndx2].sequence))
            e = etree.SubElement(ele, 'ChainId')
            e.text = self.chains1[chIndx1].name.__str__()
            e = etree.SubElement(ele, 'ChainId2')
            e.text = self.chains2[chIndx2].name.__str__()
            e = etree.SubElement(ele, 'Alignment')
            e.text = align.formattedAlignment(str(self.chains1[chIndx1].name), str(self.chains2[chIndx2].name))
        return root

#===========================================================================================================
import unittest

def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testPdbData)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testRange))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testSeqData))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)


class testPdbData(unittest.TestCase):

    def test1(self):
        p = CPdbData()
        p.loadFile(CPdbDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1df7.pdb'))
        self.assertEqual(p.composition.chains,['A'], 'Error loading CPDBData - wrong chains')
        self.assertEqual(p.composition.moleculeType, ['PROTEIN', 'MONOMER'], 'Error loading CPDBData - wrong moleculeType')

    def test2(self):
        p = CPdbData()
        resNum, altLoc = p.interpretResidueInput('123.1')
        self.assertEqual(resNum, 123, 'Error in CPdbData.interpretResidueInput - wrong resNum')
        self.assertEqual(altLoc, '1', 'Error in CPdbData.interpretResidueInput - wrong altLoc')

    def test3(self):
        from core import CCP4Utils
        p = CPdbDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1df7.pdb')
        p.selection = 'A/10-20'
        self.assertEqual(p.isSelectionSet(), True, 'Error - CPdbDataFile.isSelectionSet() wrong')
        p.getSelectedAtomsPdbFile(fileName=os.path.join(CCP4Utils.getTestTmpDir(), '1df7_selection.pdb'))


class testRange(unittest.TestCase):

    def test1(self):
        from core import CCP4Container
        c = CCP4Container.CContainer()
        c.addContent(name='XYZIN', cls=CPdbDataFile)
        c.addContent(name='DOMAINLIST', cls=CCP4Data.CList, subItem={'class':CResidueRangeList } )
        c.DOMAINLIST.append([])
        c.DOMAINLIST[0].append({'chainId':'A', 'firstRes':'1', 'lastRes':'20'})
        c.DOMAINLIST[0].append({'chainId':'B', 'firstRes':'1', 'lastRes':'40'})
        #print 'testRange.test1',c.DOMAINLIST[0].subItemClass()
        #print 'testRange.test1',c.DOMAINLIST[0],type(c.DOMAINLIST[0])
        #print 'testRange.test1',c.DOMAINLIST[0][1],type(c.DOMAINLIST[0][1])
        self.assertEqual(c.DOMAINLIST[0][1].chainId, 'B', 'Failed to set CResidueRangeList data')


class testSeqData(unittest.TestCase):

    def test1(self):
        seq = CSequence()
        seq.loadFile(CSeqDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1mzr.fasta'))
        print(seq.__dict__['_value'])
        seq.loadFile(CSeqDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1mzr.pir'))
        print(seq.__dict__['_value'])

class testCAsuContent(unittest.TestCase):

    def testMolecularWeight(self):
        asuFile = CAsuDataFile(project='CCP4I2_TOP', relPath='demo_data/gamma', baseName='gamma.asu.xml')
        asuFile.loadFile()

    def testCreateCAsuContent(self):
        sequence = CAsuContentSeq({'sequence':'QWERTY', 'nCopies':1, 'polymerType':'PROTEIN', 'name':'DUMMY', 'desccription':'', 'source':''})
        asuContent = CAsuContent()
        asuContent.seqList.addItem(sequence)
        print(asuContent)
        self.assertAlmostEqual(asuContent.molecularWeight(), second=881.93, places=2)
