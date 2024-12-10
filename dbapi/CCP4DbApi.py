from __future__ import print_function

"""
    CCP4DbApi.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

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
   Initially based on code from by George Pelios at Diamond - version from 15 July 2010
   Liz Potterton Mar 2011 - Imported into CCP4i2 and radically rewritten
"""

import sqlite3
import os
import re
import types
import sys
import traceback
import copy
if sys.version_info >= (3,7):
    from collections.abc import Callable
else:
    from collections import Callable

from PySide2 import QtCore, QtSql

from core import CCP4Config
from core.CCP4ErrorHandling import *
import collections
from core.CCP4QtObject import CObject

def isAlive(qobj):
  return True


USE_PERFORMANCE_CLASSES = True

UUIDTYPE = str

PRIVILEGE_NONE = 0
PRIVILEGE_READ = 1
PRIVILEGE_COMMENT = 2
PRIVILEGE_WRITE = 3
PRIVILEGE_DELETE = 4
PRIVILEGE_PROJECT = 5
PRIVILEGE_EXTEND = 6

JOB_STATUS_UNKNOWN = 0
JOB_STATUS_PENDING = 1
JOB_STATUS_QUEUED = 2
JOB_STATUS_RUNNING = 3
JOB_STATUS_INTERRUPTED = 4
JOB_STATUS_FAILED = 5
JOB_STATUS_FINISHED = 6
JOB_STATUS_REMOTE = 7
JOB_STATUS_FILE_HOLDER = 8
JOB_STATUS_TO_DELETE = 9
JOB_STATUS_UNSATISFACTORY = 10

USER_ROLE_UNKNOWN = 0
USER_ROLE_MANAGER = 1
USER_ROLE_OWNER = 2
USER_ROLE_USER = 3
USER_ROLE_REMOVED = 4

JOB_STATUS_TEXT = ['Unknown','Pending','Queued','Running','Interrupted','Failed','Finished','Running remotely','File holder','To delete','Unsatisfactory']

FINISHED_JOB_STATUS = ['Finished','Interrupted','To delete','Unsatisfactory']
JOB_EVALUATION_TEXT = ['Unknown','Best','Good','Rejected']
USER_AGENT_TEXT = ['Unknown','CCP4i2','CCP4mg','Coot']
FILETYPES_TEXT = ['Unknown', 'application/CCP4-seq', 'chemical/x-pdb', 'MultiPDB', 'application/CCP4-mtz', 'application/CCP4-unmerged-mtz',
                  'application/CCP4-unmerged-experimental', 'application/CCP4-map', 'application/refmac-dictionary', 'application/refmac-TLS',
                  'application/CCP4-mtz-freerflag', 'application/CCP4-mtz-observed', 'application/CCP4-mtz-phases', 'application/CCP4-mtz-map', '',
                  'application/CCP4-seqalign', 'application/CCP4-mtz-mini', 'application/coot-script', 'application/refmac-external-restraints',
                  'application/CCP4-scene', 'application/CCP4-shelx-FA', 'application/phaser-sol', 'chemical/x-mdl-molfile',
                  'application/iMosflm-xml', 'application/CCP4-image', 'application/CCP4-generic-reflections', 'application/HHPred-alignments',
                  'application/Blast-alignments', 'chemical/x-pdb-ensemble', 'application/CCP4-asu-content', 'application/dials-jfile',
                  'application/dials-pfile', 'application/phaser-rfile', 'application/refmac-keywords',
                  'application/x-pdf', 'application/postscript', 'application/EBI-validation-xml', 'chemical/x-cif' ]
FILETYPES_CLASS = ['DataFile', 'SeqDataFile', 'PdbDataFile', '', 'MtzDataFile', 'MtzDataFile',
                   'UnmergedDataFile', 'MapDataFile', 'DictDataFile', 'TLSDataFile',
                   'FreeRDataFile', 'ObsDataFile', 'PhsDataFile', 'MapCoeffsDataFile', '',
                   'SeqAlignDataFile', 'MiniMtzDataFile', 'CootHistoryDataFile', 'RefmacRestraintsDataFile',
                   'SceneDataFile', 'ShelxFADataFile', 'PhaserSolDataFile', 'MDLMolDataFile',
                   'ImosflmXmlDataFile', 'ImageFile', 'GenericReflDataFile', 'HhpredDataFile',
                   'BlastDataFile', 'EnsemblePdbDataFile', 'AsuDataFile', 'DialsJsonFile',
                   'DialsPickleFile', 'PhaserRFileDataFile', 'RefmacKeywordFile', 'PDFDataFile', 'PostscriptDataFile', 'EBIValidationXMLDataFile', 'MmcifReflDataFile' ]
MINIMTZFILETYPES = [10, 11, 12, 13]
FILE_ROLE_OUT = 0
FILE_ROLE_IN = 1
FILE_ROLE_IMPORT = 2

PATH_FLAG_JOB_DIR = 1
PATH_FLAG_IMPORT_DIR = 2

"""
Do not ever be tempted to renumber or delete anything here. The position of a type in this list is critically important.
"""
FILETYPELIST = [ (0,'Unknown','File type unknown'),
                    (1,'application/CCP4-seq','Model sequence'),
                    (2,'chemical/x-pdb','Model coordinates'),
                    (3,'MultiPDB','Multiple model coordinates'),
                    (4,'application/CCP4-mtz','Merged experimental data'),
                    (5,'application/CCP4-mtz-unmerged','Unmerged experimental data'),
                    (6,'application/CCP4-unmerged-experimental','Unmerged experimental data any format'),
                    (7,'application/CCP4-map','Electron density map'),
                    (8,'application/refmac-dictionary','Refmac dictionary'),
                    (9,'application/refmac-TLS','Refmac TLS'),
                    (10,'application/CCP4-mtz-freerflag','FreeR flag'),
                    (11,'application/CCP4-mtz-observed','Observed intensities and structure factors'),
                    (12,'application/CCP4-mtz-phases','Phases'),
                    (13,'application/CCP4-mtz-map','Map coefficients'),
                    (14,'Dummy','Dummy'),
                    (15,'application/CCP4-seqalign', 'Sequence alignment' ),
                    (16,'application/CCP4-mtz-mini','Experimental data object'),
                    (17,'application/coot-script','Coot script'),
                    (18,'application/refmac-external-restraints','Refmac external restraints'),
                    (19,'application/CCP4-scene','CCP4mg scene file'),
                    (20,'application/CCP4-shelx-FA','Shelx FA'),
                    (21,'application/phaser-sol','Phaser solutions'),
                    (22,'chemical/x-mdl-molfile','MDL Molfile'),
                    (23,'application/iMosflm-xml','iMosflm data'),
                    (24,'application/CCP4-image','Image file'),
                    (25,'application/CCP4-generic-reflections','Merged reflection data'),
                    (26,'application/CCP4-HHPred-alignments','HHPred sequence search results'),
                    (27,'application/CCP4-Blast-alignments','Blast sequence search results'),
                    (28,'chemical/x-pdb-ensemble','Ensemble model coordinates'),
                    (29,'application/CCP4-asu-content','Asu content'),
                    (30,'application/dials-jfile', 'Dials json data file'),
                    (31,'application/dials-pfile', 'Dials pickle data file'),
                    (32,'application/phaser-rfile', 'Phaser rotation solutions'),
                    (33,'application/refmac-keywords', 'Refmac5 keyword file'),
                    (34,'application/x-pdf', 'PDF File'),
                    (35,'application/postscript', 'Postscript file'),
                    (36,'application/EBI-validation-xml', 'Validation XML'),
                    (37,'chemical/x-cif', 'mmCif reflection data'),
                ]

KEYTYPELIST =  [ (0,'Unknown','Key type unknown'),
                    (1,'RFactor','R Factor'),
                    (2,'RFree','Free R Factor'),
                    (3,'completeness','model completeness'),
                    (4, 'spaceGroup','space group'),
                    (5, 'highResLimit', 'high resolution limit'),
                    (6, 'rMeas', 'Rmeas data consistency'),
                    (7,'FOM','figure of merit of phases'),
                    (8,'CFOM','correlation FOM'),
                    (9,'Hand1Score','Hand 1 score'),
                    (10,'Hand2Score','Hand 2 score'),
                    (11,'CC','correlation coefficient between Fo and Fc'),
                    (12,'nAtoms','number of atoms in model'),
                    (13,'nResidues','number of residues in model'),
                    (14,'phaseError','phase error'),
                    (15,'weightedPhaseError','weighted phase error'),
                    (16,'reflectionCorrelation','reflection correlation'),
                    (17,'RMSxyz','RMS displacement'),
                    (18,'cutoff','Pairef cutoff'),
                    (19,'ccHalf','correlation coefficient between two half datasets'),
                    (20,'R','R value'),
                    (21,'R1Factor','R1 value'),
                    (22,'R1Free','Free R1 value'),
                    (23,'R1','R1 value'),
                    (24,'FSCaverage','Weighted average of map-model Fourier shell correlation'),
                    (25,'CCFwork_avg','Weighted average correlation between observed and calculated amplitudes within resolution bins'),
                    (26,'CCFfree_avg','Free weighted average correlation between observed and calculated amplitudes within resolution bins'),
                    (27,'CCF_avg','Weighted average correlation between observed and calculated amplitudes within resolution bins'),
                    (28,'CCIwork_avg','Weighted average correlation between observed and calculated intensities within resolution bins'),
                    (29,'CCIfree_avg','Free weighted average correlation between observed and calculated intensities within resolution bins'),
                    (30,'CCI_avg','Weighted average correlation between observed and calculated intensities within resolution bins'),]

FILEASSOCIATIONTYPELIST =  [(0,'Unknown','File association type unknown'),
                            (1,'Observed-Free','Observed data and FreeR set') ]

FILEASSOCIATIONROLELIST =  [ (0,'Unknown','File association role unknown',0),
                             (1,'Observed data','Observed data',1),
                             (2,'FreeR set','FreeR set',1) ]

class CDbApi(CObject):
    '''The main API class for the CCP4DB database'''
    projectReset = QtCore.Signal(dict)
    jobDeleted = QtCore.Signal(dict)
    jobCreated = QtCore.Signal(dict)
    jobUpdated = QtCore.Signal(dict)
    jobStatusUpdated = QtCore.Signal(dict)
    jobStarted = QtCore.Signal(dict)
    jobFinished = QtCore.Signal(dict)
    jobToDelete = QtCore.Signal(dict)
    importFileCreated = QtCore.Signal(dict)
    importFileUpdated = QtCore.Signal(dict)
    importFileDeleted = QtCore.Signal(dict)
    exportFileCreated = QtCore.Signal(dict)
    commentEdited = QtCore.Signal(dict)
    commentDeleted = QtCore.Signal(dict)
    setJobToImportSignal = QtCore.Signal(dict)
    fileUseCreated = QtCore.Signal(dict)
    fileCreated = QtCore.Signal(dict)
    fileUpdated = QtCore.Signal(dict)
    projectsListChanged = QtCore.Signal()
    followFromJobChanged = QtCore.Signal(tuple)
    projectUpdated = QtCore.Signal(dict)
    projectDeleted = QtCore.Signal(dict)
    projectTagsChanged = QtCore.Signal(str)
    tagCreated = QtCore.Signal(dict)
    tagDeleted = QtCore.Signal(dict)
    tagUpdated = QtCore.Signal(dict)
    unusedTagsDeleted = QtCore.Signal()
    projectCommentEdited = QtCore.Signal(dict)
    projectCommentDeleted = QtCore.Signal(dict)

    ERROR_CODES = {101 : { 'description' : 'Error connecting to database' },
                   102 : { 'description' : 'Error creating database cursor' },
                   103 : { 'description' : 'Database commit called when there is no open connection' },
                   104 : { 'description' : 'Error in database commit', 'severity' : SEVERITY_CRITICAL },
                   105 : { 'description' : 'Error reading SQL file' },
                   106 : { 'description' : 'Attempting operation for which do not have pivilege' },
                   107 : { 'description' : 'User name known but no longer has access privilege' },
                   108 : { 'description' : 'Unable to retreive permissions for project and user' },
                   109 : { 'description' : 'Error enabling foreign keys in sqlite' },
                   110 : { 'description' : 'There is already project with same name in database' },
                   111 : { 'description' : 'There is no project with this name in database' },
                   112 : { 'description' : 'User does not have permission to write to project' },
                   113 : { 'description' : 'Attempting to create project without giving directory' },
                   114 : { 'description' : 'Attempting to access non-existent jobId/fileId' },
                   115 : { 'description' : 'Attempting to access non-existent job data type' },
                   116 : { 'description' : 'No suitable job status provided' },
                   117 : { 'description' : 'Attempting to create project with same directory as existing project' },
                   118 : { 'description' : 'Attempting to access non-existent projectId' },
                   119 : { 'description' : 'Attempting to access non-existent project data type' },
                   120 : { 'description' : 'There is already user with same name in database' },
                   121 : { 'description' : 'There is no user with this name in database' },
                   122 : { 'description' : 'Attempting to set inappropriate user role' },
                   123 : { 'description' : 'Error creating unique next job number for project' },
                   124 : { 'description' : 'Error attempting to set inappropriate followFromJobId' },
                   125 : { 'description' : 'Attempting to retrieve user info with invalid user id' },
                   126 : { 'description' : 'Error retrieve user info with user id' },
                   127 : { 'description' : 'Attempting to create file record for file not in job directory' },
                   128 : { 'description' : 'Error attempting to test if file is in correct directory' },
                   131 : { 'description' : 'Unable to retrieve the project for jobID' },
                   132 : { 'description' : 'Can not queue a job if the controlFile is not set or does not exist' },
                   133 : { 'description' : 'Error creating/updating directory alias' },
                   134 : { 'description' : 'Attempting to access non-existant projectName/jobNumber' },
                   135 : { 'description' : 'Attempting to set invalid job control file' },
                   136 : { 'description' : 'User does not have permission to view project' },
                   137 : { 'description' : 'No job information recovered for this jobId' },
                   140 : { 'description' : 'Failed to import QtSql - is Qt installed?' },
                   141 : { 'description' : 'Failed to create Qt database connection' },
                   142 : { 'description' : 'Failed to connect to database file' },
                   143 : { 'description' : 'Created database connection not valid' },
                   150 : { 'description' : 'Error executing database query'  , 'severity' : SEVERITY_CRITICAL},
                   151 : { 'description' : 'Error commiting database query' , 'severity' : SEVERITY_CRITICAL},
                   152 : { 'description' : 'Error opening connection to database for command' },
                   153 : { 'description' : 'Error creating database transaction for command' },
                   154 : { 'description' : 'Error creating database using Qt interface - failed to start Qt application' },
                   160 : { 'description' : 'There is already fileType with same name in database' },
                   161 : { 'description' : 'Unable to enter filname in database - bad job id' },
                   162 : { 'description' : 'Unable to enter filname in database - file does not exist' },
                   163 : { 'description' : 'Attempting to use non-existant jobId' },
                   164 : { 'description' : 'There is no file type with this name in database' },
                   165 : { 'description' : 'There is already a file type with this name in database' },
                   166 : { 'description' : 'Invalid value for RoleID when retrieving job file information' },
                   167 : { 'description' : 'Invalid value for mode when retrieving job file information' },
                   168 : { 'description' : 'Invalid fileId in getFullPath' },
                   169 : { 'description' : 'Invalid key name in updateJob' },
                   170 : { 'description' : 'Job params file not found' },
                   171 : { 'description' : 'Error gleaning file info for job - error creating file record' },
                   172 : { 'description' : 'Error gleaning file info for job - error creating fileUse record' },
                   173 : { 'severity': SEVERITY_WARNING, 'description' : 'Possible job output file does not exist' },
                   174 : { 'description' : 'Calling getProjectInfo with no valid project' },
                   175 : { 'description' : 'Attempting to set preceedingJobId with inappropriate jobId' },
                   180 : { 'description' : 'Calling getFileInfo with invalid fileId' },
                   190 : { 'description' : 'Attempting to set invalid data in a job record' },
                   200 : { 'description' : 'Error executing database command' },
                   201 : { 'severity': SEVERITY_WARNING, 'description' : 'Invalid key name in updateProject' },
                   202 : { 'description' : 'Invalid value data type in updateProject' },
                   203 : { 'description' : 'In updateProject setting parent project is already child of this project' },
                   211 : { 'description' : 'Error attempting to convert project data to XML' },
                   212 : { 'description' : 'Error attempting to write project data to XML file' },
                   220 : { 'description' : 'Database has no database version information' },
                   221 : { 'description' : 'Database has no schema version information' },
                   222 : { 'description' : 'Database is based on out-of-date schema' },
                   230 : { 'description' : 'Invalid key name in updateFile' },
                   231 : { 'description' : 'Unable to retrieve the project for fileID' },
                   240 : { 'description' : 'Attempting to create a Comment with no/invalid jobId' },
                   241 : { 'description' : 'Attempting to update a Comment with no/invalid commentId' },
                   242 : { 'description' : 'Attempting to delete a Comment with no/invalid commentId' },
                   250 : { 'description' : 'Attempting to create file import record with invalid source fileId' },
                   251 : { 'description' : 'Attempting to create file import record without source information' },
                   252 : { 'description' : 'Attempting to update file import record with invalid key' },
                   253 : { 'description' : 'Attempting to create file export record with invalid source fileId' },
                   254 : { 'description' : 'Attempting to create file export record with invalid (non-existant?) filename' },
                   255 : { 'description' : 'Error attempting to write import file to xml' },
                   256 : { 'description' : 'Error attempting to write export file to xml' },
                   257 : { 'description' : 'Error opening project xml file' },
                   258 : { 'description' : 'Error attempting to get instances of file exports' },
                   270 : { 'description' : 'Calling getImportFileInfo with invalid mode' },
                   271 : { 'description' : 'Calling getImportFileInfo with invalid importFileId' },
                   272 : { 'description' : 'getImportFileInfo can not find imported file in the File table' },
                   280 : { 'description' : 'getDatabaseInfo did not find a single database record' },
                   290 : { 'description' : 'There is no key type with this name in database' },
                   291 : { 'description' : 'There is no file association type with this name in database' },
                   292 : { 'description' : 'Project-tag already exists in database' },
                   293 : { 'description' : 'Invalid tag update - bad key' },
                   294 : { 'description' : 'Atempting to create project tag that already exists' }}

    DATABASEITEMS = ['databaseid','creationtime','creationhostname','creationusername','schemaversion','schemadate']
    DATABASETYPES = [UUIDTYPE,float,str,str,str,str]
    TIMEFORMAT = '%H:%M %d-%b-%y'
    USERITEMS = ['userid','username','userroleid']
    USERITEMTYPES = [UUIDTYPE,str,int]
    PROJECTITEMS = ['projectid','projectname','projectcreated','userid','parentprojectid','lastjobnumber','projectdirectory','followfromjobid','lastcleanuptime',
                    'i1projectname','i1projectdirectory','lastaccess']
    PROJECTTYPES = [ UUIDTYPE,str,float,UUIDTYPE,UUIDTYPE,int,str,UUIDTYPE,float,str,str,float]
    # used in getTablesEtree..
    PROJECTITEMS0 = ['projectid','projectname','projectcreated','lastjobnumber','projectdirectory']
    PROJECTTYPES0 = [ UUIDTYPE,str,float,int,str]

    JOBITEMS = ['jobid','jobnumber','creationtime','finishtime','status','evaluation',
                  'useragent','jobtitle','projectid','taskname','taskversion',
                  'parentjobid','preceedingjobid','processid']
    JOBTYPES = [ UUIDTYPE,str,float,float,int,int,
                 int,str,UUIDTYPE,str,str,
                 UUIDTYPE,UUIDTYPE,int]

    FILEITEMS = [ 'fileid','filename','jobid','jobparamname','filetypeid','annotation','filesubtype','filecontent','pathflag']
    FILETYPES = [UUIDTYPE,str,UUIDTYPE,str,int,str,int,int,int]

    IMPORTFILEITEMS = ['importid','fileid','creationtime','lastmodifiedtime','checksum','sourcefilename','sourcefileid','exportfileid','annotation','importnumber','reference']
    IMPORTFILETYPES = [UUIDTYPE,UUIDTYPE,float,float,str,str,UUIDTYPE,UUIDTYPE,str,int,str]

    EXPORTFILEITEMS = ['exportid','creationtime','fileid','exportfilename','annotation']
    EXPORTFILETYPES = [UUIDTYPE,float,UUIDTYPE,str,str]

    FILEUSEITEMS = ['fileid','jobid','roleid','jobparamname']
    FILEUSETYPES = [ UUIDTYPE,UUIDTYPE,int,str]

    XDATAITEMS = ['xdataid','xdataclass','jobid','xdataxml']
    XDATATYPES = [ UUIDTYPE,str,UUIDTYPE,str]

    COMMENTITEMS = ['commentid','jobid','userid','timeofcomment','comment']
    COMMENTTYPES = [ UUIDTYPE,UUIDTYPE,UUIDTYPE,float,str]

    PROJECTCOMMENTITEMS = ['projectcommentid','projectid','userid','timeofcomment','comment']
    PROJECTCOMMENTTYPES = [ UUIDTYPE,UUIDTYPE,UUIDTYPE,float,str]

    JOBKEYVALUEITEMS = [ 'jobid','keytypeid','value']
    JOBKEYVALUETYPES = [ UUIDTYPE, int, float ]
    JOBKEYCHARVALUEITEMS = [ 'jobid','keytypeid','value']
    JOBKEYCHARVALUETYPES = [ UUIDTYPE, int, str ]

    FILEASSOCIATIONITEMS = ['fileassociationid', 'fileassociationtypeid','fileassociationname' ]
    FILEASSOCIATIONTYPES = [ UUIDTYPE, int, str ]

    FILEASSOCIATIONMEMBERITEMS = [ 'fileassociationid','fileid','fileassociationroleid']
    FILEASSOCIATIONMEMBERTYPES = [ UUIDTYPE, UUIDTYPE, int ]

    PROJECTIMPORTITEMS = ['projectimportid','projectid','projectimporttime','projectexportid','projectexportdatabaseid','projectexporttime','projectexportafter']
    PROJECTIMPORTTYPES = [ UUIDTYPE,UUIDTYPE,float,UUIDTYPE,UUIDTYPE,float,float ]

    TAGITEMS = [ 'tagid', 'parenttagid','text']
    TAGTYPES = [UUIDTYPE,UUIDTYPE,str]

    PROJECTTAGITEMS = [ 'tagid' ,'projectid' ]
    PROJECTTAGTYPES = [UUIDTYPE,UUIDTYPE]

    insts = None

    def __init__(self,parent=None,mode=None,fileName=None,hostName='localhost',userName=None,userPassword=None,createDb=False,**kw):
        '''initialise CCP4DbApi class'''
        #print 'CDbApi.__init__ fileName',fileName
        CObject.__init__(self,parent)
        # Save te first instance of this class
        #print 'CDbApi.__init__ insts',CDbApi.insts,'self',self
        if CDbApi.insts is None:
            CDbApi.insts = self
        if mode is None or fileName is None:
            from core import CCP4Config
            mode = CCP4Config.DBMODE()
            fileName = CCP4Config.DBFILE()
            if userName is None: userName =  CCP4Config.DBUSER()
        self._schemaVersion = self.getSchemaVersion()
        print('Current schema version:', self._schemaVersion)
        if fileName is  None:
            if mode in ['sqlite','qt_sqlite']:
                from core import CCP4Utils
                fileName = os.path.join(CCP4Utils.getDotDirectory(),'db','database.sqlite')
        if mode in ['qt_sqlite']:
            # Ensure there is a QApplication - otherwise database drivers don't load
            from core import CCP4Modules
            try:
                qApp = CCP4Modules.QTAPPLICATION()
            except:
                raise CException(self.__class__,154)
        if userName is None:
            from core import CCP4Utils
            userName = CCP4Utils.getUserId()
        self._projectList = []
        self._dbMode = mode
        self._fileName=fileName
        self._hostName=hostName
        self._userName=userName
        self._userPassword=userPassword
        self._userRole = USER_ROLE_UNKNOWN
        self._preferences = {'testFileExists':True, 'useTransaction':False, 'useImportFileDir':True }
        self._conn=None
        self._cur=None
        self._open = False
        # Save record of commands before a commit
        self._commands = []
        self._diagnostic =kw.get('diagnostic',False)
        # Cache regularly accessed data
        self._projectPermissions = {}
        self._projectDirectories = {}
        self._followFromJob = {}    
        self.openDb(fileName=fileName,createDb=createDb,diagnostic=kw.get('loadDiagnostic',True))
        #self.setupBackupJobToXml()
        #self.connect(self,QtCore.SIGNAL('followFromJobChanged'),self.bleep)

    def loadFileTypes(self,diagnostic=True):
        self.execute('SELECT MAX(FiletypeID) FROM Filetypes')
        fmaxList = self.fetchAll2Py(int)
        if fmaxList[0] is None:
            fmax = 0
        else:
            fmax = fmaxList[0] + 1
        if fmax<len(FILETYPELIST):
            if diagnostic: print('Loading new file type into database',FILETYPELIST[fmax:])
            for item in FILETYPELIST[fmax:]:
                self.execute(' INSERT INTO Filetypes (FiletypeID,FiletypeName,FiletypeDescription) VALUES (?,?,?)',item)

    def loadJobStatus(self,diagnostic=True):
        self.execute('SELECT MAX(StatusID) FROM JobStatus')
        fmaxList = self.fetchAll2Py(int)
        if fmaxList[0] is None:
            fmax = 0
        else:
            fmax = fmaxList[0]+1
        if fmax<=len(JOB_STATUS_TEXT):
            for idx in range(fmax,len(JOB_STATUS_TEXT)):
                if diagnostic:
                    print('Loading new job status into database',idx,JOB_STATUS_TEXT[idx])
                self.execute(' INSERT INTO JobStatus (StatusID,StatusText) VALUES (?,?)',(idx,JOB_STATUS_TEXT[idx]))

    def loadKeyTypes(self,diagnostic=True):
        self.execute('SELECT MAX(KeyTypeID) FROM KeyTypes')
        fmaxList = self.fetchAll2Py(int)
        if fmaxList[0] is None:
            fmax = 0
        else:
            fmax = fmaxList[0] + 1
        #print 'loadFileTypes',fmaxList,fmax
        if fmax<len(KEYTYPELIST):
            if diagnostic: print('Loading new key type into database',KEYTYPELIST[fmax:])
            for item in KEYTYPELIST[fmax:]:
                self.execute(' INSERT INTO KeyTypes (KeyTypeID,KeyTypeName,KeyTypeDescription) VALUES (?,?,?)',item)

    def loadFileAssociationTypes(self):
        self.execute('SELECT MAX(FileAssociationTypeID) FROM FileAssociationTypes')
        fmaxList = self.fetchAll2Py(int)
        if fmaxList[0] is None:
            fmax = 0
        else:
            fmax = fmaxList[0] + 1
        #print 'loadFileTypes',fmaxList,fmax
        if fmax<len(FILEASSOCIATIONTYPELIST):
            print('Loading new file association type into database',FILEASSOCIATIONTYPELIST[fmax:])
            for item in FILEASSOCIATIONTYPELIST[fmax:]:
                self.execute(' INSERT INTO FileAssociationTypes (FileAssociationTypeID,FileAssociationTypeName,FileAssociationTypeDescription) VALUES (?,?,?)', item)
        self.execute('SELECT MAX(FileAssociationRoleID) FROM FileAssociationRoles')
        fmaxList = self.fetchAll2Py(int)
        if fmaxList[0] is None:
            fmax = 0
        else:
            fmax = fmaxList[0] + 1
        #print 'loadFileTypes',fmaxList,fmax
        if fmax<len(FILEASSOCIATIONROLELIST):
            print('Loading new file association roles into database',FILEASSOCIATIONROLELIST[fmax:])
            for item in FILEASSOCIATIONROLELIST[fmax:]:
                self.execute(' INSERT INTO FileAssociationRoles (FileAssociationRoleID,FileAssociationRoleName,FileAssociationRoleDescription,FileAssociationTypeID) VALUES (?,?,?,?)',item)

    def setPreference(self,key=None,value=None):
        if key in self._preferences:
            self._preferences[key] = value
        else:
            pass

    def bleep(self,args):
        print('DbApi.bleep', args)

    def cacheTextLabels(self):
        self.execute("SELECT StatusID,StatusText FROM JobStatus ORDER BY StatusID")
        self.fetchAll2PyList(self, toType=[int, str])
        # ?? What about the JOB_STATUS_UNKNOWN etc.


    def iAmUserAgent(self):
        #from core.CCP4Config import VERSION
        #return 'CCP4i2 '+str(VERSION())
        return 1  # ccp4i2

    def loadSchema(self,fileName=None):
        if fileName is  None:
            from core import CCP4Utils
            fileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'dbapi','database_schema.sql')
        print('Loading database schema from:', fileName, 'Version:', self._schemaVersion)
        self.read(fileName=fileName)
        print('Finished loading database schema.')

    def getSchemaVersion(self, fileName=None):
        from core import CCP4Utils
        if fileName is  None:
            fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), 'dbapi', 'database_schema.sql')
        text = CCP4Utils.readFile(fileName)
        hit = re.search(r'\n-- VERSION(.*)\n',text)
        if hit is not None:
            version = hit.groups()[0].strip()
        else:
            version = None
        hit = re.search(r'\n-- DATE(.*)\n',text)
        if hit is not None:
            date = hit.groups()[0].strip()
        else:
            date = None
        #print 'CDbApi.getSchemaVersion',version,date
        return (version, date)

    def getDbSchemaVersion(self):
      try:
        self.execute('SELECT SchemaVersion,SchemaDate FROM Databases WHERE ThisIsMe > 0 ORDER BY CreationTime DESC')
        rv = self.fetchAll2PyList([str,str])
      except:
        raise CException(self.__class__,220)
      #print 'CDbApi.validateSchema',rv
      if len(rv) != 1:
        raise CException(self.__class__,221)
      return rv[0]

    def validateSchemaVersion(self):
      version,date = self.getDbSchemaVersion()
      if version != self._schemaVersion[0]:
        raise CException(self.__class__,222,'Db version: '+str(version)+' current schema version: '+str(self._schemaVersion[0]))

    def updateDbSchema(self):
      version,date = self.getDbSchemaVersion()
      print('updateDbSchema',version,date)
      if version == '0.1.10':
        version,date = self.update_0_1_10()
      if version == '0.1.11':
       version,date = self.update_0_1_11()
      if version == '0.1.12':
       version,date = self.update_0_1_12()
      if version == '0.1.13':
       version,date = self.update_0_1_13()
      if version == '0.1.14':
       version,date = self.update_0_1_14()
      if version == '0.1.15':
       version,date = self.update_0_1_15()
      if version == '0.1.16':
       version,date = self.update_0_1_16()
      if version == '0.1.17':
       version,date = self.update_0_1_17()
      if version == '0.1.18':
       version,date = self.update_0_1_18()
      if version == '0.1.19':
       version,date = self.update_0_1_19()
      if version == '0.1.20':
       version,date = self.update_0_1_20()
      if version == '0.1.21':
       version,date = self.update_0_1_21()

      self._schemaVersion = version,date

    def update_0_1_10(self):
      print('Updating database from 0.1.10 to 0.1.11 ..')
      self.execute('ALTER TABLE Files ADD JobParamName VARCHAR(32)')
      self.execute('ALTER TABLE FileUses ADD JobParamName VARCHAR(32)')
      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',('0.1.11','23-07-2013'))
      self.commit()
      print('DONE updating database from 0.1.10 to 0.1.11')
      return ('0.1.11','23-07-2013')

    def update_0_1_11(self):
      last = '0.1.11'
      next =  '0.1.12'
      date = '14-12-2013'
      print('Updating database from '+last+' to '+next)

      # The important bit..
      self.execute('ALTER TABLE ImportFiles ADD LastModifiedTime REAL')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_12(self):
      last = '0.1.12'
      next = '0.1.13'
      date = '06-02-2014'
      print('Updating database from '+last+' to '+next)

      # The important bit..
      self.execute('ALTER TABLE ImportFiles ADD Checksum VARCHAR(32)')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_13(self):
      last = '0.1.13'
      next = '0.1.14'
      date = '07-02-2014'
      print('Updating database from '+last+' to '+next)

      # The important bit..
      self.execute('ALTER TABLE ImportFiles ADD ImportNumber MEDIUMINT')
      self.execute('ALTER TABLE ImportFiles ADD Reference TEXT')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_14(self):
      last = '0.1.14'
      next = '0.1.15'
      date = '15-05-2014'
      print('Updating database from '+last+' to '+next)

      # Oh dear! Here trying to rollback the ADD Userid if the subsequent add constraint fails.
      # Rollback is not possible due to Python sqlite bug that automatically commits DDL transactions
      # http://bugs.python.org/issue10740
      self.execute('ALTER TABLE Jobs ADD UserId VARCHAR(32) REFERENCES Users (UserId)')
      self.execute('''CREATE TABLE KeyTypes (
                         KeyTypeID MEDIUMINT NOT NULL,
                         KeyTypeName VARCHAR(50) NOT NULL,
	                 KeyTypeDescription TEXT,
	                 PRIMARY KEY (KeyTypeID),
	                 UNIQUE (KeyTypeID),
	                 UNIQUE (KeyTypeName) )''')

      self.execute('''CREATE TABLE JobKeyValues (
 	                   JobID VARCHAR(32) NOT NULL,
                           KeyTypeID MEDIUMINT NOT NULL,
                           Value REAL,
                           PRIMARY KEY (JobID,KeyTypeID),
                           FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
                           FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) )''')

      self.execute('''CREATE TABLE FileAssociationTypes (
        FileAssociationTypeID MEDIUMINT NOT NULL,
        FileAssociationTypeName VARCHAR(50) NOT NULL,
	FileAssociationTypeDescription TEXT,
	PRIMARY KEY (FileAssociationTypeID),
	UNIQUE (FileAssociationTypeID),
	UNIQUE (FileAssociationTypeName) )''')

      self.execute('''CREATE TABLE FileAssociationRoles (
        FileAssociationRoleID MEDIUMINT NOT NULL,
	FileAssociationTypeID MEDIUMINT NOT NULL,
        FileAssociationRoleName VARCHAR(50) NOT NULL,
	FileAssociationRoleDescription TEXT,
	PRIMARY KEY (FileAssociationRoleID),
	FOREIGN KEY (FileAssociationTypeID) REFERENCES FileAssociationTypes (FileAssociationTypeID),
	UNIQUE (FileAssociationRoleID),
	UNIQUE (FileAssociationRoleName) )''')

      self.execute('''CREATE TABLE FileAssociations (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileAssociationTypeId MEDIUMINT NOT NULL,
       FileAssociationName VARCHAR(50),
       PRIMARY KEY (FileAssociationId),
       UNIQUE (FileAssociationId) )''')

      self.execute('''CREATE TABLE FileAssociationMembers (
       FileAssociationID VARCHAR(32) NOT NULL,
       FileID VARCHAR(32) NOT NULL,
       FileAssociationRoleID MEDIUMINT NOT NULL,
       PRIMARY KEY (FileID,FileAssociationID,FileAssociationRoleID),
       FOREIGN KEY (FileID) REFERENCES Files (FileID),
       FOREIGN KEY (FileAssociationID) REFERENCES FileAssociations (FileAssociationID) )''')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_15(self):
      last = '0.1.15'
      next = '0.1.16'
      date = '19-08-2014'
      print('Updating database from '+last+' to '+next)


      self.execute('''CREATE TABLE JobKeyCharValues (
 	                   JobID VARCHAR(32) NOT NULL,
                           KeyTypeID MEDIUMINT NOT NULL,
                           Value VARCHAR(255),
                           PRIMARY KEY (JobID,KeyTypeID),
                           FOREIGN KEY (JobId)  REFERENCES Jobs (JobId),
                           FOREIGN KEY (KeyTypeId)  REFERENCES KeyTypes (KeyTypeId) )''')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_16(self):
      last = '0.1.16'
      next = '0.1.17'
      date = '30-01-15'
      print('Updating database from '+last+' to '+next)

      self.execute('ALTER TABLE Files ADD PathFlag MEDIUMINT')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_17(self):
      last = '0.1.17'
      next = '0.1.18'
      date = '01-03-16'
      print('Updating database from '+last+' to '+next)

      self.execute('ALTER TABLE Projects ADD I1ProjectName VARCHAR(100)')
      self.execute('ALTER TABLE Projects ADD I1ProjectDirectory TEXT')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_18(self):
      last = '0.1.18'
      next = '0.1.19'
      date = '09-06-16'
      print('Updating database from '+last+' to '+next)

      self.execute('ALTER TABLE Projects ADD LastAccess REAL')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))

      self.commit()

      self.setAllProjectsLastAccess()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_19(self):
      last = '0.1.19'
      next = '0.1.20'
      date = '02-08-16'
      print('Updating database from '+last+' to '+next)

      self.execute('ALTER TABLE Jobs ADD ProcessId INT')

      self.execute('''CREATE TABLE ServerJobs
    (
    JobId VARCHAR(32) NOT NULL,
	ServerProcessId INT,
	Machine VARCHAR(255),
	Username  VARCHAR(100),
	Mechanism VARCHAR(32),
	RemotePath VARCHAR(255),
	CustomCodeFile  VARCHAR(255),
	Validate VARCHAR(32),
	KeyFilename  VARCHAR(255),
    PRIMARY KEY(JobId),
	UNIQUE (JobId),
    FOREIGN KEY (JobId) REFERENCES Jobs  (JobId)
        ) ;
    ''' )

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))

      self.commit()

      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_20(self):
      last = '0.1.20'
      next = '0.1.21'
      date = '18-08-16'
      print('Updating database from '+last+' to '+next)

      self.execute('ALTER TABLE ServerJobs ADD ServerGroup VARCHAR(32)')

      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))

      self.commit()

      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def update_0_1_21(self):
      last = '0.1.21'
      next = '0.1.22'
      date = '23-09-16'
      print('Updating database from '+last+' to '+next)

      self.execute('''
  CREATE TABLE Tags (
       TagID VARCHAR(32) NOT NULL,
       ParentTagID VARCHAR(32),
       Text VARCHAR(50) NOT NULL,
       PRIMARY KEY (TagID),
       FOREIGN KEY (ParentTagID) REFERENCES Tags (TagID) );''')

      self.execute('''
    CREATE TABLE ProjectTags (
       TagID VARCHAR(32) NOT NULL,
       ProjectID VARCHAR(32) NOT NULL,
       PRIMARY KEY ( ProjectID, TagID ),
       FOREIGN KEY (TagID) REFERENCES Tags (TagID) ,
       FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID) );''')

      self.execute('''
     CREATE TABLE ProjectComments
      (
	ProjectCommentID VARCHAR(32) NOT NULL,
	ProjectID VARCHAR(32) NOT NULL,
	UserID VARCHAR(255),
	TimeOfComment FLOAT NOT NULL,
	Comment TEXT,
	PRIMARY KEY (ProjectID, ProjectCommentID),
        FOREIGN KEY (ProjectID) REFERENCES Projects (ProjectID)
) ;''')
      self.execute('UPDATE Databases SET SchemaVersion=? , SchemaDate=? WHERE ThisIsMe>0',(next,date))
      self.commit()
      print('DONE Updating database from '+last+' to '+next)
      return (next,date)

    def validateUser(self, userName, userPassword):
        try:
            uid = self.getUserId(userName)
        except CException as e:
            # What to do if invalid user!
            self.close()
            raise e
        role = self.getUserRole(userName)
        if role in [USER_ROLE_UNKNOWN,USER_ROLE_REMOVED]:
            raise CException(self.__class__,107,userName)


    def close(self):
        #print 'cDBApi.close',self._dbMode
        if self._dbMode == 'sqlite':
            if self._cur is not None:
                self._cur.close()
            self._cur = None
        else:
            if self._open:
                self._conn.close()
                self._open = False
            if self._conn is not None:
                QtSql.QSqlDatabase.removeDatabase('CCP4DB')
            self._conn = None

    def connection(self):
      if self._dbMode == 'sqlite':
        if self._conn is None:
          try:
            # MN added check_same_thread=False....A SERIOUS FIXME IF WE DEVELOP WEB ACCESS TO PROJECT DATA
            self._conn = sqlite3.connect(self._fileName,5.0,1, check_same_thread = False)
          except Exception as e:
            raise CException(self.__class__,101,self._fileName+'\nError:'+str(e))
          try:
            self.execute('PRAGMA foreign_keys = ON')
          except Exception as e:
            raise CException(self.__class__,109,'Error:'+str(e))
          try:
            self.execute('PRAGMA synchronous=OFF')
#            self.execute('PRAGMA locking_mode = EXCLUSIVE')
            self.execute('PRAGMA page_size=4096')
            self.execute('PRAGMA temp_store=MEMORY')
            self.execute('PRAGMA journal_mode=WAL')
            self.execute('PRAGMA count_changes=OFF')
            self.execute('PRAGMA cache_size=80000')
            self.execute('PRAGMA default_cache_size=80000')
          except Exception as e:
            raise CException(self.__class__,109,'Error:'+str(e))

        return self._conn

      else:
        self._conn = QtSql.QSqlDatabase.database('CCP4DB')
        if not self._conn.isValid():
          try:
            self._conn = QtSql.QSqlDatabase.addDatabase('QSQLITE','CCP4DB')
          except Exception as e:
            raise CException(self.__class__,141,'Error:'+str(e))
          try:
            self._conn.setDatabaseName(self._fileName)
          except Exception as e:
            raise CException(self.__class__,142,self._fileName+'\nError:'+str(e))
          if not self._conn.isValid():
            raise CException(self.__class__,143,self._fileName)
          #self._preferences['useTransaction'] = self._conn.driver().hasFeature(QtSql.QSqlDriver.Transactions)
        return self._conn

    def openDb(self,fileName=None,createDb=False,diagnostic=True):

      if self._conn is not None:
        self._conn.close()
        self._conn = None
        self._fileName = None

      if fileName is None: return

      if os.path.exists(fileName):
        print('CCP4i2 opening database file',fileName)
      elif createDb:
        print('CCP4i2 creating database file',fileName)
      else:
        print('CCP4i2 can not open non-existant database',fileName)
        return
      self._fileName = fileName
      if not os.path.exists(fileName):
        self.loadSchema()
        self._userRole = USER_ROLE_OWNER
        # Need to side-step the usual checks
        uid = self.uniqueId(table='Users',identifier='UserID')
        dbid = self.uniqueId(table='Databases',identifier='DatabaseID')
        try:
          from core import CCP4Utils
          hostName = CCP4Utils.getHostName()
        except:
          hostName=None
        self.execute("INSERT INTO Users (UserId,UserName,UserPassword,UserRoleID) VALUES( ?,?,?,?)",
                       (uid,self._userName,self.encodePassword(self._userPassword),self._userRole) )
        self.execute("INSERT INTO Databases (DatabaseId,SchemaVersion,SchemaDate,CreationTime,CreationUserName,CreationHostName,ThisIsMe) VALUES(?,?,?,?,?,?,?)",
                     (dbid,self._schemaVersion[0],self._schemaVersion[1],self.currentTime(),self._userName,hostName,1))
        self.commit()
      else:
        self.updateDbSchema()
        self.validateSchemaVersion()
        self.validateUser(self._userName,self._userPassword)
        self._userRole = self.getUserRole(self._userName)

      self.loadFileTypes(diagnostic=diagnostic)
      self.loadKeyTypes(diagnostic=diagnostic)
      self.loadJobStatus(diagnostic=diagnostic)
      #Add a commit here, since any or all of the above might have caused execution of SQL commands
      self.commit()

    def commit(self):
      if self._diagnostic:
        import traceback
        print('Commit from',traceback.extract_stack(limit=2)[-2][2])

      if self._conn is None:
        raise CException(self.__class__,103)

      if self._dbMode == 'sqlite':
        try:
          self._conn.commit()
          self.close()
        except Exception as e:
          raise CException(self.__class__,104,str(self._commands)+'\nError:'+str(e))
        self._commands = []
        self._cur = None
      else:
        if self._preferences['useTransaction'] and len(self._commands)>0:
          rv = self._conn.commit()
          if not rv:
            raise CException(self.__class__,151,str(self._commands))
        self._commands = []
        self._conn.close()
        self._open = False

    def rollback(self):
      if self._dbMode == 'sqlite':
        try:
          self._conn.rollback()
        except:
          print('Error in CDbApi.rollback')


    def read(self,fileName):
      from core import CCP4Utils
      #print 'CDbApi.read',fileName
      try:
        text = CCP4Utils.readFile(fileName)
      except Exception as e:
        raise CException(self.__class__,105,fileName+'\nError:'+str(e))

      if self._dbMode == 'sqlite':
        self.connection().executescript(text)
        self.commit()
      else:
        comList = text.split(';')
        for com0 in  comList[0:-1]:
          #print com0
          com = re.sub('\n','',com0)
          self.execute(com)
        self.commit()

    def setDiagnostic(self,mode=False):
      if mode in (True,False):
        self._diagnostic = mode


    def backupDB(self,f2):
        if self._dbMode == 'sqlite':
            try:
                self.backupDBInternal(f2)
            except Exception as e:
                if str(e).count('database') and str(e).count('locked'):
                    print('****CDbApi ERROR - trying to backup locked database')
                    self.commit()
                    self.close()
                    try:
                        self.backupDBInternal(f2)
                    except:
                        print('****CDbApi ERROR on retry:',e)
                        print('Backing up to:',f2)
                        if self._diagnostic: traceback.print_stack()
                        raise e
                else:
                    print('CDbApi ERROR',e)
                    print('Backing up to:',f2)
                    if self._diagnostic: traceback.print_stack()
                    raise e
        else:
            print('CDbApi ERROR',e)
            print('Backing up to:',f2)
            print('with database mode',self._dbMode)

    def backupDBInternal(self,f2):
        import tempfile
        import shutil

        con = self.connection()
        sql = "".join([s+"\n" for s in con.iterdump()])

        tfile = tempfile.NamedTemporaryFile(delete=False)
        tfn = tfile.name
        tfile.close()

        f2trunc = open(tfn,"w")
        f2trunc.close()
        print("Writing to",tfn)

        conbak = sqlite3.connect(tfn)
        cur = conbak.cursor()

        #I am doing a simpler thing. Splitting based on semi-colon + newline is not safe. Strings may contain
        #this sequence. e.g. in comments. The simpler way may be (theoretically?) slower, but should be safer.
        try:
            cur.executescript(sql)
        except:
            print("Fail",com)
            conbak.close()
            raise
        """
        p = re.compile(";[\r\n]+")

        sql_commands = p.split(sql)

        for com in sql_commands:
            if com == "COMMIT":
                continue
            try:
                cur.execute(com)
            except:
                print "Fail",com
                conbak.close()
                raise
        """

        conbak.commit()

        shutil.copy(tfn,f2)

        print("Saved backup database to",f2)
        conbak.close()

    def execute(self,com,args=[]):

      if self._diagnostic:
        print('CDbApi.execute:',com,args)
      self._commands.append((com,args))

      if self._dbMode == 'sqlite':
        try:
          #if self._diagnostic: print 'connection().execute(com,args)'
          self._cur=self.connection().execute(com,args)
          #if self._diagnostic: print 'DONE connection().execute(com,args)'
        except Exception as e:
          if str(e).count('database') and str(e).count('locked'):
            # This probably wont work as it is a different process that has the lock
            print('****CDbApi ERROR - trying to handle locked database')
            self.commit()
            self.close()
            #self.openDb(fileName=self.fileName)
            try:
              self._cur=self.connection().execute(com,args)
              print('****CDbApi apparently OK on retry')
            except:
              print('****CDbApi ERROR on retry:',e)
              print('Executing:',com,args)
              if self._diagnostic: traceback.print_stack()
              raise e
          else:
            print('CDbApi ERROR',e)
            print('Executing:',com,args)
            if self._diagnostic: traceback.print_stack()
            raise e
      else:
        if not self._open:
          rv = self.connection().open()
          if not rv:
            raise CException(self.__class__,152,com+' '+str(args),stack=False)
          if self._preferences['useTransaction']:
            rv = self.connection().transaction()
            if not rv:
              raise CException(self.__class__,153,com+' '+str(args),stack=False)
          self._open = True

        if len(args) == 0:
          self._query = self.connection().exec_(com)
        else:
          self._query = QtSql.QSqlQuery(com,self._conn)
          for item in args:
            if item is None:
#FIXME PYQT - No idea what this is about
              #self._query.addBindValue(QtCore.QVariant(QtCore.QVariant.String))
              self._query.addBindValue(str) # MAybe this ?
            else:
              self._query.addBindValue(item)
          self._query.exec_()

        lastError = self._query.lastError()
        if lastError.type() != QtSql.QSqlError.NoError:
          # TODO better error reporting
          print('execute isValid',self._query.isValid())
          raise CException(self.__class__,150,'Code: '+str(lastError.type())+'  '+com+' '+str(args))

    def fetchAll(self):
      if self._dbMode == 'sqlite':
        if self._cur is None:
          return []
        else:
          sqlList = self._cur.fetchall()
          self.close()
        if self._diagnostic:
          print('CDbApi.fetchAll',sqlList)
        return sqlList
      else:
        pass

    def fetchAll2Py(self,toType=str):
      pyList = []
      if self._dbMode == 'sqlite':
        if self._cur is None:
          return []
        else:
          sqlList = self._cur.fetchall()
          self.close()
        for row in sqlList:
          if row[0] is None:
            pyList.append(None)
          else:
            pyList.append(toType(row[0]))
        #print('CDbApi.fetchAll2Py',pyList)
        return pyList
      else:
        while next(self._query):
          qVal =  self._query.value(0)
          if qVal.isValid():
            try:
              pyVal = toType(qVal)
            except:
               pyVal = None
          else:
            pyVal = None
          pyList.append(pyVal)
        #print('CDbApi.fetchAll2Py',pyList)
        return pyList

    def fetchAll2PyList(self,toType=[str]):
      pyList = []
      if self._dbMode == 'sqlite':
        if self._cur is None:
          return []
        else:
          sqlList = self._cur.fetchall()
          self.close()
        pyList = []
        for row in sqlList:
          pyRow = []
          for ii in range(len(row)):
            if row[ii] is None:
              pyRow.append(None)
            else:
              try:
                pyRow.append(toType[ii](row[ii]))
              except:
                print('ERROR in CCP4DbApi.fetchAll2PyList',row,toType,ii)
                import traceback
                traceback.print_stack(limit=5)
                pyRow.append(None)
          pyList.append(pyRow)
        return pyList
      else:
        while next(self._query):
          pyRow = []
          for ii in range(len(toType)):
            qVal = self._query.value(ii)
            if qVal.isValid():
              try:
                pyVal = toType[ii](qVal)
              except:
                pyVal = None
              pyRow.append(pyVal)
            else:
              pyRow.append(None)
          pyList.append(pyRow)
        return pyList

    def idExists(self,idValue=None,table=None,idName=None):
      if idName is None:
        idName = { 'ImportFiles' : 'ImportID',
                   'ExportFiles' : 'ExportID' }.get(table,None)
        if idName is None:
          if table[-1] == 's':
            idName = table[0:-1]+'ID'
          else:
            idName = table+'ID'
      self.execute('SELECT '+idName+' FROM '+table+' WHERE '+idName+'=?',(idValue,))
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)==1:
        return True
      else:
        return False


    def setupBackupJobToXml(self):
      import functools
      self.jobCreated.connect(functools.partial(self.backupJobToXml,'jobCreated'))
      self.jobUpdated.connect(functools.partial(self.backupJobToXml,'jobUpdated'))
      self.importFileCreated.connect(functools.partial(self.backupJobToXml,'importFileCreated'))
      self.importFileUpdated.connect(functools.partial(self.backupJobToXml,'importFileUpdated'))
      self.exportFileCreated.connect(functools.partial(self.backupJobToXml,'exportFileCreated'))
      self.commentEdited.connect(functools.partial(self.backupJobToXml,'commentEdited'))
      self.commentDeleted.connect(functools.partial(self.backupJobToXml,'commentDeleted'))

    @QtCore.Slot(str,dict)
    def backupJobToXml(self,signal,args):
      from dbapi import CCP4DbUtils
      #print 'CDbApi.backupJobToXml',signal,args
      if signal == 'jobUpdated':
        if not args['key'] in ['preceedingjob','title','evaluation']: return
      jobId = args['jobId']

      jobInfo = self.getJobInfo(jobId=jobId,mode=['jobnumber','projectname','taskname'])
      jobDirectory = self.jobDirectory(jobId=jobId)
      backup = CCP4DbUtils.CJobDbBackup(jobId=jobId,jobNumber=jobInfo['jobnumber'],taskName=jobInfo['taskname'],
                              jobDirectory=jobDirectory,projectName=jobInfo['projectname'])

      if signal == 'jobUpdated':
        backup.updateJob(key=args['key'],value=args['value'])
        backup.save()

      if signal in [ 'importFileCreated','importFileUpdated']:
        importInfo = self.getImportFileInfo(importId=args['importId'])
        backup.editImportFile(fileId=importInfo.get('fileid',None),importId=args['importId'],
                    sourceFileName=importInfo.get('sourcefilename',None),annotation=importInfo.get('annotation',None),
                            creationTime=importInfo.get('creationtime',None))
        backup.save()

      if signal in [ 'exportFileCreated']:
        exportInfo = self.getExportFileInfo(exportId=args['exportId'])
        #print 'backupJobToXml exportInfo',exportInfo
        backup.editExportFile(fileId=exportInfo.get('fileid',None),exportId=args['exportId'],
               exportFileName=exportInfo.get('exportfilename',None),creationTime=exportInfo.get('creationtime',None))
        backup.save()

      if signal in ['commentEdited']:
        commentInfo = self.getCommentInfo(commentId=args['commentId'])
        backup.editComment(commentId=args['commentId'],userName=commentInfo.get('username',None),
                 timeOfComment=commentInfo.get('timeofcomment',None),comment=commentInfo.get('comment',None))
        backup.save()

    def uniqueId(self,table,identifier):
      #print 'uniqueId',table,identifier
      import uuid
      return uuid.uuid1().hex

    def integerId(self,table,identifier):
      self.execute('SELECT LastID FROM LastUniqueIds WHERE TableName=?',(table,))
      rv = self.fetchAll2Py(int)
      if len(rv) == 1:
        newValue = rv[0] + 1
        self.execute('UPDATE LastUniqueIds SET LastID=? WHERE TableName=?',(newValue,table,))
        return newValue
      else:
        args = (identifier,table)
        self.execute('SELECT '+identifier+' FROM '+table+' ORDER BY '+identifier)
        l = self.fetchAll2Py(toType=int)
        if len(l)==0:
          newValue = 1
        else:
          newValue = l[-1]+1
        self.execute('INSERT INTO LastUniqueIds (TableName,LastID) VALUES (?,?)',(table,newValue))
        return newValue

    def createUser(self, userName, userPassword=None, userRole=USER_ROLE_USER):
        if not self._userRole in [USER_ROLE_MANAGER, USER_ROLE_OWNER]:
            raise CException(self.__class__, 106, 'You are ' + self._userName + ' calling createUser')
        try:
            self.getUserId(userName)
        except:
            pass
        else:
            # Already id for this userName - bad
            raise CException(self.__class__,120,userName)
        uniqueID = self.uniqueId(table='Users',identifier='UserID')
        if not userRole in [USER_ROLE_MANAGER,USER_ROLE_OWNER,USER_ROLE_USER]:
            raise CException(self.__class__,122,userName)
        self.execute("INSERT INTO Users (UserId,UserName,UserPassword,UserRoleID) VALUES( ?,?,?,?)",
                       (uniqueID,userName,self.encodePassword(userPassword),userRole) )
        self.commit()
        return uniqueID

    def encodePassword(self,userPassword=None):
      if userPassword is None or userPassword=='':
        return None
      else:
        '''
        import hashlib
        m = hashlib.md5()
        m.update(userPassword.encode())
        userPasswordhexdigest=m.hexdigest()
        #print ("%s" % userPasswordhexdigest)
        return userPasswordhexdigest
        '''
        return userPassword

    def getUserInfo(self,userId=None,mode='all',userName=None):
      if userId is None and userName is not None:
        userId = self.getUserId(userName)
      if userId is None:
        raise CException(self.__class__,124)

      # permission???

      if mode == 'all':
          itemList = []
          itemList.extend(self.USERITEMS)
      elif not isinstance(mode,list):
          itemList = [mode.lower()]
      else:
        itemList = []
        for item in mode: itemList.append(item.lower())

      selcom = "SELECT "
      typeConv = []
      for item in itemList:
        try:
          ii = self.USERITEMS.index(item)
          selcom = selcom + item +', '
          typeConv.append(self.USERITEMTYPES[ii])
        except:
          pass
      selcom = selcom[0:-2] +  " FROM Users WHERE UserID=?"
      self.execute(selcom,(userId,))
      rv = self.fetchAll2PyList(typeConv)
      #print 'getUserInfo rv',rv
      if len(rv) != 1:
        raise CException(self.__class__,125,str(userId))
      else:
        if len(itemList)==1:
          # Just return the value if only one item requested
          return rv[0][0]
        else:
          # Return multiple items as a dictionary
          ret = {}
          i = -1
          for itemName in itemList:
            i = i + 1
            ret[itemName] = rv[0][i]
          return ret


    def getUserId(self, userName):
        u = (userName,)
        self.execute('SELECT UserID FROM Users WHERE UserName = ?',u)
        rv = self.fetchAll2Py(UUIDTYPE)
        names = self.listUsers()
        if len(rv) != 1:
            # nb. Do not change einfo without checking utils/startup.py first.
            einfo =  {'c_user' : userName, 'p_users' : names}
            raise CException(self.__class__, 121, einfo, stack=False)
        else:
            return rv[0]

    def getUserRole(self, userName):
        u = (userName,)
        self.execute('SELECT UserRoleID FROM Users WHERE UserName = ?',u)
        rv = self.fetchAll2Py(int)
        if len(rv)==0:
            raise CException(self.__class__, 121, userName)
        else:
            return rv[0]

    def updateUser(self, userName, newUserName=None, userPassword=None, userRole=None):
        '''run an UPDATE query to update an info on a user in the Users table.'''
        if not self._userRole in (USER_ROLE_MANAGER,USER_ROLE_OWNER):
            raise CException(self.__class__,106)
            #'You are '+str(self._userName)+' calling updateUser')
        userId = self.getUserId(userName)
        if newUserName is not None:
            try:
                newId = self.getUserId(newUserName)
            except:
                pass
            else:
                raise CException(self.__class__, 120, newUserName)
            args = (newUserName,userId,)
            self.execute("UPDATE Users SET userName= ? WHERE userID= ?", args)
        if userPassword is not None:
            args = (self.encodePassword(userPassword),userId,)
            self.execute("UPDATE Users SET userPassword= ? WHERE userID= ?", args)
        if userRole is not None:
            # It must be one of allowed roles
            if not userRole in  [USER_ROLE_MANAGER,USER_ROLE_OWNER,USER_ROLE_USER,USER_ROLE_REMOVED]:
                raise CException(self.__class__, 122, userName)
            args =  (userRole,userId,)
            self.execute("UPDATE Users SET userRoleID= ? WHERE userID= ?",args)
        self.commit()
        self._userName = newUserName

    def listUsers(self, toTerm=False):
        self.execute("SELECT UserID, UserName, UserRoleID  FROM Users")
        userList = self.fetchAll2PyList([UUIDTYPE,str,int])
        if toTerm:
            print('listUsers')
            for item in userList:
                print (item)
        return userList

    def deleteUser(self,userName):
        '''Set userRole to USER_ROLE_REMOVED - safer than removing the record completely '''
        if not self._userRole in [USER_ROLE_MANAGER,USER_ROLE_OWNER]:
            raise CException(self.__class__,106,'You are '+self._userName+' calling deleteUser')
        userId = self.getUserId(userName)
        args =  (USER_ROLE_REMOVED,userId,)
        self.execute("UPDATE Users SET userRoleID= ? WHERE userID= ?",args)
        self.commit()

    def timeString(self,machineTime):
        return time.strftime(self.TIMEFORMAT,time.localtime(machineTime))

    def timeIsToday(self,machineTime):
        # Convert to 'year day' to
        return time.strftime('%y %j',time.localtime(machineTime))  == time.strftime('%y %j',time.localtime(time.time()))

    def matchProjectDirectory(self,directory=None,testAliases=False):
      '''Is there a  project with this directory already in the database'''
      from core import CCP4Utils
      if directory is None: return None,None,None
      self.execute('SELECT ProjectID, ProjectName, ProjectDirectory FROM Projects')
      dirList = self.fetchAll2PyList([UUIDTYPE,str,str])

      dirIn = os.path.abspath(directory)
      bestResult = [None,None,None,999]
      for  ProjectID, ProjectName, ProjectDirectory in dirList:
        try:
          relPath = os.path.relpath(dirIn,ProjectDirectory)
        except:
          pass
        else:
          if relPath == '.':
            return ProjectID,ProjectName,None
          elif relPath[0:2] != '..':
            nSubDir = CCP4Utils.splitPath(relPath)
            if nSubDir<bestResult[3]:
              bestResult = ProjectID,ProjectName,relPath, nSubDir

      '''
      if testAliases:
        #for name,path in self.directories.items():
        self.execute('SELECT DirectoryAlias, Directory FROM DirectoryAliases')
        dirList = self.fetchAll2PyList([str,str])
        for name,path in dirList:
          relPath = os.path.relpath(dirIn,str(path))
          if relPath == '.':
            return None,name
          elif relPath[0:2] != '..':
            nSubDir = CCP4Utils.splitPath(relPath)
            if nSubDir<bestResult[3]:
            bestResult = None,name,relPath,nSubDir
      '''

      return bestResult[0:3]


    def createProject(self,projectName,userName=None,parentProjectId=None,projectDirectory=None,projectId=None):
      '''Create a new project with the name "projectName" for user "userName".'''

      if not self._userRole in [USER_ROLE_MANAGER,USER_ROLE_OWNER]:
        raise CException(self.__class__,106,'You are '+self._userName+' calling createProject')

      if projectDirectory is None:
        raise CException(self.__class__,113,projectName)
      projectDirectory = os.path.normpath(projectDirectory)
      matchID,matchName,relpath = self.matchProjectDirectory(projectDirectory)
      if matchID is not None and (relpath is None or len(relpath)==0):
        raise CException(self.__class__,117,'Matching: '+str(matchName)+' Relative path: '+str(relpath),stack=False)

      try:
        pid = self.getProjectId(projectName)
      except:
        pass
      else:
        # Succeeded in getting projectId - baaad!
        raise CException(self.__class__,110,projectName)

      if userName is None:
        userId =self.getUserId(self._userName)
      else:
        userId = self.getUserId(userName)

      #print 'createProject',projectName,'userName',userName,userId
      if projectId is None:
        projectId = self.uniqueId(table='Projects',identifier='ProjectID')
      time = self.currentTime()

      try:
        args = (projectId,projectName,time,userId,parentProjectId,0,projectDirectory,time)
        self.execute("INSERT INTO Projects (ProjectID,ProjectName,ProjectCreated,UserID,ParentProjectID,LastJobNumber,ProjectDirectory,LastAccess) VALUES (?,?,?,?,?,?,?,?);",args)

        args = (projectId,userId,PRIVILEGE_EXTEND,userId)
        self.execute("INSERT INTO ProjectsUsersPermissions (ProjectID,UserID,PrivilegeID,CreatorID) VALUES (?,?,?,?);",args)
      except Exception as e:
        raise CException(self.__class__,200,'args: '+str(args)+' '+str(e))

      self.commit()
      self.projectsListChanged.emit()
      return projectId

    def getProjectInfo(self,projectId=None,mode='all',projectName=None,checkPermission=True):
      if projectId is None:
        if projectName is None: raise CException(self.__class__,174)
        projectId = self.getProjectId(projectName=projectName)


      if checkPermission:
        permission = self.getProjectPermission(projectId)
        if permission < PRIVILEGE_READ:
          raise CException(self.__class__,136)

      if mode == 'lastjobaccesstime':
        self.execute('SELECT MAX (FinishTime) FROM Jobs WHERE ProjectId=?',(projectId,))
        lastFinish = self.fetchAll2Py(float)
        self.execute('SELECT MAX (CreationTime) FROM Jobs WHERE ProjectId=?',(projectId,))
        lastCreation = self.fetchAll2Py(float)
        #print 'CDbApi.getProjectInfo lastjobaccesstime',lastFinish,lastCreation
        if len(lastFinish)>0 and lastFinish[0] is not None:
          if len(lastCreation)>0 and lastCreation[0] is not None:
            return max(lastFinish[0],lastCreation[0])
          else:
            return lastFinish[0]
        elif len(lastCreation)>0 and lastCreation[0] is not None:
          return lastCreation[0]
        else:
          return 0.0

      itemList = []
      extraItemList = []
      if mode == 'all':
          itemList.extend(self.PROJECTITEMS)
      elif not isinstance(mode,list):
          itemList = [mode.lower()]
      else:
        for item in mode:
          item0 =  item.lower()
          if item0 in self.PROJECTITEMS:
            itemList.append(item0)
          else:
            extraItemList.append(item0)

      ret = {}
      selcom = "SELECT "
      conv = []
      for item in itemList:
        selcom = selcom + item +', '
        conv.append(self.PROJECTTYPES[self.PROJECTITEMS.index(item)])
      selcom = selcom[0:-2] +  " FROM Projects WHERE ProjectID =?"

      #print 'getProjectInfo selcom',selcom
      self.execute(selcom,(projectId,))
      rv = self.fetchAll2PyList(conv)
      if len(rv) != 1:
         raise CException(self.__class__,118,str(projectId))
      else:
        ii = -1
        for item in itemList:
          ii = ii + 1
          ret[item] = rv[0][ii]

      if len(extraItemList)==0:
        if len(itemList) == 1:
          return ret[itemList[0]]
        else:
          return ret

      if 'parentprojectname' in extraItemList:
        self.execute('SELECT parentprojectid FROM Projects WHERE ProjectID =?',(projectId,))
        rv = self.fetchAll2Py(UUIDTYPE)
        if len(rv)!=1:
          raise CException(self.__class__,118,str(projectId))
        elif rv[0] is None:
          ret['parentprojectname'] = None
        else:
          ret['parentprojectname'] =  self.getProjectInfo(rv[0],'projectname')
      elif 'childprojects' in extraItemList :
        ret['childprojects'] = self._getChildProjects(projectId,False)
      elif 'decendantprojects' in extraItemList:
        ret['decendantprojects'] = self._getChildProjects(projectId,True)

      if len(itemList)==0 and len(extraItemList)==1:
        return ret[extraItemList[0]]
      else:
        return ret

    def getProjectJobSearchInfo(self,projectId,topOnly=True):
      com = 'SELECT DISTINCT TaskName FROM Jobs WHERE ProjectId = ?'
      if topOnly: com = com + ' AND ParentJobId IS NULL'
      self.execute(com,(projectId,))
      taskNameList = self.fetchAll2Py(str)
      self.execute('SELECT MIN(FinishTime) FROM Jobs WHERE Jobs.ProjectId = ? AND ParentJobId IS NULL',(projectId,))
      minTime = self.fetchAll2Py(int)
      self.execute('SELECT MAX(FinishTime) FROM Jobs WHERE Jobs.ProjectId = ? AND ParentJobId IS NULL',(projectId,))
      maxTime = self.fetchAll2Py(int)
      print('getProjectJobSearchInfo minTime',minTime,type(minTime))
      return { 'taskNameList' : taskNameList, 'minTime' : minTime[0], 'maxTime' : maxTime[0] }

    def getProjectJobsByFileType(self,projectId,fileTypeId=None):
      self.execute('''SELECT DISTINCT Jobs.JobNumber, Jobs.JobId, Jobs.TaskName, Jobs.CreationTime FROM Jobs
      INNER JOIN Files ON Files.JobId = Jobs.JobId
      WHERE Jobs.ProjectId = ? AND Jobs.ParentJobId IS NULL AND Files.FileTypeId = ?
      UNION
      SELECT DISTINCT Jobs.JobNumber, Jobs.JobId, Jobs.TaskName,Jobs.CreationTime  FROM Jobs
      INNER JOIN FileUses ON FileUses.JobId = Jobs.JobId
      INNER JOIN Files ON FileUses.FileId = Files.FileId
      WHERE Jobs.ProjectId = ? AND Jobs.ParentJobId IS NULL AND Files.FileTypeId = ? AND FileUses.RoleId = 1
      ORDER BY Jobs.CreationTime''',
      (projectId,fileTypeId,projectId,fileTypeId))
      jobList = self.fetchAll2PyList([str,str,str,float])

      return jobList

    def loadTaskTitleTable(self,inputList):
      '''Load task shorttitles to temporary table for searching'''
      self.execute('''CREATE TABLE IF NOT EXISTS  TempTaskTitles
( TaskName VARCHAR(100) NOT NULL,
TaskTitle TEXT );''')
      self.execute('DELETE FROM TempTaskTitles;')
      self.commit()

      for taskName,title in inputList:
        self.execute('INSERT INTO TempTaskTitles (TaskName,TaskTitle) VALUES (?,?)',(taskName,title))
      self.commit()


    def projectJobSearch(self,projectId,topOnly=True,searchParams={}):
      jobList = []
      args = [projectId]
      #print 'projectJobSearch searchParams',searchParams

      if searchParams.get('searchText',None) is None:
        com = 'SELECT Jobs.JobId,Jobs.JobNumber FROM Jobs WHERE Jobs.ProjectId = ?'
      else:
        com = '''SELECT Jobs.JobId,Jobs.JobNumber FROM Jobs
                 LEFT OUTER JOIN Comments ON Comments.JobId=Jobs.JobId
                 LEFT OUTER JOIN TempTaskTitles ON TempTaskTitles.TaskName = Jobs.TaskName
                 WHERE Jobs.ProjectId = ? '''

      if topOnly: com = com + ' AND Jobs.ParentJobId IS NULL'


      if searchParams.get('taskName',None) is not None:
        if not isinstance(searchParams['taskName'],list):
          com = com + ' AND Jobs.TaskName = ?'
          args.append(searchParams['taskName'])
        else:
          com = com + ' AND Jobs.TaskName IN (?'
          for item in searchParams['taskName'][1:]: com = com+',?'
          com = com + ')'
          args.extend(searchParams['taskName'])

      if searchParams.get('minTime',None) is not None:
        com = com + ' AND Jobs.FinishTime > ?'
        args.append(searchParams['minTime'])

      if searchParams.get('maxTime',None) is not None:
        com = com + ' AND Jobs.FinishTime < ?'
        args.append(searchParams['maxTime'])

      if searchParams.get('searchText',None) is not None:
        # textSearchMode values 0,1,2 == is, startwith, contains
        # may need wildcard adding
        if searchParams.get('textSearchMode',0)>0:
          if searchParams['searchText'][-1] != '%':
            searchParams['searchText'] = searchParams['searchText'] + '%'
          if searchParams['textSearchMode'] == 2 and searchParams['searchText'][0] != '%':
            searchParams['searchText'] = '%' + searchParams['searchText']
        com = com + ' AND ( Jobs.JobTitle LIKE ? OR Comments.Comment LIKE ? OR (Jobs.JobTitle IS NULL AND TempTaskTitles.TaskTitle LIKE ?))'
        args.append(searchParams['searchText'])
        args.append(searchParams['searchText'])
        args.append(searchParams['searchText'])

      #print 'projectJobSearch',com,args
      self.execute(com,args)
      jobList = self.fetchAll2PyList([str,str])
      #print 'projectJobSearch',jobList
      return jobList

    def projectSearch(self,searchParams={}):

      #print 'projectSearch searchParams',searchParams
      args = []
      joinTagsComments = """
      LEFT OUTER JOIN ProjectComments ON ProjectComments.ProjectId=Projects.ProjectId
      LEFT OUTER JOIN ProjectTags ON ProjectTags.ProjectId = Projects.ProjectId
      LEFT OUTER JOIN Tags ON Tags.TagId = ProjectTags.TagId
      """
      joinImportFiles = """
      LEFT OUTER JOIN Jobs ON Jobs.ProjectId = Projects.ProjectId
      LEFT OUTER JOIN Files ON Files.JobId = Jobs.JobId
      LEFT OUTER JOIN ImportFiles ON ImportFiles.FileId = Files.FileId
      """

      if searchParams.get('parent',None) is not None:
        com = """WITH RECURSIVE subProjects(n) AS (
              VALUES ( '"""+str(searchParams['parent'])+"""' )
              UNION
              SELECT Projects.ProjectId FROM Projects, subProjects
              WHERE Projects.ParentProjectId=subProjects.n
            )
            SELECT Projects.ProjectId,Projects.ProjectName FROM Projects"""
        if searchParams.get('searchText',None) is not None: com += joinTagsComments
        if searchParams.get('importedFile',None) is not None: com += joinImportFiles
        com +=   """\nWHERE Projects.ProjectId IN subProjects AND Projects.ProjectId != ?"""
        args.append(str(searchParams['parent']))
        and0 = ' AND'

      else:
        com = 'SELECT Projects.ProjectId,Projects.ProjectName FROM Projects'
        if searchParams.get('searchText',None) is not None: com += joinTagsComments
        if searchParams.get('importedFile',None) is not None: com +=joinImportFiles
        com += ' WHERE '
        and0 = ''


      if searchParams.get('useDate',False) and searchParams.get('minTime',None) is not None:
        if searchParams['dateMode'] in ['active', 'lastAccess']:
          com = com + and0 + ' Projects.LastAccess > ?'
        elif searchParams['dateMode'] == 'created':
          com = com + and0 + ' Projects.ProjectCreated > ?'
        and0 = ' AND'
        args.append(searchParams['minTime'])

      if searchParams.get('useDate',False) and searchParams.get('maxTime',None) is not None:
        if searchParams['dateMode'] == 'lastAccess':
          com = com + and0 + ' Projects.LastAccess < ?'
        elif searchParams['dateMode'] in ['active','created']:
          com = com + and0 + ' Projects.ProjectCreated < ?'
        and0 = ' AND'
        args.append(searchParams['maxTime'])

      if searchParams.get('importedFile',None) is not None:
         com = com + and0 + ' ImportFiles.SourceFileName = ?'
         args.append(searchParams['importedFile'])
         and0 = ' AND'

      if searchParams.get('projectDir',None) is not None:
         com = com + and0 + ' Projects.ProjectDirectory = ?'
         args.append(searchParams['projectDir'])
         and0 = ' AND'

      if searchParams.get('searchText',None) is not None:
        # textSearchMode = 0,1,2 == is, startwith, contains
        # may need to add wildcard at start/end
        if searchParams.get('textSearchMode',0)>0:
          if searchParams['searchText'][-1] != '%':
            searchParams['searchText'] = searchParams['searchText'] + '%'
          if searchParams['textSearchMode'] == 2 and searchParams['searchText'][0] != '%':
              searchParams['searchText'] = '%' + searchParams['searchText']
        com1 = ''
        and1 = ''
        if searchParams.get('name',True):
          com1 = 'Projects.ProjectName LIKE ?'
          and1 = ' OR'
          args.append(searchParams['searchText'])
        if searchParams.get('annotation',True):
          com1 = com1 + and1 + ' ProjectComments.Comment LIKE ?'
          and1 = ' OR'
          args.append(searchParams['searchText'])
        if searchParams.get('tags',True):
          com1 = com1 + and1 + ' Tags.Text LIKE ?'
          args.append(searchParams['searchText'])
        com = com + and0 +' ('+com1+')'


      #print 'projectSearch',com,args
      self.execute(com,args)
      retList = self.fetchAll2PyList([str,str])
      #print 'projectSearch',retList
      return retList


    def getProjectSearchInfo(self):
      self.execute('SELECT DISTINCT p1.ProjectId, p1.ProjectName FROM Projects p1 LEFT OUTER JOIN Projects p2 ON p1.ProjectId = p2.parentProjectId WHERE p2.parentProjectId IS NOT NULL ORDER BY p1.ProjectName')
      info = { 'parentProjects' : self.fetchAll2PyList([str,str]) }

      self.execute('SELECT MIN (ProjectCreated), MAX(ProjectCreated), MIN(LastAccess),MAX(LastAccess) FROM Projects')
      ret = self.fetchAll2PyList([float,float,float,float])
      info['projectCreated'] = ret[0][0:2]
      info['lastAccessed'] = ret[0][2:]
      print('getProjectSearchInfo',info)
      self.execute('SELECT DISTINCT TaskName FROM Jobs ORDER BY TaskName')
      info['taskName'] = self.fetchAll2Py(str)
      return info


    def setAllProjectsLastAccess(self):
      ''' Set the Projects.LastAccess for existing projects when updating to schema 0.1.19'''
      self.execute('SELECT ProjectId FROM Projects WHERE LastAccess IS NULL')
      projectIdList = self.fetchAll2Py(str)
      for pid in projectIdList:
        self.execute('UPDATE Projects SET LastAccess = ? WHERE ProjectId = ?',(self.getProjectInfo(pid,'lastjobaccesstime'),pid))
      self.commit()

    def updateProjectLastAcess(self,closingProjectId=None,openProjectIdList=[]):
      t = self.currentTime()
      self.execute('UPDATE Projects SET LastAccess = ? WHERE ProjectId = ?',(t,closingProjectId))
      for pid in openProjectIdList:
        self.execute('UPDATE Projects SET LastAccess = ? WHERE ProjectId = ?',(t+1.0,pid))
      self.commit()


    def resetLastJobNumber(self,projectId,commit=True):
      self.execute('SELECT JobNumber FROM Jobs WHERE ProjectId = ? AND ParentJobId IS NULL',(projectId,))
      # Misusing fetchAll2Py to return the job numbers as ints
      # MN: THis fails if we have a doubly nested job (117.8.1)
      #PRopose we retrieve as string, and then see if we need to split
      jobNumListAsStr = self.fetchAll2Py(str)
      jobNumList = [int(jobNum.split('.')[0]) for jobNum in jobNumListAsStr]
      if len(jobNumList) == 0:
        maxJobNum = 0
      else:
        maxJobNum = max(jobNumList)
      #print 'resetLastJobNumber',jobNumList,maxJobNum
      self.execute('UPDATE Projects SET LastJobNumber = ? WHERE ProjectId = ?',(str(maxJobNum),projectId))
      if commit: self.commit()

    def _getChildProjects(self,projectId,recurse=False):
      self.execute('SELECT ProjectID FROM Projects WHERE ParentProjectID =?',(projectId,))
      rv = self.fetchAll2Py(UUIDTYPE)
      if not recurse:
        return rv
      else:
        retval = []
        for pid in rv:
          retval.append(pid)
          children = self._getChildProjects(pid,True)
          retval.extend(children)
        return retval



    def getProjectId(self,projectName=None):
      #print 'CDbApi.getProjectId',projectName,type(projectName)
      if projectName is None:
        raise CException(self.__class__,111,str(projectName))
      u = (str(projectName),)
      self.execute('SELECT ProjectID FROM Projects WHERE ProjectName = ?',u)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)==0:
        raise CException(self.__class__,111,str(projectName))
      else:
        return rv[0]


    def getProjectUnfinishedJobs(self,projectId=None):
      self.execute('SELECT JobId FROM JOBS WHERE ProjectId=? AND Status != ?',(str(projectId),JOB_STATUS_FINISHED))
      rv = self.fetchAll2Py(UUIDTYPE)
      return rv

    def getProjectPermission(self,projectId=None,userName=None,projectName=None,jobId=None):

      return 4

      if projectName is None and projectId is None and jobId is not None:
        projectId = self.getJobInfo(jobId,'projectid')

      # Is it cached?
      if projectName is not None:
        permission = self._projectPermissions.get(userName,{}).get(projectName,None)
      else:
        permission = self._projectPermissions.get(userName,{}).get(projectId,None)
      if permission is not None: return permission

      # Retrieve from database
      if projectId is None:
        if projectName is not None:
          projectId = self.getProjectId(projectName)
      if userName is None: userName = self._userName
      userId = self.getUserId(userName)
      #print 'getProjectPermission userName,userId',userName,userId,type(userId),'projectId',projectId,type(projectId)
      arg = (UUIDTYPE(projectId),userId)
      self.execute('SELECT PrivilegeID FROM ProjectsUsersPermissions  WHERE ProjectID = ? AND UserID = ?',arg)
      rv = self.fetchAll2Py(int)
      if len(rv)==0:
        raise CException(self.__class__,108,'For '+str(projectId)+':'+str(userName))
      else:
        if userName not in self._projectPermissions: self._projectPermissions[userName] = {}
        if projectName is not None:
          self._projectPermissions[userName][projectName] = rv[0]
        else:
          self._projectPermissions[userName][projectId] = rv[0]
        return rv[0]

    def getProjectJobFileName(self,projectId=None,fileName=None,jobNumber="1",subJobNumber=""):
        if subJobNumber != "":
            fname = os.path.join(self.getProjectDirectory(projectId=projectId),"CCP4_JOBS","job_"+jobNumber,"job_"+subJobNumber,fileName)
        else:
            if "." in jobNumber:
                subJobNumber = jobNumber.split(".")[1]
                jobNumber = jobNumber.split(".")[0]
                fname = os.path.join(self.getProjectDirectory(projectId=projectId),"CCP4_JOBS","job_"+jobNumber,"job_"+subJobNumber,fileName)
            else:
                fname = os.path.join(self.getProjectDirectory(projectId=projectId),"CCP4_JOBS","job_"+jobNumber,fileName)
        return fname

    def getProjectJobFile(self,projectId=None,fileName=None,jobNumber="1",subJobNumber=""):
        fname = self.getProjectJobFileName(projectId,fileName,jobNumber,subJobNumber)
        d = fname
        if fname.endswith(".png"):
            with open(fname,"rb") as f:
               d = f.read()
        else:
            with open(fname) as f:
               d = f.read()
        return d

    def getProjectDirectory(self,projectId=None,projectName=None,jobId=None):

      #print 'getProjectDirectory input',projectId,projectName,jobId
      #print 'getProjectDirectory cache',self._projectDirectories
      if projectName is None and projectId is None and jobId is not None:
        projectId = self.getJobInfo(jobId,'projectid')

      # Is it cached?
      diry = None
      if projectId is not None:
        diry = self._projectDirectories.get(projectId,None)
      elif projectName is not None:
        diry = self._projectDirectories.get(projectName,None)
      if diry is not None: return diry

      # Retrieve from database
      if projectId is None:
        if projectName is not None: projectId = self.getProjectId(projectName)
        if projectId is None: raise CException(self.__class__,118)

      #print 'getProjectDirectory',projectId,type(projectId)
      arg = (projectId,)
      self.execute('SELECT ProjectDirectory FROM Projects  WHERE ProjectID = ?',arg)
      rv = self.fetchAll2Py(str)
      #print 'getProjectDirectory execute rv',rv
      if len(rv)==0:
        raise CException(self.__class__,118,'For '+str(projectId))
      else:
        diry = os.path.normpath(rv[0])
        if projectId is not None:
          self._projectDirectories[projectId] = diry
        else:
          self._projectDirectories[projectName] = diry
        return   diry

    def getProjectDirectoryList(self,order='ASC'):
      if not order in ['ASC','DESC']: order = 'ASC'
      self.execute('SELECT ProjectID,ProjectName,ProjectDirectory,ParentProjectID,ProjectCreated,LastAccess FROM Projects ORDER BY ProjectName '+order)
      return self.fetchAll2PyList([UUIDTYPE,str,str,UUIDTYPE,float,float])

    def getRecentProjects(self,order='name',limit=10):
      # LastAccess is last time that project window was open
      # Get most recent for projects menu
      self.execute('SELECT ProjectID,ProjectName,LastAccess FROM Projects ORDER BY LastAccess DESC LIMIT ?',(limit,))
      retList = self.fetchAll2PyList([UUIDTYPE,str,float])
      return retList

    def getProjectsToCleanup(self,period=86400.0):
      # get projects accessed after LastCleanupTime (- period which is probably 86400 (1 day) to allow for cleanup leaving recent files)
      self.execute('SELECT ProjectID,ProjectDirectory,ProjectName FROM Projects WHERE LastCleanupTime IS NULL OR LastAccess > LastCleanupTime - ?',(period,))
      retList = self.fetchAll2PyList([UUIDTYPE,str,str])
      return retList

    def lastProjectCleanup(self):
      self.execute('SELECT MAX (LastCleanupTime) FROM Projects')
      return self.db.fetchAll2Py(float)[0]

    def getProjectFollowFromJobId(self,projectId=None):
      if projectId in self._followFromJob:
        return self._followFromJob[projectId]
      arg = (projectId,)
      self.execute('SELECT followFromJobId FROM Projects  WHERE ProjectID = ?',arg)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)==0:
        return None
      else:
        self._followFromJob[projectId] = rv[0]
        return rv[0]

    def setProjectFollowFromJobId(self,projectId=None,jobId=None,clear=False):
      permission = self.getProjectPermission(projectId=projectId)
      if permission <  PRIVILEGE_WRITE:
        raise CException(self.__class__,112,'setProjectFollowFromJobId to project: '+str(projectId))

      # Check the job is in this project
      u = (str(jobId),)
      self.execute('SELECT ProjectID FROM Jobs WHERE JobID = ?',u)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)==0:
        raise CException(self.__class__,124,str(jobId))
      elif rv[0] != projectId:
        raise CException(self.__class__,124,str(jobId))

      # Get the old followFromJobId so can add to emited signal
      previousFollowFromJobId = self.getProjectFollowFromJobId(projectId=projectId)

      if clear:
        if jobId == previousFollowFromJobId:
          args = (None,projectId)
          self.execute('UPDATE Projects SET followFromJobId= ? WHERE ProjectID= ?',args)
          jobId = None
        else:
          return
      else:
        args = (jobId,projectId)
        self.execute('UPDATE Projects SET followFromJobId= ? WHERE ProjectID= ?',args)
      self.commit()
      self._followFromJob[projectId] = jobId
      #print 'setProjectFollowFromJobId',projectId,previousFollowFromJobId,jobId
      self.followFromJobChanged.emit((projectId,previousFollowFromJobId,jobId))
      return

    def considerUpdatingFollowFrom(self,projectId,finishedJobId):
      jobInfo = self.getJobInfo(jobId=finishedJobId,mode=['status','parentjobid','preceedingjobid'])
      if jobInfo['status'] != 'Finished' or jobInfo['parentjobid'] is not None: return
      followJobId = self.getProjectFollowFromJobId(projectId=projectId)
      #print 'CDbApi.considerUpdatingFollowFrom followJobId',finishedJobId,jobInfo,followJobId
      if followJobId is None:
        # Theres no follow from job so set it if this job has output mtz/pdb
        fileTypeList = self.getJobFiles(jobId=finishedJobId,mode='fileTypeId')
        #print 'considerUpdatingFollowFrom fileInfo',fileTypeList
        nF = 0
        for fType in [2,3,4,5,6,11,12,13]: nF = nF + fileTypeList.count(fType)
        if nF>0:
          self.setProjectFollowFromJobId(projectId=projectId,jobId=finishedJobId)
      else:
        #print 'considerUpdatingFollowFrom preceedingjobid followJobId',jobInfo['preceedingjobid'],followJobId,type(jobInfo['preceedingjobid']),type(followJobId)
        if jobInfo['preceedingjobid'] == followJobId:
          self.setProjectFollowFromJobId(projectId=projectId,jobId=finishedJobId)

    def getJobPermission(self,jobId,userName=None):
      '''Return the permission this user has for this job'''

      # Find the project for the JobID
      u = (str(jobId),)
      #self.setDiagnostic(True)
      self.execute('SELECT ProjectID FROM Jobs WHERE JobID = ?',u)
      rv = self.fetchAll2Py(UUIDTYPE)
      #self.setDiagnostic(False)
      if len(rv)==0:
        raise CException(self.__class__,131,str(jobId))
      else:
        projectId = rv[0]

      if userName is None: userName = self._userName

      return projectId,self.getProjectPermission(projectId=projectId,userName=userName)

    def getFilePermission(self,fileId,userName=None):
      '''Return the permission this user has for this job'''

      # Find the project for the FileID
      u = (str(fileId),)
      self.execute('SELECT Jobs.JobID,Jobs.ProjectID FROM Jobs INNER JOIN Files ON Jobs.JobId = Files.JobId WHERE Files.FileID = ?',u)
      rv = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE])
      if len(rv)==0:
        raise CException(self.__class__,231,str(fileId))

      if userName is None: userName = self._userName

      return rv[0][0],rv[0][1],self.getProjectPermission(projectId=rv[0][1],userName=userName)

    def setProjectPrivileges(self,userNameList=[],projectNameList=[],privilege=PRIVILEGE_READ):
      '''Enable another user access to projects'''
      pass

    def closeProject(self,projectName):
      '''Prevent any further entries - remove all privileges'''
      pass

    def copyProject(self,projectName):
      '''Copy project to in-memory db and return'''
      pass

    def archiveProject(self,projectName):
      '''Make XML copy? and delete'''
      pass

    def exportJobSelectionCommand(self,projectId=None,after=None,status=None,jobList=None):
      if jobList is not None:
        # Expand jobList to include child jobs
        childJobs = []
        com = 'SELECT jobnumber FROM jobs WHERE jobid in ('
        for job in jobList: com += '?,'
        com = com[0:-1] + ')'
        self.execute(com,jobList)
        jobNumberList = self.fetchAll2Py(str)
        for jobNo in jobNumberList:
          self.execute("SELECT jobid FROM jobs where projectId = ? AND jobnumber LIKE '"+jobNo+".%'",(projectId,))
          childJobs.extend(self.fetchAll2Py(UUIDTYPE))
        #print 'exportJobSelectionCommand',jobNumberList,childJobs
        jobList.extend(childJobs)

        # Create a standard selection for the list of jobs
        jobSelectCom = ' WHERE Jobs.JobId IN ('
        jobSelectArgs = jobList
        for job in jobList: jobSelectCom +='?,'
        jobSelectCom = jobSelectCom[0:-1] +')'
      else:
        # Create a standard selection based on the alternative criteria
        jobSelectArgs = [projectId]
        jobSelectCom = ' WHERE Jobs.ProjectId = ?'
        if after is not None:
          jobSelectArgs.append(after)
          jobSelectCom +=' AND  Jobs.FinishTime >= ?'
        if status is not None:
          if status == 'exportable':
            status = [JOB_STATUS_PENDING,JOB_STATUS_INTERRUPTED,JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_UNSATISFACTORY,JOB_STATUS_FILE_HOLDER]
          elif not isinstance(status,list):
            status = [status]
          jobSelectCom +=' AND Jobs.Status IN ('
          for item in status:
            jobSelectArgs.append(item)
            jobSelectCom +='?,'
          jobSelectCom = jobSelectCom[0:-1] +')'
      jobSelectCom = jobSelectCom + ' ORDER BY Jobs.CreationTime'
      return jobSelectCom,jobSelectArgs

    def getJobsByTime(self,projectId=None,before=None,after=None,mode='JobID',topLevel=False,status=None,jobList=None):
      jobSelectCom,jobSelectArgs = self.exportJobSelectionCommand(projectId=projectId,after=after,status=status,jobList=jobList)
      com = 'SELECT '+mode+' FROM Jobs ' +  re.sub('Jobs\.','',jobSelectCom)

      if topLevel: com = com + ' AND parentJobId IS NULL'
      #self.setDiagnostic(True)
      self.execute(com,jobSelectArgs)
      #self.setDiagnostic(False)
      return self.fetchAll2Py(UUIDTYPE)


    def getUnfinishedExportJobs(self,projectId=None,after=None,status=None,jobList=None):
      status = [JOB_STATUS_PENDING,JOB_STATUS_INTERRUPTED]
      jobSelectCom,jobSelectArgs = self.exportJobSelectionCommand(projectId=projectId,after=after,status=status,jobList=jobList)
      com = 'SELECT JobId FROM Jobs ' + re.sub('Jobs\.','',jobSelectCom)
      self.execute(com,jobSelectArgs)
      jobList = self.fetchAll2Py(UUIDTYPE)
      return jobList

    def getFilesUsedInJobList(self,jobList=None):
      # Utility to aid exporting a limited number of job and ensuring that we get the input files
      # required by the job(s) that are not in the exported job directory(s)
      # Get the input files to a list of jobs that are not output by a job in that list
      # and also include all imported files
      if jobList is None or len(jobList)==0:
        return []

      com0 = '(?'
      for item in jobList[1:]: com0 = com0 + ',?'
      com0 = com0 + ')'
      com1 = 'SELECT DISTINCT Files.FileId,Files.Filename,Jobs.JobId,Jobs.JobNumber,ImportFiles.ImportId,ImportFiles.SourceFilename FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId INNER JOIN FileUses ON Files.FileId = FileUses.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId'
      com = com1 + ' WHERE FileUses.JobId IN '+com0+' AND Files.JobId NOT IN '+com0
      #print 'CDbApi.getFilesUsedInJobList',com
      args = []
      args.extend(jobList)
      args.extend(jobList)
      self.execute(com,args)
      ret = self.fetchAll2PyList([UUIDTYPE,str,UUIDTYPE,str,UUIDTYPE,str])
      #print 'getFilesUsedInJobList ret',ret

      com = com1 + ' WHERE ImportFiles.ImportId IS NOT NULL AND Jobs.JobId IN '+com0
      args = jobList
      self.execute(com,args)
      reti = self.fetchAll2PyList([UUIDTYPE,str,UUIDTYPE,str,UUIDTYPE,str])
      #print 'CDbApi.getFilesUsedInJobList reti',reti
      ret.extend(reti)

      return ret

    def getTablesEtree(self,projectId=None,after=None,status=None,jobList=None,inputFileList=None,inputFileFromJobList=None):
      from lxml import etree

      #print 'getTablesEtree',jobList,inputFileList,inputFileFromJobList

      errReport = CErrorReport()
      jobNumberList = []
      root = etree.Element('project')

      table = self.getDatabaseInfo(returnList=True)
      databaseTableEle=self.tableEtree('database',self.DATABASEITEMS,table)
      root.append(databaseTableEle)

      com = 'SELECT '+self.PROJECTITEMS0[0]
      for item in self.PROJECTITEMS0[1:]: com = com + ', '+item
      com = com + ' FROM Projects WHERE ProjectID =?'

      self.execute(com,(projectId,))
      table = self.fetchAll2PyList(self.PROJECTTYPES0)
      projectTableEle = self.tableEtree('project',self.PROJECTITEMS0,table)
      root.append(projectTableEle)

      info = self.getProjectInfo(projectId=projectId,mode=['userid','parentprojectname'])
      self.execute('SELECT UserName FROM Users WHERE UserId=?',(info['userid'],))
      names = self.fetchAll2Py(str)
#FIXME - The above may return empty list. Hmm. Seems ProjectComments are not imported properly. Fix here.
      if len(names) == 0:
         self.execute('SELECT UserName FROM Users')
         allNames = self.fetchAll2Py(str)
         names = allNames
      projectTableEle[0].set('username',names[0])
      if info['parentprojectname'] is not None: projectTableEle[0].set('parentprojectname',info['parentprojectname'])

      jobSelectCom,jobSelectArgs = self.exportJobSelectionCommand(projectId=projectId,after=after,status=status,jobList=jobList)

      # Get list of jobs and data to populate xml table
      com = 'SELECT '+self.JOBITEMS[0]
      for item in self.JOBITEMS[1:]: com = com + ', '+item
      com = com + ' FROM Jobs'
      com += re.sub('Jobs\.','',jobSelectCom)
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.JOBTYPES)
      # Make the list of jobs before adding the jobs that provide the input files
      idx = self.JOBITEMS.index('jobnumber')
      for row in table:
        if not row[idx].count('.'): jobNumberList.append(row[idx])


      # If there are jobs which provide input files then add them to the table
      # Beware changing the job status to 'file holder'
      if inputFileFromJobList is not None and len(inputFileFromJobList)>0:
        com = 'SELECT '+self.JOBITEMS[0]
        for item in self.JOBITEMS[1:]: com = com + ', '+item
        com = com + ' FROM Jobs WHERE JobId IN ('
        for item in inputFileFromJobList: com = com + '?,'
        com = com[0:-1]+')'
        print('getTablesEtree',com,inputFileFromJobList)
        self.execute(com,inputFileFromJobList)
        inputJobTable =  self.fetchAll2PyList(self.JOBTYPES)
        idxStatus = self.JOBITEMS.index('status')
        for item in inputJobTable:
          item[idxStatus] = JOB_STATUS_FILE_HOLDER
          if item not in table: table.append(item)

      # Remove preceedingjobid if it is not in the jobList
      # Likely necessary for exporting jobs
      if jobList is not None:
        iPreJob = self.JOBITEMS.index('preceedingjobid')
        for row in table:
          if not jobList.count(row[iPreJob]): row[iPreJob] = None

      root.append(self.tableEtree('job',self.JOBITEMS,table))

      com = 'SELECT Files.'+self.FILEITEMS[0]
      for item in self.FILEITEMS[1:]: com = com + ', Files.'+item
      com = com + ' FROM Files INNER JOIN Jobs ON Files.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.FILETYPES)

      if inputFileList is not None and len(inputFileList)>0:
        com = 'SELECT '+self.FILEITEMS[0]
        for item in self.FILEITEMS[1:]: com = com + ', '+item
        com = com + ' FROM Files WHERE FileId IN ('
        for item in inputFileList: com = com + '?,'
        com = com[0:-1]+')'
        self.execute(com,inputFileList)
        fileTable =  self.fetchAll2PyList(self.FILETYPES)
        #print 'getTablesEtree fileTable from fileList', fileTable

        iJob = self.FILEITEMS.index('jobid')
        for row in fileTable:
          if row[iJob] is not None and  row[iJob] not in jobList:
            #row[iJob] = None
            table.append(row)


      root.append(self.tableEtree('file',self.FILEITEMS,table))

      com = 'SELECT FileUses.'+self.FILEUSEITEMS[0]
      for item in self.FILEUSEITEMS[1:]: com = com + ', FileUses.'+item
      com = com + ' FROM FileUses INNER JOIN Jobs ON FileUses.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.FILEUSETYPES)
      root.append(self.tableEtree('fileuse',self.FILEUSEITEMS,table))

      com = 'SELECT ImportFiles.'+self.IMPORTFILEITEMS[0]
      for item in self.IMPORTFILEITEMS[1:]: com = com + ', ImportFiles.'+item
      com = com + ' FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId = Files.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.IMPORTFILETYPES)

      if inputFileList is not None and len(inputFileList)>0:
        #print 'CDbApi.getTablesEtree getting import files from fileList',inputFileList
        com = 'SELECT '+self.IMPORTFILEITEMS[0]
        for item in self.IMPORTFILEITEMS[1:]: com = com + ', '+item
        com = com + ' FROM ImportFiles WHERE FileId IN ('
        for item in inputFileList: com = com + '?,'
        com = com[0:-1]+')'
        self.execute(com,inputFileList)
        importTable =  self.fetchAll2PyList(self.IMPORTFILETYPES)
        #print 'CDbApi.getTablesEtree importTable from specified files',importTable
        for item in importTable:
          if item not in table: table.append(item)



      root.append(self.tableEtree('importfile',self.IMPORTFILEITEMS,table))


      com = 'SELECT ExportFiles.'+self.EXPORTFILEITEMS[0]
      for item in self.EXPORTFILEITEMS[1:]: com = com + ', ExportFiles.'+item
      com = com + ' FROM ExportFiles INNER JOIN Files ON ExportFiles.FileId = Files.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.EXPORTFILETYPES)
      root.append(self.tableEtree('exportfile',self.EXPORTFILEITEMS,table))

      com = 'SELECT XData.'+self.XDATAITEMS[0]
      for item in self.XDATAITEMS[1:]: com = com + ', XData.'+item
      com = com + ' FROM XData INNER JOIN Jobs ON XData.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.XDATATYPES)
      tableEtree = self.tableEtree('xdata',self.XDATAITEMS[0:-1],table)
      # Add the actual xml representation in as etree
      for ii in range(len(table)):
        dataTree = etree.fromstring(table[ii][-1])
        dataTree.tag = 'xdataxml'
        tableEtree[ii].append(dataTree)
      root.append(tableEtree)

      com = 'SELECT Comments.'+self.COMMENTITEMS[0]
      for item in self.COMMENTITEMS[1:]: com = com + ', Comments.'+item
      com = com + ' FROM Comments INNER JOIN Jobs ON Comments.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.COMMENTTYPES)
      root.append(self.tableEtree('comment',self.COMMENTITEMS,table))

      com = 'SELECT ProjectComments.'+self.PROJECTCOMMENTITEMS[0]
      for item in self.PROJECTCOMMENTITEMS[1:]: com = com + ', ProjectComments.'+item
      com = com + ' FROM ProjectComments WHERE ProjectComments.ProjectId = ?'
      self.execute(com,[projectId,])
      table = self.fetchAll2PyList(self.PROJECTCOMMENTTYPES)
      self.execute("SELECT * FROM Users",())
      uList = self.fetchAll2Py(str)
      try:
          for t in table:
              commentUser = t[self.PROJECTCOMMENTITEMS.index('userid')]
              if not commentUser in uList:
                  t[self.PROJECTCOMMENTITEMS.index('userid')] = uList[0]
      except:
          exc_type, exc_value,exc_tb = sys.exc_info()[:3]
          sys.stderr.write(str(exc_type)+'\n')
          sys.stderr.write(str(exc_value)+'\n')
          
      root.append(self.tableEtree('projectcomment',self.PROJECTCOMMENTITEMS,table))

      com = 'SELECT JobKeyValues.'+self.JOBKEYVALUEITEMS[0]
      for item in self.JOBKEYVALUEITEMS[1:]: com = com + ', JobKeyValues.'+item
      com = com + ' FROM JobKeyValues INNER JOIN Jobs ON JobKeyValues.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.JOBKEYVALUETYPES)
      root.append(self.tableEtree('jobkeyvalue',self.JOBKEYVALUEITEMS,table))

      com = 'SELECT JobKeyCharValues.'+self.JOBKEYCHARVALUEITEMS[0]
      for item in self.JOBKEYCHARVALUEITEMS[1:]: com = com + ', JobKeyCharValues.'+item
      com = com + ' FROM JobKeyCharValues INNER JOIN Jobs ON JobKeyCharValues.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.JOBKEYCHARVALUETYPES)
      root.append(self.tableEtree('jobkeycharvalue',self.JOBKEYCHARVALUEITEMS,table))

      com = 'SELECT DISTINCT FileAssociationMembers.'+self.FILEASSOCIATIONMEMBERITEMS[0]
      for item in self.FILEASSOCIATIONMEMBERITEMS[1:]: com = com + ', FileAssociationMembers.'+item
      com = com + ' FROM FileAssociationMembers INNER JOIN Files ON FileAssociationMembers.FileId = Files.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.FILEASSOCIATIONMEMBERTYPES)
      root.append(self.tableEtree('fileassociationmember',self.FILEASSOCIATIONMEMBERITEMS,table))

      com = 'SELECT DISTINCT FileAssociations.'+self.FILEASSOCIATIONITEMS[0]
      for item in self.FILEASSOCIATIONITEMS[1:]: com = com + ', FileAssociations.'+item
      com = com + ' FROM FileAssociations INNER JOIN FileAssociationMembers ON FileAssociations.FileAssociationId = FileAssociationMembers.FileAssociationId INNER JOIN Files ON FileAssociationMembers.FileId = Files.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId '+jobSelectCom
      self.execute(com,jobSelectArgs)
      table = self.fetchAll2PyList(self.FILEASSOCIATIONTYPES)
      root.append(self.tableEtree('fileassociation',self.FILEASSOCIATIONITEMS,table))

      com = 'SELECT DISTINCT ProjectTags.'+self.PROJECTTAGITEMS[0]
      for item in self.PROJECTTAGITEMS[1:]: com = com + ', ProjectTags.'+item
      com = com + ' FROM ProjectTags WHERE ProjectTags.ProjectId = ?'
      self.execute(com,[projectId,])
      table = self.fetchAll2PyList(self.PROJECTTAGTYPES)
      root.append(self.tableEtree('projecttag',self.PROJECTTAGITEMS,table))

      com = 'SELECT DISTINCT Tags.'+self.TAGITEMS[0]
      for item in self.TAGITEMS[1:]: com = com + ', Tags.'+item
      com = com + ' FROM Tags INNER JOIN ProjectTags ON Tags.TagId = ProjectTags.TagId WHERE ProjectTags.ProjectId = ?'
      self.execute(com,[projectId,])
      table = self.fetchAll2PyList(self.TAGTYPES)
      root.append(self.tableEtree('tag',self.TAGITEMS,table))

      #print 'getTablesEtree jobNumberList',jobNumberList
      return root,jobNumberList,errReport

    def tableEtree(self,itemName,attributeList,table):
      from lxml import etree
      tableEle = etree.Element(itemName+'Table')
      for row in table:
        ele = etree.Element(itemName)
        ii = 0
        for attrib in attributeList:
          if row[ii] is not None: ele.set(attrib,str(row[ii]))
          ii = ii + 1
        tableEle.append(ele)
      return tableEle


    def getProjectEtree(self,projectId):
      from lxml import etree
      errReport = CErrorReport()
      root = etree.Element('project')
      info = self.getProjectInfo(projectId)
      for itemName in ['projectname','projectdirectory','projectcreated']:
        item = etree.Element(itemName)
        item.text = str(info[itemName])
        root.append(item)
      item = etree.Element('userid')
      userName = self.getUserInfo(info['userid'],'username')
      item = etree.Element('username')
      item.text = str(userName)
      root.append(item)

      jobBranch =  etree.Element('jobs')
      allJobsInfo = self.getProjectJobListInfo(projectId=projectId,order='ASC')
      if len(allJobsInfo)>0:
        jobBranch =  etree.Element('jobs')
        for jobInfo in allJobsInfo:
          jobItem,err = self.getJobEtree(jobId=jobInfo['jobid'],jobInfo=jobInfo)
          errReport.extend(err)
          jobBranch.append(jobItem)


      root.append(jobBranch)

      return root,errReport


    def getJobEtree(self,jobId=None,jobInfo={}):
      from lxml import etree
      errReport = CErrorReport()

      if len(jobInfo) == 0:
        jobInfo = self.getJobInfo(jobId)
      root = etree.Element('job')
      for itemName in self.JOBITEMS:
        if jobInfo[itemName] is not None:
          if itemName not in ['projectid','preceedingjobid']:
            item = etree.Element(itemName)
            item.text = str(jobInfo[itemName])
            root.append(item)
          elif itemName == 'preceedingjobid' and jobInfo['preceedingjobid'] is not None:
            preceedingNumber = self.getJobInfo(jobInfo['preceedingjobid'],'jobnumber')
            if preceedingNumber is not None:
              item = etree.Element('preceedingjobnumber')
              item.text = str(preceedingNumber)
              root.append(item)

      projectId = jobInfo['projectid']

      outputFileIds = self.getJobFiles(jobId=jobInfo['jobid'])
      if len(outputFileIds)>0:
        fileBranch =  etree.Element('outputFiles')
        for fileId in outputFileIds:
          fileItem,err = self.getFileEtree(fileId=fileId,role=FILE_ROLE_OUT,projectId=projectId)
          errReport.extend(err)
          fileBranch.append(fileItem)
        root.append(fileBranch)

      inputFileIds = self.getJobFiles(jobId=jobInfo['jobid'],role=FILE_ROLE_IN)
      if len(inputFileIds)>0:
        fileBranch =  etree.Element('inputFiles')
        for fileId in inputFileIds:
          try:
            fileItem,err = self.getFileEtree(fileId=fileId,role=FILE_ROLE_IN,projectId=projectId)
            errReport.extend(err)
            fileBranch.append(fileItem)
          except:
            print('Error creating entry for fileId',fileId)
        root.append(fileBranch)

      xDataIds = self.getXData(jobId=jobInfo['jobid'])
      #print 'getJobEtree xDataIds',xDataIds
      if len(xDataIds)>0:
        dataBranch = etree.Element('outputData')
        for xDataId in xDataIds:
          xDataItem = self.getXDataEtree(xDataId)
          if xDataItem is not None: dataBranch.append(xDataItem)
        root.append(dataBranch)

      return root,errReport

    def getJobOutputHtmlEtree(self,jobId=None,projectId=None):
      from lxml import etree
      root = etree.Element('div')
      root.set('id','result_data')
      outputFileIds = self.getJobFiles(jobId=jobId)
      for fileId in outputFileIds:
        root.append(self.getFileHtmlEtree(fileId,projectId=projectId))

      return root


    def getFileEtree(self,fileId,fileInfo={},role=FILE_ROLE_OUT,projectId=None):
      # !!!! Does not include jobId !! Assume this is only used in heirarchical representation below jobs
      # Create xml compatible with CDataFile xml
      from lxml import etree
      err = CErrorReport()
      if len(fileInfo) == 0: fileInfo = self.getFileInfo(fileId,mode= [ 'filename','relpath','annotation','fileclass','importid','filecontent','filesubtype'])
      #print 'getFileEtree fileInfo',fileInfo
      root = etree.Element(fileInfo['fileclass'])
      for itemName,xmlTag in  [ ['filename','baseName'],['relpath','relPath'],['annotation','annotation'],['filecontent','contentFlag'],['filesubtype','subType']]:
        if fileInfo[itemName] is not None:
          item = etree.Element(xmlTag)
          item.text = str(fileInfo[itemName])
          root.append(item)
      '''
      if  fileInfo['projectid'] is None:
        projectname = 'FULLPATH'
      else:
        projectname = self.getProjectInfo(fileInfo['projectid'],'projectname')
      item = etree.Element('project')
      item.text = str(projectname)
      root.append(item)
      '''
      item = etree.Element('dbFileId')
      item.text = str(fileId)
      root.append(item)
      if fileInfo['importid'] is not None:
        try:
          importInfo = self.getImportFileInfo(importId=fileInfo['importid'])
          importItem = etree.Element('importFile')
          root.append(importItem)
          for itemName,xmlTag in [['sourcefilename','sourceFilename'],['sourcefileid','sourceFileId'],['exportfileid','exportFileId'],['annotation','annotation']]:
            if importInfo[itemName] is not None:
              item = etree.Element(xmlTag)
              item.text = str(importInfo[itemName])
              importItem.append(item)
        except:
          err.append(self.__class__,255,str(fileInfo['importid']))
      try:
        exportInstances = self.getExportFileInstances(fileId=fileId)
      except:
        err.append(self.__class__,258,'FileID: '+str(fileId))
        exportInstances = []

      for exportInfo in exportInstances:
        try:
          exportItem = etree.Element('fileExport')
          root.append(exportItem)
          for itemName,xmlTag in [['exportid','exportFileId'],['exportfilename','exportFilename']]:
            if exportInfo[itemName] is not None:
                item = etree.Element(xmlTag)
                item.text = str(exportInfo[itemName])
                exportItem.append(item)
        except:
           err.append(self.__class__,256,str(exportInfo['exportid']))
      return root,err



    def getFileHtmlEtree(self,fileId,fileInfo={},role=FILE_ROLE_OUT,projectId=None):
      # Create xml compatible with CDataFile xml
      if projectId is None:
        print('NEED projectId in getFileHtmlEtree')
      from lxml import etree
      if len(fileInfo) == 0: fileInfo = self.getFileInfo(fileId,mode= [ 'filename','relpath','annotation','fileclass'])
      # Aiming for something like..
      # <object type="x-ccp4-widget/CPdbDataFile" id="result_file_pdb" width="600" height="30">
      # <param name="mimeTypeName" value="text/pdb" />  --- not necessary?
      # <param name="fullPath" value="/Users/lizp/Desktop/dev/ccp4i2/test/data/1df7.pdb" />
      # <param name="relPath" value="CCP4_JOBS/job_16" />
      # <param name="baseName" value="buccaneer_16_XYZOUT.pdb" />
      # </object>

      #print 'getFileEtree fileInfo',fileInfo
      root = etree.Element('object')
      root.set('type','x-ccp4-widget/C'+fileInfo['fileclass'])
      root.set('id','file_'+str(fileId))
      root.set('width','600')
      root.set('height','300')
      for itemName,xmlTag in  [ ['filename','baseName'],['relpath','relPath'],['annotation','annotation']]:
        if fileInfo[itemName] is not None:
          item = etree.Element('param')
          item.set('name',xmlTag)
          item.set('value',str(fileInfo[itemName]))
          root.append(item)
      #if fileInfo['projectid'] != projectId:

      projectname = self.getProjectInfo(projectId,'projectname')
      item = etree.Element('param')
      item.set('name','project')
      item.set('value',str(projectname))
      root.append(item)
      item = etree.Element('param')
      item.set('name','dbFileId')
      item.set('value',str(fileId))
      root.append(item)
      #print 'getFileHtmlEtree',etree.tostring(root,pretty_print=True)
      return root

    def getXDataEtree(self,xDataId,xDataInfo={}):
      from lxml import etree
      if len(xDataInfo)==0: xDataInfo = self.getXDataInfo(xDataId)
      if len(xDataInfo)==0:
        #print 'getXDataEtree no data',xDataId
        return None
      item = etree.Element(xDataInfo['xdataclass'])
      item.text = str(xDataInfo['xdataxml'])
      return item


    def exportProjectXml(self,projectId=None,fileName=None,recordExport=False,after=None,status=None,jobList=None,inputFileList=None,inputFileFromJobList=None):
      from core import CCP4File
      print('exportProjectXml',projectId,fileName)
      errReport = CErrorReport()
      try:
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        f.header.setCurrent()
        f.header.function.set('PROJECTDATABASE')
        f.header.projectId.set(projectId)
        f.header.projectName.set(self.getProjectInfo(projectId,'projectname'))
      except:
        errReport.append(self.__class__,257,fileName)
      #tree,err = self.getProjectEtree(projectId)
      tree,jobNumberList,err = self.getTablesEtree(projectId,after=after,status=status,jobList=jobList,inputFileList=inputFileList,inputFileFromJobList=inputFileFromJobList)
      #print 'CDbApi.exportProjectXml',tree,err
      errReport.extend(err)
      if recordExport:
        exportId = self.createProjectExport(projectId=projectId,projectExportAfter=after)
        table = self.getProjectExportInfo(projectExportId=exportId,returnList=True)
        exportTableEle=self.tableEtree('projectexport',('projectexportid','projectexporttime','projectexportafter'),table)
        #print 'exportProjectXml',exportId,table,exportTableEle
        tree.append(exportTableEle)
      try:
        f.saveFile(bodyEtree=tree)
      except CException as err:
        errReport.extend(err)
      except:
        err = CException(self.__class__,212)
        errReport.extend(err)
      #print 'CDbApi.exportProjectXml',errReport.report()
      return jobNumberList,errReport


    def updateProject(self,projectId=None,key=None,value=None):
      '''Update the Project'''
      if not self._userRole in [USER_ROLE_MANAGER,USER_ROLE_OWNER]:
        raise CException(self.__class__,106,'You are '+self._userName+' calling updateProject')

      key = key.lower()
      if not key in ['projectname','parentprojectid','projectdirectory','projectcreated','lastjobnumber','lastcleanuptime',
                     'i1projectname','i1projectdirectory','lastaccess']:
        raise CException(self.__class__,201,key)

      permission = self.getProjectPermission(projectId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      # Sanity check
      if key == 'parentprojectid':
        if value is not None and not isinstance(value,UUIDTYPE):
          raise CException(self.__class__,202,key)
        #The new parent must not be a child of this project
        children = self._getChildProjects(projectId,True)
        #print 'updateProject childProjects',children
        if children.count(value):
          raise CException(self.__class__,203,key)
      if key == 'projectcreated':
        if value is None: return
      if key == 'lastaccess' and value is None:
        value = self.currentTime()

      #self.setDiagnostic(True)
      self.execute("UPDATE Projects SET "+key+" = ? WHERE ProjectID = ?",(value,projectId,))
      #self.setDiagnostic(False)
      self.commit()
      self.projectUpdated.emit({'projectId':projectId,'key':key,'value':value})

    def updateProjectOwner(self,userName):
      pass

    def deleteProject(self,projectId=None,deleteChildren=False):
      '''Delete project and all jobs & permissions'''
      #self.setDiagnostic(True)
      self.removeTempTables()

      self.execute("DELETE FROM ExportFiles WHERE FileId IN (SELECT FileId FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectID = ?)",(projectId,))
      self.execute("DELETE FROM ImportFiles WHERE FileId IN (SELECT FileId FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectID = ?)",(projectId,))


     # Delete all jobs/files/comments/XData in the project
      print('Start deleting files etc')
      self.execute('DELETE FROM FileUses WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM Files WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM XData WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM Comments WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM JobKeyValues WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM JobKeyCharValues WHERE JobId IN (SELECT JobId FROM Jobs WHERE Jobs.ProjectId = ?)',(projectId,))
      self.execute('DELETE FROM ProjectImports WHERE ProjectId = ?',(projectId,))
      self.execute('DELETE FROM ProjectExports WHERE ProjectId = ?',(projectId,))
      self.execute('DELETE FROM ProjectTags WHERE ProjectId = ?',(projectId,))
      self.execute('DELETE FROM ProjectComments WHERE ProjectId = ?',(projectId,))

      # Delete Jobs and ProjectsUsersPermissions
      print('Start deleting jobs etc')
      self.execute("UPDATE Projects SET FollowFromJobId=NULL WHERE ProjectId = ?",(projectId,))
      self.execute("UPDATE Projects SET ParentProjectId=NULL WHERE ParentProjectId = ?",(projectId,))
      self.execute("DELETE FROM Jobs WHERE ProjectId = ?",(projectId,))
      self.execute("DELETE FROM ProjectsUsersPermissions WHERE ProjectID = ?",(projectId,))


      # Delete from Projects table too
      #print 'CDbApi.deleteProject',projectId,type(projectId)
      print('Start deleting project')
      self.execute("DELETE FROM Projects WHERE ProjectID = ?",(str(projectId),))
      print('Done deleting project')
      #self.setDiagnostic(False)

      self.commit()
      self.projectsListChanged.emit()
      self.projectDeleted.emit({'projectId':projectId})
      pass

    def getAllInfo(self,toTerm=False):

       retval = {}
       retval["pdbs"] = {}
       retval["mtzs"] = {}
       retval["mtzdiffs"] = {}
       projects = self.listProjects(toTerm);
       retval["projects"] = projects
       for project in projects:
         projectId,projectName,projectDir,dummy = project
         pdbs     = self.getJobsWithOutputFiles(projectId=projectId,fileType='chemical/x-pdb')
         mtzs     = self.getJobsWithOutputFiles(projectId=projectId,fileType='application/CCP4-mtz-map',subType=1)
         mtzdiffs = self.getJobsWithOutputFiles(projectId=projectId,fileType='application/CCP4-mtz-map',subType=2)
         retval["pdbs"][projectId] = pdbs
         retval["mtzs"][projectId] = mtzs
         retval["mtzdiffs"][projectId] = mtzdiffs
       return retval


    def listProjectNames(self,toTerm=False):
       self.execute("SELECT  ProjectName  FROM Projects ORDER BY ProjectCreated")
       projectList = self.fetchAll2Py(str)
       if toTerm:
         print('listProjectNames')
         for item in projectList: print (item)
       return projectList

    def listProjects(self,toTerm=False,order='date'):
      try:
        if order == 'date':
          self.execute("SELECT ProjectID, ProjectName, ProjectDirectory, ParentProjectID  FROM Projects ORDER BY ProjectCreated")
        elif order == 'name':
          self.execute("SELECT ProjectID, ProjectName, ProjectDirectory, ParentProjectID  FROM Projects ORDER BY ProjectName")
        elif order == 'access':
          self.execute("SELECT ProjectID, ProjectName, ProjectDirectory, ParentProjectID  FROM Projects ORDER BY LastAccess DESC")
        self._projectList = self.fetchAll2PyList([UUIDTYPE,str,str,UUIDTYPE])
      except:
        pass
      if toTerm:
         print('listProjects')
         for item in self._projectList: print (item)
      return self._projectList

    def nextJobNumber(self,projectId,childOf=None,reserve=1):
      if childOf is None:
        # Use Projects.LastJobNumber which hold last top level job number
        self.execute("SELECT LastJobNumber FROM Projects WHERE ProjectID = ? ",(projectId,))
        jobNoList = self.fetchAll2Py(int)
        #print 'nextJobNumber',projectId,jobNoList
        if len(jobNoList) != 1:
          raise CException(self.__class__,123,str(projectId))
        if reserve>0:
          jobNo = jobNoList[0] + reserve
          self.execute("UPDATE Projects SET LastJobNumber = ? WHERE ProjectID = ? ",(jobNo,projectId))
          self.commit()
          return  jobNoList[0] + 1
        else:
          return  jobNoList[0]
      else:
        # Find the subJobs and get the maximum subJobNumber
        self.execute("SELECT jobNumber from Jobs WHERE ProjectID = ? and JobID = ?",(projectId,childOf))
        jobNoList = self.fetchAll2Py(str)
        if len(jobNoList) != 1:
          raise CException(self.__class__,123,str(projectId))
        jobNoStr = jobNoList[0]
        self.execute("SELECT jobNumber from Jobs WHERE ProjectID = ? and ParentJobID = ?",(projectId,childOf))
        jobNoList = self.fetchAll2Py(str)
        maxNo = 0
        for item in jobNoList:
          lastNo = item.split('.')[-1]
          maxNo = max(maxNo,int(lastNo))
        return jobNoStr+'.'+str(maxNo+1)

    def getTaskNameLookup(self,projectId=None,jobId=None,extras=False):
      if extras:
        if jobId is not None:
          #Finding subjobs follows code example from
          #https://www.sqlite.org/lang_with.html  (works_for_alice example)
          # Using bindings did not work - hum?
          self.execute("""
            WITH RECURSIVE subJobs(n) AS (
              VALUES ( '"""+jobId+"""' )
              UNION
              SELECT Jobs.JobId FROM Jobs, subJobs
              WHERE Jobs.ParentJobId=subJobs.n
            )
            SELECT Jobs.JobNumber,Jobs.TaskName,Jobs.JobId,Jobs.Status,Jobs.FinishTime FROM Jobs WHERE Jobs.JobId IN subJobs;""")
        else:
          self.execute('SELECT JobNumber,TaskName,JobId,Status,FinishTime FROM Jobs WHERE ProjectID = ?',(projectId,))
        rv = self.fetchAll2PyList([str,str,UUIDTYPE,int,float])
      else:
        self.execute('SELECT JobNumber,TaskName FROM Jobs WHERE ProjectID = ?',(projectId,))
        rv = self.fetchAll2PyList([str,str])
      return rv


    def getJobId(self,projectName=None,jobNumber=None,projectId=None):
      if projectId is None: projectId = self.getProjectId(projectName)
      if jobNumber is None:
        raise CException(self.__class__,134,str(projectName)+':'+str(jobNumber))
      u = (projectId,jobNumber,)
      self.execute('SELECT JobID FROM Jobs WHERE ProjectID = ? AND JobNumber = ?',u)
      rv = self.fetchAll2Py(UUIDTYPE)
      #print('CDbApi.getJobId',projectId,jobNumber,rv)
      if len(rv)==0:
        raise CException(self.__class__,134,str(projectId)+':'+str(jobNumber))
      else:
        return rv[0]

    def getJobsInRange(self,projectId=None,jobIdRange=None,jobNumberRange=None):
      #print 'CDbApi.getJobsInRange',jobIdRange,jobNumberRange
      args = [projectId]
      if jobIdRange is not None:
        args.extend(jobIdRange)
        self.execute('SELECT CreationTime FROM Jobs WHERE ProjectId =? AND JobID IN (?,?)',args)
      elif jobNumberRange is not None:
        args.extend(jobNumberRange)
        self.execute('SELECT CreationTime FROM Jobs WHERE  ProjectId =? AND JobNumber IN (?,?)',args)

      timeRange = self.fetchAll2Py(float)
      #print 'CDbApi.getJobsInRange timeRange',timeRange
      if len(timeRange)!=2: return []
      if timeRange[0]<timeRange[1]:
        self.execute('SELECT JobId FROM Jobs WHERE CreationTime >= ? AND CreationTime <= ? AND ParentJobId IS NULL',timeRange)
      else:
        self.execute('SELECT JobId FROM Jobs WHERE CreationTime <= ? AND CreationTime >= ?  AND ParentJobId IS NULL',timeRange)
      rv = self.fetchAll2Py(UUIDTYPE)
      return rv


    def createJob(self,projectId,taskName,parentJobId=None,jobTitle=None,status=JOB_STATUS_PENDING,jobNumber=None,taskVersion=None,userAgent=None,userName=None):
      '''Create a job record for project projectId'''

      #print 'CDbApi.createJob jobTitle',jobTitle


      permission = self.getProjectPermission(projectId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112,str(projectId))


      if jobNumber is None:
        jobNumber = self.nextJobNumber(projectId,childOf=parentJobId)
      time = self.currentTime()
      if userAgent is None: userAgent = self.iAmUserAgent()
      if userName is not None:
        userId = self.getUserId(userName)
      else:
        userId = self.getUserId(self._userName)
      #print 'createJob userId',userId

      # Does it already exist?
      try:
        jobId = self.getJobId(projectId=projectId,jobNumber=jobNumber)
      except:
        jobId = None

      if jobId is not None:
        # Job already exists - liable to be a developer rerunning jobs
        args = [projectId,parentJobId,jobNumber,userAgent,taskName,jobTitle,time,JOB_STATUS_PENDING,userId,jobId]
        com = "UPDATE Jobs SET ProjectID=?, parentJobId=?, JobNumber=?, userAgent=?, TaskName=?, JobTitle=?, CreationTime=?, Status=?, UserId= ? WHERE JobId= ?"

      else:
        jobId = self.uniqueId(table='Jobs',identifier='JobID')
        args = [jobId,projectId,parentJobId,jobNumber,userAgent,taskName,jobTitle,time,JOB_STATUS_PENDING,userId]
        com = "INSERT INTO Jobs (JobID,ProjectID,parentJobId,JobNumber,userAgent,TaskName,JobTitle,CreationTime,Status,UserId) VALUES (?,?,?,?,?,?,?,?,?,?)"
      self.execute(com,args)
      self.commit()
      self.jobCreated.emit({'jobId':jobId,'projectId':projectId,'parentJobId':parentJobId})

      # Must commit before emiting! And better to signal job creation
      # before job finish
      if status != JOB_STATUS_PENDING:
        self.updateJobStatus(jobId,status)

      return jobId


    def updateJobStatus(self,jobId=None,status=None,jobNumber=None,projectName=None):
      '''Update job status '''
      #print 'CDbApi.updateJobStatus',jobId,status,type(jobId)
      if status is None:
          raise CException(self.__class__,116)
      if status in JOB_STATUS_TEXT:
        status = JOB_STATUS_TEXT.index(status)

      if jobId is None:
        jobId = self.getJobId(jobNumber=jobNumber,projectName=projectName)

      projectId,permission = self.getJobPermission(jobId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      if status == JOB_STATUS_QUEUED:
        # A control file needs to exist
        controlFile = self._makeJobFileName(jobId=jobId,mode='JOB_INPUT')
        if controlFile is None:
          raise CException(self.__class__,132,str(controlFile))
        elif not os.path.exists(controlFile):
          raise CException(self.__class__,132,str(controlFile))

      self.execute("UPDATE Jobs SET status = ? WHERE JobID = ?",(status,jobId,))

      #if status in [JOB_STATUS_RUNNING]:
      #  self.execute("UPDATE Jobs SET CreationTime = ? WHERE JobID = ?", (self.currentTime(),jobId,))
      #  self.commit()
      if status in (JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_INTERRUPTED,JOB_STATUS_TO_DELETE,JOB_STATUS_UNSATISFACTORY):
        self.execute("UPDATE Jobs SET FinishTime = ? WHERE JobID = ?", (self.currentTime(),jobId,))
        #print 'CDbApi.updateJobStatus finished',jobId,status
        self.commit()
        self.jobFinished.emit({'jobId':jobId,'projectId':projectId,'status':status})
      else:
        self.commit()

      self.jobUpdated.emit({'jobId':jobId,'projectId':projectId,'key':'status','value':status})
      self.jobStatusUpdated.emit({'jobId':jobId,'status':JOB_STATUS_TEXT[status]})

    def updateJob(self,jobId=None,key=None,value=None):

      if key is None: raise CException(self.__class__,169,'None')
      key = key.lower()
      if not key in self.JOBITEMS:
        raise CException(self.__class__,169,key)

      projectId,permission = self.getJobPermission(jobId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      # Sensible value?
      if key == 'preceedingjobid':
        #print 'CDbApi.updateJob preceedingjobid',value,type(value)
        # Should be in same project..
        if value is None or len(value)==0:
          return
        else:
          self.execute('SELECT ProjectId from Jobs WHERE JobId = ?',(UUIDTYPE(value),))
          pList = self.fetchAll2Py(UUIDTYPE)
          if len(pList)>0:
            preceedingJobProject = pList[0]
          else:
            #raise CException(self.__class__,175,'No project')
            return
          if preceedingJobProject != projectId or value == jobId:
            raise CException(self.__class__,175,str(value))
      elif key == 'evaluation':
        if isinstance(value,int) and value>=0 and value <len(JOB_EVALUATION_TEXT):
          pass
        elif JOB_EVALUATION_TEXT.count(value):
          value = JOB_EVALUATION_TEXT.index(value)
        else:
          raise CException(self.__class__,190,'Key: '+str(key)+'  Value: '+str(value))
      if value is None or value == 0 or value == '':
        self.execute("UPDATE Jobs SET "+key+" = NULL WHERE JobID = ?",(jobId,))
      else:
        self.execute("UPDATE Jobs SET "+key+" = ? WHERE JobID = ?",(value,jobId,))
      self.commit()

      self.jobUpdated.emit({'jobId':jobId,'projectId':projectId,'key':key,'value':value})


    def deleteJob(self,jobId=None,projectName=None,jobNumber=None,deleteChildren=False):
      '''Delete a record in the Jobs table'''
      #print 'CDbApi.deleteJob',jobId
      if jobId is None:
        jobId = self.getJobId(projectName=projectName,jobNumber=jobNumber)

      projectId,permission = self.getJobPermission(jobId)
      if permission < PRIVILEGE_DELETE:
        raise CException(self.__class__,112)

      jobInfo = self.getJobInfo(jobId=jobId,mode=['projectid','jobnumber'])

      #if deleteChildren:
      self._deleteJobTree(jobId)
      #else:
      #self._deleteOneJob(jobId)
      self.resetLastJobNumber(projectId=projectId,commit=False)
      self.commit()
      self.jobDeleted.emit({'jobId':jobId,'projectId':jobInfo['projectid'],'jobNumber':jobInfo['jobnumber'],'deleteChildren':deleteChildren})

    def _deleteJobTree(self,jobId):
      self.execute("SELECT JobID FROM Jobs WHERE ParentJobID = ?",(jobId,))
      childJobs = self.fetchAll2Py(UUIDTYPE)
      for job in childJobs:
        self._deleteJobTree(job)
      self._deleteOneJob(jobId)

    def _deleteOneJob(self,jobId,saveJob=False,saveImportFiles=False):
      #print '_deleteOneJob',jobId
      #self.setDiagnostic(True)
      self.execute("UPDATE Jobs SET PreceedingJobId = ? WHERE PreceedingJobId = ?",(None,jobId))
      self.execute("UPDATE Projects SET FollowFromJobId = ? WHERE FollowFromJobId = ?",(None,jobId))

      self.execute("DELETE FROM FileUses WHERE JobID = ?",(jobId,))
      self.execute("DELETE FROM XData WHERE JobID = ?",(jobId,))
      self.execute("DELETE FROM Comments WHERE JobID = ?",(jobId,))
      self.execute("DELETE FROM JobKeyValues WHERE JobID = ?",(jobId,))
      self.execute("DELETE FROM JobKeyCharValues WHERE JobID = ?",(jobId,))
      if not saveImportFiles:
        self.execute("DELETE FROM ImportFiles WHERE FileId IN (SELECT FileId FROM Files WHERE Files.JobId= ? )",(jobId,))
      self.execute("DELETE FROM ExportFiles WHERE FileId IN (SELECT FileId FROM Files WHERE Files.JobId= ? )",(jobId,))
      self.execute("DELETE FROM FileUses WHERE FileId IN (SELECT FileId FROM Files WHERE Files.JobId= ? )",(jobId,))
      # If we zap one of the files in a file association do we zap the whole file association???
      # Maybe count number of members - could also change association type if one memeber lost?
      #self.execute("DELETE FROM FileAssociations WHERE FileAssociationId IN (SELECT FileAssociationId FROM FileAssociationMembers WHERE FileId IN (SELECT FileId FROM Files WHERE Files.JobId= ? ))",(jobId,))
      self.execute("DELETE FROM FileAssociationMembers WHERE FileId IN (SELECT FileId FROM Files WHERE Files.JobId= ? )",(jobId,))

      if saveImportFiles:
        self.execute("DELETE FROM Files WHERE JobID = ? AND FileId NOT IN (SELECT Files.FileId FROM Files INNER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId AND Files.JobId = ?)",(jobId,jobId,))
      else:
        self.execute("DELETE FROM Files WHERE JobID = ?",(jobId,))
      #self.execute("UPDATE Files SET JobId = NULL WHERE JobID = ?",(jobId,))
      self.execute("DELETE FROM ServerJobs WHERE JobID = ?",(jobId,))

      if not saveJob:
        self.execute("DELETE FROM Jobs WHERE JobID = ?",(jobId,))
      #print 'done _deleteOneJob'
      #self.setDiagnostic(False)


    def getJobAncestors(self,jobId=None):
      jobPath = [jobId]
      while 1:
        self.execute('SELECT parentjobid FROM Jobs WHERE jobid=?',(jobPath[0],))
        rv = self.fetchAll2Py(UUIDTYPE)
        if len(rv)==0 or rv[0] is None:
          return jobPath
        else:
          jobPath.insert(0,rv[0])

    def getJobInfo(self,jobId=None,mode='all',projectName=None,jobNumber=None,returnType=None):
      # Retrieves details for one job

      if jobId is None:
        jobId = self.getJobId(projectName=projectName,jobNumber=jobNumber)

      try:
        projectId,permission = self.getJobPermission(jobId)
        if permission < PRIVILEGE_READ:
          raise CException(self.__class__,135)
      except:
        pass

      xtrasList =[]

      if mode == 'all':
          itemList = []
          itemList.extend(self.JOBITEMS)
      elif not isinstance(mode,list):
          itemList = [mode.lower()]
      else:
        itemList = []
        for item in mode: itemList.append(item.lower())
      # 'projectname' needs to be last in list to avoid problems in jobInfoDict
      retList = []
      for item in itemList:
        if item == 'projectname':
          if not itemList.count('projectid'): retList.append('projectid')
          xtrasList.append('projectname')
        elif item == 'runtime':
          if not itemList.count('creationtime'): retList.append('creationtime')
          if not itemList.count('finishtime'): retList.append('finishtime')
          xtrasList.append('runtime')
        elif item == 'parentjobnumber':
          if not itemList.count('parentjobid'): retList.append('parentjobid')
          xtrasList.append('parentjobnumber')
        elif item in ['childjobs','descendentjobs']:
          xtrasList.append(item)
        elif item in ['performanceclass']:
          if not itemList.count('taskname'): retList.append('taskname')
          xtrasList.append(item)
        elif item in ['performance']:
          xtrasList.append(item)
        elif item in self.JOBITEMS:
          retList.append(item)

      if len(retList)>0:
        selcom = "SELECT "
        typeConv = []
        for item in retList:
          selcom = selcom + item +', '
          typeConv.append(self.JOBTYPES[self.JOBITEMS.index(item)])
        selcom = selcom[0:-2] +  " FROM Jobs WHERE JobID=?"
        self.execute(selcom,(jobId,))
        rv = self.fetchAll2PyList(typeConv)
        #print 'getProjectJobListInfo rv',rv
        if len(rv) != 1:
          raise CException(self.__class__,137,'JobID:'+str(jobId))
        ret = self.jobInfoDict(retList,rv[0],extras=xtrasList)
      else:
        ret = {}


      if xtrasList.count('childjobs'): ret['childjobs'] = self.getChildJobs(jobId)
      if xtrasList.count('descendentjobs'): ret['descendentjobs'] = self.getChildJobs(jobId,descendents=True)
      if xtrasList.count('performanceclass'):
        ret['performanceclass'] = self.getJobPerformanceClass(jobId,taskName=ret.get('taskname'))
      if xtrasList.count('performance'):
        ret['performance'] = self.getJobPerformance(jobId=jobId)
      if (len(retList)+len(xtrasList))==1 and returnType != dict:
        # Just return the value if only one item requested
        if len(retList)==1:
          return ret[retList[0]]
        else:
          return ret[xtrasList[0]]
      else:
        # Return multiple items as a dictionary
        return ret

    def getJobPerformanceClass(self,jobId=None,taskName=None):
      from core import CCP4PerformanceData
      selcom = 'SELECT XDataClass,XDataXml FROM XData WHERE JobId=? AND XDataClass in '+CCP4PerformanceData.performanceIndicatorClasses()
      self.execute(selcom,(jobId,))
      rv = self.fetchAll2PyList([str,str])
      if len(rv)>0:
        try:
          from lxml import etree
          from core import CCP4DataManager
          obj = CCP4DataManager.DATAMANAGER().getClass(rv[0][0])()
          obj.setEtree(etree.fromstring(rv[0][1]))
          return obj
        except Exception as e:
          print(e)
      else:
        from core import CCP4TaskManager
        cls = CCP4TaskManager.TASKMANAGER().getPerformanceClass(taskName)
        if cls is None: return None
        perfDict = self.getJobPerformance(jobId=jobId)
        try:
          obj = cls()
          obj.set(perfDict)
          #print 'CDbApi.getJobPerformanceClass',perfDict,cls,obj
          return obj
        except Exception as e:
          print(e)
        return None

    def getJobPerformance(self,jobId=None):
      self.execute('SELECT KeyTypeId,Value FROM JobKeyValues WHERE jobId = ?',(jobId,))
      rv = self.fetchAll2PyList([int,float])
      self.execute('SELECT KeyTypeId,Value FROM JobKeyCharValues WHERE jobId = ?',(jobId,))
      rv.extend(self.fetchAll2PyList([int,str]))
      if len(rv)==0:
        return None
      else:
        ret = {}
        for KeyTypeId,Value in rv:
          ret[KEYTYPELIST[KeyTypeId][1]] = Value
        #print 'CDbApi.getJobPerformance',ret
        return ret

    def getJobFilesInfo(self,jobId=None,jobParamName=None,fileTypeId=None,input=False):

      if input:
        com = '''SELECT Files.FileId,Files.Filename,Files.Annotation,Files.JobParamName,Files.FileTypeId,Files.FileSubType,
                  Files.FileContent,Files.JobId,ImportFiles.ImportId from Files
                INNER JOIN FileUses ON Files.FileId = FileUses.FileId LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId
                WHERE FileUses.JobId= ?'''
        args=[jobId]
        if jobParamName is not None:
          com += 'AND FileUses.JobParamName LIKE ?'
          args.append(jobParamName+'%')
        if fileTypeId is not None:
          com += 'AND Files.FileTypeId=?'
          args.append(fileTypeId)
      else:
        com = '''SELECT Files.FileId,Files.Filename,Files.Annotation,Files.JobParamName,Files.FileTypeId,Files.FileSubType,
               Files.FileContent,Files.JobId,ImportFiles.ImportId from Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId
               WHERE JobId= ? '''
        args=[jobId]
        if jobParamName is not None:
          com += 'AND JobParamName LIKE ?'
          args.append(jobParamName+'%')
        if fileTypeId is not None:
          com += 'AND FileTypeId=?'
          args.append(fileTypeId)


      self.execute(com,args)
      rv = self.fetchAll2PyList([str,str,str,str,int,int,int,str,str])
      ret = []

      for item in rv:
        ret.append( { 'fileId' : item[0], 'fileName' : item[1], 'annotation' : item[2], 'jobParamName' : item[3],
                      'fileTypeId' : item[4], 'subType': item[5], 'fileContent' :item[6], 'jobId' :item[7],'importId' :  item[8]   } )
        ret[-1]['fullPath'] = self.getFullPath(fileId=ret[-1]['fileId'])
        # Since this method is utility for export_mtz code put some useful extra info here
        ret[-1].update(self.getJobInfo(ret[-1]['jobId'],['jobnumber','taskname']))

      return ret



    def getProjectJobListInfo(self,mode='all',projectId=None,projectName=None,topLevelOnly=False,order=None,
                                    firstJob=None,lastJob=None,maxJobs=None,jobStatus=[]):
      # Return list of dicts containing data specified by mode

      #print 'CDbApi.getProjectJobListInfo',projectId,projectName,mode
      if order is not None:
        order = order.upper()
        if not order in ['ASC','DESC']: order = 'ASC'

      if projectId is None: projectId = self.getProjectId(projectName)
      permission = self.getProjectPermission(projectId)
      if permission < PRIVILEGE_READ:
        raise CException(self.__class__,135)

      if mode == 'all':
          itemList = []
          itemList.extend(self.JOBITEMS)
          itemList.append('childjobs')
      elif not isinstance(mode,list):
          itemList = [mode.lower()]
      else:
        itemList = []
        for item in mode: itemList.append(item.lower())
      if ('childjobs' in mode or order is not None or \
          firstJob is not None or lastJob is not None) and not 'jobid' in itemList : itemList.append('jobid')
      if 'projectname' in itemList and not 'projectid' in itemList: itemList.append('projectid')
      if ('performance' in itemList or 'performanceclass' in itemList) and not 'jobid' in itemList : itemList.append('jobid')
      if 'performanceclass' in itemList and not 'taskname' in itemList : itemList.append('taskname')


      selcom = "SELECT "
      typeConv = []
      for item in itemList:
        if item == 'parentjobnumber':
          selcom = selcom + 'parentjobid, '
          typeConv.append(int)
        elif item == 'childjobs':
          pass
        else:
          try:
            ii = self.JOBITEMS.index(item)
            selcom = selcom + item +', '
            typeConv.append(self.JOBTYPES[ii])
          except:
            pass
      selcom = selcom[0:-2] + " FROM Jobs"

      # For one job or all jobs? (??Ranges)
      selcom = selcom + " WHERE ProjectID=?"

      if len(jobStatus)>0:
        selcom = selcom + " AND (Status = "+str(jobStatus[0])
        for item in jobStatus[1:]:
          selcom = selcom + ' OR Status ='+str(item)
        selcom = selcom +')'
        #print 'getProjectJobListInfo',selcom

      if topLevelOnly:
        selcom = selcom + ' AND ParentJobId IS NULL'
      if order is not None:
        selcom = selcom + ' ORDER BY CreationTime '+order
      if firstJob is not None:
        pass
      if lastJob is not None:
        pass


      #print 'getJobListInfo selcom',selcom,typeConv
      self.execute(selcom,(projectId,))

      rv = self.fetchAll2PyList(typeConv)
      #print 'getProjectJobListInfo rv',rv
      if maxJobs is not None and len(rv)>maxJobs: rv = rv[0:maxJobs]

      """
      if 'performanceclass' in itemList:
        from core import CCP4PerformanceData,CCP4DataManager
        from lxml import etree
        DM = CCP4DataManager.DATAMANAGER()
        self.execute('SELECT JobId,XDataClass,XDataXml FROM XData WHERE JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=?) AND XDataClass IN '+CCP4PerformanceData.performanceIndicatorClasses(),(projectId,))
        perfList = self.fetchAll2PyList([UUIDTYPE,str,str])
        perfDict = {}
        for jid,cls,data in perfList:
          try:
            perfDict[jid] = DM.getClass(cls)()
            perfDict[jid].setEtree(etree.fromstring(data))
          except:
            print 'Error in CDbApi.getProjectJobListInfo creating Performace Indicator data object'
        """


      perfDict = {}
      if 'performance' in itemList or 'performanceclass' in itemList:
          self.execute('SELECT JobId,KeyTypeId,Value FROM JobKeyValues WHERE JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=?  )',(projectId,))
          perfList = self.fetchAll2PyList([UUIDTYPE,int,float])
          self.execute('SELECT JobId,KeyTypeId,Value FROM JobKeyCharValues WHERE JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=? )',(projectId,))
          perfList.extend(self.fetchAll2PyList([UUIDTYPE,int,str]))

          for jid,keyId,value in perfList:
            try:
              if jid not in perfDict: perfDict[jid] = {}
              perfDict[jid][KEYTYPELIST[keyId][1]] = value
            except:
              print('Error in CDbApi.getProjectJobListInfo creating Performace Indicator data object')
          #print 'getProjectJobListInfo.perfDict', perfList, perfDict

      retList = []
      for jobValue in rv:
        retList.append(self.jobInfoDict(itemList,jobValue,performanceDict=perfDict))
      #print 'getJobListInfo retList',retList

      if 'childjobs' in itemList:
        for item in retList: item['childjobs'] = []
        maxIdx = len(retList)-1
        #self.execute("SELECT JobID,ParentJobId FROM Jobs WHERE ProjectID=? AND ParentJobId IS NOT NULL ORDER BY ParentJobId",(projectId,))
        self.execute("SELECT j1.JobID,j1.ParentJobId FROM Jobs j1 INNER JOIN Jobs j2 ON j1.ParentJobId = j2.jobId WHERE j1.ProjectID=? AND j1.ParentJobId IS NOT NULL ORDER BY j2.CreationTime",(projectId,))
        rv = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE])
        idx = 0
        for job,parent in rv:
          #print 'getProjectJobListInfo  job,parent', job,parent,idx,retList[idx]['jobid']
          while (retList[idx]['jobid'] != parent) and idx<maxIdx:
            idx = idx + 1
          retList[idx]['childjobs'].append(job)

        #print 'getProjectJobListInfo',retList

      return retList


    def getProjectJobList(self,projectId=None,topLevelOnly=False,maxLength=None):
      permission = self.getProjectPermission(projectId)
      if permission < PRIVILEGE_READ:
        raise CException(self.__class__,135)

      com = 'SELECT JobId FROM Jobs WHERE ProjectID = ?'
      if topLevelOnly: com = com + ' AND ParentJobId IS NULL'
      com = com +  ' ORDER BY CreationTime DESC'

      self.execute(com,[projectId])
      jobIdList= self.fetchAll2Py(UUIDTYPE)
      if maxLength is None or len(jobIdList)<maxLength:
        return jobIdList
      else:
        return jobIdList[0:maxLength]

    def jobInfoDict(self,itemList,jobValue,extras=[],performanceDict={}):
      # Utility to convert job info to Python dict and nicer forms
      #print 'CDbApi.jobInfoDict',itemList,jobValue
      from core import CCP4Annotation,CCP4TaskManager
      ii = -1
      ret = {}
      for item in itemList:
        ii = ii + 1
        # Skip converting time/data format - leave it to gui preferences
        #if item in ['creationtime','finishtime']:
        #  if jobValue[ii] is None:
        #    ret[item] = None
        #  else:
        #    ret[item] = str(CCP4Annotation.CTime(value=jobValue[ii]))
        if item == 'status':
          try:
            ret[item] = JOB_STATUS_TEXT[jobValue[ii]]
          except:
            ret[item] = 'Unknown'
        elif item == 'evaluation':
          if jobValue[ii] is None:
            ret[item] = 'Unknown'
          else:
            ret[item] = JOB_EVALUATION_TEXT[jobValue[ii]]
        elif item == 'useragent':
          ret[item] = USER_AGENT_TEXT[jobValue[ii]]
        elif item == 'parentjobnumber':
          ret[item] = self.getJobInfo(jobValue[ii],'jobnumber')
        elif item in ['performance']:
          ret[item] = performanceDict.get(jobValue[itemList.index('jobid')],None)
        elif item in ['performanceclass']:
          ret[item] = self.getJobPerformanceClass(jobValue[itemList.index('jobid')],jobValue[itemList.index('taskname')])
          cls = CCP4TaskManager.TASKMANAGER().getPerformanceClass(jobValue[itemList.index('taskname')])
          if cls is None:
            ret[item] = None
          else:
            try:
              ret[item] = cls()
            except Exception as e:
              print('Error in CDbApi.jobInfoDict creating Performance data object',e)
              ret[item] = None

        elif item in ['childjobs']:
          pass
        else:
          try:
            ret[item] = jobValue[ii]
          except:
            print('Error in jobInfoDict',item,ii,jobValue)

      for item in extras:
        if item == 'runtime':
          finishtime = jobValue[itemList.index('finishtime')]
          if finishtime is not None:
            ret[item] = finishtime-jobValue[itemList.index('creationtime')]
          else:
            ret[item] = 0.0
        elif item == 'projectname':
          ret[item] = self.getProjectInfo(projectId=jobValue[itemList.index('projectid')],mode='projectname')
        elif item == 'parentjobnumber':
          if jobValue[itemList.index('parentjobid')] is not None:
            ret[item] = self.getJobInfo(jobId=jobValue[itemList.index('parentjobid')],mode='jobnumber')
          else:
            ret[item] = None

      #print 'jobInfoDict',extras,ret
      return ret

    def getChildJobs(self,jobId,descendents=False,details=False):
      self.execute("SELECT JobId FROM Jobs WHERE ParentJobId = ? ORDER BY CreationTime",(jobId,))
      childList = self.fetchAll2Py(UUIDTYPE)
      if not descendents:
        if details:
          detailsList = []
          for child in childList:
            ret = self.getJobInfo(jobId=child,mode=['taskname','jobnumber'])
            detailsList.append([ret['jobnumber'],child,ret['taskname']])
          detailsList.sort()
          return detailsList
        else:
          return childList
      retList = []
      for item in childList:
        cL = self.getChildJobs(item,descendents=True,details=details)
        retList.append( [item,cL] )
      return retList


    def fileTypeCommand(self,selcom='',args=[],subType=None,contentFlag=None):
      if subType is not None and subType is not NotImplemented:
        if not isinstance(subType,list):
          if subType == 0:
            selcom += ' AND (Files.FileSubType IS NULL OR Files.FileSubType = ?)'
          else:
            selcom += ' AND Files.FileSubType = ?'
          args.append(subType)
        else:
          if 0 in subType:
            selcom += ' AND (Files.FileSubType IS NULL OR Files.FileSubType IN ('
          else:
            selcom += ' AND Files.FileSubType IN ('
          for subTypeItem in subType:
            args.append(subTypeItem)
            selcom +='?,'
          selcom = selcom[0:-1]+')'
          if 0 in subType: selcom += ')'
      if contentFlag is not None and contentFlag is not NotImplemented:
        if not isinstance(contentFlag,list):
          selcom += ' AND Files.FileContent = ?'
          args.append(contentFlag)
        else:
          selcom += ' AND Files.FileContent IN ('
          for contentFlagItem in contentFlag:
            args.append(contentFlagItem)
            selcom +='?,'
          selcom=selcom[0:-1]+ ')'
      return selcom,args

    def getFileByJobContext(self,contextJobId=None,fileType=None,projectId=None,subType=None,contentFlag=None):

      #print 'CDbApi.getFileByContext',contextJobId,fileType,subType,contentFlag

      if fileType in FILETYPES_TEXT: fileType = FILETYPES_TEXT.index(fileType)

      if fileType in ['application/CCP4-mtz',4]:
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE  ImportFiles.ImportId IS NULL AND Files.JobId = ? AND Files.FiletypeID IN (10,11,12,13)'
        args = [contextJobId]
      elif fileType in ['application/CCP4-mtz-mini',16]:
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE  ImportFiles.ImportId IS NULL AND Files.JobId = ? AND Files.FiletypeID IN (10,11,12,13)'
        args = [contextJobId]
      else:
        # Is there an output file from the job?
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE ImportFiles.ImportId IS NULL AND Files.JobId = ? AND Files.FiletypeID = ?'
        args = [contextJobId,fileType]
      selcom,args = self.fileTypeCommand(selcom,args,subType,contentFlag)
      #print 'getFileByJobContext',selcom,args
      self.execute(selcom,args)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)>1 and  self.getJobInfo(contextJobId,'taskname') == 'coot_rebuild':
        #print 'getFileByJobContext reversing coot output'
        rv.reverse()
      if len(rv)>0: return rv

      # Look for imported files
      if fileType in ['application/CCP4-mtz',4]:
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE  ImportFiles.ImportId IS NOT NULL AND Files.JobId = ? AND Files.FiletypeID IN (10,11,12,13)'
        args = [contextJobId]
      elif fileType in ['application/CCP4-mtz-mini',16]:
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE  ImportFiles.ImportId IS NOT NULL AND Files.JobId = ? AND Files.FiletypeID IN (10,11,12,13)'
        args = [contextJobId]
      else:
        selcom = 'SELECT Files.FileId FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE ImportFiles.ImportId IS NOT NULL AND Files.JobId = ? AND Files.FiletypeID = ?'
        args = [contextJobId,fileType]
      selcom,args = self.fileTypeCommand(selcom,args,subType,contentFlag)
      #print 'getFileByJobContext imported files',selcom,args
      self.execute(selcom,args)
      rv = self.fetchAll2Py(UUIDTYPE)
      #print 'getFileByJobContext imported files',rv
      if len(rv)>0: return rv

      # Search FileUses for either an output or input file
      selcom = 'SELECT FileUses.FileId,FileUses.RoleID FROM FileUses INNER JOIN Files ON FileUses.FileId=Files.FileId WHERE FileUses.JobId = ? AND Files.FiletypeID = ?'
      args = [contextJobId,fileType]
      selcom,args = self.fileTypeCommand(selcom,args,subType,contentFlag)
      self.execute(selcom,args)
      rv = self.fetchAll2PyList((UUIDTYPE,int))
      # 'CDbApi.getFileByContext FileUses',rv
      if len(rv)==1:
        # One hit - don't care if it was input or output
        return rv[0]
      elif len(rv)>1:
        # More than one file of appropriate type - precedence to output files
        outList = []
        inList = []
        for fileId,role in rv:
          if role == FILE_ROLE_IN:
            inList.append(fileId)
          else:
            outList.append(fileId)
        # Beware for coot_rebuild want output files in reverse order
        if len(outList)>1 and self.getJobInfo(contextJobId,'taskname') == 'coot_rebuild':
          #print 'getFileByJobContext reversing outlist'
          outList.reverse()

        if len(outList)>0: return outList
        if len(inList)>0: return inList

      # The job has neither input nor output file of the appropriate type
      # We need to find preceeding jobs
      args = (contextJobId,)
      self.execute('SELECT PreceedingJobID FROM Jobs WHERE JobId = ?',args)
      rv = self.fetchAll2Py(UUIDTYPE)
      #print 'CDbApi.getFileByContextPreceedingJobID',contextJobId,rv
      if len(rv)==1 and rv[0] is not None:
        return self.getFileByJobContext(contextJobId=rv[0],fileType=fileType,projectId=projectId,subType=subType,contentFlag=contentFlag)

      return []


    def getJobsByStatus(self,status=JOB_STATUS_QUEUED,projectId=None):
      # Wicked use of TreeLeft to hold remoteMode
      com = 'SELECT JobID FROM Jobs WHERE status = ?'
      args = [str(status)]
      if projectId is not None:
        com = com + ' AND ProjectId=?'
        args.append(projectId)
      com = com + ' ORDER BY CreationTime'
      self.execute(com,args)
      return self.fetchAll2Py(UUIDTYPE)

    '''
    def getJobsByStatus(self,status=JOB_STATUS_QUEUED,projectId=None):
      if projectId is None:
        self.execute('SELECT JobID FROM Jobs WHERE status = ? ORDER BY CreationTime',(status,))
        return self.fetchAll2Py(UUIDTYPE)
      else:
        self.execute('SELECT JobID FROM Jobs WHERE status = ? AND ProjectId=? ORDER BY CreationTime',(status,projectId))
        return self.fetchAll2Py(UUIDTYPE)
    '''

    def createFileType(self,fileTypeName,fileTypeDescription=None):
      try:
        self.getFileTypeId(fileTypeName)
      except:
        pass
      else:
        # Already id for this fileTypeName - bad
        raise CException(self.__class__,165,fileTypeName)

      uniqueID = self.integerId(table='FileTypes',identifier='FileTypeID')

      self.execute("INSERT INTO FileTypes (FileTypeID,FileTypeName,fileTypeDescription) VALUES( ?,?,?)",
                   (uniqueID,fileTypeName,fileTypeDescription ) )
      self.commit()

    def getFileTypeId(self,fileTypeName=None,mimeType=None):
      if fileTypeName is not None:
        u = (fileTypeName,)
        self.execute('SELECT FileTypeID FROM FileTypes WHERE FileTypeName = ?',u)
        rv = self.fetchAll2Py(int)
        if len(rv)==0:
          raise CException(self.__class__,164,fileTypeName)
        else:
          return rv[0]
      elif mimeType is not None:
        if mimeType in FILETYPES_TEXT:
          return FILETYPES_TEXT.index(mimeType)
        else:
          return 0
      else:
        return 0

    def getParamsContainer(self,jobId):
      from core import CCP4Container,CCP4Modules
      paramsFile = self._makeJobFileName(jobId=jobId,mode='PARAMS')
      #print 'CDbApi.getParamsContainer',paramsFile
      if not os.path.exists(paramsFile): raise CException(self.__class__,170,paramsFile)
      if self.parent is not None:
#SJM - I had to make the reference to parent a function call. Otherwise deep-down things try to __init__ with a builtin as
#the parent. (Not sure where this new parent stuff is from anyway).
          if isinstance(self.parent, Callable):
              container = CCP4Container.CContainer(parent=self.parent())
          else:
              container = CCP4Container.CContainer(parent=self.parent)
      else:
          container = CCP4Container.CContainer(parent=CCP4Modules.QTAPPLICATION())
      container.loadDataFromXml(paramsFile)
      return container


    def deleteFileUses(self,jobId=None):
      self.execute("DELETE FROM FileUses WHERE JobId=? AND RoleId=?",(jobId,FILE_ROLE_IN))
      self.commit()


    def gleanJobFiles(self,jobId=None,container=None,projectId=None,dbOutputData=None,roleList=[FILE_ROLE_IN,FILE_ROLE_OUT],unSetMissingFiles=False):

      errorReport = CErrorReport()
      from core import CCP4Data,CCP4File
      if container is None:
        container = self.getParamsContainer(jobId=jobId)

      if projectId is None:
        projectId = self.getJobInfo(jobId=jobId,mode='projectId')
      projectName = self.getProjectInfo(projectId=projectId,mode='projectname')
      #print 'Database loading file data for job',jobId,projectId,projectName,dbOutputData,roleList
      preceedingJobs = []

      #if container.header.pluginTitle.isSet():
        #print 'gleanJobFiles update jobtitle',container.header.pluginTitle.__str__()
      #  jInfo = self.getJobInfo(jobId,'jobtitle')
      #  if jInfo is not None and jInfo != "":
      #    print 'CDbApi.gleanJobFiles jobtitle already set to',jInfo,'not updating'
      #  else:
      #    self.updateJob(jobId=jobId,key='jobtitle',value=container.header.pluginTitle.__str__())
      if container.header.pluginVersion.isSet():
        self.updateJob(jobId=jobId,key='taskversion',value=container.header.pluginVersion.__str__())


      for file_role in roleList:
        subcontainer = (container.inputData,container.outputData)[[FILE_ROLE_IN,FILE_ROLE_OUT].index(file_role)]
        if file_role == FILE_ROLE_OUT and dbOutputData is not None:
          keyList = dbOutputData
        else:
          keyList = subcontainer.dataOrder()
        for key in keyList:
          obj0 = subcontainer.__getattr__(key)
          try:
            objList,xmlText,keyValues = obj0.saveToDb()
            jobParamName = obj0.objectName()
          except:
            print('ERROR in gleanJobFiles for',key)
            objList,xmlText,keyValues = [],None,{}
            jobParamName = ''
          #print 'gleanJobFiles saveToDb',key,objList,xmlText,keyValues,jobParamName
          # Create list of existing files to import - this converts CDataFile and CList to same form for
          # procesing
          importList = []
          idx = -1
          for obj in objList:
            #print 'gleanJobFiles obj',obj,type(obj),isinstance(obj,CCP4File.CDataFile), obj.isSet()
            if isinstance(obj,CCP4File.CDataFile) and obj.isSet():
              idx += 1
              if not obj.exists():
                errorReport.append(self.__class__,173,str(obj.objectName())+' '+str(obj),stack=False)
                if unSetMissingFiles: obj.unSet()
              else:
                if len(objList)>1:
                  importList.append([obj,jobParamName+'.'+str(idx)])
                else:
                  importList.append([obj,jobParamName])
            elif isinstance(obj,CCP4Data.CList) and isinstance (obj.subItemObject(),CCP4File.CDataFile):
              idx=-1
              for fileObj in obj:
                idx += 1
                if fileObj.exists():
                  importList.append([fileObj,jobParamName+'.'+str(idx)])
                elif fileObj.isSet():
                  errorReport.append(self.__class__,173,str(fileObj.objectName())+' '+str(fileObj),stack=False)
                  if unSetMissingFiles: fileObj.unSet()
          #print 'gleanJobFiles importList',self.getJobInfo(jobId=jobId,mode='jobnumber'),importList
          for obj,jobParamName in importList:
              fileType = obj.qualifiers('mimeTypeName')
              fileId,fromJobId = self.matchFileName(fileObject=obj)
              if fileId is not None:
                  print('Recording file in database matches existing',key,str(obj),fileId,fromJobId)
              if fileId is None:
                if file_role == FILE_ROLE_OUT:
                  try:
                    #print 'gleanJobFiles output file',obj.objectName()
                    fileId = self.createFile(jobId=jobId,fileObject=obj,projectId=projectId,jobParamName=jobParamName)
                    obj.dbFileId = fileId
                    #print 'gleanJobFile setting dbFileId',obj.objectName(),obj,fileId
                  except CException as e:
                    errorReport.extend(e,stack=False)
                  except Exception as e:
                    errorReport.append(self.__class__,171,str(e),stack=False)
                else:
                  print('ERROR in gleanJobFiles - no fileId for an input file',obj)
                  try:
                    fileId = self.createFile(jobId=None,fileObject=obj,projectId=projectId,jobParamName=jobParamName)
                    obj.dbFileId = fileId
                    #print 'gleanJobFile setting dbFileId',obj.objectName(),obj,fileId
                  except CException as e:
                    errorReport.extend(e,stack=False)
                  except Exception as e:
                    errorReport.append(self.__class__,171,str(e),stack=False)
                  try:
                    self.createFileUse(jobId=jobId,fileId=fileId,role=file_role,jobParamName=jobParamName)
                  except CException as e:
                    errorReport.extend(e,stack=False)
                  except Exception as e:
                    errorReport.append(self.__class__,172,str(e),stack=False)
              else:
                if file_role == FILE_ROLE_OUT:
                  #print 'CDbApi.gleanJobFiles calling createFileUse'
                  # This could be a parent job that is 'passing on' the output from a sub-job
                  # The pipeline script ought to have copied the file to the parent jobDir
                  # and we should not be here! But we'll record asa file use
                  #self.updateFile(fileId=fileId,key='jobId',value=jobId)
                  if not obj.dbFileId.isSet(): obj.dbFileId = fileId
                  #print 'gleanJobFile setting dbFileId',obj.objectName(),obj,fileId
                  try:
                    self.createFileUse(jobId=jobId,fileId=fileId,role=file_role,jobParamName=jobParamName)
                    #print 'Done createFileUse'
                  except CException as e:
                    errorReport.extend(e,stack=False)
                  except Exception as e:
                    errorReport.append(self.__class__,172,str(e),stack=False)
                  '''
                  try:
                    self.createFileUse(jobId=jobId,fileId=fileId,role=file_role)
                  except CException as e:
                    errorReport.extend(e)
                  except Exception as e:
                    errorReport.append(self.__class__,172,str(e))
                  '''

                else:
                  # File already in in db so add a FileUse record
                  if fromJobId is not None:
                    if not (fileType,fileId,fromJobId) in preceedingJobs: preceedingJobs.append([fileType,fileId,fromJobId])
                  try:
                    self.createFileUse(jobId=jobId,fileId=fileId,role=file_role,jobParamName=jobParamName)
                    if not obj.dbFileId.isSet():  obj.dbFileId = fileId
                    print('Recording file in database',jobParamName,obj,fileId)
                  except CException as e:
                    errorReport.extend(e,stack=False)
                  except Exception as e:
                    errorReport.append(self.__class__,172,str(e),stack=False)
          if file_role == FILE_ROLE_OUT and xmlText is not None:
            # Not a file but is output data - save as xdata
            #xmlText =obj.xmlText(pretty_print=False)
            dataClass = obj0.__class__.__name__
            #print 'gleanJobFiles',dataClass,xmlText
            self.createXData(jobId=jobId,dataClass=dataClass,dataXml=xmlText)
          if file_role == FILE_ROLE_OUT and len(keyValues)>0:
            print('Saving key:value pairs',keyValues)
            for keyTypeName,value in list(keyValues.items()):
              self.createJobKeyValue(jobId=jobId,keyTypeName=keyTypeName,value=value)
            print('Done saving key:value pairs')


      #print 'gleanJobFiles preceedingJobs',preceedingJobs
      for fileType,fileId,fromJobId in preceedingJobs:
        if fromJobId != jobId and fileType == 'chemical/x-pdb':
          try:
            self.updateJob(jobId,'preceedingjobid',fromJobId)
          except CException as e:
            errorReport.extend(e)
      return errorReport


    def _makeJobFileName(self,jobId=None,mode='PARAMS'):

      '''
      exts = { 'PARAMS' : '.params.xml',
               'JOB_INPUT' : '.input_params.xml',
               'PROGRAMXML' : '.program.xml',
               'LOG' : '.log',
               'STDOUT' : '.stdout.txt',
               'STDERR' : '.stderr.txt' }
      '''
      defFiles = { 'PARAMS' : 'params.xml',
               'JOB_INPUT' : 'input_params.xml',
               'PROGRAMXML' : 'program.xml',
               'LOG' : 'log.txt',
               'COM' : 'com.txt',
               'STDOUT' : 'stdout.txt',
               'STDERR' : 'stderr.txt' }


      import os
      path = os.path.join(self.jobDirectory(jobId=jobId),defFiles.get(mode,'unknown.unk'))
      #print 'CDbApi._makeJobFileName',mode,path
      return path

    def createImportFile(self,fileId=None,sourceFileName=None,sourceFileId=None,jobId=None,annotation=None,
                         checksum=None,importNumber=None,reference=None):
      # Record the 'outside world' name of file OR the fileId in another project
      importId = self.uniqueId(table='ImportFiles',identifier='ImportID')
      if sourceFileId is not None:
        try:
          jobId = self.getFileInfo(fileId=sourceFileId,mode='jobid')
        except:
          raise CException(self.__class__,251,str(sourceFileId))
        sourceFileName = None
      elif sourceFileName is None:
        raise CException(self.__class__,250)

      time = self.currentTime()
      try:
        lastModifiedTime = os.path.getmtime(sourceFileName)
      except:
        lastModifiedTime = None

      if checksum is None and sourceFileName is not None and os.path.exists(sourceFileName):
        try:
          blockSize = 256*128
          import hashlib
          md5 = hashlib.md5()
          with open(sourceFileName,'rb') as f:
            for chunk in iter(lambda: f.read(blockSize), b''):
              md5.update(chunk)
          checksum = md5.hexdigest()
          print('CDbApi.createImportFile calculated checksum',sourceFileName,checksum)
        except:
          print('Failed attempting to calculate checksum',sourceFileName)

      args = (importId,fileId,time,lastModifiedTime,checksum,sourceFileName,sourceFileId,annotation,importNumber,reference)
      #print 'CDbApi.createImportFile',args
      self.execute("INSERT INTO ImportFiles (ImportID,FileId,CreationTime,lastModifiedTime,checksum,SourceFilename,SourceFileId,Annotation,ImportNumber,Reference) VALUES (?,?,?,?,?,?,?,?,?,?)",args)
      self.commit()

      self.importFileCreated.emit({'jobId':jobId,'importId':importId})

      return importId

    def updateImportFile(self,importId=None,key=None,value=None,jobId=None):
      if key is None: raise CException(self.__class__,252)
      key = key.lower()
      if not key in ['annotation','exportfileid','creationtime']:
        raise CException(self.__class__,252,key)

      #self.setDiagnostic(True)
      self.execute("UPDATE ImportFiles SET "+key+" = ? WHERE ImportID = ?",(value,importId,))
      #self.setDiagnostic(False)
      self.commit()

      if jobId is None:
        self.execute('SELECT Files.JobId FROM Files INNER JOIN ImportFiles ON Files.FileId = ImportFiles.FileId WHERE ImportFiles.ImportId=?',(importId,))
        jobId = self.fetchAll2Py(UUIDTYPE)[0]

      self.importFileUpdated.emit({'jobId':jobId,'importId':importId,'key':key,'value':value})

    def deleteImportFile(self,importId=None):
      self.execute("DELETE FROM ImportFiles WHERE ImportID = ?",(importId,))
      self.importFileDeleted.emit({'importId':importId})

    def getImportFileInfo(self,importId=None,mode='all',fileId=None):
      itemList = []
      fileList = []
      if mode == 'all':
        itemList.extend(self.IMPORTFILEITEMS)
      elif (not isinstance(mode,list)):
        if mode.lower() in self.IMPORTFILEITEMS:
          itemList = [mode.lower()]
        elif mode.lower() in  ['fileid','jobid']:
          fileList = [mode.lower()]
      else:
        for item in mode:
          if item.lower() in  self.IMPORTFILEITEMS: itemList.append(item.lower())
          if item.lower() in ['fileid','jobid'] : fileList.append(item.lower())
      if len(itemList)<0 and len(fileList)<0:
        raise CException(self.__class__,270,mode)

      selcom = 'SELECT '
      typeConv = []
      ret = {}
      if len(itemList)>0:
        for item in itemList:
          selcom = selcom + item +', '
          typeConv.append(self.IMPORTFILETYPES[self.IMPORTFILEITEMS.index(item)])
        if importId is not None:
          selcom = selcom[0:-2] +  " FROM ImportFiles WHERE ImportID=?"
          self.execute(selcom,(importId,))
        elif fileId is not None:
          selcom = selcom[0:-2] +  " FROM ImportFiles WHERE FileID=?"
          self.execute(selcom,(fileId,))
        rv = self.fetchAll2PyList(typeConv)
        #print 'getProjectJobListInfo rv',rv
        if len(rv) != 1: raise CException(self.__class__,271)
        i = 0
        for item in itemList:
          ret[item] = rv[0][i]
          i = i + 1
      if len(fileList)>0:
        selcom = 'SELECT '
        typeConv = []
        for item in fileList:
          selcom = selcom + item +', '
          typeConv.append(UUIDTYPE)
        selcom = selcom[0:-2] +  " FROM Files WHERE ImportID=?"
        self.execute(selcom,(importId,))
        rv1 = self.fetchAll2PyList(typeConv)
        i = 0
        for item in itemList:
          ret[item] = rv1[0][i]
          i = i + 1
      return ret

    def getImportedFile(self,sourceFileName,projectId=None,lastModifiedTime=None,mimeTypeName=None,fileContent=None,reference=None,fileType=None):
      # Beware different mini-MTZ types may have been imported from the same source file
      # The importAnnotation is the list of input columns and should match if exactly the same data
      # was imported from a monster-MTZ

      cmd = 'SELECT ImportFiles.ImportID,Files.FileId,ImportFiles.checksum,Files.Annotation, Jobs.JobId, Jobs.JobNumber FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId =Files.FileId INNER JOIN Jobs ON Files.JobId = Jobs.JobId WHERE ImportFiles.SourceFilename=?'
      args = [sourceFileName]
      if projectId is not None:
        cmd = cmd + ' AND Jobs.ProjectId = ?'
        args.append(projectId)
      if lastModifiedTime is not None:
        cmd = cmd + ' AND ImportFiles.LastModifiedTime BETWEEN ? AND ?'
        args.extend(lastModifiedTime-0.1,lastModifiedTime+0.1)
      if mimeTypeName  is not None and mimeTypeName in FILETYPES_TEXT:
        fileType = FILETYPES_TEXT.index(mimeTypeName)
        cmd = cmd + ' AND Files.FileTypeId = ?'
        args.append(fileType)
      if fileContent is not None:
         cmd = cmd + ' AND Files.FileContent = ?'
         args.append(fileContent)
      if reference is not None:
        cmd = cmd + ' AND ImportFiles.reference = ?'
        args.append(reference)
      if fileType is not None:
        cmd = cmd + ' AND Files.FiletypeID = ?'
        args.append(fileType)
      #print 'CDbApi.getImportedFile',cmd,args
      self.execute(cmd,args)
      rv = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,str,str,str,str])

      #print 'getImportedFile',sourceFileName,rv
      return rv

    def getImportFileInstances(self,jobId=None,brief=True):
      if brief:
        self.execute('SELECT ImportFiles.ImportID, ImportFiles.CreationTime, ImportFiles.LastModifiedTime, ImportFiles.SourceFilename, Files.FileId, ImportFiles.Annotation FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId =Files.FileId WHERE Files.JobId = ?',(jobId,))
        rv = self.fetchAll2PyList([UUIDTYPE,float,float,str,UUIDTYPE,str])
        retList = []
        for importid,time,lastmodifiedtime,sourcefilename,fileid,annotation in rv:
          retList.append( { 'importid' : importid, 'sourcefilename': sourcefilename, 'creationtime' : time, 'lastmodifiedtime': lastmodifiedtime, 'fileid':fileid,'annotation':annotation } )
        return retList
      else:
        # Return complete importfile and file info
        com = 'SELECT'
        for item in self.IMPORTFILEITEMS:
          com = com + ' ImportFiles.'+item+','
        for item in self.FILEITEMS:
          com = com + ' Files.'+item+','
        com = com[0:-1] + ' FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId =Files.FileId WHERE Files.JobId = ?'
        #print 'getImportFileInstances',com
        self.execute(com, (jobId,))
        rv = self.fetchAll2PyList(self.IMPORTFILETYPES+self.FILETYPES)
        #print 'getImportFileInstances rv',rv
        infoList = []
        for f in rv:
          infoList.append({})
          for i in range(len(self.IMPORTFILEITEMS)):
            if self.IMPORTFILEITEMS[i] not in ['fileid','exportfileid','lastmodifiedtime']:
              infoList[-1][self.IMPORTFILEITEMS[i]] = f[i]
          for i in range(len(self.FILEITEMS)):
            if self.FILEITEMS[i] not in ['jobid']:
              infoList[-1][self.FILEITEMS[i]] = f[i+len(self.IMPORTFILEITEMS)]
          infoList[-1]['filetypename'] = FILETYPES_TEXT[infoList[-1]['filetypeid']][1]

        return infoList


    def createExportFile(self,fileId=None,exportFilename=None):
      # Record the 'outside world' name of exported file
      # Ensure fileId exists
      if fileId is None: raise CException(self.__class__,253,str(fileId))
      if exportFilename is None or not os.path.exists(exportFilename):
        raise CException(self.__class__,254,str(exportFilename))
      try:
        jobId = self.getFileInfo(fileId=fileId,mode='jobid')
      except:
        raise CException(self.__class__,253,str(fileId))

      exportId = self.uniqueId(table='ExportFiles',identifier='ExportID')
      time = self.currentTime()

      args = (exportId,time,exportFilename,fileId)
      #print 'CDbApi.createEportFile',args
      self.execute("INSERT INTO ExportFiles (ExportID,CreationTime,ExportFilename,FileId) VALUES (?,?,?,?)",args)
      self.commit()

      self.exportFileCreated.emit({'jobId':jobId,'exportId':exportId})

      return exportId

    def getExportFileInstances(self,fileId=None,jobId=None):
      #print 'getExportFileInstances',fileId,type(fileId)
      if fileId is not None:
        self.execute('SELECT ExportID, CreationTime, ExportFilename, FileId FROM ExportFiles WHERE FileId = ?',(fileId,))
      elif jobId is not None:
        self.execute('SELECT ExportFiles.ExportID, ExportFiles.CreationTime, ExportFiles.ExportFilename, ExportFiles.FileId FROM ExportFiles INNER JOIN Files ON ExportFiles.FileId =Files.FileId WHERE Files.JobId = ?',(jobId,))
      else:
        return []
      rv = self.fetchAll2PyList([UUIDTYPE,float,str,UUIDTYPE])
      retList = []
      for exportid,time,exportfilename,fileid in rv:
        retList.append( { 'exportid' : exportid, 'exportfilename': exportfilename, 'creationtime' : time, 'fileid':fileid } )
      return retList

    def getExportFileInfo(self,exportId=None,mode='all'):
      self.execute('SELECT ExportFiles.ExportID, ExportFiles.CreationTime, ExportFiles.ExportFilename, Files.JobId, Jobs.JobNumber, Jobs.TaskName FROM ExportFiles INNER JOIN Files ON ExportFiles.FileId=Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE ExportFiles.ExportID = ?',(exportId,))
      rv = self.fetchAll2PyList([UUIDTYPE,float,str,UUIDTYPE,str,str])
      retList = []
      for exportid,creationtime,exportfilename,jobid,jobnumber,taskname in rv:
          retList.append( { 'exportid' : exportid, 'exportfilename': exportfilename, 'creationtime' : creationtime, 'jobid': jobid, 'jobnumber' : jobnumber, 'taskname': taskname } )
      if len(retList)==0:
        return {}
      else:
        return retList[0]



    def getExportFilesByFileType(self,fileTypeClass=None,before=None,projectId=None):
      com = 'SELECT ExportFiles.ExportID, ExportFiles.CreationTime, ExportFiles.ExportFilename, Files.JobId, Jobs.JobNumber, Jobs.TaskName FROM ExportFiles INNER JOIN Files ON ExportFiles.FileId=Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId'
      args = []
      where = ' WHERE '
      if fileTypeClass in FILETYPES_CLASS:
        args.append(FILETYPES_CLASS.index(fileTypeClass))
        com = com + where +'Files.FileType = ?'
        where = ' AND '
      if before is not None:
        com = com + where +'ExportFiles.CreationTime < ?'
        args.append(before)
        where = ' AND '
      if projectId is not None:
        com = com + where +'Jobs.ProjectId=?'
        args.append(projectId)
      com = com + ' ORDER BY Jobs.CreationTime DESC'
      self.execute(com,args)
      rv = self.fetchAll2PyList([UUIDTYPE,float,str,UUIDTYPE,str,str])
      retList = []
      for exportid,creationtime,exportfilename,jobid,jobnumber,taskname in rv:
        retList.append( { 'exportid' : exportid, 'exportfilename': exportfilename, 'creationtime' : time, 'jobid': jobid, 'jobnumber' : jobnumber, 'taskname': taskname } )
      return retList

    def createFile(self,jobId=None,projectId=None,fileObject=None,fileTypeId=None,fileTypeName=None,
             sourceFileName=None,sourceFileId=None,sourceFileAnnotation=None,importNumber=None,reference=None,
                   **kw):

      # fileObject is a CDataFile
      # File is output by jobId or
      # imported from sourceFileName or sourceFileId (file from another project)
      #print 'createFile',jobId,projectId,fileObject,fileObject.objectName(),sourceFileName,kw
      if projectId is not None:
        permission = self.getProjectPermission(projectId=projectId)
      else:
        projectId,permission = self.getJobPermission(jobId=jobId)
      if permission < PRIVILEGE_WRITE: raise CException(self.__class__,112)

      # The file should exist! - but gleanJobFiles has checked this
      if self._preferences['testFileExists']:
        if fileObject is not None and not fileObject.exists():
          raise CException(self.__class__,162,fileObject.fullPath.__str__())
        # And the path should be correct for the job (if it is not imported file)
        #print 'createFile jobdir',fileObject,os.path.split(fileObject.fullPath.__str__())[0],self.jobDirectory(jobId=jobId)

        if sys.platform[0:3] != 'win' and fileObject is not None and sourceFileName is None and sourceFileId is None:
          from core import CCP4Utils
          try:
            isSame = CCP4Utils.samefile ( os.path.normpath(os.path.split(fileObject.fullPath.__str__())[0]),
                                                 self.jobDirectory(jobId=jobId) )
          except:
            raise CException(self.__class__,128,fileObject.fullPath.__str__()+' '+self.jobDirectory(jobId=jobId),stack=False)
          if not isSame:
            raise CException(self.__class__,127,fileObject.fullPath.__str__()+' expected in '+str(self.jobDirectory(jobId=jobId)),stack=False)
        else:
          pass

      annotation = None
      if fileObject is not None:
        baseName = fileObject.baseName.__str__()
        relPath = fileObject.relPath.__str__()
        if relPath.count('CCP4_IMPORTED_FILES'):
          filePathFlag = PATH_FLAG_IMPORT_DIR
        else:
          filePathFlag = PATH_FLAG_JOB_DIR
        if fileObject.annotation.isSet():
          annotation = str(fileObject.annotation)
        if fileObject.contentFlag.isSet():
          fileContent = int(fileObject.contentFlag)
        else:
          fileContent = None
        jobParamName = fileObject.objectName()
        if jobParamName is None or len(jobParamName)==0: jobParamName = kw.get('jobParamName',None)
        if fileObject.subType.isSet():
          subType = int(fileObject.subType)
        else:
          subType = None
      else:
        baseName = kw.get('baseName',None)
        relPath = kw.get('relPath',None)
        filePathFlag = kw.get('pathFlag',None)
        annotation =  kw.get('annotation',None)
        fileContent = kw.get('fileContent',None)
        jobParamName = kw.get('jobParamName',None)
        subType =  kw.get('subType',None)

      if fileTypeId is None:
        if fileTypeName is None: fileTypeName = fileObject.qualifiers('mimeTypeName')
        fileTypeId = self.getFileTypeId(mimeType=fileTypeName)
        #print 'CDbApi.createFile',fileObject.objectName(),fileObject.qualifiers('mimeTypeName'),fileTypeId

      fileId = self.uniqueId(table='Files',identifier='FileID')

      args = [fileId,jobId,jobParamName,baseName,fileTypeId,annotation,fileContent,subType,filePathFlag]
      #print 'CDbApi.createFile args',args
      self.execute("INSERT INTO Files (FileID,JobID,JobParamName,Filename,FiletypeID,Annotation,FileContent,FileSubType,PathFlag) VALUES (?,?,?,?,?,?,?,?,?)",args)

      if sourceFileName is not None or sourceFileId is not None:
        importId = self.createImportFile(fileId=fileId,sourceFileName=sourceFileName,sourceFileId=sourceFileId,jobId=jobId,
                                  annotation=sourceFileAnnotation,importNumber=importNumber,reference=reference)


      self.commit()

      self.fileCreated.emit({'fileId':fileId,'jobId':jobId,'projectId':projectId})
      return fileId


    # Why would you edit/delete a File or FileUse record other than deleting a job???
    # Answer: cos someones deleted the file

    def createFileUse(self,jobId=None,fileId=None,role=FILE_ROLE_IN,jobParamName=None):
      if jobId is None:
        raise CException(self.__class__,114,'jobId')
      if fileId is None:
        raise CException(self.__class__,114,'fileId')

      projectId,permission = self.getJobPermission(jobId=jobId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      args = [fileId,jobId,role,jobParamName]
      #print 'CDbApi.createFileUse',args
      self.execute("INSERT INTO FileUses (FileID,JobID,RoleID,JobParamName) VALUES (?,?,?,?)",args)
      self.commit()
      self.fileUseCreated.emit({'fileId':fileId,'jobId':jobId,'projectId':projectId,'role':role})

    def matchFileName(self,fileName=None,fileObject=None):
      if fileObject is not None:
        if fileObject.dbFileId.isSet():
          try:
            jobId = self.getFileInfo(fileId=fileObject.dbFileId.get(),mode='jobid')
          except:
            pass
          else:
            #print 'matchFileName from dbFileId',fileObject.dbFileId.get(),jobId
            return fileObject.dbFileId.get(),jobId
        baseName = str(fileObject.baseName)
        relPath = str(fileObject.relPath)
        fileProjectId = None
        if fileObject.project.isSet():
          try:
            fileProjectId = self.getProjectId(projectName=str(fileObject.project))
          except CException as e:
            # Not recognised as project - is it an alias?
            try:
              aliasPath = self.getAliasDirectory(str(fileObject.project))
            except:
              return None,None
            relPath = os.path.join(aliasPath,relPath)

      elif fileName is not None:
        path,baseName = os.path.split(os.path.abspath(fileName))
        fileProjectId,fileProjectName,fileProjectPath  = self.matchProjectDirectory(path)
        #if projectId is None: return None
        #print 'CDbApi.matchFileName',path,fileProjectId
        if fileProjectId is None:
          relPath = path
        else:
          relPath = fileProjectPath
        if len(relPath)==0: relPath = None
      else:
        return None,None

      if fileProjectId is not None:
        jobNumber = ''
        testPath = relPath
        while len(testPath)>0:
          testPath,term = os.path.split(testPath)
          if term[0:4] == 'job_':
            jobNumber =  term[4:] + '.' + jobNumber
          elif term == 'CCP4_JOBS':
            pass
          else:
            jobNumber = ''
            testPath = ''
        #print 'matchFileName',relPath,'jobNumber',jobNumber
        if len(jobNumber)>0:
          jobNumber = jobNumber[0:-1]
          args = (baseName,jobNumber,fileProjectId,)
          self.execute('SELECT  Files.FileId,Files.JobId FROM Files INNER JOIN Jobs ON Files.JobId = Jobs.JobId WHERE Files.FileName= ? AND Jobs.JobNumber = ? AND Jobs.ProjectId = ?',args)
          fileIdList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE])
          if len(fileIdList)==1:
            return fileIdList[0][0],fileIdList[0][1]
        else:
          # We don't have a standard relPath - keh?
          args = (baseName,fileProjectId)
          self.execute('SELECT Files.FileId,Files.JobId FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Files.Filename = ?  AND Jobs.ProjectID = ?',args)
          fileIdList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE])
          #print 'CDbApi.matchFileName fileIdList',args,fileIdList
          if len(fileIdList) == 1:
            return fileIdList[0][0],fileIdList[0][1]

      return None,None

    def getFileInfo(self,fileId=None,mode='all'):

      fileId = UUIDTYPE(fileId)
      if mode == 'all':
          itemList = []
          itemList.extend(self.FILEITEMS)
      elif not isinstance(mode,list):
          itemList = [mode.lower()]
      else:
        itemList = []
        for item in mode:
          itemList.append(item.lower())

      selcom = "SELECT "
      typeConv = []
      retList = []
      jobItems = []
      for jobItem in ['jobnumber','taskname','projectname','projectid']:
        if jobItem in itemList:
          if not 'jobid' in itemList: itemList.append('jobid')
          jobItems.append(jobItem)
      if 'sourcefilename' in itemList:
        if not 'importid' in itemList: itemList.append('importid')
      if 'relpath' in itemList:
        if not 'jobnumber' in jobItems: jobItems.append('jobnumber')
        if not 'jobid' in itemList: itemList.append('jobid')
        if not 'importid' in itemList: itemList.append('importid')
        if not 'filetypeid' in itemList: itemList.append('filetypeid')
        if not 'pathflag' in itemList: itemList.append('pathflag')
        ifRelPath = True
      else:
        ifRelPath = False
      for item in itemList:
        try:
          if item in ['filetypeid','filetype','fileclass']:
            if selcom.count('filetypeid') == 0:
              selcom = selcom + 'Files.filetypeid, '
              typeConv.append(int)
              retList.append('filetypeid')
          elif item in ['importid','sourcefilename']:
            selcom = selcom + 'ImportFiles.'+item +', '
            typeConv.append(UUIDTYPE)
            retList.append(item)
          else:
            ii = self.FILEITEMS.index(item)
            if ii>0:
              selcom = selcom + 'Files.'+item +', '
              typeConv.append(self.FILETYPES[ii])
              retList.append(item)
        except:
          pass
      if 'importid' in retList:
        selcom = selcom[0:-2] +  " FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE Files.FileID=?"
      else:
        selcom = selcom[0:-2] +  " FROM Files WHERE Files.FileID=?"
      #print 'getFileInfo selcom',selcom,jobItems
      self.execute(selcom,(fileId,))
      rv = self.fetchAll2PyList(typeConv)
      if len(rv) != 1:
        #import traceback
        #traceback.print_stack(limit=7)
        raise CException(self.__class__,180,str(fileId))
      else:
        ret = self.fileInfoDict(retList,rv[0])
        #print 'getFileInfo after fileInfoDict',ret
        # Need to append the 'derived' file type info?
        if 'fileid' in itemList: ret['fileid'] = fileId
        if 'filetype' in itemList:
          ret['filetype'] = FILETYPES_TEXT[ret['filetypeid']]
        if 'fileclass' in itemList:
          ret['fileclass'] = FILETYPES_CLASS[ret['filetypeid']]
        if len(jobItems)>0:
          if ret['jobid'] is None:
            jobInfo = { 'jobnumber' : None, 'taskname':None,'projectid':None,'projectname': None, 'relpath' : None}
          else:
            jobInfo = self.getJobInfo(jobId=ret['jobid'],mode=jobItems,returnType=dict)
            if isinstance(jobInfo,dict):
              ret.update(jobInfo)
            else:
              ret[jobItems[0]] = jobInfo

            if ifRelPath:
              ret['relpath'] = self.getRelPath(ret['pathflag'],ret['filetypeid'],jobInfo['jobnumber'],ret['importid'])
        #print 'getFileInfo with derived',jobItems,ret
        if len(itemList)==1:
          # Just return the value if only one item requested
          return ret[itemList[0]]
        else:
          # Return multiple items as a dictionary
          #print 'getFileInfo ret',ret
          return ret


    def fileInfoDict(self,itemList=[],fileValue=[]):
      #print 'fileInfoDict',itemList,fileValue
      ii = -1
      ret = {}
      for item in itemList:
        ii = ii + 1
        if item == 'filetype':
          if fileValue[ii]>0 and fileValue[ii]<len(FILETYPES_TEXT):
            ret[item]= FILETYPES_TEXT[fileValue[ii]]
          else:
            # FILETYPES_TEXT[0] is 'Unknown'
            ret[item] = FILETYPES_TEXT[0]
        elif item == 'fileclass':
          if fileValue[ii]>0 and fileValue[ii]<len(FILETYPES_CLASS):
            ret[item]= FILETYPES_CLASS[fileValue[ii]]
          else:
            # FILETYPES_TEXT[0] is 'Unknown'
            ret[item] = FILETYPES_CLASS[0]
        else:
          ret[item] = fileValue[ii]
      #print 'jobInfoDict',ret
      return ret


    def getFullPath(self,fileId=None,projectDirectory=None):
      self.execute('SELECT Files.JobId,ImportFiles.ImportId,Files.Filename, Files.FileTypeId, Files.PathFlag FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE Files.FileId = ?', (fileId,))
      fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,str,int,int])
      #print 'CDbApi.getFullPath',fileId,fileList
      if len(fileList) != 1:  raise CException(self.__class__,168)
      jobId,importId,fileName,fileTypeId,pathFlag = fileList[0]

      if projectDirectory is None or importId is None:
        if jobId is None:
          raise CException(self.__class__,168)
        self.execute('SELECT JobNumber,ProjectId FROM Jobs WHERE JobId = ?',(jobId,))
        jobList = self.fetchAll2PyList([str,UUIDTYPE])
        if len(jobList)!=1:  raise CException(self.__class__,168)
        jobNumber,projectId = jobList[0]
        if projectDirectory is None: projectDirectory = self.getProjectDirectory(projectId)

      relPath=self.getRelPath(pathFlag,fileTypeId,jobNumber,importId)
      path = os.path.join(projectDirectory,relPath,fileName)
      print('CDbApi.getFullPath',path)
      return path

    def getRelPath(self,pathFlag,fileTypeId,jobNumber,importId=None):
      if pathFlag is not None:
        if pathFlag == PATH_FLAG_IMPORT_DIR:
          relPath= 'CCP4_IMPORTED_FILES'
        else:
          relPath=self.jobRelPath(jobNumber=jobNumber)
      else:
        # pathFlag not set in old data
        if importId is not None and fileTypeId not in MINIMTZFILETYPES:
          relPath= 'CCP4_IMPORTED_FILES'
        else:
          relPath=self.jobRelPath(jobNumber=jobNumber)
      return relPath


    def getJobImportFiles(self,jobId=None):
      args = (jobId,)
      self.execute('SELECT Files.JobId,Files.FileID,ImportFiles.ImportId,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files INNER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE Files.JobId = ?',args)
      fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,UUIDTYPE,int,str,str])
      return fileList


    def getJobFiles(self,jobId=None,role=FILE_ROLE_OUT,mode='fileId',searchFileUses=True,fileTypes=None):

      outputList = []
      diry = self.getProjectDirectory(jobId=jobId)
      relPath =None
      if fileTypes is not None:
        fileTypeSele = ' AND Files.FileTypeID IN ('+str(fileTypes[0])
        for item in fileTypes[1:]:
          fileTypeSele = fileTypeSele + ',' + str(item)
        fileTypeSele += ')'
      else:
        fileTypeSele = ''

      if role == FILE_ROLE_OUT:
        # Search the Files table for output files
        args = [jobId,]
        seleList = ' LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE ImportFiles.ImportID is NULL AND JobID = ?'+fileTypeSele

        if mode == 'all':
          self.execute('SELECT Files.JobId,Files.FileID,ImportFiles.ImportId,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files'+seleList,args)
          outputList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,UUIDTYPE,int,str,str])
        elif mode == 'fileId':
          self.execute('SELECT Files.FileID FROM Files'+seleList,args)
          outputList = self.fetchAll2Py(UUIDTYPE)
        elif mode == 'fullPath':
          self.execute('SELECT Files.Filename FROM Files'+seleList,args)
          fileList = self.fetchAll2Py(str)
          #print 'getJobFiles fileList',fileList
          for item in fileList:
            if relPath is None: relPath = self.jobRelPath(jobId=jobId)
            outputList.append(os.path.join(diry,relPath,item))
        elif mode == 'fileName':
          self.execute('SELECT Files.Filename FROM Files'+seleList,args)
          outputList = self.fetchAll2Py(str)
        elif mode == 'annotation':
          self.execute('SELECT Files.Annotation FROM Files'+seleList,args)
          outputList = self.fetchAll2Py(str)
        elif mode in ['fileType','fileTypeId']:
          self.execute('SELECT Files.FileTypeID FROM Files'+seleList,args)
          typesList = self.fetchAll2Py(int)
          #print 'getJobFiles FileTypeID',typesList
          if mode == 'fileTypeId':
            outputList = typesList
          else:
            for item in typesList:
              if item is None:
                outputList.append(FILETYPES_TEXT[0])
              else:
                outputList.append(FILETYPES_TEXT[item])
        else:
          raise CException(self.__class__,167)

      if not searchFileUses: return outputList

      # Search the FileUses table for input or output files
      args = [jobId,role]
      if mode == 'all':
        self.execute('SELECT Files.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobID = ? AND FileUses.RoleID = ?'+fileTypeSele,args)
        fileUseList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,int,str,str])
        for fileList in fileUseList:
          outputList.append([fileList[0],fileList[1],None,fileList[2],fileList[3],fileList[4]])
      elif mode == 'fileId':
        self.execute('SELECT FileID FROM FileUses WHERE JobID = ? AND RoleID = ?',args)
        rv = self.fetchAll2Py(UUIDTYPE)
        #print 'CDbApi.getJobFiles FileUses',args,rv
        outputList.extend(rv)
      elif mode == 'fullPath':
        self.execute('SELECT Files.JobId,Files.Filename FROM Files INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobID = ? AND FileUses.RoleID = ?'+fileTypeSele,args)
        fileList = self.fetchAll2PyList([UUIDTYPE,str])
        for jobId,fileName in fileList:
          if jobId  is not None:
            relPath = self.jobRelPath(jobId=jobId)
            outputList.append(os.path.join(diry,relPath,fileName))
          else:
            outputList.append(os.path.join(diry,'CCP4_IMPORTED_FILES',fileName))
      elif mode == 'fileName':
        self.execute('SELECT Files.Filename FROM Files INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobID = ? AND FileUses.RoleID = ?'+fileTypeSele,args)
        outputList.extend(self.fetchAll2Py(str))
      elif mode == 'annotation':
        self.execute('SELECT Files.Annotation FROM Files INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobID = ? AND FileUses.RoleID = ?'+fileTypeSele,args)
        outputList.extend(self.fetchAll2Py(str))
      elif mode == 'fileType':
        self.execute('SELECT Files.FileTypeID FROM Files INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobID = ? AND FileUses.RoleID = ?'+fileTypeSele,args)
        typesList = self.fetchAll2Py(int)
        #print 'getJobFiles join FileTypeID',typesList
        for item in typesList:
          if item is None:
            outputList.append(FILETYPES_TEXT[0])
          else:
            outputList.append(FILETYPES_TEXT[item])
      return outputList

    def getProjectFiles(self,projectId=None,fileType=None,subType=None,contentFlag=None,topLevelOnly=True):
      args = [projectId]
      com = 'SELECT Files.JobId,Files.FileID,ImportFiles.ImportId,Files.FileTypeId,Files.Filename,Files.Annotation,Jobs.JobNumber FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ?'
      if topLevelOnly:
        com = com + ' AND Jobs.ParentJobId IS NULL'
      if fileType is not None:
        if fileType == 'application/CCP4-mtz-mini':
          com = com + ' AND Files.FileTypeID IN  (10,11,12,13)'
        elif fileType in FILETYPES_TEXT:
          fileType = FILETYPES_TEXT.index(fileType)
          com = com + ' AND Files.FileTypeID = ?'
          args.append(fileType)
        elif isinstance(fileType,int):
          com = com + ' AND Files.FileTypeID = ?'
          args.append(fileType)
      if subType is not None or contentFlag is not None:
        com,args = self.fileTypeCommand(com,args,subType=subType,contentFlag=contentFlag)
      self.execute(com,args)
      fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,UUIDTYPE,int,str,str,str])
      return fileList

    def deleteFilesOnJobIdAndParamName(self,jobIdParamList=[]):
      #com = 'DELETE FROM Files WHERE (JobId,JobParamName) IN ( ( * , * ) )'
      com = 'DELETE FROM Files WHERE JobId = ?  AND JobParamName LIKE ?'
      self.setDiagnostic(True)
      for ii in range(len(jobIdParamList)):
        #print 'deleteFilesOnJobIdAndParamName',jobIdParamList[ii]
        try:
          self.execute(com,(jobIdParamList[ii][0],jobIdParamList[ii][1]+'%'))
        except:
          print('Failed delete in deleteFilesOnJobIdAndParamName',jobIdParamList[ii])
      self.setDiagnostic(False)
      self.commit()

    def deleteFilesOnJobNumberAndParamName(self,projectId=None,jobNumberParamList=[]):
      #com = 'DELETE FROM Files WHERE (JobId,JobParamName) IN ( ( * , * ) )'
      com = 'DELETE FROM Files WHERE FileId IN (SELECT Files.FileId FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? AND Jobs.JobNumber = ? AND Files.JobParamName LIKE ?)'
      for ii in range(len(jobNumberParamList)):
        #print 'deleteFilesOnJobIdAndParamName',jobNumberParamList[ii]
        try:
          self.execute(com,(projectId,jobNumberParamList[ii][0],jobNumberParamList[ii][1]+'%'))
        except:
          print('Failed delete in deleteFilesOnJobNumberAndParamName',jobNumberParamList[ii])
      self.commit()

    def getProjectOutputFiles(self,projectId=None):
      args = (projectId,)
      #self.execute('SELECT Files.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? AND Files.ImportID IS NULL',args)
      self.execute('SELECT Files.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? AND ImportFiles.ImportId IS NULL',args)
      fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,int,str,str])
      #print 'CDbApi.getProjectOutputFiles',fileList
      '''
      # The following would return sub-jobs reference to files that have been 'passed on' to parent
      args = (projectId,FILE_ROLE_OUT)
      self.execute('SELECT FileUses.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation FROM Files INNER JOIN Jobs ON Files.JobId=Jobs.JobId INNER JOIN FileUses ON Files.FileID=FileUses.FileID WHERE Jobs.ProjectId = ? AND FileUses.RoleID = ?',args)
      fileUsesList = self.fetchAll2PyList([int,int,int,str,str])
      fileList.extend(fileUsesList)
      '''
      return fileList

    def getMultiFollowOnJobs(self,jobIdList,traceImportFiles=True):
      jobTreeList = []
      for jobId in jobIdList:
        jobTreeList.append(self.getFollowOnJobs(jobId,traceImportFiles=traceImportFiles))
      return jobTreeList

    def getPreceedingJobs(self,jobId=None,fileType=None):
      args = [jobId,jobId]
      com = 'SELECT DISTINCT Files.JobId FROM Files INNER JOIN FileUses ON Files.FileId = FileUses.FileId INNER JOIN Jobs ON FileUses.JobID = Jobs.JobID WHERE Jobs.JobID = ? AND Jobs.parentJobId IS NULL AND Files.JobId != ?'
      if fileType is None:
        pass
      elif not isinstance(fileType,list):
        args.append(fileType)
        com = com + ' AND Files.FileTypeId = ?'
      elif len(fileType) == 1:
        args.append(fileType[0])
        com = com + ' AND Files.FileTypeId = ?'
      else:
        args.extend(fileType)
        com = com + ' AND Files.FileTypeId IN (?,?'
        for ii in range(len(fileType)-2): com = com + ',?'
        com = com + ')'
      com = com + ' ORDER BY Jobs.CreationTime'
      self.execute(com,args)
      pJobList = self.fetchAll2Py(UUIDTYPE)
      #print 'getPreceedingJobs',jobId,pJobList
      preceedJobList = []
      for jid in pJobList:
        preceedJobList.append(self.getPreceedingJobs(jid,fileType))
      jobNum = self.getJobInfo(jobId,'jobnumber')
      return [jobId,jobNum,preceedJobList]

    def getSucceedingJobs(self,jobId=None,fileType=None):
      args = [jobId,jobId]
      com = 'SELECT DISTINCT FileUses.JobID FROM FileUses INNER JOIN Files ON FileUses.FileId = Files.FileId INNER JOIN Jobs ON FileUses.JobID= Jobs.JobID WHERE Files.JobId = ? AND FileUses.JobID != ? AND Jobs.parentJobId IS NULL'
      if fileType is None:
        pass
      elif not isinstance(fileType,list):
        args.append(fileType)
        com = com + ' AND Files.FileTypeId = ?'
      elif len(fileType) == 1:
        args.append(fileType[0])
        com = com + ' AND Files.FileTypeId = ?'
      else:
        args.extend(fileType)
        com = com + ' AND Files.FileTypeId IN (?,?'
        for ii in range(len(fileType)-2): com = com + ',?'
        com = com + ')'

      com = com + ' ORDER BY Jobs.CreationTime'
      self.execute(com,args)
      myList = self.fetchAll2Py(UUIDTYPE)
      retList = []
      for childJob in myList:
        retList.append(self.getSucceedingJobs(childJob,fileType))
      jobNum = self.getJobInfo(jobId,'jobnumber')
      #print 'CDbApi.get
      return jobId,jobNum,retList


    def getFollowOnJobs(self,jobId=None,traceImportFiles=True,fileId=None):
      # Get import and output files associated with job and test if subsequently used
      # This is prelude to deleting job
      # Alternatively get jobs which subsequently use a given file of fileId
      if fileId is not None:
        args = (fileId,)
        self.execute('SELECT DISTINCT JobID FROM FileUses WHERE FileId = ?')
        childJobList = self.fetchAll2Py(UUIDTYPE)
        args = (fileId,)
        self.execute('SELECT ImportFiles.ImportId,Files.Filename FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId=Files.FileId WHERE =Files.FileId= ?',args)
        importFileList = self.fetchAll2PyList([UUIDTYPE,str])
      else:
        if traceImportFiles:
          args = (jobId,jobId,)
          self.execute('SELECT DISTINCT FileUses.JobID FROM FileUses INNER JOIN Files ON FileUses.FileId = Files.FileId INNER JOIN Jobs ON FileUses.JobID= Jobs.JobID WHERE Files.JobId = ? AND FileUses.JobID != ? AND Jobs.parentJobId IS NULL ORDER BY Jobs.CreationTime',args)
        else:
          args = (jobId,)
          # Exclude tracing follow on jobs that use imported files
          self.execute('SELECT DISTINCT FileUses.JobID FROM FileUses INNER JOIN Files ON FileUses.FileId = Files.FileId INNER JOIN Jobs ON FileUses.JobID= Jobs.JobID LEFT OUTER JOIN ImportFiles ON ImportFiles.FileId =  Files.FileId WHERE Files.JobId = ? AND Jobs.parentJobId IS NULL AND ImportFiles.ImportId IS NULL ORDER BY Jobs.CreationTime',args)
        """
        # This ought to be redundant due to 'SELECT DISTINCT'
        retList = self.fetchAll2Py(UUIDTYPE)
        childJobList = []
        for item in retList:
          if item not in childJobList: childJobList.append(item)
        """
        childJobList = self.fetchAll2Py(UUIDTYPE)

        args = (jobId,)
        self.execute('SELECT ImportFiles.ImportId,Files.Filename  FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId=Files.FileId INNER JOIN Jobs ON Files.JobID = Jobs.JobID WHERE Files.JobId = ? AND Jobs.parentJobId IS NULL',args)
        importFileList = self.fetchAll2PyList([UUIDTYPE,str])

      #print 'CDbApi.getFollowOnJobs',jobId,childJobList
      jobTree = []
      for childJob in childJobList:
        jobTree.append(self.getFollowOnJobs(childJob))
      #print 'CDbApi.getFollowOnJobs retVal',retVal
      return jobId,importFileList,jobTree

    def getJobSourceOfFiles(self,fileIdList=[]):
      '''
      Get the list of distinct jobs that are the source of the files in fileIdList
      '''
      com = 'SELECT DISTINCT JobId FROM Files WHERE FileId IN (?'
      for j in fileIdList[1:]:
          com = com + ',?'
      self.execute(com+')',fileIdList)
      jobIdList = self.fetchAll2Py(UUIDTYPE)
      return jobIdList

    def setJobToImport(self,jobId=None,projectId=None,importFiles=[]):
      # Delete child jobs
      self.execute("SELECT JobID FROM Jobs WHERE ParentJobID = ?",(jobId,))
      childJobs = self.fetchAll2Py(UUIDTYPE)
      for job in childJobs:
        self._deleteJobTree(job)
      #Delete file and fileuses
      self._deleteOneJob(jobId,saveJob=True,saveImportFiles=True)
      #self.updateJob(jobId,'taskName','import_files')
      self.updateJob(jobId,'status',JOB_STATUS_FILE_HOLDER)
      self.updateJob(jobId,'evaluation',0)
      jobTitle = self.getJobInfo(jobId=jobId,mode='jobtitle')
      #self.updateJob(jobId,'jobTitle','Imported files from deleted: '+jobTitle)

      if projectId is None: projectId = self.getJobInfo(jobId,'projectid')
      self.setJobToImportSignal.emit({'projectId':projectId,'jobId':jobId})

    def getProjectImportFiles(self,projectId=None,ifJobInfo=False):
      args = (projectId,)
      if ifJobInfo:
        self.execute('SELECT ImportFiles.ImportID, Files.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation, Jobs.JobNumber, Jobs.TaskName,ImportFiles.SourceFileName,ImportFiles.CreationTime FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId=Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? ORDER BY Jobs.CreationTime',args)
        fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,UUIDTYPE,int,str,str,str,str,str,float])
      else:
        self.execute('SELECT  ImportFiles.ImportID,Files.JobId,Files.FileID,Files.FileTypeId,Files.Filename,Files.Annotation,ImportFiles.SourceFileName,ImportFiles.CreationTime FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId=Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE Jobs.ProjectId = ?',args)
        fileList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,UUIDTYPE,int,str,str,str,float])
      #print 'CDbApi.getProjectImportFiles',fileList
      # Convert time to readable string
      if len(fileList)>0:
        t = len(fileList[0])-1
        for i in range(len(fileList)):
          fileList[i][t] = self.timeString(fileList[i][t])
      return fileList

    def refineSelectionCommand(self,selcom0,args0,fileTypeId,subType,contentFlag):
        allowAny = False
        if fileTypeId in [4,16]:
          selcom0 +=' AND Files.FileTypeId IN (10,11,12,13)'
        else:
          selcom0 +=' AND Files.FileTypeId= ?'
          args0.append(fileTypeId)
        subType0 = copy.deepcopy(subType)
        #print 'refineSelectionCommand',subType0,type(subType0)
        if subType0 is not None:
          # Check if '0' is in the list - implies allow any but list after the
          # preferred selection
          if isinstance(subType0,list):
            try:
              subType0.remove(0)
              allowAny = True
              if len(subType0) == 1: subType0 = subType0[0]
            except:
              subType0 = subType
            else:
              subType0 = subType

          if not isinstance(subType0,list):
            if subType0 == 0:
              selcom0 += ' AND ( Files.FileSubType IS NULL OR Files.FileSubType = ?)'
            elif subType0 < 0 :
              selcom0 += ' AND Files.FileSubType != ?'
              subType0 = -subType0
            else:
              selcom0 += ' AND Files.FileSubType = ?'
            args0.append(subType0)
          else:
            if 0 in subType0:
              selcom0 += ' AND (Files.FileSubType IS NULL OR Files.FileSubType IN ('
            else:
              selcom0 += ' AND Files.FileSubType IN ('
            for subTypeItem in subType0:
              args0.append(subTypeItem)
              selcom0 +='?,'
            selcom0 = selcom0[0:-1]+')'
            if 0 in subType0: selcom0 += ')'
        if contentFlag is not None:
          if not isinstance(contentFlag,list):
            if contentFlag==0:
              selcom0 += ' AND ( Files.FileContent IS NULL OR Files.FilContent = ?)'
            elif contentFlag < 0:
              selcom0 += ' AND Files.FileContent != ?'
              contentFlag = -contentFlag
            else:
              selcom0 += ' AND Files.FileContent = ?'
            args0.append(contentFlag)
          else:
            if 0 in contentFlag:
              selcom0 += ' AND (Files.FileContent IS NULL OR Files.FileContent IN ('
            else:
              selcom0 += ' AND Files.FileContent IN ('
            for contentFlagItem in contentFlag:
              args0.append(contentFlagItem)
              selcom0 +='?,'
            selcom0=selcom0[0:-1]+ ')'
            if 0 in contentFlag: selcom0 += ')'
        return selcom0,args0,allowAny

    def getJobsWithOutputFiles(self,projectId=None,fileTypeId=None,projectName=None,fileType=None,fileTypeClass=None,
                               subType=None,contentFlag=None,topLevelOnly=True,importFiles=False,jobId=None):

      jobInfo = []
      if projectId is None:
        if projectName is not None:
          projectId = self.getProjectId(projectName=projectName)
      if fileTypeId is None:
        if fileTypeClass is not None and fileTypeClass in FILETYPES_CLASS:
          fileTypeId = FILETYPES_CLASS.index(fileTypeClass)
        elif fileType is not None and fileType in FILETYPES_TEXT:
          fileTypeId = FILETYPES_TEXT.index(fileType)


      # Get list of files generated by jobs in the project (no imported files)
      selcom1 = 'SELECT Jobs.TaskName,Jobs.JobNumber, Files.FileId, Files.Annotation,ImportFiles.SourceFilename,Files.FiletypeID,Files.FileSubType,Files.JobParamName, Jobs.CreationTime FROM Jobs INNER JOIN Files ON Files.JobId=Jobs.JobId LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE '
      if jobId is not None:
        args1 =  [jobId,]
        selcom1 = selcom1 + 'Jobs.JobId=?'
      else:
        args1 =  [projectId,]
        selcom1 = selcom1 + 'Jobs.ProjectId=?'


      selcom,args,allowAny = self.refineSelectionCommand(selcom1,args1,fileTypeId,subType,contentFlag)

      selcom += ' AND ImportFiles.ImportId IS NULL'
      if topLevelOnly: selcom = selcom + ' AND Jobs.ParentJobId IS NULL'
      selcom +=' ORDER BY Jobs.CreationTime DESC'
      #print 'CDbApi.getJobsWithOutputFiles',selcom,args
      self.execute(selcom,args)
      jobList = self.fetchAll2PyList([str,str,UUIDTYPE,str,str,int,int,str,float])
      #print 'getJobsWithOutputFiles',args,jobList

      #print 'getJobsWithOutputFiles args',args
      #print 'getJobsWithOutputFiles first',jobInfo


      # Add imported files to the list - after any created files
      if importFiles:
        args =  [projectId,]
        selcom = 'SELECT Jobs.TaskName,Jobs.JobNumber, Files.FileId, Files.Annotation,ImportFiles.SourceFilename,Files.FiletypeID,Files.FileSubType,Files.JobParamName, Jobs.CreationTime FROM Jobs INNER JOIN Files ON Files.JobId=Jobs.JobId LEFT OUTER JOIN ImportFiles ON Files.FileId=ImportFiles.FileId WHERE '
        if jobId is not None:
          args =  [jobId,]
          selcom = selcom + 'Jobs.JobId=?'
        else:
          args =  [projectId,]
          selcom = selcom + 'Jobs.ProjectId=?'


        selcom,args,allowAny = self.refineSelectionCommand(selcom,args,fileTypeId,subType,contentFlag)

        selcom += ' AND ImportFiles.ImportId IS NOT NULL'
        if topLevelOnly: selcom = selcom + ' AND Jobs.ParentJobId IS NULL'
        selcom +=' ORDER BY Jobs.CreationTime DESC'
        self.execute(selcom,args)
        impList = self.fetchAll2PyList([str,str,UUIDTYPE,str,str,int,int,str,float])
        if len(impList)>0:
          jobList.extend(impList)
          def get_key(item): return item[8]
          jobList = sorted(jobList,key=get_key,reverse=True)

      # Do another search for any files of right type but wrong subtype
      #print 'getJobsWithOutputFiles allowAny',allowAny,subType
      #print 'CDbApi.getJobsWithOutputFiles allowAny selcom1',selcom1,args1
      if subType is not None and allowAny:
        selcom2,args2,allowAny = self.refineSelectionCommand(selcom1,args1[0:1],fileTypeId,None,contentFlag)
        selcom2 += ' AND ImportFiles.ImportId IS NULL'
        if topLevelOnly: selcom2 = selcom2 + ' AND Jobs.ParentJobId IS NULL'
        selcom2 +=' ORDER BY Jobs.CreationTime DESC'

        #print 'CDbApi.getJobsWithOutputFiles allowAny selcom2',selcom2,args2
        self.execute(selcom2,args2)
        jobList.extend(self.fetchAll2PyList([str,str,UUIDTYPE,str,str,int,int,str,float]))


      # Sort it all into a list of dicts
      ll = 0
      insPt = None
      insJob = None
      fileIdList = []
      while ll<len(jobList):
        if jobList[ll][2] in fileIdList:
          ll = ll + 1
        else:
          fileIdList.append(jobList[ll][2])
          if jobList[ll][0] in ['coot_rebuild']:
            if insPt is None or jobList[ll][1] != insJob:
              insPt = ll
              insJob = jobList[ll][1]
          else:
            insPt = None
            insJob = None
          if insPt is not None:
            jobInfo.insert(insPt,{ 'taskname' :  jobList[ll][0], 'jobnumber' :  jobList[ll][1], 'fileid' : jobList[ll][2] , 'annotation' : jobList[ll][3],
                           'importSource' : jobList[ll][4] , 'filetype' : FILETYPES_TEXT[jobList[ll][5]], 'filesubtype':jobList[ll][6] ,
                           'jobparamname' : jobList[ll][7] } )
          else:

            jobInfo.append( {  'taskname' :  jobList[ll][0], 'jobnumber' :  jobList[ll][1], 'fileid' : jobList[ll][2] , 'annotation' : jobList[ll][3],
                           'importSource' : jobList[ll][4] , 'filetype' : FILETYPES_TEXT[jobList[ll][5]], 'filesubtype':jobList[ll][6] ,
                           'jobparamname' : jobList[ll][7] } )
          ll = ll + 1


      #if contentFlag is not None: print 'CDbApi.getJobsWithOutputFiles',jobInfo
      return jobInfo

    def updateFile(self,fileId=None,key=None,value=None):

      if key is None: raise CException(self.__class__,230,'None')
      key = key.lower()
      if not key in ['annotation','jobid']:
        raise CException(self.__class__,230,key)

      jobId,projectId,permission = self.getFilePermission(fileId)
      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      #self.setDiagnostic(True)
      self.execute("UPDATE Files SET "+key+" = ? WHERE FileID = ?",(value,fileId,))
      #self.setDiagnostic(False)
      self.commit()
      self.fileUpdated.emit({'fileId':fileId,'jobId':jobId,'projectId':projectId,'key':key,'value':value})

    '''
    BEWARE Comments (for jobs) and ProjectComments are two different tables but share the following api
    '''

    def getComments(self,jobId=None,projectId=None,userName=None):
      if jobId is not None:
        com = 'SELECT CommentID,UserId,TimeOfComment,Comment FROM Comments WHERE JobID = ?'
        args = [jobId]
      else:
        com = 'SELECT ProjectCommentID,UserId,TimeOfComment,Comment FROM ProjectComments WHERE ProjectID = ?'
        args = [projectId]
      if  userName is not None:
        userId=self.getUserId(userName=userName)
        if userId is not None:
          com += ' AND UserId = ?'
          args.append(userId)
      com+= ' ORDER BY TimeOfComment DESC'
      self.execute(com,args)
      tmpList = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,int,str])
      commentList = []
      for cid,uid,machTime,comment in tmpList:
#FIXME - This may fail, inexplicably and intermittently. Don't know why but tweaking userids of ProjectComments is not working.
        try:
            uName = self.getUserInfo(uid,mode='UserName')
        except:
            self.execute('SELECT UserName FROM Users')
            allNames = self.fetchAll2Py(str)
            uName = allNames[0]
        #timeString = self.timeString(machTime)
        commentList.append([cid,uName,machTime,comment])
      return commentList

    def getCommentInfo(self,jobId=None,projectId=None):
      commentList = self.getComments(jobId=jobId,projectId=projectId)
      if len(commentList) == 0:
        return {}
      else:
        return { 'commentid' : commentList[0][0],
                 'userid' : commentList[0][1],
                 'timeofcomment' : commentList[0][2],
                 'comment' : commentList[0][3] }

    def createComment(self,jobId=None,projectId=None,userName=None,comment=None):

      if jobId is None and projectId is None: raise CException(self.__class__,240,'None')

      if userName is None:
        userId =self.getUserId(self._userName)
      else:
        userId = self.getUserId(userName)
      time = self.currentTime()

      if jobId is not None:
        commentId = self.uniqueId(table='Comments',identifier='CommentID')
        args = [commentId,jobId,userId,time,comment]
        self.execute("INSERT INTO Comments (CommentId,JobID,UserID,TimeOfComment,Comment) VALUES (?,?,?,?,?);",args)
        self.commit()
        self.commentEdited.emit({'jobId':jobId,'commentId':commentId,'comment':comment})
      else:
        commentId = self.uniqueId(table='ProjectComments',identifier='ProjectCommentID')
        args = [commentId,projectId,userId,time,comment]
        self.execute("INSERT INTO ProjectComments (projectCommentId,projectID,UserID,TimeOfComment,Comment) VALUES (?,?,?,?,?);",args)
        self.commit()
        self.projectCommentEdited.emit({'projectId':projectId,'commentId':commentId,'comment':comment})

      return commentId


    def updateComment(self,commentId,comment=None):

      time = self.currentTime()
      self.execute('SELECT JobId FROM Comments WHERE CommentId = ?',(commentId,))
      jobIdList = self.fetchAll2Py(UUIDTYPE)
      if len(jobIdList)==1:
        self.execute("UPDATE Comments SET Comment = ?, TimeOfComment = ? WHERE CommentID = ?",(comment,time,commentId))
        self.commit()
        self.commentEdited.emit({'jobId':jobIdList[0],'commentId':commentId,'comment':comment})
      else:
        self.execute('SELECT ProjectId FROM ProjectComments WHERE ProjectCommentId = ?',(commentId,))
        projectIdList = self.fetchAll2Py(UUIDTYPE)
        if len(projectIdList)==1:
          self.execute("UPDATE ProjectComments SET Comment = ?, TimeOfComment = ? WHERE ProjectCommentID = ?",(comment,time,commentId))
          self.commit()
          self.projectCommentEdited.emit({'projectId':projectIdList[0],'commentId':commentId,'comment':comment})

    def deleteComment(self,commentId):

      self.execute('SELECT JobId FROM Comments WHERE CommentId = ?',(commentId))
      jobIdList = self.fetchAll2Py(UUIDTYPE)
      if len(jobIdList)==1:
        self.execute("DELETE FROM Comments WHERE CommentID = ?",(commentId,))
        self.commit()
        self.commentDeleted.emit({'jobId':jobIdList[0],'commentId':commentId})
      else:
        self.execute('SELECT ProjectId FROM ProjectComments WHERE ProjectCommentId = ?',(commentId))
        projectIdList = self.fetchAll2Py(UUIDTYPE)
        print('deleteComment',projectIdList)
        if len(projectIdList)==1:
          self.execute("DELETE FROM ProjectComments WHERE ProjectCommentID = ?",(commentId,))
          self.commit()
          self.projectCommentDeleted.emit({'projectId':projectIdList[0],'commentId':commentId})

    def currentTime(self):
      '''Return machine time as a float'''
      import time
      return time.time()

    def resetProjectTags(self,projectId=None,tagIdList=[]):
      currentProjectTags = self.getProjectTags(projectId)
      delList = []
      addList = []
      for tid in currentProjectTags:
        if not tid in tagIdList: delList.append(tid)
      for tid in tagIdList:
        if not tid in currentProjectTags: addList.append(tid)
      #print 'resetProjectTags del',delList,'add',addList
      if len(delList)>0:
        com = "DELETE FROM ProjectTags WHERE ProjectID = ? AND TagId IN (?"
        for ii in range(1,len(delList)):com = com + ',?'
        com = com +')'
        delList.insert(0,projectId)
        self.execute(com,delList)
      if len(addList)>0:
        for tid in addList:
          self.createProjectTag(projectId,tid,False)
      if len(delList)>0 or len(addList)>0:
        self.commit()
        self.projectTagsChanged.emit(projectId)

    def getProjectTags(self,projectId=None,tagText=False):
      if tagText:
        self.execute('SELECT ProjectTags.TagId, Tags.Text FROM ProjectTags INNER JOIN Tags ON ProjectTags.TagId=Tags.TagId WHERE ProjectTags.ProjectId = ? ORDER BY Tags.Text',(projectId,))
        ret = self.fetchAll2PyList([UUIDTYPE,str])
        return ret
      else:
        self.execute('SELECT TagId FROM ProjectTags WHERE ProjectId = ?',(projectId,))
        ret = self.fetchAll2Py(UUIDTYPE)
        return ret

    def getProjectsWithTag(self,tagId):
      self.execute('''SELECT ProjectTags.ProjectId, Projects.ProjectName
                    FROM ProjectTags
                    INNER JOIN Projects ON ProjectTags.ProjectId = Projects.ProjectId
                    WHERE ProjectTags.TagId = ?''',(tagId,))
      return self.fetchAll2PyList([UUIDTYPE,str])

    def createProjectTag(self,projectId=None,tagId=None,commit=True):
      if self.projectTagExists(projectId,tagId):
        raise CException(self.__class__,292)
      else:
        self.execute("INSERT INTO ProjectTags (ProjectId,TagId) VALUES (?,?)",(projectId,tagId))
        if commit:
          self.commit()
          print('createProjectTag emit projectTagsChanged',projectId)
          self.projectTagsChanged.emit(projectId)

    def projectTagExists(self,projectId=None,tagId=None):
      self.execute('SELECT ProjectId,TagId FROM ProjectTags WHERE ProjectId = ? AND TagId = ?',(projectId,tagId))
      ret = self.fetchAll2PyList([UUIDTYPE,UUIDTYPE])
      #print 'projectTagExists',ret
      return len(ret)>0

    def deleteProjectTag(self,projectId=None,tagId=None):
       self.execute("DELETE FROM ProjectTags WHERE ProjectID = ? AND TagId = ?",(projectId,tagId))
       self.commit()
       self.projectTagsChanged.emit(projectId)

    def createTag(self,text,parentTagId=None):
      self.execute("SELECT tagId FROM Tags WHERE text = ?",[text,])
      ret = self.fetchAll2Py(UUIDTYPE)
      #print 'createTag',text,ret
      if len(ret)>0:
        raise CException(self.__class__,294,ret)
      tagId = self.uniqueId(table='Tags',identifier='TagID')
      self.execute("INSERT INTO Tags (TagId,text,parentTagId) VALUES (?,?,?)",[tagId,text,parentTagId])
      self.commit()
      self.tagCreated.emit({'tagId':tagId,'parentTagId':parentTagId,'text':text})
      return tagId

    def updateTag(self,tagId,key=None,value=None):
      key = key.lower()
      if not key in ['parenttagid','text']:
        raise CException(self.__class__,293,key)
      # Check is unique text
      if key == 'text':
        self.execute("SELECT tagId FROM Tags WHERE text = ?",[value,])
        ret = self.fetchAll2Py(UUIDTYPE)
        #print 'updateTag',value,ret
        if len(ret)>0:
          raise CException(self.__class__,294,ret)
      self.execute("UPDATE Tags SET "+key+" = ? WHERE TagID = ?",(value,tagId,))
      self.commit()
      self.tagUpdated.emit({'tagId':tagId,'key':key,'value':value})
      for pid,pname in self.getProjectsWithTag(tagId):
        self.projectTagsChanged.emit(pid)

    def deleteTag(self,tagId,deleteChildren=False):
      self.execute("SELECT ParentTagId FROM Tags WHERE TagId = ?",[tagId,])
      parentId = self.fetchAll2Py(UUIDTYPE)
      #print 'CDbApi.deleteTag parentId',parentId
      if parentId is not None:
        if deleteChildren:
          pass
        else:
          self.execute("UPDATE Tags SET ParentTagId = NULL  WHERE ParentTagID = ?",(parentId))
      projectsWithTag = self.getProjectsWithTag(tagId)
      self.execute("DELETE FROM ProjectTags WHERE TagID = ?",(tagId,))
      self.execute("DELETE FROM Tags WHERE TagID = ?",(tagId,))
      self.commit()
      self.tagDeleted.emit({'tagId' : tagId})
      for pid,pname in projectsWithTag:
        self.projectTagsChanged.emit(pid)

    def deleteUnusedTags(self):
      self.execute('DELETE FROM Tags WHERE Tags.TagId NOT IN (SELECT DISTINCT TagId FROM ProjectTags)',())
      self.commit()
      self.unusedTagsDeleted.emit()

    def getTagList(self):
      #self.execute('SELECT Tags.TagId, Tags.ParentTagId, Tags.Text, COUNT(ProjectTags.TagId) FROM Tags LEFT OUTER JOIN ProjectTags ON ProjectTags.TagId = Tags.TagId')
      #self.execute('SELECT Tags.TagId, Tags.ParentTagId, Tags.Text, COUNT(ProjectTags.TagId) FROM ProjectTags LEFT OUTER JOIN Tags ON ProjectTags.TagId = Tags.TagId')
      self.execute('SELECT Tags.TagId, Tags.ParentTagId, Tags.Text FROM Tags ORDER BY Tags.Text')
      return self.fetchAll2PyList([UUIDTYPE,UUIDTYPE,str])

    def getProjectTagList(self):
      self.execute('SELECT DISTINCT ProjectTags.ProjectId,Tags.Text FROM ProjectTags INNER JOIN Tags ON ProjectTags.TagId = Tags.TagId INNER JOIN Projects ON ProjectTags.ProjectId = Projects.ProjectId ORDER BY Projects.ProjectName, Tags.Text')
      return self.fetchAll2PyList([UUIDTYPE,str])

    def newTag(self,projectId=None,text=None,parentTagId=None):
      tagId = self.createTag(text=text,parentTagId=parentTagId)
      #print 'CDbApi.newTag',tagId
      if projectId is not None:
        self.createProjectTag(projectId,tagId)
      return tagId

    def getProjectCommentsList(self):
      self.execute('SELECT DISTINCT ProjectId, Comment FROM ProjectComments')
      return self.fetchAll2PyList([UUIDTYPE,str])

    def setDirectoryAlias(self,alias=None,directory=None):

      if not self._userRole in [USER_ROLE_OWNER]:
        raise CException(self.__class__,106,'You are '+self._userName+' calling createDirectoryAlias')

      if directory is None:
        raise CException(self.__class__,133,alias)
      if alias != 'CCP4I2_TOP':
        matchId,matchDir,relPath = self.matchProjectDirectory(directory)
        if matchDir is not None and len(relPath)==0:
          raise CException(self.__class__,133,alias)

      try:
        diry = self.getAliasDirectory(alias)
      except:
        diry = None

      if diry is None:
        args = (alias,directory)
        self.execute("INSERT INTO DirectoryAliases (DirectoryAlias,Directory) VALUES (?,?);",args)
      else:
        args = (directory,alias)
        self.execute("UPDATE DirectoryAliases SET Directory= ? WHERE DirectoryAlias= ?",args)
      self.commit()

    def deleteDirectoryAlias(self,alias):
      # TBD
      pass

    def getAliasDirectory(self,alias=None):
      if alias is None:
        raise CException(self.__class__,134,str(alias))
      args = (str(alias),)
      self.execute('SELECT Directory FROM DirectoryAliases WHERE DirectoryAlias = ?',args)
      rv = self.fetchAll2Py(str)
      if len(rv)==0:
        raise CException(self.__class__,134,alias)
      else:
        return rv[0]

    def listDirectoryAliases(self,toTerm=False):
       self.execute("SELECT DirectoryAlias FROM DirectoryAliases")
       aliasList = self.fetchAll2Py(str)
       if toTerm:
         print('listDirectoryAliases')
         for item in aliasList: print (item)
       return aliasList

    def getRecentlyFinishedJobs(self,after=0.0):
      #self.execute('SELECT JobId FROM Jobs WHERE (status='+str(JOB_STATUS_FAILED)+' OR status='+str(JOB_STATUS_FINISHED)+') AND finishtime>'+str(after))
      #print 'CDbApi.getRecentlyFinishedJobs'
      args = (str(after),)
      self.execute('SELECT JobId,ProjectId,Status,JobTitle,FinishTime,ParentJobId FROM Jobs WHERE finishtime> ? ORDER BY finishtime',args)
      finishedJobs = self.fetchAll2PyList((UUIDTYPE,UUIDTYPE,int,str,float,UUIDTYPE))
      #print 'getRecentlyFinishedJobs jobList',finishedJobs

      for jobId,projectId,status,jobTitle,finishTime,parentJobId in finishedJobs:
        #print 'CDbApi.getRecentlyFinishedJobs',jobId
        if status == JOB_STATUS_TO_DELETE:
          self.jobToDelete.emit({'jobId':jobId,'projectId':projectId,'status':status,'jobTitle':jobTitle,'finishTime':finishTime,'parentJobId':parentJobId})
        else:
          self.jobFinished.emit({'jobId':jobId,'projectId':projectId,'status':status,'jobTitle':jobTitle,'finishTime':finishTime,'parentJobId':parentJobId})
          self.considerUpdatingFollowFrom(projectId,jobId)

      return finishedJobs

    def getRecentlyStartedJobs(self,after=0.0):
      #self.execute('SELECT JobId FROM Jobs WHERE (status='+str(JOB_STATUS_FAILED)+' OR status='+str(JOB_STATUS_FINISHED)+') AND finishtime>'+str(after))
      #print 'getRecentlyFinishedJobs'
      args = (str(after),str(JOB_STATUS_RUNNING))
      self.execute('SELECT JobId,ProjectId,Status,ParentJobId,TaskName FROM Jobs WHERE creationtime> ? AND STATUS= ? ORDER BY creationtime',args)
      startedJobs = self.fetchAll2PyList((UUIDTYPE,UUIDTYPE,int,UUIDTYPE,str))
      #print 'getRecentlyStartedJobs jobList',startedJobs

      for jobId,projectId,status,parentJobId,taskName in startedJobs:
        #print 'CDbApi.getRecentlyFinishedJobs',jobId
        self.jobStarted.emit({'jobId':jobId,'projectId':projectId,'status':status,'parentJobId':parentJobId,'taskName':taskName})

      return startedJobs


    def getRunningSubJob(self,jobId):
      # Return the jobId and taskName of a running subJob.  Return nothing if there is >1 running sub-jobs
      # This is currently used in support of interupting jobs and multiple running sub-jobs not supported
      # in that context
      self.execute('SELECT JobId,TaskName FROM  Jobs WHERE parentJobId = ? AND STATUS = ?',(jobId,JOB_STATUS_RUNNING) )
      subJobs = self.fetchAll2PyList((UUIDTYPE,str))
      #print 'getRunningSubJob',subJobs
      if len(subJobs) == 1:
        return { 'jobId' : subJobs[0][0] , 'taskName' : subJobs[0][1] }
      else:
        return None

    def getRunningJobs(self,remote=False):
      if not remote:
        self.execute('SELECT JobId,JobNumber,TaskName,ProjectId,ProcessId,parentJobId FROM Jobs WHERE STATUS = ?',(str(JOB_STATUS_RUNNING)))
      else:
        self.execute('SELECT Jobs.JobId,Jobs.JobNumber,Jobs.TaskName,Jobs.ProjectId,ServerJobs.ServerProcessId,Jobs.parentJobId FROM Jobs INNER JOIN ServerJobs ON Jobs.JobId=ServerJobs.JobId WHERE Jobs.STATUS = ?',(str(JOB_STATUS_REMOTE)))
      ret = self.fetchAll2PyList((UUIDTYPE,str,str,UUIDTYPE,int,UUIDTYPE))
      return ret

    def jobDirectory(self,jobId=None,jobNumber=None,projectId=None,projectDirectory=None):
      import os

      if jobNumber is None or (projectId is None and projectDirectory is None):
        jobInfo = self.getJobInfo(jobId=jobId,mode=['projectid','jobnumber'])
        jobNumber = jobInfo['jobnumber']
        if projectDirectory is None:
          projectDirectory = self.getProjectInfo(projectId=jobInfo['projectid'],mode='projectdirectory')
      elif projectDirectory is None:
        projectDirectory = self.getProjectInfo(projectId=projectId,mode='projectdirectory')

      directory =  os.path.join(projectDirectory,'CCP4_JOBS')

      for job in jobNumber.split('.'):
       directory = os.path.join(directory,'job_'+job)
      #print 'CDbApi.jobDirectory',jobId,jobNumber,directory

      return os.path.normpath(directory)

    def jobRelPath(self,jobId=None,jobNumber=None):
      if jobNumber is None:
        jobNumber = self.getJobInfo(jobId=jobId,mode='jobnumber')
      directory = 'CCP4_JOBS'
      for job in jobNumber.split('.'):
       directory = os.path.normpath(os.path.join(directory,'job_'+job))
      return directory


    def createXData(self,jobId=None,dataClass=None,dataXml=None,projectId=None):

      if projectId is not None:
        permission = self.getProjectPermission(projectId=projectId)
      elif jobId is not None:
        projectId,permission = self.getJobPermission(jobId=jobId)
      else:
        raise CException(self.__class__,161,str(jobId))

      if permission < PRIVILEGE_WRITE:
        raise CException(self.__class__,112)

      dataId = self.uniqueId(table='XData',identifier='XDataID')
      args = [dataId,jobId,dataClass,dataXml]
      self.execute("INSERT INTO XData (XDataID,JobID,XDataClass,XDataXml) VALUES (?,?,?,?)",args)

      self.commit()
      return dataId

    def getXData(self,dataClass=None,projectId=None,jobId=None):
      # Is there an xdata from the job?
      if jobId is None:
        if dataClass is None:
          args = (projectId,)
          self.execute('SELECT XData.XDataId FROM XData INNER JOIN Jobs ON XData.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? ',args)
        else:
          args = (dataClass,projectId)
          self.execute('SELECT XData.XDataId FROM XData WHERE XDataClass = ? INNER JOIN Jobs ON XData.JobId=Jobs.JobId WHERE Jobs.ProjectId = ? ',args)
      else:
        if dataClass is None:
          args = (jobId,)
          self.execute('SELECT XDataId FROM XData WHERE JobId = ? ',args)
        else:
          args = (dataClass,jobId)
          self.execute('SELECT XDataId FROM XData WHERE XDataClass = ? AND JobId = ? ',args)
      rv = self.fetchAll2Py(UUIDTYPE)
      return rv


    def getXDataInfo(self,xDataId=None,mode='all'):
      #self.setDiagnostic(True)
      args = (xDataId,)
      if mode.lower() in ['xdataclass','class','xdataxml','xml']:
        if mode.lower() in ['xdataclass','class']:
          mode = 'xdataclass'
        else:
          mode = 'xdataxml'
        self.execute('SELECT '+mode+' FROM XData WHERE XDataID = ?',args)
        rv = self.fetchAll2Py(str)
      else:
        self.execute('SELECT XDataClass,XDataXml FROM XData WHERE XDataID = ?',args)
        rv = self.fetchAll2PyList([str,str])
      rv = self.fetchAll2Py(str)
      #self.setDiagnostic(False)
      if len(rv)>0:
        return rv[0]
      else:
        return []


    def getXDataByJobContext(self,contextJobId=None,dataClass=None,projectId=None):

      #print 'CDbApi.getXDataByContext',contextJobId,dataClass,projectId

      # Is there an xdata from the job?
      args = (contextJobId,dataClass)
      self.execute('SELECT XDataId FROM XData WHERE JobId = ? AND XDataClass = ?',args)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)>0:
        return rv

      # We need to find preceeding jobs
      args = (contextJobId,)
      self.execute('SELECT PreceedingJobID FROM Jobs WHERE JobId = ?',args)
      rv = self.fetchAll2Py(UUIDTYPE)
      if len(rv)==1:
        return self.getXDataByJobContext(contextJobId=rv[0],dataClass=dataClass,projectId=projectId)

      return []



    def getJobXData(self,jobId=None,role=FILE_ROLE_OUT,mode='xdataid'):
      outputList = []

      if role == FILE_ROLE_OUT:
        # Search the Files table for output files
        args = [jobId,]
        seleList = 'JobID = ?'

        if mode == 'xdataid':
          self.execute('SELECT xdataid FROM XData WHERE '+seleList,args)
          outputList = self.fetchAll2Py(UUIDTYPE)

          return outputList
      return []

    def getKeyTypeId(self,keyTypeName):
      u = (keyTypeName,)
      self.execute('SELECT KeyTypeID FROM KeyTypes WHERE KeyTypeName = ?',u)
      rv = self.fetchAll2Py(int)
      #print 'getKeyTypeId',keyTypeName,rv
      if len(rv)==0:
        raise CException(self.__class__,290,keyTypeName)
      else:
        return rv[0]

    def createJobKeyValue(self,jobId=None,keyTypeId=None,keyTypeName=None,value=None):
      if keyTypeId is None and keyTypeName is not None:
        keyTypeId = self.getKeyTypeId(keyTypeName)

      if isinstance(value,str):
        args = (jobId,keyTypeId,value)
        self.execute("INSERT INTO JobKeyCharValues (JobId,KeyTypeId,Value) VALUES (?,?,?)",args)
      else:
        args = (jobId,keyTypeId,float(value))
        #print 'createJobKeyValues args',args
        self.execute("INSERT INTO JobKeyValues (JobId,KeyTypeId,Value) VALUES (?,?,?)",args)
      self.commit()

    def updateJobKeyValue(self,jobId=None,keyTypeId=None,keyTypeName=None,value=None):
      if keyTypeId is None and keyTypeName is not None:
        keyTypeId = self.getKeyTypeId(keyTypeName)
      if isinstance(value,str):
        self.db.execute("UPDATE JobKeyCharValues SET Value = ? WHERE JobID = ? AND KeyTypeId = ?",(value,jobId,keyTypeId))
      else:
        self.db.execute("UPDATE JobKeyValues SET Value = ? WHERE JobID = ? AND KeyTypeId = ?",(float(value),jobId,keyTypeId))
      self.commit()

    def deleteJobKeyValue(self,jobId=None,keyTypeId=None,keyTypeName=None):
      if keyTypeId is None and keyTypeName is not None:
        keyTypeId = self.getKeyTypeId(keyTypeName)
      self.db.execute('DELETE FROM JobKeyValues WHERE JobId = ? AND KeyTypeId = ?',(jobId,keyTypeId,))
      self.db.execute('DELETE FROM JobKeyCharValues WHERE JobId = ? AND KeyTypeId = ?',(jobId,keyTypeId,))

    def getFileAssociationTypeId(self,fileAssociationTypeName):
      u = (fileAssociationTypeName,)
      self.execute('SELECT FileAssociationTypeID FROM FileAssociationTypes WHERE FileAssociationTypeName = ?',u)
      rv = self.fetchAll2Py(int)
      if len(rv)==0:
        raise CException(self.__class__,291,fileAssociationTypeName)
      else:
        return rv[0]

    def createFileAssociation(self,typeId=None,typeName=None,associationName=None,fileList=[]):
      if typeId is None and typeName is not None:
        typeId = self.getFileAssociationTypeId(typeName)
      pid = self.uniqueId(table='FileAssociations',identifier='FileAssociationID')

      args = [pid,typeId,associationName]
      self.execute("INSERT INTO FileAssociations (FileAssociationID,FileAssociationTypeID,FileAssociationName) VALUES (?,?,?)",args)

      if len(fileList)>0:
        for fileItem in fileList:
          if not isinstance(fileItem,(list,tuple)): fileItem = [fileItem,0]
          self.createFileAssociationMember(self,fileAssociationId=pid,fileId=fileItem[0],roleId=fileItem[1],commit=False)
      self.commit()
      return pid

    def deleteFileAssociation(self,fileAssociationId):
      self.db.execute('DELETE FROM FileAssociations WHERE FileAssociationId = ?',(fileAssociationId,))
      self.db.execute('DELETE FROM FileAssociationMembers WHERE FileAssociationId = ?',(fileAssociationId,))
      self.commit()

    def addFileAssociationMember(self,fileAssociationId=None,fileId=None,roleId=0,commit=True):
      self.execute("INSERT INTO FileAssociationMembers (FileID,FileAssociationID,RoleId) VALUES (?,?,?)",
                   (fileId,fileAssociationId,roleId) )
      if commit: self.commit()

    def removeFileAssociationMember(self, fileAssociationId=None,fileId=None,commit=True):
      self.db.execute('DELETE FROM FileAssociationMembers WHERE FileId = ? AND FileAssociationId = ?',(fileId,fileAssociationId))
      if commit: self.commit()

    def getFileAssociationInfo(self,fileAssociationId=None,fileId=None,typeId=None,roleId=None):

      com = 'SELECT DISTINCT FileAssociations.FileAssociationId,FileAssociations.FileAssociationTypeId,FileAssociation.FileAssociationName,FileAssociationMembers.FileId,FileAssociationMembers.FileAssociationRoleID FROM FileAssociations INNER JOIN FileAssociationMembers ON FileAssociations.FileAssociationID = FileAssociationMembers.FileAssociationID'
      if fileAssociationId is not None:
         com = com + ' WHERE FileAssociations.FileAssociationId=?'
         args = (fileAssociationId,)

      else:
        com = com + ' WHERE FileAssociationMembers.FileID=?'
        args = [fileId,]
        if roleId is not None:
          com = com + ' AND FileAssociationMembers.FileAssociationRoleID = ?'
          args.append(roleId)
        if typeId is not None:
          com = com + ' AND FileAssociations.FileAssociationTypeID = ?'
          args.append(typeId)

      self.db.execute(com,args)
      rv = self.fetchAll2PyList(UUIDTYPE,int,str,UUIDTYPE,int)
      #print 'getFileAssociationInfo rv',rv

      info = {}
      for assocId,assocTypeId,assocName,fileId,roleId in rv:
        if assocId not in info:
          info[assocId] = { 'associationtypeid': assocTypeId, 'associationname':assocName, 'members' : [] }
        info[assocId]['members'].append( { 'fileId' : fileId, 'roleId': roleId } )

      return info




    def createProjectExport(self,projectId=None,projectExportAfter=None,projectExportSelection=None):
      pid = self.uniqueId(table='ProjectExports',identifier='ProjectExportID')
      exportTime = self.currentTime()
      #args = [pid,projectId,exportTime,projectExportAfter,projectExportSelection]
      #self.execute("INSERT INTO ProjectExports (ProjectExportID,ProjectID,ProjectExportTime,ProjectExportAfter,ProjectExportSelection) VALUES (?,?,?,?,?)",args)
      args = [pid,projectId,exportTime,projectExportAfter]
      self.execute("INSERT INTO ProjectExports (ProjectExportID,ProjectID,ProjectExportTime,ProjectExportAfter) VALUES (?,?,?,?)",args)
      self.commit()
      return pid

    def getProjectExportInfo(self,projectExportId=None,projectId=None,returnList=False):
      if projectExportId is not None:
        self.execute('SELECT ProjectExportID,ProjectExportTime,ProjectExportAfter FROM ProjectExports WHERE ProjectExportId=?',(projectExportId,))
      elif projectId is not None:
        self.execute('SELECT ProjectExportID,ProjectExportTime,ProjectExportAfter FROM ProjectExports WHERE ProjectId=? ORDER BY ProjectExportTime DESC',(projectId,))
      rv = self.fetchAll2PyList((UUIDTYPE,float,float))
      if returnList or projectId is not None:
        return rv
      else:
        ret={}
        ii=0
        for item in ('projectexportid','projectexporttime','projectexportafter'):
          ret[item.lower()] = rv[0][ii]
          ii += 1
        return ret

    def getProjectImportInfo(self,projectImportId=None,projectId=None,returnList=False):
      databaseHostNames = {}
      if projectImportId is not None:
        self.execute('SELECT ProjectImportID,ProjectImportTime,ProjectExportTime,ProjectExportDatabaseId FROM ProjectImports WHERE ProjectImportId=?',(projectImportId,))
      elif projectId is not None:
        self.execute('SELECT ProjectImportID,ProjectImportTime,ProjectExportTime,ProjectExportDatabaseId FROM ProjectImports WHERE ProjectId=? ORDER BY ProjectImportTime DESC',(projectId,))
      rv = self.fetchAll2PyList((UUIDTYPE,float,float,UUIDTYPE))
      if returnList or projectId is not None:
        # Substitute database hostname for the databaseid
        databaseHostNames = {}
        for ii in range(len(rv)):
          if databaseHostNames.get(rv[ii][3],None) is None:
            databaseHostNames[rv[ii][3]] = self.getDatabaseInfo(databaseId=rv[ii][3])['creationhostname']
          rv[ii][3] = databaseHostNames[rv[ii][3]]
        return rv
      else:
        ret={}
        ii=0
        for item in ('projectimportid','projectimporttime','projectexporttime'):
          ret[item.lower()] = rv[0][ii]
          ii += 1
        ret['projectexportdatabase'] = self.getDatabaseInfo(databaseId= rv[0][ii])['creationhostname']
        return ret

    def createDatabase(self, databaseId=None, creationTime=None, creationHostName=None, 
                       createUserName=None, schemaVersion=None, schemaDate=None):
      # Record details of another database from which we have imported data
      loaded = self.idExists(idValue=databaseId,table='Databases',idName='DatabaseID')
      if loaded: return True
      args = [databaseId,creationTime,creationHostName,createUserName,schemaVersion,schemaDate]
      self.execute("INSERT INTO Databases (DatabaseId,CreationTime,CreationHostName,CreationUserName,SchemaVersion,SchemaDate) VALUES (?,?,?,?,?,?)",args)
      self.commit()
      return False

    def getDatabaseInfo(self,databaseId=None,returnList=False):
      if databaseId is None:
        # Assume we want info for the current db
        self.execute('SELECT DatabaseId,CreationTime,CreationHostName,CreationUserName,SchemaVersion,SchemaDate FROM Databases WHERE ThisIsMe=?',(1,))
      else:
        self.execute('SELECT DatabaseId,CreationTime,CreationHostName,CreationUserName,SchemaVersion,SchemaDate FROM Databases WHERE DatabaseId=?',(databaseId,))

      rv = self.fetchAll2PyList(self.DATABASETYPES)
      if len(rv)!=1:
        raise CException(self.__class__,280)
      if returnList:
        return rv
      else:
        ret={}
        ii=0
        for item in self.DATABASEITEMS:
          ret[item.lower()] = rv[0][ii]
          ii += 1
        return ret

    def createProjectImport (self,projectId=None,projectExportId=None,projectExportDatabaseId=None,projectExportTime=None,projectExportAfter=None):
      if projectId is None or not self.idExists(projectId,'Projects'):
        pass

      if projectExportId is None or projectExportDatabaseId is None or not self.idExists(projectExportDatabaseId,'Databases'):
        #print 'CDbApi.createProjectImport invalid projectExportDatabaseId',projectExportDatabaseId
        return
      pid = self.uniqueId(table='ProjectImports',identifier='ProjectImportID')
      importTime = self.currentTime()

      args = [pid,projectId,importTime,projectExportId,projectExportTime,projectExportDatabaseId]
      self.execute("INSERT INTO ProjectImports (ProjectImportID,ProjectID,ProjectImportTime,ProjectExportId,projectExportTime,projectExportDatabaseId) VALUES (?,?,?,?,?,?)",args)
      self.commit()
      return pid

    def removeTempTables(self):
      for table in ['TempJobs','TempFiles','TempImportFiles','TempExportFiles','TempFileUses','TempXData','TempComments','TempJobKeyValues','TempJobKeyCharValues','TempFileAssociations','TempFileAssociationMembers','TempTags','TempProjectTags','TempProjectComments']:
        self.execute('DROP TABLE IF EXISTS '+table)
      self.commit()

    def getJobTree(self,jobId=None,projectId=None,ifBack=True):
      # Trace back or forward from a given job -
      # beware does not deal with 'export file - run non-ccp4 program - reimport a file' scenario
      # From johndpope's reply on http://stackoverflow.com/questions/7456957/basic-recursive-query-on-sqlite3
      # (also at http://dje.me/2011/03/26/sqlite-data-trees.html)
      # This query be better done by a WITH statement (CommonTableExpressions) which is supported in later sqlite (after 3.8.3) but this is not yet in Python
      # The code here is a fudge only likely to work for sqlite and needs to be sqlite 3.6.18 or later

      """
      self.execute('SELECT Jobs.JobId, Jobs.JobNumber, Jobs.TaskName, Files.FileId, Files.JobParamName FROM Jobs INNER JOIN Files ON Jobs.JobId = Files.JobId INNER JOIN FileUses ON Files.FileId=FileUses.FileId WHERE FileUses.JobId = ?',(jobId,))
      ret = self.fetchAll2PyList([str,str,str,str,str])
      print 'getJobTree',ret
      """

      self.execute('PRAGMA recursive_triggers = TRUE;')

      # The critical selection command depends on wether going back
      selcom = ''' SELECT Jobs.JobId, Jobs.JobNumber, Jobs.TaskName, Jobs.CreationTime,Jobs.ProjectId,
                      Files.FileId, Files.JobParamName,Files.FiletypeID, FileUses.JobId  FROM Jobs
                      '''
      if ifBack:
        # tracing back to preceding jobs
        selcom += '''INNER JOIN Files ON Jobs.JobId = Files.JobId
                     INNER JOIN FileUses ON Files.FileId=FileUses.FileId
                     WHERE FileUses.JobId=new.JobId;'''
      else:
        # tracing subsequent jobs
        selcom += '''INNER JOIN FileUses ON FileUses.JobId = Jobs.JobId
                     INNER JOIN Files ON  FileUses.FileId = Files.FileID
                     WHERE Files.JobId=new.JobId;'''


      # Create a temporary table
      self.execute('DROP TABLE IF EXISTS TempJobTree')
      self.execute('''CREATE TEMP TABLE TempJobTree ( JobId VARCHAR(32), JobNumber VARCHAR(50),TaskName VARCHAR(100), CreationTime REAL, ProjectId  VARCHAR(32),
                                                  FileID VARCHAR(32), JobParamName VARCHAR(32),FiletypeID VARCHAR(32), ChildJobId VARCHAR(32), UNIQUE (JobId)); ''')

      # Create trigger
      self.execute('''CREATE TRIGGER getJobTree AFTER INSERT ON TempJobTree BEGIN
                      INSERT OR IGNORE INTO TempJobTree  (JobId,JobNumber,TaskName,CreationTime,ProjectId,FileID,JobParamName,FiletypeID,ChildJobId) ''' +
                     selcom + 'END;' )

      # Insert the fist jobId into the temporary table - this kicks off the trigger
      self.execute('''INSERT INTO TempJobTree (JobId,JobNumber,TaskName,CreationTime,ProjectId)
                      SELECT JobId,JobNumber,TaskName,CreationTime,ProjectId FROM Jobs WHERE JobId=?''',(jobId,))

      #self.execute('SELECT JobID,JobNumber,TaskName,CreationTime,FileID,JobParamName,ChildJobId FROM TempJobTree ORDER BY CreationTime')
      self.execute('SELECT JobID,JobNumber,TaskName,CreationTime,ProjectId,FileID,JobParamName,FiletypeID,ChildJobId FROM TempJobTree')
      ret = self.fetchAll2PyList([str,str,str,float,str,str,str,int,str])
      #for item in ret: print 'getJobTree',item
      self.execute('DROP TABLE TempJobTree')

      return ret

    def createServerJob(self,jobId,params={}):
      args = [jobId]
      #print 'CDbApi.createServerJob',jobId,params
      for key in ['machine','username','mechanism','remotePath','customCodeFile','validate','keyFilename','serverGroup']:
        args.append(getattr(params,key,None))
      self.execute('INSERT INTO ServerJobs (JobId,Machine,Username,Mechanism,RemotePath,CustomCodeFile,Validate,KeyFilename,ServerGroup) VALUES (?,?,?,?,?,?,?,?,?)',args)
      self.commit()

    def deleteServerJob(self,jobId,commit=True):
      self.execute('DELETE FROM ServerJobs  WHERE JobId = ?',(jobId,))
      if commit: self.commit()

    def updateServerJob(self,jobId,name,value):
      self.execute('UPDATE ServerJobs SET '+name+'= ? WHERE JobId = ?',(value,jobId))
      self.commit()

    def getServerJobs(self):

      self.execute('SELECT JobId,Machine,Username,Mechanism,RemotePath,CustomCodeFile,Validate,KeyFilename,ServerProcessId,ServerGroup FROM ServerJobs')
      ret = self.fetchAll2PyList([UUIDTYPE,str,str,str,str,str,str,str,int,str])
      return ret




class CDbXml(QtCore.QObject):

  fileDeleted = QtCore.Signal(str)

  Instances = set()

  ERROR_CODES = { 1 : { 'description' : 'Failed to read project file for import' },
                  2 : { 'description' : 'Attempting to import project with existing project name' },
                  3 : { 'description' : 'Error creating new project when attempting to import project' },
                  4 : { 'description' : 'Error importing job - job with same number already exists in project' },
                  5 : { 'description' : 'Error retreiving jobid/jobnumber/taskname data for job' },
                  6 : { 'description' : 'Error importing job - job with same jobId already exists in project' },
                  7 : { 'description' : 'Error importing ImportFile record - FileId not found in database','severity' : SEVERITY_WARNING  },
                  8 : { 'description' : 'Error importing ExportFile record - FileId not found in database','severity' : SEVERITY_WARNING  },
                  9 : { 'description' : 'Warning importing job - new job number assigned to avoid clash', 'severity' : SEVERITY_WARNING},

                 10 : { 'description' : 'Warning - attempting to import job already in database','severity' : SEVERITY_WARNING },
                100 : { 'description' : 'No project data found in the database file' },
                101 : { 'description' : 'Unknown error creating project in database' },
                102 : { 'description' : 'Project successfully imported into database', 'severity' : SEVERITY_OK },
                103 : { 'description' : 'Project import to database is NOT commited','severity' : SEVERITY_WARNING },
                104 : { 'description' : ' ' },
                105 : { 'description' : ' ' },
                106 : { 'description' : ' ' },
                200 : { 'description' : 'Expected table not found in import file', 'severity' : SEVERITY_WARNING },
                201 : { 'description' : 'Loading file failed to find mapping for jobId ' },
                202 : { 'description' : 'Loading file failed to find mapping for importId ' },
                203 : { 'description' : 'Loading file failed to find mapping for fileId ' },
                204 : { 'description' : 'Loading file failed to find mapping for exportId' },
                205 : { 'description' : 'Loading file failed to find/interpret xdataxml' },
                210 : { 'description' : 'Error entering job data to database' },
                211 : { 'description' : 'Error entering file data to database' },
                212 : { 'description' : 'Error interpreting data - wrong datatype' },
                213 : { 'description' : 'Error inserting row into database' },
                214 : { 'description' : 'Error inserting comment with unrecognised UserId', 'severity' : SEVERITY_WARNING },
                215 : { 'description' : 'Error inserting comment with unrecognised JobId' },
                216 : { 'description' : 'Error inserting file with unrecognised JobId' },
                217 : { 'description' : 'Error inserting fileUse with unrecognised JobId' },
                218 : { 'description' : 'Error inserting fileUse with unrecognised FileId' },
                219 : { 'description' : 'Error inserting XData with unrecognised JobId' },
                120 : { 'description' : 'Error opening project database xml file' },
                121 : { 'description' : 'Error reading project database xml file' },
                122 : { 'description' : 'Error reading database/export info from project database xml file' },
                123 : { 'description' : 'Error processing job data from project database xml file' }
                      }
  COMMIT_POLICY_NO_ERRORS = 0
  COMMIT_POLICY_FORCE = 1
  COMMIT_POLICY_NO_COMMIT = 2

  @staticmethod
  def updateInstances():
    CDbXml.Instances = set([obj for obj in CDbXml.Instances if isAlive(obj)])
    print('CDbXml.updateInstances',CDbXml.Instances)


  def __init__(self,db=None,xmlFile=None,diagnostic=False):
    self.xmlFile = xmlFile
    self.db=db
    self.projectDirectory = None
    self.projectName = None
    self.projectId = None
    self.projectInfo = {}

    self.errReport = CErrorReport()
    self._diagnostic = diagnostic
    self.stats = {}   # Collect stats to report to user
    self.stats['incrJobNumber'] = 0

    self.PROJECTITEMS = []
    self.PROJECTITEMS.extend(self.db.PROJECTITEMS0)
    self.PROJECTITEMS.extend(['username','parentprojectname'])
    self.PROJECTTYPES = []
    self.PROJECTTYPES.extend(self.db.PROJECTTYPES0)
    self.PROJECTTYPES.extend([str,str])
    self.makeAllSql()

  def setDiagnostic(self,mode=False):
    if mode in (True,False):
      self._diagnostic = mode

  def loadFile(self):
    from core import CCP4File
    f = CCP4File.CI2XmlDataFile(fullPath=self.xmlFile)
    f.loadFile()
    root = f.getEtreeRoot().xpath('./ccp4i2_body')[0]
    return root

  def headerInfo(self,load=False):
    from core import CCP4File
    f = CCP4File.CI2XmlDataFile(fullPath=self.xmlFile)
    f.loadFile()
    f.loadHeader()
    if load:
      self.projectName =  f.header.projectName.__str__()
      self.projectId = f.header.projectId.__str__()
    return f.header.get()

  def loadExportInfo(self):
    try:
      root = self.loadFile()
    except:
      raise CException(self.__class__,120,self.xmlFile)
    eleList = root.xpath('databaseTable/database')
    if len(eleList)!= 1:
      raise CException(self.__class__,122,self.xmlFile)
    dbInfo = {}
    exportInfo = {}
    for key,value in list(eleList[0].items()):
      dbInfo[key] = value
    eleList = root.xpath('projectexportTable/projectexport')
    if len(eleList)!= 1:
      raise CException(self.__class__,122,self.xmlFile,)
    for key,value in list(eleList[0].items()):
      exportInfo[key] = value

    return dbInfo,exportInfo

  def setProjectImport(self):
    dbInfo,exportInfo = self.loadExportInfo()
    #print 'CDbXml.setProjectImport',dbInfo,exportInfo
    if not self.db.idExists(dbInfo['databaseid'],'Databases'):
      loaded = self.db.createDatabase(databaseId=dbInfo['databaseid'],creationTime=dbInfo['creationtime'],
             creationHostName=dbInfo['creationhostname'],createUserName=dbInfo['creationusername'],
                             schemaVersion=dbInfo.get('schemaversion',None),schemaDate=dbInfo.get('schemadate',None) )
      #print 'CDbXml.setProjectImport database already loaded',loaded
    self.db.createProjectImport(projectId=self.projectId,projectExportId=exportInfo['projectexportid'],
       projectExportDatabaseId=dbInfo['databaseid'],projectExportTime=exportInfo['projectexporttime'])

  def loadProjectInfo(self):
    try:
      root = self.loadFile()
    except:
      raise CException(self.__class__, 120, self.xmlFile)

    pEleList = root.xpath('projectTable/project')
    if len(pEleList)!= 1:
      raise CException(self.__class__, 121, self.xmlFile, stack=False)

    self.projectInfo = {}
    for key,value in list(pEleList[0].items()):
      self.projectInfo[key] = value
    self.projectId = self.projectInfo.get('projectid',None)
    self.projectName = self.projectInfo.get('projectname',None)
    self.projectDirectory = self.projectInfo.get('projectdirectory',None)
    del self.projectInfo['projectid']
    del self.projectInfo['projectname']
    del self.projectInfo['projectdirectory']
    return self.projectInfo


  def createProject(self,projectDirectory=None):
    errReport = CErrorReport()
    if projectDirectory is not None: self.projectDirectory = projectDirectory
    print('Attempting to create database project:',self.projectName,'with directory',self.projectDirectory,'and projectId',self.projectId)
    try:
      pid = self.getProjectId(projectName=self.projectName)
      return CErrorReport(self.__class__,2,str(self.projectName))
    except:
      pass
    try:
      pid = self.db.createProject(projectName=self.projectName,projectDirectory=self.projectDirectory,projectId=self.projectId)
    except CException as e:
      return e
    except:
      return CErrorReport(self.__class__,3,str(self.projectName))
    for key,value in list(self.projectInfo.items()):
      try:
        self.db.updateProject(self.projectId,key,value)
      except CException as e:
        errReport.extend(e)
    return CErrorReport()


  def loadJob(self,root,parentJobId=None):
    import time
    try:
      jobId = UUIDTYPE(root.findtext('jobid'))
      jobNumber = root.findtext('jobnumber')
      taskName = root.findtext('taskname')
    except:
      self.errReport.append(self.__class__,5,'Project: '+str(self.projectName))
      return False
    creationTime = root.findtext('creationtime')
    if creationTime is not None:
      creationTime = time.mktime(time.strptime(creationTime,CDbApi.TIMEFORMAT))
    finishTime =root.findtext('finishtime')
    if finishTime is not None:
      finishTime = time.mktime(time.strptime(finishTime,CDbApi.TIMEFORMAT))
    try:
      statusId = JOB_STATUS_TEXT.index(root.findtext('status'))
    except:
      statusId = 0
    try:
      evaluation = JOB_EVALUATION_TEXT.index(root.findtext('evaluation'))
    except:
      evaluation = 0
    userAgent = root.findtext('useragent')
    jobTitle = root.findtext('jobtitle')
    #print 'CDbXml.createJob',oldJobId,taskName

    try:
      # Expect this to throw error because job should not exist -
      # If it does not fail then throw error
      testJobId = self.db.getJobId(projectId=self.projectId,jobNumber=jobNumber)
      self.errReport.append(self.__class__,4,'Project: '+str(self.projectName)+' JobNumber: '+str(jobNumber),stack=False)
      return False
    except:
      pass


    args = [jobId,self.projectId,parentJobId,jobNumber,userAgent,taskName,jobTitle,creationTime,finishTime,statusId,evaluation]
    if self._diagnostic: print('CDbXml.createJob',args)
    self.db.execute("INSERT INTO Jobs (JobID,ProjectID,parentJobId,JobNumber,userAgent,TaskName,JobTitle,CreationTime,FinishTime,Status,Evaluation) VALUES (?,?,?,?,?,?,?,?,?,?,?)",args)

    self.db.commit()

    outputFilesEleList = root.find('outputFiles').getchildren()
    if self._diagnostic: print('Loading',len(outputFilesEleList),'output files')
    for outputFileEle in outputFilesEleList:
      self.createFile(outputFileEle,jobId=jobId)

    inputFilesEleList = root.find('inputFiles').getchildren()
    if self._diagnostic: print('Loading',len(inputFilesEleList),'input files')
    for inputFileEle in inputFilesEleList:
      # Trickier cos some of these should be fileUses not files
      pass


    return True


  def createFile(self,root=None,jobId=None):
    from core import CCP4File
    fileClass = str(root.tag)
    try:
      fileTypeId = FILETYPES_CLASS.index(fileClass)
    except:
      fileTypeId = 0

    '''
    fileId = self.db.uniqueId(table='Files',identifier='FileID')
    try:
      oldFileId = int(root.findtext('dbFileId'))
    except:
      oldFileId = None
    '''
    # !!!!!!!!!!!!!!! This use of importfiles is broken !!!!!!!!!!!!!!
    impRoot = root.find('importFile')
    if impRoot is None:
      importId = None
    else:
      sourceFilename =  impRoot.findtext('sourceFilename')
      annotation = impRoot.findtext('annotation')
      importId = self.db.createImportFile(sourceFileName=sourceFilename,annotation=annotation,jobId=jobId)

    fileId = UUIDTYPE(root.findtext('fileId'))
    baseName = root.findtext('baseName')
    relPath = root.findtext('relPath')
    annotation = root.findtext('annotation')

    args = [fileId,jobId,importId,baseName,relPath,fileTypeId,annotation,]
    if self._diagnostic: print('CDbXml.createFile args',args)
    self.db.execute("INSERT INTO Files (FileID,JobID,ImportId,Filename,RelPath,FiletypeID,Annotation) VALUES (?,?,?,?,?,?,?)",args)

    self.db.commit()

  def deleteFile(self,fileId=None):
    fileId = UUIDTYPE(fileId)
    #print 'CDbApi.deleteFile',fileId
    self.execute("DELETE FROM Files WHERE FileID = ?",(fileId,))
    # ON DELETE CASCADE should deal with these - also FileUse and ExportFile
    self.execute("DELETE FROM ImportFiles WHERE FileID = ?",(fileId,))
    self.execute("DELETE FROM ExportFiles WHERE FileID = ?",(fileId,))
    self.execute("DELETE FROM FileUses WHERE FileID = ?",(fileId,))
    self.fileDeleted.emit(fileId)

  #def createFileUse(self,root=None,jobId=None):    # KJS : This function is a mess. Need to revise later.
  #  args = [fileId,jobId,role,]
  #  if self._diagnostic: print 'CDbApi.createFileUse',args
  #  self.execute("INSERT INTO FileUses (FileID,JobID,RoleID) VALUES (?,?,?)",args)
  #  self.commit()

  def makeAllSql(self):
    self.sql = {}
    self.makeSql('Jobs',self.db.JOBITEMS)
    self.makeSql('Files',self.db.FILEITEMS)
    self.makeSql('ImportFiles',self.db.IMPORTFILEITEMS)
    self.makeSql('ExportFiles',self.db.EXPORTFILEITEMS)
    self.makeSql('FileUses',self.db.FILEUSEITEMS)
    self.makeSql('XData',self.db.XDATAITEMS)
    self.makeSql('Comments',self.db.COMMENTITEMS)
    self.makeSql('JobKeyValues',self.db.JOBKEYVALUEITEMS)
    self.makeSql('JobKeyCharValues',self.db.JOBKEYCHARVALUEITEMS)
    self.makeSql('FileAssociations',self.db.FILEASSOCIATIONITEMS)
    self.makeSql('FileAssociationMembers',self.db.FILEASSOCIATIONMEMBERITEMS)

    self.makeSql('TempJobs',self.db.JOBITEMS)
    self.makeSql('TempFiles',self.db.FILEITEMS)
    self.makeSql('TempImportFiles',self.db.IMPORTFILEITEMS)
    self.makeSql('TempExportFiles',self.db.EXPORTFILEITEMS)
    self.makeSql('TempFileUses',self.db.FILEUSEITEMS)
    self.makeSql('TempXData',self.db.XDATAITEMS)
    self.makeSql('TempComments',self.db.COMMENTITEMS)
    self.makeSql('TempJobKeyValues',self.db.JOBKEYVALUEITEMS)
    self.makeSql('TempJobKeyCharValues',self.db.JOBKEYCHARVALUEITEMS)
    self.makeSql('TempFileAssociations',self.db.FILEASSOCIATIONITEMS)
    self.makeSql('TempFileAssociationMembers',self.db.FILEASSOCIATIONMEMBERITEMS)
    self.makeSql('TempTags',self.db.TAGITEMS)
    self.makeSql('TempProjectTags',self.db.PROJECTTAGITEMS)
    self.makeSql('TempProjectComments',self.db.PROJECTCOMMENTITEMS)

    self.makeSql('ImportProjects',self.db.PROJECTIMPORTITEMS)

    #print 'CDbXml',self.sql


  def makeSql(self,tableName,attributeList):
    com = "INSERT INTO "+tableName + " (" + attributeList[0]
    val = '?'
    for item in attributeList[1:]:
      com = com + ','+item
      val = val + ',?'
    self.sql[tableName] = com + ")  VALUES ("+val+")"


  def loadJobsFromTable(self,commitPolicy=None):
    if commitPolicy is None: commitPolicy = CDbXml.COMMIT_POLICY_NO_ERRORS
    from lxml import etree
    root = self.loadFile()
    if root is None:
      self.errReport.append(self.__class__,1,str(self.xmlFile))
      return False
    #print 'loadTable root',root.tag

    if commitPolicy == CDbXml.COMMIT_POLICY_NO_COMMIT:
      return False
    elif commitPolicy == CDbXml.COMMIT_POLICY_NO_ERRORS and  self.errReport.maxSeverity()>SEVERITY_WARNING:
      return False
    else:
      self.db.commit()
      return True


  def loadTable(self,commitPolicy=None,newJobNumber=None,selectJobIdList=None,updateJobStatus=False):
    if commitPolicy is None: commitPolicy = CDbXml.COMMIT_POLICY_NO_ERRORS
    from lxml import etree
    self.loadedJobs = {}
    root = self.loadFile()
    if root is None:
      self.errReport.append(self.__class__,1,str(self.xmlFile))
      return False
    #self.db.setDiagnostic(True)

    # Jobs
    try:
      table = root.xpath('./jobTable')[0]
    except:
      self.errReport.append(self.__class__,200,'jobTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.JOBITEMS,self.db.JOBTYPES)
        jobId = rowData[self.db.JOBITEMS.index('jobid')]
        loadedJob = None
        if selectJobIdList is None or jobId in selectJobIdList:
          if self.db.idExists(jobId,'Jobs'):
            if not updateJobStatus:
              self.errReport.append(self.__class__,10,'JobId '+str(jobId))
            else:
              # Assume this is from a remotely run job - update status
              for key in ['status','finishtime']:
                value = rowData[self.db.JOBITEMS.index(key)]
                self.db.execute("UPDATE Jobs SET "+key+" = ? WHERE JobID = ?",(value,jobId,))
              loadedJob=(rowData[self.db.JOBITEMS.index('jobnumber')],None)
          else:
            # Ensure the job number is unique - we could be merging diverged projects
            try:
              # getJobId will raise error if there is no job with this jobNumber - thats good
              oldJobNumber = rowData[self.db.JOBITEMS.index('jobnumber')]
              parentJobId =  rowData[self.db.JOBITEMS.index('parentjobid')]
              if newJobNumber is None:
                self.loadedJobs[jobId]=(oldJobNumber,None)
              elif parentJobId is None:
                # use next number if it is top level job
                rowData[self.db.JOBITEMS.index('jobnumber')] = str(newJobNumber)
                self.errReport.append(self.__class__,9,'Imported number '+ str(oldJobNumber) + ' changed to '+ str(newJobNumber))
                self.loadedJobs[jobId]=(oldJobNumber,newJobNumber)
                newJobNumber+=1
              else:
                if parentJobId in self.loadedJobs:
                  subJobNumber = str(self.loadedJobs[parentJobId][1]) + rowData[self.db.JOBITEMS.index('jobnumber')].split('.')[-1]
                  rowData[self.db.JOBITEMS.index('jobnumber')] = subJobNumber
                  loadedJob=(oldJobNumber,subJobNumber)
                else:
                  loadedJob=(oldJobNumber,None)
            except:
              print('ERROR in loadTable processing jobId',jobId)
              self.errReport.append(self.__class__,123,str(jobId))
            if self._diagnostic: print('CDbXml create Job',rowData)
            try:
              self.db.execute(self.sql['Jobs'],rowData)
            except Exception as e:
              self.errReport.append(self.__class__,210,str(rowData)+' '+str(e))
            else:
              if loadedJob is not None: self.loadedJobs[jobId] = loadedJob

    #print 'CDbXml.loadTable loadedJobs',self.loadedJobs
    if len(self.loadedJobs)>0:
      # Update the project lastJobNumber
      if newJobNumber is not None:
        lastJobNumber = newJobNumber
      else:
        lastJobNumber = oldJobNumber.split('.')[0]
      #print 'loadTable Updating LastJobNumber',lastJobNumber
      self.db.execute('UPDATE Projects SET LastJobNumber=? WHERE ProjectId = ?',(int(lastJobNumber),self.projectId))

    # Files
    loadedFiles = []
    try:
      table = root.xpath('./fileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEITEMS,self.db.FILETYPES)
        # test that job has been imported
        idx = self.db.FILEITEMS.index('jobid')
        #print 'CDbXml.loadTable attemting file',self.loadedJobs.has_key(rowData[idx]), self.db.idExists(rowData[idx],'Jobs','JobId')
        if (rowData[idx] not in self.loadedJobs) and (not self.db.idExists(rowData[idx],'Jobs','JobId')):
          self.errReport.append(self.__class__,216,idx)
        else:
          if self._diagnostic: print('CDbXml create File',rowData)
          try:
            self.db.execute(self.sql['Files'],rowData)
            loadedFiles.append(rowData[0])
          except:
            self.errReport.append(self.__class__,211,self.sql['Files']+' '+str(rowData))

    # ImportFiles
    # There is issue that file may be recorded with an importID that is not yet valid
    try:
      table = root.xpath('./importfileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'importfileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.IMPORTFILEITEMS,self.db.IMPORTFILETYPES)
        #Test if we have the corresponding row in the file table
        idx = self.db.IMPORTFILEITEMS.index('fileid')
        if (not rowData[idx] in loadedFiles) and (not self.db.idExists(rowData[idx],'Files','FileId')):
          self.errReport.append(self.__class__,7,str(rowData))
        else:
          try:
            self.db.execute(self.sql['ImportFiles'],rowData)
          except Exception as e:
            self.errReport.append(self.__class__,212,str(rowData)+' '+str(e))

    # ExportFiles
    try:
      table = root.xpath('./exportfileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'exportfileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.EXPORTFILEITEMS,self.db.EXPORTFILETYPES)
        # Test if fileId of exported file is in the database
        idx = self.db.EXPORTFILEITEMS.index('fileid')
        if (not not rowData[idx] in loadedFiles) and (not self.db.idExists(rowData[idx],'Files','FileId')):
          self.errReport.append(self.__class__,8,str(rowData))
        else:
          if self._diagnostic: print('CDbXml create ExportFile',rowData)
          self.db.execute(self.sql['ExportFiles'],rowData)

    # FileUses
    try:
      table = root.xpath('./fileuseTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileuseTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEUSEITEMS,self.db.FILEUSETYPES)
        # Test that job has been imported and that FileId exisits
        idx = self.db.FILEUSEITEMS.index('jobid')
        if (rowData[idx] not in self.loadedJobs) and (not self.db.idExists(rowData[idx],'Jobs','JobId')):
          self.errReport.append(self.__class__,217,'FileUses '+str(rowData))
        else:
          idx = self.db.FILEUSEITEMS.index('fileid')
          if (not rowData[idx] in loadedFiles)  and  (not self.db.idExists(rowData[idx],'Files')):
            self.errReport.append(self.__class__,218,'FileUses '+str(rowData))
          else:
            try:
              self.db.execute(self.sql['FileUses'],rowData)
            except:
              self.errReport.append(self.__class__,213,'FileUses '+str(rowData))

    # XData
    try:
      table = root.xpath('./xdataTable')[0]
    except:
      self.errReport.append(self.__class__,200,'xdataTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.XDATAITEMS[0:-1],self.db.XDATATYPES[0:-1])
        #print 'CImportDb.loadTable xdata rowData',rowData
        idx = self.db.XDATAITEMS.index('jobid')
        if (rowData[idx] not in self.loadedJobs) and (not self.db.idExists(rowData[idx],'Jobs','JobId')):
          self.errReport.append(self.__class__,219,'XData: '+str(rowData))
        else:
          try:
            xdataEtree = rowTree.xpath('./xdataxml')[0]
            cls = rowData[self.db.XDATAITEMS.index('xdataclass')]
            xdataEtree.tag = rowData[self.db.XDATAITEMS.index('xdataclass')]
            xdataString = etree.tostring(xdataEtree)
          except Exception as e:
            self.errReport.append(self.__class__,205,str(rowData)+' '+str(e))
          else:
            rowData.append(xdataString)
            if self._diagnostic: print('CDbXml create XData',rowData)
            try:
              self.db.execute(self.sql['XData'],rowData)
            except:
              self.errReport.append(self.__class__,213,'XData '+str(rowData))

    # Comments
    try:
      table = root.xpath('./commentTable')[0]
    except:
      self.errReport.append(self.__class__,200,'commentTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.COMMENTITEMS,self.db.COMMENTTYPES)
        if self.db.idExists(rowData[self.db.COMMENTITEMS.index('commentid')],'Comments'):
          pass
        else:
          idx = self.db.COMMENTITEMS.index('jobid')
          if (rowData[idx] not in self.loadedJobs) and (not self.db.idExists(rowData[idx],'Jobs','JobId')):
            # job not recognised in this database !!!!
            self.errReport.append(self.__class__,215,'Comment '+str(rowData))
          else:
            if not self.db.idExists(rowData[self.db.COMMENTITEMS.index('userid')],'Users'):
              # user not recognised in this database !!!!
              rowData[self.db.COMMENTITEMS.index('userid')] = None
              self.errReport.append(self.__class__,214,'Comment '+str(rowData))
            if self._diagnostic: print('CDbXml.create comment',rowData)
            try:
              self.db.execute(self.sql['Comments'],rowData)
            except:
              self.errReport.append(self.__class__,213,'Comment '+str(rowData))

    # ProjectComments
    try:
      table = root.xpath('./projectcommentTable')[0]
    except:
      self.errReport.append(self.__class__,200,'projectcommentTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.PROJECTCOMMENTITEMS,self.db.PROJECTCOMMENTTYPES)
        if self.db.idExists(rowData[self.db.PROJECTCOMMENTITEMS.index('projectcommentid')],'ProjectComments'):
          pass
        else:
          idx = self.db.PROJECTCOMMENTITEMS.index('projectid')
          if (not self.projectId == rowData[idx]) and (not self.db.idExists(rowData[idx],'Projects','ProjectId')):
            # job not recognised in this database !!!!
            self.errReport.append(self.__class__,215,'ProjectComment '+str(rowData))
          else:
            if not self.db.idExists(rowData[self.db.PROJECTCOMMENTITEMS.index('userid')],'Users'):
              # user not recognised in this database !!!!
              rowData[self.db.PROJECTCOMMENTITEMS.index('userid')] = None
              self.errReport.append(self.__class__,214,'ProjectComment '+str(rowData))
            if self._diagnostic: print('CDbXml.create projectcomment',rowData)
            try:
              self.db.execute(self.sql['ProjectComments'],rowData)
            except:
              self.errReport.append(self.__class__,213,'ProjectComment '+str(rowData))

    # JobKeyValue
    for tableName,items,types,label in [ [ 'jobkeyvalueTable',self.db.JOBKEYVALUEITEMS,self.db.JOBKEYVALUETYPES,'JobKeyValues' ],
                                     [ 'jobkeycharvalueTable',self.db.JOBKEYCHARVALUEITEMS,self.db.JOBKEYCHARVALUETYPES,'JobKeyCharValues' ] ]:
      try:
        table = root.xpath('./'+tableName)[0]
      except:
        self.errReport.append(self.__class__,200,tableName)
      else:
        for rowTree in table.getchildren():
          rowData = self.getRowFromEtree(rowTree,items,types)
          if 0:   # Need test if item already in db
            pass
          else:
            idx = items.index('jobid')
            if (rowData[idx] not in self.loadedJobs) and (not self.db.idExists(rowData[idx],'Jobs','JobId')):
              # job not recognised in this database !!!!
              self.errReport.append(self.__class__,215,label+' '+str(rowData))
            else:
              try:
                self.db.execute(self.sql[label],rowData)
              except:
                self.errReport.append(self.__class__,213,label+' '+str(rowData))
   # FileAssociations
    try:
      table = root.xpath('./fileassociationTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileassociationTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEASSOCIATIONITEMS,self.db.FILEASSOCIATIONTYPES)
        if self.db.idExists(rowData[self.db.FILEASSOCIATIONITEMS.index('fileassociationid')],'FileAssociations'):
          pass
        else:
          try:
            self.db.execute(self.sql['FileAssociations'],rowData)
          except:
            self.errReport.append(self.__class__,213,'FileAssociations '+str(rowData))

   # FileAssociationMembers
    try:
      table = root.xpath('./fileassociationmemberTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileassociationmemberTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEASSOCIATIONMEMBERITEMS,self.db.FILEASSOCIATIONMEMBERTYPES)
        if 0:  # need to test if alrady loaded
          pass
        else:
          try:
            self.db.execute(self.sql['FileAssociationMembers'],rowData)
          except:
            self.errReport.append(self.__class__,213,'FileAssociationMembers '+str(rowData))


    self.db.setDiagnostic(False)
    if commitPolicy == CDbXml.COMMIT_POLICY_NO_COMMIT:
      return False
    elif commitPolicy == CDbXml.COMMIT_POLICY_NO_ERRORS and  self.errReport.maxSeverity()>SEVERITY_WARNING:
      return False
    else:
      self.db.commit()
      self.db.projectReset.emit({'projectId' : self.projectId})
      return True


  def removeTempTables(self):
    self.db.removeTempTables()

  def createTempTables(self,fileName=None):
    self.db.removeTempTables()
    # Create temporary tables - beware may exist from previous run
    # Should have a mechanism to block two processes trying to do this at once!
    if fileName is  None:
        from core import CCP4Utils
        fileName = os.path.join(CCP4Utils.getCCP4I2Dir(),'dbapi','temp_database_schema.sql')
    #print 'Loading temp database schema from:',fileName
    self.db.read(fileName=fileName)
    #print 'Finished loading temp database schema.'



  def loadTempTable(self,commitPolicy=None,resetProjectId=None):
    if commitPolicy is None: commitPolicy = CDbXml.COMMIT_POLICY_NO_ERRORS
    from lxml import etree
    self.loadedJobs = {}
    root = self.loadFile()
    if root is None:
      self.errReport.append(self.__class__,1,str(self.xmlFile))
      return False
    #self.db.setDiagnostic(True)

    # Jobs
    try:
      table = root.xpath('./jobTable')[0]
    except:
      self.errReport.append(self.__class__,200,'jobTable')
    else:
      doneLoaded = []
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.JOBITEMS,self.db.JOBTYPES)
        if resetProjectId is not None: rowData[self.db.JOBITEMS.index('projectid')] = resetProjectId
        jobId = rowData[self.db.JOBITEMS.index('jobid')]
        #print 'loadTempTable',jobId,doneLoaded
        if not jobId in doneLoaded:
          try:
            self.db.execute(self.sql['TempJobs'],rowData)
          except Exception as e:
            print('loadTempTable:Error loading job:',e)
            self.errReport.append(self.__class__,210,str(rowData)+' '+str(e))
          else:
            doneLoaded.append(jobId)


    # Files
    try:
      table = root.xpath('./fileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEITEMS,self.db.FILETYPES)
        try:
          self.db.execute(self.sql['TempFiles'],rowData)
        except:
          self.errReport.append(self.__class__,211,str(rowData))

    # ImportFiles
    # There is issue that file may be recorded with an importID that is not yet valid
    try:
      table = root.xpath('./importfileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'importfileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.IMPORTFILEITEMS,self.db.IMPORTFILETYPES)
        try:
          self.db.execute(self.sql['TempImportFiles'],rowData)
        except Exception as e:
          self.errReport.append(self.__class__,212,str(rowData)+' '+str(e))

    # ExportFiles
    try:
      table = root.xpath('./exportfileTable')[0]
    except:
      self.errReport.append(self.__class__,200,'exportfileTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.EXPORTFILEITEMS,self.db.EXPORTFILETYPES)
        try:
          self.db.execute(self.sql['TempExportFiles'],rowData)
        except Exception as e:
          self.errReport.append(self.__class__,212,str(rowData)+' '+str(e))

    # FileUses
    try:
      table = root.xpath('./fileuseTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileuseTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEUSEITEMS,self.db.FILEUSETYPES)
        try:
          self.db.execute(self.sql['TempFileUses'],rowData)
        except:
          self.errReport.append(self.__class__,213,'FileUses '+str(rowData))

    # XData
    try:
      table = root.xpath('./xdataTable')[0]
    except:
      self.errReport.append(self.__class__,200,'xdataTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.XDATAITEMS[0:-1],self.db.XDATATYPES[0:-1])
        #print 'CImportDb.loadTable xdata rowData',rowData
        try:
          xdataEtree = rowTree.xpath('./xdataxml')[0]
          cls = rowData[self.db.XDATAITEMS.index('xdataclass')]
          xdataEtree.tag = rowData[self.db.XDATAITEMS.index('xdataclass')]
          xdataString = etree.tostring(xdataEtree)
        except Exception as e:
          self.errReport.append(self.__class__,205,str(rowData)+' '+str(e))
        else:
          rowData.append(xdataString)
          if self._diagnostic: print('CDbXml create XData',rowData)
          try:
            self.db.execute(self.sql['TempXData'],rowData)
          except:
            self.errReport.append(self.__class__,213,'XData '+str(rowData))

    # JobKeyValues
    for tableName,tempTableName,items,types in [ [ 'jobkeyvalueTable','TempJobKeyValues',self.db.JOBKEYVALUEITEMS,self.db.JOBKEYVALUETYPES],
                                                 [ 'jobkeycharvalueTable','TempJobKeyCharValues',self.db.JOBKEYCHARVALUEITEMS,self.db.JOBKEYCHARVALUETYPES] ]:
      try:
        table = root.xpath('./'+tableName)[0]
      except:
        self.errReport.append(self.__class__,200,tableName)
      else:
        for rowTree in table.getchildren():
          rowData = self.getRowFromEtree(rowTree,items,types)
          if self._diagnostic: print('CDbXml create JobKeyValue',rowData,self.sql[tempTableName])
          try:
            self.db.execute(self.sql[tempTableName],rowData)
          except:
            self.errReport.append(self.__class__,213,tableName+' '+str(rowData))


    # Comments
    try:
      table = root.xpath('./commentTable')[0]
    except:
      self.errReport.append(self.__class__,200,'commentTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.COMMENTITEMS,self.db.COMMENTTYPES)
        # user not recognised in this database !!!!
        if not self.db.idExists(rowData[self.db.COMMENTITEMS.index('userid')],'Users'):
            self.db.execute('SELECT * FROM Users',())
            userlist = self.db.fetchAll2Py(str)
            if len(userlist)>0:
                rowData[self.db.COMMENTITEMS.index('userid')] = userlist[0]
            else:
                continue
        if self.db.idExists(rowData[self.db.COMMENTITEMS.index('commentid')],'Comments'):
          pass
        else:
          try:
            self.db.execute(self.sql['TempComments'],rowData)
          except:
            self.errReport.append(self.__class__,213,'Comment '+str(rowData))

    # ProjectComments
    try:
      table = root.xpath('./projectcommentTable')[0]
    except:
      self.errReport.append(self.__class__,200,'projectcommentTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.PROJECTCOMMENTITEMS,self.db.PROJECTCOMMENTTYPES)
        # user not recognised in this database !!!!
        if not self.db.idExists(rowData[self.db.PROJECTCOMMENTITEMS.index('userid')],'Users'):
            self.db.execute('SELECT * FROM Users',())
            userlist = self.db.fetchAll2Py(str)
            if len(userlist)>0:
                rowData[self.db.PROJECTCOMMENTITEMS.index('userid')] = userlist[0]
            else:
                continue
        if self.db.idExists(rowData[self.db.PROJECTCOMMENTITEMS.index('projectcommentid')],'ProjectComments'):
          pass
        else:
          try:
            self.db.execute(self.sql['TempProjectComments'],rowData)
          except:
            self.errReport.append(self.__class__,213,'ProjectComment '+str(rowData))


    # FileAssociations
    try:
      table = root.xpath('./fileassociationTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileassociationTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEASSOCIATIONITEMS,self.db.FILEASSOCIATIONTYPES)
        if self._diagnostic: print('CDbXml create FileAssociations')
        try:
          self.db.execute(self.sql['TempFileAssociations'],rowData)
        except:
          self.errReport.append(self.__class__,213,'FileAssociation '+str(rowData))

    # FileAssociationmembers
    try:
      table = root.xpath('./fileassociationmemberTable')[0]
    except:
      self.errReport.append(self.__class__,200,'fileassociationmemberTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.FILEASSOCIATIONMEMBERITEMS,self.db.FILEASSOCIATIONMEMBERTYPES)
        if self._diagnostic: print('CDbXml create FileAssociationMembers')
        try:
          self.db.execute(self.sql['TempFileAssociationMembers'],rowData)
        except:
          self.errReport.append(self.__class__,213,'FileAssociationMember '+str(rowData))


    # Tags
    try:
      table = root.xpath('./tagTable')[0]
    except:
      self.errReport.append(self.__class__,200,'tagTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.TAGITEMS,self.db.TAGTYPES)
        if self.db.idExists(rowData[self.db.TAGITEMS.index('tagid')],'Tags'):
          pass
        else:
          try:
            self.db.execute(self.sql['Tags'],rowData)
          except:
            self.errReport.append(self.__class__,213,'Tags '+str(rowData))
    try:
      table = root.xpath('./projecttagTable')[0]
    except:
      self.errReport.append(self.__class__,200,'projecttagTable')
    else:
      for rowTree in table.getchildren():
        rowData = self.getRowFromEtree(rowTree,self.db.PROJECTTAGITEMS,self.db.PROJECTTAGTYPES)

        if self.db.projectTagExists(rowData[self.db.PROJECTTAGITEMS.index('projectid')],rowData[self.db.PROJECTTAGITEMS.index('tagid')]):
          pass
        else:
          try:
            self.db.execute(self.sql['ProjectTags'],rowData)
          except:
            self.errReport.append(self.__class__,213,'ProjectTags '+str(rowData))



    self.db.setDiagnostic(False)
    if commitPolicy == CDbXml.COMMIT_POLICY_NO_COMMIT:
      return False
    elif commitPolicy == CDbXml.COMMIT_POLICY_NO_ERRORS and  self.errReport.maxSeverity()>SEVERITY_WARNING:
      return False
    else:
      self.db.commit()
      return True

  def setExclInTempTables(self):
    # Exclude jobs that are already in db
    # TempJobs.Excl= 2 Job in db finishd and no further changes allowed - exclude from copying
    # TempJobs.Excl= 1 Job in db - unfinished - allow update status and copying files
    # TempJobs.Excl= 0 Job not in db - allow all copying
    # TempJobs.Excl= -1 - job has been copied in project dir so can be copied  into db
    if self._diagnostic:
      self.db.execute('SELECT JobId FROM Jobs WHERE ProjectId=? AND Status IN (?,?,?,?)',[self.projectId,JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_UNSATISFACTORY,JOB_STATUS_FILE_HOLDER])
      jobIdList = self.db.fetchAll2Py(str)
      print('setExclInTempTables TempJobs finished/failed excl=2',jobIdList)
      self.db.execute('SELECT JobId FROM Jobs WHERE ProjectId=? AND Status NOT IN (?,?,?,?)',[self.projectId,JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_UNSATISFACTORY,JOB_STATUS_FILE_HOLDER])
      jobIdList = self.db.fetchAll2Py(str)
      print('setExclInTempTables TempJobs not finished/failed excl=1',jobIdList)

    self.db.execute('UPDATE TempJobs SET Excl=2 WHERE JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=? AND Status IN (?,?,?,?))',[self.projectId,JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_UNSATISFACTORY,JOB_STATUS_FILE_HOLDER])

    self.db.execute('UPDATE TempJobs SET Excl=1 WHERE JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=? AND Status NOT IN (?,?,?,?))',[self.projectId,JOB_STATUS_FINISHED,JOB_STATUS_FAILED,JOB_STATUS_UNSATISFACTORY,JOB_STATUS_FILE_HOLDER])

    self.db.commit()

    # Do we need to revise job numbers
    self.db.execute('SELECT COUNT(*) FROM TempJobs WHERE Excl=0 and parentJobId IS NULL')
    nImportJobs = self.db.fetchAll2Py(int)[0]
    if self._diagnostic: print('setExclInTempTables nImportJobs',nImportJobs)
    if nImportJobs>0:
      self.db.execute('SELECT MIN(CAST(JobNumber AS integer)) FROM TempJobs WHERE Excl=0 and parentJobId IS NULL' )
      importMin = self.db.fetchAll2Py(int)[0]
      self.db.execute('SELECT MAX(CAST(JobNumber AS integer)) FROM TempJobs WHERE Excl=0 and parentJobId IS NULL' )
      importMax = self.db.fetchAll2Py(int)[0]

      self.db.execute('SELECT COUNT(*) FROM Jobs WHERE parentJobId IS NULL AND ProjectId = ?',(self.projectId,))
      nCurrentJobs = self.db.fetchAll2Py(int)[0]
      if nCurrentJobs>0:
        self.db.execute('SELECT MAX(CAST(JobNumber AS integer)) FROM Jobs WHERE parentJobId IS NULL AND ProjectId = ?',(self.projectId,) )
        currentMax = self.db.fetchAll2Py(int)[0]
      else:
        currentMax = 0
      #print 'excludeDuplicates job number nImportJobs,nCurrentJobs,importMin,importMax,currentMax',nImportJobs,nCurrentJobs,importMin,importMax,currentMax

      self.stats['importMin'] = importMin
      self.stats['importMax'] = importMax
      if importMin<=currentMax:
        incrJobNumber = currentMax - importMin + 1
        self.stats['incrJobNumber'] = incrJobNumber
        # Custom function to reset job number
        # See http://stackoverflow.com/questions/2108870/how-to-create-custom-functions-in-sqlite
        def updateJobNumber(number):
          #print 'updateJobNumber',number,incrJobNumber,
          jList = number.split('.')
          ret = str(int(jList[0])+incrJobNumber)
          for n in jList[1:]: ret = ret + '.' + n
          #print ret
          return ret

        self.db.connection().create_function('updateJobNumber',1,updateJobNumber)
        self.db.execute('UPDATE TempJobs SET NewJobNumber=updateJobNumber(JobNumber) WHERE Excl<=0')



    # Delete duplicate files/importfiles/fileuses
    # Bother - sqlite does not do inner joins in update and the alternative 'JobId IN' business is probably not v. efficient
    if self._diagnostic:
      self.listTempJobs('Jobs after setExclInTempTables')
      self.db.execute('SELECT FileId,FileName FROM TempFiles  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl = 2)')
      fileList = self.db.fetchAll2PyList([str,str])
      print('setExclInTempTables TempFiles from Jobs not finished/failed excl=2',fileList)
      self.db.execute('SELECT FileId,FileName FROM TempFiles WHERE FileId IN (SELECT FileId FROM Files)')
      fileList = self.db.fetchAll2PyList([str,str])
      print('setExclInTempTables TempFiles fileId already in Files excl=1',fileList)
    self.db.execute('UPDATE TempFiles SET Excl=1  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl = 2) OR FileId IN (SELECT FileId FROM Files)')

  def setExclImportedFiles(self):
    # mark the imported files as 'NOT excluded'
    if self._diagnostic:
      self.db.execute('SELECT FileId,FileName FROM TempFiles WHERE Excl=0 AND FileId IN (SELECT FileId FROM TempImportFiles)')
      fileList = self.db.fetchAll2PyList([str,str])
      #print 'setExclImportedFiles imported files',fileList
    self.db.execute('UPDATE TempFiles SET Excl=-1  WHERE Excl=0 AND FileId IN (SELECT FileId FROM TempImportFiles)')

  def cleanupTempTables(self):
    # Exclude TempImportFiles that are already recorded in db or that refer to a TempFile that has not been imported
    self.db.execute('UPDATE TempImportFiles SET Excl=1  WHERE FileId IN (SELECT FileId FROM TempFiles WHERE Excl != -1) OR ImportId IN (SELECT ImportFiles.ImportId FROM ImportFiles INNER JOIN Files ON ImportFiles.FileId =Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE JOBS.ProjectId = ?)',(self.projectId,))
    #Exclude TempFileUses that refer to a TempJob that has not been imported
    self.db.execute('UPDATE TempFileUses SET Excl=1  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl > 0)')
    #Exclude TempXData that refer to a TempJob that has not been imported or that are already in Db
    self.db.execute('UPDATE TempXData SET Excl=1  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl > 0) OR XDataId IN (SELECT XData.XDataId FROM XData INNER JOIN Jobs ON XData.JobId=Jobs.JobId WHERE Jobs.ProjectId=? )',(self.projectId,))
    #Exclude TempJobKeyValues that refer to a TempJob that has not been imported or that are already in Db
    #self.db.execute('UPDATE TempJobKeyValues SET Excl=1  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl > 0) OR JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=? )',(self.projectId,))
    #self.db.execute('UPDATE TempJobKeyCharValues SET Excl=1  WHERE JobId IN (SELECT JobId FROM TempJobs WHERE Excl > 0) OR JobId IN (SELECT JobId FROM Jobs WHERE ProjectId=? )',(self.projectId,))
    #Exclude TempExportFiles and TempComments that are already in Db
    self.db.execute('UPDATE TempExportFiles SET Excl=1 WHERE ExportId IN (SELECT ExportFiles.ExportId FROM ExportFiles INNER JOIN Files ON ExportFiles.FileId =Files.FileId INNER JOIN Jobs ON Files.JobId=Jobs.JobId WHERE JOBS.ProjectId = ?)',(self.projectId,))
    self.db.execute('UPDATE TempComments SET Excl=1 WHERE CommentId IN (SELECT Comments.CommentId FROM Comments INNER JOIN Jobs ON Comments.JobId=Jobs.JobId WHERE Jobs.ProjectId=? )',(self.projectId,))

    self.db.commit()

  def importStats(self):
    stats = {}
    self.db.execute('SELECT COUNT (*) FROM TempJobs')
    stats['jobsTotal'] = self.db.fetchAll2Py(int)[0]
    self.db.execute('SELECT COUNT (*) FROM TempJobs WHERE Excl=-1')
    stats['jobsImported'] = self.db.fetchAll2Py(int)[0]
    self.db.execute('SELECT COUNT (*) FROM TempFiles')
    stats['filesTotal'] = self.db.fetchAll2Py(int)[0]
    self.db.execute('SELECT COUNT (*) FROM TempFiles WHERE Excl=-1')
    stats['filesImported'] = self.db.fetchAll2Py(int)[0]
    self.db.execute('SELECT TempFiles.FileId,TempFiles.JobId,TempFiles.FileName,TempFiles.Excl,TempJobs.JobNumber from TempFiles INNER JOIN TempJobs ON TempFiles.JobId=TempJobs.JobId  WHERE TempFiles.Excl!=-1')
    stats['failedFiles'] = self.db.fetchAll2PyList([str,str,str,int,str])
    stats.update( self.stats )
    return stats

  def importProjectCommentsTempTables(self):
    """
    SJM 2/6/2020 - Apologies for the existence of this method. My sqlite3 skills are too weak to risk breaking an existing function. But I think this is sane.
    """

    self.db.execute('PRAGMA foreign_keys = OFF')

    tableList = ['ProjectComments',self.db.PROJECTCOMMENTITEMS,self.db.PROJECTCOMMENTTYPES,0]

    name,columns,typeList,excl = tableList
    columnText = columns[0]
    for item in columns[1:]:columnText  = columnText  + ','+item
    com = 'INSERT OR REPLACE INTO '+name+' ('+columnText + ')'
    selectCom =  ' SELECT ' +columnText+' FROM Temp'+name
    self.db.execute(com+selectCom)
    self.db.commit()
    self.db.execute('PRAGMA foreign_keys = ON')

  def importTempTables(self,testExl=True,includeCode=-1,ignoreList=[]):
    """
    SJM 2/6/2020 - Added ignoreList, because things with userids (Comments, ProjectComments) currently cause failure.
    """
    if self._diagnostic:
      self.db.execute('SELECT JobId,JobNumber,Excl FROM TempJobs')
      jobList = self.db.fetchAll2PyList([str,str,int])
      for job in  jobList: print('importTempTables job',job)

    self.db.execute('PRAGMA foreign_keys = OFF')

    tableList = [ ['Jobs',self.db.JOBITEMS,self.db.JOBTYPES,includeCode],
                  ['JobKeyValues',self.db.JOBKEYVALUEITEMS,self.db.JOBKEYVALUETYPES,0],
                  ['JobKeyCharValues',self.db.JOBKEYCHARVALUEITEMS,self.db.JOBKEYCHARVALUETYPES,0],
                  ['Files',self.db.FILEITEMS,self.db.FILETYPES,includeCode],
                  ['FileUses',self.db.FILEUSEITEMS,self.db.FILEUSETYPES,0],
                  ['ImportFiles',self.db.IMPORTFILEITEMS,self.db.IMPORTFILETYPES,0],
                  ['ExportFiles',self.db.EXPORTFILEITEMS,self.db.EXPORTFILETYPES,0],
                  ['Comments',self.db.COMMENTITEMS,self.db.COMMENTTYPES,0],
                  ['XData',self.db.XDATAITEMS,self.db.XDATATYPES,0],
                  ['FileAssociations',self.db.FILEASSOCIATIONITEMS,self.db.FILEASSOCIATIONTYPES,0],
                  ['FileAssociationMembers',self.db.FILEASSOCIATIONMEMBERITEMS,self.db.FILEASSOCIATIONMEMBERTYPES,0]
                  ]
    for name,columns,typeList,excl in tableList:
        if name in ignoreList: continue
        columnText = columns[0]
        for item in columns[1:]:columnText  = columnText  + ','+item
        if name == 'Jobs' and self.stats['incrJobNumber']>0 :
          # Copy jobNumber from TempJobs.NewJobNumber
          com = 'INSERT OR IGNORE INTO '+name+' ('+columnText + ')'
          selectCom = ' SELECT '
          selectCom = selectCom + re.sub('jobnumber','newjobnumber',columnText) +' FROM Temp'+name
        else:
          # Copy everything from Tempx to x
          if name in ['JobKeyValues','JobKeyCharValues','Comments']:
            com = 'INSERT OR REPLACE INTO '+name+' ('+columnText + ')'
          else:
            com = 'INSERT OR IGNORE INTO '+name+' ('+columnText + ')'
          selectCom =  ' SELECT ' +columnText+' FROM Temp'+name
        if testExl: selectCom = selectCom + ' WHERE Excl='+str(excl)
        if self._diagnostic:
          valueCom = 'VALUES('
          for item in typeList: valueCom = valueCom + '?,'
          valueCom = valueCom[0:-1] + ')'
          self.db.execute(selectCom)
          rv = self.db.fetchAll2PyList(typeList)
          for item in rv:
            print('importTempTables',name,item)
            self.db.execute(com+valueCom,item)
        else:
          self.db.execute(com+selectCom)
        if name in ['Jobs','Files']: self.db.commit()
    self.db.commit()


    # if TempJob.Excl=-1 or -2 then the job is in db but may need updating
    # sqlite wil not do INNER JOIN in UPDATE
    # Just loop in Python and trust there are not too many jobs with changed status
    com = 'SELECT '+self.db.JOBITEMS[0]
    for item in self.db.JOBITEMS[1:]: com = com + ',' + item
    com = com + ' FROM TempJobs WHERE Excl IN (-1,1,-2)'
    #print 'importTempTables',com
    self.db.execute(com)
    ret = self.db.fetchAll2PyList(self.db.JOBTYPES)
    #print 'importTempTables',ret
    for row in ret:
      #print 'importTempTables updating job',row
      for column in ['jobnumber','finishtime','status','evaluation','jobtitle']:
        self.db.execute('UPDATE Jobs SET '+column+'=? WHERE JobId =?',(row[self.db.JOBITEMS.index(column)],row[0]))

    '''
    if self.stats.get('importMax',None) is not None:
      lastJobNumber = str(int(self.stats['importMax']) + self.stats['incrJobNumber'])
      print 'importTempTables lastJobNumber',lastJobNumber
      self.db.execute('Update Projects SET LastJobNumber=? WHERE ProjectId=?',(lastJobNumber,self.projectId))
    else:
    '''
    self.db.resetLastJobNumber(self.projectId,commit=False)
    self.db.commit()
    self.db.execute('PRAGMA foreign_keys = ON')

  def listTempJobs(self,label=''):
    # Dignostic
    self.db.execute('SELECT JobId,JobNumber,TaskName,Excl,NewJobNumber FROM TempJobs')
    rv = self.db.fetchAll2PyList([str,str,str,int,str])
    print(label)
    print('JobId,JobNumber,TaskName,Excl,NewJobNumber')
    for item in rv: print(item)
 
  def listTempFiles(self,label=''):
    # Dignostic
    self.db.execute('SELECT TempJobs.JobNumber,TempFiles.FileId,TempFiles.FileName,TempFiles.Excl FROM TempFiles INNER JOIN TempJobs ON TempFiles.JobId = TempJobs.JobId')
    rv = self.db.fetchAll2PyList([str,str,str,int])
    print(label)
    print('TempJobs.JobNumber,TempFiles.FileId,TempFiles.FileName,TempFiles.Excl')
    for item in rv: print(item)

  def listTempFileUses(self,label=''):
    # Dignostic
    self.db.execute('SELECT FileId,JobId,Excl FROM TempFileUses')
    rv = self.db.fetchAll2PyList([str,str,int])
    print(label)
    for item in rv: print(item)


  def loadTempTables(self):

    for table,colList in [['Jobs',self.db.JOBITEMS],
                          ['Files',self.db.FILEITEMS],
                          ['ImportFiles',self.db.IMPORTFILEITEMS],
                          ['ExportFiles',self.db.EXPORTFILEITEMS],
                          ['FileUses',self.db.FILEUSEITEMS],
                          ['Comments',self.db.COMMENTITEMS],
                          ['XData',self.db.XDATAITEMS],
                          ['JobKeyValues',self.db.JOBKEYVALUEITEMS],
                          ['JobKeyCharValues',self.db.JOBKEYCHARVALUEITEMS],
                          ['FileAssociations',self.db.FILEASSOCIATIONITEMS],
                          ['FileAssociationMembers',self.db.FILEASSOCIATIONMEMBERITEMS]]:
      com = 'SELECT '+colList[0]
      for item in colList[1:]: com = com + ','+item
      com = ' INTO '+table+' FROM Temp'+table+' WHERE Excl = 0'
      #print 'loadTempTables com',com
      self.db.execute(com)
    self.db.commit()

  def importThisFile(self,jobNumber,fileName=None):
    # Invoked from CCP4Export.ImportProjectThread to return if directory/file should be copied into project
    # and if it needs renaming due to changed job number
    # Also resets TempJobs.Excl and  TempFiles.Excl to -1 or -2 to flag that job/file have been copied
    if fileName is None:
      self.db.execute('SELECT JobId,Excl,NewJobNumber FROM TempJobs WHERE JobNumber = ?',(jobNumber,))
      ret = self.db.fetchAll2PyList([str,int,str])
      #print 'CDbXml.importThisFile',jobNumber,ret
      if len(ret)==0:
        return False,None
      elif ret[0][1] ==0:
        self.db.execute('UPDATE TempJobs SET Excl=-1 WHERE JobNumber = ?',(jobNumber,))
        if self._diagnostic: print('CDbXml.importThisFile TempJobs',jobNumber,-1)
        return True,ret[0][2]
      elif ret[0][1] ==1:
        # Job is in db but is unfinished
        self.db.execute('UPDATE TempJobs SET Excl=-2 WHERE JobNumber = ?',(jobNumber,))
        if self._diagnostic: print('CDbXml.importThisFile TempJobs',jobNumber,-2)
        return True,ret[0][2]
      elif ret[0][1] == -1:
        return True,ret[0][2]
      else:
        return False,None
    else:
      # Importing a file
      self.db.execute('SELECT TempJobs.JobId,TempJobs.Excl,TempJobs.NewJobNumber,TempFiles.FileId,TempFiles.Excl FROM TempJobs INNER JOIN TempFiles ON TempJobs.JobId=TempFiles.JobId LEFT OUTER JOIN TempImportFiles ON TempFiles.FileId=TempImportFiles.FileId WHERE TempJobs.JobNumber = ? AND TempFiles.FileName=? AND TempImportFiles.ImportId IS NULL',(jobNumber,fileName))
      ret = self.db.fetchAll2PyList([str,int,str,str,int])
      #print 'CDbXml.importThisFile',jobNumber,fileName,ret
      if len(ret)==0:
        # It could be a job that we want to import but a file that is not in db
        self.db.execute('SELECT JobId,Excl,NewJobNumber FROM TempJobs WHERE JobNumber = ?',(jobNumber,))
        ret = self.db.fetchAll2PyList([str,int,str])
        if len(ret)==0:
          return False,None
        elif ret[0][1] <= 0:
          return True,ret[0][2]
        else:
          return False,None
      elif ret[0][4] == 0:
        # Required job and a file in the db
        self.db.execute('UPDATE TempFiles SET Excl=-1 WHERE FileId = ?',(ret[0][3],))
        if self._diagnostic:
           self.db.execute('SELECT FileId,FileName,Excl FROM TempFiles WHERE FileId = ?',(ret[0][3],))
           fileList = self.db.fetchAll2PyList([str,str,int])
           if self._diagnostic: print('CDbXml.importThisFile TempFiles',fileList)
        return True,ret[0][2]
      elif ret[0][4] == -1:
        return True,ret[0][2]
      else:
        # Explicitly excluded file/job
        return False,None

  def setAllToImport(self):
    self.db.execute('UPDATE TempJobs SET Excl=-1')
    self.db.execute('UPDATE TempFiles SET Excl=-1')

  def getRowFromEtree(self,ele,attributeList,attributeTypes):
    rv = []
    idx = -1
    for att  in attributeList:
      idx = idx + 1
      '''
      ele = tree.get(att)
      if ele is None or ele.text is None:
        rv.append(None)
      else:
        rv.append(attributeTypes[idx](ele.text))
      '''
      val = ele.get(att,None)
      if val is None or val=='None':
        rv.append(None)
      else:
        try:
          rv.append(attributeTypes[idx](val))
        except:
          rv.append(None)
          self.errReport.append(self.__class__,212,'Data type:'+str(att)+' is: '+str(val))
    return rv


class CDatabaseDef:
  def __init__(self,fileName):
    self.headerArray = {}
    self.dataArray = {}
    self.project_name =''
    self.project_dir = ''
    self.jobs = {}
    self.headerArray,typeArray,self.dataArray = self.read(fileName)
    #print 'dataArray',self.dataArray

    # Get project name and dir from header
    #if self.headerArray.has_key('PROJECT'):
    if 'PROJECT' in self.headerArray:
      self.project_name = self.headerArray['PROJECT'][0]
      if len(self.headerArray['PROJECT'])>1:
        self.project_dir = self.headerArray['PROJECT'][1]
    #print ('Read project',self.project_name,self.project_dir)

    self.nJobs = int(self.dataArray.get('NJOBS',0))

    self.jobs = self.extractJobs(self.dataArray)

  def read (self,fileName='',headerOnly=0):
      '''
      Read the def file into dicts for header, types and data
      A database.def should not have any type information
      This is left in the code in case want to reuse elsewhere
      '''
      headerArray = {}
      typeArray = {}
      dataArray = {}
      try:
        f = open(fileName)
      except:
        print(('ERROR opening',fileName))
        return headerArray,typeArray,dataArray
      if not f:
        print(('ERROR opening',fileName))
        return headerArray,typeArray,dataArray
      try:
        content = f.readlines()
        f.close()
      except:
        print(('ERROR reading ',fileName))
        return headerArray,typeArray,dataArray


      for line in content:
        words = self.splitDefLine(line)
        #print 'CCP4Data.parseDefFile',words
        if len(words)>=2:
          if words[0] == '#CCP4I' and len(words)>=3:
            #print 'header',words
            headerArray[words[1]] = words[2:]
          elif words[0][0] == '#':
            pass
          elif not headerOnly:
           if len(words)>2:
             if  words[1][0]=='_':
               typeArray[words[0]] = words[1][1:]
             else:
               typeArray[words[0]] = words[1]
             idata = 2
           else:
             idata = 1
           if words[idata].strip('"')=='':
             dataArray[words[0]] = ''
           else:
             dataArray[words[0]]=words[idata]

      return  headerArray,typeArray,dataArray

  def splitDefLine(self,line=''):
    m = re.search(r'(.*)\"(.*)\"(.*)',line)
    if not m:
      return line.split()
    else:
      a,b,c = m.groups()
      a = a.strip()
      c = c.strip()
      rv = []
      if a:rv.extend(a.split())
      rv.append(b)
      if c: rv.append(c.split())
      return rv

  def splitDefList(self,line=''):
    '''
    Parse the INPUT_FILE and OUTPUT_FILE from database.def
    This is a space-separated list
    If the file name includes spaces then the name is enclosed in curly brace
    '''
    rv = []
    while len(line)>0:
      m = re.search(r'(.*?)\{(.*?)\}(.*)',line)
      if not m:
        rv.extend(line.split())
        line = ''
      else:
        a,b,c = m.groups()
        #print 'splitDefList',a,'*',b,'*',c
        a = a.strip()
        c = c.strip()
        if a:rv.extend(a.split())
        rv.append(b)
        line = c.strip()
        #print 'new line',line
    return rv


  def extractJobs(self,dataArray):

    jobs = {}

    for key,value in list(dataArray.items()):
      try:
        if key.count(',') == 1:
          dataType,jobId = key.split(',')
          try:
            jobId = int(jobId)
          except:
            jobId = -1
          if jobId>0:
            # Initialise the data for one job
            if not jobId in jobs:
              jobs[jobId] = { 'STATUS' : '',
                              'DATE' : 0,
                              'LOGFILE' : '',
                              'TASKNAME' : '',
                              'TITLE' : '',
                              'INPUT_FILES' : [],
                              'INPUT_FILES_DIR' : [],
                              'OUTPUT_FILES' : [],
                              'OUTPUT_FILES_DIR' : []   }
            # Check that the dataType is one of the recognised properties for a job
            if dataType in jobs[jobId]:
              if ['INPUT_FILES','OUTPUT_FILES'].count(dataType):
                jobs[jobId][dataType] = self.splitDefList(value)
              elif ['INPUT_FILES_DIR','OUTPUT_FILES_DIR'].count(dataType):
                jobs[jobId][dataType] = value.split(' ')
              else:
                jobs[jobId][dataType] = value
      except:
         print(('ERROR interpreting ',key))

    return jobs



#=========================================================================================================
import unittest
from core import CCP4Utils

def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testSqliteDb)
  #unittest.TestLoader.testMethodPrefix = 'testx'
  #suite = unittest.defaultTestLoader.loadTestsFromTestCase(testQtDb)
  return suite

def testModule():

  from core import CCP4File
  CCP4File.CDataFile.BLOCK_LOAD_FILES = True

  from core import CCP4Utils
  if os.path.exists(testQtDb.TESTDBFILE): os.remove(testQtDb.TESTDBFILE)
  if os.path.exists(testSqliteDb.TESTDBFILE): os.remove(testSqliteDb.TESTDBFILE)
  for testDir in [testDb.TESTDBDIR]:
    if not os.path.exists(testDir):
      os.mkdir(testDir)
      jobsDir = os.path.join(testDir,'CCP4_JOBS')
      os.mkdir(jobsDir)

      os.mkdir(os.path.join(jobsDir,'job_1'))
      os.mkdir(os.path.join(jobsDir,'job_2'))
      os.mkdir(os.path.join(jobsDir,'job_3'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_1'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_2'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_3'))
      os.mkdir(os.path.join(jobsDir,'job_5'))
      # Emulating imported data
      for name,job in [['starting_data.mtz','job_2'],['starting_data.seq','job_1']]:
        fileName = os.path.join(testDir,name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
        fileName = os.path.join(testDir,'CCP4_JOBS',job,name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      # Emulating pipeline
      for name in ['model_1.pdb','model_2.pdb']:
        fileName = os.path.join(jobsDir,'job_3','job_1',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      for name in ['built.pdb','built.mtz']:
        fileName = os.path.join(jobsDir,'job_3','job_2',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      for name in ['built.pdb','built.mtz']:
        fileName = os.path.join(jobsDir,'job_3','job_3',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      fileName = os.path.join(jobsDir,'job_5','XYZOUT.pdb')
      CCP4Utils.saveFile(fileName=fileName,text='Dummy file')

  from core import CCP4Modules
  pm = CCP4Modules.PROJECTSMANAGER()

  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

  CCP4File.CDataFile.BLOCK_LOAD_FILES = False

class testDb(unittest.TestCase):
  TESTDBDIR = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.test.dir')


  def test0User(self):
    uid0 = self.db.createUser('metoo')
    uid1 = self.db.getUserId('metoo')
    self.assertEqual(uid0,uid1,'Failed to create and retrieve user id')
    try:
      self.db.createUser('metoo')
    except CException as e:
      self.assertEqual( e[0]['code'], 120,'Wrong error when attempting to create user with same name')
    except:
      self.fail('Unknown error when attempting to create user with same name')
    else:
      self.fail('No error when attempting to create user with same name')


    self.db._userRole = USER_ROLE_USER
    try:
      self.db.createUser('another')
    except CException as e:
      self.assertEqual( e[0]['code'], 106,'Wrong error when attempting to create user with wrong permission')
    except:
      self.fail('Unknown error when attempting to create user with wrong permission')
    else:
      self.fail('No error when attempting to create user with wrong permission')

    self.db._userRole = USER_ROLE_OWNER
    self.db.updateUser('metoo','another')

    self.db.commit()
    self.db.close()

    self.db = CDbApi(mode=self.mode,fileName=self.TESTDBFILE,userName='me')
    uList = self.db.listUsers(True)
    print('test0User uList',uList)
    self.assertEqual(uList[1][1],'another','Failed changing user name and listing')

    info = self.db.getUserInfo(uid0)
    #print 'after getUserInfo',info
    self.assertEqual(info['username'],'another','Error calling getUserInfo')


  def test1Project(self):

    pid0 = self.db.createProject('myproject',projectDirectory=self.TESTDBDIR)
    pid1 = self.db.getProjectId('myproject')
    self.assertEqual(pid0,pid1,'Error creating project and retrieving project ID')

    baby = self.db.createProject('babyproject',userName='me',parentProjectId=pid1,projectDirectory=os.path.join(self.TESTDBDIR,'baby'))
    self.db.commit()
    self.db.close()

    self.db = CDbApi(mode=self.mode,fileName=self.TESTDBFILE,userName='me')
    #print 'test1Project _userName',self.db._userName
    pList = self.db.listProjects(True)
    self.assertEqual(pList[1][3],pid0,'Failed creating child project and listing')

    self.db.updateProject(baby,key='projectDirectory',value=os.path.join(self.TESTDBDIR,'baby2'))
    pList = self.db.listProjects(True)
    self.assertEqual(os.path.split(pList[1][2])[-1],'baby2','Failed updating project directory')


    '''
    # Recent projects concept currently dropped
    #self.db.listContexts(True)
    recent = self.db.getRecentProjects()
    #print 'test1Project getRecentProjects',recent
    self.assertEqual(recent,[[1,'myproject']],'Failed retrieving recent projects')
    '''


  def test2Aliases(self):
     from core import CCP4Utils
     self.db.setDirectoryAlias('CCP4I2_TOP',CCP4Utils.getCCP4I2Dir())
     self.db.setDirectoryAlias('HOME',CCP4Utils.getHOME())
     self.db.commit()

     aliasList = self.db.listDirectoryAliases(toTerm=True)
     self.assertEqual(len(aliasList),2,'Alias list wrong length')
     self.assertTrue(('CCP4I2_TOP' in aliasList),'Alias list does not have CCP4I2_TOP')



  def test3Jobs(self):
    #pid = self.db.createProject('myproject',projectDirectory=self.TESTDBDIR)
    pid = self.db.getProjectId('myproject')
    print('test3Jobs pid',pid)
    self.db.listProjects(True)
    from core import CCP4ModelData,CCP4XtalData

    # Two jobs loading contents and exptal data
    jid1 =  self.db.createJob(pid,'contents',status=JOB_STATUS_FINISHED)
    fo1 = CCP4ModelData.CSeqDataFile(project='myproject',relPath='CCP4_JOBS/job_1',baseName='starting_data.seq')
    f1=self.db.createFile(jobId=jid1,fileObject=fo1,sourceFileName=os.path.join(CCP4Utils.getTestTmpDir(),'starting_data.seq'))

    jid2 =  self.db.createJob(pid,'import_xtal',status=JOB_STATUS_FINISHED)
    fo2 = CCP4XtalData.CMtzDataFile(project='myproject',relPath='CCP4_JOBS/job_2',baseName='starting_data.mtz')
    f2=self.db.createFile(jobId=jid2,fileObject=fo2,sourceFileName=os.path.join(CCP4Utils.getTestTmpDir(),'starting_data.mtz'))


    # An 'auto mr' pipeline
    jid3 =  self.db.createJob(pid,'auto_mr')
    self.db.createFileUse(jobId=jid3,fileId=f1)
    self.db.createFileUse(jobId=jid3,fileId=f2)
    jobNo = self.db.getJobInfo(jid3,'jobNumber')
    self.assertEqual(jobNo,'3','Error retrieving jobNumber')

    # mr run creating two pdbs files
    jid4 =  self.db.createJob(pid,'mr',parentJobId=jid3,jobTitle='testing a parentJobId and jobTitle')
    self.db.createFileUse(jobId=jid4,fileId=f1)
    self.db.createFileUse(jobId=jid4,fileId=f2)
    fo3 = CCP4ModelData.CPdbDataFile(project='myproject',relPath='CCP4_JOBS/job_3/job_1',baseName='model_1.pdb')
    fo4 = CCP4ModelData.CPdbDataFile(project='myproject',relPath='CCP4_JOBS/job_3/job_1',baseName='model_2.pdb')
    f3=self.db.createFile(jobId=jid4,fileObject=fo3)
    f4=self.db.createFile(jobId=jid4,fileObject=fo4)
    self.assertEqual(self.db.getJobInfo(jid4,'ParentJobId'),jid3,'Error setting and retrieving parent job number')
    self.db.updateJobStatus(jid4,JOB_STATUS_FINISHED)

    # Two runs of model building
    jid5 =  self.db.createJob(pid,'build',parentJobId=jid3,jobTitle='build')
    self.db.createFileUse(jobId=jid5,fileId=f2)
    self.db.createFileUse(jobId=jid5,fileId=f3)
    fo5 = CCP4ModelData.CPdbDataFile(project='myproject',baseName='built.pdb',relPath='CCP4_JOBS/job_3/job_2')
    f5=self.db.createFile(jobId=jid5,fileObject=fo5,fileTypeName='chemical/x-pdb')
    self.db.updateJobStatus(jid5,JOB_STATUS_FINISHED)

    fileIds = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_IN,mode='fileId')
    self.assertTrue(fileIds==[f2,f3] or fileIds==[f3,f2],'Failed restoring input file ids')
    fileNames = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_OUT,mode='fileName')
    self.assertEqual(fileNames,['built.pdb'],'Failed restoring input file ids')

    # Lets say second model building failed
    jid6 =  self.db.createJob(pid,'build',parentJobId=jid3,jobTitle='build')
    self.db.createFileUse(jobId=jid6,fileId=f2)
    self.db.createFileUse(jobId=jid6,fileId=f4)
    #f6=self.db.createFile(jobId=jid6,fileName='job_6_built.pdb',relPath='CCP4_JOBS/job_6',fileTypeName='text/pdb')
    self.db.updateJobStatus(jid6,JOB_STATUS_FAILED)

    # Finish off the pipeline

    self.db.updateJobStatus(status='Finished',projectName='myproject',jobNumber='3')
    self.db.createFileUse(jobId=jid3,fileId=f5,role=FILE_ROLE_OUT)
    info = self.db.getJobInfo(jid3)
    self.assertEqual(info['status'],'Finished','Error changing job status')
    self.assertEqual(info['taskname'],'auto_mr','Error retrieving taskName')
    fileTypeList = self.db.getJobFiles(jobId=jid3,role=FILE_ROLE_OUT,mode='fileType')
    # Output should be one PDB file
    #print 'test3Jobs fileTypeList',fileTypeList
    self.assertEqual(fileTypeList,['chemical/x-pdb'])

    fileList = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_OUT,mode='fullPath')
    self.assertEqual(fileList,[os.path.join(self.TESTDBDIR,'CCP4_JOBS','job_3','job_2','built.pdb')],'Failed running getJobFiles')
    fileId,jobId = self.db.matchFileName(fileName = os.path.join(self.TESTDBDIR,'CCP4_JOBS','job_3','job_2','built.pdb'))
    #print 'test3Jobs matchFileName',fileId,jobId ,'f5',f5
    self.assertEqual(fileId,f5,'Error using matchFileName')

    info = self.db.getProjectJobListInfo(projectName='myproject')
    #print 'testJob getJobListInfo',info


  def test4Task(self):
    pid = self.db.getProjectId('myproject')
    jid = self.db.createJob(pid,'pdbset')
    self.db.commit()
    """
    # make config file with reference to test database file
    from core import CCP4Modules,CCP4Config,CCP4Utils
    config = CCP4Config.CConfig()
    config.set('dbMode',self.mode)
    config.set('dbFile',self.TESTDBFILE)
    config.set('dbUser','me')
    import os
    from core import CCP4Utils

    configFile = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi_test_config.data.xml')
    print 'test4Task configFile',configFile
    config.saveDataToXml(configFile)
    # Set the configFile in JOBCONTROLLER so it is part of arguments to subprocess
    JC = CCP4Modules.JOBCONTROLLER()
    JC.setConfigFile(configFile)
    JC.startTimer()

    controlFile= os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi.params.xml')
    self.db.updateJobControlFile(jobId=jid,controlFile=controlFile)
    self.db.updateJobStatus(jobId=jid,status=JOB_STATUS_QUEUED)
    self.db.commit()

    import time
    time.sleep(30)
    status = self.db.getJobInfo(jid,'status')
    print 'test4Task status',status
    self.assertEqual(status,'Finished','Failed to run task and update database')
    """


  def test5glean(self):
    # This will e jobNumber=5
    pid = self.db.getProjectId('myproject')
    import os
    from core import CCP4Utils
    jid7 =  self.db.createJob(projectId=pid,taskName='pdbset',jobTitle='test5glean')
    import shutil
    from core import CCP4Modules
    shutil.copyfile(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi.params.xml'),CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jid7,mode='PARAMS'))
    self.db.updateJobStatus(jid7,status='Finished')
    #self.db.setDiagnostic(True)
    #print 'test5glean jid7',jid7,CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jid7,mode='PARAMS')
    rv = self.db.gleanJobFiles(jobId=jid7)
    #print 'from gleanJobFiles',rv.report()

    fileList = self.db.getJobFiles(jobId=jid7,role=FILE_ROLE_IN,mode='fullPath')
    #print 'test5glean fileList',fileList
    self.assertEqual(len(fileList),1,'Error gleaning files from params - wrong number of input files')
    self.assertEqual(os.path.split(fileList[0])[1],'model_1.pdb','Error gleaning files from params - wrong file name')

    fileList = self.db.getJobFiles(jobId=jid7,role=FILE_ROLE_OUT)
    self.assertEqual(len(fileList),1,'Error gleaning files from params - wrong number of output files')
    filePath = self.db.getFullPath(fileId=fileList[0])
    self.assertEqual(os.path.split(filePath)[1],'XYZOUT.pdb','Error gleaning files from params - wrong output file name')

  def test6glean(self):
    pid = self.db.getProjectId('myproject')
    import os
    from core import CCP4Utils
    pf =  os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi_2.params.xml')
    import shutil
    from core import CCP4Modules
    jid8 =  self.db.createJob(pid,'newProject',jobTitle='test6glean')
    shutil.copyfile(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi_2.params.xml'),CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jid8,mode='PARAMS'))

    CCP4Modules.PROJECTSMANAGER().importFiles(jobId=jid8)
    self.db.updateJobStatus(jid8,status='Finished')
    rv = self.db.gleanJobFiles(jobId=jid8)

    containsSe = self.db.getXDataByJobContext(contextJobId=jid8,dataClass='CContainsSeMet',projectId=pid)
    #print 'CCP4DbApi.test6glean containsSe',containsSe

    import os
    fileName=os.path.join(CCP4Utils.getTestTmpDir(),'export_db.xml')
    if os.path.exists(fileName): os.remove(fileName)
    self.db.exportProjectXml(projectId=pid,fileName=fileName)
    self.assertTrue(os.path.exists(fileName),'Export xml file not created')

  def test7comment(self):

    jid3 =  self.db.getJobId(projectName='myproject',jobNumber='3')
    cid = self.db.createComment(jobId=jid3,comment='Well whatever really')

    commentList = self.db.getComments(jobId=jid3)
    #print 'test7comment commentList',commentList
    self.assertEqual('Well whatever really',commentList[0][3],'Error creating and retrieving comment')
  '''
  def test8ImportProject(self):
    projectXml = '/Users/lizp/Desktop/t4b.ccp4i2db.xml'
    projectDirectory = '/Users/lizp/Desktop/t4'
    dbImport = CDbXml(db=self.db,projectDirectory=projectDirectory)
    err = dbImport.load(projectXml)
    print 'test8ImportProject errorReport',err


  def test9ImportProject(self):
    projectXml = '/Users/lizp/Desktop/t4b.ccp4i2db.xml'
    projectDirectory = '/Users/lizp/Desktop/t4'
    dbImport = CDbXml(db=self.db,projectDirectory=projectDirectory)
    err = dbImport.loadTable(projectXml)
    print 'test9ImportProject errorReport',err
  '''

class testSqliteDb(testDb):

  TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testsqlite.db')

  def setUp(self):

    self.db = CDbApi(mode='sqlite',userName='me',userPassword='foo',fileName=testSqliteDb.TESTDBFILE)
    self.mode = 'sqlite'
    from core import CCP4Modules
    CCP4Modules.PROJECTSMANAGER().setDatabase(self.db)

class testQtDb(testDb):

  TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testqt.db')

  def setUp(self):
    self.db = CDbApi(mode='qt_sqlite',fileName=testQtDb.TESTDBFILE,userName='me',userPassword='foo')
    self.mode = 'qt_sqlite'
    from core import CCP4Modules
    CCP4Modules.PROJECTSMANAGER().setDatabase(self.db)
