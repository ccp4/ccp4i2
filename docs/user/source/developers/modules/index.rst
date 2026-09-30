###################
CCP4i2 Data Classes
###################

.. toctree::
    :maxdepth: 2

    CCP4Annotation.rst
    CCP4ComFilePatchManager.rst
    CCP4CootData.rst
    CCP4CustomTaskManager.rst
    CCP4Data.rst
    CCP4File.rst
    CCP4ImportedJobManager.rst
    CCP4MathsData.rst
    CCP4ModelData.rst
    CCP4PerformanceData.rst
    CCP4Preferences.rst
    CCP4RefmacData.rst
    CCP4WorkflowManager.rst
    CCP4XtalData.rst

Current Data Classes
~~~~~~~~~~~~~~~~~~~~

CCP4Annotation
^^^^^^^^^^^^^^

| `CAnnotation <./CCP4Annotation.html#CAnnotation>`__ Annotation text
  with user id and time
| `CAnnotationList <./CCP4Annotation.html#CAnnotationList>`__ A list of
  annotation
| `CAuthor <./CCP4Annotation.html#CAuthor>`__ Placeholder for
  bibliographic author
| `CBibReference <./CCP4Annotation.html#CBibReference>`__ Bibliographic
  reference
| `CBibReferenceGroup <./CCP4Annotation.html#CBibReferenceGroup>`__ Set
  of bibliographic references for a task
| `CFont <./CCP4Annotation.html#CFont>`__ Simplified Qt font options
| `CHostName <./CCP4Annotation.html#CHostName>`__ Computer name
| `CHostname <./CCP4Annotation.html#CHostname>`__
| `CMetaDataTag <./CCP4Annotation.html#CMetaDataTag>`__ This class will
  extend list of enumerators if new value for string is entered
| `CMetaDataTagList <./CCP4Annotation.html#CMetaDataTagList>`__
| `CServerGroup <./CCP4Annotation.html#CServerGroup>`__ One or more
  compute servers used in "remote" running
| `CTime <./CCP4Annotation.html#CTime>`__ The time. Uses Python time
  module
| `CUserAddress <./CCP4Annotation.html#CUserAddress>`__ User id and
  platform node
| `CUserId <./CCP4Annotation.html#CUserId>`__ A user ID

CCP4CootData
^^^^^^^^^^^^

| `CCootHistoryDataFile <./CCP4CootData.html#CCootHistoryDataFile>`__

CCP4Data
^^^^^^^^

| `CBoolean <./CCP4Data.html#CBoolean>`__ A Boolean
| `CFloat <./CCP4Data.html#CFloat>`__ A float
| `CFloatRange <./CCP4Data.html#CFloatRange>`__ Two floats defining
  start and end of range
| `CFollowFromJob <./CCP4Data.html#CFollowFromJob>`__
| `CI2DataType <./CCP4Data.html#CI2DataType>`__
| `CInt <./CCP4Data.html#CInt>`__ An integer
| `CIntRange <./CCP4Data.html#CIntRange>`__ Two integers defining start
  and end of range
| `CJobStatus <./CCP4Data.html#CJobStatus>`__
| `CJobTitle <./CCP4Data.html#CJobTitle>`__
| `COneWord <./CCP4Data.html#COneWord>`__ A single word string - no
  white space
| `CPatchSelection <./CCP4Data.html#CPatchSelection>`__
| `CRange <./CCP4Data.html#CRange>`__ Base class for CIntRange and
  CFloatRange
| `CRangeSelection <./CCP4Data.html#CRangeSelection>`__
| `CString <./CCP4Data.html#CString>`__ A string
| `CTable <./CCP4Data.html#CTable>`__
| `CUUID <./CCP4Data.html#CUUID>`__

CCP4File
^^^^^^^^

| `CDataFile <./CCP4File.html#CDataFile>`__ A data file - expected to
  have associated class for file contents
| `CDataFileContent <./CCP4File.html#CDataFileContent>`__ Base class for
  classes holding file contents
| `CExePath <./CCP4File.html#CExePath>`__
| `CExePathList <./CCP4File.html#CExePathList>`__
| `CExportedFile <./CCP4File.html#CExportedFile>`__
| `CExportedFileList <./CCP4File.html#CExportedFileList>`__
| `CFileFunction <./CCP4File.html#CFileFunction>`__ List of recognised
  XML file functions
| `CFilePath <./CCP4File.html#CFilePath>`__ A file path
| `CI2XmlDataFile <./CCP4File.html#CI2XmlDataFile>`__ A reference to an
  XML file with CCP4i2 Header
| `CI2XmlHeader <./CCP4File.html#CI2XmlHeader>`__ Container for header
  info from XML file
| `CMmcifData <./CCP4File.html#CMmcifData>`__ Generic mmCIF data. This
  is intended to be a base class for other classes specific to
  coordinates, reflections or geometry data.
| `CMmcifDataFile <./CCP4File.html#CMmcifDataFile>`__ A generic mmCIF
  format file. This is intended to be a base class for other classes
  specific to coordinates, reflections or geometry data.
| `CPostscriptDataFile <./CCP4File.html#CPostscriptDataFile>`__ An
  postscript format file
| `CProjectId <./CCP4File.html#CProjectId>`__ The CCP4i2 database
  project id - a global unique id
| `CProjectName <./CCP4File.html#CProjectName>`__ The name of a CCP4i
  project or directory alias
| `CSceneDataFile <./CCP4File.html#CSceneDataFile>`__ An xml format file
  for defining scene in CCP4mg.
| `CSearchPath <./CCP4File.html#CSearchPath>`__
| `CSearchPathList <./CCP4File.html#CSearchPathList>`__
| `CTextDataFile <./CCP4File.html#CTextDataFile>`__ A text data file
| `CVersion <./CCP4File.html#CVersion>`__ A (string) version number of
  the form n.m.i
| `CXmgrDataFile <./CCP4File.html#CXmgrDataFile>`__ An xmgr format file.
  This is the input format for xmgrace, as output by scala or aimless
| `CXmlDataFile <./CCP4File.html#CXmlDataFile>`__ A reference to an XML
  file

CCP4MathsData
^^^^^^^^^^^^^

| `CAngle <./CCP4MathsData.html#CAngle>`__ An angle
| `CEulerRotation <./CCP4MathsData.html#CEulerRotation>`__
| `CMatrix33 <./CCP4MathsData.html#CMatrix33>`__
| `CTransformation <./CCP4MathsData.html#CTransformation>`__
| `CXyz <./CCP4MathsData.html#CXyz>`__
| `CXyzBox <./CCP4MathsData.html#CXyzBox>`__

CCP4ModelData
^^^^^^^^^^^^^

| `CAtomSelection <./CCP4ModelData.html#CAtomSelection>`__
| `CBlastData <./CCP4ModelData.html#CBlastData>`__
| `CBlastDataFile <./CCP4ModelData.html#CBlastDataFile>`__
| `CBlastItem <./CCP4ModelData.html#CBlastItem>`__
| `CChemComp <./CCP4ModelData.html#CChemComp>`__
| `CContainsSeMet <./CCP4ModelData.html#CContainsSeMet>`__
| `CDictData <./CCP4ModelData.html#CDictData>`__
| `CDictDataFile <./CCP4ModelData.html#CDictDataFile>`__ A refmac
  dictionary file
| `CElement <./CCP4ModelData.html#CElement>`__ Chemical element
| `CEnsemble <./CCP4ModelData.html#CEnsemble>`__ An ensemble of models.
  Typically, this would be a set of related PDB files, but models could
  also be xtal or EM maps. This should be indicated by the types entry.
  A single ensemble is a CList of structures.
| `CEnsembleList <./CCP4ModelData.html#CEnsembleList>`__
| `CEnsemblePdbDataFile <./CCP4ModelData.html#CEnsemblePdbDataFile>`__ A
  PDB coordinate file containing ensemble of structures as 'NMR' models
| `CHhpredData <./CCP4ModelData.html#CHhpredData>`__
| `CHhpredDataFile <./CCP4ModelData.html#CHhpredDataFile>`__
| `CHhpredItem <./CCP4ModelData.html#CHhpredItem>`__
| `CHomolog <./CCP4ModelData.html#CHomolog>`__
| `CHomologList <./CCP4ModelData.html#CHomologList>`__
| `CMDLMolDataFile <./CCP4ModelData.html#CMDLMolDataFile>`__ A molecule
  definition file (MDL)
| `CMonomer <./CCP4ModelData.html#CMonomer>`__ A monomer compound.
  ?smiles
| `CPdbData <./CCP4ModelData.html#CPdbData>`__ Contents of a PDB file -
  a subset with functionality for GUI
| `CPdbDataFile <./CCP4ModelData.html#CPdbDataFile>`__ A PDB coordinate
  file
| `CPdbDataFileList <./CCP4ModelData.html#CPdbDataFileList>`__
| `CPdbEnsembleItem <./CCP4ModelData.html#CPdbEnsembleItem>`__
| `CRefmacNcs <./CCP4ModelData.html#CRefmacNcs>`__ Definition of a NCS
  for Refmac - not implemented yet
| `CRefmacRigidDomain <./CCP4ModelData.html#CRefmacRigidDomain>`__
  Definition of a rigid domain for Refmac
| `CResidueRange <./CCP4ModelData.html#CResidueRange>`__ A residue range
  selection
| `CResidueRangeList <./CCP4ModelData.html#CResidueRangeList>`__ A list
  of residue range selections
| `CSeqAlignDataFile <./CCP4ModelData.html#CSeqAlignDataFile>`__ A
  (multiple) sequence alignment file
| `CSeqDataFile <./CCP4ModelData.html#CSeqDataFile>`__ A sequence file
| `CSeqDataFileList <./CCP4ModelData.html#CSeqDataFileList>`__
| `CSequence <./CCP4ModelData.html#CSequence>`__ A string of sequence
  one-letter codes Need to be able to parse common seq file formats Do
  we need to support alternative residues What about nucleic/polysach?
| `CSequenceAlignment <./CCP4ModelData.html#CSequenceAlignment>`__ An
  alignment of two or more sequences. Each sequence is obviously related
  to class CSequence, but will also contain gaps relevant to the
  alignment. We could implement the contents as a list of CSequence
  objects? The alignment is typically formatted in a file as consecutive
  or interleaved sequences.
| `CSequenceMeta <./CCP4ModelData.html#CSequenceMeta>`__
| `CTLSDataFile <./CCP4ModelData.html#CTLSDataFile>`__ A refmac TLS file

CCP4PerformanceData
^^^^^^^^^^^^^^^^^^^

| `CAtomCountPerformance <./CCP4PerformanceData.html#CAtomCountPerformance>`__
| `CDataReductionPerformance <./CCP4PerformanceData.html#CDataReductionPerformance>`__
| `CExpPhasPerformance <./CCP4PerformanceData.html#CExpPhasPerformance>`__
| `CModelBuildPerformance <./CCP4PerformanceData.html#CModelBuildPerformance>`__
| `CPerformanceIndicator <./CCP4PerformanceData.html#CPerformanceIndicator>`__
| `CPhaseErrorPerformance <./CCP4PerformanceData.html#CPhaseErrorPerformance>`__
| `CRefinementPerformance <./CCP4PerformanceData.html#CRefinementPerformance>`__
| `CSuperposePerformance <./CCP4PerformanceData.html#CSuperposePerformance>`__
| `CTestObsConversionsPerformance <./CCP4PerformanceData.html#CTestObsConversionsPerformance>`__

CCP4RefmacData
^^^^^^^^^^^^^^

| `CRefmacAnomalousAtom <./CCP4RefmacData.html#CRefmacAnomalousAtom>`__
| `CRefmacRestraintsDataFile <./CCP4RefmacData.html#CRefmacRestraintsDataFile>`__
| `CRefmacRigidGroupItem <./CCP4RefmacData.html#CRefmacRigidGroupItem>`__
| `CRefmacRigidGroupList <./CCP4RefmacData.html#CRefmacRigidGroupList>`__
| `CRefmacRigidGroupSegment <./CCP4RefmacData.html#CRefmacRigidGroupSegment>`__

CCP4XtalData
^^^^^^^^^^^^

| `CAltSpaceGroup <./CCP4XtalData.html#CAltSpaceGroup>`__
| `CAltSpaceGroupList <./CCP4XtalData.html#CAltSpaceGroupList>`__
| `CAnomalousColumnGroup <./CCP4XtalData.html#CAnomalousColumnGroup>`__
  Selection of F/I and AnomF/I columns from MTZ. Expected to be part of
  ab initio phasing dataset ( CDataset)
| `CAnomalousIntensityColumnGroup <./CCP4XtalData.html#CAnomalousIntensityColumnGroup>`__
  Selection of I and AnomI columns from MTZ. Expected to be part of ab
  initio phasing dataset ( CDataset)
| `CAnomalousScatteringElement <./CCP4XtalData.html#CAnomalousScatteringElement>`__
  Definition of a anomalous scattering element
| `CAsuComponent <./CCP4XtalData.html#CAsuComponent>`__ A component of
  the asymmetric unit. This is for use in MR, defining what we are
  searching for. There are similarities to CCrystalComponents and it
  should maybe be merged.
| `CAsuComponentList <./CCP4XtalData.html#CAsuComponentList>`__
| `CAtomicFormFactors <./CCP4XtalData.html#CAtomicFormFactors>`__ Table
  of form factors for element v wavelength
| `CCell <./CCP4XtalData.html#CCell>`__ A unit cell
| `CCellAngle <./CCP4XtalData.html#CCellAngle>`__ A cell angle
| `CCellLength <./CCP4XtalData.html#CCellLength>`__ A cell length
| `CColumnGroup <./CCP4XtalData.html#CColumnGroup>`__ Groups of columns
  in MTZ - probably from analysis by hklfile
| `CColumnGroupItem <./CCP4XtalData.html#CColumnGroupItem>`__ Definition
  of set of columns that form a 'group'
| `CColumnGroupList <./CCP4XtalData.html#CColumnGroupList>`__
| `CColumnType <./CCP4XtalData.html#CColumnType>`__ A list of recognised
  MTZ column types
| `CColumnTypeList <./CCP4XtalData.html#CColumnTypeList>`__ A list of
  acceptable MTZ column types
| `CCrystalComponents <./CCP4XtalData.html#CCrystalComponents>`__ A list
  of sequences, monomers and anomalous scatterers expected in a crystal
| `CCrystalCompositionLabel <./CCP4XtalData.html#CCrystalCompositionLabel>`__
  Serves as column header for CCrystalComposition - is a name for a
  composition model The composition model can be for one of three units
| `CCrystalCompositionTable <./CCP4XtalData.html#CCrystalCompositionTable>`__
  A table of crystal components v. composition models
| `CCrystalName <./CCP4XtalData.html#CCrystalName>`__
| `CDataset <./CCP4XtalData.html#CDataset>`__ The experimental data
  model for ab initio phasing
| `CDatasetList <./CCP4XtalData.html#CDatasetList>`__
| `CDatasetName <./CCP4XtalData.html#CDatasetName>`__
| `CExperimentalDataType <./CCP4XtalData.html#CExperimentalDataType>`__
  Experimental data type e.g. native or peak
| `CFPairColumnGroup <./CCP4XtalData.html#CFPairColumnGroup>`__
| `CFSigFColumnGroup <./CCP4XtalData.html#CFSigFColumnGroup>`__
| `CFormFactor <./CCP4XtalData.html#CFormFactor>`__ The for factor (Fp
  and Fpp) for a giving element and wavelength
| `CFreeRColumnGroup <./CCP4XtalData.html#CFreeRColumnGroup>`__
| `CFreeRDataFile <./CCP4XtalData.html#CFreeRDataFile>`__
| `CGenericReflDataFile <./CCP4XtalData.html#CGenericReflDataFile>`__
| `CHLColumnGroup <./CCP4XtalData.html#CHLColumnGroup>`__
| `CIPairColumnGroup <./CCP4XtalData.html#CIPairColumnGroup>`__
| `CISigIColumnGroup <./CCP4XtalData.html#CISigIColumnGroup>`__
| `CImageFile <./CCP4XtalData.html#CImageFile>`__
| `CImosflmXmlDataFile <./CCP4XtalData.html#CImosflmXmlDataFile>`__ An
  iMosflm data file
| `CImportUnmerged <./CCP4XtalData.html#CImportUnmerged>`__
| `CImportUnmergedList <./CCP4XtalData.html#CImportUnmergedList>`__
| `CMapCoeffsDataFile <./CCP4XtalData.html#CMapCoeffsDataFile>`__
| `CMapColumnGroup <./CCP4XtalData.html#CMapColumnGroup>`__
| `CMapDataFile <./CCP4XtalData.html#CMapDataFile>`__ A CCP4 Map file
| `CMergeMiniMtz <./CCP4XtalData.html#CMergeMiniMtz>`__
| `CMergeMiniMtzList <./CCP4XtalData.html#CMergeMiniMtzList>`__
| `CMiniMtzDataFile <./CCP4XtalData.html#CMiniMtzDataFile>`__
| `CMiniMtzDataFileList <./CCP4XtalData.html#CMiniMtzDataFileList>`__
| `CMmcifReflData <./CCP4XtalData.html#CMmcifReflData>`__ Reflection
  data in mmCIF format
| `CMmcifReflDataFile <./CCP4XtalData.html#CMmcifReflDataFile>`__ A
  reflection file in mmCIF format
| `CMtzColumn <./CCP4XtalData.html#CMtzColumn>`__ An MTZ column with
  column label and column type
| `CMtzColumnGroup <./CCP4XtalData.html#CMtzColumnGroup>`__
| `CMtzColumnGroupType <./CCP4XtalData.html#CMtzColumnGroupType>`__
| `CMtzData <./CCP4XtalData.html#CMtzData>`__ Some of the data contents
  of an MTZ file
| `CMtzDataFile <./CCP4XtalData.html#CMtzDataFile>`__ An MTZ
  experimental data file
| `CMtzDataset <./CCP4XtalData.html#CMtzDataset>`__
| `CObsDataFile <./CCP4XtalData.html#CObsDataFile>`__
| `CPhaserSolDataFile <./CCP4XtalData.html#CPhaserSolDataFile>`__
| `CPhasingGroup <./CCP4XtalData.html#CPhasingGroup>`__
| `CPhiFomColumnGroup <./CCP4XtalData.html#CPhiFomColumnGroup>`__
| `CPhsDataFile <./CCP4XtalData.html#CPhsDataFile>`__
| `CProgramColumnGroup <./CCP4XtalData.html#CProgramColumnGroup>`__ A
  group of MTZ columns required for program input
| `CProgramColumnGroup0 <./CCP4XtalData.html#CProgramColumnGroup0>`__
| `CReindexOperator <./CCP4XtalData.html#CReindexOperator>`__
| `CResolutionRange <./CCP4XtalData.html#CResolutionRange>`__
| `CRunBatchRange <./CCP4XtalData.html#CRunBatchRange>`__
| `CRunBatchRangeList <./CCP4XtalData.html#CRunBatchRangeList>`__
| `CShelxFADataFile <./CCP4XtalData.html#CShelxFADataFile>`__
| `CShelxLabel <./CCP4XtalData.html#CShelxLabel>`__
| `CSpaceGroup <./CCP4XtalData.html#CSpaceGroup>`__ A string holding the
  space group
| `CSpaceGroupCell <./CCP4XtalData.html#CSpaceGroupCell>`__ Cell space
  group and parameters
| `CUnmergedDataContent <./CCP4XtalData.html#CUnmergedDataContent>`__
| `CUnmergedDataFile <./CCP4XtalData.html#CUnmergedDataFile>`__ Handle
  MTZ, XDS and scalepack files. Allow wildcard filename
| `CUnmergedDataFileList <./CCP4XtalData.html#CUnmergedDataFileList>`__
| `CUnmergedMtzDataFile <./CCP4XtalData.html#CUnmergedMtzDataFile>`__
| `CWavelength <./CCP4XtalData.html#CWavelength>`__ Wavelength in
  Angstrom
| 

Introduction
~~~~~~~~~~~~

These Python data classes are used by by CCP4i2 and by compatible
pipelines and wrappers. It is the objective that these classes cover all
the crystallographic data and that CCP4i2 provides widgets suitable for
each data class. The classes are in modules named *CCP4Whatever* and the
classes are named *CWhatever* (the C is for CCP4 as distinct from Q for
Qt). The key fuctionality of each class is:

-  Read/write XML definition of the class to a *.def.xml* file
-  Read/write XML data contained in an instance of the class to a
   *.data.xml* file
-  Validate the data. Eg. check that data is within allowed ranges
-  Method to 'fix' bad data - this is particular to each class and never
   applied by the data class library itself

It is intended that the classes should also provide useful scientific
functionality.

All data classes are ultimately derived from *CCP4Data.CData*. There are
'simple' classes: *CBoolean*, *CFloat*, *CInt*, *CString* etc. which are
derived from *CCP4Data.CBaseData*. These can be sub-classed but the ues
of *qualifiers* is intended to minimise the need to sub-class. For
example *CInt* has the optional qualifiers *min* and *max* which can be
define the limits of allowed values. These limits are then used in the
*CInt.validate()* method. For example in Python code qualifiers for
*CCellLength* are set:

::

     class CCellLength(CFloat):
       QUALIFIERS = { 'min' : 0.0,
                      'toolTip' : 'Cell length in Angstrom' }

In this example the *toolTip* shows the other important use of
qualifiers: to provide helpful information for the GUI.

'Complex' classes contain other classes - the contents are defined by a
*CONTENTS* statement in the class:

::

   class CCell(CData):
       CONTENTS = 
                 {  'a' :     { 'class' : CCellLength },
                    'b' :     { 'class' : CCellLength },
                    'c' :     { 'class' : CCellLength },
                    'alpha' : { 'class' : CCellAngle },
                    'beta'  : { 'class' : CCellAngle },
                    'gamma' : { 'class' : CCellAngle }  }

Nothing else would be required to define a functional *CCell* class. The
CData.build() method will build the required data structure for the
class from this definition.

There is also support for lists and tables as *CList* and *CTable* which
(unlike Python lists) must have all elements of the same data type which
can be any CData class other than CList or CTable. These are both
derived from *CCollection* which handles a *subItem* which is the
definition of the type of the elements of the CList or CTable.

Containers
~~~~~~~~~~

*CCP4Container.CContainer* is a sub-class of *CCP4Data.CData* that can
hold a set of *CCP4Data.CData* objects that are defined at run time; the
definition usually comes from an XML *def* file. Typically a CContainer
will hold all of the data for one pipeline, wrapper or GUI. CContainers
can contain other CContainers and it is recommended that the container
for a pipeline or wrapper has three sub-containers called 'inputData',
'outputData' and 'parameters' to distinguish the different functions of
the data. The content of the *def* file can be defined most easily using
the *defEd* program.

--------------

Class list last updated: 15:51 20/Oct/16
