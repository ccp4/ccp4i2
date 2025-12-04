type CInt = Number;
type CFloat = Number;
type CBoolean = Boolean;
type CString = string;
type CDict = any;

export interface CAltSpaceGroup extends CSpaceGroup {}
export type CAltSpaceGroupList = CAltSpaceGroup[];
export type CAngle = CFloat;
export interface CAnnotation {
  text: CString;
  time: CTime;
  author: CUserId;
}
export type CAnnotationList = CAnnotation[];
export interface CAnomalousColumnGroup extends CProgramColumnGroup {}
export interface CAnomalousIntensityColumnGroup extends CProgramColumnGroup {}
export interface CAnomalousScatteringElement extends CElement {}
export interface CAsuComponent {
  moleculeType: CString;
  seqFile: CSeqDataFile;
  numberOfCopies: CInt;
}
export type CAsuComponentList = CAsuComponent[];
export interface CAsuContent extends CDataFileContent {
  seqList: CAsuContentSeqList;
}
export interface CAsuContentSeq {
  sequence: CSequenceString;
  nCopies: CInt;
  polymerType: CString;
  name: CString;
  description: CString;
  source: CDataFile;
}
export type CAsuContentSeqList = CAsuContentSeq[];
export interface CAsuDataFile extends CI2XmlDataFile {
  selection: CDict;
}
export interface CAtomCountPerformance extends CPerformanceIndicator {
  nAtoms: CInt;
  nResidues: CInt;
}
export interface CAtomRefmacSelection {
  groupId: CInt;
  chainId: COneWord;
  firstRes: CInt;
  lastRes: CInt;
}
export interface CAtomRefmacSelectionGroups {
  groupIds: CString;
}
export type CAtomRefmacSelectionList = CAtomRefmacSelection[];
export interface CAtomRefmacSelectionOccupancy {
  groupId: CInt;
  chainIds: CString;
  firstRes: CInt;
  lastRes: CInt;
  atoms: CString;
  alt: COneWord;
}
export interface CAtomSelection {
  text: CString;
}
export type CAuthor = CString;
export interface CBaseData {}
export interface CBibReference {
  pmid: CInt;
  title: CString;
  authorList: CAuthor[];
  source: CString;
  url: CString;
  selected: CBoolean;
}
export interface CBibReferenceGroup {
  taskName: CString;
  version: CString;
  title: CString;
  references: CBibReference[];
}
export interface CBlastData extends CDataFileContent {
  queryId: CString;
  alignmentList: CBlastItem[];
}
type CBlastDataFile = CDataFile;
export interface CBlastItem {
  hitId: CString;
  querySequence: CString;
  hitSequence: CString;
}
export interface CCell {
  a: CCellLength;
  b: CCellLength;
  c: CCellLength;
  alpha: CCellAngle;
  beta: CCellAngle;
  gamma: CCellAngle;
}
export type CCellAngle = CFloat;
export type CCellLength = CFloat;
export interface CChemComp {
  id: COneWord;
  three_letter_code: COneWord;
  name: CString;
  group: CString;
  number_atoms_all: CInt;
  number_atoms_nh: CInt;
  desc_level: CInt;
}
export interface CCollection {}
export interface CColumnGroup {
  columnGroupType: COneWord;
  contentFlag: CInt;
  dataset: CString;
  columnList: CMtzColumn[];
  selected: CBoolean;
}
export interface CColumnGroupItem {
  columnName: COneWord;
  defaultList: CString;
  columnType: CColumnTypeList;
  partnerTo: COneWord;
  partnerOffset: CInt;
}
export type CColumnGroupList = CColumnGroup[];
export type CColumnType = CString;
export type CColumnTypeList = CColumnType[];
export interface CContainer {}
type CCootHistoryDataFile = CDataFile;
export type CCrystalName = CString;
export interface CCustomComFile {
  text: CString;
  name: CString;
}
export type CCustomComFileList = CCustomComFile[];
export interface CCustomTaskDefinition extends CContainer {
  name: COneWord;
  title: CString;
  tagCharacter: CString;
  comLine: CString;
  comFileList: CCustomComFileList;
  paramList: CCustomTaskParamList;
}
export type CCustomTaskFileFunction = CString;
export interface CCustomTaskParam {
  name: COneWord;
  dataType: CI2DataType;
  label: CString;
  obligatory: CBoolean;
  saveDataToDb: CBoolean;
  function: CCustomTaskFileFunction;
  mergeTo: CString;
  splitColumns: CString;
  requiredContentType: CBoolean[];
  outputFilePath: CString;
}
export type CCustomTaskParamList = CCustomTaskParam[];
export interface CData extends CDataQualifiers {}
export interface CDataFile {
  project: CProjectId;
  baseName: CFilePath;
  relPath: CFilePath;
  annotation: CString;
  dbFileId: CUUID;
  subType: CInt;
  contentFlag: CInt;
}
export interface CDataFileContent {}
export interface CDataReductionCCPerformance extends CPerformanceIndicator {
  spaceGroup: CSpaceGroup;
  highResLimit: CFloat;
  ccHalf: CFloat;
}
export interface CDataReductionPerformance extends CPerformanceIndicator {
  spaceGroup: CSpaceGroup;
  highResLimit: CFloat;
  rMeas: CFloat;
}
type CDataReflFile = CDataFile;
export interface CDataset {
  selected: CBoolean;
  obsDataFile: CObsDataFile;
  crystalName: CCrystalName;
  datasetName: CDatasetName;
  formFactors: CFormFactor;
  formFactorSource: CString;
}
export type CDatasetList = CDataset[];
export type CDatasetName = CString;
export interface CDateRange {
  year: CInt;
  month: CString;
  day: CInt;
  yearRange: CInt;
  monthRange: CInt;
  dayRange: CInt;
}
type CDialsJsonFile = CDataFile;
type CDialsPickleFile = CDataFile;
export interface CDictData {
  monomerList: CChemComp[];
}
type CDictDataFile = CDataFile;
export interface CEBIValidationXMLDataFile extends CXmlDataFile {}
export interface CElement extends COneWord {}
export interface CEnsemble {
  label: COneWord;
  number: CInt;
  use: CBoolean;
  pdbItemList: CPdbEnsembleItem[];
}
export type CEnsembleList = CEnsemble[];
export interface CEnsemblePdbDataFile extends CPdbDataFile {
  selection: CAtomSelection;
}
export interface CEulerRotation {
  alpha: CAngle;
  beta: CAngle;
  gamma: CAngle;
}
export interface CExePath {
  exeName: CString;
  exePath: CDataFile;
}
export type CExePathList = CExePath[];
export interface CExpPhasPerformance extends CPerformanceIndicator {
  FOM: CFloat;
  CFOM: CFloat;
  Hand1Score: CFloat;
  Hand2Score: CFloat;
  CC: CFloat;
  RFactor: CFloat;
  RFree: CFloat;
  annotation: CString;
}
export type CExperimentalDataType = CString;
export interface CExportedFile {
  exportId: CUUID;
}
export type CExportedFileList = CExportedFile[];
export interface CFPairColumnGroup extends CProgramColumnGroup {}
export interface CFSigFColumnGroup extends CProgramColumnGroup {}
export type CFileFunction = CString;
export type CFilePath = CString;
export interface CFloatRange extends CRange {
  start: CFloat;
  end: CFloat;
}
export interface CFollowFromJob extends CUUID {}
export interface CFont {
  family: CString;
  style: CInt;
  pointSize: CInt;
  weight: CInt;
}
export interface CFormFactor {
  Fp: CFloat;
  Fpp: CFloat;
}
export interface CFreeRColumnGroup extends CProgramColumnGroup {}
export interface CFreeRDataFile extends CMiniMtzDataFile {}
type CGenericReflDataFile = CDataFile;
export interface CHLColumnGroup extends CProgramColumnGroup {}
export interface CHhpredData extends CDataFileContent {
  alignmentList: CHhpredItem[];
}
type CHhpredDataFile = CDataFile;
export interface CHhpredItem {
  annotation: CString;
  identifier: CString;
  chain: CString;
}
export type CHostName = CString;
export interface CHostname extends CHostName {}
export type CI2DataType = CString;
export interface CI2XmlDataFile extends CXmlDataFile {
  header: CI2XmlHeader;
}
export interface CI2XmlHeader {
  function: CFileFunction;
  userId: CUserId;
  hostName: CHostName;
  creationTime: CTime;
  ccp4iVersion: CVersion;
  pluginName: CString;
  pluginVersion: CVersion;
  pluginTitle: CString;
  projectName: CProjectName;
  projectId: CProjectId;
  jobId: CUUID;
  jobNumber: CString;
  comment: CString;
  OS: CString;
}
export interface CIPairColumnGroup extends CProgramColumnGroup {}
export interface CISigIColumnGroup extends CProgramColumnGroup {}
type CImageFile = CDataFile;
export type CImageFileList = CImageFile[];
type CImosflmXmlDataFile = CDataFile;
export interface CImportUnmerged {
  file: CUnmergedDataFile;
  cell: CCell;
  wavelength: CWavelength;
  crystalName: CString;
  dataset: CString;
  excludeSelection: CRangeSelection;
}
export type CImportUnmergedList = CImportUnmerged[];
export interface CImportedJobData {
  name: COneWord;
  dataType: CI2DataType;
  label: CString;
  fileName: CDataFile;
}
export type CImportedJobDataList = CImportedJobData[];
export interface CImportedJobDefinition extends CContainer {
  name: COneWord;
  title: CString;
  ifImportDir: CBoolean;
  commandFile: CDataFile;
  logFile: CDataFile;
  inputFileDefinitionList: CImportedJobDataList;
  outputFileDefinitionList: CImportedJobDataList;
}
export interface CIntRange extends CRange {
  start: CInt;
  end: CInt;
}
export type CJobStatus = CInt;
export type CJobTitle = CString;
export interface CList extends CCollection {}
type CMDLMolDataFile = CDataFile;
export interface CMapCoeffsDataFile extends CMiniMtzDataFile {}
export interface CMapColumnGroup extends CProgramColumnGroup {}
type CMapDataFile = CDataFile;
export interface CMatrix33 {}
export interface CMergeMiniMtz {
  fileName: CMiniMtzDataFile;
  columnTag: CString;
  columnNames: CString;
}
export type CMergeMiniMtzList = CMergeMiniMtz[];
export interface CMetaDataTag {
  tag: CString;
}
export type CMetaDataTagList = CMetaDataTag[];
export interface CMiniMtzDataFile extends CMtzDataFile {}
export type CMiniMtzDataFileList = CMiniMtzDataFile[];
export interface CMmcifData extends CDataFileContent {}
type CMmcifDataFile = CDataFile;
export interface CMmcifReflData extends CMmcifData {
  cell: CCell;
  spaceGroup: CSpaceGroup;
  wavelength: CWavelength;
  haveFreeRColumn: CBoolean;
  haveFobsColumn: CBoolean;
  haveFpmObsColumn: CBoolean;
  haveIobsColumn: CBoolean;
  haveIpmObsColumn: CBoolean;
}
export interface CMmcifReflDataFile extends CMmcifDataFile {}
export interface CModelBuildPerformance extends CPerformanceIndicator {
  RFactor: CFloat;
  completeness: CFloat;
  annotation: CString;
}
export interface CMonomer {
  identifier: CString;
  formula: CString;
  dictionaryName: CString;
  smiles: CString;
}
export interface CMtzColumn {
  columnLabel: COneWord;
  columnType: CColumnType;
  dataset: COneWord;
  groupIndex: CInt;
}
export interface CMtzColumnGroup {
  groupType: CMtzColumnGroupType;
  columns: CMtzColumn[];
}
export interface CMtzColumnGroupType extends CColumnType {}
export interface CMtzData extends CDataFileContent {
  cell: CCell;
  spaceGroup: CSpaceGroup;
  resolutionRange: CResolutionRange;
  listOfColumns: CMtzColumn[];
  datasets: CString[];
  crystalNames: CString[];
  wavelengths: CWavelength[];
  datasetCells: CCell[];
  merged: CBoolean;
}
type CMtzDataFile = CDataFile;
export interface CMtzDataset {
  name: CString;
  columnGroups: CMtzColumnGroup[];
}
export interface CObsDataFile extends CMiniMtzDataFile {}
export type COccRefmacSelectionList = CAtomRefmacSelectionOccupancy[];
export type COccRelationRefmacList = CAtomRefmacSelectionGroups[];
export type COneWord = CString;
export type COutputFileList = CString[];
type CPDFDataFile = CDataFile;
export interface CPairefPerformance extends CPerformanceIndicator {
  cutoff: CFloat;
}
export interface CPatchDefinition extends CContainer {
  taskNameList: any[];
  projectId: CUUID;
  jobId: CUUID;
  text1: CString;
  text2: CString;
  diffs: CString;
  controlParameters: CContainer;
}
export interface CPatchSelection {
  taskName: CString;
  patch: CString;
}
export interface CPdbData extends CDataFileContent {}
export interface CPdbDataFile extends CDataFile {
  selection: CAtomSelection;
}
export type CPdbDataFileList = CPdbDataFile[];
export interface CPdbEnsembleItem {
  structure: CPdbDataFile;
  identity_to_target: CFloat;
  rms_to_target: CFloat;
}
export interface CPerformanceIndicator {
  value: CFloat;
  annotation: CString;
}
export interface CPhaseErrorPerformance extends CPerformanceIndicator {
  phaseError: CFloat;
  weightedPhaseError: CFloat;
  reflectionCorrelation: CFloat;
}
type CPhaserRFileDataFile = CDataFile;
type CPhaserSolDataFile = CDataFile;
export interface CPhiFomColumnGroup extends CProgramColumnGroup {}
export interface CPhsDataFile extends CMiniMtzDataFile {}
type CPostscriptDataFile = CDataFile;
export interface CPreferences extends CContainer {
  TASK_WINDOW_LAYOUT: CString;
  JOB_LIST_DATE_TIME: CBoolean;
  TABLES_ALTERNATING_COLOR: CBoolean;
  AUTO_INFO_ON_FILE_IMPORT: CBoolean;
  EXTERNAL_FILES_IN_EXTERNAL_BROWSER: CBoolean;
  EXTERNAL_FILES_IN_IFRAME: CBoolean;
  HD_ICONS: CBoolean;
  COMPACT_TASK_MENU: CBoolean;
  GUI_FONT_SIZE: CInt;
  REPORT_ZOOM_FACTOR: CFloat;
  BROWSER_ZOOM_FACTOR: CFloat;
  REPORT_FONT_SIZE: CInt;
  INVALID_FRAME_MODE: CInt;
  INVALID_FRAME_WIDTH: CInt;
  INVALID_FRAME_COLOUR: CString;
  WINDOWS_STYLE: CString;
  PROJECT_VIEWER_BUTTONS: CString;
  TOOLBARBUTTONSSTYLE: CInt;
  TEXT_VIEW_LINE_LIMIT: CInt;
  COOT_EXECUTABLE: CDataFile;
  CCP4MG_EXECUTABLE: CDataFile;
  SHELXDIR: CDataFile;
  DIALSDIR: CDataFile;
  BUSTERDIR: CDataFile;
  EXEPATHLIST: CExePathList;
  BZR_DOWNLOAD: CBoolean;
  DBLOCAL_QUIT_RUNNING: CBoolean;
  DELETE_INTERACTIVE_JOBS: CBoolean;
  SHOW_DELETE_INTERACTIVE_JOBS: CBoolean;
  FILESYSTEMWATCHERPOLLER: CBoolean;
  RETAIN_DIAGNOSTIC_FILES: CBoolean;
  NATIVEFILEBROWSER: CBoolean;
  SHOW_TASK_MENU_TOOLBUTTON: CBoolean;
  SHOW_JOB_SEARCH_TOOLBUTTON: CBoolean;
  SHOW_EXPORT_PROJECT_TOOLBUTTON: CBoolean;
  SHOW_RUN_TOOLBUTTON: CBoolean;
  SHOW_RUN_REMOTE_TOOLBUTTON: CBoolean;
  SHOW_CLONE_TOOLBUTTON: CBoolean;
  SHOW_TASK_HELP_TOOLBUTTON: CBoolean;
  SHOW_REFERENCES_TOOLBUTTON: CBoolean;
  SHOW_EXPORT_MTZ_TOOLBUTTON: CBoolean;
  SHOW_VIEW_COOT_TOOLBUTTON: CBoolean;
  SHOW_VIEW_CCP4MG_TOOLBUTTON: CBoolean;
  SHOW_SHOW_LOG_TOOLBUTTON: CBoolean;
  SHOW_SHOW_I2RUN_TOOLBUTTON: CBoolean;
  NEW_PROJECT_TOOLBUTTON: CBoolean;
  SHOW_TIPS_OF_THE_DAY: CBoolean;
  SHOW_TASK_MODE_BUTTONS: CBoolean;
  RESTORE_TO_TASKLIST: CBoolean;
  AUTO_UPDATE_REPORT80: CBoolean;
  SHOW_WRAPPERS: CBoolean;
  DISABLE_WEBGL: CBoolean;
  PDB_REDO_TOKEN_ID: CString;
  PDB_REDO_TOKEN_SECRET: CString;
  PDB_REDO_TOKEN_EXPIRES: CString;
}
export interface CProgramColumnGroup {}
export interface CProgramColumnGroup0 {
  columnGroup: CMtzColumnGroup;
  datasetName: CString;
}
export interface CProjectId extends CUUID {}
export type CProjectName = CString;
export interface CRange {}
export type CRangeSelection = CString;
export interface CRefinementPerformance extends CPerformanceIndicator {
  RFactor: CFloat;
  RFree: CFloat;
  RMSBond: CFloat;
  RMSAngle: CFloat;
  weightUsed: CFloat;
  annotation: CString;
}
export interface CRefmacAnomalousAtom {
  atomType: CString;
  Fp: CFloat;
  Fpp: CFloat;
}
type CRefmacKeywordFile = CDataFile;
type CRefmacRestraintsDataFile = CDataFile;
export interface CRefmacRigidGroupItem {
  rigid_group_id: CString;
  segmentList: CRefmacRigidGroupSegment[];
}
export type CRefmacRigidGroupList = CRefmacRigidGroupItem[];
export interface CRefmacRigidGroupSegment {
  chain_id: CString;
  residue_1: CInt;
  residue_2: CInt;
}
export interface CReindexOperator {
  h: CString;
  k: CString;
  l: CString;
}
export interface CResidueRange {
  chainId: COneWord;
  firstRes: COneWord;
  lastRes: COneWord;
}
export type CResidueRangeList = CResidueRange[];
export interface CResolutionRange {
  low: CFloat;
  high: CFloat;
}
export interface CRunBatchRange {
  runNumber: CInt;
  batchRange0: CInt;
  batchRange1: CInt;
  resolution: CFloat;
  fileNumber: CInt;
}
export type CRunBatchRangeList = CRunBatchRange[];
type CSceneDataFile = CDataFile;
export interface CSearchPath {
  name: CString;
  path: CDataFile;
}
export type CSearchPathList = CSearchPath[];
type CSeqAlignDataFile = CDataFile;
type CSeqDataFile = CDataFile;
export type CSeqDataFileList = CSeqDataFile[];
export interface CSequence extends CBioPythonSeqInterface {
  identifier: CString;
  referenceDb: CString;
  reference: CString;
  name: CString;
  description: CString;
  sequence: CString;
  moleculeType: CString;
}
export interface CSequenceAlignment extends CBioPythonSeqInterface {
  identifier: CString;
  moleculeType: CString;
}
export interface CSequenceMeta {
  uniprotId: CString;
  organism: CString;
  expressionSystem: CString;
}
export type CSequenceString = CString;
export interface CServerGroup {
  name: CString;
  mechanism: CString;
  serverList: CHostname[];
  userExtensible: CBoolean;
  customCodeFile: CDataFile;
  queueOptionsFile: CDataFile;
  ccp4Dir: CString;
  tempDir: CString;
  sge_root: CString;
  keyFilename: CString;
  validate: CString;
  timeout: CFloat;
  maxTries: CInt;
}
type CShelxFADataFile = CDataFile;
export type CShelxLabel = CString;
export type CSpaceGroup = CString;
export interface CSpaceGroupCell {
  spaceGroup: CSpaceGroup;
  cell: CCell;
}
export interface CSuperposePerformance extends CPerformanceIndicator {
  RMSxyz: CFloat;
  nResidues: CInt;
}
type CTLSDataFile = CDataFile;
export interface CTestObsConversionsPerformance extends CPerformanceIndicator {
  columnLabelsString: CString;
}
type CTextDataFile = CDataFile;
export type CTime = CInt;
export interface CTransformation {
  translation: CXyz;
  rotation: CEulerRotation;
}
export type CUUID = CString;
export interface CUnmergedDataContent extends CDataFileContent {
  format: CString;
  merged: CString;
  crystalName: CCrystalName;
  datasetName: CDatasetName;
  cell: CCell;
  spaceGroup: CSpaceGroup;
  batchs: CString;
  lowRes: CFloat;
  highRes: CFloat;
  knowncell: CBoolean;
  knownwavelength: CBoolean;
  numberLattices: CInt;
  wavelength: CWavelength;
  numberofdatasets: CInt;
}
type CUnmergedDataFile = CDataFile;
export type CUnmergedDataFileList = CUnmergedDataFile[];
export interface CUnmergedMtzDataFile extends CMtzDataFile {}
export interface CUserAddress {
  platformNode: CString;
  userId: CUserId;
}
export type CUserId = CString;
export type CVersion = CString;
export type CWavelength = CFloat;
export type CWorkflowContainerList = CDict;
export interface CWorkflowDataFlow {
  fromJob: CString;
  fromKey: CString;
  toKey: CString;
  annotation: CString;
}
export type CWorkflowDataFlowList = CWorkflowDataFlow[];
export interface CWorkflowDefinition extends CContainer {
  jobDef: CWorkflowJobDefinitionDict;
  title: CString;
}
export interface CWorkflowFileOut {
  key: CString;
  className: CString;
}
export interface CWorkflowJobDefinition {
  taskName: CString;
  input: CWorkflowDataFlowList;
  allOutputFiles: CWorkflowFileOut[];
  output: CWorkflowDataFlowList;
}
export type CWorkflowJobDefinitionDict = CDict;
export interface CXia2ImageSelection {
  imageFile: CImageFile;
  imageStart: CInt;
  imageEnd: CInt;
}
export type CXia2ImageSelectionList = CXia2ImageSelection[];
type CXmgrDataFile = CDataFile;
type CXmlDataFile = CDataFile;
export interface CXyz {
  x: CFloat;
  y: CFloat;
  z: CFloat;
}
export interface CXyzBox {
  xMin: CFloat;
  yMin: CFloat;
  zMin: CFloat;
  xMax: CFloat;
  yMax: CFloat;
  zMax: CFloat;
}
type CYmlFile = CDataFile;
