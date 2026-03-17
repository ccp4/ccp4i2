/**
 * Type definitions for CCP4i2 item class names.
 *
 * These types provide compile-time safety when working with item._class values.
 * They help catch typos and provide IDE autocomplete support.
 */

/**
 * Simple value types (CInt, CFloat, CString, CBoolean, etc.)
 */
export type SimpleItemClass =
  | "CInt"
  | "CFloat"
  | "CCellLength"
  | "CCellAngle"
  | "CWavelength"
  | "CString"
  | "CSequenceString"
  | "CSMILESString"
  | "CFilePath"
  | "COneWord"
  | "CCrystalName"
  | "CRangeSelection"
  | "CDatasetName"
  | "CAtomSelection"
  | "CBoolean";

/**
 * File data types (various file containers)
 */
export type FileItemClass =
  | "CPdbDataFile"
  | "CAsuDataFile"
  | "CImportUnmerged"
  | "CDictDataFile"
  | "CTLSDataFile"
  | "CPhaserSolDataFile"
  | "CPhaserRFileDataFile"
  | "CRefmacRestraintsDataFile"
  | "CSeqDataFile"
  | "CSeqAlignDataFile"
  | "CHhpredDataFile"
  | "CBlastDataFile"
  | "CDataFile"
  | "CUnmergedDataFile"
  | "CCootHistoryDataFile"
  | "CGenericReflDataFile"
  | "CDialsJsonFile"
  | "CDialsPickleFile"
  | "CMDLMolDataFile"
  | "CRefmacKeywordFile"
  | "CMol2DataFile"
  | "CMapDataFile"
  | "CUnmergedMtzDataFile"
  | "CPhaserTngDagFile"
  | "CMmcifDataFile"
  | "CMmcifReflDataFile";

/**
 * MTZ-related file types
 */
export type MtzItemClass =
  | "CObsDataFile"
  | "CMapCoeffsDataFile"
  | "CXmlDataFile"
  | "CPhsDataFile"
  | "CMiniMtzDataFile"
  | "CMtzDataFile"
  | "CFreeRDataFile";

/**
 * List/collection types
 */
export type ListItemClass =
  | "CList"
  | "CImportUnmergedList"
  | "CAltSpaceGroupList"
  | "CColumnGroupList"
  | "CEnsembleList"
  | "CRunBatchRangeList"
  | "CAsuContentSeqList"
  | "CAtomRefmacSelectionList"
  | "COccRefmacSelectionList"
  | "COccRelationRefmacList"
  | "CTLSRangeList";

/**
 * Container/composite types
 */
export type ContainerItemClass =
  | "CContainer"
  | "CSpaceGroupCell"
  | "CCell"
  | "CEnsemble"
  | "CFloatRange"
  | "CAsuContentSeq"
  | "CPdbEnsembleItem";

/**
 * Space group types
 */
export type SpaceGroupItemClass = "CSpaceGroup" | "CAltSpaceGroup";

/**
 * Specialized types
 */
export type SpecializedItemClass = "CReindexOperator" | "CColumnGroup" | "CRunBatchRange";

/**
 * All known item class names
 */
export type ItemClass =
  | SimpleItemClass
  | FileItemClass
  | MtzItemClass
  | ListItemClass
  | ContainerItemClass
  | SpaceGroupItemClass
  | SpecializedItemClass;

/**
 * Type guard to check if a string is a known item class
 */
export function isKnownItemClass(className: string): className is ItemClass {
  const knownClasses: ItemClass[] = [
    // Simple types
    "CInt",
    "CFloat",
    "CCellLength",
    "CCellAngle",
    "CWavelength",
    "CString",
    "CSequenceString",
    "CSMILESString",
    "CFilePath",
    "COneWord",
    "CCrystalName",
    "CRangeSelection",
    "CDatasetName",
    "CAtomSelection",
    "CBoolean",
    // File types
    "CPdbDataFile",
    "CAsuDataFile",
    "CImportUnmerged",
    "CDictDataFile",
    "CTLSDataFile",
    "CPhaserSolDataFile",
    "CPhaserRFileDataFile",
    "CRefmacRestraintsDataFile",
    "CSeqDataFile",
    "CSeqAlignDataFile",
    "CHhpredDataFile",
    "CBlastDataFile",
    "CDataFile",
    "CUnmergedDataFile",
    "CCootHistoryDataFile",
    "CGenericReflDataFile",
    "CDialsJsonFile",
    "CDialsPickleFile",
    "CMDLMolDataFile",
    "CRefmacKeywordFile",
    "CMol2DataFile",
    "CMapDataFile",
    "CUnmergedMtzDataFile",
    "CPhaserTngDagFile",
    "CMmcifDataFile",
    "CMmcifReflDataFile",
    // MTZ types
    "CObsDataFile",
    "CMapCoeffsDataFile",
    "CXmlDataFile",
    "CPhsDataFile",
    "CMiniMtzDataFile",
    "CMtzDataFile",
    "CFreeRDataFile",
    // List types
    "CList",
    "CImportUnmergedList",
    "CAltSpaceGroupList",
    "CColumnGroupList",
    "CEnsembleList",
    "CRunBatchRangeList",
    "CAsuContentSeqList",
    "CAtomRefmacSelectionList",
    "COccRefmacSelectionList",
    "COccRelationRefmacList",
    "CTLSRangeList",
    // Container types
    "CContainer",
    "CSpaceGroupCell",
    "CCell",
    "CEnsemble",
    "CFloatRange",
    "CAsuContentSeq",
    "CPdbEnsembleItem",
    // Space group types
    "CSpaceGroup",
    "CAltSpaceGroup",
    // Specialized types
    "CReindexOperator",
    "CColumnGroup",
    "CRunBatchRange",
  ];

  return knownClasses.includes(className as ItemClass);
}

/**
 * The authoritative class→component registry now lives in task-element.tsx
 * (COMPONENT_REGISTRY).  This file defines only the ItemClass types used
 * for compile-time safety across the codebase.
 */
