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
  | "CMol2DataFile";

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
  | "CAsuContentSeqList";

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
export type SpecializedItemClass = "CReindexOperator" | "CColumnGroup";

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
    "CAsuContentSeqList",
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
  ];

  return knownClasses.includes(className as ItemClass);
}

/**
 * Mapping of item classes to their corresponding React component element type.
 * This helps document which component handles which class.
 */
export const ITEM_CLASS_COMPONENT_MAP = {
  // Simple types -> CSimpleElement
  CInt: "CIntElement",
  CFloat: "CFloatElement",
  CCellLength: "CFloatElement",
  CCellAngle: "CFloatElement",
  CWavelength: "CFloatElement",
  CString: "CStringElement",
  CSequenceString: "CStringElement",
  CSMILESString: "CSMILESStringElement",
  CFilePath: "CStringElement",
  COneWord: "CStringElement",
  CCrystalName: "CStringElement",
  CRangeSelection: "CStringElement",
  CDatasetName: "CStringElement",
  CAtomSelection: "CStringElement",
  CBoolean: "CBooleanElement",

  // File types -> CSimpleDataFileElement
  CPdbDataFile: "CPdbDataFileElement",
  CAsuDataFile: "CAsuDataFileElement",
  CImportUnmerged: "CImportUnmergedElement",
  CDictDataFile: "CSimpleDataFileElement",
  CTLSDataFile: "CSimpleDataFileElement",
  // ... (truncated for brevity - maps all file types)

  // MTZ types -> CMiniMtzDataFileElement
  CObsDataFile: "CMiniMtzDataFileElement",
  CMapCoeffsDataFile: "CMiniMtzDataFileElement",
  CMiniMtzDataFile: "CMiniMtzDataFileElement",
  CMtzDataFile: "CMiniMtzDataFileElement",
  CFreeRDataFile: "CMiniMtzDataFileElement",

  // List types -> CListElement
  CList: "CListElement",
  CColumnGroupList: "CListElement",
  CAltSpaceGroupList: "CListElement",
  CEnsembleList: "CListElement",
  CAsuContentSeqList: "CAsuContentSeqListElement",

  // Container types
  CContainer: "CCP4i2ContainerElement",
  CSpaceGroupCell: "CCP4i2ContainerElement",
  CCell: "CCellElement",
  CEnsemble: "CEnsembleElement",

  // Space group types
  CSpaceGroup: "CAltSpaceGroupElement",
  CAltSpaceGroup: "CAltSpaceGroupElement",

  // Specialized types
  CReindexOperator: "CReindexOperatorElement",
  CColumnGroup: "CColumnGroupElement",
} as const;

export type ItemClassComponentMap = typeof ITEM_CLASS_COMPONENT_MAP;
