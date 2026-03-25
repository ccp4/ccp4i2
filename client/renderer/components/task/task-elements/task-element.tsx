import { Job } from "../../../types/models";
import { CIntElement } from "./cint";
import { PropsWithChildren, useMemo } from "react";
import { SxProps, Theme, Typography } from "@mui/material";
import { CStringElement } from "./cstring";
import { useJob } from "../../../utils";
import { CFloatElement } from "./cfloat";
import { CPdbDataFileElement } from "./cpdbdatafile";
import { CAsuDataFileElement } from "./casudatafile";
import { CFreeRDataFileElement } from "./cfreerfile";
import { CMiniMtzDataFileElement } from "./cminimtzdatafile";
import { CBooleanElement } from "./cboolean";
import { CListElement } from "./clist";
// Lazy import to break circular dependency: ccontainer.tsx imports task-element.tsx
let _CCP4i2ContainerElement: React.FC<any> | null = null;
const getLazyContainerElement = (): React.FC<any> => {
  if (!_CCP4i2ContainerElement) {
    _CCP4i2ContainerElement =
      require("./ccontainer").CCP4i2ContainerElement;
  }
  return _CCP4i2ContainerElement!;
};
import { CImportUnmergedElement } from "./cimportunmerged";
import { CCellElement } from "./ccell";
import { CEnsembleElement } from "./censemble";
import { CAltSpaceGroupElement } from "./caltspacegroupelement";
import { CSimpleDataFileElement } from "./csimpledatafile";
import { CReindexOperatorElement } from "./creindexoperator";
import { CRunBatchRangeElement } from "./crunbatchrange";
import { CRangeElement } from "./crange";
import { v4 as uuid4 } from "uuid";
import { CAsuContentSeqElement } from "./casucontentseq";
import { CColumnGroupElement } from "./ccolumngroup";
import { CPdbEnsembleItemElement } from "./cpdbensembleitem";
import { CSMILESStringElement } from "./csmilesstring";
import { CAsuContentSeqListElement } from "./casucontentseqlist";
import {
  CAtomRefmacSelectionListElement,
  COccRefmacSelectionListElement,
  COccRelationRefmacListElement,
} from "./coccupancygroups";
import { CTLSRangeListElement } from "./ctlsranges";
import { CAtomSelectionElement } from "./catomselection";
import { useInferredVisibility } from "./hooks/useInferredVisibility";
import type { ItemClass } from "./types/item-classes";

export interface CCP4i2TaskElementProps extends PropsWithChildren {
  job: Job;
  itemName: string;
  sx?: SxProps<Theme>;
  pathOfItem?: (item: HTMLElement) => string;
  visibility?: boolean | (() => boolean);
  disabled?: boolean | (() => boolean);
  qualifiers?: any;
  onChange?: (updatedItem: any) => void;
  suppressMutations?: boolean;
}

/**
 * Registry entry: a React component and optional qualifier overrides.
 * When extraQualifiers is set, it is merged on top of the resolved qualifiers.
 */
interface RegistryEntry {
  component: React.FC<any>;
  extraQualifiers?: Record<string, any>;
}

/**
 * Maps every known ItemClass to its rendering component (and optional extra
 * qualifiers).  Adding a new item type is a one-liner here — no switch/case
 * needed.
 */
const COMPONENT_REGISTRY: Record<string, RegistryEntry> = {
  // Simple value types
  CInt: { component: CIntElement },
  CFloat: { component: CFloatElement },
  CCellLength: { component: CFloatElement },
  CCellAngle: { component: CFloatElement },
  CWavelength: { component: CFloatElement },
  CString: { component: CStringElement },
  CSequenceString: { component: CStringElement },
  CFilePath: { component: CStringElement },
  COneWord: { component: CStringElement },
  CCrystalName: { component: CStringElement },
  CRangeSelection: { component: CStringElement },
  CDatasetName: { component: CStringElement },
  CAtomSelection: { component: CAtomSelectionElement },
  CSMILESString: { component: CSMILESStringElement },
  CBoolean: { component: CBooleanElement },

  // File types → CSimpleDataFileElement
  CPdbDataFile: { component: CPdbDataFileElement },
  CImportUnmerged: { component: CImportUnmergedElement },
  CAsuDataFile: { component: CAsuDataFileElement },
  CDictDataFile: { component: CSimpleDataFileElement },
  CTLSDataFile: { component: CSimpleDataFileElement },
  CPhaserSolDataFile: { component: CSimpleDataFileElement },
  CPhaserRFileDataFile: { component: CSimpleDataFileElement },
  CRefmacRestraintsDataFile: { component: CSimpleDataFileElement },
  CSeqDataFile: { component: CSimpleDataFileElement },
  CSeqAlignDataFile: { component: CSimpleDataFileElement },
  CHhpredDataFile: { component: CSimpleDataFileElement },
  CBlastDataFile: { component: CSimpleDataFileElement },
  CDataFile: { component: CSimpleDataFileElement },
  CUnmergedDataFile: { component: CSimpleDataFileElement },
  CCootHistoryDataFile: { component: CSimpleDataFileElement },
  CDialsJsonFile: { component: CSimpleDataFileElement },
  CDialsPickleFile: { component: CSimpleDataFileElement },
  CMDLMolDataFile: { component: CSimpleDataFileElement },
  CRefmacKeywordFile: { component: CSimpleDataFileElement },
  CMol2DataFile: { component: CSimpleDataFileElement },
  CMapDataFile: { component: CSimpleDataFileElement },
  CUnmergedMtzDataFile: { component: CSimpleDataFileElement },
  CPhaserTngDagFile: { component: CSimpleDataFileElement },
  CMmcifDataFile: { component: CSimpleDataFileElement },
  CMmcifReflDataFile: { component: CSimpleDataFileElement },
  // CGenericReflDataFile was duplicated in the old switch (first match won → CSimpleDataFileElement)
  CGenericReflDataFile: { component: CSimpleDataFileElement },

  // MTZ reflection file types
  CObsDataFile: { component: CMiniMtzDataFileElement },
  CMapCoeffsDataFile: { component: CMiniMtzDataFileElement },
  CXmlDataFile: { component: CMiniMtzDataFileElement },
  CPhsDataFile: { component: CMiniMtzDataFileElement },
  CMiniMtzDataFile: {
    component: CMiniMtzDataFileElement,
    extraQualifiers: {
      mimeTypeName: [
        "application/CCP4-mtz",
        "application/CCP4-mtz-observed",
        "application/CCP4-mtz-mini",
        "application/CCP4-mtz-phases",
        "application/CCP4-mtz-freerflag",
      ],
    },
  },
  CMtzDataFile: {
    component: CMiniMtzDataFileElement,
    extraQualifiers: {
      mimeTypeName: [
        "application/CCP4-mtz",
        "application/CCP4-mtz-observed",
        "application/CCP4-mtz-mini",
        "application/CCP4-mtz-phases",
        "application/CCP4-mtz-freerflag",
      ],
    },
  },
  CFreeRDataFile: { component: CFreeRDataFileElement },

  // List / collection types
  CList: { component: CListElement },
  CImportUnmergedList: { component: CListElement },
  CAltSpaceGroupList: { component: CListElement },
  CColumnGroupList: { component: CListElement },
  CEnsembleList: { component: CListElement },
  CRunBatchRangeList: { component: CListElement },
  CAsuContentSeqList: { component: CAsuContentSeqListElement },
  CAtomRefmacSelectionList: { component: CAtomRefmacSelectionListElement },
  COccRefmacSelectionList: { component: COccRefmacSelectionListElement },
  COccRelationRefmacList: { component: COccRelationRefmacListElement },
  CTLSRangeList: { component: CTLSRangeListElement },

  // Container / composite types (use lazy getter to break circular dependency)
  CSpaceGroupCell: { get component() { return getLazyContainerElement(); } },
  CContainer: { get component() { return getLazyContainerElement(); } },
  CCell: { component: CCellElement },
  CEnsemble: { component: CEnsembleElement },
  CFloatRange: { component: CRangeElement },
  CAsuContentSeq: { component: CAsuContentSeqElement },
  CPdbEnsembleItem: { component: CPdbEnsembleItemElement },

  // Space group types
  CSpaceGroup: { component: CAltSpaceGroupElement },
  CAltSpaceGroup: { component: CAltSpaceGroupElement },

  // Specialized types
  CReindexOperator: { component: CReindexOperatorElement },
  CColumnGroup: { component: CColumnGroupElement },
  CRunBatchRange: { component: CRunBatchRangeElement },
};

/**
 * CCP4i2TaskElement renders the appropriate widget for a given task parameter
 * based on its item class, looked up in COMPONENT_REGISTRY.
 */
export const CCP4i2TaskElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(props.itemName);

  const isVisible = useInferredVisibility(props.visibility);
  // Use a stable key based on the item's object path instead of generating a new UUID on each render
  // A new UUID on each render causes the component to unmount/remount, losing all state
  const the_uuid = useMemo(
    () => item?._objectPath || props.itemName || uuid4(),
    [item?._objectPath, props.itemName]
  );

  const qualifiers = useMemo<any>(() => {
    if (item?._qualifiers) {
      try {
        const overriddenQualifiers = props.qualifiers
          ? { ...item._qualifiers, ...props.qualifiers }
          : item._qualifiers;
        return overriddenQualifiers;
      } catch (err) {
        console.warn(`Error getting qualifiers on ${props.itemName}`);
      }
    }
    return props.qualifiers;
  }, [item, props.qualifiers]);

  const interfaceElement = useMemo(() => {
    const itemClass = item?._class as ItemClass | undefined;
    const entry = itemClass ? COMPONENT_REGISTRY[itemClass] : undefined;

    if (!entry) {
      return (
        <Typography key={the_uuid}>
          {item ? item._class : "No item"}
        </Typography>
      );
    }

    const Component = entry.component;
    const mergedQualifiers = entry.extraQualifiers
      ? { ...qualifiers, ...entry.extraQualifiers }
      : qualifiers;

    return <Component key={the_uuid} {...props} qualifiers={mergedQualifiers} />;
  }, [item, qualifiers]);

  return isVisible ? <>{interfaceElement}</> : null;
};
