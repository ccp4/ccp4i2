import { useCallback, useEffect, useMemo, useState } from "react";
import { LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { apiJson } from "../../../api-fetch";

// Layout constants
const LABEL_WIDTH = "14rem";
const SHORT_FIELD = "10rem";
const DROPDOWN_FIELD = "16rem";

// Types
interface MonomerBond {
  atom1: string;
  atom2: string;
  type: string;
}

interface MonomerInfo {
  atoms: string[];
  bonds: MonomerBond[];
}

/** Dict keyed by monomer code -> { atoms, bonds } */
type MonomerDict = Record<string, MonomerInfo>;

const EMPTY_MONOMER: MonomerInfo = { atoms: [], bonds: [] };

/**
 * Fetch atom/bond info for a monomer from the CCP4 library.
 * Returns { atoms, bonds } or the empty default on failure.
 */
function useMonomerInfo(resName: string | undefined): MonomerInfo {
  const [info, setInfo] = useState<MonomerInfo>(EMPTY_MONOMER);

  useEffect(() => {
    if (!resName || resName.trim().length === 0) {
      setInfo(EMPTY_MONOMER);
      return;
    }

    let cancelled = false;
    const code = resName.trim().toUpperCase();
    apiJson<{ success: boolean; data?: { atoms: string[]; bonds: MonomerBond[] } }>(
      `monomer-info/${code}`
    )
      .then((res) => {
        if (!cancelled && res.success && res.data) {
          setInfo({ atoms: res.data.atoms, bonds: res.data.bonds });
        }
      })
      .catch((err) => {
        console.warn(`[MakeLink] Fetch failed for "${code}":`, err);
        if (!cancelled) setInfo(EMPTY_MONOMER);
      });

    return () => {
      cancelled = true;
    };
  }, [resName]);

  return info;
}

/** Format bonds as "ATOM1-ATOM2" labels for the bond-change autocomplete. */
function useBondLabels(bonds: MonomerBond[] | undefined): string[] {
  return useMemo(
    () => (bonds ?? []).map((b) => `${b.atom1}-${b.atom2}`),
    [bonds]
  );
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container, setParameter, fetchDigest } = useJob(props.job.id);
  const { value: MON_1_TYPE } = useTaskItem("MON_1_TYPE");
  const { value: MON_2_TYPE } = useTaskItem("MON_2_TYPE");
  const { value: RES_NAME_1_TLC } = useTaskItem("RES_NAME_1_TLC");
  const { value: RES_NAME_2_TLC } = useTaskItem("RES_NAME_2_TLC");
  const { value: RES_NAME_1_CIF, item: resNameItem1CIF } = useTaskItem("RES_NAME_1_CIF");
  const { value: RES_NAME_2_CIF, item: resNameItem2CIF } = useTaskItem("RES_NAME_2_CIF");
  const { value: ATOM_NAME_1_TLC, item: atomItem1TLC } = useTaskItem("ATOM_NAME_1_TLC");
  const { value: ATOM_NAME_1_CIF, item: atomItem1CIF } = useTaskItem("ATOM_NAME_1_CIF");
  const { value: ATOM_NAME_2_TLC, item: atomItem2TLC } = useTaskItem("ATOM_NAME_2_TLC");
  const { value: ATOM_NAME_2_CIF, item: atomItem2CIF } = useTaskItem("ATOM_NAME_2_CIF");
  const { value: ATOM_NAME_1, item: atomItem1Plain } = useTaskItem("ATOM_NAME_1");
  const { value: ATOM_NAME_2, item: atomItem2Plain } = useTaskItem("ATOM_NAME_2");
  const { item: dict1Item } = useTaskItem("DICT_1");
  const { item: dict2Item } = useTaskItem("DICT_2");
  const toggleDelete1 = useBoolToggle(useTaskItem, "TOGGLE_DELETE_1");
  const toggleDelete2 = useBoolToggle(useTaskItem, "TOGGLE_DELETE_2");
  const toggleChange1 = useBoolToggle(useTaskItem, "TOGGLE_CHANGE_1");
  const toggleChange2 = useBoolToggle(useTaskItem, "TOGGLE_CHANGE_2");
  const toggleCharge1 = useBoolToggle(useTaskItem, "TOGGLE_CHARGE_1");
  const toggleCharge2 = useBoolToggle(useTaskItem, "TOGGLE_CHARGE_2");
  const toggleLink = useBoolToggle(useTaskItem, "TOGGLE_LINK");
  const { value: LINK_MODE } = useTaskItem("LINK_MODE");

  // _LIST / _TYPE UI fields and their plain backend counterparts
  const { value: DELETE_1_LIST } = useTaskItem("DELETE_1_LIST");
  const { value: DELETE_2_LIST } = useTaskItem("DELETE_2_LIST");
  const { value: CHARGE_1_LIST } = useTaskItem("CHARGE_1_LIST");
  const { value: CHARGE_2_LIST } = useTaskItem("CHARGE_2_LIST");
  const { value: CHANGE_BOND_1_LIST } = useTaskItem("CHANGE_BOND_1_LIST");
  const { value: CHANGE_BOND_2_LIST } = useTaskItem("CHANGE_BOND_2_LIST");
  const { value: CHANGE_BOND_1_TYPE } = useTaskItem("CHANGE_BOND_1_TYPE");
  const { value: CHANGE_BOND_2_TYPE } = useTaskItem("CHANGE_BOND_2_TYPE");
  const { value: DELETE_1, item: delete1Item } = useTaskItem("DELETE_1");
  const { value: DELETE_2, item: delete2Item } = useTaskItem("DELETE_2");
  const { value: CHARGE_1, item: charge1Item } = useTaskItem("CHARGE_1");
  const { value: CHARGE_2, item: charge2Item } = useTaskItem("CHARGE_2");
  const { value: CHANGE_BOND_1, item: changeBond1Item } = useTaskItem("CHANGE_BOND_1");
  const { value: CHANGE_BOND_2, item: changeBond2Item } = useTaskItem("CHANGE_BOND_2");
  const { value: CHANGE_1_TYPE, item: change1TypeItem } = useTaskItem("CHANGE_1_TYPE");
  const { value: CHANGE_2_TYPE, item: change2TypeItem } = useTaskItem("CHANGE_2_TYPE");

  // --- TLC mode: fetch monomer info from CCP4 library ---
  const mon1ResName = MON_1_TYPE === "TLC" ? RES_NAME_1_TLC : undefined;
  const mon2ResName = MON_2_TYPE === "TLC" ? RES_NAME_2_TLC : undefined;
  const mon1TLC = useMonomerInfo(mon1ResName);
  const mon2TLC = useMonomerInfo(mon2ResName);

  // --- CIF mode: per-monomer dictionary from digest ---
  const [dict1Monomers, setDict1Monomers] = useState<MonomerDict>({});
  const [dict2Monomers, setDict2Monomers] = useState<MonomerDict>({});

  const handleDict1Change = useCallback(async () => {
    if (!dict1Item?._objectPath) return;
    const digest = await fetchDigest(dict1Item._objectPath);
    if (digest?.monomers && typeof digest.monomers === "object") {
      const codes = Object.keys(digest.monomers);
      setDict1Monomers(digest.monomers);
      // Auto-select first code only if current value isn't already valid
      if (codes.length > 0 && resNameItem1CIF?._objectPath && !codes.includes(RES_NAME_1_CIF ?? "")) {
        setParameter({ object_path: resNameItem1CIF._objectPath, value: codes[0] });
      }
    }
  }, [dict1Item?._objectPath, fetchDigest, resNameItem1CIF?._objectPath, RES_NAME_1_CIF, setParameter]);

  const handleDict2Change = useCallback(async () => {
    if (!dict2Item?._objectPath) return;
    const digest = await fetchDigest(dict2Item._objectPath);
    if (digest?.monomers && typeof digest.monomers === "object") {
      const codes = Object.keys(digest.monomers);
      setDict2Monomers(digest.monomers);
      // Auto-select first code only if current value isn't already valid
      if (codes.length > 0 && resNameItem2CIF?._objectPath && !codes.includes(RES_NAME_2_CIF ?? "")) {
        setParameter({ object_path: resNameItem2CIF._objectPath, value: codes[0] });
      }
    }
  }, [dict2Item?._objectPath, fetchDigest, resNameItem2CIF?._objectPath, RES_NAME_2_CIF, setParameter]);

  // Fetch CIF digest on mount if dictionary files are already set
  useEffect(() => {
    if (MON_1_TYPE === "CIF" && dict1Item?.dbFileId) handleDict1Change();
  }, [MON_1_TYPE, dict1Item?.dbFileId, handleDict1Change]);

  useEffect(() => {
    if (MON_2_TYPE === "CIF" && dict2Item?.dbFileId) handleDict2Change();
  }, [MON_2_TYPE, dict2Item?.dbFileId, handleDict2Change]);

  // Derive active CIF monomer info from selected residue name
  const dict1Codes = useMemo(() => Object.keys(dict1Monomers), [dict1Monomers]);
  const dict2Codes = useMemo(() => Object.keys(dict2Monomers), [dict2Monomers]);
  const mon1CIF = (RES_NAME_1_CIF && dict1Monomers[RES_NAME_1_CIF]) ?? EMPTY_MONOMER;
  const mon2CIF = (RES_NAME_2_CIF && dict2Monomers[RES_NAME_2_CIF]) ?? EMPTY_MONOMER;

  // --- Merge: pick TLC or CIF monomer info based on mode ---
  const mon1 = MON_1_TYPE === "CIF" ? mon1CIF : mon1TLC;
  const mon2 = MON_2_TYPE === "CIF" ? mon2CIF : mon2TLC;
  const mon1BondLabels = useBondLabels(mon1?.bonds);
  const mon2BondLabels = useBondLabels(mon2?.bonds);

  // Auto-select the first atom when monomer info arrives, but only if the
  // current value is empty or not in the atom list (avoids redundant writes on load).
  const atomItem1 = MON_1_TYPE === "CIF" ? atomItem1CIF : atomItem1TLC;
  const atomItem2 = MON_2_TYPE === "CIF" ? atomItem2CIF : atomItem2TLC;
  const currentAtom1 = MON_1_TYPE === "CIF" ? ATOM_NAME_1_CIF : ATOM_NAME_1_TLC;
  const currentAtom2 = MON_2_TYPE === "CIF" ? ATOM_NAME_2_CIF : ATOM_NAME_2_TLC;

  const mon1Atoms = mon1?.atoms ?? [];
  const mon2Atoms = mon2?.atoms ?? [];

  useEffect(() => {
    if (mon1Atoms.length > 0 && atomItem1?._objectPath && !mon1Atoms.includes(currentAtom1 ?? "")) {
      setParameter({ object_path: atomItem1._objectPath, value: mon1Atoms[0] });
    }
  }, [mon1Atoms, currentAtom1, atomItem1?._objectPath, setParameter]);

  useEffect(() => {
    if (mon2Atoms.length > 0 && atomItem2?._objectPath && !mon2Atoms.includes(currentAtom2 ?? "")) {
      setParameter({ object_path: atomItem2._objectPath, value: mon2Atoms[0] });
    }
  }, [mon2Atoms, currentAtom2, atomItem2?._objectPath, setParameter]);

  // Sync the active mode's atom value into the plain ATOM_NAME_1/ATOM_NAME_2
  // fields that the backend script reads to build AceDRG instructions.
  // Only write when the value actually differs to avoid redundant API calls.
  useEffect(() => {
    if (currentAtom1 && atomItem1Plain?._objectPath && currentAtom1 !== ATOM_NAME_1) {
      setParameter({ object_path: atomItem1Plain._objectPath, value: currentAtom1 });
    }
  }, [currentAtom1, ATOM_NAME_1, atomItem1Plain?._objectPath, setParameter]);

  useEffect(() => {
    if (currentAtom2 && atomItem2Plain?._objectPath && currentAtom2 !== ATOM_NAME_2) {
      setParameter({ object_path: atomItem2Plain._objectPath, value: currentAtom2 });
    }
  }, [currentAtom2, ATOM_NAME_2, atomItem2Plain?._objectPath, setParameter]);

  // Sync _LIST / _TYPE UI fields → plain backend fields.
  // The backend reads DELETE_1, CHARGE_1, CHANGE_BOND_1, CHANGE_1_TYPE etc.
  // but the GUI populates DELETE_1_LIST, CHARGE_1_LIST, CHANGE_BOND_1_LIST, CHANGE_BOND_1_TYPE.
  // Bond format: UI uses "C1-C2", backend expects "C1 -- C2".
  useEffect(() => {
    if (DELETE_1_LIST && delete1Item?._objectPath && DELETE_1_LIST !== DELETE_1)
      setParameter({ object_path: delete1Item._objectPath, value: DELETE_1_LIST });
  }, [DELETE_1_LIST, DELETE_1, delete1Item?._objectPath, setParameter]);

  useEffect(() => {
    if (DELETE_2_LIST && delete2Item?._objectPath && DELETE_2_LIST !== DELETE_2)
      setParameter({ object_path: delete2Item._objectPath, value: DELETE_2_LIST });
  }, [DELETE_2_LIST, DELETE_2, delete2Item?._objectPath, setParameter]);

  useEffect(() => {
    if (CHARGE_1_LIST && charge1Item?._objectPath && CHARGE_1_LIST !== CHARGE_1)
      setParameter({ object_path: charge1Item._objectPath, value: CHARGE_1_LIST });
  }, [CHARGE_1_LIST, CHARGE_1, charge1Item?._objectPath, setParameter]);

  useEffect(() => {
    if (CHARGE_2_LIST && charge2Item?._objectPath && CHARGE_2_LIST !== CHARGE_2)
      setParameter({ object_path: charge2Item._objectPath, value: CHARGE_2_LIST });
  }, [CHARGE_2_LIST, CHARGE_2, charge2Item?._objectPath, setParameter]);

  useEffect(() => {
    if (CHANGE_BOND_1_LIST && changeBond1Item?._objectPath) {
      const backendVal = CHANGE_BOND_1_LIST.replace("-", " -- ");
      if (backendVal !== CHANGE_BOND_1)
        setParameter({ object_path: changeBond1Item._objectPath, value: backendVal });
    }
  }, [CHANGE_BOND_1_LIST, CHANGE_BOND_1, changeBond1Item?._objectPath, setParameter]);

  useEffect(() => {
    if (CHANGE_BOND_2_LIST && changeBond2Item?._objectPath) {
      const backendVal = CHANGE_BOND_2_LIST.replace("-", " -- ");
      if (backendVal !== CHANGE_BOND_2)
        setParameter({ object_path: changeBond2Item._objectPath, value: backendVal });
    }
  }, [CHANGE_BOND_2_LIST, CHANGE_BOND_2, changeBond2Item?._objectPath, setParameter]);

  useEffect(() => {
    if (CHANGE_BOND_1_TYPE && change1TypeItem?._objectPath && CHANGE_BOND_1_TYPE !== CHANGE_1_TYPE)
      setParameter({ object_path: change1TypeItem._objectPath, value: CHANGE_BOND_1_TYPE });
  }, [CHANGE_BOND_1_TYPE, CHANGE_1_TYPE, change1TypeItem?._objectPath, setParameter]);

  useEffect(() => {
    if (CHANGE_BOND_2_TYPE && change2TypeItem?._objectPath && CHANGE_BOND_2_TYPE !== CHANGE_2_TYPE)
      setParameter({ object_path: change2TypeItem._objectPath, value: CHANGE_BOND_2_TYPE });
  }, [CHANGE_BOND_2_TYPE, CHANGE_2_TYPE, change2TypeItem?._objectPath, setParameter]);

  if (!container) return <LinearProgress />;

  // Build qualifier overrides for monomer-dependent fields.
  // Always pass enumerators (even empty) to suppress the "tmp" placeholder from def.xml.
  const mon1AtomQuals = { guiLabel: " ", enumerators: mon1Atoms };
  const mon2AtomQuals = { guiLabel: " ", enumerators: mon2Atoms };
  const mon1DeleteQuals = { guiLabel: "Atom to delete", enumerators: mon1Atoms };
  const mon2DeleteQuals = { guiLabel: "Atom to delete", enumerators: mon2Atoms };
  const mon1BondQuals = { guiLabel: " ", enumerators: mon1BondLabels };
  const mon2BondQuals = { guiLabel: " ", enumerators: mon2BondLabels };
  const mon1ChargeQuals = { guiLabel: " ", enumerators: mon1Atoms };
  const mon2ChargeQuals = { guiLabel: " ", enumerators: mon2Atoms };
  const resName1CIFQuals = { guiLabel: " ", enumerators: dict1Codes };
  const resName2CIFQuals = { guiLabel: " ", enumerators: dict2Codes };

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* --- First monomer to be linked --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "First monomer to be linked" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Get ligand description from" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                <CCP4i2TaskElement itemName="MON_1_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              {MON_1_TYPE === "CIF" && (
                <CCP4i2TaskElement itemName="DICT_1" {...props} onChange={handleDict1Change} />
              )}
              {MON_1_TYPE === "CIF" && (
                <InlineField label="Residue name" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="RES_NAME_1_CIF" {...props} qualifiers={resName1CIFQuals} />
                </InlineField>
              )}
              {MON_1_TYPE === "CIF" && (
                <InlineField label="Linking atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="ATOM_NAME_1_CIF" {...props} qualifiers={mon1AtomQuals} />
                </InlineField>
              )}
              {MON_1_TYPE === "TLC" && (
                <InlineField label="Residue name" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement
                    itemName="RES_NAME_1_TLC"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              )}
              {MON_1_TYPE === "TLC" && (
                <InlineField label="Linking atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="ATOM_NAME_1_TLC" {...props} qualifiers={mon1AtomQuals} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_DELETE_1" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
              {toggleDelete1.value && (
                <CCP4i2TaskElement itemName="DELETE_1_LIST" {...props} qualifiers={mon1DeleteQuals} />
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHANGE_1" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
              {toggleChange1.value && (
                <InlineField label="Bond" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHANGE_BOND_1_LIST" {...props} qualifiers={mon1BondQuals} />
                </InlineField>
              )}
              {toggleChange1.value && (
                <InlineField label="New bond order" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHANGE_BOND_1_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHARGE_1" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
              {toggleCharge1.value && (
                <InlineField label="Atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHARGE_1_LIST" {...props} qualifiers={mon1ChargeQuals} />
                </InlineField>
              )}
              {toggleCharge1.value && (
                <InlineField label="New charge" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHARGE_1_VALUE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Second monomer to be linked --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Second monomer to be linked" }}
            containerHint="FolderLevel"
          >
            <CCP4i2ContainerElement
              {...props}
              itemName=""
              qualifiers={{ initiallyOpen: true }}
              containerHint="BlockLevel"
            >
              <InlineField label="Get ligand description from" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                <CCP4i2TaskElement itemName="MON_2_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
              </InlineField>
              {MON_2_TYPE === "CIF" && (
                <CCP4i2TaskElement itemName="DICT_2" {...props} onChange={handleDict2Change} />
              )}
              {MON_2_TYPE === "CIF" && (
                <InlineField label="Residue name" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="RES_NAME_2_CIF" {...props} qualifiers={resName2CIFQuals} />
                </InlineField>
              )}
              {MON_2_TYPE === "CIF" && (
                <InlineField label="Linking atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="ATOM_NAME_2_CIF" {...props} qualifiers={mon2AtomQuals} />
                </InlineField>
              )}
              {MON_2_TYPE === "TLC" && (
                <InlineField label="Residue name" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement
                    itemName="RES_NAME_2_TLC"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              )}
              {MON_2_TYPE === "TLC" && (
                <InlineField label="Linking atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="ATOM_NAME_2_TLC" {...props} qualifiers={mon2AtomQuals} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_DELETE_2" {...props} qualifiers={{ guiLabel: "Delete atom" }} />
              {toggleDelete2.value && (
                <CCP4i2TaskElement itemName="DELETE_2_LIST" {...props} qualifiers={mon2DeleteQuals} />
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHANGE_2" {...props} qualifiers={{ guiLabel: "Change order of bond" }} />
              {toggleChange2.value && (
                <InlineField label="Bond" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHANGE_BOND_2_LIST" {...props} qualifiers={mon2BondQuals} />
                </InlineField>
              )}
              {toggleChange2.value && (
                <InlineField label="New bond order" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHANGE_BOND_2_TYPE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
              <CCP4i2TaskElement itemName="TOGGLE_CHARGE_2" {...props} qualifiers={{ guiLabel: "Change the formal charge of an atom" }} />
              {toggleCharge2.value && (
                <InlineField label="Atom" width={DROPDOWN_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHARGE_2_LIST" {...props} qualifiers={mon2ChargeQuals} />
                </InlineField>
              )}
              {toggleCharge2.value && (
                <InlineField label="New charge" width={SHORT_FIELD} labelWidth={LABEL_WIDTH}>
                  <CCP4i2TaskElement itemName="CHARGE_2_VALUE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
              )}
            </CCP4i2ContainerElement>
          </CCP4i2ContainerElement>

          {/* --- Bond order --- */}
          <InlineField label="Order of the bond between linked atoms" width={SHORT_FIELD}>
            <CCP4i2TaskElement itemName="BOND_ORDER" {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>

          {/* --- Apply links to model (optional) --- */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Apply links to model (optional)" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="TOGGLE_LINK" {...props} qualifiers={{ guiLabel: "Apply links to model" }} />
            {toggleLink.value && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="XYZIN" {...props} />
                <InlineField label="Apply links" width={DROPDOWN_FIELD}>
                  <CCP4i2TaskElement itemName="LINK_MODE" {...props} qualifiers={{ guiLabel: " " }} />
                </InlineField>
                {LINK_MODE === "AUTO" && (
                  <InlineField label="within" width={SHORT_FIELD} hint="times the dictionary value for this bond">
                    <CCP4i2TaskElement itemName="LINK_DISTANCE" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
                {LINK_MODE === "MANUAL" && (
                  <InlineField label="Create link between residues:" width={DROPDOWN_FIELD}>
                    <CCP4i2TaskElement itemName="MODEL_RES_LIST" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
                {LINK_MODE === "MANUAL" && (
                  <InlineField label="Filter list by atom proximity:" width={SHORT_FIELD}>
                    <CCP4i2TaskElement itemName="MODEL_LINK_DISTANCE" {...props} qualifiers={{ guiLabel: " " }} />
                  </InlineField>
                )}
              </CCP4i2ContainerElement>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab key="controlParameters" label="Advanced">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Advanced AceDRG options" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="EXTRA_ACEDRG_INSTRUCTIONS" {...props} />
            <CCP4i2TaskElement itemName="EXTRA_ACEDRG_KEYWORDS" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
