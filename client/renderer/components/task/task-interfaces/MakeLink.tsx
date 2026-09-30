import { useCallback, useEffect, useMemo, useState } from "react";
import { Box, LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";
import { apiJson } from "../../../api-fetch";
import { MonomerEditor } from "../../monomer/monomer-editor";
import type { MonomerAtomDetail } from "../../../lib/monomer-molblock";
import {
  NO_EDITS,
  dictionaryOrder,
  editedCharge,
  editedOrder,
  type BondOrder,
  type MonomerEdits,
} from "../../../lib/monomer-edits";

// Layout constants
// Inside a monomer section, which may have half the pane.
const MONOMER_LABEL_WIDTH = "7.5rem";
const MONOMER_FIELD = "14rem";
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
  /** The same atoms in the same order, with element, charge and hydrogens. */
  atom_details?: MonomerAtomDetail[];
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
    apiJson<{
      success: boolean;
      data?: { atoms: string[]; bonds: MonomerBond[]; atom_details?: MonomerAtomDetail[] };
    }>(`monomer-info/${code}`)
      .then((res) => {
        if (!cancelled && res.success && res.data) {
          setInfo({
            atoms: res.data.atoms,
            bonds: res.data.bonds,
            atom_details: res.data.atom_details ?? [],
          });
        } else if (!cancelled) {
          // Not the previous monomer's atoms under a new name.
          setInfo(EMPTY_MONOMER);
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

const isTruthy = (value: unknown) => value === true || value === "True" || value === "true";

/**
 * A monomer's edits as the task holds them: the DELETE_ATOMS_n, BOND_ORDERS_n
 * and CHARGES_n lists, with the single-edit fields of older jobs folded in
 * when ticked -- the same fold the plugin makes (MakeLink.monomerEdits), so
 * what is drawn is what will run.
 */
function readEdits(values: {
  deletes: unknown;
  bondOrders: unknown;
  charges: unknown;
  toggleDelete: unknown;
  deleteAtom: unknown;
  toggleChange: unknown;
  changeBond: unknown;
  changeType: unknown;
  toggleCharge: unknown;
  chargeAtom: unknown;
  chargeValue: unknown;
}): MonomerEdits {
  let edits: MonomerEdits = {
    deletes: Array.isArray(values.deletes) ? values.deletes.filter(Boolean).map(String) : [],
    bondOrders: Array.isArray(values.bondOrders)
      ? values.bondOrders
          .filter((b: any) => b?.ATOM_1 && b?.ATOM_2 && b?.ORDER)
          .map((b: any) => ({ atom1: b.ATOM_1, atom2: b.ATOM_2, order: String(b.ORDER).toUpperCase() as BondOrder }))
      : [],
    charges: Array.isArray(values.charges)
      ? values.charges
          .filter((c: any) => c?.ATOM)
          .map((c: any) => ({ atom: c.ATOM, charge: Number(c.CHARGE ?? 0) }))
      : [],
  };
  if (isTruthy(values.toggleDelete) && values.deleteAtom && !edits.deletes.includes(String(values.deleteAtom))) {
    edits = { ...edits, deletes: [...edits.deletes, String(values.deleteAtom)] };
  }
  if (isTruthy(values.toggleChange) && values.changeBond && values.changeType) {
    const [atom1, atom2] = String(values.changeBond).split(" -- ").map((s) => s.trim());
    if (atom1 && atom2 && !editedOrder(edits, atom1, atom2)) {
      edits = { ...edits, bondOrders: [...edits.bondOrders, { atom1, atom2, order: String(values.changeType).toUpperCase() as BondOrder }] };
    }
  }
  if (isTruthy(values.toggleCharge) && values.chargeAtom && editedCharge(edits, String(values.chargeAtom)) === null) {
    edits = { ...edits, charges: [...edits.charges, { atom: String(values.chargeAtom), charge: Number(values.chargeValue ?? 0) }] };
  }
  return edits;
}

/**
 * Everything one monomer section needs: where its description comes from,
 * its linking atom, and its edits, with the writes that change them.
 */
function useMonomer(
  n: 1 | 2,
  useTaskItem: ReturnType<typeof useJob>["useTaskItem"],
  setParameter: ReturnType<typeof useJob>["setParameter"],
  fetchDigest: ReturnType<typeof useJob>["fetchDigest"]
) {
  const { value: type } = useTaskItem(`MON_${n}_TYPE`);
  const { value: tlc } = useTaskItem(`RES_NAME_${n}_TLC`);
  const { value: cifName, item: cifNameItem } = useTaskItem(`RES_NAME_${n}_CIF`);
  const { value: atomTLC, item: atomTLCItem } = useTaskItem(`ATOM_NAME_${n}_TLC`);
  const { value: atomCIF, item: atomCIFItem } = useTaskItem(`ATOM_NAME_${n}_CIF`);
  const { value: atomPlain, item: atomPlainItem } = useTaskItem(`ATOM_NAME_${n}`);
  const { item: dictItem } = useTaskItem(`DICT_${n}`);

  const { value: deletes, item: deletesItem } = useTaskItem(`DELETE_ATOMS_${n}`);
  const { value: bondOrders, item: bondOrdersItem } = useTaskItem(`BOND_ORDERS_${n}`);
  const { value: charges, item: chargesItem } = useTaskItem(`CHARGES_${n}`);
  const { value: toggleDeleteValue, item: toggleDeleteItem } = useTaskItem(`TOGGLE_DELETE_${n}`);
  const { value: deleteAtom } = useTaskItem(`DELETE_${n}`);
  const { value: toggleChangeValue, item: toggleChangeItem } = useTaskItem(`TOGGLE_CHANGE_${n}`);
  const { value: changeBond } = useTaskItem(`CHANGE_BOND_${n}`);
  const { value: changeType } = useTaskItem(`CHANGE_${n}_TYPE`);
  const { value: toggleChargeValue, item: toggleChargeItem } = useTaskItem(`TOGGLE_CHARGE_${n}`);
  const { value: chargeAtom } = useTaskItem(`CHARGE_${n}`);
  const { value: chargeValue } = useTaskItem(`CHARGE_${n}_VALUE`);

  const isCIF = type === "CIF";
  const libraryMonomer = useMonomerInfo(isCIF ? undefined : tlc);

  // CIF mode: the monomers in the chosen dictionary, from its digest.
  const [dictMonomers, setDictMonomers] = useState<MonomerDict>({});
  const handleDictChange = useCallback(async () => {
    if (!dictItem?._objectPath) return;
    const digest = await fetchDigest(dictItem._objectPath);
    if (digest?.monomers && typeof digest.monomers === "object") {
      const codes = Object.keys(digest.monomers);
      setDictMonomers(digest.monomers);
      // Auto-select only when the dictionary holds exactly one monomer.
      if (codes.length === 1 && cifNameItem?._objectPath && cifName !== codes[0]) {
        setParameter({ object_path: cifNameItem._objectPath, value: codes[0] });
      }
    }
  }, [dictItem?._objectPath, fetchDigest, cifNameItem?._objectPath, cifName, setParameter]);

  useEffect(() => {
    if (isCIF && dictItem?.dbFileId) handleDictChange();
  }, [isCIF, dictItem?.dbFileId, handleDictChange]);

  const dictCodes = useMemo(() => Object.keys(dictMonomers), [dictMonomers]);
  const monomer: MonomerInfo = isCIF ? (cifName && dictMonomers[cifName]) || EMPTY_MONOMER : libraryMonomer;

  const linkAtomItem = isCIF ? atomCIFItem : atomTLCItem;
  const linkAtom: string | undefined = (isCIF ? atomCIF : atomTLC) || undefined;

  // The plugin reads ATOM_NAME_n; keep it in step with the mode's own field.
  useEffect(() => {
    if (linkAtom && atomPlainItem?._objectPath && linkAtom !== atomPlain) {
      setParameter({ object_path: atomPlainItem._objectPath, value: linkAtom });
    }
  }, [linkAtom, atomPlain, atomPlainItem?._objectPath, setParameter]);

  const edits = useMemo(
    () =>
      readEdits({
        deletes, bondOrders, charges,
        toggleDelete: toggleDeleteValue, deleteAtom,
        toggleChange: toggleChangeValue, changeBond, changeType,
        toggleCharge: toggleChargeValue, chargeAtom, chargeValue,
      }),
    [deletes, bondOrders, charges, toggleDeleteValue, deleteAtom, toggleChangeValue,
     changeBond, changeType, toggleChargeValue, chargeAtom, chargeValue]
  );

  const pickLink = useCallback(
    (atom: string) => {
      if (linkAtomItem?._objectPath) setParameter({ object_path: linkAtomItem._objectPath, value: atom });
    },
    [linkAtomItem?._objectPath, setParameter]
  );

  const writeEdits = useCallback(
    (next: MonomerEdits) => {
      // Only the lists that changed; usually a click changes one.
      const legacy = [toggleDeleteValue, toggleChangeValue, toggleChargeValue].some(isTruthy);
      const changed = (a: unknown, b: unknown) => legacy || JSON.stringify(a) !== JSON.stringify(b);
      if (changed(next.deletes, edits.deletes) && deletesItem?._objectPath) {
        setParameter({ object_path: deletesItem._objectPath, value: next.deletes });
      }
      if (changed(next.bondOrders, edits.bondOrders) && bondOrdersItem?._objectPath) {
        setParameter({
          object_path: bondOrdersItem._objectPath,
          value: next.bondOrders.map((b) => ({ ATOM_1: b.atom1, ATOM_2: b.atom2, ORDER: b.order })),
        });
      }
      if (changed(next.charges, edits.charges) && chargesItem?._objectPath) {
        setParameter({
          object_path: chargesItem._objectPath,
          value: next.charges.map((c) => ({ ATOM: c.atom, CHARGE: c.charge })),
        });
      }
      // An older job's single edits are now in the lists: untick them, so
      // the lists are the only place the edits live.
      for (const [on, item] of [
        [toggleDeleteValue, toggleDeleteItem],
        [toggleChangeValue, toggleChangeItem],
        [toggleChargeValue, toggleChargeItem],
      ] as const) {
        if (isTruthy(on) && item?._objectPath) setParameter({ object_path: item._objectPath, value: false });
      }
    },
    [edits, deletesItem?._objectPath, bondOrdersItem?._objectPath, chargesItem?._objectPath, setParameter,
     toggleDeleteValue, toggleDeleteItem, toggleChangeValue, toggleChangeItem, toggleChargeValue, toggleChargeItem]
  );

  const code: string | undefined = (isCIF ? cifName : tlc?.trim().toUpperCase()) || undefined;
  const emptyMessage = isCIF
    ? dictCodes.length === 0
      ? "Choose a dictionary file to draw its monomer"
      : "Choose a residue name"
    : tlc?.trim()
      ? `"${tlc.trim().toUpperCase()}" is not in the CCP4 monomer library`
      : "Type a residue name to draw it";

  return {
    isCIF, monomer, code, emptyMessage, dictCodes, handleDictChange,
    linkAtom, pickLink, edits, writeEdits,
  };
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container, setParameter, fetchDigest } = useJob(props.job.id);
  const toggleLink = useBoolToggle(useTaskItem, "TOGGLE_LINK");
  // A cloned or autofilled job can carry a model with the toggle off, which
  // validity() warns about; the model stays in view so the warning can be acted on.
  const { item: xyzinItem } = useTaskItem("XYZIN");
  const showModel = Boolean(toggleLink.value || xyzinItem?.dbFileId);
  const { value: LINK_MODE } = useTaskItem("LINK_MODE");
  const { value: BOND_ORDER } = useTaskItem("BOND_ORDER");
  const linkOrder = (dictionaryOrder(BOND_ORDER) ?? "SINGLE") as BondOrder;

  const mon1 = useMonomer(1, useTaskItem, setParameter, fetchDigest);
  const mon2 = useMonomer(2, useTaskItem, setParameter, fetchDigest);

  if (!container) return <LinearProgress />;

  const monomerSection = (n: 1 | 2, mon: typeof mon1, label: string) => (
    <CCP4i2ContainerElement
      {...props}
      itemName=""
      qualifiers={{ guiLabel: label }}
      containerHint="FolderLevel"
    >
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ initiallyOpen: true }}
        containerHint="BlockLevel"
      >
        <InlineField label="Source" width={MONOMER_FIELD} labelWidth={MONOMER_LABEL_WIDTH}>
          <CCP4i2TaskElement itemName={`MON_${n}_TYPE`} {...props} qualifiers={{ guiLabel: " " }} />
        </InlineField>
        {mon.isCIF && (
          <CCP4i2TaskElement itemName={`DICT_${n}`} {...props} onChange={mon.handleDictChange} />
        )}
        {mon.isCIF ? (
          <InlineField label="Residue name" width={MONOMER_FIELD} labelWidth={MONOMER_LABEL_WIDTH}>
            <CCP4i2TaskElement
              itemName={`RES_NAME_${n}_CIF`}
              {...props}
              qualifiers={{ guiLabel: " ", enumerators: mon.dictCodes }}
            />
          </InlineField>
        ) : (
          <InlineField label="Residue name" width={SHORT_FIELD} labelWidth={MONOMER_LABEL_WIDTH}>
            <CCP4i2TaskElement itemName={`RES_NAME_${n}_TLC`} {...props} qualifiers={{ guiLabel: " " }} />
          </InlineField>
        )}
        <MonomerEditor
          monomer={mon.monomer}
          code={mon.code}
          emptyMessage={mon.emptyMessage}
          linkAtom={mon.linkAtom}
          linkOrder={linkOrder}
          edits={mon.edits ?? NO_EDITS}
          onPickLink={mon.pickLink}
          onEdits={mon.writeEdits}
        />
      </CCP4i2ContainerElement>
    </CCP4i2ContainerElement>
  );

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        <CCP4i2Tab key="inputData" label="Input data">
          {/* The two halves of one link, side by side when the pane has room
              for both; the grid follows the pane's width, not the window's. */}
          <Box
            sx={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fit, minmax(380px, 1fr))",
              gap: 1,
              alignItems: "start",
            }}
          >
            {monomerSection(1, mon1, "First monomer to be linked")}
            {monomerSection(2, mon2, "Second monomer to be linked")}
          </Box>

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
            {showModel && (
              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ initiallyOpen: true }}
                containerHint="BlockLevel"
              >
                <CCP4i2TaskElement itemName="XYZIN" {...props} />
                {toggleLink.value && (
                  <>
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
                  </>
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
