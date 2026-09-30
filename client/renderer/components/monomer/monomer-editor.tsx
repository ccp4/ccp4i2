"use client";

/**
 * One monomer, drawn from its dictionary with its edits laid over it, and
 * edited by clicking: choose the linking atom, delete atoms, set bond orders
 * and formal charges.
 *
 * The edits are a description of the modified monomer (lib/monomer-edits.ts),
 * not commands: the server writes AceDRG's instructions from them. So there
 * is no order to get wrong and no way to say two things about one atom or
 * bond. Nothing here knows about jobs; the caller stores what it is given.
 */
import { useMemo, useState } from "react";
import {
  Alert,
  Box,
  Button,
  Chip,
  Stack,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";

import { MonomerPicker } from "./monomer-picker";
import type { MonomerAtomDetail, MonomerBond } from "../../lib/monomer-molblock";
import {
  BOND_ORDERS,
  dictionaryOrder,
  editedCharge,
  editedOrder,
  setBondOrder,
  setCharge,
  toggleDelete,
  valenceAdvisories,
  type BondOrder,
  type MonomerEdits,
} from "../../lib/monomer-edits";

/** A monomer as the server describes it: atom names, bonds, atom details. */
export interface EditableMonomer {
  atoms: string[];
  bonds: MonomerBond[];
  atom_details?: MonomerAtomDetail[];
}

/** What a click on the drawing does. */
type PickMode = "link" | "delete" | "bond" | "charge";

const PICK_HINTS: Record<PickMode, string> = {
  link: "Click the atom that bonds to the other monomer",
  delete: "Click an atom to delete it; click again to restore it",
  bond: "Click a bond to choose its order",
  charge: "Click an atom to change its formal charge",
};

const chargeText = (charge: number) => (charge > 0 ? `+${charge}` : `${charge}`);

export interface MonomerEditorProps {
  monomer: EditableMonomer;
  code?: string;
  emptyMessage: string;
  linkAtom?: string;
  linkOrder: BondOrder;
  edits: MonomerEdits;
  onPickLink: (atom: string) => void;
  onEdits: (edits: MonomerEdits) => void;
}

export function MonomerEditor({
  monomer,
  code,
  emptyMessage,
  linkAtom,
  linkOrder,
  edits,
  onPickLink,
  onEdits,
}: MonomerEditorProps) {
  const [mode, setMode] = useState<PickMode>("link");
  const [bond, setBond] = useState<{ atom1: string; atom2: string } | null>(null);
  const [atom, setAtom] = useState<string | null>(null);
  const [notice, setNotice] = useState<string | null>(null);

  const atomDetails = monomer.atom_details ?? [];
  const known = new Set(monomer.atoms);
  const deleted = new Set(edits.deletes);

  // A selection from another monomer, or one since deleted, is no selection.
  const selectedBond =
    bond && known.has(bond.atom1) && known.has(bond.atom2) ? bond : null;
  const selectedAtom = atom && known.has(atom) ? atom : null;

  const advisories = useMemo(
    () => valenceAdvisories(atomDetails, monomer.bonds, edits, linkAtom, linkOrder),
    [atomDetails, monomer.bonds, edits, linkAtom, linkOrder]
  );

  const pickAtom = (name: string) => {
    setNotice(null);
    if (mode === "link") {
      if (deleted.has(name)) {
        setNotice(`${name} is deleted, so it cannot be the linking atom`);
        return;
      }
      onPickLink(name);
    } else if (mode === "delete") {
      if (name === linkAtom) {
        setNotice(`${name} is the linking atom, so it cannot be deleted`);
        return;
      }
      onEdits(toggleDelete(edits, name));
    } else if (mode === "charge") {
      setAtom(name);
    }
  };

  const dictionaryBond = selectedBond
    ? monomer.bonds.find(
        (b) =>
          (b.atom1 === selectedBond.atom1 && b.atom2 === selectedBond.atom2) ||
          (b.atom1 === selectedBond.atom2 && b.atom2 === selectedBond.atom1)
      )
    : undefined;
  const bondBaseline = dictionaryOrder(dictionaryBond?.type);
  const bondGone = selectedBond ? deleted.has(selectedBond.atom1) || deleted.has(selectedBond.atom2) : false;

  const atomDetail = selectedAtom ? atomDetails.find((a) => a.name === selectedAtom) : undefined;
  const atomBaseline = atomDetail?.charge ?? 0;
  const atomCharge = selectedAtom ? editedCharge(edits, selectedAtom) ?? atomBaseline : 0;

  const hasEdits = edits.deletes.length + edits.bondOrders.length + edits.charges.length > 0;

  return (
    <Box sx={{ my: 1, maxWidth: 420 }}>
      <ToggleButtonGroup
        size="small"
        exclusive
        value={mode}
        onChange={(_, value: PickMode | null) => {
          if (value) {
            setMode(value);
            setNotice(null);
          }
        }}
        aria-label="What a click on the drawing does"
        sx={{ mb: 0.5, flexWrap: "wrap" }}
      >
        <ToggleButton value="link">Linking atom</ToggleButton>
        <ToggleButton value="delete">Delete</ToggleButton>
        <ToggleButton value="bond">Bond order</ToggleButton>
        <ToggleButton value="charge">Charge</ToggleButton>
      </ToggleButtonGroup>
      <Typography variant="caption" color="text.secondary" component="div" sx={{ mb: 0.5 }}>
        {PICK_HINTS[mode]}
      </Typography>
      <MonomerPicker
        atomDetails={atomDetails}
        bonds={monomer.bonds}
        mode={mode === "bond" ? "bond" : "atom"}
        selectedAtom={mode === "link" ? linkAtom ?? null : mode === "charge" ? selectedAtom : null}
        selectedBond={mode === "bond" ? selectedBond : null}
        markedAtoms={advisories.map((a) => a.atom)}
        edits={edits}
        onPickAtom={pickAtom}
        onPickBond={(atom1, atom2) => {
          setNotice(null);
          setBond({ atom1, atom2 });
        }}
        title={code}
        emptyMessage={emptyMessage}
      />

      <Typography variant="body2" sx={{ mt: 0.5 }} data-testid="link-atom-line">
        Linking atom:{" "}
        {linkAtom && known.has(linkAtom) ? (
          <b>{linkAtom}</b>
        ) : linkAtom && monomer.atoms.length > 0 ? (
          <Typography component="span" variant="body2" color="warning.main">
            {linkAtom} is not an atom of {code ?? "this monomer"}; click one
          </Typography>
        ) : (
          <Typography component="span" variant="body2" color="text.secondary">
            not chosen yet
          </Typography>
        )}
      </Typography>

      {notice && (
        <Typography variant="body2" color="warning.main" sx={{ mt: 0.5 }}>
          {notice}
        </Typography>
      )}

      {mode === "bond" && selectedBond && (
        <Stack direction="row" spacing={1} alignItems="center" sx={{ mt: 1 }} flexWrap="wrap">
          <Typography variant="body2">
            Bond {selectedBond.atom1}–{selectedBond.atom2}:
          </Typography>
          {bondGone ? (
            <Typography variant="body2" color="text.secondary">
              one of its atoms is deleted
            </Typography>
          ) : bondBaseline === null ? (
            <Typography variant="body2" color="text.secondary">
              {dictionaryBond?.type ?? "this"} bonds cannot be changed here
            </Typography>
          ) : (
            <>
              <ToggleButtonGroup
                size="small"
                exclusive
                value={editedOrder(edits, selectedBond.atom1, selectedBond.atom2) ?? bondBaseline}
                onChange={(_, value: BondOrder | null) =>
                  value && onEdits(setBondOrder(edits, selectedBond.atom1, selectedBond.atom2, value, bondBaseline))
                }
                aria-label={`Order of bond ${selectedBond.atom1}-${selectedBond.atom2}`}
              >
                {BOND_ORDERS.map((order) => (
                  <ToggleButton key={order} value={order}>
                    {order.toLowerCase()}
                  </ToggleButton>
                ))}
              </ToggleButtonGroup>
              <Typography variant="caption" color="text.secondary">
                dictionary: {bondBaseline.toLowerCase()}
              </Typography>
            </>
          )}
        </Stack>
      )}

      {mode === "charge" && selectedAtom && (
        <Stack direction="row" spacing={1} alignItems="center" sx={{ mt: 1 }}>
          <Typography variant="body2">Charge on {selectedAtom}:</Typography>
          {deleted.has(selectedAtom) ? (
            <Typography variant="body2" color="text.secondary">
              it is deleted
            </Typography>
          ) : (
            <>
              <Button
                size="small"
                variant="outlined"
                sx={{ minWidth: 32 }}
                disabled={atomCharge <= -3}
                aria-label="Lower the charge"
                onClick={() => onEdits(setCharge(edits, selectedAtom, atomCharge - 1, atomBaseline))}
              >
                −
              </Button>
              <Typography variant="body2" sx={{ minWidth: 24, textAlign: "center" }} data-testid="charge-value">
                {chargeText(atomCharge)}
              </Typography>
              <Button
                size="small"
                variant="outlined"
                sx={{ minWidth: 32 }}
                disabled={atomCharge >= 3}
                aria-label="Raise the charge"
                onClick={() => onEdits(setCharge(edits, selectedAtom, atomCharge + 1, atomBaseline))}
              >
                +
              </Button>
              <Typography variant="caption" color="text.secondary">
                dictionary: {chargeText(atomBaseline)}
              </Typography>
            </>
          )}
        </Stack>
      )}

      {hasEdits && (
        <Stack direction="row" spacing={0.5} useFlexGap flexWrap="wrap" sx={{ mt: 1 }} aria-label="Changes to this monomer">
          {edits.deletes.map((name) => (
            <Chip key={`d-${name}`} size="small" label={`delete ${name}`} onDelete={() => onEdits(toggleDelete(edits, name))} />
          ))}
          {edits.bondOrders.map((b) => (
            <Chip
              key={`b-${b.atom1}-${b.atom2}`}
              size="small"
              label={`${b.atom1}–${b.atom2} ${b.order.toLowerCase()}`}
              onDelete={() => onEdits(setBondOrder(edits, b.atom1, b.atom2, null, null))}
            />
          ))}
          {edits.charges.map((c) => (
            <Chip
              key={`c-${c.atom}`}
              size="small"
              label={`${c.atom} ${chargeText(c.charge)}`}
              onDelete={() => onEdits(setCharge(edits, c.atom, null, 0))}
            />
          ))}
        </Stack>
      )}

      {advisories.length > 0 && (
        <Alert severity="warning" sx={{ mt: 1, py: 0 }}>
          {advisories.map((a) => (
            <div key={a.atom}>{a.message}</div>
          ))}
          <Typography variant="caption" component="div">
            AceDRG has the final say; it will refuse a monomer it cannot make.
          </Typography>
        </Alert>
      )}
    </Box>
  );
}

export default MonomerEditor;
