"use client";

/**
 * "Add ligand here": put the dataset's fragment at the view centre in one
 * click, instead of the seven-interaction Get monomer / Merge / tidy-up
 * sequence (docs/campaign-place-ligand-design.md).
 *
 * The button's whole claim is that it already knows which ligand and which
 * molecule. When it does not, it is disabled and says why; it never asks for
 * a code, because a version that asks is the Get monomer dialog, which is
 * one menu away. With several candidate codes it offers them all, since two
 * fragments in one dictionary is exactly the case where guessing produces a
 * wrong model that looks right.
 */

import React, { useCallback, useMemo, useState } from "react";
import { Button, Menu, MenuItem, Tooltip } from "@mui/material";
import { Add as AddIcon, ArrowDropDown as ArrowDropDownIcon } from "@mui/icons-material";
import { moleculeForFile } from "../../lib/ligand-codes";

export const NO_DICTIONARY_REASON = "No ligand dictionary for this job";
export const NO_MOLECULE_REASON = "The dataset's coordinates are not loaded";

interface AddLigandButtonProps {
  /** The codes the dataset's dictionaries define; empty disables the button. */
  ligandCodes: string[];
  /** Everything loaded in Moorhen; the target is found by file, never by position. */
  molecules: { uniqueId?: string }[];
  /** The file the dataset's coordinates were loaded from. */
  memberCoordFileId: number | null;
  /** Place one code. Reporting the outcome is the caller's job. */
  onAddLigand: (code: string) => Promise<void>;
}

export const AddLigandButton: React.FC<AddLigandButtonProps> = ({
  ligandCodes,
  molecules,
  memberCoordFileId,
  onAddLigand,
}) => {
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const [busy, setBusy] = useState(false);

  const target = useMemo(
    () => moleculeForFile(molecules, memberCoordFileId),
    [molecules, memberCoordFileId],
  );
  const disabledReason =
    ligandCodes.length === 0 ? NO_DICTIONARY_REASON : !target ? NO_MOLECULE_REASON : null;

  const place = useCallback(
    async (code: string) => {
      setAnchorEl(null);
      setBusy(true);
      try {
        await onAddLigand(code);
      } finally {
        setBusy(false);
      }
    },
    [onAddLigand],
  );

  const several = ligandCodes.length > 1;
  const button = (
    <Button
      size="small"
      variant="outlined"
      startIcon={<AddIcon />}
      endIcon={several ? <ArrowDropDownIcon /> : undefined}
      disabled={!!disabledReason || busy}
      onClick={(e) => (several ? setAnchorEl(e.currentTarget) : place(ligandCodes[0]))}
    >
      Add ligand here
    </Button>
  );

  return (
    <>
      <Tooltip
        title={disabledReason ?? (several ? "Choose which ligand to add" : `Add ${ligandCodes[0]} at the view centre`)}
      >
        {/* A disabled button fires no events, so the tooltip needs this wrapper. */}
        <span>{button}</span>
      </Tooltip>
      {several && (
        <Menu anchorEl={anchorEl} open={anchorEl !== null} onClose={() => setAnchorEl(null)}>
          {ligandCodes.map((code) => (
            <MenuItem key={code} onClick={() => place(code)}>
              {code}
            </MenuItem>
          ))}
        </Menu>
      )}
    </>
  );
};
