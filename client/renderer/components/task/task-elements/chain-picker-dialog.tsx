/**
 * Dialog for selecting a chain from a multi-chain PDB entry.
 * Used when fetching coordinates for sequence extraction and the
 * structure contains multiple polymer chains.
 */

import {
  Button,
  Chip,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  List,
  ListItemButton,
  ListItemIcon,
  ListItemText,
  Radio,
  Typography,
} from "@mui/material";
import { useState } from "react";
import type { ChainSequenceInfo } from "./mmcif-sequence-parser";

interface ChainPickerDialogProps {
  open: boolean;
  onClose: () => void;
  onSelect: (chain: ChainSequenceInfo) => void;
  chains: ChainSequenceInfo[];
  pdbId: string;
}

const polymerColor: Record<string, "primary" | "info" | "success" | "warning"> = {
  PROTEIN: "primary",
  DNA: "info",
  RNA: "success",
  OTHER: "warning",
};

export const ChainPickerDialog: React.FC<ChainPickerDialogProps> = ({
  open,
  onClose,
  onSelect,
  chains,
  pdbId,
}) => {
  const [selected, setSelected] = useState<number>(0);

  return (
    <Dialog open={open} onClose={onClose} maxWidth="sm" fullWidth>
      <DialogTitle>
        Select chain from {pdbId.toUpperCase()}
      </DialogTitle>
      <DialogContent>
        <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
          This entry contains {chains.length} polymer chain{chains.length > 1 ? "s" : ""}.
          Select the chain to use as sequence:
        </Typography>
        <List dense>
          {chains.map((chain, idx) => (
            <ListItemButton
              key={`${chain.chainId}-${idx}`}
              selected={selected === idx}
              onClick={() => setSelected(idx)}
            >
              <ListItemIcon sx={{ minWidth: 36 }}>
                <Radio
                  checked={selected === idx}
                  size="small"
                />
              </ListItemIcon>
              <ListItemText
                primary={
                  <>
                    Chain {chain.chainId}
                    {" "}
                    <Chip
                      label={chain.polymerType}
                      size="small"
                      color={polymerColor[chain.polymerType] || "default"}
                      variant="outlined"
                      sx={{ ml: 1 }}
                    />
                  </>
                }
                secondary={`${chain.length} residues`}
              />
            </ListItemButton>
          ))}
        </List>
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose}>Cancel</Button>
        <Button
          variant="contained"
          onClick={() => {
            onSelect(chains[selected]);
            onClose();
          }}
        >
          Use chain {chains[selected]?.chainId}
        </Button>
      </DialogActions>
    </Dialog>
  );
};
