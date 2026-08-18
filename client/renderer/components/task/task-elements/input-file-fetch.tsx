import { IconButton, SxProps, Tooltip } from "@mui/material";
import { useCallback } from "react";
import { useTaskInterface } from "../../../providers/task-provider";
import { Language } from "@mui/icons-material";

interface InputFileFetchProps {
  disabled: boolean;
  modes?: string[];
  item?: any;
  sx?: SxProps;
  onChange?: (updatedItem: any) => void;
}

/**
 * "Choose a file" source button: fetch from an online resource (PDB, EBI, ...).
 * Only offered for file types whose def.xml declares downloadModes, so this is
 * the one source button that is not always present.
 */
export const InputFileFetch: React.FC<InputFileFetchProps> = ({
  disabled,
  sx,
  item,
  modes,
  onChange,
}) => {
  const { setDownloadDialogOpen, setFetchItemParams } = useTaskInterface();

  const handleFetchClick = useCallback(
    (ev: any) => {
      ev.stopPropagation();
      if (setFetchItemParams) setFetchItemParams({ item, modes, onChange });
      if (setDownloadDialogOpen) setDownloadDialogOpen(true);
    },
    [item, modes, onChange, setDownloadDialogOpen, setFetchItemParams]
  );

  return (
    <Tooltip title="Fetch a file from an online resource">
      <span>
        <IconButton
          disabled={disabled}
          size="small"
          tabIndex={-1}
          aria-label="Fetch a file from an online resource"
          sx={sx}
          onClick={handleFetchClick}
        >
          <Language fontSize="small" />
        </IconButton>
      </span>
    </Tooltip>
  );
};
