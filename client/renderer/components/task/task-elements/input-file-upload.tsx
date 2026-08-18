import { IconButton, styled, SxProps, Tooltip } from "@mui/material";
import { ChangeEvent } from "react";
import { Folder } from "@mui/icons-material";

interface InputFileUploadProps {
  handleFileChange: (ev: ChangeEvent<HTMLInputElement>) => void;
  disabled: boolean;
  accept: string;
  sx?: SxProps;
}

/**
 * "Choose a file" source button: browse the local file system.
 * One of the three file-source buttons grouped together in CDataFileElement
 * (file system / project database / internet) — keep its look and feel in
 * step with input-file-fetch.tsx and the project-database button.
 */
export const InputFileUpload: React.FC<InputFileUploadProps> = ({
  handleFileChange,
  disabled,
  accept,
  sx,
}) => {
  return (
    <Tooltip title="Browse files on this computer">
      <span>
        <IconButton
          disabled={disabled}
          component="label"
          size="small"
          tabIndex={-1}
          aria-label="Browse files on this computer"
          sx={sx}
        >
          <Folder fontSize="small" />
          <VisuallyHiddenInput
            disabled={disabled}
            type="file"
            onChange={handleFileChange}
            accept={accept}
          />
        </IconButton>
      </span>
    </Tooltip>
  );
};

export const VisuallyHiddenInput = styled("input")({
  clip: "rect(0 0 0 0)",
  clipPath: "inset(50%)",
  height: 1,
  overflow: "hidden",
  position: "absolute",
  bottom: 0,
  left: 0,
  whiteSpace: "nowrap",
  width: 1,
});
