import { Button, styled, SxProps } from "@mui/material";
import { ChangeEvent, useContext } from "react";
import { TaskInterfaceContext } from "../../../providers/task-provider";
import { Folder } from "@mui/icons-material";

interface InputFileUploadProps {
  handleFileChange: (ev: ChangeEvent<HTMLInputElement>) => void;
  disabled: boolean;
  accept: string;
  sx?: SxProps;
}
export const InputFileUpload: React.FC<InputFileUploadProps> = ({
  handleFileChange,
  disabled,
  accept,
  sx,
}) => {
  return (
    <Button
      disabled={disabled}
      component="label"
      role={undefined}
      variant="outlined"
      tabIndex={-1}
      size="small"
      startIcon={<Folder fontSize="small" />}
      sx={sx}
    >
      <VisuallyHiddenInput
        disabled={disabled}
        type="file"
        onChange={handleFileChange}
        accept={accept}
        sx={sx}
      />
    </Button>
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
