import { ChangeEventHandler } from "react";
import { Button, styled } from "@mui/material";
import { Upload } from "@mui/icons-material";

const VisuallyHiddenInput = styled("input")({
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

export default function FileUpload({
  text,
  onChange,
}: {
  text: string;
  onChange: ChangeEventHandler<HTMLInputElement>;
}) {
  return (
    <Button component="label" variant="contained" startIcon={<Upload />}>
      {text}
      <VisuallyHiddenInput type="file" multiple onChange={onChange} />
    </Button>
  );
}
