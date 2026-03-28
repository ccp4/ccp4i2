/*
 * Copyright (C) 2025 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
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
