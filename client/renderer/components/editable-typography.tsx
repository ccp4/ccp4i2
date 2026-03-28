/*
 * Copyright (C) 2025-2026 Newcastle University
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
import { Typography } from "@mui/material";
import { KeyboardEvent, useRef } from "react";
import type { TypographyProps, SxProps, Theme } from "@mui/material";

type TypographyVariant = TypographyProps["variant"];

export default function EditableTypography(props: {
  variant: TypographyVariant;
  text: string;
  onDelay?: (text: string) => void;
  sx?: SxProps<Theme>;
}) {
  const timeout = useRef<NodeJS.Timeout | undefined>(undefined);

  function handleKeyUp(e: KeyboardEvent<HTMLSpanElement>) {
    const text = e.currentTarget.innerText.trim();
    clearTimeout(timeout.current);
    timeout.current = setTimeout(() => props.onDelay?.(text), 1000);
  }

  return (
    <Typography
      variant={props.variant}
      contentEditable={Boolean(props.onDelay)}
      suppressContentEditableWarning
      onKeyUp={handleKeyUp}
      sx={props.sx}
    >
      {props.text}
    </Typography>
  );
}
