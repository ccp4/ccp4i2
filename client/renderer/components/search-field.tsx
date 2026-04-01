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
import {
  IconButton,
  InputAdornment,
  SxProps,
  TextField,
  Theme,
} from "@mui/material";
import { Close, Search } from "@mui/icons-material";

export default function SearchField(props: {
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
  size?: "small" | "medium";
  sx?: SxProps<Theme>;
}) {
  return (
    <TextField
      value={props.value}
      onChange={(e) => props.onChange(e.target.value)}
      placeholder={props.placeholder || "Search..."}
      slotProps={{
        input: {
          startAdornment: (
            <InputAdornment position="start">
              <Search />
            </InputAdornment>
          ),
          endAdornment: (
            <InputAdornment position="end">
              <IconButton
                onClick={() => props.onChange("")}
                sx={{ visibility: props.value ? "visible" : "hidden" }}
              >
                <Close />
              </IconButton>
            </InputAdornment>
          ),
        },
      }}
      fullWidth
      size={props.size}
      sx={props.sx}
    />
  );
}
