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
import { ChangeEvent, useEffect, useState } from "react";
import { IconButton, InputAdornment, TextField } from "@mui/material";
import { Close } from "@mui/icons-material";

export default function SearchField(props: {
  what: string;
  onDelay: (query: string) => void;
}) {
  const [query, setQuery] = useState("");

  const onChange = (event: ChangeEvent<HTMLInputElement>) => {
    setQuery(event.target.value);
  };

  useEffect(() => {
    const timeout = setTimeout(() => props.onDelay(query), 500);
    return () => clearTimeout(timeout);
  }, [query]);

  return (
    <TextField
      value={query}
      label={`Search ${props.what}`}
      onChange={onChange}
      InputProps={{
        endAdornment:
          query.length > 0 ? (
            <InputAdornment position="end">
              <IconButton onClick={() => setQuery("")}>
                <Close />
              </IconButton>
            </InputAdornment>
          ) : undefined,
      }}
      sx={{ flex: "auto" }}
    />
  );
}
