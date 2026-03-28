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
import { FormControlLabel, Switch } from "@mui/material";
import { useCCP4i2Window } from "../app-context";

export const DevModeToggle: React.FC = () => {
  const { devMode, setDevMode } = useCCP4i2Window();

  const onToggleDevMode = async (
    ev: React.ChangeEvent<HTMLInputElement>
  ): Promise<void> => {
    setDevMode(!devMode);
    if (!window.electronAPI) {
      console.error("Electron API is not available");
      return;
    }
    window.electronAPI.sendMessage("toggle-dev-mode", {});
    ev.preventDefault();
    ev.stopPropagation();
  };

  return (
    <FormControlLabel
      control={
        <Switch
          checked={devMode}
          onChange={onToggleDevMode}
          name="devModeToggle"
          color="warning"
        />
      }
      label="Dev Mode"
    />
  );
};
