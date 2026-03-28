/*
 * Copyright (C) 2026 Newcastle University
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
"use client";
import { ImportProjectContent } from "../../../components/import-project-content";
import { Paper } from "@mui/material";
import { NavigationShortcutsProvider } from "../../../providers/navigation-shortcuts-provider";
import CCP4i2TopBar from "../../../components/ccp4i2-topbar";

export default function ImportProjectPage() {
  return (
    <NavigationShortcutsProvider>
      <CCP4i2TopBar title="Import Project" showBackButton backPath="/ccp4i2" />
      <Paper
        sx={{
          display: "flex",
          justifyContent: "center",
          alignItems: "center",
          height: "calc(100vh - 80px)",
        }}
      >
        <ImportProjectContent />
      </Paper>
    </NavigationShortcutsProvider>
  );
}
