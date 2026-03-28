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
"use client";
import React from "react";
import { useRouter } from "next/navigation";
import { useMsal } from "@azure/msal-react";
import Toolbar from "@mui/material/Toolbar";
import IconButton from "@mui/material/IconButton";
import ArrowBackIcon from "@mui/icons-material/ArrowBack";
import ArrowForwardIcon from "@mui/icons-material/ArrowForward";
import LogoutIcon from "@mui/icons-material/Logout";
import Stack from "@mui/material/Stack";

const REQUIRE_AUTH = process.env.NEXT_PUBLIC_REQUIRE_AUTH === "true";

export default function HistoryToolbar(props) {
  const router = useRouter();
  const { instance } = useMsal();

  const handleLogout = () => {
    instance.logoutRedirect();
  };

  return (
    <Toolbar>
      <Stack direction="row" spacing={1} sx={{ mr: 2 }}>
        <IconButton
          aria-label="Back"
          onClick={() => router.back()}
          size="large"
        >
          <ArrowBackIcon />
        </IconButton>
        <IconButton
          aria-label="Forward"
          onClick={() => router.forward()}
          size="large"
        >
          <ArrowForwardIcon />
        </IconButton>
        {REQUIRE_AUTH && (
          <IconButton aria-label="Logout" onClick={handleLogout} size="large">
            <LogoutIcon />
          </IconButton>
        )}
      </Stack>
      {props.children}
    </Toolbar>
  );
}
