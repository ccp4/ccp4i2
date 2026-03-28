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
import { Editor, EditorProps } from "@monaco-editor/react";
import { useApi } from "../api";
import { useJob } from "../utils";
import { Stack } from "@mui/material";
import { use, useEffect, useRef, useState } from "react";
import { useTheme } from "../theme/theme-provider";

export const JobCommentEditor: React.FC<{
  jobId: number;
}> = ({ jobId }) => {
  const api = useApi();
  const { job } = useJob(jobId);
  const timerRef = useRef<any | null>(null);
  const [inFlight, setInFlight] = useState<boolean>(false);
  const { mode } = useTheme();

  useEffect(() => {
    return () => {
      if (timerRef.current) {
        clearTimeout(timerRef.current);
      }
    };
  }, []);

  const handleChange = (value: string, event: any) => {
    if (timerRef.current) {
      clearTimeout(timerRef.current);
    }

    timerRef.current = setTimeout(() => {
      saveComments(value);
    }, 2000);
  };

  const saveComments = async (value: string) => {
    setInFlight(true);

    try {
      const formData = new FormData();
      formData.append("comments", value);

      const response = await api.patch(`jobs/${jobId}`, formData);
    } catch (error: any) {
      if (error?.name === "Canceled") {
        console.warn("Monaco operation was canceled; ignoring.");
      } else {
        console.error("Error updating comments", error);
      }
    } finally {
      setInFlight(false);
    }
  };

  return (
    <Stack direction={"column"} spacing={1} sx={{ height: "100%" }}>
      <Editor
        height="calc(100vh - 15rem)"
        value={job?.comments || ""}
        language="text"
        onChange={handleChange}
        theme={mode === "dark" ? "vs-dark" : "light"}
        options={{ readOnly: inFlight }}
      />
    </Stack>
  );
};
