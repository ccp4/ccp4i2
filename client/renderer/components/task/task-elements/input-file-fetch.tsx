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
import { Button, SxProps } from "@mui/material";
import { ChangeEvent, useCallback } from "react";
import { useTaskInterface } from "../../../providers/task-provider";
import { Language } from "@mui/icons-material";

interface InputFileFetchProps {
  handleFileChange: (ev: ChangeEvent<HTMLInputElement>) => void;
  disabled: boolean;
  modes?: string[];
  item?: any;
  sx?: SxProps;
  onChange?: (updatedItem: any) => void;
}
export const InputFileFetch: React.FC<InputFileFetchProps> = ({
  disabled,
  sx,
  item,
  modes,
  onChange,
}) => {
  const { setDownloadDialogOpen, setFetchItemParams } = useTaskInterface();

  const handleFetchClick = useCallback(
    (ev: any) => {
      ev.stopPropagation();
      if (setFetchItemParams) setFetchItemParams({ item, modes, onChange });
      if (setDownloadDialogOpen) setDownloadDialogOpen(true);
    },
    [item, modes, onChange, setDownloadDialogOpen, setFetchItemParams]
  );

  return (
    <Button
      disabled={disabled}
      component="label"
      role={undefined}
      variant="outlined"
      tabIndex={-1}
      size="small"
      startIcon={<Language fontSize="small" />}
      sx={sx}
      onClick={handleFetchClick}
    />
  );
};
