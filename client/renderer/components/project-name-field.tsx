"use client";
import { ChangeEvent } from "react";
import { TextField } from "@mui/material";
import { Project } from "../types/models";

export function getProjectNameError(
  name: string,
  projects?: Project[],
  excludeId?: number
): string {
  if (name.length === 0) return "Name is required";
  if (!name.match("^[A-z0-9_-]+$"))
    return "Name can only contain letters, numbers, underscores, and hyphens";
  if (projects?.find((p) => p.name === name && p.id !== excludeId))
    return "Name is already taken";
  return "";
}

interface ProjectNameFieldProps {
  value: string;
  onChange: (value: string) => void;
  /** Called when Enter is pressed and the field is valid */
  onSubmit?: () => void;
  error: string;
  autoFocus?: boolean;
}

export function ProjectNameField({
  value,
  onChange,
  onSubmit,
  error,
  autoFocus,
}: ProjectNameFieldProps) {
  function handleChange(e: ChangeEvent<HTMLInputElement>) {
    onChange(e.target.value);
  }

  function handleKeyDown(e: React.KeyboardEvent<HTMLInputElement>) {
    if (e.key === "Enter" && onSubmit && error.length === 0) {
      onSubmit();
    }
  }

  return (
    <TextField
      label="Name"
      value={value}
      onChange={handleChange}
      onKeyDown={handleKeyDown}
      required
      error={error.length > 0}
      helperText={error}
      autoFocus={autoFocus}
    />
  );
}
