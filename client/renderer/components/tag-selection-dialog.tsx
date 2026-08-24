"use client";
import React, { useMemo, useState } from "react";
import {
  Autocomplete,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  TextField,
  Typography,
} from "@mui/material";

import { useApi } from "../api";
import { ProjectTag } from "../types/models";

interface TagOption {
  id?: number;
  label: string;
  isNew?: boolean;
}

/**
 * Apply one tag to a whole selection of projects.
 *
 * Goes through the tag's bulk `add_projects` action rather than one request
 * per project, so tagging forty projects is one round trip.
 */
export default function TagSelectionDialog({
  open,
  projectIds,
  onClose,
  onApplied,
}: {
  open: boolean;
  projectIds: number[];
  onClose: () => void;
  onApplied: (label: string, count: number) => void;
}) {
  const api = useApi();
  const { data: allTags = [], mutate: mutateTags } =
    api.get<ProjectTag[]>("projecttags/");
  const [value, setValue] = useState<TagOption | null>(null);
  const [inputValue, setInputValue] = useState("");
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const options: TagOption[] = useMemo(
    () =>
      allTags.map((tag) => ({
        id: tag.id,
        label: tag.display_path ?? tag.text,
      })),
    [allTags]
  );

  const apply = async () => {
    if (!value) return;
    setBusy(true);
    setError(null);
    try {
      let tagId = value.id;
      if (value.isNew) {
        const created = await api.post<ProjectTag>("projecttags/", {
          text: value.label.trim(),
          parent: null,
          projects: [],
        });
        tagId = created.id;
        await mutateTags();
      }
      await api.post(`projecttags/${tagId}/add_projects/`, {
        project_ids: projectIds,
      });
      onApplied(value.label.trim(), projectIds.length);
      handleClose();
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    } finally {
      setBusy(false);
    }
  };

  const handleClose = () => {
    setValue(null);
    setInputValue("");
    setError(null);
    onClose();
  };

  return (
    <Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
      <DialogTitle>
        Tag {projectIds.length} project{projectIds.length === 1 ? "" : "s"}
      </DialogTitle>
      <DialogContent>
        <Autocomplete
          autoFocus
          value={value}
          onChange={(_event, option) => setValue(option)}
          inputValue={inputValue}
          onInputChange={(_event, next) => setInputValue(next)}
          options={options}
          getOptionLabel={(option) => option.label}
          isOptionEqualToValue={(option, selected) =>
            option.id === selected.id && option.label === selected.label
          }
          filterOptions={(available, state) => {
            const query = state.inputValue.trim();
            const filtered = available.filter((option) =>
              option.label.toLowerCase().includes(query.toLowerCase())
            );
            const exists = available.some(
              (option) => option.label.toLowerCase() === query.toLowerCase()
            );
            if (query && !exists) {
              filtered.push({ label: query, isNew: true });
            }
            return filtered;
          }}
          renderOption={(props, option) => (
            <li {...props} key={`${option.id ?? "new"}-${option.label}`}>
              {option.isNew ? `Create tag "${option.label}"` : option.label}
            </li>
          )}
          renderInput={(params) => (
            <TextField
              {...params}
              label="Tag"
              placeholder="Choose or type a new tag"
              margin="dense"
            />
          )}
        />
        <Typography variant="caption" color="text.secondary">
          A project keeps any tags it already has — this adds one more.
        </Typography>
        {error && (
          <Typography variant="caption" color="error" sx={{ display: "block" }}>
            {error}
          </Typography>
        )}
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose}>Cancel</Button>
        <Button onClick={apply} disabled={!value || busy} variant="contained">
          Apply
        </Button>
      </DialogActions>
    </Dialog>
  );
}
