import React, { useCallback, useState } from "react";
import {
  Autocomplete,
  Box,
  Chip,
  Paper,
  TextField,
} from "@mui/material";
import { Add as AddIcon } from "@mui/icons-material";
import { useApi } from "../api";
import { ProjectTag } from "../types/models";

interface TagOption {
  id?: number;
  text: string;
  isNew?: boolean;
}

export default function EditTags(props: {
  tags: number[];
  onChange: (tags: number[]) => void;
}) {
  const api = useApi();
  const { data: allTags = [], mutate: mutateAllTags } =
    api.get<ProjectTag[]>("projecttags/");
  const [inputValue, setInputValue] = useState("");
  const [isCreating, setIsCreating] = useState(false);

  // Convert available tags to options
  const availableOptions: TagOption[] = React.useMemo(
    () => allTags.map((tag) => ({ id: tag.id, text: tag.text })),
    [allTags]
  );

  // Current selected tags as TagOptions
  const selectedTags: TagOption[] = React.useMemo(
    () =>
      props.tags
        .map((id) => allTags.find((t) => t.id === id))
        .filter((t): t is ProjectTag => t !== undefined)
        .map((t) => ({ id: t.id, text: t.text })),
    [props.tags, allTags]
  );

  const handleChange = useCallback(
    async (
      _event: React.SyntheticEvent,
      _newValue: TagOption[],
      reason: string,
      details?: any
    ) => {
      if (reason === "selectOption" && details?.option) {
        const added = details.option as TagOption;

        if (added.isNew) {
          // Create the tag via API, then add its ID
          setIsCreating(true);
          try {
            const newTag = await api.post<ProjectTag>("projecttags/", {
              text: added.text.trim(),
              parent: null,
              projects: [],
            });
            mutateAllTags();
            props.onChange([...props.tags, newTag.id]);
          } catch (err) {
            console.error("Failed to create tag:", err);
          } finally {
            setIsCreating(false);
          }
        } else if (added.id && !props.tags.includes(added.id)) {
          props.onChange([...props.tags, added.id]);
        }
      } else if (reason === "removeOption" && details?.option) {
        const removed = details.option as TagOption;
        if (removed.id) {
          props.onChange(props.tags.filter((id) => id !== removed.id));
        }
      }
    },
    [props.tags, props.onChange, api, mutateAllTags]
  );

  const filterOptions = useCallback(
    (options: TagOption[], { inputValue }: any) => {
      const trimmed = inputValue.trim();
      const inputLower = trimmed.toLowerCase();

      const filtered = options.filter(
        (option) =>
          !selectedTags.some((s) => s.id === option.id) &&
          option.text.toLowerCase().includes(inputLower)
      );

      // Offer to create a new tag if the input doesn't match any existing
      // tag (case-insensitive, whitespace-trimmed)
      if (trimmed) {
        const exactMatch =
          allTags.some(
            (t) => t.text.trim().toLowerCase() === inputLower
          ) ||
          selectedTags.some(
            (s) => s.text.trim().toLowerCase() === inputLower
          );

        if (!exactMatch) {
          filtered.push({ text: trimmed, isNew: true });
        }
      }

      return filtered;
    },
    [selectedTags, allTags]
  );

  const renderOption = (renderProps: any, option: TagOption) => (
    <Box
      component="li"
      {...renderProps}
      sx={
        option.isNew
          ? {
              borderTop: "1px solid",
              borderColor: "divider",
              color: "primary.main",
              fontWeight: 500,
            }
          : undefined
      }
    >
      {option.isNew && <AddIcon sx={{ mr: 1, fontSize: 16 }} />}
      {option.isNew ? `Create "${option.text}"` : option.text}
    </Box>
  );

  return (
    <Autocomplete
      multiple
      value={selectedTags}
      onChange={handleChange}
      inputValue={inputValue}
      onInputChange={(_event, value) => setInputValue(value)}
      options={availableOptions}
      filterOptions={filterOptions}
      getOptionLabel={(option) => option.text}
      isOptionEqualToValue={(a, b) => a.id === b.id}
      renderOption={renderOption}
      renderTags={(value, getTagProps) =>
        value.map((option, index) => (
          <Chip
            {...getTagProps({ index })}
            key={option.id}
            label={option.text}
            size="small"
          />
        ))
      }
      loading={isCreating}
      disabled={isCreating}
      selectOnFocus
      clearOnBlur
      handleHomeEndKeys
      size="small"
      renderInput={(params) => (
        <TextField
          {...params}
          label="Tags"
          placeholder={
            selectedTags.length === 0 ? "Add tags..." : "Add more..."
          }
          size="small"
        />
      )}
      PaperComponent={(paperProps) => <Paper {...paperProps} sx={{ mt: 1 }} />}
    />
  );
}
