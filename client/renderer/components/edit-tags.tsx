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
  const { data: allTags = [] } = api.get<ProjectTag[]>("projecttags/");
  const [inputValue, setInputValue] = useState("");

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
    (
      _event: React.SyntheticEvent,
      newValue: TagOption[],
      reason: string,
      details?: any
    ) => {
      if (reason === "selectOption" && details?.option) {
        const added = details.option as TagOption;
        if (added.id && !props.tags.includes(added.id)) {
          props.onChange([...props.tags, added.id]);
        }
        // Note: creating new tags inline is not supported here — tags must
        // exist before they can be applied to a new project. New tags can
        // be created from the project settings page after creation.
      } else if (reason === "removeOption" && details?.option) {
        const removed = details.option as TagOption;
        if (removed.id) {
          props.onChange(props.tags.filter((id) => id !== removed.id));
        }
      }
    },
    [props.tags, props.onChange]
  );

  const filterOptions = useCallback(
    (options: TagOption[], { inputValue }: any) =>
      options.filter(
        (option) =>
          !selectedTags.some((s) => s.id === option.id) &&
          option.text.toLowerCase().includes(inputValue.toLowerCase())
      ),
    [selectedTags]
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
