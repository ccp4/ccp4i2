import React, { useState, useCallback } from "react";
import {
  Autocomplete,
  TextField,
  Chip,
  Box,
  Typography,
  CircularProgress,
  Paper,
} from "@mui/material";
import { Add as AddIcon, LocalOffer as TagIcon } from "@mui/icons-material";
import { useApi } from "../api";
import { usePopcorn } from "../providers/popcorn-provider";
import { ProjectTag } from "../types/models";
import { useProject } from "../utils";

interface TagOption {
  id?: number;
  text: string;
  isNew?: boolean;
  parent?: number | null;
}

export const TagsOfProject: React.FC<{
  projectId: number;
}> = ({ projectId }) => {
  const { project } = useProject(projectId);
  const api = useApi();
  const { setMessage } = usePopcorn();

  // Get current project tags
  const {
    data: projectTags = [],
    mutate: mutateProjectTags,
    isLoading: isLoadingProjectTags,
  } = api.get<ProjectTag[]>(`projects/${projectId}/tags/`);

  // Get all available tags for autocomplete
  const {
    data: allTags = [],
    mutate: mutateAllTags,
    isLoading: isLoadingAllTags,
  } = api.get<ProjectTag[]>("projecttags/");

  const [inputValue, setInputValue] = useState("");
  const [isSubmitting, setIsSubmitting] = useState(false);

  // Convert ProjectTags to TagOptions for the autocomplete
  const availableOptions: TagOption[] = React.useMemo(() => {
    return allTags.map((tag) => ({
      id: tag.id,
      text: tag.text,
      parent: tag.parent,
      isNew: false,
    }));
  }, [allTags]);

  // Current selected tags as TagOptions
  const selectedTags: TagOption[] = React.useMemo(() => {
    return projectTags.map((tag) => ({
      id: tag.id,
      text: tag.text,
      parent: tag.parent,
      isNew: false,
    }));
  }, [projectTags]);

  // Handle adding a tag to the project
  const handleAddTag = useCallback(
    async (tagOption: TagOption) => {
      if (!tagOption) return;

      setIsSubmitting(true);
      try {
        let tagId = tagOption.id;

        // If it's a new tag, create it first
        if (tagOption.isNew || !tagId) {
          console.log("Creating new tag:", {
            text: tagOption.text,
            parent: tagOption.parent || null,
          });

          try {
            const newTag = await api.post<ProjectTag>("projecttags/", {
              text: tagOption.text,
              parent: tagOption.parent || null,
              projects: [], // Start with empty projects array
            });
            console.log("New tag created:", newTag);
            tagId = newTag.id;
            // Refresh the all tags list
            mutateAllTags();
          } catch (createError) {
            console.error("Failed to create new tag:", createError);
            // Check if it's a uniqueness error - maybe the tag already exists
            if (
              createError instanceof Error &&
              createError.message.includes("unique")
            ) {
              // Try to find existing tag with same text
              const existingTag = allTags.find(
                (tag) =>
                  tag.text === tagOption.text &&
                  tag.parent === (tagOption.parent || null)
              );
              if (existingTag) {
                console.log("Found existing tag:", existingTag);
                tagId = existingTag.id;
              } else {
                throw createError;
              }
            } else {
              throw createError;
            }
          }
        }

        // Add the tag to the project
        console.log("Adding tag to project:", { tag_id: tagId });
        await api.post(`projects/${projectId}/tags/`, {
          tag_id: tagId,
        });

        // Refresh project tags
        mutateProjectTags();
        setMessage(`✅ Tag "${tagOption.text}" added to project`);
      } catch (error) {
        console.error("Failed to add tag:", error);
        let errorMessage = "Unknown error";
        if (error instanceof Error) {
          errorMessage = error.message;
        } else if (typeof error === "string") {
          errorMessage = error;
        } else if (error && typeof error === "object" && "message" in error) {
          errorMessage = (error as any).message;
        }
        setMessage(`❌ Failed to add tag: ${errorMessage}`);
      } finally {
        setIsSubmitting(false);
      }
    },
    [projectId, api, mutateProjectTags, mutateAllTags, setMessage]
  );

  // Handle removing a tag from the project
  const handleRemoveTag = useCallback(
    async (tagOption: TagOption) => {
      if (!tagOption.id) return;

      setIsSubmitting(true);
      try {
        await api.delete(`projects/${projectId}/tags/${tagOption.id}/`);
        mutateProjectTags();
        setMessage(`Tag "${tagOption.text}" removed from project`);
      } catch (error) {
        console.error("Failed to remove tag:", error);
        setMessage(`Failed to remove tag: ${error}`);
      } finally {
        setIsSubmitting(false);
      }
    },
    [projectId, api, mutateProjectTags, setMessage]
  );

  // Handle tag selection change
  const handleTagChange = useCallback(
    (
      event: React.SyntheticEvent,
      newValue: TagOption[],
      reason: string,
      details?: any
    ) => {
      if (reason === "selectOption" && details?.option) {
        // Adding a tag (existing or new via "+ Create" option)
        const addedTag = details.option;
        // Case-insensitive check to prevent duplicates
        if (
          !selectedTags.find(
            (tag) => tag.text.toLowerCase() === addedTag.text.toLowerCase()
          )
        ) {
          handleAddTag(addedTag);
        }
      } else if (reason === "removeOption" && details?.option) {
        // Removing a tag
        handleRemoveTag(details.option);
      }
    },
    [selectedTags, handleAddTag, handleRemoveTag]
  );

  // Filter options to show available tags that aren't already selected
  const filterOptions = useCallback(
    (options: TagOption[], { inputValue }: any) => {
      const filtered = options.filter(
        (option) =>
          !selectedTags.some((selected) => selected.text === option.text) &&
          option.text.toLowerCase().includes(inputValue.toLowerCase())
      );

      // Add "+ Create new tag" at the bottom when there's input text
      // and it's not already an exact match for an existing tag or selected tag
      const inputLower = inputValue.toLowerCase();
      const exactMatchExists =
        filtered.some((option) => option.text.toLowerCase() === inputLower) ||
        selectedTags.some((selected) => selected.text.toLowerCase() === inputLower);

      if (inputValue && !exactMatchExists) {
        filtered.push({
          text: inputValue,
          isNew: true,
        });
      }

      return filtered;
    },
    [selectedTags]
  );

  // Custom option rendering
  const renderOption = (props: any, option: TagOption) => (
    <Box
      component="li"
      {...props}
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

  // Custom tag rendering
  const renderTags = (tagValue: TagOption[], getTagProps: any) =>
    tagValue.map((option, index) => (
      <Chip
        {...getTagProps({ index })}
        key={option.id || option.text}
        label={option.text}
        size="small"
        disabled={isSubmitting}
      />
    ));

  if (isLoadingProjectTags || isLoadingAllTags) {
    return (
      <Box display="flex" alignItems="center" gap={1}>
        <TagIcon sx={{ fontSize: 20, color: "text.secondary" }} />
        <CircularProgress size={16} />
        <Typography variant="body2">Loading...</Typography>
      </Box>
    );
  }

  return (
    <Box display="flex" alignItems="center" gap={1} sx={{ minHeight: 40 }}>
      {/* Tag icon */}
      <TagIcon sx={{ fontSize: 20, color: "text.secondary" }} />

      {/* Count chip */}
      <Chip
        label={selectedTags.length}
        size="small"
        variant="outlined"
        sx={{
          minWidth: 32,
          height: 24,
          fontSize: "0.75rem",
          color: "text.secondary",
          borderColor: "divider",
        }}
      />

      {/* Autocomplete field */}
      <Box sx={{ flex: 1, minWidth: 200 }}>
        <Autocomplete
          multiple
          value={selectedTags}
          onChange={handleTagChange}
          inputValue={inputValue}
          onInputChange={(event, newInputValue) => setInputValue(newInputValue)}
          options={availableOptions}
          filterOptions={filterOptions}
          getOptionLabel={(option) => option.text}
          isOptionEqualToValue={(option, value) => option.text === value.text}
          renderOption={renderOption}
          renderTags={renderTags}
          loading={isSubmitting}
          disabled={isSubmitting}
          selectOnFocus
          clearOnBlur
          handleHomeEndKeys
          size="small"
          renderInput={(params) => (
            <TextField
              {...params}
              placeholder={
                selectedTags.length === 0 ? "Add tags..." : "Add more tags..."
              }
              variant="outlined"
              size="small"
              InputProps={{
                ...params.InputProps,
                sx: { height: 32 }, // Consistent height
                endAdornment: (
                  <>
                    {isSubmitting && <CircularProgress size={16} />}
                    {params.InputProps.endAdornment}
                  </>
                ),
              }}
            />
          )}
          PaperComponent={(props) => <Paper {...props} sx={{ mt: 1 }} />}
        />
      </Box>
    </Box>
  );
};
