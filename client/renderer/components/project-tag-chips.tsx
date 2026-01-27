"use client";

import React from "react";
import { Box, Chip, Tooltip, Typography } from "@mui/material";
import { ProjectTag } from "../types/models";

export interface ProjectTagChipsProps {
  /** Tags to display - can be array of tag objects or tag IDs */
  tags: (number | ProjectTag)[] | undefined;
  /** Maximum number of visible tags before showing "+N" overflow (default: 3) */
  maxVisible?: number;
  /** Chip size (default: "small") */
  size?: "small" | "medium";
  /** Hide "No tags" text when empty (default: false) */
  hideEmpty?: boolean;
}

/**
 * Displays project tags as styled chips.
 * Handles both old format (number[]) and new format (ProjectTag[]).
 * Shows overflow indicator when tags exceed maxVisible.
 */
export const ProjectTagChips = React.memo(
  ({
    tags,
    maxVisible = 3,
    size = "small",
    hideEmpty = false,
  }: ProjectTagChipsProps) => {
    // Handle both old format (number[]) and new format (ProjectTag[])
    const projectTagsData = Array.isArray(tags)
      ? tags.filter(
          (tag): tag is ProjectTag => typeof tag === "object" && tag !== null
        )
      : [];

    if (projectTagsData.length === 0) {
      if (hideEmpty) {
        return null;
      }
      // Return empty container to maintain consistent layout
      return (
        <Box
          sx={{
            minHeight: size === "small" ? 20 : 24,
            display: "flex",
            alignItems: "center",
          }}
        >
          <Typography
            variant="caption"
            color="text.disabled"
            sx={{ fontSize: size === "small" ? "0.65rem" : "0.7rem" }}
          >
            No tags
          </Typography>
        </Box>
      );
    }

    const visibleTags = projectTagsData.slice(0, maxVisible);
    const hiddenTags = projectTagsData.slice(maxVisible);
    const hiddenCount = hiddenTags.length;

    return (
      <Box
        sx={{
          display: "flex",
          flexWrap: "wrap",
          gap: 0.5,
          minHeight: size === "small" ? 20 : 24,
        }}
      >
        {visibleTags.map((tag) => (
          <Chip
            key={tag.id}
            label={tag.text}
            size={size}
            variant="outlined"
            sx={{
              height: size === "small" ? 20 : 24,
              fontSize: size === "small" ? "0.7rem" : "0.75rem",
              bgcolor: "primary.50",
              borderColor: "primary.200",
              color: "primary.700",
              fontWeight: 500,
              "&:hover": {
                bgcolor: "primary.100",
                borderColor: "primary.300",
              },
            }}
          />
        ))}
        {hiddenCount > 0 && (
          <Tooltip
            title={
              <Box>
                {hiddenTags.map((tag) => (
                  <Typography key={tag.id} variant="caption" display="block">
                    {tag.text}
                  </Typography>
                ))}
              </Box>
            }
            placement="top"
          >
            <Chip
              label={`+${hiddenCount}`}
              size={size}
              variant="outlined"
              sx={{
                height: size === "small" ? 20 : 24,
                fontSize: size === "small" ? "0.7rem" : "0.75rem",
                bgcolor: "grey.100",
                borderColor: "grey.300",
                color: "grey.600",
                fontWeight: 500,
                cursor: "help",
              }}
            />
          </Tooltip>
        )}
      </Box>
    );
  }
);

ProjectTagChips.displayName = "ProjectTagChips";
