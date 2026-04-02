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
import {
  Box,
  Card,
  CardContent,
  CardHeader,
  Collapse,
  Divider,
  Stack,
  SxProps,
  Typography,
} from "@mui/material";
import React, { PropsWithChildren, useMemo, useState } from "react";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { ErrorInfo } from "./error-info";
import { MyExpandMore } from "../../expand-more";
import { ExpandMore } from "@mui/icons-material";
import { FIELD_SPACING } from "./field-sizes";
import { useExpertLevel } from "./expert-level-context";

interface CCP4i2ContainerElementProps extends CCP4i2TaskElementProps {
  initiallyOpen?: boolean;
  containerHint?: "FolderLevel" | "BlockLevel" | "RowLevel";
  /**
   * @deprecated Use dedicated composite widgets (CCellElement, CReindexOperatorElement)
   * or wrap children in FieldRow with size prop instead.
   *
   * Base widgets are now full-width by default; containers control sizing.
   */
  elementSx?: SxProps;
  excludeItems?: string[];
}

/**
 * Vertical stacking container - each field gets its own line.
 * This is the default layout for FolderLevel and BlockLevel containers.
 */
const COLUMN_CONTAINER_SX = {
  display: "flex",
  flexDirection: "column",
  gap: FIELD_SPACING.rowGap,
  alignItems: "stretch",
  pl: FIELD_SPACING.marginLeft,
  py: 1,
} as const;

/**
 * Horizontal row container - fields flow in a row with wrapping.
 * Used for RowLevel containers where fields should be side-by-side.
 */
const ROW_CONTAINER_SX = {
  display: "flex",
  flexDirection: "row",
  flexWrap: "wrap",
  gap: FIELD_SPACING.columnGap,
  alignItems: "center",
  pl: FIELD_SPACING.marginLeft,
  py: 1,
} as const;

export const CCP4i2ContainerElement: React.FC<
  PropsWithChildren<CCP4i2ContainerElementProps>
> = (props) => {
  const {
    job,
    itemName,
    children,
    containerHint = "FolderLevel",
    initiallyOpen = true,
    visibility,
    qualifiers,
    elementSx,
  } = props;

  const { useTaskItem, getValidationColor } = useJob(job.id);
  const { item } = useTaskItem(itemName);
  const [open, setOpen] = useState(initiallyOpen);
  const maxExpertLevel = useExpertLevel();

  const inferredVisibility = useMemo(() => {
    if (!visibility) return true;
    if (typeof visibility === "function") {
      return visibility();
    }
    return visibility;
  }, [visibility]);

  // Get validation color for border when itemName is provided
  const validationBorderColor = useMemo(() => {
    if (itemName && item) {
      return getValidationColor(item);
    }
    return "divider";
  }, [itemName, item, getValidationColor]);

  const childNames = useMemo(() => {
    if (item) {
      let names: string[] = [];

      if (
        Array.isArray(item?._CONTENTS_ORDER) &&
        item._CONTENTS_ORDER.length > 0
      ) {
        names = item._CONTENTS_ORDER;
      } else if (item._value && item._value.constructor == Object) {
        names = Object.keys(item._value);
      }

      // Filter out excluded items if excludeItems prop is provided
      if (props.excludeItems && props.excludeItems.length > 0) {
        names = names.filter((name) => !props.excludeItems?.includes(name));
      }

      return names;
    }
    return [];
  }, [item, props.excludeItems]);

  // Generate content for column layout (default - each field on its own line)
  const columnContent = useMemo(() => {
    return item ? (
      <Box sx={COLUMN_CONTAINER_SX} key={item._objectPath}>
        {childNames.map((childName: string) => {
          const childObjectPath = `${item._objectPath}.${childName}`;
          const { item: childItem } = useTaskItem(childObjectPath);

          // Expert level filtering: hide children above the threshold
          if (maxExpertLevel !== undefined) {
            const childLevel = childItem?._qualifiers?.expertLevel;
            if (childLevel !== undefined && childLevel !== null && childLevel > maxExpertLevel) {
              return null;
            }
          }

          return (
            <CCP4i2TaskElement
              key={childObjectPath}
              {...props}
              sx={elementSx}
              itemName={childObjectPath}
              qualifiers={{ ...childItem._qualifiers }}
            />
          );
        })}
      </Box>
    ) : null;
  }, [item, elementSx, childNames, useTaskItem, props, maxExpertLevel]);

  // Generate content for row layout (RowLevel - fields side-by-side)
  const rowContent = useMemo(() => {
    return item ? (
      <Box sx={ROW_CONTAINER_SX} key={item._objectPath}>
        {childNames.map((childName: string) => {
          const childObjectPath = `${item._objectPath}.${childName}`;
          const { item: childItem } = useTaskItem(childObjectPath);

          // Expert level filtering: hide children above the threshold
          if (maxExpertLevel !== undefined) {
            const childLevel = childItem?._qualifiers?.expertLevel;
            if (childLevel !== undefined && childLevel !== null && childLevel > maxExpertLevel) {
              return null;
            }
          }

          return (
            <CCP4i2TaskElement
              key={childObjectPath}
              {...props}
              sx={elementSx}
              itemName={childObjectPath}
              qualifiers={{ ...childItem._qualifiers }}
            />
          );
        })}
      </Box>
    ) : null;
  }, [item, elementSx, childNames, useTaskItem, props, maxExpertLevel]);

  const columnChildren = useMemo(() => {
    if (children) {
      return <Box sx={COLUMN_CONTAINER_SX}>{children}</Box>;
    }
    return null;
  }, [children]);

  const rowChildren = useMemo(() => {
    if (children) {
      return <Box sx={ROW_CONTAINER_SX}>{children}</Box>;
    }
    return null;
  }, [children]);

  // Row-level styling - horizontal flow with tighter spacing
  const rowSx = useMemo(
    () => ({
      display: "flex",
      flexWrap: "wrap",
      alignItems: "center",
      gap: 1,
      mx: 2,
      px: 2,
      py: 1,
      border: 2,
      borderColor: validationBorderColor,
      borderRadius: 1,
    }),
    [validationBorderColor]
  );

  if (!inferredVisibility) return null;

  if (containerHint === "FolderLevel") {
    return (
      <Card sx={{ m: 2 }}>
        <CardHeader
          sx={{ py: 1 }}
          title={
            <Typography variant="subtitle2" sx={{ color: "primary.main", fontWeight: 700 }}>
              {qualifiers?.guiLabel}
            </Typography>
          }
          onClick={(ev) => {
            ev.stopPropagation();
            setOpen(!open);
          }}
          variant="lightGrey"
          action={
            <Stack direction="row">
              <MyExpandMore expand={open} aria-expanded={open} aria-label="show more">
                <ExpandMore sx={{ color: "text.primary" }} />
              </MyExpandMore>
              {item && <ErrorInfo {...props} />}
            </Stack>
          }
        />
        <CardContent sx={{ px: 0, pt: 0, "&:last-child": { pb: 1 } }}>
          <Collapse in={open} timeout="auto" unmountOnExit>
            {columnContent}
            {columnChildren}
          </Collapse>
        </CardContent>
      </Card>
    );
  }

  if (containerHint === "BlockLevel") {
    return (
      <Box sx={{ mx: 2, my: 1 }}>
        <Typography
          variant="subtitle2"
          sx={{
            pl: FIELD_SPACING.marginLeft,
            mb: 0.5,
            fontWeight: 700,
            color: "primary.main",
          }}
        >
          {qualifiers?.guiLabel}
        </Typography>
        <Divider sx={{ borderColor: "primary.main", mb: 0.5 }} />
        {columnContent}
        {columnChildren}
      </Box>
    );
  }

  if (containerHint === "RowLevel") {
    return (
      <Box sx={rowSx}>
        {rowContent}
        {rowChildren}
      </Box>
    );
  }

  // Default: no container chrome, use column layout
  return (
    <>
      {columnContent}
      {columnChildren}
    </>
  );
};
