import React, { useMemo } from "react";
import { Box } from "@mui/material";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { FieldRow } from "./field-row";
import { FIELD_SPACING } from "./field-sizes";

/**
 * CCellElement - displays unit cell parameters (a, b, c, alpha, beta, gamma) in a row.
 *
 * DESIGN: Base widgets are full-width. This container constrains them to 'xs' size.
 * The FieldRow with size="xs" wraps each child in a container with maxWidth: 8rem.
 *
 * Before (old approach):
 *   - CCellElement passed elementSx={{ width: FIELD_SIZES.xs }} to children
 *   - Each base widget had to accept and apply elementSx
 *   - Width logic was duplicated in base widgets and containers
 *
 * After (new approach):
 *   - Base widgets are full-width by default
 *   - FieldRow controls children's width via size="xs"
 *   - Single point of width control in the container
 */
export const CCellElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { job, itemName, visibility, qualifiers } = props;
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  const isVisible = useMemo(() => {
    if (!visibility) return true;
    if (typeof visibility === "function") return visibility();
    return visibility;
  }, [visibility]);

  // Get child names from _CONTENTS_ORDER or _value keys
  const childNames = useMemo(() => {
    if (!item) return [];

    if (Array.isArray(item._CONTENTS_ORDER) && item._CONTENTS_ORDER.length > 0) {
      return item._CONTENTS_ORDER;
    }

    if (item._value && typeof item._value === "object") {
      return Object.keys(item._value);
    }

    return [];
  }, [item]);

  if (!isVisible || !item) return null;

  return (
    <Box
      sx={{
        mx: 2,
        px: 2,
        py: 1,
        border: 2,
        borderColor: "divider",
        borderRadius: 1,
      }}
    >
      {/* Container controls width - children are full-width and get constrained here */}
      <FieldRow equalWidth={false} size="xs">
        {childNames.map((childName: string) => {
          const childObjectPath = `${item._objectPath}.${childName}`;
          return (
            <CCP4i2TaskElement
              key={childObjectPath}
              {...props}
              itemName={childObjectPath}
            />
          );
        })}
      </FieldRow>
    </Box>
  );
};
