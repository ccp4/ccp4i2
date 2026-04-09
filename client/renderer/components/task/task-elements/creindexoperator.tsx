import React, { useMemo } from "react";
import { Box } from "@mui/material";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob } from "../../../utils";
import { FieldRow } from "./field-row";

/**
 * CReindexOperatorElement - displays reindex operator matrix values in a row.
 *
 * DESIGN: Base widgets are full-width. This container constrains them to 'xs' size.
 */
export const CReindexOperatorElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { job, itemName, visibility } = props;
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
      {/* Container controls width - 'xs' gives ~8rem per field, good for small integers */}
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
