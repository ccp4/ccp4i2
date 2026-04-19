import React, { useMemo } from "react";
import {
  Box,
  Card,
  CardContent,
  CardHeader,
  IconButton,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import { Add, Delete } from "@mui/icons-material";

import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useProject } from "../../../utils";
import { useCCP4i2Window } from "../../../app-context";
import { ErrorTrigger } from "./error-info";
import { useContainerList } from "./hooks/useContainerList";

interface CListElementProps extends CCP4i2TaskElementProps {
  initiallyOpen?: boolean;
}

export const CListElement: React.FC<CListElementProps> = ({
  itemName,
  job,
  qualifiers,
  onChange,
  visibility,
  ...restProps
}) => {
  const { projectId } = useCCP4i2Window();
  const { project } = projectId
    ? useProject(projectId)
    : { project: undefined };

  const {
    item,
    items,
    isVisible,
    isEditable,
    validationColor,
    addItem,
    deleteAt,
  } = useContainerList({
    job,
    itemName,
    project,
    visibility,
    onChange,
  });

  const guiLabel =
    qualifiers?.guiLabel ||
    item?._objectPath?.split(".").at(-1) ||
    "Unnamed List";

  const borderColor = itemName && item ? validationColor : "divider";

  const cardSx = useMemo(
    () => ({
      mx: 1,
      border: 1,
      borderColor,
      borderRadius: 2,
      boxShadow: "none",
      "&:hover": {
        borderColor:
          borderColor === "divider" ? "primary.light" : borderColor,
      },
    }),
    [borderColor]
  );

  if (!isVisible) return null;

  return (
    <Card sx={cardSx}>
      <CardHeader
        title={<Typography variant="body2">{guiLabel}</Typography>}
        action={
          <Stack direction="row" alignItems="center">
            <IconButton
              disabled={!isEditable}
              onClick={() => addItem()}
              size="small"
              sx={{ color: "primary.text" }}
              aria-label="Add new item to list"
            >
              <Add />
            </IconButton>
            <ErrorTrigger item={item} job={job} />
          </Stack>
        }
      />
      <CardContent>
        {items.length === 0 ? (
          <Typography
            variant="caption"
            color="text.secondary"
            sx={{ fontStyle: "italic", textAlign: "center" }}
          >
            No elements in this list
          </Typography>
        ) : (
          <Stack spacing={1}>
            {items.map((content: any, index: number) => {
              // During add/delete round-trips the cache may briefly expose
              // raw values (string/null) before the server response is
              // patched in.  Skip those rather than crashing; they'll be
              // replaced on the next render.
              const itemPath =
                content && typeof content === "object"
                  ? content._objectPath
                  : null;
              if (!itemPath) return null;
              return (
                <Stack
                  key={itemPath}
                  direction="row"
                  alignItems="center"
                  spacing={0.5}
                >
                  <Box sx={{ flex: 1, minWidth: 0 }}>
                    <CCP4i2TaskElement
                      {...restProps}
                      itemName={itemPath}
                      job={job}
                      qualifiers={qualifiers}
                      onChange={onChange}
                    />
                  </Box>
                  <Tooltip title={`Delete item ${index + 1}`} placement="left">
                    <IconButton
                      disabled={!isEditable}
                      onClick={() => deleteAt(index)}
                      size="small"
                      color="error"
                      aria-label={`Delete item ${index + 1}`}
                    >
                      <Delete />
                    </IconButton>
                  </Tooltip>
                </Stack>
              );
            })}
          </Stack>
        )}
      </CardContent>
    </Card>
  );
};

CListElement.displayName = "CListElement";

export default CListElement;
