import React from "react";
import {
  Box,
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
import { FieldShell } from "./field-shell";
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

  // An interface that gives the list no label gets its items' label
  // ("Atomic model") before the parameter's name, which is never meant for
  // users ("XYZIN_LIST", "DICT_LIST" and "UNMERGEDFILES" all showed).
  const guiLabel =
    qualifiers?.guiLabel ||
    item?._subItem?._qualifiers?.guiLabel ||
    item?._value?.[0]?._qualifiers?.guiLabel ||
    item?._objectPath?.split(".").at(-1) ||
    "Unnamed List";

  // The list's label names the list. An item that has a label of its own
  // ("Atomic model") keeps it, rather than every row repeating the list's
  // ("Models to superpose"); an item without one still takes the list's, so
  // it never falls back to its bare name ("[0]").
  const { guiLabel: _listLabel, ...qualifiersWithoutLabel } = qualifiers || {};
  const itemQualifiers = (content: any) =>
    content?._qualifiers?.guiLabel ? qualifiersWithoutLabel : qualifiers;

  const borderColor = itemName && item ? validationColor : "divider";

  if (!isVisible) return null;

  return (
    <FieldShell
      title={guiLabel}
      borderColor={borderColor}
      action={
        <Tooltip title="Add item">
          <span>
            <IconButton
              disabled={!isEditable}
              onClick={() => addItem()}
              size="small"
              color="primary"
              aria-label="Add new item to list"
            >
              <Add fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>
      }
      errorTrigger={<ErrorTrigger item={item} job={job} />}
    >
      <Box>
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
                      qualifiers={itemQualifiers(content)}
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
                      <Delete fontSize="small" />
                    </IconButton>
                  </Tooltip>
                </Stack>
              );
            })}
          </Stack>
        )}
      </Box>
    </FieldShell>
  );
};

CListElement.displayName = "CListElement";

export default CListElement;
