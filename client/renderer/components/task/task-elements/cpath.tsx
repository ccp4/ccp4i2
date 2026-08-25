import React, { useCallback, useEffect, useMemo, useState } from "react";
import {
  Avatar,
  Box,
  IconButton,
  Stack,
  TextField,
  Tooltip,
} from "@mui/material";
import { DeleteOutline, FolderOpen } from "@mui/icons-material";

import { CCP4i2TaskElementProps } from "./task-element";
import { useContainerField } from "./hooks/useContainerField";
import { ErrorTrigger } from "./error-info";
import { FIELD_SPACING } from "./field-sizes";
import {
  BrowseMode,
  browsePath,
  canBrowsePath,
} from "../../../utils/browse-path";

export interface CPathElementProps extends CCP4i2TaskElementProps {
  /** "directory" (the default) or "file" — which native chooser to open. */
  mode?: BrowseMode;
}

/**
 * A file parameter that names a path on the machine running the job, rather
 * than a file to upload into the project.
 *
 * Used for CDataFiles qualified `isDirectory` — a folder of diffraction images
 * for xia2, a BORGES library for ARCIMBOLDO, an existing run directory — and
 * for diffraction images themselves. None of the ordinary file sources apply:
 * these cannot be uploaded through the browser, are never project-database
 * files, and have nothing to fetch from the internet. They must also stay
 * where they are, because the program will read them in place — uploading one
 * frame of a sweep would separate it from the siblings that make it useful.
 *
 * So this offers a path: the native chooser in the desktop build, and a typed
 * absolute path in either build.
 */
export const CPathElement: React.FC<CPathElementProps> = (props) => {
  const { job, itemName, qualifiers, visibility, disabled, onChange } = props;
  const mode: BrowseMode = props.mode ?? "directory";

  const {
    item,
    unwrappedValue,
    isVisible,
    isDisabled,
    validationColor,
    commit,
  } = useContainerField({
    job,
    itemName,
    visibility,
    disabled,
    onChange,
    suppressMutations: props.suppressMutations,
  });

  // An external absolute path is stored whole in baseName with relPath
  // cleared; a path inside the project keeps the two separate.
  const serverPath = useMemo(() => {
    const baseName = unwrappedValue?.baseName ?? "";
    const relPath = unwrappedValue?.relPath ?? "";
    if (!baseName) return "";
    return relPath ? `${relPath}/${baseName}` : `${baseName}`;
  }, [unwrappedValue?.baseName, unwrappedValue?.relPath]);

  const [draft, setDraft] = useState(serverPath);
  useEffect(() => setDraft(serverPath), [serverPath]);

  const guiLabel =
    qualifiers?.guiLabel ||
    item?._objectPath?.split(".").at(-1) ||
    (mode === "file" ? "File" : "Directory");

  const commitDraft = useCallback(() => {
    const trimmed = draft.trim();
    if (trimmed === serverPath) return;
    commit(trimmed === "" ? null : trimmed);
  }, [draft, serverPath, commit]);

  const handleBrowse = useCallback(async () => {
    const picked = await browsePath({
      mode,
      title: `Select ${guiLabel}`,
      message: qualifiers?.toolTip || `Select ${guiLabel}`,
    });
    if (picked) {
      setDraft(picked);
      commit(picked);
    }
  }, [commit, guiLabel, mode, qualifiers?.toolTip]);

  if (!isVisible) return null;

  const browseEnabled = canBrowsePath() && !isDisabled;

  return (
    <Box sx={{ mx: FIELD_SPACING.marginLeft, my: 0 }}>
      <Stack direction="row" alignItems="center">
        <Avatar
          src="/svgicons/DataFile.svg"
          alt={mode === "file" ? "File" : "Directory"}
          sx={{
            width: 32,
            height: 32,
            mr: 1,
            flexShrink: 0,
            bgcolor: serverPath ? "primary.light" : "action.hover",
          }}
        >
          <FolderOpen fontSize="small" />
        </Avatar>

        <Box sx={{ flex: 1, minWidth: 0 }}>
          <Tooltip title={qualifiers?.toolTip || ""}>
            <TextField
              fullWidth
              size="small"
              label={guiLabel}
              value={draft}
              disabled={isDisabled}
              error={validationColor === "error.light"}
              onChange={(event) => setDraft(event.target.value)}
              onBlur={commitDraft}
              onKeyDown={(event) => {
                if (event.key === "Enter") commitDraft();
              }}
              placeholder={
                mode === "file"
                  ? "Absolute path to a file on the machine running the job"
                  : "Absolute path to a directory on the machine running the job"
              }
              slotProps={{
                inputLabel: { shrink: true, disableAnimation: true },
              }}
            />
          </Tooltip>
        </Box>

        <Stack
          direction="row"
          spacing={0.5}
          alignItems="center"
          sx={{ ml: 1, flexShrink: 0 }}
        >
          <Tooltip
            title={
              canBrowsePath()
                ? "Browse for a directory"
                : "The desktop app can open a directory chooser here; in the browser, type the path"
            }
          >
            <span>
              <IconButton
                size="small"
                onClick={handleBrowse}
                disabled={!browseEnabled}
                aria-label={`Browse for a ${mode}`}
              >
                <FolderOpen fontSize="small" />
              </IconButton>
            </span>
          </Tooltip>

          {serverPath && (
            <Tooltip title="Clear">
              <span>
                <IconButton
                  size="small"
                  onClick={() => {
                    setDraft("");
                    commit(null);
                  }}
                  disabled={isDisabled}
                  aria-label={`Clear ${mode}`}
                >
                  <DeleteOutline fontSize="small" />
                </IconButton>
              </span>
            </Tooltip>
          )}

          <ErrorTrigger item={item} job={job} />
        </Stack>
      </Stack>
    </Box>
  );
};

CPathElement.displayName = "CPathElement";
