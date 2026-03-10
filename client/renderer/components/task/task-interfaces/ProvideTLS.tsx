import React, { useCallback, useEffect, useState } from "react";
import {
  Box,
  Button,
  LinearProgress,
  Paper,
  ToggleButton,
  ToggleButtonGroup,
  Typography,
} from "@mui/material";
import { AutoFixHigh, Upload } from "@mui/icons-material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

type EditMode = "table" | "text";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, useTaskItem, setParameter, callPluginMethod } =
    useJob(props.job.id);
  const { item: tlsGroupsItem } = useTaskItem("TLSGROUPS");
  const { item: editModeItem } = useTaskItem("EDIT_MODE");
  const [suggesting, setSuggesting] = useState(false);
  const [importing, setImporting] = useState(false);
  const [statusMsg, setStatusMsg] = useState<string | null>(null);
  const [editMode, setEditMode] = useState<EditMode>("table");

  // Sync local state from server parameter
  useEffect(() => {
    const serverVal = editModeItem?._value;
    if (serverVal === "text" || serverVal === "table") {
      setEditMode(serverVal);
    }
  }, [editModeItem?._value]);

  const isEditable = props.job.status === 1;

  /** Toggle edit mode — persists to server via guiParameters.EDIT_MODE */
  const handleEditModeChange = useCallback(
    async (_e: React.MouseEvent, newMode: EditMode | null) => {
      if (!newMode || !editModeItem) return;
      setEditMode(newMode);
      await setParameter({
        object_path: editModeItem._objectPath,
        value: newMode,
      });
    },
    [editModeItem, setParameter]
  );

  /** Suggest TLS groups from XYZIN coordinates */
  const handleSuggest = useCallback(async () => {
    if (!isEditable || suggesting) return;
    setSuggesting(true);
    setStatusMsg(null);
    try {
      const result = await callPluginMethod("suggest_tls_groups");
      if (Array.isArray(result) && result.length > 0 && tlsGroupsItem) {
        await setParameter({
          object_path: tlsGroupsItem._objectPath,
          value: result,
        });
        setStatusMsg(`Suggested ${result.length} TLS range(s) from coordinates.`);
      } else if (result?.error) {
        setStatusMsg(result.error);
      } else {
        setStatusMsg("No polymer chains found in coordinate file.");
      }
    } catch (e) {
      console.error("Error suggesting TLS groups:", e);
      setStatusMsg("Failed to suggest TLS groups.");
    } finally {
      setSuggesting(false);
    }
  }, [isEditable, suggesting, callPluginMethod, tlsGroupsItem, setParameter]);

  /** Import from TLSTEXT into structured table */
  const handleImportFromText = useCallback(async () => {
    if (!isEditable || importing) return;
    setImporting(true);
    setStatusMsg(null);
    try {
      const result = await callPluginMethod("import_tls_from_text");
      if (Array.isArray(result) && result.length > 0 && tlsGroupsItem) {
        await setParameter({
          object_path: tlsGroupsItem._objectPath,
          value: result,
        });
        setStatusMsg(`Imported ${result.length} TLS range(s) from text.`);
      } else if (result?.error) {
        setStatusMsg(result.error);
      } else {
        setStatusMsg("No TLS definitions found in text.");
      }
    } catch (e) {
      console.error("Error importing TLS text:", e);
      setStatusMsg("Failed to import TLS text.");
    } finally {
      setImporting(false);
    }
  }, [isEditable, importing, callPluginMethod, tlsGroupsItem, setParameter]);

  /** Import from TLSIN file into structured table */
  const handleImportFromFile = useCallback(async () => {
    if (!isEditable || importing) return;
    setImporting(true);
    setStatusMsg(null);
    try {
      const result = await callPluginMethod("import_tls_from_file");
      if (Array.isArray(result) && result.length > 0 && tlsGroupsItem) {
        await setParameter({
          object_path: tlsGroupsItem._objectPath,
          value: result,
        });
        setStatusMsg(`Imported ${result.length} TLS range(s) from file.`);
      } else if (result?.error) {
        setStatusMsg(result.error);
      } else {
        setStatusMsg("No TLS definitions found in file.");
      }
    } catch (e) {
      console.error("Error importing TLS file:", e);
      setStatusMsg("Failed to import TLS file.");
    } finally {
      setImporting(false);
    }
  }, [isEditable, importing, callPluginMethod, tlsGroupsItem, setParameter]);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input files for import/suggestion */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" color="text.secondary">
          Provide a coordinate file to suggest TLS groups, or a TLS file to import existing definitions.
        </Typography>
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <CCP4i2TaskElement
          itemName="TLSIN"
          {...props}
          qualifiers={{ guiLabel: "Starting TLS set" }}
        />
      </CCP4i2ContainerElement>

      {/* TLS definitions — toggle between table and raw text */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "TLS group definitions" }}
        containerHint="FolderLevel"
      >
        <Box sx={{ display: "flex", alignItems: "center", gap: 2, mb: 1 }}>
          <ToggleButtonGroup
            value={editMode}
            exclusive
            onChange={handleEditModeChange}
            size="small"
          >
            <ToggleButton value="table">Table</ToggleButton>
            <ToggleButton value="text">Free text</ToggleButton>
          </ToggleButtonGroup>

          {editMode === "table" && (
            <>
              <Button
                variant="outlined"
                size="small"
                startIcon={<AutoFixHigh />}
                onClick={handleSuggest}
                disabled={!isEditable || suggesting}
              >
                {suggesting ? "Suggesting..." : "Suggest from coordinates"}
              </Button>
              <Button
                variant="outlined"
                size="small"
                startIcon={<Upload />}
                onClick={handleImportFromFile}
                disabled={!isEditable || importing}
              >
                Import from TLS file
              </Button>
              <Button
                variant="outlined"
                size="small"
                startIcon={<Upload />}
                onClick={handleImportFromText}
                disabled={!isEditable || importing}
              >
                Import from TLS text
              </Button>
            </>
          )}
        </Box>

        {statusMsg && (
          <Typography variant="body2" color="text.secondary" sx={{ mb: 1 }}>
            {statusMsg}
          </Typography>
        )}

        {editMode === "table" ? (
          <CCP4i2TaskElement itemName="TLSGROUPS" {...props} />
        ) : (
          <CCP4i2TaskElement itemName="TLSTEXT" {...props} />
        )}
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
