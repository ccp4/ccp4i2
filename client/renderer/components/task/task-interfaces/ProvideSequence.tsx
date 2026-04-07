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
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback } from "react";
import { apiText, apiGet } from "../../../api-fetch";
import { useApi } from "../../../api";
import {
  Table, TableBody, TableCell, TableContainer, TableHead, TableRow,
  Paper, Typography, Alert, Tooltip,
} from "@mui/material";
import { Warning as WarningIcon } from "@mui/icons-material";

interface SequencePreview {
  id: string;
  name: string;
  description: string;
  sequence: string;
  length: number;
}

interface PreviewResult {
  sequences: SequencePreview[];
  format: string | null;
  commentary: string;
}

/** Return a list of warnings for a parsed sequence entry. */
function sequenceWarnings(seq: SequencePreview): string[] {
  const warnings: string[] = [];
  if (!seq.id || seq.id === "<unknown id>")
    warnings.push("No sequence identifier");
  if (!seq.description || seq.description === "<unknown description>")
    warnings.push("No description");
  if (seq.length < 10) warnings.push(`Very short (${seq.length} residues)`);
  return warnings;
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const api = useApi();
  const { useTaskItem } = useJob(job.id);

  const { value: sequenceText, forceUpdate: forceSetSEQUENCETEXT } = useTaskItem("SEQUENCETEXT");

  // Handle XYZIN coordinate file change - digest file and populate SEQUENCETEXT
  const handleXyzInChange = useCallback(
    async (updatedItem: any) => {
      const dbFileId =
        updatedItem?._value?.dbFileId?._value?.trim() ||
        updatedItem?.dbFileId;
      if (!dbFileId) return;

      try {
        const response = await apiGet(`files/${dbFileId}/digest_by_uuid`);
        const sequences: Record<string, string> = response?.data?.sequences;
        if (sequences && Object.keys(sequences).length > 0) {
          const fastaText = Object.entries(sequences)
            .map(([chainId, seq]) => `>${chainId}\n${seq}`)
            .join("\n");
          await forceSetSEQUENCETEXT(fastaText);
        }
      } catch (error) {
        console.error("Error extracting sequences from coordinate file:", error);
      }
    },
    [forceSetSEQUENCETEXT]
  );

  // Handle SEQIN file change - read raw file content and populate SEQUENCETEXT
  // Matches Qt behavior: open file, read content, set SEQUENCETEXT
  const handleSeqInChange = useCallback(
    async (updatedItem: any) => {
      // dbFileId is a CData object with nested _value, not a plain string
      const dbFileId =
        updatedItem?._value?.dbFileId?._value?.trim() ||
        updatedItem?.dbFileId;
      if (!dbFileId) return;

      try {
        const content = await apiText(`files/${dbFileId}/download_by_uuid`);
        if (content) {
          await forceSetSEQUENCETEXT(content);
        }
      } catch (error) {
        console.error("Error reading sequence file:", error);
      }
    },
    [forceSetSEQUENCETEXT]
  );

  // Preview: parse SEQUENCETEXT on the server and show what will be produced
  const hasText = typeof sequenceText === "string" && sequenceText.trim().length > 0;
  const { data: preview } = api.objectMethod<PreviewResult>(
    job.id,
    "ProvideSequence",
    "previewSequences",
    {},
    [sequenceText],
    hasText,
  );
  const previewData = (preview as any)?.data?.result as PreviewResult | undefined;

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Sequence" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="SEQUENCETEXT"
            qualifiers={{ guiLabel: "Sequence text" }}
            sx={{ minWidth: "100%", minHeight: "10rem" }}
          />

          {hasText && previewData && (
            <>
              {previewData.format && (
                <Typography variant="body2" sx={{ mb: 1 }}>
                  Detected format: <strong>{previewData.format}</strong>
                </Typography>
              )}

              {previewData.sequences.length === 0 ? (
                <Alert severity="warning" sx={{ mb: 1 }}>
                  Could not parse any sequences from the text provided.
                </Alert>
              ) : (
                <TableContainer component={Paper} variant="outlined" sx={{ mb: 1 }}>
                  <Table size="small">
                    <TableHead>
                      <TableRow>
                        <TableCell>#</TableCell>
                        <TableCell>ID</TableCell>
                        <TableCell>Description</TableCell>
                        <TableCell align="right">Length</TableCell>
                        <TableCell>Sequence</TableCell>
                      </TableRow>
                    </TableHead>
                    <TableBody>
                      {previewData.sequences.map((seq, i) => {
                        const warnings = sequenceWarnings(seq);
                        return (
                          <TableRow
                            key={i}
                            sx={warnings.length > 0 ? { bgcolor: "warning.light", "& td": { color: "warning.contrastText" } } : undefined}
                          >
                            <TableCell>
                              {warnings.length > 0 ? (
                                <Tooltip title={warnings.join("; ")}>
                                  <WarningIcon fontSize="small" color="warning" sx={{ verticalAlign: "middle" }} />
                                </Tooltip>
                              ) : (
                                i + 1
                              )}
                            </TableCell>
                            <TableCell>{seq.id}</TableCell>
                            <TableCell>{seq.description}</TableCell>
                            <TableCell align="right">{seq.length}</TableCell>
                            <TableCell
                              sx={{
                                fontFamily: "monospace",
                                fontSize: "0.75rem",
                                maxWidth: 300,
                                overflow: "hidden",
                                textOverflow: "ellipsis",
                                whiteSpace: "nowrap",
                              }}
                            >
                              {seq.sequence}
                            </TableCell>
                          </TableRow>
                        );
                      })}
                    </TableBody>
                  </Table>
                </TableContainer>
              )}
            </>
          )}
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Input files" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="SEQIN"
            qualifiers={{ guiLabel: "File from which to extract sequence" }}
            onChange={handleSeqInChange}
          />

          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN"
            qualifiers={{ guiLabel: "Coordinate file (for extracting sequence or Matthews calc)" }}
            onChange={handleXyzInChange}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
