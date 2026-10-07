import React, { useEffect, useMemo, useCallback, useState } from "react";
import {
  alpha,
  Box,
  Typography,
  Chip,
  IconButton,
  Stack,
  ToggleButton,
  ToggleButtonGroup,
  useTheme as useMuiTheme,
} from "@mui/material";
import {
  ExpandLess,
  ExpandMore,
  CheckCircle,
  Code,
  ViewList,
} from "@mui/icons-material";
import { Editor } from "@monaco-editor/react";
import { useJob } from "../utils";
import { useApi } from "../api";
import { Job } from "../types/models";
import { useCCP4i2Window } from "../app-context";
import { useTheme } from "../theme/theme-provider";

interface ValidationViewerProps {
  job?: Job;
}

interface ValidationError {
  key: string;
  messages: string[];
  maxSeverity: number;
  displayName: string;
}

interface SeverityGroup {
  severity: number;
  label: string;
  color: string;
  icon: React.ReactNode;
  errors: ValidationError[];
}

type ViewMode = "formatted" | "xml";

// Severities as the validation parser in api.ts assigns them: 2 for ERROR,
// 1 for WARNING (and UNDEFINED_ERROR), 0 for anything else. This table used
// to read 2/3/1/0 as Errors/Warnings/Info/Debug, so warnings were filed under
// a collapsed "Info" group and the "Warnings" group could never appear.
// Colours are palette keys, resolved against the active theme (#613).
const SEVERITY_CONFIG = {
  2: { label: "Errors", palette: "error", icon: "⚠️" },
  1: { label: "Warnings", palette: "warning", icon: "⚡" },
  0: { label: "Other", palette: "info", icon: "ℹ️" },
} as const;
const SEVERITY_ORDER = [2, 1, 0] as const;

// Lightweight helper functions
const formatFieldName = (key: string): string => {
  const fieldName = key.split(".").pop() || key;
  // Handle underscores: F_SIGF -> F SIGF
  // Handle camelCase: myVariableName -> My Variable Name
  // But don't split all-caps sequences: XYZIN stays as XYZIN
  return fieldName
    .replace(/_/g, " ")
    // Only add space before uppercase if preceded by lowercase
    .replace(/([a-z])([A-Z])/g, "$1 $2")
    .replace(/\b\w/g, (l) => l.toUpperCase())
    .trim();
};

const cleanErrorMessage = (message: string): string => {
  return message
    .replace(/^[^:]+:\s*/, "")
    .replace(/Data has undefined value/, "Not defined");
};

export const ValidationViewer: React.FC<ValidationViewerProps> = ({ job }) => {
  const { customColors } = useTheme();
  const { mode } = useTheme();
  const muiTheme = useMuiTheme();
  const { validation } = useJob(job?.id);
  const { devMode } = useCCP4i2Window();
  const api = useApi();

  const [viewMode, setViewMode] = useState<ViewMode>("formatted");
  const [expandedGroups, setExpandedGroups] = useState<Set<number>>(
    new Set([2, 1]) // Expand errors and warnings by default
  );

  // Fetch raw XML for the XML view (with error handling)
  const {
    data: validationXml,
    error: xmlError,
    mutate: mutateValidationXml,
  } = api.get_pretty_endpoint_xml({
    type: "jobs",
    id: job?.id,
    endpoint: "validation",
  });

  // The XML is cached apart from the parsed validation (see get_validation),
  // so refresh it whenever the parsed copy is re-fetched.
  useEffect(() => {
    mutateValidationXml();
  }, [validation, mutateValidationXml]);

  const compiledErrors = useMemo(() => {
    return validation || {};
  }, [validation]);

  const processedData = useMemo((): SeverityGroup[] => {
    if (!compiledErrors || Object.keys(compiledErrors).length === 0) {
      return [];
    }

    // Quick transformation without heavy processing
    const errorsByseverity: Record<number, ValidationError[]> = {};

    Object.entries(compiledErrors).forEach(([key, errorData]) => {
      // Type assertion to ValidationError
      const err = errorData as ValidationError;
      const severity = err.maxSeverity;
      if (!errorsByseverity[severity]) {
        errorsByseverity[severity] = [];
      }

      errorsByseverity[severity].push({
        key,
        messages: err.messages.map(cleanErrorMessage),
        maxSeverity: severity,
        displayName: formatFieldName(key),
      });
    });

    // Build groups only for severities that have errors
    const groups: SeverityGroup[] = [];
    SEVERITY_ORDER.forEach((severity) => {
      const errors = errorsByseverity[severity];
      if (errors && errors.length > 0) {
        const config = SEVERITY_CONFIG[severity];
        groups.push({
          severity,
          label: config.label,
          color: muiTheme.palette[config.palette].main,
          icon: config.icon,
          errors: errors.sort((a, b) =>
            a.displayName.localeCompare(b.displayName)
          ),
        });
      }
    });

    return groups;
  }, [compiledErrors, muiTheme]);

  const totalErrors = processedData.reduce(
    (total, group) => total + group.errors.length,
    0
  );

  const toggleGroup = useCallback((severity: number) => {
    setExpandedGroups((prev) => {
      const newSet = new Set(prev);
      if (newSet.has(severity)) {
        newSet.delete(severity);
      } else {
        newSet.add(severity);
      }
      return newSet;
    });
  }, []);

  const handleViewModeChange = (
    _event: React.MouseEvent<HTMLElement>,
    newMode: ViewMode | null
  ) => {
    if (newMode !== null) {
      setViewMode(newMode);
    }
  };

  // Prettified JSON for dev mode
  const prettifiedValidation = useMemo(() => {
    return JSON.stringify(compiledErrors, null, 2);
  }, [compiledErrors]);

  // If in dev mode, show Monaco editor with JSON
  if (devMode) {
    return (
      <Box sx={{ height: "calc(100vh - 15rem)" }}>
        <Editor
          height="100%"
          defaultLanguage="json"
          value={prettifiedValidation}
          theme={mode === "dark" ? "vs-dark" : "light"}
          options={{
            readOnly: true,
            minimap: { enabled: false },
            scrollBeyondLastLine: false,
            wordWrap: "on",
            lineNumbers: "on",
            folding: true,
            automaticLayout: true,
          }}
        />
      </Box>
    );
  }

  // View toggle toolbar
  const ViewToggle = () => (
    <Box sx={{ display: "flex", justifyContent: "flex-end", mb: 1 }}>
      <ToggleButtonGroup
        value={viewMode}
        exclusive
        onChange={handleViewModeChange}
        size="small"
        aria-label="validation view mode"
      >
        <ToggleButton value="formatted" aria-label="formatted view">
          <ViewList sx={{ mr: 0.5 }} fontSize="small" />
          Formatted
        </ToggleButton>
        <ToggleButton value="xml" aria-label="XML view">
          <Code sx={{ mr: 0.5 }} fontSize="small" />
          XML
        </ToggleButton>
      </ToggleButtonGroup>
    </Box>
  );

  // XML view
  if (viewMode === "xml") {
    // Determine what to show in the editor
    let xmlContent = "<!-- Loading... -->";
    if (xmlError) {
      xmlContent = `<!-- Error loading validation XML: ${xmlError.message || "Unknown error"} -->`;
    } else if (validationXml) {
      xmlContent = validationXml;
    }

    return (
      <Box sx={{ p: 2 }}>
        <ViewToggle />
        <Box sx={{ height: "calc(100vh - 18rem)" }}>
          <Editor
            height="100%"
            defaultLanguage="xml"
            value={xmlContent}
            theme={mode === "dark" ? "vs-dark" : "light"}
            options={{
              readOnly: true,
              minimap: { enabled: false },
              scrollBeyondLastLine: false,
              wordWrap: "on",
              lineNumbers: "on",
              folding: true,
              automaticLayout: true,
            }}
          />
        </Box>
      </Box>
    );
  }

  // Regular beautified view for non-dev mode
  if (totalErrors === 0) {
    return (
      <Box sx={{ p: 2 }}>
        <ViewToggle />
        <Box
          sx={{
            p: 3,
            border: `1px solid ${customColors.ui.mediumGray}`,
            borderRadius: 1,
            backgroundColor: customColors.ui.veryLightGray,
            textAlign: "center",
          }}
        >
          <CheckCircle sx={{ color: "success.main", fontSize: 40, mb: 1 }} />
          <Typography variant="h6" sx={{ color: "success.main", mb: 0.5 }}>
            No Validation Issues
          </Typography>
          <Typography variant="body2" color="text.secondary">
            All required fields are properly configured.
          </Typography>
        </Box>
      </Box>
    );
  }

  return (
    <Box sx={{ p: 2 }}>
      <ViewToggle />
      <Box sx={{ maxHeight: "70vh", overflow: "auto" }}>
        {/* Lightweight Summary */}
        <Box
          sx={{
            p: 2,
            mb: 2,
            border: `1px solid ${customColors.ui.mediumGray}`,
            borderRadius: 1,
            backgroundColor: customColors.ui.veryLightGray,
          }}
        >
          <Stack
            direction="row"
            spacing={2}
            alignItems="center"
            flexWrap="wrap"
          >
            <Typography variant="h6">
              📋 {totalErrors} validation issues
            </Typography>
            {processedData.map((group) => (
              <Chip
                key={group.severity}
                size="small"
                label={`${group.errors.length} ${group.label.toLowerCase()}`}
                sx={{
                  backgroundColor: alpha(group.color, 0.12),
                  color: group.color,
                  border: `1px solid ${alpha(group.color, 0.25)}`,
                }}
              />
            ))}
          </Stack>
        </Box>

        {/* Lightweight Groups */}
        {processedData.map((group) => (
          <Box
            key={group.severity}
            sx={{
              mb: 2,
              border: `1px solid ${customColors.ui.mediumGray}`,
              borderRadius: 1,
              overflow: "hidden",
            }}
          >
            {/* Group Header */}
            <Box
              onClick={() => toggleGroup(group.severity)}
              sx={{
                p: 2,
                backgroundColor: alpha(group.color, 0.06),
                borderBottom: expandedGroups.has(group.severity)
                  ? `1px solid ${customColors.ui.mediumGray}`
                  : "none",
                cursor: "pointer",
                display: "flex",
                alignItems: "center",
                justifyContent: "space-between",
                "&:hover": {
                  backgroundColor: alpha(group.color, 0.12),
                },
              }}
            >
              <Stack direction="row" alignItems="center" spacing={1}>
                <span style={{ fontSize: "1.2em" }}>{group.icon}</span>
                <Typography variant="h6" sx={{ color: group.color }}>
                  {group.label}
                </Typography>
                <Chip
                  size="small"
                  label={group.errors.length}
                  sx={{
                    backgroundColor: group.color,
                    color: muiTheme.palette.getContrastText(group.color),
                    fontSize: "0.75rem",
                  }}
                />
              </Stack>
              <IconButton size="small" sx={{ color: group.color }}>
                {expandedGroups.has(group.severity) ? (
                  <ExpandLess />
                ) : (
                  <ExpandMore />
                )}
              </IconButton>
            </Box>

            {/* Group Content */}
            {expandedGroups.has(group.severity) && (
              <Box sx={{ p: 0 }}>
                {group.errors.map((error, index) => (
                  <Box
                    key={error.key}
                    sx={{
                      p: 2,
                      borderBottom:
                        index < group.errors.length - 1
                          ? `1px solid ${muiTheme.palette.divider}`
                          : "none",
                      "&:hover": {
                        backgroundColor: "action.hover",
                      },
                    }}
                  >
                    <Stack spacing={1}>
                      <Typography variant="subtitle2" sx={{ fontWeight: 600 }}>
                        {error.displayName}
                      </Typography>
                      {error.messages.map((message, msgIndex) => (
                        <Typography
                          key={msgIndex}
                          variant="body2"
                          sx={{ color: "text.secondary", pl: 1 }}
                        >
                          • {message}
                        </Typography>
                      ))}
                      <Typography
                        variant="caption"
                        sx={{
                          fontFamily: "monospace",
                          color: "text.disabled",
                          fontSize: "0.7rem",
                          pl: 1,
                        }}
                      >
                        {error.key}
                      </Typography>
                    </Stack>
                  </Box>
                ))}
              </Box>
            )}
          </Box>
        ))}
      </Box>
    </Box>
  );
};
