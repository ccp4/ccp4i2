import React, { useState, useMemo } from "react";
import {
  Card,
  CardHeader,
  CardContent,
  Collapse,
  Typography,
  Box,
  IconButton,
  Accordion,
  AccordionSummary,
  AccordionDetails,
  Alert,
  Chip,
  Button,
} from "@mui/material";
import {
  ExpandMore as ExpandMoreIcon,
  Error as ErrorIcon,
  Warning as WarningIcon,
  Info as InfoIcon,
  KeyboardArrowDown as KeyboardArrowDownIcon,
  KeyboardArrowRight as KeyboardArrowRightIcon,
} from "@mui/icons-material";

interface DiagnosticProps {
  xmlDocument: string;
}

interface ErrorReport {
  className: string;
  code: string;
  description: string;
  severity: string;
  details: string;
  stack: string;
}

interface HeaderInfo {
  function: string;
  userId: string;
  hostName: string;
  creationTime: string;
  ccp4iVersion: string;
  pluginName: string;
  projectName: string;
  projectId: string;
  jobId: string;
  jobNumber: string;
}

const Diagnostic: React.FC<DiagnosticProps> = ({ xmlDocument }) => {
  const [expandedCard, setExpandedCard] = useState<number>(-1); // Initially nothing expanded
  const [expandedStacks, setExpandedStacks] = useState<Set<number>>(new Set()); // Track which stack traces are open

  const parseXML = (xmlString: string) => {
    const parser = new DOMParser();
    const doc = parser.parseFromString(xmlString, "text/xml");

    // Parse header information
    const header: Partial<HeaderInfo> = {};
    const headerElement = doc.querySelector("ccp4i2_header");
    if (headerElement) {
      [
        "function",
        "userId",
        "hostName",
        "creationTime",
        "ccp4iVersion",
        "pluginName",
        "projectName",
        "projectId",
        "jobId",
        "jobNumber",
      ].forEach((field) => {
        const element = headerElement.querySelector(field);
        if (element) {
          header[field as keyof HeaderInfo] = element.textContent || "";
        }
      });
    }

    // Parse error reports
    const errorReports: ErrorReport[] = [];
    const errorElements = doc.querySelectorAll("errorReport");

    errorElements.forEach((errorElement) => {
      const getTextContent = (selector: string) =>
        errorElement.querySelector(selector)?.textContent?.trim() || "";

      errorReports.push({
        className: getTextContent("className"),
        code: getTextContent("code"),
        description: getTextContent("description"),
        severity: getTextContent("severityName") || getTextContent("severity"),
        details: getTextContent("details"),
        stack: getTextContent("stack"),
      });
    });

    return { header, errorReports };
  };

  const { header, errorReports } = useMemo(
    () => parseXML(xmlDocument),
    [xmlDocument]
  );

  // Auto-expand the first ERROR on mount
  React.useEffect(() => {
    const firstErrorIndex = errorReports.findIndex(
      (report) => report.severity.toUpperCase() === "ERROR"
    );
    if (firstErrorIndex !== -1) {
      setExpandedCard(firstErrorIndex);
    }
  }, [errorReports]);

  const getSeverityIcon = (severity: string) => {
    switch (severity.toUpperCase()) {
      case "ERROR":
        return <ErrorIcon color="error" />;
      case "WARNING":
        return <WarningIcon color="warning" />;
      default:
        return <InfoIcon color="info" />;
    }
  };

  const getSeverityColor = (severity: string): "error" | "warning" | "info" => {
    switch (severity.toUpperCase()) {
      case "ERROR":
        return "error";
      case "WARNING":
        return "warning";
      default:
        return "info";
    }
  };

  const handleCardToggle = (index: number) => {
    setExpandedCard(expandedCard === index ? -1 : index);
  };

  const handleStackToggle = (index: number) => {
    setExpandedStacks((prev) => {
      const newSet = new Set(prev);
      if (newSet.has(index)) {
        newSet.delete(index);
      } else {
        newSet.add(index);
      }
      return newSet;
    });
  };

  if (!errorReports.length) {
    return (
      <Box
        sx={{
          height: "calc(100vh - 19rem)",
          overflowY: "auto",
          p: 2,
        }}
      >
        <Alert severity="info">
          No error reports found in the diagnostic data.
        </Alert>
      </Box>
    );
  }

  return (
    <Box
      sx={{
        height: "calc(100vh - 21rem)",
        overflowY: "auto",
        width: "100%",
        p: 2,
      }}
    >
      {/* Header Information */}
      {header.pluginName && (
        <Card sx={{ mb: 2 }}>
          <CardHeader
            title="Job Information"
            slotProps={{ title: { typography: { variant: "h6" } } }}
          />
          <CardContent>
            <Box
              sx={{
                display: "grid",
                gridTemplateColumns: "repeat(auto-fit, minmax(200px, 1fr))",
                gap: 1,
              }}
            >
              {header.pluginName && (
                <Typography variant="body2">
                  <strong>Plugin:</strong> {header.pluginName}
                </Typography>
              )}
              {header.projectName && (
                <Typography variant="body2">
                  <strong>Project:</strong> {header.projectName}
                </Typography>
              )}
              {header.jobNumber && (
                <Typography variant="body2">
                  <strong>Job #:</strong> {header.jobNumber}
                </Typography>
              )}
              {header.creationTime && (
                <Typography variant="body2">
                  <strong>Created:</strong> {header.creationTime}
                </Typography>
              )}
              {header.userId && (
                <Typography variant="body2">
                  <strong>User:</strong> {header.userId}
                </Typography>
              )}
              {header.hostName && (
                <Typography variant="body2">
                  <strong>Host:</strong> {header.hostName}
                </Typography>
              )}
            </Box>
          </CardContent>
        </Card>
      )}

      {/* Error Reports */}
      <Typography variant="h6" gutterBottom>
        Error Reports ({errorReports.length})
      </Typography>

      {errorReports.map((report, index) => {
        const isError = report.severity.toUpperCase() === "ERROR";
        return (
          <Accordion
            key={index}
            expanded={expandedCard === index}
            onChange={() => handleCardToggle(index)}
            sx={{
              mb: 1,
              ...(isError && {
                border: "1px solid",
                borderColor: "error.main",
                borderLeft: "4px solid",
                borderLeftColor: "error.main",
              }),
            }}
          >
            <AccordionSummary
              expandIcon={<ExpandMoreIcon />}
              sx={{
                backgroundColor: isError
                  ? expandedCard === index
                    ? "error.light"
                    : "error.lighter"
                  : expandedCard === index
                  ? "action.selected"
                  : "background.paper",
                "&:hover": {
                  backgroundColor: isError ? "error.light" : "action.hover",
                },
              }}
            >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                width: "100%",
              }}
            >
              {getSeverityIcon(report.severity)}
              <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
                {report.className} - {report.code}
              </Typography>
              <Chip
                label={report.severity}
                color={getSeverityColor(report.severity)}
                size="small"
              />
            </Box>
          </AccordionSummary>
          <AccordionDetails>
            <Box sx={{ display: "flex", flexDirection: "column", gap: 2 }}>
              {/* Description */}
              <Box>
                <Typography
                  variant="subtitle2"
                  color="text.secondary"
                  gutterBottom
                >
                  Description
                </Typography>
                <Typography variant="body2">{report.description}</Typography>
              </Box>

              {/* Details */}
              {report.details && (
                <Box>
                  <Typography
                    variant="subtitle2"
                    color="text.secondary"
                    gutterBottom
                  >
                    Details
                  </Typography>
                  <Box
                    component="pre"
                    sx={{
                      backgroundColor: "grey.100",
                      p: 2,
                      borderRadius: 1,
                      fontSize: "0.875rem",
                      fontFamily: "monospace",
                      whiteSpace: "pre-wrap",
                      wordBreak: "break-word",
                      maxHeight: "200px",
                      overflow: "auto",
                    }}
                  >
                    {report.details}
                  </Box>
                </Box>
              )}

              {/* Stack Trace - Now Collapsible */}
              {report.stack && (
                <Box>
                  <Button
                    onClick={() => handleStackToggle(index)}
                    startIcon={
                      expandedStacks.has(index) ? (
                        <KeyboardArrowDownIcon />
                      ) : (
                        <KeyboardArrowRightIcon />
                      )
                    }
                    sx={{
                      justifyContent: "flex-start",
                      p: 0,
                      minWidth: "auto",
                      textTransform: "none",
                      color: "text.secondary",
                      "&:hover": {
                        backgroundColor: "transparent",
                        color: "primary.main",
                      },
                    }}
                  >
                    <Typography variant="subtitle2">Stack Trace</Typography>
                  </Button>
                  <Collapse in={expandedStacks.has(index)}>
                    <Box
                      component="pre"
                      sx={{
                        backgroundColor: "grey.100",
                        p: 2,
                        borderRadius: 1,
                        fontSize: "0.75rem",
                        fontFamily: "monospace",
                        whiteSpace: "pre-wrap",
                        wordBreak: "break-word",
                        maxHeight: "300px",
                        overflow: "auto",
                        mt: 1,
                      }}
                    >
                      {report.stack}
                    </Box>
                  </Collapse>
                </Box>
              )}
            </Box>
          </AccordionDetails>
        </Accordion>
        );
      })}
    </Box>
  );
};

export default Diagnostic;
