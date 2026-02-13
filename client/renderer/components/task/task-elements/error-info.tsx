import React, {
  memo,
  SyntheticEvent,
  useCallback,
  useEffect,
  useMemo,
  useState,
} from "react";
import {
  Box,
  Button,
  Card,
  CardContent,
  CardHeader,
  ClickAwayListener,
  Collapse,
  LinearProgress,
  Popper,
  Stack,
  Typography,
} from "@mui/material";
import {
  ExpandLess,
  ExpandMore,
  Info,
  Warning,
  CheckCircle,
  Error as ErrorIcon,
} from "@mui/icons-material";

import { CCP4i2TaskElementProps } from "./task-element";
import { useJob, ValidationError, valueOfItem } from "../../../utils";
import { useTaskInterface } from "../../../providers/task-provider";
import { Job } from "../../../types/models";
import { SimpleObjectTable } from "../../simple-object-table";
import { useCCP4i2Window } from "../../../app-context";

// ===== TYPES =====
interface ErrorTriggerProps {
  item: any;
  job: Job;
}

interface ErrorPopperProps {}

interface ValidationMessage {
  path: string;
  messages: string[];
  severity: number;
}

interface ProcessedErrorInfo {
  hasErrors: boolean;
  hasWarnings: boolean;
  messages: ValidationMessage[];
  backgroundColor: string;
  icon: React.ReactNode;
}

interface CollapsibleState {
  isOpen: boolean;
  toggle: () => void;
  close: () => void;
  open: () => void;
}

// ===== CONSTANTS =====
const VALIDATION_COLORS = {
  SUCCESS: "success.light",
  WARNING: "warning.light",
  ERROR: "error.light",
} as const;

const SEVERITY_LEVELS = {
  INFO: 0,
  WARNING: 1,
  ERROR: 2,
} as const;

const POPPER_STYLES = {
  zIndex: 1000,
  maxWidth: "min(40rem, 90vw)",
  maxHeight: "min(calc(30 * 1.5em), 80vh)", // 30 lines * 1.5em line-height, or 80vh max
  bgcolor: "background.paper",
  border: 1,
  borderColor: "divider",
  borderRadius: 1,
  boxShadow: 3,
  overflow: "hidden", // Ensure the popper itself doesn't overflow
} as const;

const ERROR_CONTAINER_STYLES = {
  p: 2,
  borderRadius: 1,
  mb: 2,
} as const;

const POPPER_MODIFIERS = [
  {
    name: "preventOverflow",
    options: {
      boundary: "viewport",
      padding: 8,
    },
  },
  {
    name: "flip",
    options: {
      fallbackPlacements: ["auto", "auto-start", "auto-end", "top", "bottom"],
    },
  },
  {
    name: "computeStyles",
    options: {
      adaptive: true,
    },
  },
];

// ===== UTILITY FUNCTIONS =====
const formatObjectPath = (objectPath?: string): string => {
  return objectPath?.split(".").pop() || objectPath || "Unknown item";
};

const getValidationSummary = (messages: ValidationMessage[]): string => {
  const errorCount = messages.filter(
    (msg) => msg.severity >= SEVERITY_LEVELS.ERROR
  ).length;
  const warningCount = messages.filter(
    (msg) => msg.severity === SEVERITY_LEVELS.WARNING
  ).length;

  if (errorCount > 0 && warningCount > 0) {
    return `${errorCount} error(s) and ${warningCount} warning(s)`;
  }
  if (errorCount > 0) {
    return `${errorCount} error(s)`;
  }
  if (warningCount > 0) {
    return `${warningCount} warning(s)`;
  }
  return "No issues";
};

const getValidationIcon = (validationColor: string) => {
  switch (validationColor) {
    case VALIDATION_COLORS.SUCCESS:
      return <CheckCircle fontSize="small" color="success" />;
    case VALIDATION_COLORS.WARNING:
      return <Warning fontSize="small" color="warning" />;
    case VALIDATION_COLORS.ERROR:
      return <ErrorIcon fontSize="small" color="error" />;
    default:
      return <Info fontSize="small" />;
  }
};

const processItemValue = (item: any) => {
  if (!item) return null;

  const result = valueOfItem(item);

  if (result === null || result === undefined) {
    return null;
  }

  if (typeof result === "boolean") {
    return { value: result ? "true" : "false" };
  }

  if (typeof result === "object") {
    return result;
  }

  return { value: String(result) };
};

// ===== CUSTOM HOOKS =====
const useValidationIcon = (validationColor: string) => {
  return useMemo(() => getValidationIcon(validationColor), [validationColor]);
};

const useProcessedErrorInfo = (
  item: any,
  fieldErrors: ValidationError[],
  validationColor: string
): ProcessedErrorInfo => {
  return useMemo(() => {
    if (!item || !fieldErrors?.length) {
      return {
        hasErrors: false,
        hasWarnings: false,
        messages: [],
        backgroundColor: VALIDATION_COLORS.SUCCESS,
        icon: getValidationIcon(VALIDATION_COLORS.SUCCESS),
      };
    }

    const messages: ValidationMessage[] = fieldErrors.map(
      (validationError) => ({
        path: validationError.path,
        messages: validationError.error?.messages || [],
        severity: validationError.error?.maxSeverity || 0,
      })
    );

    const hasErrors = messages.some(
      (msg) => msg.severity >= SEVERITY_LEVELS.ERROR
    );
    const hasWarnings = messages.some(
      (msg) => msg.severity === SEVERITY_LEVELS.WARNING
    );

    let backgroundColor: string = VALIDATION_COLORS.SUCCESS;
    if (hasErrors) {
      backgroundColor = VALIDATION_COLORS.ERROR;
    } else if (hasWarnings) {
      backgroundColor = VALIDATION_COLORS.WARNING;
    }

    return {
      hasErrors,
      hasWarnings,
      messages,
      backgroundColor,
      icon: getValidationIcon(backgroundColor),
    };
  }, [item, fieldErrors, validationColor]);
};

const useItemValue = (item: any) => {
  return useMemo(() => processItemValue(item), [item]);
};

const useCollapsibleState = (initialState = false): CollapsibleState => {
  const [isOpen, setIsOpen] = useState(initialState);

  const toggle = useCallback(() => setIsOpen((prev) => !prev), []);
  const close = useCallback(() => setIsOpen(false), []);
  const open = useCallback(() => setIsOpen(true), []);

  return { isOpen, toggle, close, open };
};

// ===== COMPONENT IMPLEMENTATIONS =====

// Memoized validation messages component
const ValidationMessages = memo<{
  processedErrorInfo: ProcessedErrorInfo;
  objectPath?: string;
}>(({ processedErrorInfo, objectPath }) => {
  if (!processedErrorInfo.hasErrors && !processedErrorInfo.hasWarnings) {
    return (
      <Box
        sx={{
          ...ERROR_CONTAINER_STYLES,
          bgcolor: processedErrorInfo.backgroundColor,
        }}
        title={objectPath}
      >
        <Stack direction="row" alignItems="center" spacing={1}>
          {processedErrorInfo.icon}
          <Typography variant="subtitle2">
            No validation issues for {formatObjectPath(objectPath)}
          </Typography>
        </Stack>
      </Box>
    );
  }

  return (
    <Box
      sx={{
        ...ERROR_CONTAINER_STYLES,
        bgcolor: processedErrorInfo.backgroundColor,
      }}
      title={objectPath}
    >
      <Stack spacing={1}>
        <Stack direction="row" alignItems="center" spacing={1}>
          {processedErrorInfo.icon}
          <Typography variant="subtitle2" fontWeight="medium">
            {getValidationSummary(processedErrorInfo.messages)} in{" "}
            {formatObjectPath(objectPath)}
          </Typography>
        </Stack>

        {processedErrorInfo.messages.map((validationMessage, index) =>
          validationMessage.messages.map((message, messageIndex) => (
            <Typography
              key={`${validationMessage.path}_${index}_${messageIndex}`}
              variant="body2"
              sx={{
                ml: 3,
                wordWrap: "break-word",
                maxWidth: "100%",
              }}
            >
              • {message}
            </Typography>
          ))
        )}
      </Stack>
    </Box>
  );
});

ValidationMessages.displayName = "ValidationMessages";

// Memoized collapsible card component
const CollapsibleCard = memo<{
  title: string;
  children: React.ReactNode;
  collapsibleState: CollapsibleState;
  variant?: "outlined" | "elevation";
}>(({ title, children, collapsibleState, variant = "outlined" }) => (
  <Card variant={variant} sx={{ mb: 1 }}>
    <CardHeader
      avatar={collapsibleState.isOpen ? <ExpandLess /> : <ExpandMore />}
      title={title}
      titleTypographyProps={{ variant: "subtitle2" }}
      onClick={collapsibleState.toggle}
      sx={{
        cursor: "pointer",
        "&:hover": {
          backgroundColor: "action.hover",
        },
      }}
    />
    <Collapse in={collapsibleState.isOpen} timeout="auto" unmountOnExit>
      <CardContent sx={{ pt: 0 }}>{children}</CardContent>
    </Collapse>
  </Card>
));

CollapsibleCard.displayName = "CollapsibleCard";

// Memoized simple card component
const SimpleCard = memo<{
  title: string;
  children: React.ReactNode;
  variant?: "outlined" | "elevation";
}>(({ title, children, variant = "outlined" }) => (
  <Card variant={variant}>
    <CardHeader title={title} titleTypographyProps={{ variant: "subtitle2" }} />
    <CardContent sx={{ pt: 0 }}>{children}</CardContent>
  </Card>
));

SimpleCard.displayName = "SimpleCard";

// Main error trigger component
// Note: Not using memo() here because this component depends on processedErrors
// from context (via useJob -> useRunCheck), which wouldn't trigger re-renders
// if we memoized based only on props.
export const ErrorTrigger: React.FC<ErrorTriggerProps> = ({ item, job }) => {
  const { setErrorInfoAnchor, setErrorInfoItem } = useTaskInterface();
  const { getValidationColor } = useJob(job.id);
  const { devMode } = useCCP4i2Window();

  const validationColor = useMemo(
    () => getValidationColor(item),
    [getValidationColor, item]
  );

  const hasValidationIssue =
    validationColor === VALIDATION_COLORS.ERROR ||
    validationColor === VALIDATION_COLORS.WARNING;

  const icon = useValidationIcon(validationColor);

  const handleClick = useCallback(
    (event: SyntheticEvent) => {
      event.stopPropagation();
      event.preventDefault();
      setErrorInfoAnchor(event.currentTarget);
      setErrorInfoItem(item);
    },
    [setErrorInfoAnchor, setErrorInfoItem, item]
  );

  // Only show when in developer mode or when the field has validation issues
  if (!devMode && !hasValidationIssue) {
    return null;
  }

  return (
    <Button
      size="small"
      variant="text"
      onClick={handleClick}
      sx={{
        minWidth: "auto",
        p: 0.5,
        color: validationColor,
        "&:hover": {
          backgroundColor: `${validationColor}20`,
        },
      }}
      aria-label={`Show validation information for ${formatObjectPath(
        item?._objectPath
      )}`}
    >
      {icon}
    </Button>
  );
};

ErrorTrigger.displayName = "ErrorTrigger";

// Main error popper component
export const ErrorPopper: React.FC<ErrorPopperProps> = memo(() => {
  const { jobId } = useCCP4i2Window();
  const { job } = useJob(jobId);
  const {
    setErrorInfoAnchor,
    errorInfoAnchor,
    setErrorInfoItem,
    errorInfoItem,
  } = useTaskInterface();

  const { getValidationColor, getErrors } = useJob(job?.id || 0);

  const fieldErrors = useMemo(
    () => (errorInfoItem ? getErrors(errorInfoItem) : null),
    [errorInfoItem, getErrors]
  );

  const validationColor = useMemo(
    () => (errorInfoItem ? getValidationColor(errorInfoItem) : ""),
    [errorInfoItem, getValidationColor]
  );

  const qualifiersState = useCollapsibleState(false);
  const processedErrorInfo = useProcessedErrorInfo(
    errorInfoItem,
    fieldErrors ?? [],
    validationColor
  );
  const itemValue = useItemValue(errorInfoItem);

  const handleClickAway = useCallback(
    (event: MouseEvent | TouchEvent) => {
      event.stopPropagation();
      setErrorInfoAnchor(null);
      setErrorInfoItem(null);
    },
    [setErrorInfoAnchor, setErrorInfoItem]
  );

  // Debug logging in development only
  useEffect(() => {
    if (process.env.NODE_ENV === "development") {
      console.log("ErrorPopper state:", {
        errorInfoItem: errorInfoItem?._objectPath,
        hasAnchor: Boolean(errorInfoAnchor),
        fieldErrors,
      });
    }
  }, [errorInfoItem, errorInfoAnchor, fieldErrors]);

  if (!errorInfoItem || !errorInfoAnchor) {
    return null;
  }

  return (
    <Popper
      anchorEl={errorInfoAnchor}
      placement="auto-end"
      open={Boolean(errorInfoAnchor)}
      sx={POPPER_STYLES}
      modifiers={POPPER_MODIFIERS}
    >
      <ClickAwayListener onClickAway={handleClickAway}>
        <Box
          sx={{
            p: 2,
            maxHeight: "calc(30 * 1.5em + 2rem)", // 30 lines + padding
            overflow: "auto",
            lineHeight: 1.5, // Ensure consistent line height
            fontSize: "0.875rem", // Match MUI body2 font size for consistency
            "&::-webkit-scrollbar": {
              width: "8px",
            },
            "&::-webkit-scrollbar-track": {
              backgroundColor: "rgba(0,0,0,0.05)",
              borderRadius: "4px",
            },
            "&::-webkit-scrollbar-thumb": {
              backgroundColor: "rgba(0,0,0,0.2)",
              borderRadius: "4px",
              "&:hover": {
                backgroundColor: "rgba(0,0,0,0.3)",
              },
            },
            // Firefox scrollbar styling
            scrollbarWidth: "thin",
            scrollbarColor: "rgba(0,0,0,0.2) rgba(0,0,0,0.05)",
          }}
        >
          <Stack spacing={2}>
            <ValidationMessages
              processedErrorInfo={processedErrorInfo}
              objectPath={errorInfoItem?._objectPath}
            />

            {errorInfoItem?._qualifiers && (
              <CollapsibleCard
                title="Item Qualifiers"
                collapsibleState={qualifiersState}
              >
                <SimpleObjectTable object={errorInfoItem._qualifiers} />
              </CollapsibleCard>
            )}

            {itemValue && (
              <SimpleCard title="Item Value">
                <SimpleObjectTable object={itemValue} />
              </SimpleCard>
            )}
          </Stack>
        </Box>
      </ClickAwayListener>
    </Popper>
  );
});

ErrorPopper.displayName = "ErrorPopper";

// Simple loading component
export const ErrorInfo: React.FC<CCP4i2TaskElementProps> = memo(() => (
  <LinearProgress
    variant="indeterminate"
    sx={{
      height: 2,
      borderRadius: 1,
    }}
  />
));

ErrorInfo.displayName = "ErrorInfo";
