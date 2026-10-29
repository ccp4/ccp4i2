"use client";
import React, { ErrorInfo, PropsWithChildren } from "react";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Stack,
  Typography,
} from "@mui/material";

/**
 * Keeps a render error inside the panel it happened in (#614).
 *
 * React unmounts everything up to the nearest error boundary when a render
 * throws, and with no boundary at all that was the whole page: Next.js's
 * "Application error: a client-side exception has occurred" with nothing to
 * do but restart the app. A boundary around each panel (job list, project
 * directory, the job's tabs) confines the damage to that panel: the rest of
 * the window keeps working, and the panel offers to try again, or clears
 * itself when the user moves to another job or tab (`resetKeys`).
 *
 * Only render errors are caught. Errors in event handlers and async code do
 * not reach a boundary; those keep going to the snackbar and the console.
 */
interface PanelErrorBoundaryProps extends PropsWithChildren {
  /** What the panel shows, for the message: "the job list". */
  name: string;
  /** Values whose change clears the error: a job id, a tab index. */
  resetKeys?: readonly unknown[];
}

interface PanelErrorBoundaryState {
  error: Error | null;
  componentStack: string | null;
}

function sameKeys(a: readonly unknown[] = [], b: readonly unknown[] = []) {
  return a.length === b.length && a.every((value, i) => Object.is(value, b[i]));
}

/** The error as plain text for a bug report: where, what, and the stacks. */
export function describePanelError(
  name: string,
  error: Error,
  componentStack: string | null,
  userAgent: string,
  url: string,
): string {
  const lines = [
    `CCP4i2 panel error: ${name}`,
    `Client: ${userAgent}`,
    `Page: ${url}`,
    "",
    `${error.name}: ${error.message}`,
  ];
  if (error.stack) lines.push("", error.stack);
  if (componentStack) lines.push("", "Component stack:", componentStack.trim());
  return lines.join("\n") + "\n";
}

export class PanelErrorBoundary extends React.Component<
  PanelErrorBoundaryProps,
  PanelErrorBoundaryState
> {
  state: PanelErrorBoundaryState = { error: null, componentStack: null };

  static getDerivedStateFromError(
    error: Error,
  ): Partial<PanelErrorBoundaryState> {
    return { error };
  }

  componentDidCatch(error: Error, info: ErrorInfo) {
    this.setState({ componentStack: info.componentStack ?? null });
    console.error(
      `Error in the ${this.props.name} panel`,
      error,
      info.componentStack,
    );
  }

  componentDidUpdate(prevProps: PanelErrorBoundaryProps) {
    if (
      this.state.error &&
      !sameKeys(prevProps.resetKeys, this.props.resetKeys)
    ) {
      this.reset();
    }
  }

  reset = () => this.setState({ error: null, componentStack: null });

  copyDetails = async () => {
    const { error, componentStack } = this.state;
    if (!error) return;
    const text = describePanelError(
      this.props.name,
      error,
      componentStack,
      navigator.userAgent,
      window.location.href,
    );
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      /* no clipboard access: the details stay on screen to select by hand */
    }
  };

  render() {
    const { error } = this.state;
    if (!error) return this.props.children;
    return (
      <Box sx={{ p: 2, overflow: "auto" }}>
        <Alert
          severity="error"
          action={
            <Stack direction="row" spacing={1}>
              <Button color="inherit" size="small" onClick={this.copyDetails}>
                Copy details
              </Button>
              <Button
                color="inherit"
                size="small"
                variant="outlined"
                onClick={this.reset}
              >
                Try again
              </Button>
            </Stack>
          }
        >
          <AlertTitle>
            Something went wrong showing the {this.props.name}
          </AlertTitle>
          <Typography
            variant="body2"
            sx={{ fontFamily: "monospace", wordBreak: "break-word" }}
          >
            {error.message || String(error)}
          </Typography>
          <Typography variant="body2" sx={{ mt: 1 }}>
            The rest of the window still works. Try again redraws this panel; if
            it fails the same way, copy the details into a bug report.
          </Typography>
        </Alert>
      </Box>
    );
  }
}
