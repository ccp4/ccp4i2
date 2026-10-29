"use client";
import { useEffect } from "react";
import { useRouter } from "next/navigation";
import {
  Alert,
  AlertTitle,
  Box,
  Button,
  Stack,
  Typography,
} from "@mui/material";
import { describePanelError } from "@/components/panel-error-boundary";

/**
 * The last resort for a render error nothing nearer caught (#614).
 *
 * Next.js shows "Application error: a client-side exception has occurred"
 * for an uncaught render error, and that page offers no way on but
 * restarting the app. This one keeps the user in the app: try the page
 * again, or go to the projects list. Errors inside a panel never reach here,
 * since each panel has a boundary of its own (PanelErrorBoundary).
 */
export default function CCP4i2Error({
  error,
  reset,
}: {
  error: Error & { digest?: string };
  reset: () => void;
}) {
  const router = useRouter();

  useEffect(() => {
    console.error("Unhandled render error", error);
  }, [error]);

  const copyDetails = async () => {
    const text = describePanelError(
      "page",
      error,
      null,
      navigator.userAgent,
      window.location.href,
    );
    try {
      await navigator.clipboard.writeText(text);
    } catch {
      /* no clipboard access: the message stays on screen to select by hand */
    }
  };

  return (
    <Box sx={{ p: 3, maxWidth: 720, mx: "auto" }}>
      <Alert severity="error">
        <AlertTitle>Something went wrong showing this page</AlertTitle>
        <Typography
          variant="body2"
          sx={{ fontFamily: "monospace", wordBreak: "break-word" }}
        >
          {error.message || String(error)}
        </Typography>
        <Typography variant="body2" sx={{ mt: 1 }}>
          Your projects and jobs are not affected. Try the page again, or go
          back to the projects list. If it fails the same way, copy the details
          into a bug report.
        </Typography>
        <Stack direction="row" spacing={1} sx={{ mt: 2 }}>
          <Button variant="contained" size="small" onClick={reset}>
            Try again
          </Button>
          <Button
            variant="outlined"
            size="small"
            onClick={() => router.push("/ccp4i2")}
          >
            Projects list
          </Button>
          <Button size="small" onClick={copyDetails}>
            Copy details
          </Button>
        </Stack>
      </Alert>
    </Box>
  );
}
