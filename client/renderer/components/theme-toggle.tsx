"use client";
import React from "react";
import { IconButton, Tooltip } from "@mui/material";
import { Brightness4, Brightness7 } from "@mui/icons-material";
import { useTheme } from "../theme/theme-provider";

export const ThemeToggle = React.forwardRef<HTMLButtonElement>((props, ref) => {
  const { mode, toggleTheme } = useTheme();

  const handleToggleTheme = () => {
    const newMode = mode === "light" ? "dark" : "light";

    // Send message to Electron API if available
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("set-theme-mode", { theme: newMode });
    }

    // Toggle the theme
    toggleTheme();
  };

  return (
    <Tooltip title={`Switch to ${mode === "light" ? "dark" : "light"} mode`}>
      <IconButton
        ref={ref}
        onClick={handleToggleTheme}
        color="inherit"
        size="small"
        sx={{
          ml: 1,
          "&:hover": {
            backgroundColor: "rgba(255, 255, 255, 0.1)",
          },
        }}
      >
        {mode === "light" ? <Brightness4 /> : <Brightness7 />}
      </IconButton>
    </Tooltip>
  );
});
