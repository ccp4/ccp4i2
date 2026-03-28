/*
 * Copyright (C) 2026 Newcastle University
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
"use client";
import {
  createContext,
  PropsWithChildren,
  useCallback,
  useContext,
  useEffect,
  useRef,
  useState,
} from "react";
import {
  IconButton,
  InputBase,
  Paper,
  Typography,
} from "@mui/material";
import {
  Close,
  KeyboardArrowUp,
  KeyboardArrowDown,
} from "@mui/icons-material";
import { isElectron } from "../utils/platform";

interface FindInPageContextValue {
  open: () => void;
}

const FindInPageContext = createContext<FindInPageContextValue | null>(null);

export function useFindInPage() {
  return useContext(FindInPageContext);
}

export function FindInPageProvider({ children }: PropsWithChildren) {
  const [visible, setVisible] = useState(false);
  const [query, setQuery] = useState("");
  const [activeMatch, setActiveMatch] = useState(0);
  const [totalMatches, setTotalMatches] = useState(0);
  const inputRef = useRef<HTMLInputElement>(null);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const electron = isElectron();

  const openFindBar = useCallback(() => {
    setVisible(true);
    // Focus input after render
    setTimeout(() => inputRef.current?.focus(), 50);
  }, []);

  const closeFindBar = useCallback(() => {
    setVisible(false);
    setQuery("");
    setActiveMatch(0);
    setTotalMatches(0);
    if (electron && window.electronAPI) {
      window.electronAPI.sendMessage("stop-find-in-page");
    }
  }, [electron]);

  const doFind = useCallback(
    (text: string, forward: boolean, findNext: boolean) => {
      if (!electron || !window.electronAPI || !text) return;
      window.electronAPI.sendMessage("find-in-page", {
        text,
        forward,
        findNext,
      });
    },
    [electron]
  );

  const findNext = useCallback(() => {
    if (query) doFind(query, true, true);
  }, [query, doFind]);

  const findPrev = useCallback(() => {
    if (query) doFind(query, false, true);
  }, [query, doFind]);

  // Listen for Cmd/Ctrl+F keydown
  useEffect(() => {
    const handleKeyDown = (e: KeyboardEvent) => {
      const mod = e.metaKey || e.ctrlKey;
      if (mod && e.key === "f") {
        e.preventDefault();
        openFindBar();
      }
      if (e.key === "Escape" && visible) {
        e.preventDefault();
        closeFindBar();
      }
    };
    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [visible, openFindBar, closeFindBar]);

  // Listen for toggle-find-in-page IPC from main process (menu accelerator)
  useEffect(() => {
    if (!electron || !window.electronAPI) return;

    const handler = (_event: any, data: any) => {
      if (data.message === "toggle-find-in-page") {
        if (visible) {
          closeFindBar();
        } else {
          openFindBar();
        }
      } else if (data.message === "found-in-page") {
        setActiveMatch(data.activeMatchOrdinal || 0);
        setTotalMatches(data.matches || 0);
      }
    };

    window.electronAPI.onMessage("message-from-main", handler);
    return () => {
      window.electronAPI.removeMessageListener("message-from-main", handler);
    };
  }, [electron, visible, openFindBar, closeFindBar]);

  // Debounced search on query change
  useEffect(() => {
    if (!visible) return;

    if (debounceRef.current) {
      clearTimeout(debounceRef.current);
    }

    if (!query) {
      setActiveMatch(0);
      setTotalMatches(0);
      if (electron && window.electronAPI) {
        window.electronAPI.sendMessage("stop-find-in-page");
      }
      return;
    }

    debounceRef.current = setTimeout(() => {
      doFind(query, true, false);
    }, 300);

    return () => {
      if (debounceRef.current) {
        clearTimeout(debounceRef.current);
      }
    };
  }, [query, visible, doFind, electron]);

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === "Enter") {
      e.preventDefault();
      if (e.shiftKey) {
        findPrev();
      } else {
        findNext();
      }
    }
    if (e.key === "Escape") {
      e.preventDefault();
      closeFindBar();
    }
  };

  const contextValue: FindInPageContextValue = { open: openFindBar };

  return (
    <>
      {visible && (
        <Paper
          elevation={4}
          sx={{
            position: "fixed",
            top: 8,
            right: 16,
            zIndex: 9999,
            display: "flex",
            alignItems: "center",
            gap: 0.5,
            px: 1,
            py: 0.5,
            borderRadius: 1,
          }}
        >
          <InputBase
            inputRef={inputRef}
            placeholder="Find in page..."
            value={query}
            onChange={(e) => setQuery(e.target.value)}
            onKeyDown={handleKeyDown}
            size="small"
            sx={{ fontSize: 14, width: 200 }}
            autoFocus
          />
          {query && (
            <Typography
              variant="caption"
              color="text.secondary"
              sx={{ whiteSpace: "nowrap", mx: 0.5 }}
            >
              {totalMatches > 0
                ? `${activeMatch}/${totalMatches}`
                : "No matches"}
            </Typography>
          )}
          <IconButton size="small" onClick={findPrev} disabled={!query}>
            <KeyboardArrowUp fontSize="small" />
          </IconButton>
          <IconButton size="small" onClick={findNext} disabled={!query}>
            <KeyboardArrowDown fontSize="small" />
          </IconButton>
          <IconButton size="small" onClick={closeFindBar}>
            <Close fontSize="small" />
          </IconButton>
        </Paper>
      )}
      <FindInPageContext.Provider value={contextValue}>
        {children}
      </FindInPageContext.Provider>
    </>
  );
}
