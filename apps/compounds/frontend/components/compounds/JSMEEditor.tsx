'use client';

import { useEffect, useRef, useState, useCallback } from 'react';
import { Box, CircularProgress, Typography } from '@mui/material';

// Declare the global JSApplet type
declare global {
  interface Window {
    JSApplet?: {
      JSME: new (
        containerId: string,
        options: { options: string }
      ) => JSMEApplet;
    };
  }
}

interface JSMEApplet {
  setCallBack: (eventName: string, callback: (jsmeEvent: JSMEEvent) => void) => void;
  options: (optionString: string) => void;
  readGenericMolecularInput: (genericObjectString: string) => void;
  smiles: () => string;
  molFile: () => string;
  clear: () => void;
  repaint: () => void;
}

interface JSMEEvent {
  src: {
    smiles: () => string;
    molFile: () => string;
  };
}

export interface JSMEEditorProps {
  /** Unique ID for the editor container */
  id?: string;
  /** Callback fired when the structure changes */
  onChange?: (smiles: string, molFile?: string) => void;
  /** Whether the editor is editable (true) or display-only (false) */
  editable?: boolean;
  /** Enable query/SMARTS mode for substructure searches */
  query?: boolean;
  /** Initial SMILES string to display */
  initialSmiles?: string | null;
  /** Show the SMILES string below the editor */
  showPreview?: boolean;
  /** Width of the editor */
  width?: string | number;
  /** Height of the editor */
  height?: string | number;
}

/**
 * JSME Molecule Editor component for drawing chemical structures.
 *
 * JSME (JavaScript Molecule Editor) is loaded dynamically from the public folder.
 * This component handles the loading of JSME scripts and provides a React-friendly
 * interface for molecule editing.
 */
export function JSMEEditor({
  id = 'jsme-editor',
  onChange,
  editable = true,
  query = false,
  initialSmiles = '',
  showPreview = false,
  width = 400,
  height = 350,
}: JSMEEditorProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const jsmeAppletRef = useRef<JSMEApplet | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [currentSmiles, setCurrentSmiles] = useState('');
  const initialSmilesLoadedRef = useRef(false);

  // Load JSME scripts
  useEffect(() => {
    const loadJSME = async () => {
      // Check if already loaded
      if (window.JSApplet) {
        setLoading(false);
        return;
      }

      try {
        // Load CSS
        if (!document.querySelector('link[href*="jsa.css"]')) {
          const link = document.createElement('link');
          link.rel = 'stylesheet';
          link.href = '/JSME_2020-06-11/jsme/jsa.css';
          document.head.appendChild(link);
        }

        // Load JSME script
        if (!document.querySelector('script[src*="jsme.nocache.js"]')) {
          await new Promise<void>((resolve, reject) => {
            const script = document.createElement('script');
            script.src = '/JSME_2020-06-11/jsme/jsme.nocache.js';
            script.async = true;
            script.onload = () => {
              // JSME sets a global JSApplet after loading
              // We need to wait a bit for it to initialize
              const checkJSApplet = setInterval(() => {
                if (window.JSApplet) {
                  clearInterval(checkJSApplet);
                  resolve();
                }
              }, 100);
              // Timeout after 10 seconds
              setTimeout(() => {
                clearInterval(checkJSApplet);
                if (!window.JSApplet) {
                  reject(new Error('JSME failed to initialize'));
                }
              }, 10000);
            };
            script.onerror = () => reject(new Error('Failed to load JSME script'));
            document.head.appendChild(script);
          });
        } else {
          // Script tag exists but JSApplet might still be loading
          await new Promise<void>((resolve, reject) => {
            const checkJSApplet = setInterval(() => {
              if (window.JSApplet) {
                clearInterval(checkJSApplet);
                resolve();
              }
            }, 100);
            setTimeout(() => {
              clearInterval(checkJSApplet);
              if (!window.JSApplet) {
                reject(new Error('JSME failed to initialize'));
              }
            }, 10000);
          });
        }

        setLoading(false);
      } catch (err) {
        setError(err instanceof Error ? err.message : 'Failed to load JSME');
        setLoading(false);
      }
    };

    loadJSME();
  }, []);

  // Handle structure changes
  const handleStructureChange = useCallback(
    (jsmeEvent: JSMEEvent) => {
      const jsme = jsmeEvent.src;
      const smiles = jsme.smiles();
      const molFile = jsme.molFile();
      setCurrentSmiles(smiles);
      onChange?.(smiles, molFile);
    },
    [onChange]
  );

  // Initialize JSME applet
  useEffect(() => {
    if (loading || error || !window.JSApplet || !containerRef.current) {
      return;
    }

    // Don't recreate if already initialized
    if (jsmeAppletRef.current) {
      return;
    }

    try {
      const options = ['newLook'];
      options.push(editable ? 'nodepict' : 'depict');
      if (query) {
        options.push('query');
      }

      jsmeAppletRef.current = new window.JSApplet.JSME(id, {
        options: options.join(','),
      });

      jsmeAppletRef.current.setCallBack('AfterStructureModified', handleStructureChange);

      // Load initial SMILES if provided
      if (initialSmiles && initialSmiles.trim().length > 0 && !initialSmilesLoadedRef.current) {
        jsmeAppletRef.current.readGenericMolecularInput(initialSmiles);
        initialSmilesLoadedRef.current = true;
      }
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to initialize JSME');
    }

    return () => {
      // Cleanup: clear the container
      if (containerRef.current) {
        containerRef.current.innerHTML = '';
      }
      jsmeAppletRef.current = null;
      initialSmilesLoadedRef.current = false;
    };
  }, [loading, error, id, editable, query, handleStructureChange, initialSmiles]);

  // Update editable mode
  useEffect(() => {
    if (jsmeAppletRef.current) {
      jsmeAppletRef.current.options(editable ? 'newlook,nodepict' : 'newlook,depict');
    }
  }, [editable]);

  // Update initial SMILES when it changes externally
  useEffect(() => {
    if (jsmeAppletRef.current && initialSmiles !== undefined) {
      if (initialSmiles && initialSmiles.trim().length > 0) {
        jsmeAppletRef.current.readGenericMolecularInput(initialSmiles);
      } else {
        jsmeAppletRef.current.clear();
      }
      jsmeAppletRef.current.repaint();
    }
  }, [initialSmiles]);

  if (error) {
    return (
      <Box
        sx={{
          width,
          height,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: 'error.light',
          borderRadius: 1,
          p: 2,
        }}
      >
        <Typography color="error">{error}</Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ position: 'relative' }}>
      {loading && (
        <Box
          sx={{
            width,
            height,
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'center',
            bgcolor: 'grey.100',
            borderRadius: 1,
          }}
        >
          <CircularProgress size={40} />
        </Box>
      )}
      <Box
        id={id}
        ref={containerRef}
        sx={{
          width,
          height,
          bgcolor: 'common.white',
          borderRadius: 1,
          overflow: 'hidden',
          display: loading ? 'none' : 'block',
        }}
      />
      {showPreview && currentSmiles && (
        <Typography
          variant="caption"
          component="pre"
          sx={{
            mt: 1,
            p: 1,
            bgcolor: 'grey.100',
            borderRadius: 1,
            fontFamily: 'monospace',
            fontSize: '0.75rem',
            overflow: 'auto',
            maxWidth: width,
          }}
        >
          {currentSmiles}
        </Typography>
      )}
    </Box>
  );
}
