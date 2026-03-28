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
'use client';

import { useEffect, useState, useRef, useCallback } from 'react';
import { Box, Skeleton, Typography, Popper, Paper, Fade, IconButton, Tooltip, Chip } from '@mui/material';
import { ContentCopy, Check, Medication } from '@mui/icons-material';
import { useRDKit } from '@/lib/compounds/rdkit-context';

interface MoleculeViewProps {
  smiles: string;
  width?: number;
  height?: number;
  showSmiles?: boolean;
}

export const MoleculeView: React.FC<MoleculeViewProps> = ({
  smiles,
  width = 400,
  height = 400,
  showSmiles = false,
}) => {
  const { rdkitModule, isLoading: rdkitLoading } = useRDKit();
  const [svgUrl, setSvgUrl] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!rdkitModule || !smiles) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) {
        setError('Invalid SMILES');
        return;
      }

      const svg = mol.get_svg(width, height);
      mol.delete();

      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      setSvgUrl(url);
      setError(null);

      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error('RDKit error:', err);
      setError('Failed to render molecule');
    }
  }, [smiles, rdkitModule, width, height]);

  if (!smiles) {
    return (
      <Typography color="text.secondary" variant="body2">
        No SMILES
      </Typography>
    );
  }

  if (rdkitLoading) {
    return <Skeleton variant="rectangular" width={width} height={height} />;
  }

  if (error) {
    return (
      <Box sx={{ width, height, display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
        <Typography color="error" variant="body2">
          {error}
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', gap: 1 }}>
      {svgUrl && (
        <img
          src={svgUrl}
          alt={`Molecule: ${smiles}`}
          style={{ width, height, objectFit: 'contain' }}
        />
      )}
      {showSmiles && (
        <Typography
          variant="caption"
          sx={{
            fontFamily: 'monospace',
            wordBreak: 'break-all',
            maxWidth: width,
          }}
        >
          {smiles}
        </Typography>
      )}
    </Box>
  );
};

/**
 * Compact molecule view for table cells with hover-to-enlarge functionality
 */
interface MoleculeChipProps {
  smiles: string;
  size?: number;
  /** Disable hover overlay (useful when already showing large view) */
  disableHover?: boolean;
}

const HOVER_DELAY_MS = 400; // Delay before showing enlarged view
const ENLARGED_SIZE = 350; // Size of enlarged view

export const MoleculeChip: React.FC<MoleculeChipProps> = ({ smiles, size = 160, disableHover = false }) => {
  const { rdkitModule, isLoading } = useRDKit();
  const [svgUrl, setSvgUrl] = useState<string | null>(null);
  const [enlargedSvgUrl, setEnlargedSvgUrl] = useState<string | null>(null);
  const [showEnlarged, setShowEnlarged] = useState(false);
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const hoverTimeoutRef = useRef<NodeJS.Timeout | null>(null);

  // Generate thumbnail SVG
  useEffect(() => {
    if (!rdkitModule || !smiles) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) return;

      const svg = mol.get_svg(size, size);
      mol.delete();

      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      setSvgUrl(url);

      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error('RDKit error:', err);
    }
  }, [smiles, rdkitModule, size]);

  // Generate enlarged SVG (only when needed)
  useEffect(() => {
    if (!rdkitModule || !smiles || !showEnlarged) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) return;

      const svg = mol.get_svg(ENLARGED_SIZE, ENLARGED_SIZE);
      mol.delete();

      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      setEnlargedSvgUrl(url);

      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error('RDKit error:', err);
    }
  }, [smiles, rdkitModule, showEnlarged]);

  const handleMouseEnter = useCallback((event: React.MouseEvent<HTMLElement>) => {
    if (disableHover) return;
    setAnchorEl(event.currentTarget);
    hoverTimeoutRef.current = setTimeout(() => {
      setShowEnlarged(true);
    }, HOVER_DELAY_MS);
  }, [disableHover]);

  const handleMouseLeave = useCallback(() => {
    if (hoverTimeoutRef.current) {
      clearTimeout(hoverTimeoutRef.current);
      hoverTimeoutRef.current = null;
    }
    setShowEnlarged(false);
    setAnchorEl(null);
  }, []);

  // Cleanup timeout on unmount
  useEffect(() => {
    return () => {
      if (hoverTimeoutRef.current) {
        clearTimeout(hoverTimeoutRef.current);
      }
    };
  }, []);

  if (!smiles || isLoading) {
    return <Skeleton variant="rectangular" width={size} height={size} />;
  }

  if (!svgUrl) {
    return (
      <Typography variant="caption" color="text.secondary" sx={{ fontFamily: 'monospace' }}>
        {smiles.length > 20 ? smiles.slice(0, 20) + '...' : smiles}
      </Typography>
    );
  }

  return (
    <>
      <Box
        onMouseEnter={handleMouseEnter}
        onMouseLeave={handleMouseLeave}
        sx={{
          width: size,
          height: size,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          cursor: disableHover ? 'default' : 'zoom-in',
        }}
      >
        <img
          src={svgUrl}
          alt=""
          style={{ width: size, height: size, objectFit: 'contain' }}
        />
      </Box>
      {!disableHover && (
        <Popper
          open={showEnlarged}
          anchorEl={anchorEl}
          placement="right-start"
          transition
          sx={{ zIndex: 1300 }}
          modifiers={[
            {
              name: 'flip',
              enabled: true,
              options: {
                fallbackPlacements: ['left-start', 'bottom', 'top'],
              },
            },
            {
              name: 'preventOverflow',
              enabled: true,
              options: {
                boundary: 'viewport',
                padding: 8,
              },
            },
          ]}
        >
          {({ TransitionProps }) => (
            <Fade {...TransitionProps} timeout={200}>
              <Paper
                elevation={8}
                sx={{
                  p: 1,
                  bgcolor: 'background.paper',
                  border: 1,
                  borderColor: 'divider',
                }}
                onMouseEnter={() => setShowEnlarged(true)}
                onMouseLeave={handleMouseLeave}
              >
                {enlargedSvgUrl ? (
                  <img
                    src={enlargedSvgUrl}
                    alt={`Structure: ${smiles}`}
                    style={{ width: ENLARGED_SIZE, height: ENLARGED_SIZE, objectFit: 'contain' }}
                  />
                ) : (
                  <Skeleton variant="rectangular" width={ENLARGED_SIZE} height={ENLARGED_SIZE} />
                )}
              </Paper>
            </Fade>
          )}
        </Popper>
      )}
    </>
  );
};

/**
 * SMILES display with copy-to-clipboard functionality for table cells
 */
interface CopyableSmilesProps {
  smiles: string;
  maxWidth?: number;
}

/**
 * Compound name chip with hover popup showing full ID and molecular sketch.
 *
 * The default MUI Chip truncates labels with an ellipsis, which for compound
 * IDs like "NCL-00012345" shows only the common prefix ("NCL-00...").
 * This component adds a hover popup that always shows the full ID and,
 * when SMILES is available, renders a molecular structure sketch via RDKit.
 */
interface CompoundNameChipProps {
  formattedId: string;
  smiles?: string | null;
  /** MUI Chip color variant */
  chipColor?: 'default' | 'primary';
  /** Override the chip label text (hover popup still shows full formattedId) */
  label?: string;
}

const COMPOUND_HOVER_DELAY_MS = 300;
const COMPOUND_SKETCH_SIZE = 200;

export const CompoundNameChip: React.FC<CompoundNameChipProps> = ({
  formattedId,
  smiles,
  chipColor = 'default',
  label,
}) => {
  const { rdkitModule } = useRDKit();
  const [showPopup, setShowPopup] = useState(false);
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const [svgUrl, setSvgUrl] = useState<string | null>(null);
  const hoverTimeoutRef = useRef<NodeJS.Timeout | null>(null);
  const isInPopup = useRef(false);

  // Generate SVG lazily — only when popup is visible and SMILES is available
  useEffect(() => {
    if (!rdkitModule || !smiles || !showPopup) return;

    try {
      const mol = rdkitModule.get_mol(smiles);
      if (!mol) return;

      const svg = mol.get_svg(COMPOUND_SKETCH_SIZE, COMPOUND_SKETCH_SIZE);
      mol.delete();

      const blob = new Blob([svg], { type: 'image/svg+xml' });
      const url = URL.createObjectURL(blob);
      setSvgUrl(url);

      return () => URL.revokeObjectURL(url);
    } catch (err) {
      console.error('RDKit error:', err);
    }
  }, [smiles, rdkitModule, showPopup]);

  const handleMouseEnter = useCallback((event: React.MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget);
    hoverTimeoutRef.current = setTimeout(() => {
      setShowPopup(true);
    }, COMPOUND_HOVER_DELAY_MS);
  }, []);

  const hidePopup = useCallback(() => {
    // Small delay to allow mouse to travel from chip to popup
    setTimeout(() => {
      if (!isInPopup.current) {
        setShowPopup(false);
        setAnchorEl(null);
      }
    }, 100);
  }, []);

  const handleMouseLeave = useCallback(() => {
    if (hoverTimeoutRef.current) {
      clearTimeout(hoverTimeoutRef.current);
      hoverTimeoutRef.current = null;
    }
    hidePopup();
  }, [hidePopup]);

  const handlePopupEnter = useCallback(() => {
    isInPopup.current = true;
    setShowPopup(true);
  }, []);

  const handlePopupLeave = useCallback(() => {
    isInPopup.current = false;
    setShowPopup(false);
    setAnchorEl(null);
  }, []);

  // Cleanup on unmount
  useEffect(() => {
    return () => {
      if (hoverTimeoutRef.current) clearTimeout(hoverTimeoutRef.current);
    };
  }, []);

  return (
    <>
      <Chip
        icon={<Medication fontSize="small" />}
        label={label || formattedId}
        size="small"
        variant="outlined"
        color={chipColor}
        onMouseEnter={handleMouseEnter}
        onMouseLeave={handleMouseLeave}
        sx={{ '& .MuiChip-label': { fontFamily: 'monospace' } }}
      />
      <Popper
        open={showPopup}
        anchorEl={anchorEl}
        placement="bottom-start"
        transition
        sx={{ zIndex: 1300 }}
        modifiers={[
          {
            name: 'flip',
            enabled: true,
            options: { fallbackPlacements: ['top-start', 'right-start', 'left-start'] },
          },
          {
            name: 'preventOverflow',
            enabled: true,
            options: { boundary: 'viewport', padding: 8 },
          },
        ]}
      >
        {({ TransitionProps }) => (
          <Fade {...TransitionProps} timeout={200}>
            <Paper
              elevation={8}
              sx={{ p: 1.5, maxWidth: 280 }}
              onMouseEnter={handlePopupEnter}
              onMouseLeave={handlePopupLeave}
            >
              <Typography
                variant="subtitle2"
                sx={{ fontFamily: 'monospace', fontWeight: 600 }}
              >
                {formattedId}
              </Typography>
              {smiles && (
                <Box sx={{ mt: 1, display: 'flex', justifyContent: 'center' }}>
                  {svgUrl ? (
                    <img
                      src={svgUrl}
                      alt={`Structure: ${formattedId}`}
                      style={{
                        width: COMPOUND_SKETCH_SIZE,
                        height: COMPOUND_SKETCH_SIZE,
                        objectFit: 'contain',
                      }}
                    />
                  ) : (
                    <Skeleton
                      variant="rectangular"
                      width={COMPOUND_SKETCH_SIZE}
                      height={COMPOUND_SKETCH_SIZE}
                    />
                  )}
                </Box>
              )}
            </Paper>
          </Fade>
        )}
      </Popper>
    </>
  );
};

export const CopyableSmiles: React.FC<CopyableSmilesProps> = ({ smiles, maxWidth = 200 }) => {
  const [copied, setCopied] = useState(false);

  const handleCopy = useCallback(async (e: React.MouseEvent) => {
    e.stopPropagation(); // Prevent row click
    try {
      await navigator.clipboard.writeText(smiles);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    } catch (err) {
      console.error('Failed to copy:', err);
    }
  }, [smiles]);

  if (!smiles) {
    return <Typography variant="body2" color="text.secondary">-</Typography>;
  }

  return (
    <Box sx={{ display: 'flex', alignItems: 'center', gap: 0.5 }}>
      <Typography
        sx={{
          maxWidth,
          overflow: 'hidden',
          textOverflow: 'ellipsis',
          whiteSpace: 'nowrap',
          fontFamily: 'monospace',
          fontSize: '0.85rem',
        }}
        title={smiles}
      >
        {smiles}
      </Typography>
      <Tooltip title={copied ? 'Copied!' : 'Copy SMILES'}>
        <IconButton
          size="small"
          onClick={handleCopy}
          sx={{
            p: 0.5,
            color: copied ? 'success.main' : 'action.active',
            '&:hover': { bgcolor: 'action.hover' },
          }}
        >
          {copied ? <Check fontSize="small" /> : <ContentCopy fontSize="small" />}
        </IconButton>
      </Tooltip>
    </Box>
  );
};
