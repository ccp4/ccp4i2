'use client';

import { useEffect, useState } from 'react';
import { Box, Skeleton, Typography } from '@mui/material';
import { useRDKit } from '@/lib/rdkit-context';

interface MoleculeViewProps {
  smiles: string;
  width?: number;
  height?: number;
  showSmiles?: boolean;
}

export const MoleculeView: React.FC<MoleculeViewProps> = ({
  smiles,
  width = 200,
  height = 200,
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
 * Compact molecule view for table cells
 */
interface MoleculeChipProps {
  smiles: string;
  size?: number;
}

export const MoleculeChip: React.FC<MoleculeChipProps> = ({ smiles, size = 60 }) => {
  const { rdkitModule, isLoading } = useRDKit();
  const [svgUrl, setSvgUrl] = useState<string | null>(null);

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
    <Box
      sx={{
        width: size,
        height: size,
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'center',
      }}
    >
      <img
        src={svgUrl}
        alt=""
        style={{ width: size, height: size, objectFit: 'contain' }}
      />
    </Box>
  );
};
