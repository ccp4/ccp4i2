'use client';

import { useMemo } from 'react';
import { Box, Skeleton, Typography } from '@mui/material';
import { MoleculeChip } from './MoleculeView';
import { useCompoundsApi } from '@/lib/api';
import { Compound } from '@/types/models';

interface CompoundStructureCellProps {
  /** UUID of the compound if already matched */
  compoundId?: string | null;
  /** Compound name that may contain NCL-XXXXX pattern */
  compoundName?: string | null;
  /** Size of the molecule image */
  size?: number;
}

/**
 * Extract NCL-number from a string (e.g., "NCL-00027145" from compound name)
 * Returns the formatted_id pattern if found
 */
function extractNclId(name: string | null | undefined): string | null {
  if (!name) return null;

  // Match NCL- followed by digits (with optional leading zeros)
  const match = name.match(/NCL-?(\d+)/i);
  if (match) {
    // Normalize to NCL-XXXXXXXX format (8 digits with leading zeros)
    const num = parseInt(match[1], 10);
    return `NCL-${num.toString().padStart(8, '0')}`;
  }
  return null;
}

/**
 * Component that displays a compound structure for a data series row.
 *
 * Resolution order:
 * 1. If compoundId is provided, fetch compound directly by ID
 * 2. If compoundName contains NCL-XXXXX pattern, search for matching compound
 * 3. Otherwise show empty placeholder
 */
export function CompoundStructureCell({
  compoundId,
  compoundName,
  size = 80,
}: CompoundStructureCellProps) {
  const api = useCompoundsApi();

  // Extract NCL ID from compound name if no direct compound ID
  const nclId = useMemo(() => {
    if (compoundId) return null; // Don't need to extract if we have direct ID
    return extractNclId(compoundName);
  }, [compoundId, compoundName]);

  // Fetch compound by ID if available
  const { data: compoundById, isLoading: loadingById } = api.get<Compound>(
    compoundId ? `compounds/${compoundId}/` : null
  );

  // Search for compound by NCL ID if no direct match
  const { data: compoundsByNcl, isLoading: loadingByNcl } = api.get<Compound[]>(
    !compoundId && nclId ? `compounds/?formatted_id=${nclId}` : null
  );

  // Determine which compound to use
  const compound = compoundById || (compoundsByNcl && compoundsByNcl.length > 0 ? compoundsByNcl[0] : null);
  const isLoading = loadingById || loadingByNcl;
  const smiles = compound?.smiles || compound?.rdkit_smiles;

  // No compound ID and no NCL pattern found
  if (!compoundId && !nclId) {
    return (
      <Box
        sx={{
          width: size,
          height: size,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: 'grey.50',
          borderRadius: 1,
        }}
      >
        <Typography variant="caption" color="text.secondary">
          -
        </Typography>
      </Box>
    );
  }

  if (isLoading) {
    return <Skeleton variant="rectangular" width={size} height={size} sx={{ borderRadius: 1 }} />;
  }

  if (!smiles) {
    return (
      <Box
        sx={{
          width: size,
          height: size,
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center',
          bgcolor: 'grey.100',
          borderRadius: 1,
        }}
      >
        <Typography variant="caption" color="text.secondary">
          No structure
        </Typography>
      </Box>
    );
  }

  return <MoleculeChip smiles={smiles} size={size} />;
}
