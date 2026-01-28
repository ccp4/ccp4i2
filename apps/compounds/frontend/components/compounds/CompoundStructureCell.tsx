'use client';

import { useMemo } from 'react';
import { Box, Skeleton, Typography } from '@mui/material';
import { MoleculeChip } from './MoleculeView';
import { useCompoundsApi } from '@/lib/compounds/api';
import { getCompoundPattern, formatCompoundId } from '@/lib/compounds/config';
import { Compound } from '@/types/compounds/models';

interface CompoundStructureCellProps {
  /** UUID of the compound if already matched */
  compoundId?: string | null;
  /** Compound name that may contain compound ID pattern (e.g., NCL-XXXXX) */
  compoundName?: string | null;
  /** Size of the molecule image */
  size?: number;
}

/**
 * Extract compound ID from a string (e.g., "NCL-00027145" from compound name)
 * Returns the formatted_id pattern if found, using the configured prefix
 */
function extractCompoundId(name: string | null | undefined): string | null {
  if (!name) return null;

  // Match configured prefix followed by digits (with optional dash and leading zeros)
  const pattern = getCompoundPattern();
  const match = name.match(pattern);
  if (match) {
    // Normalize to configured format (e.g., PREFIX-XXXXXXXX)
    const num = parseInt(match[1], 10);
    return formatCompoundId(num);
  }
  return null;
}

/**
 * Component that displays a compound structure for a data series row.
 *
 * Resolution order:
 * 1. If compoundId is provided, fetch compound directly by ID
 * 2. If compoundName contains compound ID pattern, search for matching compound
 * 3. Otherwise show empty placeholder
 */
export function CompoundStructureCell({
  compoundId,
  compoundName,
  size = 100,
}: CompoundStructureCellProps) {
  const api = useCompoundsApi();

  // Extract compound ID from compound name if no direct compound ID
  const extractedId = useMemo(() => {
    if (compoundId) return null; // Don't need to extract if we have direct ID
    return extractCompoundId(compoundName);
  }, [compoundId, compoundName]);

  // Fetch compound by ID if available
  const { data: compoundById, isLoading: loadingById } = api.get<Compound>(
    compoundId ? `compounds/${compoundId}/` : null
  );

  // Search for compound by formatted ID if no direct match
  const { data: compoundsByFormattedId, isLoading: loadingByFormattedId } = api.get<Compound[]>(
    !compoundId && extractedId ? `compounds/?formatted_id=${extractedId}` : null
  );

  // Determine which compound to use
  const compound = compoundById || (compoundsByFormattedId && compoundsByFormattedId.length > 0 ? compoundsByFormattedId[0] : null);
  const isLoading = loadingById || loadingByFormattedId;
  const smiles = compound?.smiles || compound?.rdkit_smiles;

  // No compound ID and no compound ID pattern found
  if (!compoundId && !extractedId) {
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
