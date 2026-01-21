'use client';

import { Card, CardContent, CardActionArea, Typography, Box, Chip } from '@mui/material';
import { FiberNew } from '@mui/icons-material';
import { MoleculeChip } from './MoleculeView';
import { DashboardCompound } from '@/types/compounds/models';

interface CompoundCardProps {
  compound: DashboardCompound;
  onClick: () => void;
}

function isWithinLastWeek(dateString: string): boolean {
  const date = new Date(dateString);
  const now = new Date();
  const weekAgo = new Date(now.getTime() - 7 * 24 * 60 * 60 * 1000);
  return date >= weekAgo;
}

export function CompoundCard({ compound, onClick }: CompoundCardProps) {
  const isNew = isWithinLastWeek(compound.registered_at);

  return (
    <Card
      sx={{
        height: '100%',
        display: 'flex',
        flexDirection: 'column',
        position: 'relative',
      }}
    >
      {/* New badge */}
      {isNew && (
        <Chip
          icon={<FiberNew />}
          label="New"
          size="small"
          color="secondary"
          sx={{
            position: 'absolute',
            top: 8,
            right: 8,
            zIndex: 1,
            fontSize: '0.7rem',
            height: 24,
          }}
        />
      )}
      <CardActionArea
        onClick={onClick}
        sx={{
          height: '100%',
          display: 'flex',
          flexDirection: 'column',
          alignItems: 'stretch',
        }}
      >
        <CardContent
          sx={{
            display: 'flex',
            flexDirection: 'column',
            alignItems: 'center',
            flex: 1,
            p: 1.5,
          }}
        >
          {/* Structure image */}
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              mb: 1,
              minHeight: 120,
            }}
          >
            <MoleculeChip smiles={compound.smiles} size={110} />
          </Box>

          {/* Compound ID */}
          <Typography
            variant="subtitle2"
            fontFamily="monospace"
            align="center"
            sx={{
              fontWeight: 600,
              color: 'primary.main',
            }}
          >
            {compound.formatted_id}
          </Typography>

          {/* Molecular weight */}
          {compound.molecular_weight && (
            <Typography
              variant="caption"
              color="text.secondary"
              align="center"
            >
              MW: {compound.molecular_weight.toFixed(1)}
            </Typography>
          )}

          {/* Registration date */}
          <Typography
            variant="caption"
            color="text.secondary"
            align="center"
            sx={{ mt: 'auto' }}
          >
            {new Date(compound.registered_at).toLocaleDateString()}
          </Typography>
        </CardContent>
      </CardActionArea>
    </Card>
  );
}
