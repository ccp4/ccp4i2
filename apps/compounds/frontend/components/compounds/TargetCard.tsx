'use client';

import { Box, Card, CardActionArea, CardContent, Chip, Typography } from '@mui/material';
import { Science, Biotech, FiberNew } from '@mui/icons-material';
import { Target } from '@/types/compounds/models';
import { AuthenticatedImage } from './AuthenticatedImage';

interface TargetCardProps {
  target: Target;
  onClick?: () => void;
}

export function TargetCard({ target, onClick }: TargetCardProps) {
  const hasRecentActivity = target.has_recent_compounds || target.has_recent_assays;

  return (
    <Card
      sx={{
        height: '100%',
        display: 'flex',
        flexDirection: 'column',
        position: 'relative',
      }}
    >
      {/* Recent activity badge */}
      {hasRecentActivity && (
        <Box
          sx={{
            position: 'absolute',
            top: 8,
            right: 8,
            zIndex: 1,
            display: 'flex',
            gap: 0.5,
          }}
        >
          {target.has_recent_compounds && (
            <Chip
              icon={<FiberNew />}
              label="Compounds"
              size="small"
              color="primary"
              sx={{
                bgcolor: 'primary.main',
                color: 'white',
                '& .MuiChip-icon': { color: 'white' },
                fontSize: '0.7rem',
                height: 24,
              }}
            />
          )}
          {target.has_recent_assays && (
            <Chip
              icon={<FiberNew />}
              label="Assays"
              size="small"
              color="secondary"
              sx={{
                bgcolor: 'secondary.main',
                color: 'white',
                '& .MuiChip-icon': { color: 'white' },
                fontSize: '0.7rem',
                height: 24,
              }}
            />
          )}
        </Box>
      )}

      <CardActionArea onClick={onClick} sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column', alignItems: 'stretch' }}>
        {/* Target image or placeholder */}
        {target.image ? (
          <AuthenticatedImage
            src={target.image}
            alt={target.name}
            width="100%"
            height={120}
            objectFit="cover"
          />
        ) : (
          <Box
            sx={{
              height: 120,
              bgcolor: 'grey.100',
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
            }}
          >
            <Science sx={{ fontSize: 48, color: 'grey.400' }} />
          </Box>
        )}

        <CardContent sx={{ flexGrow: 1, display: 'flex', flexDirection: 'column' }}>
          {/* Target name */}
          <Typography
            variant="h6"
            component="div"
            sx={{
              fontWeight: 600,
              overflow: 'hidden',
              textOverflow: 'ellipsis',
              whiteSpace: 'nowrap',
            }}
          >
            {target.name}
          </Typography>

          {/* Last activity */}
          {target.latest_activity && (
            <Typography variant="caption" color="text.secondary" sx={{ mb: 1 }}>
              Last activity: {new Date(target.latest_activity).toLocaleDateString()}
            </Typography>
          )}

          {/* Count chips */}
          <Box sx={{ display: 'flex', gap: 1, mt: 'auto', flexWrap: 'wrap' }}>
            <Chip
              icon={<Science fontSize="small" />}
              label={`${target.compound_count ?? 0} compounds`}
              size="small"
              variant="outlined"
              color="primary"
            />
            <Chip
              icon={<Biotech fontSize="small" />}
              label={`${target.assay_count ?? 0} assays`}
              size="small"
              variant="outlined"
              color="secondary"
            />
          </Box>
        </CardContent>
      </CardActionArea>
    </Card>
  );
}
