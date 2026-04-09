'use client';

import {
  Card,
  CardContent,
  CardActionArea,
  Typography,
  Chip,
  Box,
} from '@mui/material';
import { Assessment, Science, FiberNew } from '@mui/icons-material';
import { DashboardAssay } from '@/types/compounds/models';

interface AssayCardProps {
  assay: DashboardAssay;
  onClick: () => void;
}

function isWithinLastWeek(dateString: string): boolean {
  const date = new Date(dateString);
  const now = new Date();
  const weekAgo = new Date(now.getTime() - 7 * 24 * 60 * 60 * 1000);
  return date >= weekAgo;
}

export function AssayCard({ assay, onClick }: AssayCardProps) {
  const isNew = isWithinLastWeek(assay.created_at);

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
          {/* Icon */}
          <Box
            sx={{
              display: 'flex',
              justifyContent: 'center',
              alignItems: 'center',
              mb: 1.5,
              width: 64,
              height: 64,
              borderRadius: 2,
              bgcolor: 'info.lighter',
            }}
          >
            <Assessment sx={{ fontSize: 36, color: 'info.main' }} />
          </Box>

          {/* Protocol name */}
          <Typography
            variant="subtitle2"
            align="center"
            noWrap
            title={assay.protocol_name}
            sx={{
              width: '100%',
              fontWeight: 600,
            }}
          >
            {assay.protocol_name}
          </Typography>

          {/* Data series count */}
          <Chip
            icon={<Science fontSize="small" />}
            label={`${assay.data_series_count} series`}
            size="small"
            variant="outlined"
            color="info"
            sx={{ mt: 1 }}
          />

          {/* Creation date */}
          <Typography
            variant="caption"
            color="text.secondary"
            align="center"
            sx={{ mt: 'auto', pt: 1 }}
          >
            {new Date(assay.created_at).toLocaleDateString()}
          </Typography>
        </CardContent>
      </CardActionArea>
    </Card>
  );
}
