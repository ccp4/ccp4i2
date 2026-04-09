'use client';

import {
  Card,
  CardContent,
  CardActionArea,
  Typography,
  Chip,
  Box,
  Stack,
  Tooltip,
} from '@mui/material';
import { Folder, WorkOutline } from '@mui/icons-material';
import { DashboardProject } from '@/types/compounds/models';

interface ProjectCardProps {
  project: DashboardProject;
  onClick: () => void;
}

export function ProjectCard({ project, onClick }: ProjectCardProps) {
  return (
    <Card
      sx={{
        height: '100%',
        display: 'flex',
        flexDirection: 'column',
      }}
    >
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
              bgcolor: 'primary.lighter',
            }}
          >
            <Folder sx={{ fontSize: 36, color: 'primary.main' }} />
          </Box>

          {/* Project name */}
          <Tooltip title={project.name} placement="top">
            <Typography
              variant="subtitle2"
              align="center"
              noWrap
              sx={{
                width: '100%',
                fontWeight: 600,
              }}
            >
              {project.name}
            </Typography>
          </Tooltip>

          {/* Stats */}
          <Stack
            direction="row"
            spacing={0.5}
            sx={{ mt: 1, flexWrap: 'wrap', justifyContent: 'center' }}
          >
            <Chip
              icon={<WorkOutline fontSize="small" />}
              label={`${project.job_count} jobs`}
              size="small"
              variant="outlined"
            />
          </Stack>

          {/* Matched compounds */}
          {project.matching_compound_ids.length > 0 && (
            <Typography
              variant="caption"
              color="text.secondary"
              align="center"
              sx={{ mt: 1 }}
              title={project.matching_compound_ids.join(', ')}
            >
              {project.matching_compound_ids.length === 1
                ? project.matching_compound_ids[0]
                : `${project.matching_compound_ids.length} compounds`}
            </Typography>
          )}

          {/* Last access date */}
          <Typography
            variant="caption"
            color="text.secondary"
            align="center"
            sx={{ mt: 'auto', pt: 1 }}
          >
            {new Date(project.last_access).toLocaleDateString()}
          </Typography>
        </CardContent>
      </CardActionArea>
    </Card>
  );
}
