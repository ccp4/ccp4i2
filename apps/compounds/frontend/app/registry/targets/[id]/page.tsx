'use client';

import { use, useRef, ChangeEvent } from 'react';
import { useRouter } from 'next/navigation';
import {
  Container,
  Typography,
  Box,
  Paper,
  Skeleton,
  Button,
  IconButton,
  Tooltip,
} from '@mui/material';
import {
  Science,
  TableChart,
  Add,
  CloudUpload,
  Delete,
  ViewList,
} from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { HorizontalCarousel } from '@/components/compounds/HorizontalCarousel';
import { CompoundCard } from '@/components/compounds/CompoundCard';
import { AssayCard } from '@/components/compounds/AssayCard';
import { ProjectCard } from '@/components/compounds/ProjectCard';
import { AddAssayMenu } from '@/components/compounds/AddAssayMenu';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import {
  TargetDashboard,
  DashboardProject,
} from '@/types/compounds/models';

interface PageProps {
  params: Promise<{ id: string }>;
}

export default function TargetDashboardPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const { canContribute } = useAuth();
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Fetch dashboard data
  const {
    data: dashboardData,
    isLoading: dashboardLoading,
    mutate: mutateDashboard,
  } = api.get<TargetDashboard>(`targets/${id}/dashboard/`);

  // Fetch recent projects
  const { data: projects, isLoading: projectsLoading } = api.get<
    DashboardProject[]
  >(`targets/${id}/recent_projects/`);

  const handleImageUpload = async (e: ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const formData = new FormData();
    formData.append('image', file);

    try {
      await api.upload(`targets/${id}/upload_image/`, formData);
      mutateDashboard();
    } catch (error) {
      console.error('Failed to upload image:', error);
    }

    // Reset file input
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  const handleDeleteImage = async () => {
    try {
      await api.delete(`targets/${id}/delete_image/`);
      mutateDashboard();
    } catch (error) {
      console.error('Failed to delete image:', error);
    }
  };

  const navigateToCompounds = () => {
    router.push(routes.registry.targetCompounds(id));
  };

  if (dashboardLoading) {
    return (
      <Container maxWidth="lg" sx={{ py: 3 }}>
        <Skeleton variant="rectangular" height={200} sx={{ mb: 3 }} />
        <Skeleton variant="rectangular" height={300} sx={{ mb: 3 }} />
        <Skeleton variant="rectangular" height={200} sx={{ mb: 3 }} />
        <Skeleton variant="rectangular" height={200} />
      </Container>
    );
  }

  const target = dashboardData;

  return (
    <Container maxWidth="lg" sx={{ py: 3 }}>
      <PageHeader
        breadcrumbs={[
          { label: 'Home', href: routes.home(), icon: 'home' },
          { label: 'Targets', href: routes.registry.targets(), icon: 'target' },
          { label: target?.name || 'Loading...', icon: 'target' },
        ]}
      />

      {/* Branding Image & Header */}
      <Paper sx={{ p: 3, mb: 3, position: 'relative' }}>
        {/* Branding image */}
        {target?.image ? (
          <Box sx={{ position: 'relative', mb: 2 }}>
            <Box
              component="img"
              src={target.image}
              alt={`${target.name} banner`}
              sx={{
                width: '100%',
                maxHeight: 200,
                objectFit: 'cover',
                borderRadius: 1,
              }}
            />
            <Tooltip title={canContribute ? "Remove image" : "Requires Contributor or Admin operating level"}>
              <span style={{ position: 'absolute', top: 8, right: 8 }}>
                <IconButton
                  onClick={handleDeleteImage}
                  size="small"
                  disabled={!canContribute}
                  sx={{
                    bgcolor: 'background.paper',
                    '&:hover': { bgcolor: canContribute ? 'error.lighter' : undefined },
                  }}
                >
                  <Delete fontSize="small" />
                </IconButton>
              </span>
            </Tooltip>
          </Box>
        ) : (
          <Box
            sx={{
              height: 120,
              bgcolor: 'grey.100',
              borderRadius: 1,
              mb: 2,
              display: 'flex',
              alignItems: 'center',
              justifyContent: 'center',
            }}
          >
              <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
              <span>
                <Button
                  variant="outlined"
                  startIcon={<CloudUpload />}
                  component="label"
                  disabled={!canContribute}
                >
                  Upload Branding Image
                  <input
                    ref={fileInputRef}
                    type="file"
                    hidden
                    accept="image/*"
                    onChange={handleImageUpload}
                    disabled={!canContribute}
                  />
                </Button>
              </span>
            </Tooltip>
          </Box>
        )}

        {/* Target info and actions */}
        <Box
          sx={{
            display: 'flex',
            alignItems: 'center',
            justifyContent: 'space-between',
            flexWrap: 'wrap',
            gap: 2,
          }}
        >
          <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
            <Science sx={{ fontSize: 48, color: 'primary.main' }} />
            <Box>
              <Typography variant="h4">{target?.name}</Typography>
              <Typography color="text.secondary">
                {target?.compound_count ?? 0} compounds registered
              </Typography>
            </Box>
          </Box>
          <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
            <Button
              component={Link}
              href={routes.registry.targetCompounds(id)}
              variant="outlined"
              startIcon={<ViewList />}
            >
              All Compounds
            </Button>
            <Button
              component={Link}
              href={routes.assays.aggregate({ target: id })}
              variant="outlined"
              startIcon={<TableChart />}
            >
              View Assay Data
            </Button>
            <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
              <span>
                <Button
                  component={Link}
                  href={`${routes.registry.new()}?target=${id}`}
                  variant="contained"
                  startIcon={<Add />}
                  disabled={!canContribute}
                >
                  New Compound
                </Button>
              </span>
            </Tooltip>
          </Box>
        </Box>
      </Paper>

      {/* Recent Compounds Carousel */}
      <HorizontalCarousel
        items={dashboardData?.recent_compounds || []}
        title="Recent Compounds"
        renderItem={(compound) => (
          <CompoundCard
            compound={compound}
            onClick={() => router.push(routes.registry.compound(compound.id))}
          />
        )}
        getItemKey={(c) => c.id}
        onBackgroundClick={navigateToCompounds}
        itemWidth={180}
        height={260}
        emptyMessage="No compounds registered yet"
        headerAction={
          <Tooltip title={canContribute ? '' : 'Requires Contributor or Admin operating level'} arrow>
            <span>
              <Button
                component={Link}
                href={`${routes.registry.new()}?target=${id}`}
                variant="outlined"
                size="small"
                startIcon={<Add />}
                disabled={!canContribute}
              >
                Register
              </Button>
            </span>
          </Tooltip>
        }
      />

      {/* Recent Assays Carousel */}
      <HorizontalCarousel
        items={dashboardData?.recent_assays || []}
        title="Recent Assays"
        renderItem={(assay) => (
          <AssayCard
            assay={assay}
            onClick={() => router.push(routes.assays.detail(assay.id))}
          />
        )}
        getItemKey={(a) => a.id}
        itemWidth={180}
        height={220}
        emptyMessage="No assays for this target"
        headerAction={<AddAssayMenu targetId={id} disabled={!canContribute} />}
      />

      {/* Related CCP4i2 Projects Carousel */}
      <HorizontalCarousel
        items={projects || []}
        title="Related CCP4i2 Projects"
        loading={projectsLoading}
        renderItem={(project) => (
          <ProjectCard
            project={project}
            onClick={() => router.push(routes.external.ccp4i2Project(project.id))}
          />
        )}
        getItemKey={(p) => String(p.id)}
        itemWidth={200}
        height={220}
        emptyMessage="No matching projects found"
      />

      {/* Lead Compounds Table (Placeholder) */}
      <Paper sx={{ p: 3 }}>
        <Typography variant="h6" gutterBottom>
          Lead Compounds
        </Typography>
        <Box
          sx={{
            p: 4,
            bgcolor: 'grey.50',
            borderRadius: 1,
            textAlign: 'center',
          }}
        >
          <Typography color="text.secondary">
            Configurable lead compounds table - coming in a future release
          </Typography>
        </Box>
      </Paper>
    </Container>
  );
}
