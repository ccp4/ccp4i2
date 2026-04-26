'use client';

import { use, useRef, useState, useEffect, ChangeEvent } from 'react';
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
  CircularProgress,
  Alert,
} from '@mui/material';
import {
  Science,
  TableChart,
  Add,
  CloudUpload,
  Delete,
  ViewList,
  Settings,
  Radar,
} from '@mui/icons-material';
import Link from 'next/link';
import { PageHeader } from '@/components/compounds/PageHeader';
import { HorizontalCarousel } from '@/components/compounds/HorizontalCarousel';
import { CompoundCard } from '@/components/compounds/CompoundCard';
import { AssayCard } from '@/components/compounds/AssayCard';
import { ProjectCard } from '@/components/compounds/ProjectCard';
import { AddAssayMenu } from '@/components/compounds/AddAssayMenu';
import { AggregationTable } from '@/components/compounds/AggregationTable';
import { AuthenticatedImage } from '@/components/compounds/AuthenticatedImage';
import { TargetGenesSection } from '@/components/compounds/TargetGenesSection';
import { useCompoundsApi } from '@/lib/compounds/api';
import { useAuth } from '@/lib/compounds/auth-context';
import { routes } from '@/lib/compounds/routes';
import {
  fetchAggregation,
  fetchProtocols,
  deleteAggregationView,
} from '@/lib/compounds/aggregation-api';
import { scorecardDataNeeds } from '@/lib/compounds/scorecard';
import {
  TargetDashboard,
  DashboardProject,
} from '@/types/compounds/models';
import { AggregationResponse } from '@/types/compounds/aggregation';

interface PageProps {
  params: Promise<{ id: string }>;
}

export default function TargetDashboardPage({ params }: PageProps) {
  const { id } = use(params);
  const router = useRouter();
  const api = useCompoundsApi();
  const { canContribute, canAdminister } = useAuth();
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

  // Saved aggregation view state
  const [aggregationData, setAggregationData] = useState<AggregationResponse | null>(null);
  const [aggregationLoading, setAggregationLoading] = useState(false);
  const [aggregationError, setAggregationError] = useState<string | null>(null);

  // Fetch aggregation data when dashboard loads and has a saved view
  useEffect(() => {
    const savedView = dashboardData?.saved_aggregation_view;
    if (!savedView || !dashboardData) return;

    const runAggregation = async () => {
      setAggregationLoading(true);
      setAggregationError(null);

      try {
        // First fetch protocols to convert names to IDs
        const allProtocols = await fetchProtocols({ target: id });
        const matchedProtocols = allProtocols.filter((p) =>
          savedView.protocol_names.some((name) => p.name.toLowerCase() === name.toLowerCase())
        );

        // Auto-inject the scorecard's required protocols and molecular
        // properties so a saved cards view doesn't show "no data" for
        // axes the scorecard depends on (Lipinski needs molecular_weight
        // / clogp / hbd / hba; protocol/ratio/worst_of axes need their
        // referenced protocol IDs to be in the request). Same logic as
        // the aggregation page so behaviour matches.
        const scorecardNeeds = scorecardDataNeeds(dashboardData.scorecard_config);
        const scorecardProtocolIds = scorecardNeeds.protocolIds.length > 0
          ? new Set([
              ...matchedProtocols.map((p) => p.id),
              ...scorecardNeeds.protocolIds,
            ])
          : new Set(matchedProtocols.map((p) => p.id));

        // Build predicates
        const predicates: any = {
          targets: [id],
        };
        if (scorecardProtocolIds.size > 0) {
          predicates.protocols = Array.from(scorecardProtocolIds);
        }
        if (savedView.compound_search) {
          predicates.compound_search = savedView.compound_search;
        }
        // Always include status (empty string means "all")
        predicates.status = savedView.status ?? 'valid';

        const effectiveProperties = Array.from(
          new Set([
            ...(savedView.include_properties ?? []),
            ...scorecardNeeds.properties,
          ]),
        );

        // Run aggregation - map pivot/cards to compact for API (they use same data structure)
        const apiOutputFormat = (savedView.output_format === 'pivot' || savedView.output_format === 'cards')
          ? 'compact'
          : savedView.output_format;
        const result = await fetchAggregation({
          predicates,
          output_format: apiOutputFormat,
          aggregations: savedView.aggregations,
          // Cast to MolecularPropertyName[] - backend validates the values
          include_properties: effectiveProperties.length > 0 ? effectiveProperties as any : undefined,
          // Default to showing identifiers when the saved view predates the flag
          include_identifiers: savedView.include_identifiers ?? true,
        });
        setAggregationData(result);
      } catch (err) {
        setAggregationError(err instanceof Error ? err.message : 'Failed to load aggregation data');
      } finally {
        setAggregationLoading(false);
      }
    };

    runAggregation();
  }, [dashboardData, id]);

  const handleDeleteSavedView = async () => {
    try {
      await deleteAggregationView(id);
      mutateDashboard();
      setAggregationData(null);
    } catch (err) {
      console.error('Failed to delete saved view:', err);
    }
  };

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
            <AuthenticatedImage
              src={target.image}
              alt={`${target.name} banner`}
              width="100%"
              height={200}
              objectFit="cover"
              sx={{ borderRadius: 1 }}
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
            <Button
              component={Link}
              href={routes.registry.targetScorecard(id)}
              variant="outlined"
              startIcon={<Radar />}
            >
              Scorecard
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

        {/* Genes section — linked Gene rows with HGNC metadata; inline edit */}
        {/* for Contributors+. Updates dashboard cache so carousels etc. see */}
        {/* the new target state without a full re-fetch. */}
        {target && (
          <TargetGenesSection
            target={target}
            canEdit={canContribute}
            onUpdated={(updated) =>
              mutateDashboard({ ...target, ...updated } as TargetDashboard, false)
            }
          />
        )}
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

      {/* Saved Aggregation View */}
      <Paper sx={{ p: 3 }}>
        <Box sx={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', mb: 2 }}>
          <Typography variant="h6">
            Lead Compounds
          </Typography>
          {canAdminister && dashboardData?.saved_aggregation_view && (
            <Box sx={{ display: 'flex', gap: 1 }}>
              <Tooltip title="Configure in Data Aggregation">
                <IconButton
                  component={Link}
                  href={routes.assays.aggregate({
                    targets: [dashboardData.name],
                    protocols: dashboardData.saved_aggregation_view.protocol_names,
                    compound: dashboardData.saved_aggregation_view.compound_search || undefined,
                    format: dashboardData.saved_aggregation_view.output_format,
                    aggregations: dashboardData.saved_aggregation_view.aggregations,
                    status: dashboardData.saved_aggregation_view.status || undefined,
                    concentrationDisplay: dashboardData.saved_aggregation_view.concentration_display,
                    properties: dashboardData.saved_aggregation_view.include_properties,
                  })}
                  size="small"
                >
                  <Settings fontSize="small" />
                </IconButton>
              </Tooltip>
              <Tooltip title="Remove saved view">
                <IconButton
                  onClick={handleDeleteSavedView}
                  size="small"
                  color="error"
                >
                  <Delete fontSize="small" />
                </IconButton>
              </Tooltip>
            </Box>
          )}
        </Box>

        {!dashboardData?.saved_aggregation_view ? (
          <Box
            sx={{
              p: 4,
              bgcolor: 'grey.50',
              borderRadius: 1,
              textAlign: 'center',
            }}
          >
            <Typography color="text.secondary" sx={{ mb: canAdminister ? 2 : 0 }}>
              No aggregation view configured for this target
            </Typography>
            {canAdminister && (
              <Button
                component={Link}
                href={routes.assays.aggregate({ targets: [dashboardData?.name || ''] })}
                variant="outlined"
                size="small"
                startIcon={<Settings />}
              >
                Configure in Data Aggregation
              </Button>
            )}
          </Box>
        ) : aggregationLoading ? (
          <Box sx={{ display: 'flex', justifyContent: 'center', p: 4 }}>
            <CircularProgress />
          </Box>
        ) : aggregationError ? (
          <Alert severity="error">{aggregationError}</Alert>
        ) : (
          <AggregationTable
            data={aggregationData}
            loading={false}
            aggregations={dashboardData.saved_aggregation_view.aggregations}
            outputFormat={dashboardData.saved_aggregation_view.output_format}
            concentrationDisplay={dashboardData.saved_aggregation_view.concentration_display}
            // Pass the target's scorecard so the cards view can render
            // spider / values / compact modes when one is configured. If
            // no scorecard exists, CardsView falls back to protocols-only
            // automatically (its `hasScorecard ? 'both' : 'protocols'`
            // default).
            scorecardConfig={dashboardData.scorecard_config}
          />
        )}
      </Paper>
    </Container>
  );
}
