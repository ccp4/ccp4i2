"use client";
import { Suspense, useMemo, useState, useCallback, useEffect } from "react";
import { useParams, useSearchParams } from "next/navigation";
import { ClientStoreProvider } from "../../../../../providers/client-store-provider";
import { useCampaignsApi } from "../../../../../lib/campaigns-api";
import { useMoorhenBreadcrumbs } from "../../../../../providers/moorhen-breadcrumb-context";
import CampaignMoorhenWrapper from "../../../../../components/moorhen/campaign-moorhen-wrapper";

// Inner component that uses useSearchParams (requires Suspense boundary)
function CampaignPageContent() {
  const { id } = useParams();
  const searchParams = useSearchParams();
  const viewParam = searchParams.get("view");
  const jobParam = searchParams.get("job"); // Optional: specific job to load
  const campaignId = id ? parseInt(id as string) : null;
  const initialJobId = jobParam ? parseInt(jobParam) : null;

  const campaignsApi = useCampaignsApi();
  const { data: campaign } = campaignsApi.useCampaign(campaignId);
  const { data: parentProject } = campaignsApi.useParentProject(campaignId);
  const { data: parentFiles } = campaignsApi.useParentFiles(campaignId);
  const { data: memberProjects } = campaignsApi.useMemberProjects(campaignId);
  const { data: sites, mutate: mutateSites } = campaignsApi.useSites(campaignId);

  // Track which member project is currently selected (null = parent)
  const [selectedMemberProjectId, setSelectedMemberProjectId] = useState<number | null>(null);
  // Track specific job to load (from URL param or selection)
  const [selectedJobId, setSelectedJobId] = useState<number | null>(initialJobId);
  // Track if we've initialized from URL params
  const [initialized, setInitialized] = useState(false);

  // When member projects load and we have a job param, find which project it belongs to
  useEffect(() => {
    if (initialized || !memberProjects || !initialJobId) return;

    // Find the project containing this job
    for (const project of memberProjects) {
      const job = project.jobs?.find((j) => j.id === initialJobId);
      if (job) {
        setSelectedMemberProjectId(project.id);
        setSelectedJobId(initialJobId);
        setInitialized(true);
        return;
      }
    }
    // Job not found in member projects - might be in parent project
    setInitialized(true);
  }, [memberProjects, initialJobId, initialized]);

  // When user selects a different project, reset to latest job
  const handleSelectMemberProject = useCallback((projectId: number | null) => {
    setSelectedMemberProjectId(projectId);
    setSelectedJobId(null); // Reset to auto-select latest job
  }, []);

  // Determine which files to load based on selection
  const fileIds = useMemo(() => {
    // If a specific job is selected, load that job
    if (selectedJobId !== null) {
      return { type: "job" as const, jobId: selectedJobId };
    }

    if (selectedMemberProjectId !== null) {
      // Find the selected member project
      const memberProject = memberProjects?.find(
        (p) => p.id === selectedMemberProjectId
      );
      if (memberProject?.jobs) {
        // Find the most recent finished refmac or dimple job
        const latestJob = [...memberProject.jobs]
          .filter((j) => j.status === 6) // FINISHED
          .filter((j) =>
            ["refmac", "i2Refmac", "i2Dimple", "dimple"].includes(j.task_name)
          )
          .sort((a, b) => b.id - a.id)[0];

        if (latestJob) {
          return { type: "job" as const, jobId: latestJob.id };
        }
      }
      return { type: "none" as const };
    }

    // Default: load parent project files (coordinates + maps)
    if (!parentFiles) {
      return { type: "none" as const };
    }

    const ids: number[] = [];
    // Add coordinate files
    if (parentFiles.coordinates) {
      ids.push(...parentFiles.coordinates.map((f) => f.id));
    }
    // Note: maps would need to be added here when parent has map files
    return { type: "files" as const, fileIds: ids };
  }, [parentFiles, selectedMemberProjectId, selectedJobId, memberProjects]);

  // Build breadcrumbs for the layout AppBar
  const { setBreadcrumbs } = useMoorhenBreadcrumbs();
  const currentJobId = fileIds.type === "job" ? fileIds.jobId : null;

  useEffect(() => {
    if (!campaign) return;

    const crumbs: { label: string; href: string }[] = [
      { label: campaign.name, href: `/ccp4i2/campaigns/${campaignId}` },
    ];

    if (selectedMemberProjectId !== null) {
      const project = memberProjects?.find(
        (p) => p.id === selectedMemberProjectId
      );
      if (project) {
        crumbs.push({
          label: project.name,
          href: `/ccp4i2/project/${project.id}`,
        });

        if (currentJobId !== null) {
          const job = project.jobs?.find((j) => j.id === currentJobId);
          if (job) {
            crumbs.push({
              label: job.title || `${job.task_name} #${job.number}`,
              href: `/ccp4i2/project/${project.id}/job/${job.id}`,
            });
          }
        }
      }
    } else if (parentProject) {
      crumbs.push({
        label: parentProject.name,
        href: `/ccp4i2/project/${parentProject.id}`,
      });
    }

    setBreadcrumbs(crumbs);
    return () => setBreadcrumbs([]);
  }, [
    campaign,
    campaignId,
    parentProject,
    selectedMemberProjectId,
    currentJobId,
    memberProjects,
    setBreadcrumbs,
  ]);

  // Handle site update
  const handleUpdateSites = useCallback(
    async (newSites: typeof sites) => {
      if (!campaignId || !newSites) return;
      await campaignsApi.updateSites(campaignId, newSites);
      mutateSites();
    },
    [campaignId, campaignsApi, mutateSites]
  );

  if (!campaign) {
    return <div>Loading campaign...</div>;
  }

  return (
    <ClientStoreProvider>
      <CampaignMoorhenWrapper
        campaign={campaign}
        fileSource={fileIds}
        viewParam={viewParam}
        sites={sites || []}
        onUpdateSites={handleUpdateSites}
        memberProjects={memberProjects || []}
        selectedMemberProjectId={selectedMemberProjectId}
        onSelectMemberProject={handleSelectMemberProject}
        parentProject={parentProject}
      />
    </ClientStoreProvider>
  );
}

// Outer component with Suspense boundary for useSearchParams
const CampaignPage = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <CampaignPageContent />
    </Suspense>
  );
};

export default CampaignPage;
