"use client";
import { Suspense, useMemo, useState, useCallback, useEffect } from "react";
import { useParams, useRouter, useSearchParams } from "next/navigation";
import { ClientStoreProvider } from "@/providers/client-store-provider";
import { useCampaignsApi } from "@/lib/campaigns-api";
import { useMoorhenBreadcrumbs } from "@/providers/moorhen-breadcrumb-context";
import { usePopcorn } from "@/providers/popcorn-provider";
import CampaignMoorhenWrapper from "@/components/moorhen/campaign-moorhen-wrapper";
import type { MoorhenScene } from "@/types/moorhen-scene";
import { parseEventParam } from "@/lib/event-evidence";
import { sceneSourceFromParams } from "@/lib/scene-source";
import { parseSiteParam } from "@/lib/site-view-link";
import { listenAsViewer } from "@/lib/viewer-channel";
import type {
  NewCampaignSite,
  SiteEvaluation,
  SiteVerdict,
} from "@/types/campaigns";

// Inner component that uses useSearchParams (requires Suspense boundary)
function CampaignPageContent() {
  const params = useParams();
  const id = params?.id as string | undefined;
  const searchParams = useSearchParams();
  const viewParam = searchParams?.get("view");
  const jobParam = searchParams?.get("job"); // Optional: specific job to load
  const summaryMode = searchParams?.get("summary") === "1"; // Campaign overview
  const siteParam = searchParams?.get("site"); // Optional: site to open on
  // Optional: a PanDDA event whose evidence to load (a filled matrix box)
  const eventParam = searchParams?.get("event");
  const initialEvent = useMemo(() => parseEventParam(eventParam), [eventParam]);
  // Optional: a scene a job builds on request (a PanDDA report's "open in
  // Moorhen" button), applied over the member's model once it is in view.
  const sceneKey = ["sceneJob", "sceneMethod", "sceneArgs"].map((k) => searchParams?.get(k) ?? "").join("\n");
  const sceneSource = useMemo(
    () => (searchParams ? sceneSourceFromParams(searchParams) : null),
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [sceneKey],
  );
  const campaignId = id ? parseInt(id as string) : null;
  const initialJobId = jobParam ? parseInt(jobParam) : null;

  // Take links for this campaign and member job from other windows
  // (lib/viewer-channel): the new URL's scene is applied over what is open.
  const router = useRouter();
  useEffect(() => listenAsViewer((href) => router.replace(href)), [router]);

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

  // Summary mode: fetch the campaign overview scene.
  //
  // The deps here must stay free of unstable values. `campaignsApi` is a fresh
  // object every render, and `summaryLoading`-style state we set in the effect
  // would also churn the deps — either makes the effect re-run mid-flight, fire
  // its cleanup, and cancel the in-flight fetch, so the scene gets fetched and
  // then silently discarded. We depend only on the stable inputs and rely on
  // the `cancelled` flag (which is StrictMode-safe: the second mount re-fetches
  // and wins).
  //
  // With a site in the URL as well, the scene is that site's: the hits
  // recorded there, fitted on its pocket, rather than the whole campaign.
  // Everything downstream (serialise, seed the Scenes panel, auto-apply) is
  // the same, so it goes through the same prop.
  // `site` names one site, or several -- a site and its subsites, or
  // neighbours -- viewed as one scene (lib/site-view-link).
  const summarySiteIds = useMemo(() => parseSiteParam(siteParam), [siteParam]);
  const summarySiteKey = summarySiteIds.join(",");
  const [summaryScene, setSummaryScene] = useState<MoorhenScene | null>(null);
  const { setMessage } = usePopcorn();
  useEffect(() => {
    if (!summaryMode || campaignId === null) return;
    let cancelled = false;
    // The server builds the scene on demand and remembers it, so the first
    // open after a change can take a while on a big campaign; say so rather
    // than leave an empty viewer.
    setMessage("Building the campaign scene; a large campaign can take a minute", "info");
    const request =
      summarySiteIds.length > 1
        ? campaignsApi.fetchSitesScene(campaignId, summarySiteIds)
        : summarySiteIds.length === 1
          ? campaignsApi.fetchSiteScene(campaignId, summarySiteIds[0])
          : campaignsApi.fetchSummaryScene(campaignId);
    request
      .then((res) => {
        if (!cancelled) setSummaryScene(res.scene);
      })
      .catch((err) => {
        if (cancelled) return;
        console.error("Failed to load campaign summary scene:", err);
        setMessage(
          `Could not build the campaign scene: ${err instanceof Error ? err.message : String(err)}`,
          "error",
        );
      });
    return () => {
      cancelled = true;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [summaryMode, campaignId, summarySiteKey]);

  // Determine which files to load based on selection
  const fileIds = useMemo(() => {
    // Summary overview: the scene is applied by the Scenes panel (via the
    // summaryScene prop), not the file loader — so load nothing here.
    if (summaryMode) {
      return { type: "none" as const };
    }

    // If a specific job is selected, load that job
    if (selectedJobId !== null) {
      return { type: "job" as const, jobId: selectedJobId };
    }

    if (selectedMemberProjectId !== null) {
      // Find the selected member project
      const memberProject = memberProjects?.find(
        (p) => p.id === selectedMemberProjectId
      );
      // The server names the dataset's current model when it can; that is
      // what the overview's site boxes open, so a plain selection agrees.
      if (memberProject?.current_model_job) {
        return { type: "job" as const, jobId: memberProject.current_model_job.id };
      }
      if (memberProject?.jobs) {
        // Older server: the most recent finished refmac or dimple job
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
  }, [parentFiles, selectedMemberProjectId, selectedJobId, memberProjects, summaryMode]);

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

  // Sites are written one at a time, addressed by id. The whole list used to
  // be PUT back on every change, so two people editing a campaign discarded
  // each other's sites, and a rename orphaned the verdicts recorded there.
  const handleAddSite = useCallback(
    async (site: NewCampaignSite) => {
      if (!campaignId) return;
      await campaignsApi.addSite(campaignId, site);
      mutateSites();
    },
    [campaignId, campaignsApi, mutateSites]
  );

  const handleUpdateSite = useCallback(
    async (siteId: number, changes: Partial<NewCampaignSite>) => {
      if (!campaignId) return;
      await campaignsApi.updateSite(campaignId, siteId, changes);
      mutateSites();
    },
    [campaignId, campaignsApi, mutateSites]
  );

  const handleDeleteSite = useCallback(
    async (siteId: number) => {
      if (!campaignId) return;
      await campaignsApi.deleteSite(campaignId, siteId);
      mutateSites();
      // A deleted site takes its verdicts with it, so what is loaded here is
      // now stale.
      setEvaluations((current) => current.filter((e) => e.site_id !== siteId));
    },
    [campaignId, campaignsApi, mutateSites]
  );

  // What was found at each site in the selected dataset.
  //
  // Fetched rather than read off `memberProjects`: that payload carries only
  // hits and unclears, and a control that records verdicts has to be able to
  // tell "looked, found nothing" from "not looked at yet" — the distinction
  // these rows exist for. The deps deliberately exclude `campaignsApi`, which
  // is a fresh object every render and would re-run this mid-flight.
  const [evaluations, setEvaluations] = useState<SiteEvaluation[]>([]);
  const refreshEvaluations = useCallback(async () => {
    if (!campaignId || !selectedMemberProjectId) {
      setEvaluations([]);
      return;
    }
    try {
      setEvaluations(
        await campaignsApi.fetchEvaluations(campaignId, selectedMemberProjectId)
      );
    } catch (err) {
      console.error("Failed to load site evaluations:", err);
      setEvaluations([]);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [campaignId, selectedMemberProjectId]);

  useEffect(() => {
    let cancelled = false;
    if (!campaignId || !selectedMemberProjectId) {
      setEvaluations([]);
      return;
    }
    campaignsApi
      .fetchEvaluations(campaignId, selectedMemberProjectId)
      .then((rows) => {
        if (!cancelled) setEvaluations(rows);
      })
      .catch((err) => {
        console.error("Failed to load site evaluations:", err);
        if (!cancelled) setEvaluations([]);
      });
    return () => {
      cancelled = true;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [campaignId, selectedMemberProjectId]);

  const handleSetVerdict = useCallback(
    async (siteId: number, verdict: SiteVerdict) => {
      if (!campaignId || !selectedMemberProjectId) return;
      await campaignsApi.setSiteEvaluation(
        campaignId,
        siteId,
        selectedMemberProjectId,
        verdict
      );
      await refreshEvaluations();
    },
    [campaignId, selectedMemberProjectId, campaignsApi, refreshEvaluations]
  );

  const handleClearVerdict = useCallback(
    async (siteId: number) => {
      if (!campaignId || !selectedMemberProjectId) return;
      await campaignsApi.clearSiteEvaluation(
        campaignId,
        siteId,
        selectedMemberProjectId
      );
      await refreshEvaluations();
    },
    [campaignId, selectedMemberProjectId, campaignsApi, refreshEvaluations]
  );

  if (!campaign) {
    return <div>Loading campaign...</div>;
  }

  return (
    <ClientStoreProvider>
      <CampaignMoorhenWrapper
        campaign={campaign}
        fileSource={fileIds}
        summaryScene={summaryMode ? summaryScene : null}
        viewParam={viewParam}
        // Going to one named site; several are framed by their scene, whose
        // camera looks at their centroid, and going to the first would undo it.
        initialSiteId={summarySiteIds.length === 1 ? summarySiteIds[0] : null}
        initialEvent={initialEvent}
        sceneSource={sceneSource}
        sites={sites || []}
        onAddSite={handleAddSite}
        onUpdateSite={handleUpdateSite}
        onDeleteSite={handleDeleteSite}
        evaluations={evaluations}
        onSetVerdict={handleSetVerdict}
        onClearVerdict={handleClearVerdict}
        memberProjects={memberProjects || []}
        selectedMemberProjectId={selectedMemberProjectId}
        onSelectMemberProject={handleSelectMemberProject}
        parentProject={parentProject}
      />
    </ClientStoreProvider>
  );
}

// Outer component with Suspense boundary for useSearchParams
const CampaignPageClient = () => {
  return (
    <Suspense fallback={<div>Loading...</div>}>
      <CampaignPageContent />
    </Suspense>
  );
};

export default CampaignPageClient;
