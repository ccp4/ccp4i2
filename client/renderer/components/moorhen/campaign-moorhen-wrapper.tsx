"use client";

/**
 * The campaign Moorhen viewer: site navigation, member project switching,
 * verdicts, and a summary scene of every discovered fragment.
 *
 * The viewer itself (refs, sizing, store wiring) comes from
 * use-moorhen-viewer-core; loading and scenes from use-moorhen-loaders, the
 * same code the project/job/file viewer runs. What is here is the campaign's
 * own: which member and job are in view and when that content is ready,
 * the sites and the camera moves between them, the evidence for a PanDDA
 * event, the ligand the member's dictionary describes, and its drawer.
 */
import {
  removeMolecule,
  removeMap,
  removeVectorsMatchingIDString,
  setOrigin,
  setQuat,
  setZoom,
  setRequestDrawScene,
} from "moorhen/react-lib";
// @ts-ignore - moorhen 0.23 type may lack .d.ts depending on build
import type { MoorhenPanel } from "moorhen/react-lib";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { moorhen } from "moorhen/types/moorhen";
import { apiText, apiGet, apiPost, apiUpload } from "../../api-fetch";
import { useMoorhenViewState } from "../../hooks/use-moorhen-view-state";
import { siteViewHref } from "../../lib/site-view-link";
import { extractFileIdFromUniqueId, readCameraState } from "../../lib/moorhen-view-state";
import {
  ProjectGroup,
  CampaignSite,
  NewCampaignSite,
  SiteEvaluation,
  SiteVerdict,
  MemberProjectWithSummary,
} from "../../types/campaigns";
import { Project } from "../../types/models";
import { serialiseScene } from "../../lib/scene";
import { preferredJobId, siteViewUrl } from "../../lib/site-verdicts";
import { fetchDictionaryTexts } from "../../lib/moorhen-dictionaries";
import { candidateLigandCodes, placeLigand } from "../../lib/ligand-codes";
import {
  evidenceLoadPlan,
  fetchEventEvidence,
  fileDownloadUrl,
  type EventRef,
} from "../../lib/event-evidence";
import type { MoorhenScene } from "../../types/moorhen-scene";
import { ANNOTATION_VECTOR_PREFIX } from "../../lib/moorhen-scene-resolver";
import { EMPTY_PROJECT_CONTEXT, ProjectContext } from "../../lib/scene-ref-urls";
import type { SceneSource } from "../../lib/scene-source";
import { CampaignMoorhenTabbedPanel } from "./campaign-moorhen-tabbed-panel";
import { MoorhenViewerFrame, useMoorhenViewerCore } from "./use-moorhen-viewer-core";
import { useMoorhenLoaders, type JobLoadInfo } from "./use-moorhen-loaders";

type FileSource =
  | { type: "none" }
  | { type: "files"; fileIds: number[] }
  | { type: "job"; jobId: number };

export interface CampaignMoorhenWrapperProps {
  campaign: ProjectGroup;
  fileSource: FileSource;
  /** Campaign summary scene (built server-side). When set, it's seeded into
   *  the Scenes panel and auto-applied through that panel's own parse/apply
   *  path — the single rendering pathway shared with hand-edited scenes. */
  summaryScene?: MoorhenScene | null;
  viewParam?: string | null;
  /** A site to move to once the scene is up — the `site` URL parameter, which
   *  is how a verdict chip in the campaign overview opens its site. */
  initialSiteId?: number | null;
  /** A PanDDA event whose evidence (event map, autobuilt pose) to load over
   *  the dataset's model once it is in view -- the `event` URL parameter,
   *  which is how a filled box in the campaign's site matrix opens. */
  initialEvent?: EventRef | null;
  /** A scene a job builds on request (a PanDDA report's button): asked for
   *  once the member's model is in view, and applied over it. The member's
   *  coordinate file is passed as ``apo_file_id`` so a scene that would
   *  load an apo model reuses the loaded one and attaches its dictionary
   *  to it. The scene carries its own camera (the event). */
  sceneSource?: SceneSource | null;
  sites: CampaignSite[];
  onAddSite: (site: NewCampaignSite) => Promise<void>;
  onUpdateSite: (
    siteId: number,
    changes: Partial<NewCampaignSite>
  ) => Promise<void>;
  onDeleteSite: (siteId: number) => Promise<void>;
  /** Verdicts recorded for the selected dataset, empties included. */
  evaluations: SiteEvaluation[];
  onSetVerdict: (siteId: number, verdict: SiteVerdict) => Promise<void>;
  onClearVerdict: (siteId: number) => Promise<void>;
  memberProjects: MemberProjectWithSummary[];
  selectedMemberProjectId: number | null;
  onSelectMemberProject: (projectId: number | null) => void;
  parentProject: Project | null | undefined;
}

const CampaignMoorhenWrapper: React.FC<CampaignMoorhenWrapperProps> = ({
  campaign,
  fileSource,
  summaryScene,
  viewParam,
  initialSiteId,
  initialEvent,
  sceneSource,
  sites,
  onAddSite,
  onUpdateSite,
  onDeleteSite,
  evaluations,
  onSetVerdict,
  onClearVerdict,
  memberProjects,
  selectedMemberProjectId,
  onSelectMemberProject,
  parentProject,
}) => {
  const core = useMoorhenViewerCore({ sidePanelKey: "campaignControls" });
  const { setMessage, dispatch, store, cootInitialized, molecules, maps } = core;

  // Representation visibility state (lifted from control panel for URL capture)
  const [visibleRepresentations, setVisibleRepresentations] = useState<string[]>(["CRs"]);

  // Apply representations to all molecules (used when restoring from URL)
  const applyRepresentationsToMolecules = useCallback(
    async (newReps: string[], currentMolecules: moorhen.Molecule[]) => {
      // Default is ["CRs"], which is what molecules load with
      const defaultReps = ["CRs"];

      // Determine what to add/remove compared to default
      const toAdd = newReps.filter((r) => !defaultReps.includes(r));
      const toRemove = defaultReps.filter((r) => !newReps.includes(r));

      for (const mol of currentMolecules) {
        for (const rep of toAdd) {
          try {
            // reps arrive as string[] from URL-restore state; they are the
            // ToggleButton representation vocabulary, so narrow to the union.
            await mol.addRepresentation(rep as moorhen.RepresentationStyles, "/*/*/*/*");
          } catch (err) {
            console.error(`Failed to add ${rep} to ${mol.name}:`, err);
          }
        }
        for (const rep of toRemove) {
          try {
            mol.clearBuffersOfStyle(rep);
          } catch (err) {
            console.error(`Failed to remove ${rep} from ${mol.name}:`, err);
          }
        }
      }
    },
    []
  );

  // View state hook for URL parameter support
  const { getViewUrl, initialRepresentations } = useMoorhenViewState({
    viewParam: viewParam ?? null,
    onViewRestored: () => {},
    representations: visibleRepresentations,
  });

  const hasInitializedReps = useRef(false);

  const loadedFileSource = useRef<FileSource | null>(null);
  /**
   * True once the current `fileSource` has finished loading into Coot.
   *
   * `cootInitialized` is not that signal, and using it as one is what left
   * the `site` URL parameter pointing at the wrong view. Coot reports itself
   * ready long before a job's coordinates have been fetched, and loading them
   * ends in `centreOn()` — so the camera was moved to the site and then had
   * the load drag it back to the molecule's centre a second later. Anything
   * that positions the view waits for this instead.
   */
  const [contentReady, setContentReady] = useState(false);
  /** Guards against an overtaken load declaring a newer one's content ready. */
  const loadToken = useRef(0);

  // Ligand dictionary file ID for 2D structure display (first dict file found)
  const [ligandDictFileId, setLigandDictFileId] = useState<number | null>(null);
  const [ligandName, setLigandName] = useState<string | null>(null);
  // What "Add ligand here" needs: the codes the member's dictionaries define
  // and the file its coordinates came from (the molecule is found by that,
  // not by which molecule is active; see lib/ligand-codes).
  const [ligandCodes, setLigandCodes] = useState<string[]>([]);
  const [memberCoordFileId, setMemberCoordFileId] = useState<number | null>(null);

  // Initialize representations from URL if available (after molecules load)
  useEffect(() => {
    if (!hasInitializedReps.current && initialRepresentations && molecules.length > 0) {
      setVisibleRepresentations(initialRepresentations);
      applyRepresentationsToMolecules(initialRepresentations, molecules);
      hasInitializedReps.current = true;
    }
  }, [initialRepresentations, molecules, applyRepresentationsToMolecules]);

  // The project a scene ref resolves against when it names none of its own:
  // the campaign's parent, where the campaign-level jobs (a PanDDA run) live.
  const projectContextRef = useRef<ProjectContext>(EMPTY_PROJECT_CONTEXT);
  useEffect(() => {
    projectContextRef.current = parentProject
      ? { pk: parentProject.id, uuid: parentProject.uuid ?? null, name: parentProject.name ?? null }
      : EMPTY_PROJECT_CONTEXT;
  }, [parentProject]);

  const {
    fetchMolecule,
    fetchMapFile,
    fetchFile,
    fetchJobFiles,
    handleMapContourLevelChange,
    handleApplyScene,
    handleCaptureScene,
    handlePromoteSceneToPortable,
  } = useMoorhenLoaders({
    core,
    projectContext: projectContextRef,
    captureProject: { projectName: campaign.name },
  });

  // Cleanup all loaded molecules and maps
  const cleanupLoadedContent = useCallback(() => {
    for (const mol of molecules) {
      dispatch(removeMolecule(mol));
      mol.delete();
    }
    for (const map of maps) {
      dispatch(removeMap(map));
      map.delete();
    }
    // The summary scene's site coins go with the content they annotate.
    dispatch(removeVectorsMatchingIDString(ANNOTATION_VECTOR_PREFIX));
    dispatch(setRequestDrawScene(true));
    // Reset ligand info
    setLigandDictFileId(null);
    setLigandName(null);
    // Reset representation state to default
    setVisibleRepresentations(["CRs"]);
    hasInitializedReps.current = false;
  }, [molecules, maps, dispatch]);

  /**
   * The selected member's job: its first dictionary drives the 2D ligand
   * panel and "Add ligand here", and the view centres on it.
   */
  const loadMemberJob = useCallback(
    (jobId: number) =>
      fetchJobFiles(jobId, {
        centre: true,
        onLoaded: ({ dictionaryFiles, dictionaries, coordFile }: JobLoadInfo) => {
          if (dictionaryFiles.length > 0) {
            const firstDict = dictionaryFiles[0];
            setLigandDictFileId(firstDict.id);
            setLigandName(
              firstDict.name?.replace(/\.cif$/i, "") || firstDict.annotation || "Ligand",
            );
            setLigandCodes(candidateLigandCodes(dictionaries.map((d) => d.text)));
          } else {
            setLigandDictFileId(null);
            setLigandName(null);
            setLigandCodes([]);
          }
          setMemberCoordFileId(coordFile ? coordFile.id : null);
        },
      }),
    [fetchJobFiles],
  );

  /** Bring an extra job into the view from the project browser: nothing
   *  about the current member is disturbed and the camera stays where it is. */
  const importJobFiles = useCallback((jobId: number) => fetchJobFiles(jobId), [fetchJobFiles]);

  // Load files when fileSource changes
  useEffect(() => {
    if (!cootInitialized) return;

    const sourceChanged =
      JSON.stringify(loadedFileSource.current) !== JSON.stringify(fileSource);
    if (!sourceChanged) return;

    // Cleanup existing content before loading new files
    if (loadedFileSource.current !== null) {
      cleanupLoadedContent();
    }

    loadedFileSource.current = fileSource;

    // Only the newest load may declare the content ready, so a source that
    // changes twice in quick succession is not called done by the first
    // load settling late.
    const token = ++loadToken.current;
    const finish = () => {
      if (loadToken.current === token) setContentReady(true);
    };

    if (fileSource.type === "files" && fileSource.fileIds.length > 0) {
      setContentReady(false);
      Promise.all(fileSource.fileIds.map((fileId) => fetchFile(fileId))).then(
        finish,
        finish,
      );
    } else if (fileSource.type === "job") {
      setContentReady(false);
      loadMemberJob(fileSource.jobId).then(finish, finish);
    } else {
      // Nothing to load here: the summary scene is applied by the Scenes
      // panel (auto-apply), so it shares one rendering pathway with
      // hand-edited scenes — and it carries its own camera.
      finish();
    }
    // The loaders are stable; this runs when the source changes.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [fileSource, cootInitialized, cleanupLoadedContent]);

  // Navigate to a site
  const handleGoToSite = useCallback(
    (site: CampaignSite) => {
      dispatch(setOrigin(site.origin));
      if (site.quat) {
        dispatch(setQuat(site.quat));
      }
      if (site.zoom) {
        dispatch(setZoom(site.zoom));
      }
      dispatch(setRequestDrawScene(true));
    },
    [dispatch]
  );

  // View a site: a fresh page on the site's scene, the way the campaign
  // overview opens the summary. In a new tab, because the current session
  // holds a dataset someone is looking at, and a site view replaces the
  // loaded molecules with the site's hits.
  const handleViewSite = useCallback(
    (site: CampaignSite) => {
      window.open(
        siteViewHref(campaign.id, [site.id]),
        "_blank"
      );
    },
    [campaign.id]
  );

  // View several sites as one scene -- a site and its subsites, or
  // neighbours -- the same way, in the order the panel lists them.
  const handleViewSites = useCallback(
    (chosen: CampaignSite[]) => {
      window.open(
        siteViewHref(campaign.id, chosen.map((s) => s.id)),
        "_blank"
      );
    },
    [campaign.id]
  );

  // Move to the site named in the URL, once there is something to move
  // around in.
  //
  // Waits for the content, not merely for coot: dispatching an origin before
  // the molecules are loaded puts the camera in the right place and then has
  // the load's own centreOn() reset it underneath us (see `contentReady`).
  // Fires once, so a user who navigates away from the site is not dragged
  // back by a later re-render.
  const appliedInitialSite = useRef(false);
  useEffect(() => {
    if (appliedInitialSite.current) return;
    if (!initialSiteId || !cootInitialized || !contentReady) return;
    if (sites.length === 0) return;
    const site = sites.find((s) => s.id === initialSiteId);
    if (!site) return;
    appliedInitialSite.current = true;
    handleGoToSite(site);
  }, [initialSiteId, cootInitialized, contentReady, sites, handleGoToSite]);

  // Load the evidence for the event named in the URL: its event map at the
  // receipt's display contour, and its autobuilt pose with the dictionary
  // PanDDA built it from, attached to the pose alone. Over the dataset's own
  // model, once that is in view, without moving the camera (the site effect
  // above has put it where the event is). Fires once per page load, and only
  // for the job the page opened on: switching dataset clears the view, and
  // the event belongs to the first one. What the receipt lacks is simply not
  // loaded -- the model is already there to look at.
  const appliedInitialEvent = useRef(false);
  useEffect(() => {
    if (appliedInitialEvent.current) return;
    if (!initialEvent || !cootInitialized || !contentReady) return;
    if (fileSource.type !== "job") return;
    appliedInitialEvent.current = true;
    const token = loadToken.current;
    const stillCurrent = () => loadToken.current === token;
    (async () => {
      const evidence = await fetchEventEvidence(initialEvent);
      if (!evidence || !stillCurrent()) return;
      const plan = evidenceLoadPlan(evidence);
      if (plan.map) {
        await fetchMapFile(fileDownloadUrl(plan.map.fileId), plan.map.name, {
          description: evidence.event_map?.annotation || plan.map.name,
          shortLabel: `Event ${evidence.event_idx}`,
          contourLevel: plan.map.contourLevel,
          colour: plan.map.colour,
          radius: plan.map.radius,
        });
      }
      if (plan.pose && stillCurrent()) {
        const dictionaries = evidence.dictionary
          ? await fetchDictionaryTexts([evidence.dictionary])
          : [];
        await fetchMolecule(fileDownloadUrl(plan.pose.fileId), plan.pose.name, {
          dictionaries,
          centre: false,
        });
      }
      dispatch(setRequestDrawScene(true));
    })().catch((err) => console.warn("[event evidence] could not load:", err));
    // The loaders are stable, and the guard makes this run once.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [initialEvent, cootInitialized, contentReady, fileSource, dispatch]);

  // Apply the scene the URL asked for, over the member's loaded model. Once
  // per page load, after the content is ready (the scene positions the
  // camera, and a load's own centreOn() would otherwise undo it), and only
  // on a job: the scene names the member's model by the file just loaded.
  const appliedSceneSource = useRef(false);
  const sceneSourceKey = sceneSource ? JSON.stringify(sceneSource) : "";
  useEffect(() => {
    if (appliedSceneSource.current) return;
    if (!sceneSourceKey || !cootInitialized || !contentReady) return;
    if (fileSource.type !== "job") return;
    appliedSceneSource.current = true;
    const source = JSON.parse(sceneSourceKey) as SceneSource;
    const token = loadToken.current;
    (async () => {
      try {
        const job = await apiGet(`jobs/${source.jobId}`);
        if (!job?.task_name) throw new Error(`job ${source.jobId} not found`);
        const response: any = await apiPost(`jobs/${source.jobId}/object_method/`, {
          object_path: job.task_name,
          method_name: source.method,
          args: [],
          kwargs: memberCoordFileId != null
            ? { ...source.kwargs, apo_file_id: memberCoordFileId }
            : source.kwargs,
        });
        const result = response?.data?.result;
        if (!result?.success || !result.scene) {
          throw new Error(result?.error || "the job built no scene");
        }
        if (loadToken.current !== token) return; // the member changed meanwhile
        // JSON is YAML: the scene goes in as the text the scenes panel would take.
        await handleApplyScene(JSON.stringify(result.scene));
      } catch (err) {
        setMessage(`Could not open the scene: ${err instanceof Error ? err.message : String(err)}`, "warning");
      }
    })();
    // memberCoordFileId is set by the load that made contentReady true.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sceneSourceKey, cootInitialized, contentReady, fileSource]);

  // Save current view as a site
  const handleSaveCurrentAsSite = useCallback(
    async (name: string) => {
      const camera = readCameraState(store.getState() as moorhen.State);
      const newSite: NewCampaignSite = {
        name,
        origin: Array.from(camera.origin).slice(0, 3) as [number, number, number],
        quat: Array.from(camera.quat).slice(0, 4) as [number, number, number, number],
        zoom: camera.zoom,
      };
      await onAddSite(newSite);
    },
    [store, onAddSite]
  );

  // Delete a site, by its id. Deleting takes that site's verdicts with it.
  const handleDeleteSite = useCallback(
    async (siteId: number) => {
      await onDeleteSite(siteId);
    },
    [onDeleteSite]
  );

  // Rename a site, and optionally move it to the current view.
  //
  // Addressed by id, not by position: a site's verdicts hang off its id, so a
  // rename has to reach the same row rather than replace a list entry.
  const handleUpdateSite = useCallback(
    async (siteId: number, name: string, updatePosition: boolean) => {
      const changes: Partial<NewCampaignSite> = { name };

      if (updatePosition) {
        const camera = readCameraState(store.getState() as moorhen.State);
        changes.origin = Array.from(camera.origin).slice(0, 3) as [number, number, number];
        changes.quat = Array.from(camera.quat).slice(0, 4) as [number, number, number, number];
        changes.zoom = camera.zoom;
      }

      await onUpdateSite(siteId, changes);
    },
    [store, onUpdateSite]
  );

  // Record what was found at a site in the selected dataset.
  //
  // This replaces tagging the project with the site's name. A tag could not
  // say that somebody looked and found nothing -- an untagged project was
  // both "empty" and "not yet looked at" -- and it broke on a rename.
  const handleSetVerdict = useCallback(
    async (siteId: number, verdict: SiteVerdict) => {
      if (!selectedMemberProjectId) {
        setMessage("Please select a member project first");
        return;
      }
      const site = sites.find((s) => s.id === siteId);
      try {
        await onSetVerdict(siteId, verdict);
        setMessage(`Recorded ${verdict} at "${site?.name ?? "site"}"`);
      } catch (err) {
        console.error("Failed to record verdict:", err);
        setMessage("Failed to record verdict");
      }
    },
    [selectedMemberProjectId, sites, onSetVerdict, setMessage]
  );

  // Withdraw a verdict: back to nobody having looked, which is not the same
  // as recording "empty".
  const handleClearVerdict = useCallback(
    async (siteId: number) => {
      if (!selectedMemberProjectId) return;
      const site = sites.find((s) => s.id === siteId);
      try {
        await onClearVerdict(siteId);
        setMessage(`Withdrew the verdict at "${site?.name ?? "site"}"`);
      } catch (err) {
        console.error("Failed to withdraw verdict:", err);
        setMessage("Failed to withdraw verdict");
      }
    },
    [selectedMemberProjectId, sites, onClearVerdict, setMessage]
  );

  // Run servalcat_pipe refinement on a molecule
  const handleRunServalcat = useCallback(
    async (mol: moorhen.Molecule) => {
      // Refinement always operates on a member project (child) - the reflections
      // and coordinates belong to the specific crystal, not the campaign parent
      if (!selectedMemberProjectId) {
        setMessage("No member project selected");
        return;
      }
      const projectId = selectedMemberProjectId;

      const mp = memberProjects.find((p) => p.id === selectedMemberProjectId);
      if (!mp?.uuid) {
        setMessage("Cannot determine project UUID");
        return;
      }
      const projectUuid = mp.uuid;
      const projectDbId = projectUuid.replace(/-/g, "");

      setMessage("Creating servalcat refinement job...");

      try {
        // Step 1: Find observation reflections and the free-R set from this
        // project's top-level jobs.
        //
        // Job outputs (directory 1) are preferred over imports (directory 2),
        // but imports are NOT excluded. A member processed from unmerged data
        // always has an F_SIGF_OUT from aimless; one processed from merged
        // data has none unless the refinement step reindexed, because the
        // observations pass through unchanged -- its only observed file is
        // the imported F_SIGF_IN. Asking for directory=1 alone reported "no
        // reflection data" for every such member.
        type MtzFile = { id: number; job: number; uuid: string; name: string;
                         directory: number; content: number | null;
                         sub_type: number | null };
        const topLevelFiles = async (mimeType: string) =>
          ((await apiGet(
            `files/?type=${mimeType}` +
            `&job__project=${projectId}&job__parent__isnull=true`
          )) as MtzFile[] | null) || [];
        // Outputs before imports, then most recent job first
        const byPreference = (a: MtzFile, b: MtzFile) =>
          a.directory - b.directory || b.job - a.job;

        const obsFiles = await topLevelFiles("application/CCP4-mtz-observed");
        if (obsFiles.length === 0) {
          setMessage("No reflection data found in project");
          return;
        }

        // Among outputs if there are any (else among imports), prefer
        // anomalous data (IPAIR content & 1, FPAIR content & 2) over IMEAN/FMEAN
        const sorted = [...obsFiles].sort(byPreference);
        const candidates = sorted.filter((f) => f.directory === sorted[0].directory);
        const hasAnomalous = (f: MtzFile) => f.content !== null && (f.content & 3) !== 0;
        const reflectionFile = candidates.find(hasAnomalous) || candidates[0];

        // The free-R set. Preferring outputs matters more here than for the
        // observations: the imported set is the campaign's shared one, from
        // the reference crystal, while FREERFLAG_OUT is that set reconciled
        // with this dataset's cell and resolution.
        const freeRFile = (
          await topLevelFiles("application/CCP4-mtz-freerflag")
        ).sort(byPreference)[0];

        // Step 2: Create servalcat_pipe job
        const jobResponse = await apiPost<{
          status: string;
          data: { new_job: { id: number; uuid: string } };
        }>(`projects/${projectId}/create_task/`, {
          task_name: "servalcat_pipe",
          title: `Servalcat refinement of ${mol.name}`,
        });
        const newJobId = jobResponse.data.new_job.id;

        // Step 3: Upload coordinates from Moorhen (preserve original format)
        setMessage("Uploading coordinates...");
        const coordText = await mol.getAtoms();
        const ismmCIF = (mol as any).coordsFormat === "mmcif";
        const mimeType = ismmCIF ? "chemical/x-cif" : "chemical/x-pdb";
        const ext = ismmCIF ? ".cif" : ".pdb";
        const coordBlob = new Blob([coordText], { type: mimeType });
        const coordFile = new File([coordBlob], `${mol.name || "coords"}${ext}`);
        const coordFormData = new FormData();
        coordFormData.append("file", coordFile);
        coordFormData.append("object_path", "servalcat_pipe.inputData.XYZIN");
        await apiUpload(`jobs/${newJobId}/upload_file_param/`, coordFormData);

        // Step 4: Link reflections via database file reference
        setMessage("Linking reflection data...");
        const reflDbFileId = reflectionFile.uuid.replace(/-/g, "");
        const hklinValue: Record<string, string | number> = {
          project: projectDbId,
          dbFileId: reflDbFileId,
        };
        if (reflectionFile.sub_type !== null) {
          hklinValue.subType = reflectionFile.sub_type;
        }
        await apiPost(`jobs/${newJobId}/set_parameter/`, {
          object_path: "servalcat_pipe.inputData.HKLIN",
          value: hklinValue,
        });

        // Without this the refinement ran with no free set at all, so its
        // R-free meant nothing and could not be compared between siblings.
        if (freeRFile) {
          await apiPost(`jobs/${newJobId}/set_parameter/`, {
            object_path: "servalcat_pipe.inputData.FREERFLAG",
            value: {
              project: projectDbId,
              dbFileId: freeRFile.uuid.replace(/-/g, ""),
            },
          });
        }

        // Step 5: Upload dictionary if available
        // Uses upload_file_param (not set_parameter) because DICT_LIST starts empty
        // and upload_file_param handles list expansion automatically
        if (ligandDictFileId) {
          setMessage("Uploading ligand dictionary...");
          const dictUrl = `/api/proxy/ccp4i2/files/${ligandDictFileId}/download/`;
          const dictText = await apiText(dictUrl);
          const dictBlob = new Blob([dictText], { type: "application/refmac-dictionary" });
          const dictFile = new File([dictBlob], "ligand.cif");
          const dictFormData = new FormData();
          dictFormData.append("file", dictFile);
          dictFormData.append("object_path", "servalcat_pipe.inputData.DICT_LIST[0]");
          await apiUpload(`jobs/${newJobId}/upload_file_param/`, dictFormData);
        }

        // Step 6: Run the job
        setMessage("Running servalcat refinement...");
        await apiPost(`jobs/${newJobId}/run/`, {});

        setMessage("Servalcat refinement job submitted successfully");
      } catch (error) {
        console.error("Failed to create servalcat job:", error);
        setMessage(
          `Failed to create servalcat job: ${error instanceof Error ? error.message : "Unknown error"}`
        );
      }
    },
    [selectedMemberProjectId, memberProjects, ligandDictFileId, setMessage]
  );

  // Place the member's ligand at the view centre. Nothing is saved: the
  // result lives in the browser until it is pushed, and it is not pushed
  // automatically because a placed, unfitted ligand in an arbitrary
  // orientation is not something anyone wants silently written into their
  // project.
  const handleAddLigand = useCallback(
    async (code: string) => {
      try {
        await placeLigand(molecules, memberCoordFileId, code);
        dispatch(setRequestDrawScene(true));
        setMessage(`Ligand ${code} added. Push to CCP4i2 to keep it.`, "success");
      } catch (err) {
        setMessage(
          `Could not add ${code}: ${err instanceof Error ? err.message : String(err)}`,
          "error",
        );
      }
    },
    [molecules, memberCoordFileId, dispatch, setMessage]
  );

  // Every loaded molecule's source file, by the fileId its uniqueId carries.
  // The scene is where a molecule's provenance lives: the loader names it
  // after the scene's file ref, and that ref says which project it came from.
  const projectUuidByFileId = useMemo(() => {
    const byFileId = new Map<number, string>();
    for (const ref of summaryScene?.files ?? []) {
      if (ref.kind === "coordinates" && ref.fileId !== undefined && ref.projectId) {
        byFileId.set(ref.fileId, ref.projectId);
      }
    }
    return byFileId;
  }, [summaryScene]);

  // The way out of a summary scene: open the dataset a drawn ligand belongs
  // to, where its maps are, and at this site when the scene is a site's. A
  // summary scene deliberately carries no density (campaign_scene emits
  // coordinates and dictionaries, never maps), so nothing drawn in it can be
  // judged on the evidence — this is the control that leads to the evidence.
  //
  // Returns a reason rather than nothing when there is no dataset to open:
  // the campaign's own reference structure is drawn in every summary scene
  // and is not a dataset. That is worth saying once, where a dead control
  // would only be puzzled over.
  const datasetLink = useCallback(
    (mol: moorhen.Molecule): { url: string; label: string } | { reason: string } => {
      const fileId = extractFileIdFromUniqueId(mol.uniqueId || "");
      const projectUuid = fileId !== null ? projectUuidByFileId.get(fileId) : undefined;
      if (!projectUuid) {
        return { reason: "Not a campaign dataset (this is the reference structure)" };
      }
      const project = memberProjects.find((p) => p.uuid === projectUuid);
      if (!project) {
        return { reason: "That dataset is no longer a member of this campaign" };
      }
      const jobId = preferredJobId(project.jobs);
      if (!jobId) return { reason: `${project.name} has no job to open` };
      // With a site, the same link a verdict chip in the campaign overview
      // follows; without one (the whole-campaign summary) the dataset opens
      // on its own job, camera wherever the job leaves it.
      const url =
        initialSiteId != null
          ? siteViewUrl(campaign.id, jobId, initialSiteId)
          : `/ccp4i2/moorhen-page/campaign/${campaign.id}?job=${jobId}`;
      if (!url) return { reason: `${project.name} has no job to open` };
      return { url, label: project.name };
    },
    [projectUuidByFileId, memberProjects, campaign.id, initialSiteId],
  );

  // When viewing the campaign summary, serialise the scene to YAML so the
  // Scenes panel can show it in the editor and auto-apply it through its own
  // parse/apply path (one shared rendering pathway).
  const summarySceneYaml = useMemo(() => {
    if (!summaryScene) return undefined;
    try {
      return serialiseScene(summaryScene);
    } catch (err) {
      console.warn("[scene] failed to serialise summary scene for editor:", err);
      return undefined;
    }
  }, [summaryScene]);

  // The drawer: campaign Controls + Scenes editor under one Moorhen
  // side-panel registration (matching the job/file viewers' tabbed panel).
  const extraSidePanels: Record<string, MoorhenPanel> = {
    campaignControls: {
      icon: "MatSymSettings",
      label: "Campaign",
      panelContent: (
        <CampaignMoorhenTabbedPanel
          controlPanelProps={{
            campaign,
            sites,
            onGoToSite: handleGoToSite,
            onViewSite: handleViewSite,
            onViewSites: handleViewSites,
            onSaveCurrentAsSite: handleSaveCurrentAsSite,
            onUpdateSite: handleUpdateSite,
            onDeleteSite: handleDeleteSite,
            memberProjects,
            selectedMemberProjectId,
            onSelectMemberProject,
            parentProject,
            getViewUrl,
            molecules,
            visibleRepresentations,
            onRepresentationsChange: setVisibleRepresentations,
            ligandDictFileId,
            ligandName,
            ligandCodes,
            memberCoordFileId,
            onAddLigand: handleAddLigand,
            maps,
            onMapContourLevelChange: handleMapContourLevelChange,
            evaluations,
            onSetVerdict: handleSetVerdict,
            onClearVerdict: handleClearVerdict,
            onFileSelect: fetchFile,
            onJobLoad: importJobFiles,
            onRunServalcat: handleRunServalcat,
            // A summary scene is the view whose molecules have no density
            // behind them and have been moved by a fit — whether it is one
            // site's or the whole campaign's.
            siteSummary: !!summaryScene,
            datasetLink: summaryScene ? datasetLink : undefined,
          }}
          onApplyScene={handleApplyScene}
          onCaptureScene={handleCaptureScene}
          onPromoteSceneToPortable={handlePromoteSceneToPortable}
          cootInitialized={cootInitialized}
          initialSceneYaml={summarySceneYaml}
          autoApplyInitialScene={!!summarySceneYaml}
        />
      ),
    },
  };

  return <MoorhenViewerFrame core={core} extraSidePanels={extraSidePanels} />;
};

export default CampaignMoorhenWrapper;
