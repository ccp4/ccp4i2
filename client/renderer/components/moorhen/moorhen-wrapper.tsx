"use client";

/**
 * The project/job/file/session Moorhen viewer.
 *
 * The viewer itself (refs, sizing, store wiring) comes from
 * use-moorhen-viewer-core; loading and scenes from use-moorhen-loaders, the
 * same code the campaign viewer runs. What is here is this page's own: which
 * project it is scoped to and how that is found, what the URL asked it to
 * open (files, a job, paths, a scene a job builds, a recorded session), the
 * natural-language scene authoring, and its drawer.
 */
import { setRequestDrawScene } from "moorhen/react-lib";
// @ts-ignore - moorhen 0.23 type may lack .d.ts depending on build
import type { MoorhenPanel } from "moorhen/react-lib";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { moorhen } from "moorhen/types/moorhen";
import { MoorhenCcp4i2TabbedPanel } from "./moorhen-ccp4i2-tabbed-panel";
import { apiGet, apiText, apiPost, apiUpload } from "../../api-fetch";
import { useMoorhenViewState } from "../../hooks/use-moorhen-view-state";
import { useMoorhenSession } from "../../hooks/use-moorhen-session";
import {
  COORDINATE_TYPES,
  DICTIONARY_TYPE,
  fetchJobDictionaryFiles,
  type DictionaryFile,
  type DictionaryToAttach,
} from "../../lib/moorhen-dictionaries";
import {
  baseName,
  classifyCif,
  isLoadableFile,
  looksLikeMask,
  planPathLoad,
} from "../../lib/moorhen-path-load";
import { projectPathUrl } from "../../lib/scene-project-path";
import { SceneSource } from "../../lib/scene-source";
import { EMPTY_PROJECT_CONTEXT, ProjectContext } from "../../lib/scene-ref-urls";
import { extractFileIdFromUniqueId } from "../../lib/moorhen-view-state";
import {
  buildContentsBlock,
  buildManifestBlock,
  buildAuthoringPrompt,
  buildGroundingBlock,
  extractPdbIds,
  fetchPdbContents,
} from "../../lib/moorhen-scene-prompt";
import { useSceneNlCapability, generateScene } from "./use-scene-nl-capability";
import { MoorhenViewerFrame, useMoorhenViewerCore } from "./use-moorhen-viewer-core";
import { useMoorhenLoaders } from "./use-moorhen-loaders";

export interface MoorhenWrapperProps {
  fileIds?: number[];
  viewParam?: string | null;
  jobId?: number | null;
  /** Project pk for a project-scoped Moorhen page: provides project context
   *  (manifest + job/param resolution) without loading any specific file/job. */
  projectId?: number | null;
  /** The job this window is a recorded Moorhen session for: load its inputs
   *  from the server's load plan, save back into it, finish it. */
  sessionJobId?: number | null;
  /** Files or folders to open by their path in the project (needs
   *  ``projectId``): what a job directory holds, registered or not. */
  filePaths?: string[];
  /** A scene the job builds on request (a report link): called once Coot is
   *  ready, and what it returns applied. */
  sceneSource?: SceneSource | null;
}

const MoorhenWrapper: React.FC<MoorhenWrapperProps> = ({ fileIds, viewParam, jobId, projectId, sessionJobId, filePaths, sceneSource }) => {
  const core = useMoorhenViewerCore({ sidePanelKey: "ccp4i2Controls" });
  const { setMessage, dispatch, cootInitialized, molecules, maps } = core;
  const session = useMoorhenSession(sessionJobId, !!sessionJobId);

  // View state hook for URL parameter support
  const { getViewUrl } = useMoorhenViewState({
    viewParam: viewParam ?? null,
    onViewRestored: () => {},
  });

  // The project this page is scoped to, mirrored into a ref so the
  // (stable-identity) scene fetchers resolve job+param and projectPath refs
  // without re-binding. Resolved below: an explicit projectId (the
  // project-scoped page); else the job's; else the first file's job's,
  // covering the file viewer. Null only on the truly contextless blank page.
  const projectContextRef = useRef<ProjectContext>(EMPTY_PROJECT_CONTEXT);
  const [projectInfo, setProjectInfo] = useState<{ id: string; name?: string } | null>(null);
  // The project's primary key, as soon as it is known. Exposed as state, not
  // just the ref, so the CCP4i2 browser can open on this project.
  const [projectPk, setProjectPk] = useState<number | null>(projectId ?? null);
  const firstFileId = fileIds?.[0];
  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        let projPk: number | null = projectId ?? null;
        if (projPk == null) {
          let jobPk: number | null = jobId ?? null;
          if (jobPk == null && firstFileId != null) {
            const fileInfo = await apiGet(`files/${firstFileId}`);
            if (cancelled) return;
            jobPk = fileInfo?.job ?? null;
          }
          if (jobPk == null) return;
          const jobInfo = await apiGet(`jobs/${jobPk}`);
          if (cancelled || !jobInfo?.project) return;
          projPk = jobInfo.project;
        }
        if (projPk == null) return;
        projectContextRef.current = { ...projectContextRef.current, pk: projPk };
        setProjectPk(projPk);
        const proj = await apiGet(`projects/${projPk}`);
        if (cancelled || !proj?.uuid) return;
        projectContextRef.current = { pk: projPk, uuid: proj.uuid, name: proj.name ?? null };
        setProjectInfo({ id: proj.uuid, name: proj.name });
      } catch (err) {
        console.warn("[wrapper] failed to resolve project for scene authoring:", err);
      }
    })();
    return () => { cancelled = true; };
  }, [jobId, firstFileId, projectId]);

  const {
    loadStructure,
    fetchMapFile,
    loadMoleculesFromDictionary,
    loadDictionaryOnlyJob,
    fetchFile,
    fetchJobFiles,
    handleMapContourLevelChange,
    handleApplyScene,
    applySceneRef,
    handleCaptureScene,
    handlePromoteSceneToPortable,
  } = useMoorhenLoaders({
    core,
    projectContext: projectContextRef,
    captureProject: projectInfo ? { projectId: projectInfo.id, projectName: projectInfo.name } : undefined,
  });

  // Run servalcat refinement using the current job as context
  const [servalcatStatus, setServalcatStatus] = useState<string | null>(null);
  const handleRunServalcat = useCallback(
    async (mol: moorhen.Molecule) => {
      if (!jobId) {
        setServalcatStatus("No source job — open from a job to use refinement");
        return;
      }
      try {
        setServalcatStatus("Creating servalcat job...");

        // Look up the job's UUID and project
        const jobInfo = await apiGet(`jobs/${jobId}`);
        if (!jobInfo?.uuid || !jobInfo?.project) {
          setServalcatStatus("Could not retrieve job information");
          return;
        }

        // Create servalcat_pipe job with context_job_uuid for auto-population
        const createResult = await apiPost(
          `projects/${jobInfo.project}/create_task/`,
          {
            task_name: "servalcat_pipe",
            title: `Servalcat refinement of ${mol.name}`,
            context_job_uuid: jobInfo.uuid,
          }
        );
        const newJobId = createResult?.data?.new_job?.id;
        if (!newJobId) {
          setServalcatStatus("Failed to create servalcat job");
          return;
        }

        // Upload Moorhen coordinates as XYZIN (preserve original format)
        setServalcatStatus("Uploading coordinates...");
        const coordText = await mol.getAtoms();
        const ismmCIF = (mol as any).coordsFormat === "mmcif";
        const mimeType = ismmCIF ? "chemical/x-cif" : "chemical/x-pdb";
        const ext = ismmCIF ? ".cif" : ".pdb";
        const coordBlob = new Blob([coordText], { type: mimeType });
        const coordFile = new File([coordBlob], `${mol.name || "coords"}${ext}`);
        const formData = new FormData();
        formData.append("file", coordFile);
        formData.append("object_path", "servalcat_pipe.inputData.XYZIN");
        await apiUpload(`jobs/${newJobId}/upload_file_param/`, formData);

        // Run the job
        setServalcatStatus("Running servalcat refinement...");
        await apiPost(`jobs/${newJobId}/run/`, {});

        setServalcatStatus("Servalcat refinement submitted");
        // Clear status after a few seconds
        setTimeout(() => setServalcatStatus(null), 4000);
      } catch (err) {
        console.error("Failed to run servalcat:", err);
        setServalcatStatus(
          `Failed: ${err instanceof Error ? err.message : "Unknown error"}`
        );
      }
    },
    [jobId]
  );

  // Assemble the per-call grounding for scene authoring: project identity +
  // referenceable-file manifest + loaded-structure (and named-PDB) contents.
  // Shared by both the copy-paste prompt and the integrated Generate path so
  // they ground the model identically. `request` is used only to spot PDB ids
  // to digest.
  const assembleGrounding = useCallback(async (request: string): Promise<{
    project?: { name?: string; id?: string };
    manifest: string;
    contents: string;
  }> => {
    const mol = molecules[0];
    const fid = mol ? extractFileIdFromUniqueId(mol.uniqueId || "") : null;

    // Project identity + manifest. Resolve the project robustly (prefer the
    // cached projectInfo, else fetch) so the prompt always pins the current
    // project, and keep the files list so we can label the loaded structure by
    // its annotation below.
    let manifestBlock = "(no project context)";
    let project: { name?: string; id?: string } | undefined =
      projectInfo ? { name: projectInfo.name, id: projectInfo.id } : undefined;
    let projectFiles: { id?: number; annotation?: string }[] = [];
    const pk = projectContextRef.current.pk;
    if (pk != null) {
      try {
        if (!project?.id) {
          const proj = await apiGet(`projects/${pk}`);
          if (proj?.uuid) project = { name: proj.name, id: proj.uuid };
        }
        const [jobs, files] = await Promise.all([
          apiGet(`jobs/?project=${pk}`),
          apiGet(`files/?project=${pk}`),
        ]);
        if (Array.isArray(jobs) && Array.isArray(files)) {
          projectFiles = files;
          manifestBlock = buildManifestBlock(jobs, files);
        }
      } catch (err) {
        console.warn("[scene-prompt] manifest failed:", err);
      }
    }

    // Contents summary. One block per source: the loaded structure (labelled by
    // the file's annotation where we have one, e.g. "CDK2-Cyclin A"), plus a
    // PDBe digest for any PDB id named in the request (so the model authoring a
    // `pdb:` ref gets accurate chains + ligand CIDs — it can't fetch them itself).
    const contentsParts: string[] = [];
    if (mol && fid != null) {
      try {
        const resp = await apiGet(`files/${fid}/digest/`);
        const digest = resp?.data ?? resp;
        const annotation = projectFiles.find((f) => f.id === fid)?.annotation;
        const label = annotation || mol.name || `file_${fid}`;
        contentsParts.push(buildContentsBlock(digest, label));
      } catch (err) {
        console.warn("[scene-prompt] digest failed:", err);
      }
    }
    try {
      const pdbBlocks = await Promise.all(
        extractPdbIds(request).map((id) => fetchPdbContents(id)),
      );
      for (const block of pdbBlocks) if (block) contentsParts.push(block);
    } catch (err) {
      console.warn("[scene-prompt] PDB digest failed:", err);
    }
    const contentsBlock = contentsParts.length
      ? contentsParts.join("\n\n")
      : "(no structure loaded)";

    return { project, manifest: manifestBlock, contents: contentsBlock };
  }, [molecules, projectInfo]);

  // Copy-paste path: the full self-contained prompt (system prompt + grounding +
  // request) for pasting into any chatbot.
  const handleBuildAuthoringPrompt = useCallback(async (request: string): Promise<string> => {
    const { project, manifest, contents } = await assembleGrounding(request);
    return buildAuthoringPrompt({ project, contents, manifest, request });
  }, [assembleGrounding]);

  // Integrated path: send request + grounding to the nlp_scene endpoint (which
  // holds the static system prompt server-side) and return the raw scene text.
  const handleGenerateScene = useCallback(async (request: string): Promise<string> => {
    const { project, manifest, contents } = await assembleGrounding(request);
    const grounding = buildGroundingBlock({ project, manifest, contents });
    return generateScene(request, grounding);
  }, [assembleGrounding]);

  const sceneNl = useSceneNlCapability();

  // The drawer: CCP4i2 controls + Scenes under one Moorhen side-panel
  // registration whose content is our own MUI Tabs container. (Registering
  // two extraSidePanels surfaces only one tab in Moorhen's switcher.)
  const extraSidePanels: Record<string, MoorhenPanel> = useMemo(() => ({
    ccp4i2Controls: {
      icon: "MatSymSettings",
      label: "CCP4i2",
      panelContent: (
        <MoorhenCcp4i2TabbedPanel
          onFileSelect={fetchFile}
          onJobLoad={fetchJobFiles}
          initialProjectPk={projectPk}
          getViewUrl={getViewUrl}
          molecules={molecules}
          maps={maps}
          onMapContourLevelChange={handleMapContourLevelChange}
          onRunServalcat={jobId ? handleRunServalcat : undefined}
          servalcatStatus={servalcatStatus}
          onApplyScene={handleApplyScene}
          onCaptureScene={handleCaptureScene}
          onPromoteSceneToPortable={handlePromoteSceneToPortable}
          onBuildAuthoringPrompt={handleBuildAuthoringPrompt}
          onGenerateScene={sceneNl.available ? handleGenerateScene : undefined}
          cootInitialized={cootInitialized}
          session={session}
        />
      ),
    },
  }), [fetchFile, fetchJobFiles, projectPk, getViewUrl, molecules, maps, handleMapContourLevelChange, jobId, handleRunServalcat, servalcatStatus, handleApplyScene, handleCaptureScene, handlePromoteSceneToPortable, handleBuildAuthoringPrompt, sceneNl.available, handleGenerateScene, cootInitialized, session]);

  useEffect(() => {
    if (!fileIds || fileIds.length === 0 || !cootInitialized) return;
    (async () => {
      // Know what each file is before loading any of them: a dictionary that
      // arrives beside coordinates belongs to those coordinates and must not
      // go into Coot's global store, which is what loading it by itself does.
      const infos = await Promise.all(
        fileIds.map((id) => apiGet(`files/${id}`).catch(() => null)),
      );
      const hasCoordinates = infos.some((info) => info && COORDINATE_TYPES.has(info.type));
      // Opened on a job: the job's dictionaries (own and inputs) are the ones.
      const jobDictionaries =
        jobId && hasCoordinates ? await fetchJobDictionaryFiles(jobId) : undefined;
      // Opened on a job with dictionaries and no coordinates: the monomers
      // the dictionaries describe are the job's result, so draw those.
      const dictionaryOnlyJob = Boolean(jobId) && !hasCoordinates;
      for (const [index, fileId] of fileIds.entries()) {
        const info = infos[index];
        if ((hasCoordinates || dictionaryOnlyJob) && info?.type === DICTIONARY_TYPE) continue;
        await fetchFile(fileId, { dictionaryFiles: jobDictionaries });
      }
      if (dictionaryOnlyJob) {
        const ownDictionaries = infos.filter((info) => info?.type === DICTIONARY_TYPE);
        if (ownDictionaries.length > 0) await loadDictionaryOnlyJob(jobId as number, ownDictionaries);
      }
    })();
    // The loaders are stable for the page's life; this runs once per load request.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [fileIds, cootInitialized]);

  // Paths in the project (from a directory browser's "Open in Moorhen"):
  // files load as they are, folders as the files directly in them.
  const filePathsKey = (filePaths ?? []).join("\n");
  useEffect(() => {
    if (!filePathsKey || !projectId || !cootInitialized) return;
    const paths = filePathsKey.split("\n");
    const urlFor = (path: string) => projectPathUrl(projectId, path);
    (async () => {
      const files: string[] = [];
      for (const path of paths) {
        if (isLoadableFile(path)) {
          files.push(path);
          continue;
        }
        const listing = await apiGet(
          `projects/${projectId}/directory?path=${encodeURIComponent(path)}`,
        ).catch(() => null);
        for (const entry of listing?.container ?? []) {
          if (entry.type === "file") files.push(path ? `${path}/${entry.name}` : entry.name);
        }
      }
      const plan = planPathLoad(files);
      const dictionaries: DictionaryToAttach[] = [];
      const models = [...plan.coordinates];
      for (const path of plan.cifs) {
        const text = await apiText(urlFor(path)).catch(() => null);
        const kind = text ? classifyCif(text) : null;
        if (kind === "coordinates") models.push(path);
        else if (kind === "dictionary") dictionaries.push({ text: text as string });
      }
      const failed: string[] = [];
      let loaded = 0;
      for (const [index, path] of models.entries()) {
        const molecule = await loadStructure(urlFor(path), baseName(path), {
          dictionaries,
          centre: index === 0,
        });
        if (molecule) loaded += 1;
        else failed.push(baseName(path));
      }
      if (models.length === 0) {
        for (const [index, dict] of dictionaries.entries()) {
          const built = await loadMoleculesFromDictionary(dict, `path-dict-${index}`, {
            centre: index === 0,
          });
          loaded += built.length;
        }
        if (dictionaries.length > 0) dispatch(setRequestDrawScene(true));
      }
      for (const path of plan.maps) {
        try {
          const map = await fetchMapFile(urlFor(path), baseName(path), {
            isMask: looksLikeMask(path),
          });
          if (map) loaded += 1;
          else failed.push(baseName(path));
        } catch (err) {
          console.warn(`[paths] could not open ${path}:`, err);
          failed.push(baseName(path));
        }
      }
      if (loaded === 0 && failed.length === 0) {
        setMessage(`Nothing in ${paths.join(", ")} that Moorhen can open`, "warning");
      } else if (failed.length > 0) {
        setMessage(`Could not open ${failed.join(", ")}`, "warning");
      }
    })();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [filePathsKey, projectId, cootInitialized]);

  // A scene a job builds on request: ask the job's plugin, apply the answer.
  const sceneSourceKey = sceneSource ? JSON.stringify(sceneSource) : "";
  useEffect(() => {
    if (!sceneSourceKey || !cootInitialized) return;
    const source = JSON.parse(sceneSourceKey) as SceneSource;
    (async () => {
      try {
        const job = await apiGet(`jobs/${source.jobId}`);
        if (!job?.task_name) throw new Error(`job ${source.jobId} not found`);
        const response: any = await apiPost(`jobs/${source.jobId}/object_method/`, {
          object_path: job.task_name,
          method_name: source.method,
          args: [],
          kwargs: source.kwargs,
        });
        const result = response?.data?.result;
        if (!result?.success || !result.scene) {
          throw new Error(result?.error || "the job built no scene");
        }
        if (!applySceneRef.current) throw new Error("the viewer is not ready to apply a scene");
        // JSON is YAML: the scene goes in as the text the scenes panel would take.
        await applySceneRef.current(JSON.stringify(result.scene));
      } catch (err) {
        setMessage(`Could not open the scene: ${err instanceof Error ? err.message : String(err)}`, "warning");
      }
    })();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [sceneSourceKey, cootInitialized]);

  // A session window loads the job's inputs from the server's load plan,
  // once, in the server's order (dictionaries first, then coordinates,
  // then maps), awaiting each so a model finds its ligand geometry.
  const sessionPlanAppliedRef = useRef(false);
  useEffect(() => {
    if (!session?.state || !cootInitialized || sessionPlanAppliedRef.current) return;
    sessionPlanAppliedRef.current = true;
    const plan = session.state.load_plan;
    // The session job's dictionary inputs belong to its coordinate inputs.
    // (The job has not been dispatched, so the database does not list its
    // inputs yet; the plan does.) With none given, each coordinate file
    // brings the dictionaries of the job that produced it.
    const sessionDictionaries: DictionaryFile[] = plan
      .filter((entry) => entry.kind === "dictionary" && entry.file_id != null)
      .map((entry) => ({ id: entry.file_id as number, name: entry.name }));
    (async () => {
      for (const entry of plan) {
        if (entry.file_id == null) {
          console.warn(`Session load plan: ${entry.label} has no file record yet; skipped`);
          continue;
        }
        if (entry.kind === "dictionary") continue; // attached to the coordinates below
        try {
          await fetchFile(
            entry.file_id,
            sessionDictionaries.length > 0 ? { dictionaryFiles: sessionDictionaries } : {},
          );
        } catch (err) {
          console.warn(`Session load plan: could not load ${entry.label}`, err);
        }
      }
    })();
  }, [session?.state, cootInitialized, fetchFile]);

  return <MoorhenViewerFrame core={core} extraSidePanels={extraSidePanels} />;
};

export default MoorhenWrapper;
