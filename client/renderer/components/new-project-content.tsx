"use client";
import { ChangeEvent, useEffect, useMemo, useState } from "react";
import path from "path";
import { useRouter } from "next/navigation";
import {
  Button,
  CircularProgress,
  Container,
  Stack,
  TextField,
  ToggleButton,
  ToggleButtonGroup,
  Tooltip,
  Typography,
} from "@mui/material";
import { Folder } from "@mui/icons-material";
import { useApi } from "../api";
import { apiGet, apiPost } from "../api-fetch";
import { Project } from "../types/models";
import EditTags from "./edit-tags";
import {
  DroppedFile,
  FileDropZone,
  TASK_FOR_TYPE,
  PARAM_FOR_TYPE,
  AUTO_RUN_FOR_TYPE,
} from "./file-drop-zone";

export const NewProjectContent: React.FC = () => {
  const api = useApi();
  const router = useRouter();
  const [name, setName] = useState("");
  const [description, setDescription] = useState("");
  const [customDirectory, setCustomDirectory] = useState(false);
  // Kept apart deliberately. One variable serving both meant that picking a
  // custom parent and then switching back to "default" left the custom path on
  // screen while the request actually said __default__ — the dialog showing one
  // location and creating the project in another. Held separately, the toggle
  // itself is the way back, with nothing to reset.
  const [proposedParent, setProposedParent] = useState<string>("");
  const [configuredParent, setConfiguredParent] = useState<string>("");
  // Set when the user overrides the proposal from this page — either by
  // browsing, or by asking for the configured root back. null means "whatever
  // the server proposes", which is what lets the request stay __default__.
  const [proposedOverride, setProposedOverride] = useState<string | null>(null);
  const [customParent, setCustomParent] = useState<string>("");

  const effectiveProposal = proposedOverride ?? proposedParent;
  const parentDirectory = customDirectory ? customParent : effectiveProposal;
  // Only worth offering when it would actually change something.
  const canUseConfiguredRoot =
    configuredParent.length > 0 && effectiveProposal !== configuredParent;
  const [directoryExists, setDirectoryExists] = useState(true);
  const [electronAPIAvailable, setElectronAPIAvailable] =
    useState<boolean>(false);
  const [tags, setTags] = useState<number[]>([]);
  const [droppedFiles, setDroppedFiles] = useState<DroppedFile[]>([]);
  const [isCreating, setIsCreating] = useState(false);
  const { data: projects } = api.get<Project[]>("projects");

  // Where the server will actually put a project created with no explicit
  // directory: the parent of the most recently created project, else the
  // configured projects directory. Asked rather than computed — reimplementing
  // that rule here would give two answers to one question.
  useEffect(() => {
    let cancelled = false;
    apiGet<{ data?: { directory?: string; configured?: string } }>(
      "config/default-project-parent/"
    )
      .then((resp) => {
        if (cancelled) return;
        if (resp?.data?.directory) setProposedParent(resp.data.directory);
        if (resp?.data?.configured) setConfiguredParent(resp.data.configured);
      })
      .catch(() => {
        /* older backend: the Electron config value below still applies */
      });
    return () => {
      cancelled = true;
    };
  }, []);

  useEffect(() => {
    // Send a message to the main process to get the config
    if (window.electronAPI) {
      setElectronAPIAvailable(true);
      window.electronAPI.sendMessage("get-config");
      // Listen for messages from the main process
      window.electronAPI.onMessage(
        "message-from-main",
        (event: any, data: any) => {
          if (data.message === "get-config") {
            // Only a fallback for the split second before the server answers,
            // and for a desktop build talking to an older backend.
            setProposedParent((current) =>
              current || data.config.CCP4I2_PROJECTS_DIR
            );
          }
          if (data.message === "check-file-exists") {
            setDirectoryExists(data.exists);
          }
          if (data.message === "choose-project-parent-directory") {
            // Local to this page only — nothing is persisted, and the default
            // above is left untouched so the toggle still reverts to it.
            setCustomParent(data.directory);
          }
        }
      );
    } else {
      setDirectoryExists(false); // Assume directory does not exist in web mode
    }
  }, []);

  const directory = useMemo(() => {
    const result = path.join(parentDirectory || "", name.toLocaleLowerCase());
    if (typeof window !== "undefined" && window.electronAPI) {
      window.electronAPI.sendMessage("check-file-exists", { path: result });
    }
    return result;
  }, [parentDirectory, name]);

  async function createProject() {
    setIsCreating(true);
    try {
      const formData = new FormData();
      formData.append("name", name);
      formData.append("description", description);
      // __default__ only while the server's proposal stands untouched; once the
      // user has overridden it, say so explicitly rather than letting the server
      // re-derive something different.
      formData.append(
        "directory",
        customDirectory || proposedOverride !== null ? directory : "__default__"
      );
      const project = await api.post<Project>("projects", formData);

      // Apply tags to the new project
      for (const tagId of tags) {
        try {
          await apiPost(`projects/${project.id}/tags/`, { tag_id: tagId });
        } catch (err) {
          console.error(`Error applying tag ${tagId}:`, err);
        }
      }

      // If files were dropped, create import jobs sequentially.
      // Sequence files are handled collectively: most consumers want a
      // CAsuDataFile rather than loose CSeqDataFiles, so every dropped
      // sequence accumulates into ONE "Define AU contents" job instead of a
      // ProvideSequence job each.
      const importableFiles = droppedFiles.filter(
        (df) => TASK_FOR_TYPE[df.detectedType] !== null
      );
      const sequenceFiles = importableFiles.filter(
        (df) => df.detectedType === "sequence"
      );
      const otherFiles = importableFiles.filter(
        (df) => df.detectedType !== "sequence"
      );

      for (const df of otherFiles) {
        try {
          await createImportJob(project.id, df);
          // Small delay between jobs to avoid DB contention (SQLite)
          if (importableFiles.length > 1) {
            await new Promise((r) => setTimeout(r, 500));
          }
        } catch (err) {
          console.error(`Error importing ${df.file.name}:`, err);
        }
      }

      if (sequenceFiles.length > 0) {
        try {
          await createAsuContentJob(project.id, sequenceFiles);
        } catch (err) {
          console.error("Error building the AU-contents job:", err);
          alert(
            "Could not build a 'Define AU contents' job from the dropped " +
              "sequence file(s): " +
              err
          );
        }
      }

      router.push(`/ccp4i2/project/${project.id}`);
    } catch (err) {
      console.error("Error creating project:", err);
      alert("Error creating project: " + err);
      setIsCreating(false);
    }
  }

  async function createImportJob(projectId: number, df: DroppedFile) {
    const taskName = TASK_FOR_TYPE[df.detectedType];
    const paramPath = PARAM_FOR_TYPE[df.detectedType];
    if (!taskName || !paramPath) return;

    // 1. Create the job
    const result = await apiPost<any>(`projects/${projectId}/create_task/`, {
      task_name: taskName,
      title: df.detectedType === "unmerged"
        ? `Merge ${df.file.name}`
        : `Import ${df.file.name}`,
    });
    if (!result?.success || !result.data?.new_job) {
      throw new Error(`Failed to create ${taskName} job`);
    }
    const jobId = result.data.new_job.id;

    // 2. Upload the file
    if (df.detectedType === "sequence") {
      // Sequences: read text content and set as parameter
      const text = await df.file.text();
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.${paramPath}`,
        value: text,
      });
    } else if (df.detectedType === "unmerged") {
      // Unmerged data: add list item, then upload to the list slot
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.${paramPath}`,
        value: [{}],
      });
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append(
        "objectPath",
        `${taskName}.container.${paramPath}[0].file`
      );
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    } else if (df.detectedType === "ligand") {
      // Ligands: the input parameter and the mode selector both depend on
      // the flavour of file. MOL mode covers .mol and .sdf alike (the task's
      // own menu calls it "a MOL or SDF file"); .smi is a SMILES *file*,
      // which is a different parameter from a typed-in SMILES string.
      const ext = df.file.name.toLowerCase().split(".").pop();
      const mode =
        ext === "mol2" ? "MOL2" : ext === "smi" ? "SMILESFILE" : "MOL";
      const actualParam =
        ext === "mol2"
          ? "inputData.MOL2IN"
          : ext === "smi"
            ? "inputData.SMILESFILEIN"
            : "inputData.MOLIN";
      await apiPost(`jobs/${jobId}/set_parameter/`, {
        object_path: `${taskName}.container.inputData.MOLSMILESORSKETCH`,
        value: mode,
      });
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append("object_path", `${taskName}.container.${actualParam}`);
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    } else {
      // Standard file upload
      const uploadForm = new FormData();
      uploadForm.append("file", df.file, df.file.name);
      uploadForm.append("object_path", `${taskName}.container.${paramPath}`);
      await apiPost(`jobs/${jobId}/upload_file_param/`, uploadForm);
    }

    // 3. Only auto-run simple import jobs; leave data reduction jobs
    //    for the user to review parameters before running
    if (AUTO_RUN_FOR_TYPE[df.detectedType]) {
      await apiPost(`jobs/${jobId}/run/`, {});
    }
  }

  /**
   * One ProvideAsuContents job accumulating every dropped sequence file.
   *
   * Each file becomes a CAsuContentSeq entry: the file itself lands on the
   * entry's `source`, and the server-side digest of that upload supplies the
   * sequence, name and polymer type. Stoichiometry is nobody's decision yet,
   * so every entry gets one copy - the job is a starting point for editing,
   * not a claim about the crystal.
   *
   * Robustness: a file the digest cannot read is reported by name and the
   * job is left PENDING (not run) so the user can fix or remove the entry;
   * files that did parse keep their entries either way.
   */
  async function createAsuContentJob(
    projectId: number,
    seqFiles: DroppedFile[]
  ) {
    const taskName = "ProvideAsuContents";
    const result = await apiPost<any>(`projects/${projectId}/create_task/`, {
      task_name: taskName,
      title:
        seqFiles.length === 1
          ? `Define AU contents (${seqFiles[0].file.name})`
          : `Define AU contents (${seqFiles.length} sequences)`,
    });
    if (!result?.success || !result.data?.new_job) {
      throw new Error(`Failed to create ${taskName} job`);
    }
    const jobId = result.data.new_job.id;
    const listPath = `${taskName}.container.inputData.ASU_CONTENT`;

    // One empty entry per dropped file
    await apiPost(`jobs/${jobId}/set_parameter/`, {
      object_path: listPath,
      value: seqFiles.map(() => ({})),
    });

    const problems: string[] = [];
    for (let i = 0; i < seqFiles.length; i++) {
      const df = seqFiles[i];
      try {
        // Upload the file onto this entry's `source` slot...
        const uploadForm = new FormData();
        uploadForm.append("file", df.file, df.file.name);
        uploadForm.append("objectPath", `${listPath}[${i}].source`);
        const up = await apiPost<any>(
          `jobs/${jobId}/upload_file_param/`,
          uploadForm
        );
        if (!up?.success) {
          throw new Error(up?.error || "upload failed");
        }

        // ...then let the server parse it
        const digestResp = await apiGet<any>(
          `jobs/${jobId}/digest?object_path=${listPath}[${i}].source`
        );
        const digest = digestResp?.data ?? digestResp;
        const sequence =
          typeof digest?.sequence === "string" ? digest.sequence.trim() : "";
        if (!sequence) {
          throw new Error(
            digest?.reason || "no sequence could be read from the file"
          );
        }

        const stem = df.file.name.replace(/\.[^.]+$/, "");
        const name =
          String(digest.name || digest.identifier || stem)
            .replace(/[^A-Za-z0-9_-]+/g, "_")
            .replace(/^_+|_+$/g, "") || `sequence_${i + 1}`;
        const polymerType = ["PROTEIN", "RNA", "DNA"].includes(
          digest.moleculeType
        )
          ? digest.moleculeType
          : "PROTEIN";
        const fields: Record<string, any> = {
          sequence,
          name,
          polymerType,
          nCopies: 1,
        };
        if (digest.description) fields.description = digest.description;
        for (const [key, value] of Object.entries(fields)) {
          await apiPost(`jobs/${jobId}/set_parameter/`, {
            object_path: `${listPath}[${i}].${key}`,
            value,
          });
        }
      } catch (err: any) {
        console.error(`Could not read a sequence from ${df.file.name}:`, err);
        problems.push(`${df.file.name}: ${err?.message ?? err}`);
      }
    }

    if (problems.length === 0) {
      await apiPost(`jobs/${jobId}/run/`, {});
    } else {
      alert(
        `A 'Define AU contents' job was created but NOT run - ` +
          `${problems.length} of ${seqFiles.length} sequence file(s) could not be read:\n\n` +
          problems.join("\n") +
          `\n\nOpen the job to fix or remove the entries, then run it. ` +
          `Copy numbers default to 1 - review the stoichiometry too.`
      );
    }
  }

  function handleNameChange(event: ChangeEvent<HTMLInputElement>) {
    setName(event.target.value);
  }

  function handleNameKeyDown(event: React.KeyboardEvent<HTMLInputElement>) {
    if (event.key === "Enter" && nameError.length === 0 && directoryError.length === 0) {
      createProject();
    }
  }

  function handleCustomDirectoryChange(
    event: React.MouseEvent<HTMLElement>,
    value: any
  ) {
    if (value !== null) {
      setCustomDirectory(value);
    }
  }

  function handleUseConfiguredRoot() {
    // Overrides the proposal for this project only. Nothing is persisted: the
    // configured root is already the preference, and a project created here
    // should not rewrite anyone's settings.
    setProposedOverride(configuredParent);
  }

  function handleDirectoryChange() {
    if (typeof window !== "undefined") {
      if (window.electronAPI) {
        // Deliberately NOT locate-ccp4i2-project-directory: that changes the
        // default projects directory for every future project. Choosing where
        // to put this one project must not move everyone else's default.
        window.electronAPI.sendMessage("choose-project-parent-directory");
      } else {
        console.error("Electron API is not available");
      }
    }
  }

  let nameError = "";
  if (name.length === 0) nameError = "Name is required";
  else if (!name.match("^[A-z0-9_-]+$"))
    nameError =
      "Name can only contain letters, numbers, underscores, and hyphens";
  else if (projects?.find((p) => p.name === name))
    nameError = "Name is already taken";

  const directoryError = useMemo(() => {
    if (customDirectory && directory.length === 0)
      return "Directory is required";
    else if (directory.length > 0 && directoryExists)
      return "Directory already exists";
    return "";
  }, [directoryExists, customDirectory, directory]);

  return (
    <Container maxWidth="sm" sx={{ my: 3 }}>
      <Stack spacing={2}>
        <Typography variant="h4">Create Project</Typography>
        <TextField
          label="Name"
          value={name}
          onChange={handleNameChange}
          onKeyDown={handleNameKeyDown}
          required
          error={nameError.length > 0}
          helperText={nameError}
          autoFocus
        />
        <TextField
          label="Description"
          value={description}
          onChange={(event) => setDescription(event.target.value)}
          multiline
          minRows={2}
          helperText="Optional. What the project is for; shown wherever it is listed."
        />
        <ToggleButtonGroup
          exclusive
          value={customDirectory}
          onChange={handleCustomDirectoryChange}
          fullWidth
        >
          {/* "Proposed", not "Default": this is where your last project went,
              which is not necessarily the configured projects directory. */}
          <ToggleButton value={false}>Proposed Directory</ToggleButton>
          <ToggleButton value={true}>Custom Directory</ToggleButton>
        </ToggleButtonGroup>
        {!customDirectory && effectiveProposal.length > 0 && (
          <Stack direction="row" spacing={2} sx={{ alignItems: "center" }}>
            <TextField
              label="Parent directory where the project directory will be created"
              value={effectiveProposal}
              disabled
              sx={{ flexGrow: 1 }}
              helperText={
                canUseConfiguredRoot
                  ? `Where your most recent project was created. Configured directory: ${configuredParent}`
                  : "Your configured projects directory."
              }
            />
            {canUseConfiguredRoot && (
              <Tooltip title={`Use ${configuredParent} instead`}>
                <Button
                  variant="outlined"
                  onClick={handleUseConfiguredRoot}
                  sx={{ whiteSpace: "nowrap" }}
                >
                  Use configured
                </Button>
              </Tooltip>
            )}
          </Stack>
        )}
        {customDirectory && electronAPIAvailable && (
          <Stack direction="row" spacing={2} sx={{ alignItems: "center" }}>
            <TextField
              label="Parent directory where the project directory will be created"
              value={customParent}
              disabled={true}
              sx={{ flexGrow: 1 }}
              required
            />
            <Button
              variant="outlined"
              startIcon={<Folder />}
              onClick={handleDirectoryChange}
            >
              Select
            </Button>
          </Stack>
        )}
        <Stack direction="row">
          <TextField
            label="Resulting name for project directory"
            value={directory}
            disabled={true}
            error={directoryError.length > 0}
            helperText={directoryError}
            sx={{ flexGrow: 1 }}
          />
        </Stack>

        <EditTags tags={tags} onChange={setTags} />

        <FileDropZone files={droppedFiles} onChange={setDroppedFiles} />

        <Stack direction="row" spacing={2} justifyContent="flex-end">
          <Button
            variant="outlined"
            onClick={() => router.push("/ccp4i2")}
            disabled={isCreating}
          >
            Cancel
          </Button>
          <Button
            variant="contained"
            onClick={createProject}
            disabled={
              nameError.length > 0 || directoryError.length > 0 || isCreating
            }
            startIcon={isCreating ? <CircularProgress size={16} /> : undefined}
          >
            {isCreating
              ? droppedFiles.length > 0
                ? "Creating & importing..."
                : "Creating..."
              : droppedFiles.length > 0
                ? `Create & import ${droppedFiles.length} file${droppedFiles.length > 1 ? "s" : ""}`
                : "Create"}
          </Button>
        </Stack>
      </Stack>
    </Container>
  );
};
