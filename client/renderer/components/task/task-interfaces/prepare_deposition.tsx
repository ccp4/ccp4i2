import { useCallback, useMemo, useRef, useState } from "react";
import { Alert, Button, LinearProgress, List, ListItemButton, ListItemText, Paper, TextField, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { useApi } from "../../../api";
import { Project } from "../../../types/models";

interface LineageResult {
  success: boolean;
  refinement_job_id?: number;
  refinement_task_name?: string;
  files: Record<string, string | string[]>;
  paths: Record<string, string>;
  params: Record<string, boolean | string>;
  warnings: string[];
  software?: { name: string; version?: string; classification?: string }[];
}

interface ReferenceHit {
  pdb_entity: string;
  uniprot: string;
  identity: number | null;
  coverage: number | null;
  description: string | null;
}

interface ReferenceResult {
  ok: boolean;
  error?: string;
  hits: ReferenceHit[];
}

interface FillResult {
  ok: boolean;
  error?: string;
  entities: string[];
  software: string[];
  wavelength: number | null;
  notes: string[];
}

type Status = {
  type: "success" | "info" | "warning" | "error";
  message: string;
} | null;

/**
 * Prepare files for deposition.
 *
 * Files: the refined model drives everything -- choosing it traces the
 * refinement and scaling jobs and fills the reflections, free set, TLS,
 * dictionaries and statistics. Entry, Molecules, Experiment and Assembly are
 * the metadata OneDep reads from the model file; the Molecules tab starts
 * from one record per AU-contents sequence and each record describes the
 * segments its construct is made of (natural, engineered, designed, a known
 * scaffold, a synthesised peptide, a tag, a linker).
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const {
    container,
    callPluginMethod,
    setParameter,
    mutateContainer,
    fileItemToParameterArg,
  } = useJob(props.job.id);
  const { jobs: projectJobs } = useProject(props.job.project);
  const { files: projectFiles } = useProjectFiles(props.job.project);
  const api = useApi();
  const { data: projects } = api.get<Project[]>("projects");
  const [status, setStatus] = useState<Status>(null);
  const [busy, setBusy] = useState(false);
  // A ref as well as state: two clicks inside one render both see busy=false.
  const busyRef = useRef(false);
  const [refQuery, setRefQuery] = useState("");
  const [refHits, setRefHits] = useState<ReferenceHit[] | null>(null);
  const taskName = props.job.task_name;

  /** Runs one action at a time; a second click while it runs is ignored. */
  const exclusive = useCallback(async (action: () => Promise<void>) => {
    if (busyRef.current) return;
    busyRef.current = true;
    setBusy(true);
    try {
      await action();
    } finally {
      busyRef.current = false;
      setBusy(false);
    }
  }, []);

  const filesByUuid = useMemo(() => {
    const map = new Map<string, any>();
    for (const f of projectFiles ?? []) map.set(f.uuid.replace(/-/g, ""), f);
    return map;
  }, [projectFiles]);

  const handleXYZINChange = useCallback(
    async (updatedItem: any) => {
      const dbFileId =
        updatedItem?._value?.dbFileId?._value?.trim() || updatedItem?.dbFileId;
      if (!dbFileId) return;
      const fileUuid = [
        dbFileId.slice(0, 8),
        dbFileId.slice(8, 12),
        dbFileId.slice(12, 16),
        dbFileId.slice(16, 20),
        dbFileId.slice(20),
      ].join("-");
      setBusy(true);
      setStatus(null);
      const result: LineageResult | null = await callPluginMethod(
        "trace_xyzin_lineage",
        { file_uuid: fileUuid }
      );
      if (!result || !result.success) {
        setBusy(false);
        setStatus({
          type: "warning",
          message: (result?.warnings ?? ["Could not trace the model's lineage"]).join(". "),
        });
        return;
      }
      for (const [param, uuid] of Object.entries(result.files)) {
        if (param === "DICT_LIST" || typeof uuid !== "string") continue;
        const file = filesByUuid.get(uuid.replace(/-/g, ""));
        const objectPath = `${taskName}.inputData.${param}`;
        if (file && projectJobs && projects) {
          await setParameter(fileItemToParameterArg(file, objectPath, projectJobs, projects));
        } else {
          await setParameter({ object_path: objectPath, value: { dbFileId: uuid.replace(/-/g, "") } });
        }
      }
      for (const [param, path] of Object.entries(result.paths)) {
        await setParameter({ object_path: `${taskName}.inputData.${param}`, value: path });
      }
      for (const [param, val] of Object.entries(result.params)) {
        await setParameter({ object_path: `${taskName}.controlParameters.${param}`, value: val });
      }
      // Program versions the scaling job recorded become software rows. The
      // server decides what is new against the saved parameters (in OneDep's
      // spelling), so tracing the same model twice adds nothing.
      if ((result.software ?? []).length) {
        await callPluginMethod("fill_form", { software: result.software });
      }
      await mutateContainer();
      setBusy(false);
      const label = result.refinement_task_name ?? "refinement";
      setStatus({
        type: result.warnings.length ? "info" : "success",
        message:
          `Filled from ${label} job #${result.refinement_job_id}` +
          (result.warnings.length ? `. Note: ${result.warnings.join(". ")}` : ""),
      });
    },
    [callPluginMethod, setParameter, mutateContainer, fileItemToParameterArg, filesByUuid, projectJobs, projects, taskName]
  );

  // One Molecule record per AU-contents sequence, with proposed segments
  // (tag, protein), and the wavelength. Decided on the server against the
  // saved parameters: a sequence that already has a record is skipped, so
  // filling twice changes nothing.
  const fillForm = useCallback(
    () =>
      exclusive(async () => {
        setStatus(null);
        await mutateContainer();
        const added: FillResult | null = await callPluginMethod("fill_form", {});
        await mutateContainer();
        if (!added || !added.ok) {
          setStatus({ type: "warning", message: added?.error ?? "Nothing could be filled in" });
          return;
        }
        const parts: string[] = [];
        if (added.entities.length) parts.push(`molecule record(s) for ${added.entities.join(", ")}`);
        if (added.wavelength) parts.push(`wavelength ${added.wavelength} A`);
        setStatus({
          type: parts.length ? "success" : "info",
          message: parts.length
            ? `Added ${parts.join(" and ")}. ${added.notes.join(". ")}`
            : "Every AU-contents sequence already has a record; nothing added.",
        });
      }),
    [callPluginMethod, exclusive, mutateContainer]
  );

  // Find UniProt accessions for an AU-contents sequence: offered, never filled in.
  const findReferences = useCallback(
    () =>
      exclusive(async () => {
        setRefHits(null);
        const result: ReferenceResult | null = await callPluginMethod("find_references", {
          sequence_name: refQuery,
        });
        if (!result || !result.ok) {
          setStatus({ type: "warning", message: result?.error ?? "The search could not be run" });
          return;
        }
        setRefHits(result.hits);
        if (!result.hits.length) setStatus({ type: "info", message: "No similar entries found." });
      }),
    [callPluginMethod, exclusive, refQuery]
  );

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {busy && <LinearProgress />}
      {status && (
        <Alert severity={status.type} onClose={() => setStatus(null)}>
          {status.message}
        </Alert>
      )}
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Files">
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Model and contents" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="XYZIN" {...props} onChange={handleXYZINChange} />
            <CCP4i2TaskElement itemName="ASUCONTENT" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Reflections and related files (filled from the model's lineage)" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="F_SIGF" {...props} qualifiers={{ toolTip: "The reflections the model was refined against" }} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />
            <CCP4i2TaskElement itemName="SCALEDUNMERGED" {...props} />
            <CCP4i2TaskElement itemName="CIFSTATSOUT" {...props} />
            <CCP4i2TaskElement itemName="AIMLESSXML" {...props} />
            <CCP4i2TaskElement itemName="REFMACINPUTPARAMSXML" {...props} />
            <CCP4i2TaskElement itemName="TLSIN" {...props} />
            <CCP4i2TaskElement itemName="DICT_LIST" {...props} qualifiers={{ toolTip: "Restraint dictionaries used in refinement" }} />
            <CCP4i2TaskElement itemName="inputData.FPHIOUT" {...props} />
            <CCP4i2TaskElement itemName="inputData.DIFFPHIOUT" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Options" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="USEANOMALOUS" {...props} />
            <CCP4i2TaskElement itemName="USE_TWIN" {...props} />
            <CCP4i2TaskElement itemName="USEAIMLESSXML" {...props} />
            <CCP4i2TaskElement itemName="INCLUDEUNMERGED" {...props} />
            <CCP4i2TaskElement itemName="KEEP_HYDROGENS" {...props} />
            <CCP4i2TaskElement itemName="KEEP_ALL_SHELLS" {...props} />
            <CCP4i2TaskElement itemName="DICTIONARY_CHECK" {...props} />
            <CCP4i2TaskElement itemName="SENDTOVALIDATIONSERVER" {...props} qualifiers={{ toolTip: "Needs internet access; takes some minutes" }} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Entry">
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Title and classification" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="TITLE" {...props} />
            <CCP4i2TaskElement itemName="KEYWORDS" {...props} />
            <CCP4i2TaskElement itemName="KEYWORD_TEXT" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Deposition authors" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="AUTHORS" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Primary citation" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="CITATION_PUBLISHED" {...props} />
            <CCP4i2TaskElement itemName="CITATION_TITLE" {...props} />
            <CCP4i2TaskElement itemName="JOURNAL" {...props} />
            <CCP4i2TaskElement itemName="VOLUME" {...props} />
            <CCP4i2TaskElement itemName="FIRST_PAGE" {...props} />
            <CCP4i2TaskElement itemName="LAST_PAGE" {...props} />
            <CCP4i2TaskElement itemName="YEAR" {...props} />
            <CCP4i2TaskElement itemName="DOI" {...props} />
            <CCP4i2TaskElement itemName="PUBMED" {...props} />
            <CCP4i2TaskElement itemName="CITATION_AUTHORS" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Contact author" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="CONTACT_SALUTATION" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_FIRST_NAME" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_LAST_NAME" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_ROLE" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_ORG_TYPE" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_ADDRESS_1" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_ADDRESS_2" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_CITY" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_STATE" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_POSTCODE" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_COUNTRY" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_EMAIL" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_PHONE" {...props} />
            <CCP4i2TaskElement itemName="CONTACT_ORCID" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Related entries" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="RELATED" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Molecules">
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Molecules in the crystal" }} containerHint="FolderLevel">
            <Typography variant="body2" sx={{ mb: 1 }}>
              One record per sequence of the AU contents, in the same order. Describe each
              construct as segments: a natural protein or an engineered variant (UniProt), a
              de novo design (no reference needed), a scaffold the PDB already holds, a
              synthesised peptide, a tag or a linker. A designed protein needs only its name
              and one designed segment.
            </Typography>
            <Button variant="outlined" size="small" onClick={fillForm} disabled={busy} sx={{ mb: 1 }}>
              Fill from the AU contents
            </Button>
            <Typography variant="body2" sx={{ mt: 1 }}>
              Find UniProt references for a sequence (RCSB search; you choose, nothing is filled in):
            </Typography>
            <TextField
              size="small"
              label="AU-contents sequence name"
              value={refQuery}
              onChange={(e) => setRefQuery(e.target.value)}
              sx={{ mr: 1, mt: 0.5 }}
            />
            <Button variant="outlined" size="small" onClick={findReferences} disabled={busy || !refQuery} sx={{ mt: 0.5 }}>
              Find UniProt reference
            </Button>
            {refHits && refHits.length > 0 && (
              <List dense>
                {refHits.map((hit) => (
                  <ListItemButton
                    key={`${hit.pdb_entity}-${hit.uniprot}`}
                    onClick={() => {
                      navigator.clipboard?.writeText(hit.uniprot);
                      setStatus({
                        type: "info",
                        message: `${hit.uniprot} copied: paste it into the UniProt field of the segment it belongs to.`,
                      });
                    }}
                  >
                    <ListItemText
                      primary={`${hit.uniprot} (via PDB ${hit.pdb_entity})`}
                      secondary={`${hit.description ?? ""} identity ${hit.identity ?? "?"} coverage ${hit.coverage ?? "?"}`}
                    />
                  </ListItemButton>
                ))}
              </List>
            )}
            <CCP4i2TaskElement itemName="ENTITIES" {...props} />
            <CCP4i2TaskElement itemName="DESIGN_DETAILS" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Experiment">
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Data collection" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="SOURCE" {...props} />
            <CCP4i2TaskElement itemName="SYNCHROTRON_SITE" {...props} />
            <CCP4i2TaskElement itemName="BEAMLINE" {...props} />
            <CCP4i2TaskElement itemName="DETECTOR" {...props} />
            <CCP4i2TaskElement itemName="DETECTOR_TYPE" {...props} />
            <CCP4i2TaskElement itemName="COLLECTION_DATE" {...props} />
            <CCP4i2TaskElement itemName="TEMPERATURE" {...props} />
            <CCP4i2TaskElement itemName="WAVELENGTH" {...props} />
            <CCP4i2TaskElement itemName="MONOCHROMATOR" {...props} />
            <CCP4i2TaskElement itemName="PROTOCOL" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Crystallisation" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="CRYSTAL_METHOD" {...props} />
            <CCP4i2TaskElement itemName="CRYSTAL_TEMP" {...props} />
            <CCP4i2TaskElement itemName="CRYSTAL_PH" {...props} />
            <CCP4i2TaskElement itemName="CRYSTAL_DETAILS" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Phasing and software" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="PHASING_METHOD" {...props} />
            <CCP4i2TaskElement itemName="STARTING_MODEL" {...props} />
            <CCP4i2TaskElement itemName="PHASING_DETAILS" {...props} />
            <CCP4i2TaskElement itemName="STARTING_MODELS" {...props} />
            <CCP4i2TaskElement itemName="SOFTWARE" {...props} />
          </CCP4i2ContainerElement>
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Ligands" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="LIGAND_OF_INTEREST" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        <CCP4i2Tab label="Assembly">
          <CCP4i2ContainerElement {...props} itemName="" qualifiers={{ guiLabel: "Biological assembly" }} containerHint="FolderLevel">
            <CCP4i2TaskElement itemName="ASSEMBLY_MODE" {...props} />
            <CCP4i2TaskElement itemName="OLIGOMERIC_COUNT" {...props} />
            <CCP4i2TaskElement itemName="ASSEMBLY_OPERATORS" {...props} />
            <CCP4i2TaskElement itemName="ASSEMBLY_EVIDENCE" {...props} />
            <CCP4i2TaskElement itemName="ASSEMBLY_DETAILS" {...props} />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
