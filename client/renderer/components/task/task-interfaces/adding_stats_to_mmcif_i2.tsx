import { useCallback, useMemo, useState } from "react";
import { Alert, LinearProgress, Paper } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob, useProject, useProjectFiles } from "../../../utils";
import { useApi } from "../../../api";
import { Project } from "../../../types/models";

interface LineageResult {
  success: boolean;
  refinement_job_id?: number;
  refinement_task_name?: string;
  files: Record<string, string | string[]>;
  paths: Record<string, string>;
  params: Record<string, boolean>;
  warnings: string[];
}

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, callPluginMethod, setParameter, mutateContainer, fileItemToParameterArg } =
    useJob(props.job.id);
  const { jobs: projectJobs } = useProject(props.job.project);
  const { files: projectFiles } = useProjectFiles(props.job.project);
  const { data: projects } = useApi().get<Project[]>("projects");
  const [status, setStatus] = useState<{
    type: "success" | "info" | "warning" | "error";
    message: string;
  } | null>(null);
  const [populating, setPopulating] = useState(false);

  // Build a UUID→File lookup from project files for fast resolution
  const filesByUuid = useMemo(() => {
    const map = new Map<string, any>();
    if (projectFiles) {
      for (const f of projectFiles) {
        map.set(f.uuid.replace(/-/g, ""), f);
      }
    }
    return map;
  }, [projectFiles]);

  const handleXYZINChange = useCallback(
    async (updatedItem: any) => {
      // Extract dbFileId from the serialized container item structure:
      // { _value: { dbFileId: { _value: "abc123..." } } }
      const dbFileId =
        updatedItem?._value?.dbFileId?._value?.trim() ||
        updatedItem?.dbFileId;
      if (!dbFileId) return;

      // Convert dashless UUID back to standard format for the server
      const fileUuid = [
        dbFileId.slice(0, 8),
        dbFileId.slice(8, 12),
        dbFileId.slice(12, 16),
        dbFileId.slice(16, 20),
        dbFileId.slice(20),
      ].join("-");

      setPopulating(true);
      setStatus(null);

      const result: LineageResult | null = await callPluginMethod(
        "trace_xyzin_lineage",
        { file_uuid: fileUuid }
      );

      if (!result || !result.success) {
        setPopulating(false);
        const warnings = result?.warnings ?? ["Unknown error"];
        setStatus({ type: "warning", message: warnings.join(". ") });
        return;
      }

      // The backend set_parameter uses skip_first=True, expecting the task name
      // as the first path element (e.g., "adding_stats_to_mmcif_i2.inputData.FREERFLAG").
      // Without the prefix, "inputData.FREERFLAG" gets "inputData" stripped → lookup fails.
      const taskName = props.job.task_name;

      // Auto-populate DB file fields using full file metadata (not just dbFileId).
      // fileItemToParameterArg includes baseName, relPath, annotation, etc. which
      // are required for the file to persist through saveDataToXml/reload cycles.
      for (const [param, uuid] of Object.entries(result.files)) {
        if (param === "DICT_LIST") continue; // TODO: handle list fields
        if (typeof uuid === "string") {
          const dashlessUuid = uuid.replace(/-/g, "");
          const file = filesByUuid.get(dashlessUuid);
          const objectPath = `${taskName}.inputData.${param}`;
          if (file && projectJobs && projects) {
            await setParameter(fileItemToParameterArg(file, objectPath, projectJobs, projects));
          } else {
            // Fallback: at minimum set dbFileId
            await setParameter({ object_path: objectPath, value: { dbFileId: dashlessUuid } });
          }
        }
      }

      // Set filesystem-path files (AIMLESSXML, REFMACINPUTPARAMSXML)
      for (const [param, path] of Object.entries(result.paths)) {
        await setParameter({
          object_path: `${taskName}.inputData.${param}`,
          value: path,
        });
      }

      // Set control parameters
      for (const [param, val] of Object.entries(result.params)) {
        await setParameter({
          object_path: `${taskName}.controlParameters.${param}`,
          value: val,
        });
      }

      // Force full container re-fetch so all CCP4i2TaskElement components re-render.
      await mutateContainer();

      setPopulating(false);

      const taskLabel = result.refinement_task_name ?? "refinement";
      setStatus({
        type: "success",
        message: `Auto-populated from ${taskLabel} job #${result.refinement_job_id}`,
      });

      // Show any non-fatal warnings
      if (result.warnings.length > 0) {
        setStatus({
          type: "info",
          message: `Populated from ${taskLabel} job #${result.refinement_job_id}. Note: ${result.warnings.join(". ")}`,
        });
      }
    },
    [callPluginMethod, setParameter, mutateContainer, fileItemToParameterArg, filesByUuid, projectJobs, projects, props.job.task_name]
  );

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* ASU content */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Asymmetric unit content (i.e. sequences)" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="ASUCONTENT"
          {...props}
          qualifiers={{ toolTip: "ASU Content" }}
        />
      </CCP4i2ContainerElement>

      {/* Coordinates */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Coordinates" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="XYZIN"
          {...props}
          qualifiers={{ toolTip: "Input model — select from a refinement job" }}
          onChange={handleXYZINChange}
        />
        {populating && <LinearProgress />}
        {status && (
          <Alert severity={status.type} sx={{ mt: 0.5 }} onClose={() => setStatus(null)}>
            {status.message}
          </Alert>
        )}
      </CCP4i2ContainerElement>

      {/* Task control */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Task control" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="SENDTOVALIDATIONSERVER"
          {...props}
          qualifiers={{
            guiLabel: "Use validation server",
            toolTip: "Send to validation server - requires internet access",
          }}
        />
        <CCP4i2TaskElement
          itemName="USEAIMLESSXML"
          {...props}
          qualifiers={{
            guiLabel: "Add statistics from Aimless",
            toolTip:
              "If you didn't run data reduction in CCP4i2 you won't have this",
          }}
        />
        <CCP4i2TaskElement
          itemName="INCLUDEUNMERGED"
          {...props}
          qualifiers={{
            guiLabel: "Include unmerged from scaling job",
            toolTip:
              "If you didn't run data reduction in CCP4i2 you won't have this",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Related files */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Related files" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="USEANOMALOUS"
          {...props}
          qualifiers={{ guiLabel: "Anomalous used in refinement job" }}
        />
        <CCP4i2TaskElement
          itemName="USE_TWIN"
          {...props}
          qualifiers={{
            guiLabel: "Twinning (i.e. intensities) used in refinement job",
          }}
        />
        <CCP4i2TaskElement
          itemName="F_SIGF"
          {...props}
          qualifiers={{ toolTip: "Input reflections" }}
        />
        <CCP4i2TaskElement itemName="SCALEDUNMERGED" {...props} />
        <CCP4i2TaskElement itemName="AIMLESSXML" {...props} />
        <CCP4i2TaskElement itemName="REFMACINPUTPARAMSXML" {...props} />
        <CCP4i2TaskElement
          itemName="FREERFLAG"
          {...props}
          qualifiers={{ toolTip: "FreeR flag" }}
        />
        <CCP4i2TaskElement
          itemName="TLSIN"
          {...props}
          qualifiers={{ toolTip: "Input TLS" }}
        />
        <CCP4i2TaskElement
          itemName="DICT_LIST"
          {...props}
          qualifiers={{ toolTip: "Dictionary list" }}
        />
        <CCP4i2TaskElement itemName="inputData.FPHIOUT" {...props} />
        <CCP4i2TaskElement itemName="inputData.DIFFPHIOUT" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
