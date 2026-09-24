import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import {
  Alert,
  Box,
  Button,
  IconButton,
  LinearProgress,
  Paper,
  Popover,
  Typography,
} from "@mui/material";
import { HelpOutline } from "@mui/icons-material";

import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../../task-elements/tabs";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { InlineField } from "../../task-elements/inline-field";
import { useContainerList } from "../../task-elements/hooks/useContainerList";
import { useJob } from "../../../../utils";

import {
  Assembly,
  EMPTY_ASSEMBLY,
  Model,
  ParsedBody,
  RawBody,
  formatAssemblyRows,
  formatBodies,
  isTerse,
  parseAssemblyRows,
  parseBodies,
  parseSegments,
} from "./dm-spec";
import { AssemblyGrid } from "./assembly-grid";
import { BodiesEditor } from "./bodies-editor";
import { CoverageStrip } from "./coverage-strip";
import { useNcsPreview } from "./use-ncs-preview";

/**
 * dm_multidomain — multi-domain NCS averaging via dm.
 *
 * The task asks two questions: which chains are copies of each other, and
 * which parts of a copy move as one unit. Both used to be asked as free text
 * in two ad-hoc grammars, cross-referenced by a role name the user invented,
 * with three load-bearing conventions (row 0 is the reference, an omitted role
 * is a partial copy, an empty list means auto-detect) that could only be
 * conveyed in prose. All of it is derivable from the model the user has
 * already chosen, so the interface reads the model instead of asking:
 *
 *   - it opens on the assembly and body the model implies, rather than blank;
 *   - the assembly is a grid of chain pickers, so roles are columns and a
 *     partial copy is a visible hole;
 *   - a body's ranges are bounded by the chain's own numbering;
 *   - the coverage strip draws the bodies on the residues, where a gap or an
 *     overlap can be seen rather than deduced;
 *   - and each body reports, per copy, its matched CA count and superposition
 *     RMSD, which is the answer to "do these residues really move as one unit"
 *     arriving before the job is run.
 *
 * The stored format is unchanged: dm-spec.ts serialises back to the same
 * "CDK=A cyclin=B" rows and "cyclin:10-95,CDK:45-60" specs that i2run writes.
 */

const HelpButton: React.FC<{ children: React.ReactNode; label: string }> = ({
  children,
  label,
}) => {
  const [anchor, setAnchor] = useState<HTMLElement | null>(null);
  return (
    <>
      <IconButton size="small" aria-label={label} onClick={(e) => setAnchor(e.currentTarget)}>
        <HelpOutline fontSize="small" />
      </IconButton>
      <Popover
        open={Boolean(anchor)}
        anchorEl={anchor}
        onClose={() => setAnchor(null)}
        anchorOrigin={{ vertical: "bottom", horizontal: "left" }}
      >
        <Box sx={{ p: 2, maxWidth: "34rem" }}>{children}</Box>
      </Popover>
    </>
  );
};

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);
  const { value: PHASE_SOURCE } = useTaskItem("PHASE_SOURCE");
  const { value: XYZIN } = useTaskItem("XYZIN");
  const fromModel = PHASE_SOURCE === "model";
  const editable = job.status === 1;

  const assemblyList = useContainerList({ job, itemName: "ASSEMBLY" });
  const domainsList = useContainerList({ job, itemName: "DOMAINS" });

  // ---- the two parameters, read as one model ------------------------------
  const assemblyRows: string[] = useMemo(
    () => (assemblyList.items || []).map((row: any) => String(row?._value ?? row ?? "")),
    [assemblyList.items]
  );
  const domainRows: RawBody[] = useMemo(
    () =>
      (domainsList.items || []).map((row: any) => ({
        segments: String(row?._value?.segments?._value ?? ""),
        mode: String(row?._value?.mode?._value ?? "average"),
      })),
    [domainsList.items]
  );

  const model: Model = useMemo(
    () => ({
      assembly: assemblyRows.length ? parseAssemblyRows(assemblyRows) : EMPTY_ASSEMBLY,
      bodies: parseBodies(domainRows),
    }),
    [assemblyRows, domainRows]
  );

  // ---- what the server can tell us about the model ------------------------
  const hasModel = Boolean(XYZIN?.dbFileId);
  const previewKey = useMemo(
    () => JSON.stringify([XYZIN?.dbFileId ?? null, assemblyRows, domainRows]),
    [XYZIN?.dbFileId, assemblyRows, domainRows]
  );
  const { preview, isLoading } = useNcsPreview(job, previewKey, hasModel);

  const chains = preview?.model?.chains ?? [];
  const entities = preview?.model?.entities ?? [];

  // ---- writing back -------------------------------------------------------
  // Roles are the one thing the two parameters share, so a change to them has
  // to reach both. Each list is written only when its own serialisation
  // actually differs, so an edit to one does not churn the other.
  const commit = useCallback(
    async (next: Model) => {
      const rows = formatAssemblyRows(next.assembly);
      const terse = isTerse(next.assembly.roles);
      const bodies = formatBodies(next.bodies, terse);
      if (JSON.stringify(rows) !== JSON.stringify(assemblyRows)) {
        await assemblyList.replaceArray(rows);
      }
      if (JSON.stringify(bodies) !== JSON.stringify(domainRows)) {
        await domainsList.replaceArray(bodies);
      }
    },
    [assemblyRows, domainRows, assemblyList, domainsList]
  );

  const commitBodies = useCallback(
    (bodies: ParsedBody[]) => commit({ ...model, bodies }),
    [commit, model]
  );

  // ---- open on something that works --------------------------------------
  // "Leave the list empty to auto-detect" is magic by absence: the user is
  // shown a blank box and told in prose what the blank means. Instead the
  // detected assembly is written in, where it can be read and corrected. The
  // guard is per model file, so clearing a row afterwards stays cleared.
  const prefilledFor = useRef<string | null>(null);
  useEffect(() => {
    if (!editable || !preview?.ok || !preview.suggestion) return;
    const fileId = String(XYZIN?.dbFileId ?? "");
    if (!fileId || prefilledFor.current === fileId) return;
    prefilledFor.current = fileId;

    const assemblyEmpty = assemblyRows.length === 0;
    const bodiesBlank =
      domainRows.length === 0 ||
      (domainRows.length === 1 && !String(domainRows[0].segments ?? "").trim());
    if (!assemblyEmpty && !bodiesBlank) return;

    const suggestedAssembly = assemblyEmpty
      ? parseAssemblyRows(preview.suggestion.assembly)
      : model.assembly;
    const suggestedBodies: ParsedBody[] = bodiesBlank
      ? [{ segments: parseSegments(preview.suggestion.segments), mode: "average" }]
      : model.bodies;
    if (assemblyEmpty && suggestedAssembly.instances.length < 2 && bodiesBlank === false) {
      return;
    }
    commit({ assembly: suggestedAssembly, bodies: suggestedBodies });
  }, [editable, preview, XYZIN?.dbFileId, assemblyRows, domainRows, model, commit]);

  const detect = useCallback(() => {
    if (!preview?.suggestion) return;
    commit({
      assembly: parseAssemblyRows(preview.suggestion.assembly),
      bodies: model.bodies,
    });
  }, [preview, model.bodies, commit]);

  const [hoveredBody, setHoveredBody] = useState<number | null>(null);

  const assemblyForDisplay: Assembly = model.assembly.instances.length
    ? model.assembly
    : parseAssemblyRows(preview?.suggestion?.assembly ?? []);

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Reflection Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections (F/SIGF or intensities)" }}
            />
            <InlineField label="Starting phases from">
              <CCP4i2TaskElement itemName="PHASE_SOURCE" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <CCP4i2TaskElement
              itemName="ABCD"
              {...props}
              qualifiers={{ guiLabel: "Starting phases" }}
              visibility={() => !fromModel}
            />
            {fromModel && (
              <Typography variant="body2" sx={{ pl: 1, color: "text.secondary" }}>
                Phases will be calculated from the model with servalcat sigmaa
                (bulk solvent + sigmaA weighting).
              </Typography>
            )}
            <CCP4i2TaskElement
              itemName="FREERFLAG"
              {...props}
              qualifiers={{ guiLabel: "Free R set (optional)" }}
            />
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "NCS model" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="XYZIN"
              {...props}
              qualifiers={{ guiLabel: "Model (the copies are read from this)" }}
            />
            {preview?.model && (
              <Typography variant="body2" color="text.secondary" sx={{ pl: 1 }}>
                {preview.model.chains.length} protein chain
                {preview.model.chains.length === 1 ? "" : "s"} in{" "}
                {preview.model.entities.length} entit
                {preview.model.entities.length === 1 ? "y" : "ies"};{" "}
                {preview.model.nCopiesDetected} cop
                {preview.model.nCopiesDetected === 1 ? "y" : "ies"} of the
                assembly detected.
              </Typography>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Domains ===== */}
        <CCP4i2Tab label="Domains">
          {!hasModel ? (
            <Alert severity="info">
              Choose a model on the Input Data tab. The copies and their residue
              ranges are read from it.
            </Alert>
          ) : (
            <>
              {isLoading && !preview && <LinearProgress />}
              {preview && !preview.ok && (
                <Alert severity="warning" sx={{ mb: 1 }}>
                  Could not read the model: {preview.error}
                </Alert>
              )}

              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: "Which chains are copies of each other?" }}
                containerHint="BlockLevel"
              >
                <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
                  <HelpButton label="About NCS copies and entities">
                    <Typography variant="subtitle2" gutterBottom>
                      Copies and entities
                    </Typography>
                    <Typography variant="body2" paragraph>
                      A <b>copy</b> is one instance of the thing that repeats in
                      the asymmetric unit; an <b>entity</b> is one of the
                      proteins it is made of. A homomer has a single entity, so
                      each copy is just a chain. A CDK/cyclin complex has two,
                      so each copy names a chain for each: CDK=A cyclin=B, then
                      CDK=C cyclin=D.
                    </Typography>
                    <Typography variant="body2" paragraph>
                      The first row is the <b>reference</b>: the averaging masks
                      are cut from it and every other copy is superposed onto
                      it.
                    </Typography>
                    <Typography variant="body2">
                      Leave a cell empty where a copy does not have that entity
                      — an A<sub>4</sub>B<sub>3</sub> complex has one copy with
                      no B, and bodies that need B simply skip it.
                    </Typography>
                  </HelpButton>
                </Box>
                <AssemblyGrid
                  model={model}
                  chains={chains}
                  entities={entities}
                  disabled={!editable}
                  onChange={commit}
                  onDetect={preview?.suggestion ? detect : undefined}
                  detectLabel="Re-detect from model"
                />
              </CCP4i2ContainerElement>

              <CCP4i2ContainerElement
                {...props}
                itemName=""
                qualifiers={{ guiLabel: "Which parts move as one unit?" }}
                containerHint="BlockLevel"
              >
                <Box sx={{ display: "flex", justifyContent: "flex-end" }}>
                  <HelpButton label="About rigid bodies">
                    <Typography variant="subtitle2" gutterBottom>
                      Rigid bodies
                    </Typography>
                    <Typography variant="body2" paragraph>
                      A <b>rigid body</b> is a set of residue ranges that move
                      together between copies. Each body is superposed on its
                      own, so different bodies can follow different NCS
                      operators — which is the whole point of this task, and
                      what plain NCS averaging cannot do.
                    </Typography>
                    <Typography variant="body2" paragraph>
                      A body can span entities: the CDK C-helix that travels
                      with the cyclin N-lobe is one body with two ranges, one on
                      each. Add a second range and pick the other entity for it.
                    </Typography>
                    <Typography variant="body2">
                      <b>average</b> holds the fitted operators fixed;{" "}
                      <b>refine</b> lets dm improve them as it goes;{" "}
                      <b>exclude</b> leaves the body out of averaging entirely.
                    </Typography>
                  </HelpButton>
                </Box>

                <CoverageStrip
                  assembly={assemblyForDisplay}
                  bodies={model.bodies}
                  chains={chains}
                  highlight={hoveredBody}
                />

                <BodiesEditor
                  bodies={model.bodies}
                  assembly={assemblyForDisplay}
                  chains={chains}
                  preview={preview?.bodies}
                  disabled={!editable}
                  onChange={commitBodies}
                  onHover={setHoveredBody}
                />

                {preview?.suggestion?.segments && (
                  <Button
                    size="small"
                    disabled={!editable}
                    sx={{ mt: 1 }}
                    onClick={() =>
                      commitBodies([
                        {
                          segments: parseSegments(preview.suggestion!.segments),
                          mode: "average",
                        },
                      ])
                    }
                  >
                    Reset to one body covering the whole copy
                  </Button>
                )}
              </CCP4i2ContainerElement>
            </>
          )}
        </CCP4i2Tab>

        {/* ===== Parameters ===== */}
        <CCP4i2Tab label="Parameters">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Density modification" }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="MODE_SOLVENT"
              {...props}
              qualifiers={{ guiLabel: "Solvent flattening" }}
            />
            <CCP4i2TaskElement
              itemName="MODE_HISTOGRAM"
              {...props}
              qualifiers={{ guiLabel: "Histogram matching" }}
            />
            <InlineField label="Number of cycles">
              <CCP4i2TaskElement itemName="NCYCLES" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
            <InlineField label="Solvent content" hint="blank = estimate from model and cell">
              <CCP4i2TaskElement
                itemName="SOLVENT_CONTENT"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Mask radius" hint="Angstrom">
              <CCP4i2TaskElement itemName="MASK_RADIUS" {...props} qualifiers={{ guiLabel: " " }} />
            </InlineField>
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
