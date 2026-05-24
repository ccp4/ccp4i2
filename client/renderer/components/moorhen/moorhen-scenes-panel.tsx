/**
 * Moorhen Scenes side-panel.
 *
 * One home for all yaml-scene operations:
 *   - Live YAML editor (Monaco) with syntax highlighting
 *   - Capture current view → editor (lifter)
 *   - Open a .scene.yaml from disk → editor
 *   - Save current editor contents → browser download
 *   - Apply now (parse + resolve + render)
 *   - Live apply (debounced, silent on parse errors)
 *   - Validation/log panel below the editor (parse errors, resolver log)
 *
 * The Apply / Download / Open actions previously lived inline next to
 * Copy View Link; they've been consolidated here so there's one place
 * for everything yaml-scene-related.
 */

import React, {
  useCallback,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import {
  Box,
  Button,
  Checkbox,
  FormControlLabel,
  IconButton,
  Stack,
  Tooltip,
  Typography,
} from "@mui/material";
import {
  Camera as CaptureIcon,
  FileOpen as OpenIcon,
  Save as SaveIcon,
  PlayArrow as ApplyIcon,
} from "@mui/icons-material";
import { Editor } from "@monaco-editor/react";
import JSZip from "jszip";

import { parseScene, serialiseSceneWithComments, SceneParseError } from "../../lib/moorhen-scene";
import type { MoorhenScene } from "../../types/moorhen-scene";
import type {
  SceneResolveResult,
  SceneResolveLogEntry,
} from "../../lib/moorhen-scene-resolver";
import type { SceneLiftHints } from "../../lib/moorhen-scene-lifter";

const DEFAULT_PLACEHOLDER = `# Type or paste a scene YAML here, or click "Capture" to
# initialise from the current Moorhen view.
#
# scene: my-scheme
# version: 1
# files:
#   - name: thing
#     fileId: 123
#     projectId: <project-uuid>
# domains:
#   - { name: domain-A, chain: A, range: 1-100, color: "#4b8bbe" }
# elements:
#   - file: thing
#     representations:
#       - { style: CRs, selection: "//A", colour: by-domain }
# resolver:
#   onMissingResidues: clamp-and-log
`;

const LIVE_APPLY_DEBOUNCE_MS = 750;

/**
 * Bundle asset map — keys are paths inside the .scene.zip (e.g.
 * "assets/complex.pdb"); values are the raw bytes. The Scenes panel
 * builds this when a .zip is opened and hands it to the apply callback
 * so the resolver's bundle: fetcher can look entries up. Lifter
 * populates a fresh map when capturing a scene that has bundle: refs.
 */
export type SceneBundleAssets = Map<string, ArrayBuffer>;

interface MoorhenScenesPanelProps {
  /** Apply YAML text to the current view. The second argument is the
   *  bundle-asset map (empty when not opened from a .scene.zip). */
  onApplyScene: (
    yamlText: string,
    assets: SceneBundleAssets,
  ) => Promise<SceneResolveResult>;
  /** Capture the current view as a scene + lift hints + bundle assets.
   *  The lifter may inline non-portable molecule coords/dicts into the
   *  assets map; the panel will save as .scene.zip when assets are
   *  present, plain .scene.yaml otherwise. */
  onCaptureScene: () => Promise<{
    scene: MoorhenScene;
    hints: SceneLiftHints;
    assets: SceneBundleAssets;
  }>;
  /** True once Moorhen / Coot is ready. Disables apply until then. */
  enabled: boolean;
}

type Severity = "info" | "success" | "warning" | "error";

interface PanelMessage {
  severity: Severity;
  text: string;
  log?: SceneResolveLogEntry[];
}

export const MoorhenScenesPanel: React.FC<MoorhenScenesPanelProps> = ({
  onApplyScene,
  onCaptureScene,
  enabled,
}) => {
  const [yamlText, setYamlText] = useState<string>("");
  // Live apply is opt-in. Loading / capturing a scene populates the editor
  // but does not touch the viewer until the user clicks Apply (or enables
  // Live apply). This lets the user inspect/edit a YAML before committing.
  const [liveApply, setLiveApply] = useState<boolean>(false);
  const [message, setMessage] = useState<PanelMessage | null>(null);
  const fileInputRef = useRef<HTMLInputElement | null>(null);
  // Bundle assets keyed by their path inside the .scene.zip. Populated
  // by Open when the user picks a .zip; consumed by Apply to satisfy
  // bundle: refs without any network or filesystem access.
  const assetsRef = useRef<SceneBundleAssets>(new Map());

  // --- actions ------------------------------------------------------------

  const handleCapture = useCallback(async () => {
    try {
      const { scene, hints, assets } = await onCaptureScene();
      const yaml = serialiseSceneWithComments(scene, hints.fileComments);
      setYamlText(yaml);
      // Take ownership of any bundle assets the lifter produced — Save
      // will bundle them up if the YAML actually references them.
      assetsRef.current = assets;
      const n = assets.size;
      setMessage({
        severity: "info",
        text:
          n > 0
            ? `Captured current view (${n} asset${n === 1 ? "" : "s"} ready for bundling)`
            : "Captured current view",
      });
    } catch (err) {
      console.error("Failed to capture scene:", err);
      setMessage({
        severity: "error",
        text: `Capture failed: ${err instanceof Error ? err.message : "unknown error"}`,
      });
    }
  }, [onCaptureScene]);

  const handleOpenClick = useCallback(() => {
    fileInputRef.current?.click();
  }, []);

  const handleOpenFile = useCallback(
    async (e: React.ChangeEvent<HTMLInputElement>) => {
      const file = e.target.files?.[0];
      e.target.value = "";
      if (!file) return;
      const isZip =
        file.name.endsWith(".zip") ||
        file.name.endsWith(".scene.zip") ||
        file.type === "application/zip";
      try {
        if (isZip) {
          // Unpack: scene.yaml at the root (required), everything else
          // becomes an asset keyed by its relative path.
          const buffer = await file.arrayBuffer();
          const zip = await JSZip.loadAsync(buffer);
          const yamlEntry = zip.file("scene.yaml") || zip.file("scene.yml");
          if (!yamlEntry) {
            setMessage({
              severity: "error",
              text: `${file.name} is missing scene.yaml at the root`,
            });
            return;
          }
          const yaml = await yamlEntry.async("string");
          const assets: SceneBundleAssets = new Map();
          // Iterate every entry; skip directories and the scene.yaml itself.
          const promises: Promise<void>[] = [];
          zip.forEach((relPath, entry) => {
            if (entry.dir) return;
            if (relPath === "scene.yaml" || relPath === "scene.yml") return;
            promises.push(
              entry.async("arraybuffer").then((buf) => {
                assets.set(relPath, buf);
              }),
            );
          });
          await Promise.all(promises);
          assetsRef.current = assets;
          setYamlText(yaml);
          setMessage({
            severity: "info",
            text: `Loaded ${file.name} (${assets.size} bundled asset${assets.size === 1 ? "" : "s"})`,
          });
        } else {
          const text = await file.text();
          // Opening a plain yaml clears any previously-loaded bundle
          // assets — they're tied to the zip they came from.
          assetsRef.current = new Map();
          setYamlText(text);
          setMessage({ severity: "info", text: `Loaded ${file.name}` });
        }
      } catch (err) {
        setMessage({
          severity: "error",
          text: `Could not read ${file.name}: ${(err as Error).message}`,
        });
      }
    },
    [],
  );

  const applyNow = useCallback(
    async (text: string, opts: { silentOnParseError: boolean }) => {
      try {
        parseScene(text); // validate first for a clean error
      } catch (err) {
        if (opts.silentOnParseError) return;
        if (err instanceof SceneParseError) {
          const first = err.errors[0];
          setMessage({
            severity: "error",
            text: `Invalid scene${first.path ? ` at ${first.path}` : ""}: ${first.message}`,
          });
        } else {
          setMessage({
            severity: "error",
            text: `Invalid scene: ${(err as Error).message}`,
          });
        }
        return;
      }

      try {
        const result = await onApplyScene(text, assetsRef.current);
        const repTotal = result.applied.reduce((sum, a) => sum + a.representations, 0);
        const severity: Severity =
          result.unresolvedFiles.length > 0 || repTotal === 0 ? "warning" : "success";
        const parts: string[] = [
          `Applied ${repTotal} representation${repTotal === 1 ? "" : "s"} across ${result.applied.length} file${result.applied.length === 1 ? "" : "s"}`,
        ];
        if (result.unresolvedFiles.length > 0) {
          parts.push(`unresolved: ${result.unresolvedFiles.join(", ")}`);
        }
        setMessage({ severity, text: parts.join(" — "), log: result.log });
      } catch (err) {
        console.error("Failed to apply scene:", err);
        setMessage({
          severity: "error",
          text: `Apply failed: ${err instanceof Error ? err.message : "unknown error"}`,
        });
      }
    },
    [onApplyScene],
  );

  const handleApplyNow = useCallback(() => {
    if (!yamlText.trim()) {
      setMessage({ severity: "warning", text: "Editor is empty" });
      return;
    }
    void applyNow(yamlText, { silentOnParseError: false });
  }, [applyNow, yamlText]);

  const handleSave = useCallback(async () => {
    if (!yamlText.trim()) {
      setMessage({ severity: "warning", text: "Nothing to save" });
      return;
    }
    let sceneName = "scene";
    let hasBundleRefs = false;
    try {
      const parsed = parseScene(yamlText);
      sceneName = parsed.scene || "scene";
      hasBundleRefs = (parsed.files ?? []).some((f) => !!f.bundle);
    } catch {
      // Save the raw text anyway so the user doesn't lose work in progress.
    }

    const safe = sceneName.replace(/[^A-Za-z0-9._-]+/g, "_");
    const a = document.createElement("a");
    document.body.appendChild(a);

    // Save as .scene.zip when the scene has any bundle: refs AND we
    // have asset bytes for them — otherwise the zip would be missing
    // its own attachments. Falls back to plain yaml if there's nothing
    // to attach.
    const shouldBundle = hasBundleRefs && assetsRef.current.size > 0;
    if (shouldBundle) {
      const zip = new JSZip();
      zip.file("scene.yaml", yamlText);
      for (const [relPath, buf] of assetsRef.current.entries()) {
        zip.file(relPath, buf);
      }
      const blob = await zip.generateAsync({ type: "blob" });
      const url = URL.createObjectURL(blob);
      a.href = url;
      a.download = `${safe}.scene.zip`;
      a.click();
      URL.revokeObjectURL(url);
    } else {
      const blob = new Blob([yamlText], { type: "application/x-yaml;charset=utf-8" });
      const url = URL.createObjectURL(blob);
      a.href = url;
      a.download = `${safe}.scene.yaml`;
      a.click();
      URL.revokeObjectURL(url);
    }
    document.body.removeChild(a);
    setMessage({ severity: "info", text: `Saved ${a.download}` });
  }, [yamlText]);

  // --- live apply ---------------------------------------------------------

  // Re-apply on edit when liveApply is on. Silent on parse errors so
  // intermediate broken yaml doesn't spam the message panel; the user's
  // next valid keystroke will re-apply cleanly.
  useEffect(() => {
    if (!liveApply || !enabled || !yamlText.trim()) return;
    const handle = setTimeout(() => {
      void applyNow(yamlText, { silentOnParseError: true });
    }, LIVE_APPLY_DEBOUNCE_MS);
    return () => clearTimeout(handle);
  }, [liveApply, enabled, yamlText, applyNow]);

  // --- rendering ----------------------------------------------------------

  const messageColor = useMemo(() => {
    if (!message) return "transparent";
    return ({
      info: "#e3f2fd",
      success: "#e8f5e9",
      warning: "#fff8e1",
      error: "#ffebee",
    } as const)[message.severity];
  }, [message]);

  return (
    <Stack
      direction="column"
      spacing={1}
      sx={{ p: 1, height: "100%", overflow: "hidden" }}
    >
      {/* Toolbar */}
      <Stack direction="row" spacing={0.5} sx={{ flexWrap: "wrap", alignItems: "center" }}>
        <Tooltip title="Capture the current view into the editor (lifter)">
          <span>
            <Button
              size="small"
              variant="outlined"
              startIcon={<CaptureIcon />}
              onClick={handleCapture}
              sx={{ fontSize: "0.75rem", textTransform: "none" }}
            >
              Capture
            </Button>
          </span>
        </Tooltip>

        <Tooltip title="Open a .scene.yaml or .scene.zip from disk (zips bring their bundled assets along)">
          <span>
            <Button
              size="small"
              variant="outlined"
              startIcon={<OpenIcon />}
              onClick={handleOpenClick}
              sx={{ fontSize: "0.75rem", textTransform: "none" }}
            >
              Open…
            </Button>
          </span>
        </Tooltip>

        <Tooltip title="Apply the editor's current YAML now">
          <span>
            <Button
              size="small"
              variant="contained"
              startIcon={<ApplyIcon />}
              onClick={handleApplyNow}
              disabled={!enabled}
              sx={{ fontSize: "0.75rem", textTransform: "none" }}
            >
              Apply
            </Button>
          </span>
        </Tooltip>

        <Tooltip title="Save the editor's YAML to disk">
          <span>
            <IconButton size="small" onClick={handleSave}>
              <SaveIcon fontSize="small" />
            </IconButton>
          </span>
        </Tooltip>

        <FormControlLabel
          control={
            <Checkbox
              size="small"
              checked={liveApply}
              onChange={(e) => setLiveApply(e.target.checked)}
            />
          }
          label={<Typography variant="caption">Live apply</Typography>}
          sx={{ ml: "auto", mr: 0 }}
        />
      </Stack>

      {/* Editor */}
      <Box sx={{ flex: "1 1 auto", minHeight: 200, border: "1px solid #ddd" }}>
        <Editor
          height="100%"
          value={yamlText}
          language="yaml"
          theme="vs"
          options={{
            minimap: { enabled: false },
            scrollBeyondLastLine: false,
            fontSize: 12,
            wordWrap: "on",
            tabSize: 2,
            renderLineHighlight: "none",
            scrollbar: { vertical: "auto", horizontal: "auto" },
            placeholder: DEFAULT_PLACEHOLDER,
          }}
          onChange={(v) => setYamlText(v ?? "")}
        />
      </Box>

      {/* Message + resolver log */}
      <Box
        sx={{
          flex: "0 0 auto",
          maxHeight: 160,
          overflow: "auto",
          p: 1,
          fontSize: "0.75rem",
          backgroundColor: messageColor,
          border: "1px solid #eee",
          borderRadius: 1,
        }}
      >
        {message ? (
          <>
            <Typography variant="caption" sx={{ display: "block", fontWeight: 600 }}>
              {message.text}
            </Typography>
            {message.log && message.log.length > 0 && (
              <Box component="ul" sx={{ m: 0, pl: 2, mt: 0.5 }}>
                {message.log.map((entry, i) => (
                  <Typography
                    key={i}
                    component="li"
                    variant="caption"
                    sx={{ display: "list-item" }}
                  >
                    <strong>{entry.domain}</strong> ({entry.file}): {entry.message}
                  </Typography>
                ))}
              </Box>
            )}
          </>
        ) : (
          <Typography variant="caption" sx={{ color: "#666" }}>
            No messages — edit YAML above and (with Live apply on) the view
            will update automatically when it parses.
          </Typography>
        )}
      </Box>

      <input
        ref={fileInputRef}
        type="file"
        accept=".yaml,.yml,.scene.yaml,.zip,.scene.zip,text/yaml,application/x-yaml,application/zip"
        style={{ display: "none" }}
        onChange={handleOpenFile}
      />
    </Stack>
  );
};
