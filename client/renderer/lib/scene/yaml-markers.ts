/**
 * Map Zod-contract validation errors onto Monaco editor markers (red
 * squiggles), located at the offending node in the YAML source.
 *
 * Isolated, dependency-light: uses the `yaml` library's Document + LineCounter
 * (already a project dep) to resolve a dotted error path to a source range,
 * and the project's own `validateScene`. No monaco-editor import — the caller
 * passes its monaco severity constant — so this stays free of the global
 * Monaco loader/worker setup that `monaco-yaml` would require.
 */
import { parseDocument, LineCounter } from "yaml";

import { validateScene } from "./index";

export interface SceneMarker {
  startLineNumber: number;
  startColumn: number;
  endLineNumber: number;
  endColumn: number;
  message: string;
  severity: number;
}

/** "elements[0].representations[1].style" → ["elements", 0, "representations", 1, "style"] */
export function pathToSegments(path: string): (string | number)[] {
  const segs: (string | number)[] = [];
  if (!path) return segs;
  for (const part of path.split(".")) {
    const key = part.replace(/\[\d+\]/g, "");
    if (key) segs.push(key);
    const idxs = part.match(/\[(\d+)\]/g);
    if (idxs) for (const i of idxs) segs.push(Number(i.slice(1, -1)));
  }
  return segs;
}

/**
 * Produce Monaco markers for a YAML scene string. `errorSeverity` is the
 * caller's `monaco.MarkerSeverity.Error`. Returns [] when the YAML itself is
 * unparseable (the editor's own YAML tokeniser / the apply toast handle that).
 */
export function sceneMarkers(text: string, errorSeverity: number): SceneMarker[] {
  const lineCounter = new LineCounter();
  let doc: ReturnType<typeof parseDocument>;
  let raw: unknown;
  try {
    doc = parseDocument(text, { lineCounter });
    raw = doc.toJS();
  } catch {
    return [];
  }
  const { errors } = validateScene(raw);
  return errors.map((e) => {
    const segs = pathToSegments(e.path);
    const node = segs.length ? (doc.getIn(segs, true) as unknown) : null;
    const range = (node as { range?: [number, number, number] } | null)?.range;
    const startOff = range ? range[0] : 0;
    const endOff = range ? range[1] : Math.min(text.length || 1, 1);
    const start = lineCounter.linePos(startOff);
    const end = lineCounter.linePos(endOff);
    return {
      startLineNumber: start.line,
      startColumn: start.col,
      endLineNumber: end.line,
      endColumn: end.col,
      message: e.message,
      severity: errorSeverity,
    };
  });
}
