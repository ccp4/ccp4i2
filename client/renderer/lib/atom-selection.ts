/**
 * The chains a selection names, when it names nothing but whole chains
 * ("A", "A/", "{A or B/}"), so the atom-selection builder's chain boxes can
 * start from a job's existing selection. null for anything else, or for a
 * chain the file does not have: the boxes then start with every chain, and
 * the job's selection is shown as typed.
 */
export function chainsFromSelection(
  selection: string | undefined | null,
  chains: string[]
): Set<string> | null {
  let text = String(selection ?? "").trim();
  if (!text) return null;
  if (text.startsWith("{") && text.endsWith("}")) text = text.slice(1, -1).trim();
  const names: string[] = [];
  for (const part of text.split(/\s+or\s+/i)) {
    const m = part.trim().match(/^([A-Za-z0-9]+)\/?$/);
    if (!m || !chains.includes(m[1])) return null;
    names.push(m[1]);
  }
  return new Set(names);
}
