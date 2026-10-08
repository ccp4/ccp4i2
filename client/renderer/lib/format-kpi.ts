/**
 * A KPI value for a chip. Every numeric KPI reaches the client through the
 * float table, counts included (the gleaner stores a CInt as a float), so an
 * integral value is printed as the integer it is; anything else to three
 * significant figures.
 */
export function formatKpiValue(value: number): string {
  if (!Number.isFinite(value)) return String(value);
  if (Number.isInteger(value)) return String(value);
  return value.toPrecision(3);
}

// A lower-case letter or digit followed by a capital; a run of capitals
// followed by a capitalised word (the run is an acronym: "XMLFile").
const WORD_BOUNDARY = /(?<=[a-z0-9])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])/g;

/**
 * A KPI key in words, for a key the server sent no label for:
 * "highResLimit" -> "High res limit", "nEvents" -> "N events". Words in
 * capitals are symbols and kept; others after the first are lower-cased.
 * Must agree with `humanise_key` in server/ccp4i2/lib/kpi_labels.py.
 */
export function humaniseKpiKey(key: string): string {
  const words = key
    .replace(/_/g, " ")
    .replace(WORD_BOUNDARY, " ")
    .split(/\s+/)
    .filter(Boolean)
    .map((word) => {
      const rest = word.slice(1);
      const allCaps =
        word === word.toUpperCase() && word !== word.toLowerCase();
      return rest === rest.toLowerCase() && !allCaps ? word.toLowerCase() : word;
    });
  if (words.length === 0) return key;
  words[0] = words[0].charAt(0).toUpperCase() + words[0].slice(1);
  return words.join(" ");
}

/**
 * What a KPI is called on screen (#596). The key is a field name
 * ("highResLimit") and is never shown: the server's label for it (job_tree
 * kpis.labels, from server/ccp4i2/lib/kpi_labels.py) is, or failing that the
 * key in words.
 */
export function kpiLabel(
  key: string,
  labels?: Record<string, string> | null
): string {
  return labels?.[key] || humaniseKpiKey(key);
}
