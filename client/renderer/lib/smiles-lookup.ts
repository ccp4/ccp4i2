/**
 * SMILES for compound registration numbers, from the compounds API.
 *
 * A campaign table asked for every member's compound in one request:
 * `reg_number__in=<1474 numbers>` is a URL of some 10 KB, past what proxies
 * and servers take, and a paginated answer delivers only its first page --
 * so on a big campaign the sketches never came. Instead the numbers go in
 * batches, each batch following the API's pages, batches a few at a time;
 * each batch's SMILES show as it arrives. What has been looked up is kept
 * (a compound's SMILES does not change), so the table's polling asks for
 * nothing twice.
 */
import { useEffect, useMemo, useState } from "react";
import { apiFetch } from "../api-fetch";

export const COMPOUNDS_ENDPOINT = "/api/proxy/compounds/compounds/";
/** Registration numbers per request: a short URL, one page of results. */
export const SMILES_BATCH_SIZE = 100;
/** Requests in flight at once. */
const CONCURRENCY = 4;
/** Pages followed per batch, a guard against a `next` that never ends. */
const MAX_PAGES = 20;

type FetchJson = (url: string) => Promise<any>;

interface Compound {
  reg_number?: number;
  smiles?: string;
}

const defaultFetchJson: FetchJson = async (url) => {
  const response = await apiFetch(url);
  if (!response.ok) throw new Error(`compounds lookup failed: HTTP ${response.status}`);
  return response.json();
};

/** One batch of registration numbers: their SMILES, following pages. Each
 *  number asked for gets an entry, null when the API has none. */
export async function fetchSmilesBatch(
  regIds: number[],
  fetchJson: FetchJson = defaultFetchJson,
): Promise<Map<number, string | null>> {
  const found = new Map<number, string | null>(regIds.map((id) => [id, null]));
  let url: string | null =
    `${COMPOUNDS_ENDPOINT}?reg_number__in=${regIds.join(",")}&page_size=${regIds.length}`;
  for (let page = 0; url && page < MAX_PAGES; page++) {
    const json = await fetchJson(url);
    const results: Compound[] = Array.isArray(json) ? json : json?.results ?? [];
    for (const compound of results) {
      if (compound.reg_number && compound.smiles) found.set(compound.reg_number, compound.smiles);
    }
    // `next` is the backend's own absolute URL: keep its query, ask the proxy.
    const next: string | null = Array.isArray(json) ? null : json?.next ?? null;
    url = next ? `${COMPOUNDS_ENDPOINT}${new URL(next, "http://x").search}` : null;
  }
  return found;
}

/** Looked-up SMILES by registration number (null: none). Module-wide, so a
 *  table that remounts or polls does not ask again. */
const cache = new Map<number, string | null>();
const inFlight = new Set<number>();

/**
 * SMILES for ``regIds``, filling in batch by batch.
 */
export function useSmilesLookup(regIds: number[]) {
  const key = useMemo(() => [...new Set(regIds)].sort((a, b) => a - b).join(","), [regIds]);
  const [version, setVersion] = useState(0);
  const [error, setError] = useState<Error | null>(null);
  const [pending, setPending] = useState(0);

  useEffect(() => {
    if (!key) return;
    const missing = key
      .split(",")
      .map(Number)
      .filter((id) => !cache.has(id) && !inFlight.has(id));
    if (missing.length === 0) return;
    // No cancellation: a batch outlives a change of key (the member list
    // polls), and whatever it brings is wanted by the next render too.
    const batches: number[][] = [];
    for (let i = 0; i < missing.length; i += SMILES_BATCH_SIZE) {
      batches.push(missing.slice(i, i + SMILES_BATCH_SIZE));
    }
    missing.forEach((id) => inFlight.add(id));
    setPending((n) => n + batches.length);
    let nextBatch = 0;
    const worker = async () => {
      while (nextBatch < batches.length) {
        const batch = batches[nextBatch++];
        try {
          const found = await fetchSmilesBatch(batch);
          found.forEach((smiles, id) => cache.set(id, smiles));
        } catch (err) {
          // Left uncached, so a later render asks again.
          setError(err instanceof Error ? err : new Error(String(err)));
        } finally {
          batch.forEach((id) => inFlight.delete(id));
          setPending((n) => n - 1);
          setVersion((v) => v + 1);
        }
      }
    };
    for (let i = 0; i < Math.min(CONCURRENCY, batches.length); i++) void worker();
  }, [key]);

  const smilesMap = useMemo(() => {
    const map: Record<number, string> = {};
    if (!key) return map;
    for (const id of key.split(",").map(Number)) {
      const smiles = cache.get(id);
      if (smiles) map[id] = smiles;
    }
    return map;
    // version: recompute as batches land
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [key, version]);

  return { smilesMap, isLoading: pending > 0, error };
}
