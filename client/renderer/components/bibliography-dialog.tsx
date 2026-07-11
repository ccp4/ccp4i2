import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  Button,
  List,
  ListItem,
  Link,
  Typography,
  Box,
  CircularProgress,
  Stack,
} from "@mui/material";
import {
  MenuBook,
  ContentCopy,
  Description,
  Code,
} from "@mui/icons-material";
import React, { useMemo, useState } from "react";
import { useApi } from "../api";

export interface BibReference {
  pmid?: string | null;
  title?: string | null;
  authors?: string[];
  source?: string | null;
  link?: string | null;
  cited_by?: string;
}

interface BibliographyResponse {
  // get_endpoint does NOT unwrap the api_success {success, data} envelope, so
  // the payload arrives as either the wrapped or the unwrapped shape — accept
  // both (same gotcha the Export MTZ menu hit).
  success?: boolean;
  data?: { result?: BibReference[] };
  result?: BibReference[];
}

/** Best-effort year extraction from a MedLine source string, e.g.
 * "Acta Crystallogr D Biol Crystallogr. 2011;67(Pt 4):235-242." -> "2011". */
function extractYear(source?: string | null): string | undefined {
  if (!source) return undefined;
  const m = source.match(/\b(19|20)\d{2}\b/);
  return m ? m[0] : undefined;
}

/** DOI from a link like https://doi.org/10.1107/... */
function extractDoi(link?: string | null): string | undefined {
  if (!link) return undefined;
  const m = link.match(/10\.\d{4,}\/\S+/);
  return m ? m[0] : undefined;
}

function referencesToText(refs: BibReference[]): string {
  return refs
    .map((r) => {
      const parts: string[] = [];
      if (r.authors && r.authors.length) parts.push(r.authors.join(", "));
      if (r.title) parts.push(r.title);
      if (r.source) parts.push(r.source);
      if (r.link) parts.push(r.link);
      return parts.join("\n");
    })
    .join("\n\n");
}

/** Generate a BibTeX @article entry per reference. Fields are best-effort from
 * the MedLine data (author/title/journal-source/year/doi/url). */
function referencesToBibtex(refs: BibReference[]): string {
  const usedKeys = new Set<string>();
  return refs
    .map((r, i) => {
      // cite key: first-author surname + year, deduped.
      const firstAuthor = (r.authors?.[0] || "anon").split(/[\s,]+/)[0];
      const year = extractYear(r.source) || "";
      let key = `${firstAuthor}${year}`.replace(/[^A-Za-z0-9]/g, "");
      if (!key) key = `ref${i + 1}`;
      let uniqueKey = key;
      let suffix = 0;
      while (usedKeys.has(uniqueKey)) {
        suffix += 1;
        uniqueKey = `${key}${String.fromCharCode(96 + suffix)}`;
      }
      usedKeys.add(uniqueKey);

      const fields: string[] = [];
      if (r.authors && r.authors.length) {
        fields.push(`  author = {${r.authors.join(" and ")}}`);
      }
      if (r.title) fields.push(`  title = {${r.title.replace(/\.$/, "")}}`);
      // The source string as the journal citation (unparsed tail kept in note).
      if (r.source) fields.push(`  journal = {${r.source.replace(/\.$/, "")}}`);
      if (year) fields.push(`  year = {${year}}`);
      const doi = extractDoi(r.link);
      if (doi) fields.push(`  doi = {${doi}}`);
      if (r.link) fields.push(`  url = {${r.link}}`);
      if (r.pmid) fields.push(`  pmid = {${r.pmid}}`);

      return `@article{${uniqueKey},\n${fields.join(",\n")}\n}`;
    })
    .join("\n\n");
}

function downloadTextFile(content: string, filename: string): void {
  const blob = new Blob([content], { type: "text/plain;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

interface BibliographyDialogProps {
  /** When set, fetch + show the bibliography for this job. null = closed. */
  jobId: number | null;
  /** When set (and jobId null), show the whole-project bibliography. */
  projectId?: number | null;
  onClose: () => void;
}

/**
 * Modal listing the compiled bibliography for a job (its task + subjobs) or a
 * project. Presentation-only: the server compiles and dedupes the references.
 */
export function BibliographyDialog({
  jobId,
  projectId,
  onClose,
}: BibliographyDialogProps) {
  const api = useApi();
  const open = jobId != null || projectId != null;

  const endpoint = useMemo(() => {
    if (jobId != null) {
      return { type: "jobs" as const, id: jobId, endpoint: "bibliography" };
    }
    if (projectId != null) {
      return {
        type: "projects" as const,
        id: projectId,
        endpoint: "bibliography",
      };
    }
    return null;
  }, [jobId, projectId]);

  const { data, isLoading } = api.get_endpoint<BibliographyResponse>(endpoint);
  const references: BibReference[] = data?.data?.result ?? data?.result ?? [];

  const scopeLabel = jobId != null ? "job" : "project";
  const filenameStem = `bibliography_${scopeLabel}_${jobId ?? projectId ?? ""}`;
  const [copied, setCopied] = useState(false);

  const handleCopy = async () => {
    try {
      await navigator.clipboard.writeText(referencesToText(references));
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    } catch {
      // clipboard unavailable (e.g. insecure context) — fall back to text file
      downloadTextFile(referencesToText(references), `${filenameStem}.txt`);
    }
  };
  const handleDownloadText = () =>
    downloadTextFile(referencesToText(references), `${filenameStem}.txt`);
  const handleDownloadBibtex = () =>
    downloadTextFile(referencesToBibtex(references), `${filenameStem}.bib`);

  const hasRefs = references.length > 0;

  return (
    <Dialog open={open} onClose={onClose} maxWidth="md" fullWidth>
      <DialogTitle>
        <Stack direction="row" spacing={1} alignItems="center">
          <MenuBook fontSize="small" />
          <span>Bibliography ({scopeLabel})</span>
        </Stack>
      </DialogTitle>
      <DialogContent dividers>
        {isLoading ? (
          <Box sx={{ display: "flex", justifyContent: "center", py: 4 }}>
            <CircularProgress />
          </Box>
        ) : references.length === 0 ? (
          <Typography color="text.secondary" sx={{ py: 2 }}>
            No references found for this {scopeLabel}.
          </Typography>
        ) : (
          <List>
            {references.map((ref, i) => (
              <ListItem
                key={ref.pmid ?? `${ref.title}-${i}`}
                alignItems="flex-start"
                sx={{ display: "block", py: 1 }}
              >
                <Typography variant="body1" sx={{ fontWeight: 500 }}>
                  {ref.link ? (
                    <Link href={ref.link} target="_blank" rel="noopener">
                      {ref.title || ref.link}
                    </Link>
                  ) : (
                    ref.title || "(untitled reference)"
                  )}
                </Typography>
                {ref.authors && ref.authors.length > 0 && (
                  <Typography variant="body2" color="text.secondary">
                    {ref.authors.join(", ")}
                  </Typography>
                )}
                {ref.source && (
                  <Typography variant="body2" color="text.secondary">
                    <em>{ref.source}</em>
                  </Typography>
                )}
              </ListItem>
            ))}
          </List>
        )}
      </DialogContent>
      <DialogActions sx={{ justifyContent: "space-between", px: 3 }}>
        <Stack direction="row" spacing={1}>
          <Button
            size="small"
            startIcon={<ContentCopy />}
            onClick={handleCopy}
            disabled={!hasRefs}
          >
            {copied ? "Copied!" : "Copy"}
          </Button>
          <Button
            size="small"
            startIcon={<Description />}
            onClick={handleDownloadText}
            disabled={!hasRefs}
          >
            Text
          </Button>
          <Button
            size="small"
            startIcon={<Code />}
            onClick={handleDownloadBibtex}
            disabled={!hasRefs}
          >
            BibTeX
          </Button>
        </Stack>
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
}
