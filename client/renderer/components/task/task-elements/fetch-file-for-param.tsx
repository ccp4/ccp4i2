/*
 * Copyright (C) 2025-2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import {
  Autocomplete,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  TextField,
  Typography,
} from "@mui/material";
import { useCallback, useEffect, useMemo, useState } from "react";
import { apiFetch, apiBlob, apiText } from "../../../api-fetch";
import { useCCP4i2Window } from "../../../app-context";
import { useJob, useProject } from "../../../utils";

// Note: cootModule is no longer used here - MTZ parsing now uses native TypeScript parser
import { usePopcorn } from "../../../providers/popcorn-provider";
import { useTaskInterface } from "../../../providers/task-provider";
import {
  parseMultiChainFasta,
  ChainSequenceInfo,
} from "./mmcif-sequence-parser";
import { ChainPickerDialog } from "./chain-picker-dialog";

/** Human-readable labels for download modes shown in the dropdown */
const MODE_LABELS: Record<string, string> = {
  ebiPdb: "PDBe (mmCIF coordinates)",
  rcsbPdb: "RCSB PDB (mmCIF coordinates)",
  uniprotAFPdb: "AlphaFold DB (predicted model)",
  ebiSFs: "PDBe (structure factors)",
  uniprotFasta: "UniProt (FASTA sequence)",
  "Uppsala-EDS": "PDB-REDO (map coefficients)",
};

/** Placeholder text for the accession code field per mode */
const MODE_PLACEHOLDERS: Record<string, string> = {
  ebiPdb: "e.g. 1cbs",
  rcsbPdb: "e.g. 1cbs",
  uniprotAFPdb: "e.g. P07550",
  ebiSFs: "e.g. 1cbs",
  uniprotFasta: "e.g. P07550 or CDK2_HUMAN",
  "Uppsala-EDS": "e.g. 1cbs",
};

/** Modes that fetch PDB coordinate files */
const COORD_MODES = new Set(["ebiPdb", "rcsbPdb"]);

/** Mime types that indicate a sequence parameter target */
const SEQ_MIME_TYPES = new Set([
  "application/CCP4-seq",
  "application/CCP4-seqalign",
]);

/** Check if the target item is a sequence parameter */
function isSequenceTarget(item: any, modes: string[] | null): boolean {
  // Check if the qualifiers indicate a sequence mime type
  const mimeType = item?._qualifiers?.mimeTypeName;
  if (mimeType && SEQ_MIME_TYPES.has(mimeType)) return true;
  // Check class name
  const cls = item?._class;
  if (cls === "CSeqDataFile" || cls === "CSeqAlignDataFile") return true;
  // Check if modes were overridden to include seq-related modes alongside coord modes
  // (e.g. casucontentseq passes ["uniprotFasta", "ebiPdb", "rcsbPdb"])
  if (modes && modes.some((m) => m === "uniprotFasta") && modes.some((m) => COORD_MODES.has(m))) {
    return true;
  }
  return false;
}

/** Build a single-chain FASTA string */
function buildFasta(pdbId: string, chain: ChainSequenceInfo): string {
  return `>${pdbId.toUpperCase()}_${chain.chainId} ${chain.polymerType} ${chain.length} residues\n${chain.sequence}\n`;
}

interface FetchFileForParamProps {
  open: boolean;
  onClose: () => void;
  onSuccess?: (updatedItem: any) => void;
}
export const FetchFileForParam: React.FC<FetchFileForParamProps> = ({
  open,
  onClose,
}) => {
  const { setMessage } = usePopcorn();
  const { fetchItemParams, setFetchItemParams, setDownloadDialogOpen } =
    useTaskInterface();

  const { item, modes, onChange } = useMemo(() => {
    //alert(JSON.stringify(itemParams));
    return fetchItemParams
      ? fetchItemParams
      : { item: null, modes: null, onChange: null };
  }, [fetchItemParams]);

  const downloadModes: string[] = useMemo(
    () => modes || item?._qualifiers?.downloadModes || [],
    [item, modes]
  );

  const { jobId } = useCCP4i2Window();
  const { job, uploadFileParam } = useJob(jobId);
  const { mutateJobs, mutateFiles } = useProject(job?.project);

  const [mode, setMode] = useState<string | null>(null);

  //Initialise identifier to empty string
  useEffect(() => {
    if (open) setIdentifier("");
  }, [open]);

  useEffect(() => {
    if (modes && modes.length > 0 && (!mode || !modes.includes(mode))) {
      setMode(modes[0]);
    }
  }, [modes]);

  const [identifier, setIdentifier] = useState<string | null>(null);
  const [inFlight, setInFlight] = useState(false);

  // Chain picker state for multi-chain PDB → sequence disambiguation
  const [chainPickerOpen, setChainPickerOpen] = useState(false);
  const [pendingChains, setPendingChains] = useState<ChainSequenceInfo[]>([]);
  const [pendingPdbId, setPendingPdbId] = useState("");

  const uploadFile = useCallback(
    async (fileBlob: Blob, fileName: string) => {
      if (job && item) {
        setMessage(`Uploading file ${fileName} for ${item._objectPath}`);

        try {
          // Use centralized uploadFileParam with local cache patching
          const uploadResult = await uploadFileParam({
            objectPath: item._objectPath,
            file: fileBlob,
            fileName,
          });

          setMessage(`File ${fileName} uploaded for ${item._objectPath}`);
          if (uploadResult?.success && uploadResult.data?.updated_item) {
            if (onChange) {
              onChange(uploadResult.data.updated_item);
            }
            // Additional mutations not handled by uploadFileParam
            mutateJobs();
            mutateFiles();
          }
        } catch (err: any) {
          setMessage(
            err.message || `Failed to upload ${fileName}`,
            "error"
          );
        }
      }
    },
    [item, job, uploadFileParam, onChange, mutateJobs, mutateFiles, setMessage]
  );

  /** Upload a selected chain as a FASTA file */
  const uploadChainAsFasta = useCallback(
    (chain: ChainSequenceInfo) => {
      const pdbId = pendingPdbId || identifier || "unknown";
      const fastaContent = buildFasta(pdbId, chain);
      const blob = new Blob([fastaContent], { type: "text/plain" });
      uploadFile(blob, `${pdbId.toUpperCase()}_${chain.chainId}.fasta`);
      onClose();
    },
    [pendingPdbId, identifier, uploadFile, onClose]
  );

  const handleEbiCoordFetch = useCallback(async () => {
    if (identifier) {
      try {
        const url = `https://www.ebi.ac.uk/pdbe/api/pdb/entry/files/${identifier.toLowerCase()}`;
        const result = await apiFetch(url);
        const data = await result.json();
        if (data && data[identifier.toLowerCase()]) {
          const file = data[identifier.toLowerCase()];
          let fetchURL = file.PDB.downloads
            .filter((item: any) => item.label === "Archive mmCIF file")
            .at(0).url;
          setMessage(`Fetching file from ${fetchURL}`);
          const content = await apiBlob(fetchURL);
          setMessage(`Fetched file from ${fetchURL}`);
          uploadFile(content, fetchURL.split("/").at(-1));
          onClose();
        } else {
          const errorText = await result.text();
          setMessage(errorText);
        }
      } catch (err: any) {
        setMessage(err.message || "Unknown error");
        console.error("FetchFileForParam handleFetch error", err);
        return;
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleEbiSFsFetch = useCallback(async () => {
    if (identifier) {
      try {
        const url = `https://www.ebi.ac.uk/pdbe/api/pdb/entry/files/${identifier.toLowerCase()}`;
        const result = await apiFetch(url);
        const data = await result.json();
        if (data && data[identifier.toLowerCase()]) {
          const file = data[identifier.toLowerCase()];
          let fetchURL = file.PDB.downloads
            .filter((item: any) => item.label === "Structure Factors")
            .at(0).url;
          setMessage(`Fetching file from ${fetchURL}`);
          const content = await apiBlob(fetchURL);
          setMessage(`Fetched file from ${fetchURL}`);
          uploadFile(content, fetchURL.split("/").at(-1));
          onClose();
        } else {
          const errorText = await result.text();
          setMessage(errorText);
        }
      } catch (err: any) {
        setMessage(err.message || "Unknown error", "error");
        console.error("FetchFileForParam handleEbiSFsFetch error", err);
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleUniprotFastaFetch = useCallback(async () => {
    if (identifier) {
      setMessage(`Fetching FASTA file for ${identifier.toUpperCase()}`);
      const data = await apiText(`/api/proxy/uniprot/uniprotkb/${identifier.toUpperCase()}.fasta`);
      setMessage(`Fetched FASTA file for ${identifier.toUpperCase()}`);
      const content = new Blob([data], {
        type: "text/plain",
      });
      uploadFile(content, `${identifier.toUpperCase()}.fasta`);
      onClose();
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleRcsbPdbFetch = useCallback(async () => {
    if (identifier) {
      const id = identifier.toLowerCase();
      const fetchURL = `https://files.rcsb.org/download/${id}.cif`;
      setMessage(`Fetching coordinates from RCSB for ${id}`);
      try {
        const content = await apiBlob(fetchURL);
        setMessage(`Fetched coordinates from RCSB for ${id}`);
        uploadFile(content, `${id}.cif`);
        onClose();
      } catch (err: any) {
        setMessage(err.message || `Failed to fetch ${id} from RCSB`);
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleUniprotAFPdbFetch = useCallback(async () => {
    if (identifier) {
      const id = identifier.toUpperCase();
      const fetchURL = `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.cif`;
      setMessage(`Fetching AlphaFold model for ${id}`);
      try {
        const content = await apiBlob(fetchURL);
        setMessage(`Fetched AlphaFold model for ${id}`);
        uploadFile(content, `AF-${id}-F1-model_v4.cif`);
        onClose();
      } catch (err: any) {
        setMessage(err.message || `Failed to fetch AlphaFold model for ${id}`);
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleUppsalaEdsFetch = useCallback(async () => {
    // Uppsala EDS was retired in 2017 — fetch from PDB-REDO instead
    if (identifier) {
      const id = identifier.toLowerCase();
      const fetchURL = `https://pdb-redo.eu/db/${id}/${id}_final.mtz`;
      setMessage(`Fetching map coefficients from PDB-REDO for ${id}`);
      try {
        const content = await apiBlob(fetchURL);
        setMessage(`Fetched map coefficients from PDB-REDO for ${id}`);
        uploadFile(content, `${id}_final.mtz`);
        onClose();
      } catch (err: any) {
        setMessage(err.message || `Failed to fetch map coefficients for ${id}`);
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  /**
   * Extract chain sequences from PDBe molecules API JSON response.
   * Returns ChainSequenceInfo[] from the /api/pdb/entry/molecules/{id} endpoint.
   */
  const parsePdbeMolecules = useCallback((data: any, pdbId: string): ChainSequenceInfo[] => {
    const entities = data[pdbId.toLowerCase()];
    if (!Array.isArray(entities)) return [];
    const chains: ChainSequenceInfo[] = [];
    for (const entity of entities) {
      if (!entity.sequence || entity.molecule_type === "water" || entity.molecule_type === "bound") continue;
      const polymerType = entity.molecule_type?.includes("polypeptide") ? "PROTEIN"
        : entity.molecule_type?.includes("polyribonucleotide") && !entity.molecule_type?.includes("deoxy") ? "RNA"
        : entity.molecule_type?.includes("polydeoxyribonucleotide") ? "DNA"
        : "OTHER" as const;
      const chainIds: string[] = entity.in_chains || [];
      const moleculeName = Array.isArray(entity.molecule_name) ? entity.molecule_name[0] : entity.molecule_name || "";
      for (const chainId of chainIds) {
        chains.push({
          chainId,
          sequence: entity.sequence,
          polymerType,
          length: entity.sequence.length,
          description: `${moleculeName} (Chain ${chainId})`,
        });
      }
    }
    return chains;
  }, []);

  /**
   * Fetch sequences for a PDB entry and extract chain information.
   * - EBI: uses /api/pdb/entry/molecules/{id} JSON endpoint
   * - RCSB: uses /fasta/entry/{ID} FASTA endpoint
   * If single chain, upload directly; if multiple, show chain picker.
   */
  const handlePdbSequenceFetch = useCallback(async (source: "ebi" | "rcsb") => {
    if (!identifier) return;
    const id = identifier.toLowerCase();
    try {
      setMessage(`Fetching sequences for ${id.toUpperCase()}...`);
      let chains: ChainSequenceInfo[];

      if (source === "ebi") {
        // PDBe molecules API returns structured JSON with sequences per entity
        const url = `https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/${id}`;
        const result = await apiFetch(url);
        const data = await result.json();
        chains = parsePdbeMolecules(data, id);
      } else {
        // RCSB provides a proper FASTA endpoint
        const fastaUrl = `https://www.rcsb.org/fasta/entry/${id.toUpperCase()}`;
        const fastaText = await apiText(fastaUrl);
        chains = parseMultiChainFasta(fastaText);
      }

      if (chains.length === 0) {
        setMessage(`No polymer chains found in ${id.toUpperCase()}`);
        return;
      }

      if (chains.length === 1) {
        // Single chain — upload directly
        const chain = chains[0];
        const fastaContent = buildFasta(id, chain);
        const blob = new Blob([fastaContent], { type: "text/plain" });
        setMessage(`Uploading sequence for chain ${chain.chainId}`);
        uploadFile(blob, `${id.toUpperCase()}_${chain.chainId}.fasta`);
        onClose();
      } else {
        // Multiple chains — show picker
        setPendingChains(chains);
        setPendingPdbId(id);
        setChainPickerOpen(true);
      }
    } catch (err: any) {
      setMessage(err.message || `Failed to fetch sequences for ${id}`);
    }
  }, [identifier, uploadFile, onClose, setMessage, parsePdbeMolecules]);

  const handleFetch = useCallback(async () => {
    if (mode) {
      setInFlight(true);
      try {
        // For coord modes targeting a sequence parameter, fetch FASTA instead
        const targetIsSeq = isSequenceTarget(item, modes);
        if (targetIsSeq && mode === "ebiPdb") {
          await handlePdbSequenceFetch("ebi");
        } else if (targetIsSeq && mode === "rcsbPdb") {
          await handlePdbSequenceFetch("rcsb");
        } else if (mode === "ebiPdb") {
          await handleEbiCoordFetch();
        } else if (mode === "ebiSFs") {
          await handleEbiSFsFetch();
        } else if (mode === "uniprotFasta") {
          await handleUniprotFastaFetch();
        } else if (mode === "rcsbPdb") {
          await handleRcsbPdbFetch();
        } else if (mode === "uniprotAFPdb") {
          await handleUniprotAFPdbFetch();
        } else if (mode === "Uppsala-EDS") {
          await handleUppsalaEdsFetch();
        }
      } finally {
        setInFlight(false);
      }
    }
  }, [mode, item, modes, handleEbiCoordFetch, handleEbiSFsFetch, handleUniprotFastaFetch, handleRcsbPdbFetch, handleUniprotAFPdbFetch, handleUppsalaEdsFetch, handlePdbSequenceFetch]);

  return (
    <>
      <Dialog open={open} onClose={onClose}>
        <DialogTitle>
          {" "}
          Fetch value for {item?._objectPath?.split(".").at(-1)} from internet
        </DialogTitle>
        <DialogContent>
          {downloadModes?.length && downloadModes?.length > 0 ? (
            <>
              <Autocomplete
                sx={{ width: "30rem", mt: 2 }}
                value={mode || ""}
                renderInput={(params) => (
                  <TextField {...params} label="Download mode" />
                )}
                options={downloadModes}
                getOptionLabel={(option) => MODE_LABELS[option] || option}
                onChange={(event, newValue) => {
                  setMode(newValue);
                }}
              />
              <TextField
                sx={{ width: "30rem", mt: 2 }}
                label="Accession code"
                placeholder={mode ? MODE_PLACEHOLDERS[mode] || "" : ""}
                value={identifier || ""}
                onChange={(event) => {
                  setIdentifier(event.target.value);
                }}
                onKeyDown={(event) => {
                  if (event.key === "Enter" && identifier && !inFlight) {
                    handleFetch();
                  }
                }}
              />
            </>
          ) : (
            <Typography>No download modes</Typography>
          )}
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose} disabled={inFlight}>
            Cancel
          </Button>
          <Button onClick={handleFetch} disabled={inFlight}>
            Fetch
          </Button>
        </DialogActions>
      </Dialog>

      <ChainPickerDialog
        open={chainPickerOpen}
        onClose={() => setChainPickerOpen(false)}
        onSelect={uploadChainAsFasta}
        chains={pendingChains}
        pdbId={pendingPdbId}
      />
    </>
  );
};
