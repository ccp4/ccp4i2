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
import { usePopcorn } from "../../../providers/popcorn-provider";
import { useTaskInterface } from "../../../providers/task-provider";
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

  const { jobId, cootModule } = useCCP4i2Window();
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

  const uploadFile = useCallback(
    async (fileBlob: Blob, fileName: string) => {
      if (job && item) {
        setMessage(`Uploading file ${fileName} for ${item._objectPath}`);

        // Use centralized uploadFileParam with intent tracking
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
      }
    },
    [item, job, uploadFileParam, onChange, mutateJobs, mutateFiles, setMessage]
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
          console.log("FetchFileForParam handleFetch result", result);
        }
      } catch (err: any) {
        setMessage(err.message || "Unknown error");
        console.log("FetchFileForParam handleFetch error", err);
        return;
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleEbiSFsFetch = useCallback(async () => {
    if (identifier) {
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
        console.log("FetchFileForParam handleFetch result", result);
      }
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleUniprotFastaFetch = useCallback(async () => {
    if (identifier) {
      setMessage(`Fetching FASTA file for ${identifier.toUpperCase()}`);
      const data = await apiText(`uniprot/${identifier.toUpperCase()}.fasta`);
      setMessage(`Fetched FASTA file for ${identifier.toUpperCase()}`);
      const content = new Blob([data], {
        type: "text/plain",
      });
      uploadFile(content, `${identifier.toUpperCase()}.fasta`);
      onClose();
    }
  }, [identifier, uploadFile, onClose, setMessage]);

  const handleFetch = useCallback(async () => {
    if (mode) {
      setInFlight(true);
      if (mode === "ebiPdb") {
        await handleEbiCoordFetch();
      } else if (mode === "ebiSFs") {
        await handleEbiSFsFetch();
      } else if (mode === "uniprotFasta") {
        await handleUniprotFastaFetch();
      }
      setInFlight(false);
    }
  }, [mode, identifier]);

  return (
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
              onChange={(event, newValue) => {
                setMode(newValue);
              }}
            />
            <TextField
              sx={{ width: "30rem", mt: 2 }}
              label="Accession code"
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
  );
};
