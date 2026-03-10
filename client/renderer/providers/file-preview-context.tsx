import {
  PropsWithChildren,
  useContext,
  useEffect,
  useMemo,
  useState,
} from "react";
import { doRetrieve, doDownload } from "../api";
import {
  Box,
  CircularProgress,
  Dialog,
  DialogContent,
  DialogTitle,
  DialogActions,
  Button,
  FormControl,
  InputLabel,
  MenuItem,
  Select,
  Typography,
} from "@mui/material";
import { apiArrayBuffer, apiJson } from "../api-fetch";
import { Editor, loader } from "@monaco-editor/react";
import { prettifyXml } from "../utils";
import { createContext } from "react";
import $ from "jquery";
import { MtzPreview } from "../components/mtz-preview";
import { CsvTable } from "../components/csv-table";
import { AlignmentViewer } from "../components/alignment-viewer";
import { MolBlockView } from "../components/campaigns/molblock-view";
import { useTheme } from "../theme/theme-provider";

// Register mmCIF language and themes with Monaco at module load time.
// This runs once before any Editor component mounts.
loader.init().then((monaco) => {
  if (monaco.languages.getLanguages().some((l: { id: string }) => l.id === "mmcif")) return;

  monaco.languages.register({ id: "mmcif" });
  monaco.languages.setMonarchTokensProvider("mmcif", {
    tokenizer: {
      root: [
        [/#.*$/, "comment"],
        [/^data_\S+/, "keyword.data"],
        [/^loop_/, "keyword.loop"],
        [/^save_\S*/, "keyword.save"],
        [/_[\w.\[\]]+/, "tag"],
        [/'[^']*'/, "string"],
        [/"[^"]*"/, "string"],
        [/^;/, { token: "string.multiline", next: "@multilineString" }],
        [/[+-]?\d+\.\d*([eE][+-]?\d+)?/, "number.float"],
        [/[+-]?\d+([eE][+-]?\d+)?/, "number"],
        [/[?.](?=\s|$)/, "keyword.missing"],
      ],
      multilineString: [
        [/^;/, { token: "string.multiline", next: "@pop" }],
        [/.*/, "string.multiline"],
      ],
    },
  });

  monaco.editor.defineTheme("mmcif-light", {
    base: "vs",
    inherit: true,
    rules: [
      { token: "keyword.data", foreground: "8B0000", fontStyle: "bold" },
      { token: "keyword.loop", foreground: "8B0000", fontStyle: "bold" },
      { token: "keyword.save", foreground: "8B0000", fontStyle: "bold" },
      { token: "keyword.missing", foreground: "999999" },
      { token: "tag", foreground: "0055AA", fontStyle: "bold" },
      { token: "string", foreground: "A31515" },
      { token: "string.multiline", foreground: "A31515" },
      { token: "number", foreground: "098658" },
      { token: "number.float", foreground: "098658" },
      { token: "comment", foreground: "808080", fontStyle: "italic" },
    ],
    colors: {},
  });

  monaco.editor.defineTheme("mmcif-dark", {
    base: "vs-dark",
    inherit: true,
    rules: [
      { token: "keyword.data", foreground: "FF6B6B", fontStyle: "bold" },
      { token: "keyword.loop", foreground: "FF6B6B", fontStyle: "bold" },
      { token: "keyword.save", foreground: "FF6B6B", fontStyle: "bold" },
      { token: "keyword.missing", foreground: "808080" },
      { token: "tag", foreground: "6BB5FF", fontStyle: "bold" },
      { token: "string", foreground: "CE9178" },
      { token: "string.multiline", foreground: "CE9178" },
      { token: "number", foreground: "B5CEA8" },
      { token: "number.float", foreground: "B5CEA8" },
      { token: "comment", foreground: "6A9955", fontStyle: "italic" },
    ],
    colors: {},
  });
});

interface DictDigest {
  ligands: Array<{ id: string; name: string | null; group: string | null }>;
  monomers: Record<string, { atoms: string[]; bonds: any[] }>;
  molblocks: Record<string, string>;
}

export interface EditorContentSpecification {
  url: string;
  title: string;
  language: string;
}

interface FilePreviewDialogProps {
  contentSpecification: EditorContentSpecification | null;
  setContentSpecification: (spec: EditorContentSpecification | null) => void;
}

export const FilePreviewContext = createContext<FilePreviewDialogProps>({
  contentSpecification: null,
  setContentSpecification: () => {},
});

export const FilePreviewProvider: React.FC<PropsWithChildren> = ({
  children,
}) => {
  const [contentSpecification, setContentSpecification] =
    useState<EditorContentSpecification | null>(null);
  const contextValue = useMemo(
    () => ({ contentSpecification, setContentSpecification }),
    [contentSpecification]
  );

  return (
    <>
      <FilePreviewContext.Provider value={contextValue}>
        {children}
        <FilePreviewDialog />
      </FilePreviewContext.Provider>
    </>
  );
};

const FilePreviewDialog: React.FC = () => {
  const { contentSpecification, setContentSpecification } =
    useFilePreviewContext();
  const [previewContent, setPreviewContent] = useState<string | null>("");
  const [mtzData, setMtzData] = useState<ArrayBuffer | null>(null);
  const [dictDigest, setDictDigest] = useState<DictDigest | null>(null);
  const { mode } = useTheme();

  const compactCifPreview = (parsed: Record<string, any>): string => {
    const compact: Record<string, any> = {};
    Object.entries(parsed).forEach(([key, value]) => {
      if (
        key.startsWith("loop_") &&
        value.values &&
        Array.isArray(value.values)
      ) {
        compact[key] = {
          keys: value.keys,
          values: value.values.slice(0, 5), // Show only first 5 rows
          ...(value.values.length > 5 && {
            more: `${value.values.length - 5} more...`,
          }),
        };
      } else {
        compact[key] = value;
      }
    });
    return JSON.stringify(compact, null, 2);
  };

  const handleCifPreview = async (fileContent: ArrayBuffer) => {
    const enc = new TextDecoder("utf-8");
    const textContent = enc.decode(fileContent);
    setPreviewContent(textContent);
  };

  const handleCsvPreview = async (fileContent: ArrayBuffer) => {
    const enc = new TextDecoder("utf-8");
    const textContent = enc.decode(fileContent);
    setPreviewContent(textContent);
  };

  useEffect(() => {
    if (contentSpecification) {
      const asyncFunc = async () => {
        if (!contentSpecification.url) {
          return;
        }

        // Image files: use the URL directly as img src, no fetch needed
        if (contentSpecification.language === "image") {
          setPreviewContent(contentSpecification.url);
          return;
        }

        // Raw molfile: download the file content directly as molblock
        if (contentSpecification.language === "molblock-raw") {
          try {
            const fileContent = await apiArrayBuffer(contentSpecification.url);
            const enc = new TextDecoder("utf-8");
            setPreviewContent(enc.decode(fileContent));
          } catch (error) {
            console.error("Failed to fetch molfile:", error);
            setPreviewContent(null);
          }
          return;
        }

        // Molblock: fetch from molblock endpoint and extract molblock string
        if (contentSpecification.language === "molblock") {
          try {
            const response = await apiJson<{
              success: boolean;
              data: { molblock: string; ligand_code: string };
            }>(contentSpecification.url);
            setPreviewContent(response.data.molblock);
          } catch (error) {
            console.error("Failed to fetch molblock:", error);
            setPreviewContent(null);
          }
          return;
        }

        // Dictionary preview: fetch digest with molblocks for all monomers
        if (contentSpecification.language === "dict-preview") {
          try {
            const response = await apiJson<{
              success: boolean;
              data: DictDigest;
            }>(contentSpecification.url);
            const data = response.data ?? (response as any);
            setDictDigest(data);
          } catch (error) {
            console.error("Failed to fetch dictionary digest:", error);
            setDictDigest(null);
          }
          return;
        }

        {
          const fileContent = await apiArrayBuffer(contentSpecification.url);
          var enc = new TextDecoder("utf-8");
          if (contentSpecification.language === "json") {
            const fileText = enc.decode(fileContent);
            setPreviewContent(JSON.stringify(JSON.parse(fileText), null, 2));
          } else if (contentSpecification.language === "clustalw") {
            const fileText = enc.decode(fileContent);
            setPreviewContent(fileText);
          } else if (contentSpecification.language === "mtz") {
            setMtzData(fileContent);
          } else if (contentSpecification.language === "cif") {
            handleCifPreview(fileContent);
          } else if (contentSpecification.language === "csv") {
            handleCsvPreview(fileContent);
          } else if (contentSpecification.language === "xml") {
            const fileText = enc.decode(fileContent);
            setPreviewContent(prettifyXml($.parseXML(fileText)));
          } else {
            const fileText = enc.decode(fileContent);
            setPreviewContent(fileText);
          }
        }
      };
      asyncFunc();
    }
  }, [contentSpecification]);

  const handleDownload = () => {
    if (!contentSpecification?.url) return;
    doDownload(contentSpecification.url, contentSpecification.title || "download");
  };

  const monacoLanguage = useMemo(() => {
    switch (contentSpecification?.language) {
      case "clustal":
        return "clustal";
      case "json":
        return "json";
      case "xml":
        return "xml";
      case "cif":
        return "mmcif";
      case "csv":
        return "csv";
      default:
        return "text";
    }
  }, [contentSpecification]);

  return (
    <Dialog
      fullWidth
      maxWidth="xl"
      open={Boolean(contentSpecification)}
      onClose={() => {
        setContentSpecification(null);
        setMtzData(null);
        setDictDigest(null);
      }}
    >
      <DialogTitle>{contentSpecification?.title}</DialogTitle>
      <DialogContent>
        {contentSpecification?.language === "image" ? (
          <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: 200 }}>
            {previewContent ? (
              <img
                src={previewContent}
                alt={contentSpecification?.title || "Image preview"}
                style={{ maxWidth: "100%", maxHeight: "calc(100vh - 20rem)", objectFit: "contain" }}
              />
            ) : (
              <CircularProgress />
            )}
          </Box>
        ) : contentSpecification?.language === "dict-preview" ? (
          dictDigest ? (
            <DictPreview digest={dictDigest} />
          ) : (
            <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: 200 }}>
              <CircularProgress />
            </Box>
          )
        ) : contentSpecification?.language === "molblock" ||
          contentSpecification?.language === "molblock-raw" ? (
          <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: 200 }}>
            {previewContent ? (
              <MolBlockView molblock={previewContent} width={500} height={400} />
            ) : (
              <CircularProgress />
            )}
          </Box>
        ) : contentSpecification?.language === "csv" ? (
          <CsvTable csvText={previewContent || ""} />
        ) : contentSpecification?.language === "clustalw" ? (
          <AlignmentViewer alignment={previewContent || ""} />
        ) : contentSpecification?.language === "mtz" && mtzData ? (
          <MtzPreview data={mtzData} />
        ) : (
          <Editor
            width="100%"
            height="calc(100vh - 20rem)"
            value={previewContent || ""}
            language={monacoLanguage}
            theme={
              monacoLanguage === "mmcif"
                ? mode === "dark" ? "mmcif-dark" : "mmcif-light"
                : mode === "dark" ? "vs-dark" : "light"
            }
            options={{ readOnly: true, wordWrap: "on" }}
          />
        )}
      </DialogContent>
      <DialogActions>
        <Button
          onClick={handleDownload}
          disabled={!contentSpecification?.url}
          variant="contained"
          color="primary"
        >
          Download File
        </Button>
      </DialogActions>
    </Dialog>
  );
};

/** Inline component for dictionary preview with monomer selector */
const DictPreview: React.FC<{ digest: DictDigest }> = ({ digest }) => {
  const codes = Object.keys(digest.molblocks);
  const [selectedCode, setSelectedCode] = useState(codes[0] ?? "");

  // Find the ligand metadata for the selected monomer
  const selectedLigand = digest.ligands.find((l) => l.id === selectedCode);
  const selectedMolblock = digest.molblocks[selectedCode];

  if (codes.length === 0) {
    return (
      <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: 200 }}>
        <Typography color="text.secondary">No structures available</Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ display: "flex", flexDirection: "column", alignItems: "center", gap: 2 }}>
      {codes.length > 1 && (
        <FormControl size="small" sx={{ minWidth: 250 }}>
          <InputLabel>Monomer</InputLabel>
          <Select
            value={selectedCode}
            label="Monomer"
            onChange={(e) => setSelectedCode(e.target.value)}
          >
            {codes.map((code) => {
              const ligand = digest.ligands.find((l) => l.id === code);
              const label = ligand?.name ? `${code} \u2013 ${ligand.name}` : code;
              return (
                <MenuItem key={code} value={code}>
                  {label}
                </MenuItem>
              );
            })}
          </Select>
        </FormControl>
      )}
      {selectedMolblock ? (
        <MolBlockView
          molblock={selectedMolblock}
          name={selectedLigand?.name ?? selectedCode}
          width={500}
          height={400}
        />
      ) : (
        <Typography color="text.secondary">
          No structure available for {selectedCode}
        </Typography>
      )}
    </Box>
  );
};

export const useFilePreviewContext = () => {
  const context = useContext(FilePreviewContext);
  if (!context) {
    throw new Error(
      "useFilePreviewContext must be used within a FilePreviewProvider"
    );
  }
  return context;
};
