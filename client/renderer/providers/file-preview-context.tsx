import {
  PropsWithChildren,
  useCallback,
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
} from "@mui/material";
import { apiArrayBuffer, apiJson } from "../api-fetch";
import { Editor } from "@monaco-editor/react";
import { prettifyXml } from "../utils";
import { createContext } from "react";
import $ from "jquery";
import { CifTableStack } from "../components/cif-table-stack";
import { parseMtzHeader } from "../lib/mtz-parser";
import { CsvTable } from "../components/csv-table";
import { AlignmentViewer } from "../components/alignment-viewer";
import { MolBlockView } from "../components/campaigns/molblock-view";
import { useTheme } from "../theme/theme-provider";

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
  const { mode } = useTheme();

  // Parse MTZ header using pure TypeScript parser (no WASM dependency)
  const handleMtzPreview = useCallback(
    async (fileContent: ArrayBuffer) => {
      try {
        const header = parseMtzHeader(fileContent);
        // Format header info for preview display
        const previewData = {
          title: header.title,
          nColumns: header.nColumns,
          nReflections: header.nReflections,
          spaceGroup: header.spaceGroup,
          cell: header.cell,
          resolution: header.resolution,
          isMerged: header.isMerged,
          columns: header.columns.map(c => ({
            label: c.label,
            type: c.type,
          })),
          datasets: header.datasets,
        };
        setPreviewContent(JSON.stringify(previewData, null, 2));
      } catch (error) {
        console.error("Failed to parse MTZ file:", error);
        setPreviewContent(JSON.stringify({ error: String(error) }, null, 2));
      }
    },
    []
  );

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
            handleMtzPreview(fileContent);
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
  }, [contentSpecification, handleMtzPreview]);

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
        return "cif";
      case "csv":
        return "csv";
      case "mtz":
        return "json";
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
        ) : contentSpecification?.language === "molblock" ||
          contentSpecification?.language === "molblock-raw" ? (
          <Box sx={{ display: "flex", justifyContent: "center", alignItems: "center", minHeight: 200 }}>
            {previewContent ? (
              <MolBlockView molblock={previewContent} width={500} height={400} />
            ) : (
              <CircularProgress />
            )}
          </Box>
        ) : contentSpecification?.language === "cif" ? (
          <CifTableStack cifText={previewContent || ""} />
        ) : contentSpecification?.language === "csv" ? (
          <CsvTable csvText={previewContent || ""} />
        ) : contentSpecification?.language === "clustalw" ? (
          <AlignmentViewer alignment={previewContent || ""} />
        ) : (
          <Editor
            width="100%"
            height="calc(100vh - 20rem)"
            value={previewContent || ""}
            language={monacoLanguage}
            theme={mode === "dark" ? "vs-dark" : "light"}
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

export const useFilePreviewContext = () => {
  const context = useContext(FilePreviewContext);
  if (!context) {
    throw new Error(
      "useFilePreviewContext must be used within a FilePreviewProvider"
    );
  }
  return context;
};
