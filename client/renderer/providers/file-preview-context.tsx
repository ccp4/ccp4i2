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
  Dialog,
  DialogContent,
  DialogTitle,
  DialogActions,
  Button,
} from "@mui/material";
import { apiArrayBuffer } from "../api-fetch";
import { Editor } from "@monaco-editor/react";
import { prettifyXml } from "../utils";
import { createContext } from "react";
import $ from "jquery";
import { useCCP4i2Window } from "../app-context";
import { CifTableStack } from "../components/cif-table-stack";
import { CsvTable } from "../components/csv-table";
import { AlignmentViewer } from "../components/alignment-viewer";
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
  const { cootModule } = useCCP4i2Window();
  const { mode } = useTheme();
  const handleMtzPreview = useCallback(
    async (fileContent: ArrayBuffer) => {
      if (cootModule) {
        const byteArray = new Uint8Array(fileContent);
        try {
          cootModule.FS_unlink("/tmp/fileName");
        } catch (e) {}
        cootModule.FS_createDataFile("/tmp", "fileName", byteArray, true, true);
        const header_info_em: any = cootModule.get_mtz_columns("/tmp/fileName");
        cootModule.FS_unlink("/tmp/fileName");
        // Convert emscripten array to regular array for easier handling
        const header_info: string[] = [];
        for (let key = 0; key < header_info_em.size(); key++) {
          header_info.push(header_info_em.get(key));
        }
        console.log(header_info);
        setPreviewContent(JSON.stringify(header_info, null, 2));
      }
    },
    [cootModule]
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
  }, [contentSpecification, cootModule]);

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
        {contentSpecification?.language === "cif" ? (
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
