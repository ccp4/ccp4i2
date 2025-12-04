import { createContext, PropsWithChildren, useContext, useState } from "react";
import { CircularProgress, Paper, Popper } from "@mui/material";
import { FetchFileForParam } from "../components/task/task-elements/fetch-file-for-param";
import { ErrorPopper } from "../components/task/task-elements/error-info";

interface TaskInterfaceContextProps {
  errorInfoAnchor: Element | null;
  setErrorInfoAnchor: (el: Element | null) => void;
  errorInfoItem: any | null;
  setErrorInfoItem: (item: any) => void;
  inFlight: any | null;
  setInFlight: (item: any) => void;
  fetchItemParams?: any | null;
  setFetchItemParams: (item: any) => void;
  downloadDialogOpen?: boolean;
  setDownloadDialogOpen?: (valie: boolean) => void;
}
export const TaskInterfaceContext = createContext<TaskInterfaceContextProps>({
  errorInfoAnchor: null,
  setErrorInfoAnchor: (item: any | null) => {},
  errorInfoItem: null,
  setErrorInfoItem: (item: any) => {},
  inFlight: null,
  setInFlight: (item: any) => {},
  fetchItemParams: null,
  setFetchItemParams: (item: any) => {},
  downloadDialogOpen: false,
  setDownloadDialogOpen: (value: boolean) => {},
});

export const TaskProvider: React.FC<PropsWithChildren> = ({ children }) => {
  const [errorInfoAnchor, setErrorInfoAnchor] = useState<Element | null>(null);
  const [errorInfoItem, setErrorInfoItem] = useState<any | null>(null);
  const [inFlight, setInFlight] = useState<boolean>(false);
  const [fetchItemParams, setFetchItemParams] = useState<any | null>(null);
  const [downloadDialogOpen, setDownloadDialogOpen] = useState<boolean>(false);

  return (
    <TaskInterfaceContext.Provider
      value={{
        errorInfoAnchor,
        setErrorInfoAnchor,
        errorInfoItem,
        setErrorInfoItem,
        inFlight,
        setInFlight,
        downloadDialogOpen,
        setDownloadDialogOpen,
        fetchItemParams,
        setFetchItemParams,
      }}
    >
      <Paper
        key="interface"
        sx={{ maxHeight: "calc(100vh - 22rem)", overflowY: "auto" }}
      >
        {children}{" "}
      </Paper>
      <Popper
        open={inFlight}
        placement="bottom-end"
        sx={{
          zIndex: 1300, // Ensure it's above other content
        }}
        style={{
          display: "grid",
          placeItems: "center" /* Centers both horizontally and vertically */,
          height: "100vh" /* Full height of the viewport */,
          width: "100vw",
        }}
      >
        <CircularProgress size={150} thickness={4} />
      </Popper>
      <FetchFileForParam
        open={downloadDialogOpen}
        onClose={() => {
          setDownloadDialogOpen(false);
        }}
      />
      <ErrorPopper key="error-popper" />
    </TaskInterfaceContext.Provider>
  );
};

export const useTaskInterface = () => {
  const context = useContext(TaskInterfaceContext);
  if (!context) {
    throw new Error("useTaskInterface must be used within a TaskProvider");
  }
  return context;
};
