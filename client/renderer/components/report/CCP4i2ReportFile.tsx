import { useCallback, useContext, useMemo } from "react";
import { CCP4i2ReportElementProps } from "./CCP4i2ReportElement";
import {
  Avatar,
  Button,
  Chip,
  LinearProgress,
  Stack,
  Typography,
} from "@mui/material";
import { Menu as MenuIcon } from "@mui/icons-material";
import { useApi } from "../../api";
import { fileTypeMapping } from "../files-table";
import EditableTypography from "../editable-typography";
import { File as DjangoFile } from "../../types/models";
import { useFileMenu } from "../../providers/file-context-menu";
import { useCCP4i2Window } from "../../app-context";
import { useProject } from "../../utils";
//import { fileTypeMapping } from "../files-table";

interface CCP4i2ReportFileProps extends CCP4i2ReportElementProps {
  uuid: string;
}
export const CCP4i2ReportFile: React.FC<CCP4i2ReportFileProps> = (props) => {
  const api = useApi();
  const { projectId } = useCCP4i2Window();
  const { files, mutateFiles } = useProject(projectId || 0);
  const file = useMemo(
    () => (files ?? []).find((f) => f.uuid === props.uuid),
    [files, props.uuid]
  );
  const { setFileMenuAnchorEl, setFile } = useFileMenu();

  const fileTypeIcon = useMemo(() => {
    if (!file?.type) return "ccp4";
    return Object.keys(fileTypeMapping).includes(file?.type)
      ? fileTypeMapping[file?.type]
      : "ccp4";
  }, [file]);

  const handleMenuClick = useCallback(
    (ev: React.MouseEvent<HTMLElement>) => {
      ev.stopPropagation();
      ev.preventDefault();
      setFileMenuAnchorEl(ev.currentTarget);
      if (file) setFile(file);
    },
    [file, setFileMenuAnchorEl, setFile]
  );

  // Handle right-click context menu
  const handleContextMenu = useCallback(
    (ev: React.MouseEvent<HTMLElement>) => {
      handleMenuClick(ev);
    },
    [handleMenuClick]
  );

  // Handle button click (left-click on menu button)
  const handleButtonClick = useCallback(
    (ev: React.MouseEvent<HTMLButtonElement>) => {
      handleMenuClick(ev);
    },
    [handleMenuClick]
  );

  if (!file) return <LinearProgress />;

  return (
    <>
      <Stack
        direction="row"
        onContextMenu={handleContextMenu}
        sx={{
          border: "3px solid",
          borderRadius: "0.5rem",
          mx: 2,
          my: 1,
          p: 1,
          cursor: "context-menu",
          "&:hover": {
            backgroundColor: "action.hover",
          },
        }}
      >
        <Avatar
          src={`/api/proxy/djangostatic/qticons/${fileTypeIcon}.png`}
          sx={{ mr: 2, width: "2rem", height: "2rem" }}
        />
        <EditableTypography
          variant="body1"
          text={
            file?.annotation
              ? file?.annotation
              : file?.job_param_name
              ? file.job_param_name
              : ""
          }
          onDelay={async (annotation) => {
            const formData = new FormData();
            formData.set("annotation", annotation);
            await api.patch(`files/${file?.id}`, formData);
            mutateFiles();
          }}
        />
        <Typography sx={{ flexGrow: 1 }} />
        <div>
          <Button
            size="small"
            sx={{ p: 0, m: 0 }}
            variant="outlined"
            onClick={handleButtonClick}
          >
            <MenuIcon fontSize="small" />
          </Button>
        </div>
      </Stack>
    </>
  );
};
