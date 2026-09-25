import {
  Avatar,
  LinearProgress,
  ListItem,
  ListItemAvatar,
  ListItemButton,
  ListItemText,
} from "@mui/material";
import { fileTypeMapping } from "../files-table";
import { useFileMenu } from "../../providers/file-context-menu";
import { useCCP4i2Window } from "../../app-context";
import { useProject } from "../../utils";
import { useApi } from "../../api";
import { File as CCP4i2File } from "../../types/models";

export function CCP4i2ReportFile(props: { uuid: string }) {
  const { projectId } = useCCP4i2Window();
  const { files } = useProject(projectId);
  const api = useApi();
  const ownFile = (files ?? []).find((f) => f.uuid === props.uuid);
  // A job's inputs are not always this project's files: a campaign-level
  // job (pandda_campaign) takes its members' outputs, and a project's file
  // list never contains those. Fetch by uuid when the local list has no
  // answer; the undefined-forever alternative rendered as a progress bar
  // per file, indefinitely.
  const { data: foreignFile } = api.get<CCP4i2File>(
    files && !ownFile ? `files_by_uuid/${props.uuid}/` : null
  );
  const file = ownFile ?? foreignFile;
  const { setFileMenuAnchorEl, setFile } = useFileMenu();

  const fileTypeIcon = (file?.type && fileTypeMapping[file.type]) || "ccp4";

  function handleMenuClick(ev: React.MouseEvent<HTMLElement>) {
    ev.stopPropagation();
    ev.preventDefault();
    setFileMenuAnchorEl(ev.currentTarget);
    if (file) setFile(file);
  }

  if (!file) return <LinearProgress />;

  return (
    <ListItem
      onClick={handleMenuClick}
      onContextMenu={handleMenuClick}
      disablePadding
    >
      <ListItemButton>
        <ListItemAvatar>
          <Avatar
            src={`/svgicons/${fileTypeIcon}.svg`}
            sx={{ width: "2rem", height: "2rem" }}
          />
        </ListItemAvatar>
        <ListItemText>
          {file.annotation || file.job_param_name || ""}
        </ListItemText>
      </ListItemButton>
    </ListItem>
  );
}
