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

export function CCP4i2ReportFile(props: { uuid: string }) {
  const { projectId } = useCCP4i2Window();
  const { files } = useProject(projectId);
  const file = (files ?? []).find((f) => f.uuid === props.uuid);
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
