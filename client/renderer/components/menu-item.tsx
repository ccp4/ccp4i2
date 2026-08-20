import Typography from "@mui/material/Typography/Typography";

import MenuItem from "@mui/material/MenuItem/MenuItem";
import ListItemIcon from "@mui/material/ListItemIcon/ListItemIcon";
import ListItemText from "@mui/material/ListItemText/ListItemText";
import { SvgIconComponent } from "@mui/icons-material";

export function CCP4i2MenuItem(props: {
  icon: SvgIconComponent;
  text: string;
  onClick: () => void;
  shortcut?: string;
}) {
  return (
    <MenuItem onClick={props.onClick}>
      <ListItemIcon>
        <props.icon fontSize="small" />
      </ListItemIcon>
      <ListItemText>{props.text}</ListItemText>
      {props.shortcut && (
        <Typography variant="body2" color="textSecondary" sx={{ ml: 2 }}>
          {props.shortcut}
        </Typography>
      )}
    </MenuItem>
  );
}
