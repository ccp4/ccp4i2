import {
  ListItemIcon,
  ListItemText,
  MenuItem,
  Typography,
} from "@mui/material";
import { SvgIconComponent } from "@mui/icons-material";

export function CCP4i2MenuItem(props: {
  icon: SvgIconComponent;
  text: string;
  onClick: () => void;
  secondary?: string;
}) {
  return (
    <MenuItem onClick={props.onClick}>
      <ListItemIcon>
        <props.icon fontSize="small" />
      </ListItemIcon>
      <ListItemText>{props.text}</ListItemText>
      {props.secondary && (
        <Typography variant="body2" color="textSecondary" sx={{ ml: 2 }}>
          {props.secondary}
        </Typography>
      )}
    </MenuItem>
  );
}
