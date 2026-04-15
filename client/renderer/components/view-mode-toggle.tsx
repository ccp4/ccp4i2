import { ViewList, ViewModule } from "@mui/icons-material";
import { ToggleButton, ToggleButtonGroup, Tooltip } from "@mui/material";

export type ViewMode = "cards" | "list";

export function ViewModeToggle(props: {
  mode: ViewMode;
  onChange: (mode: ViewMode) => void;
}) {
  return (
    <ToggleButtonGroup
      exclusive
      value={props.mode}
      onChange={(_, mode) => mode && props.onChange(mode)}
      size="small"
    >
      <ToggleButton value="cards" aria-label="card view">
        <Tooltip title="Card view">
          <ViewModule />
        </Tooltip>
      </ToggleButton>
      <ToggleButton value="list" aria-label="list view">
        <Tooltip title="List view">
          <ViewList />
        </Tooltip>
      </ToggleButton>
    </ToggleButtonGroup>
  );
}
