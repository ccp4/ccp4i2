import { FormControlLabel, Switch } from "@mui/material";
import { useCCP4i2Window } from "../app-context";

export const DevModeToggle: React.FC = () => {
  const { devMode, setDevMode } = useCCP4i2Window();

  const onToggleDevMode = async (
    ev: React.ChangeEvent<HTMLInputElement>
  ): Promise<void> => {
    setDevMode(!devMode);
    if (!window.electronAPI) {
      console.error("Electron API is not available");
      return;
    }
    window.electronAPI.sendMessage("toggle-dev-mode", {});
    ev.preventDefault();
    ev.stopPropagation();
  };

  return (
    <FormControlLabel
      control={
        <Switch
          checked={devMode}
          onChange={onToggleDevMode}
          name="devModeToggle"
          color="warning"
        />
      }
      label="Dev Mode"
    />
  );
};
