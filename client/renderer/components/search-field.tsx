import {
  IconButton,
  InputAdornment,
  SxProps,
  TextField,
  Theme,
} from "@mui/material";
import { Close, Search } from "@mui/icons-material";

export default function SearchField(props: {
  value: string;
  onChange: (value: string) => void;
  placeholder?: string;
  size?: "small" | "medium";
  sx?: SxProps<Theme>;
}) {
  return (
    <TextField
      value={props.value}
      onChange={(e) => props.onChange(e.target.value)}
      placeholder={props.placeholder || "Search..."}
      slotProps={{
        input: {
          startAdornment: (
            <InputAdornment position="start">
              <Search />
            </InputAdornment>
          ),
          endAdornment: (
            <InputAdornment position="end">
              <IconButton
                onClick={() => props.onChange("")}
                sx={{ visibility: props.value ? "visible" : "hidden" }}
              >
                <Close />
              </IconButton>
            </InputAdornment>
          ),
        },
      }}
      fullWidth
      size={props.size}
      sx={props.sx}
    />
  );
}
