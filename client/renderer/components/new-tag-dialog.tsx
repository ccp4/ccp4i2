import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Stack,
  TextField,
} from "@mui/material";

export function NewTagDialog(props: { open: boolean }) {
  return (
    <Dialog open={props.open}>
      <DialogTitle>New project tag</DialogTitle>
      <DialogContent>
        <Stack gap={1}>
          <TextField label="New tag name" required autoFocus />
          <TextField label="Parent tag" />
        </Stack>
      </DialogContent>
      <DialogActions>
        <Button>Cancel</Button>
        <Button variant="contained">Create</Button>
      </DialogActions>
    </Dialog>
  );
}
