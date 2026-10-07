import { DialogTitle, DialogTitleProps, IconButton } from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";

interface DialogTitleWithCloseProps extends DialogTitleProps {
  /** Called when the X is clicked; pass the dialog's own onClose. */
  onClose: () => void;
}

/**
 * A dialog title with a close button (an X) in the top-right corner.
 *
 * Clicking the backdrop or pressing Esc closes a dialog, but neither is
 * visible, and users asked how to get out (#675). Use this in place of
 * `DialogTitle` in any dialog that views something rather than asks a
 * question; a confirmation keeps its explicit Cancel/OK instead.
 *
 * The button is positioned against the dialog's Paper (which MUI makes
 * `position: relative`), so it stays in the corner whatever the title holds.
 */
export const DialogTitleWithClose: React.FC<DialogTitleWithCloseProps> = ({
  onClose,
  children,
  sx,
  ...props
}) => (
  <DialogTitle sx={[{ pr: 6 }, ...(Array.isArray(sx) ? sx : sx ? [sx] : [])]} {...props}>
    {children}
    <IconButton
      aria-label="Close"
      onClick={onClose}
      sx={{ position: "absolute", right: 8, top: 8 }}
    >
      <CloseIcon />
    </IconButton>
  </DialogTitle>
);
