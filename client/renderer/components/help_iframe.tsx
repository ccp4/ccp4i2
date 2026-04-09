import { Delete } from "@mui/icons-material";
import { Dialog, DialogContent, DialogTitle, IconButton } from "@mui/material";

interface HelpIframeProps {
  url: string;
  open: boolean;
  handleClose: () => void;
}

export const HelpIframe: React.FC<HelpIframeProps> = ({
  url,
  open,
  handleClose,
}) => (
  <Dialog open={open} onClose={handleClose} maxWidth="lg" fullWidth>
    <DialogTitle>
      Embedded Page
      <IconButton
        onClick={handleClose}
        style={{ position: "absolute", right: 10, top: 10 }}
      >
        <Delete />
      </IconButton>
    </DialogTitle>

    <DialogContent style={{ padding: 0 }}>
      {/* Iframe for external page */}
      <iframe
        src={url}
        style={{
          width: "100%",
          height: "500px",
          border: "none",
        }}
        sandbox="allow-same-origin allow-scripts allow-popups allow-forms"
      />
    </DialogContent>
  </Dialog>
);
