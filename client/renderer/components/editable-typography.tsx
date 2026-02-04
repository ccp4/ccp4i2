import { Typography } from "@mui/material";
import { KeyboardEvent, useRef } from "react";
import type { TypographyProps, SxProps, Theme } from "@mui/material";

type TypographyVariant = TypographyProps["variant"];

export default function EditableTypography(props: {
  variant: TypographyVariant;
  text: string;
  onDelay?: (text: string) => void;
  sx?: SxProps<Theme>;
}) {
  const timeout = useRef<NodeJS.Timeout | undefined>(undefined);

  function handleKeyUp(e: KeyboardEvent<HTMLSpanElement>) {
    const text = e.currentTarget.innerText.trim();
    clearTimeout(timeout.current);
    timeout.current = setTimeout(() => props.onDelay?.(text), 1000);
  }

  return (
    <Typography
      variant={props.variant}
      contentEditable={Boolean(props.onDelay)}
      suppressContentEditableWarning
      onKeyUp={handleKeyUp}
      sx={props.sx}
    >
      {props.text}
    </Typography>
  );
}
