import { Button, Chip, Stack } from "@mui/material";
import { Add } from "@mui/icons-material";
import { useApi } from "../api";
import { ProjectTag } from "../types/models";

export default function EditTags(props: {
  tags: number[];
  onChange: (tags: number[]) => void;
}) {
  const api = useApi();
  const { data: tags } = api.get<ProjectTag[]>("project-tags");

  function handleDelete() {}

  return (
    <Stack spacing={2} alignItems={"flex-start"}>
      {tags && tags.length > 0 && (
        <Stack direction="row" spacing={2} useFlexGap sx={{ flexWrap: "wrap" }}>
          {tags.map((tag) => (
            <Chip key={tag.id} label={tag.text} onDelete={handleDelete} />
          ))}
        </Stack>
      )}
      <Button variant="outlined" startIcon={<Add />}>
        Add Tag
      </Button>
    </Stack>
  );
}
