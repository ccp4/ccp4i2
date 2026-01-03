"use client";
import { ChangeEvent, use } from "react";
import { Container, LinearProgress, Stack, Toolbar } from "@mui/material";
import { useApi } from "../../../../api";
import { File, Project } from "../../../../types/models";
import EditableTypography from "../../../../components/editable-typography";
import FilesTable from "../../../../components/files-table";
import FileUpload from "../../../../components/file-upload";
import { useProject } from "../../../../utils";
import { usePopcorn } from "../../../../providers/popcorn-provider";

export default function FilesPage({
  params,
}: {
  params: Promise<{ id: string }>;
}) {
  const api = useApi();
  const { id } = use(params);
  const { data: project } = api.get<Project>(`projects/${id}`);
  const { files, mutateFiles } = useProject(id ? parseInt(id) : undefined);
  const { setMessage } = usePopcorn();

  async function importFiles(event: ChangeEvent<HTMLInputElement>) {
    const fileList = event.target.files;
    if (fileList && project) {
      for (let i = 0; i < fileList.length; i++) {
        const formData = new FormData();
        formData.append("project", project.id.toString());
        formData.append("file", fileList[i]);
        formData.append("name", fileList[i].name);
        formData.append("size", fileList[i].size.toString());
        try {
          const result: any = await api.post("files", formData);
          if (result?.success === false) {
            setMessage(`Failed to import ${fileList[i].name}: ${result?.error || "Unknown error"}`, "error");
          } else {
            mutateFiles();
            setMessage(`Imported ${fileList[i].name}`, "success");
          }
        } catch (error) {
          setMessage(`Error importing ${fileList[i].name}: ${error instanceof Error ? error.message : String(error)}`, "error");
        }
      }
    }
  }

  if (!project) return <LinearProgress />;
  return (
    <Stack spacing={2}>
      <Container>
        <EditableTypography
          variant="h4"
          text={project.name}
          onDelay={(name) =>
            api.patch(`projects/${project.id}`, { name: name })
          }
        />
        <Toolbar disableGutters>
          <FileUpload text="Import Files" onChange={importFiles} />
        </Toolbar>
        <FilesTable files={files} mutate={mutateFiles} />
      </Container>
    </Stack>
  );
}
