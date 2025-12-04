"use client";
import React, { PropsWithChildren, useContext, useState, useRef } from "react";
import {
  Button,
  Container,
  Paper,
  Stack,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Typography,
} from "@mui/material";
import { useApi } from "../api";
import { Cancel, Check, Folder, Upload } from "@mui/icons-material";
import { VisuallyHiddenInput } from "./task/task-elements/input-file-upload";
import { useRouter } from "next/navigation";
import { Project } from "../types/models";
import { useTheme } from "../theme/theme-provider";

export const ImportProjectContent: React.FC = () => {
  const { customColors } = useTheme();
  const api = useApi();
  const router = useRouter();
  const [files, setFiles] = useState<File[]>([]);
  const [uploading, setUploading] = useState(false);
  const { mutate: mutateProjects } = api.get<Project[]>("projects");

  // Create a ref for the hidden file input
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileUpload = (selectedFiles: FileList | null) => {
    if (selectedFiles) {
      const filesArray = Array.from(selectedFiles);
      setFiles(filesArray);

      const formData = new FormData();
      for (let i = 0; i < selectedFiles.length; i++) {
        formData.append("files", selectedFiles[i]);
      }

      setUploading(true);
      api
        .post<any>("/projects/import_project/", formData)
        .then((response) => {
          console.log("Files uploaded successfully:", response.data);
          setUploading(false);
          mutateProjects();
          router.push("/");
        })
        .catch((error) => {
          console.error("Error uploading files:", error);
          setUploading(false);
        });
    }
  };

  const onChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    handleFileUpload(event.target.files);
  };

  const handlePaperClick = () => {
    fileInputRef.current?.click();
  };

  return (
    <Container
      sx={{
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
        height: "calc(100vh - 5rem)", // Full viewport height
        margin: 0, // Remove default margin
      }}
    >
      <Paper sx={{ padding: 2, minWidth: "50rem" }}>
        <Stack spacing={2}>
          <Typography variant="h4" gutterBottom>
            Import Project(s)
          </Typography>
          <Stack spacing={2} direction="row">
            <Paper
              sx={{
                border: "2px dashed #ccc",
                padding: 4,
                textAlign: "center",
                cursor: "pointer",
                backgroundColor: customColors.ui.lightGray,
                "&:hover": {
                  backgroundColor: "#f1f1f1",
                },
              }}
              onClick={handlePaperClick}
              onDragOver={(e) => e.preventDefault()}
              onDrop={(e) => {
                e.preventDefault();
                const droppedFiles = e.dataTransfer.files;
                setFiles(droppedFiles ? Array.from(droppedFiles) : []);
                handleFileUpload(droppedFiles);
              }}
            >
              <Typography variant="body1" color="textSecondary">
                Drag and drop files here, or click here to upload
              </Typography>
            </Paper>
            <Button
              component="label"
              variant="contained"
              startIcon={<Upload />}
            >
              <VisuallyHiddenInput
                ref={fileInputRef}
                type="file"
                multiple
                onChange={onChange}
              />
            </Button>
          </Stack>
          {files?.length > 0 && (
            <Table>
              <TableHead>
                <TableRow>
                  <TableCell>File</TableCell>
                  <TableCell>Processing</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {files.map((aFile, index) => (
                  <TableRow key={index}>
                    <TableCell>
                      <Folder /> {aFile.name}
                    </TableCell>
                    <TableCell>
                      {!uploading ? (
                        <Check color="success" />
                      ) : (
                        <Cancel color="error" />
                      )}
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          )}
        </Stack>
      </Paper>
    </Container>
  );
};
