"use client";
import { SyntheticEvent, useCallback } from "react";
import { Menu, MenuItem } from "@mui/material";
import { Preview, Download } from "@mui/icons-material";
import { useCCP4i2Window } from "../app-context";
import { doDownload, useApi } from "../api";
import { Project } from "../types/models";
import { useFilePreviewContext } from "../providers/file-preview-context";
import { useFileSystemFileBrowser } from "../providers/file-system-file-browser-context";

interface FileSystemFileMenuProps {
  onClose?: () => void;
}

export const FileSystemFileMenu: React.FC<FileSystemFileMenuProps> = ({
  onClose,
}) => {
  const { anchorEl, menuNode, closeMenu } = useFileSystemFileBrowser();

  const { projectId } = useCCP4i2Window();

  const { get, post } = useApi();

  const { data: project } = get<Project>(`projects/${projectId}`);

  const isOpen = Boolean(anchorEl && menuNode);

  const { setContentSpecification } = useFilePreviewContext();

  const handleClose = () => {
    closeMenu();
    onClose?.();
  };

  const handleDownload = useCallback(
    (ev: SyntheticEvent) => {
      if (!menuNode || !project) return;
      ev.stopPropagation();
      const composite_path = `/api/proxy/ccp4i2/projects/${
        project.id
      }/project_file?path=${encodeURIComponent(
        menuNode?.path.slice(project.directory.length + 1) || ""
      )}`;

      doDownload(composite_path, menuNode.name);
      handleClose();
    },
    [project, menuNode]
  );

  const handlePreview = useCallback(
    async (ev: SyntheticEvent) => {
      if (!menuNode || !project) return;
      ev.stopPropagation();
      const composite_path = `/api/proxy/ccp4i2/projects/${
        project.id
      }/project_file?path=${encodeURIComponent(
        menuNode?.path.slice(project.directory.length + 1) || ""
      )}`;

      setContentSpecification({
        url: composite_path,
        title: menuNode.name || "Preview",
        language: menuNode.name.endsWith(".json")
          ? "json"
          : menuNode.name.endsWith(".mtz")
          ? "mtz"
          : menuNode.name.endsWith(".cif")
          ? "cif"
          : menuNode.name.endsWith(".csv")
          ? "csv"
          : menuNode.name.endsWith(".xml")
          ? "xml"
          : menuNode.name.endsWith(".aln")
          ? "clustalw"
          : "text",
      });
      handleClose();
    },
    [projectId, menuNode]
  );

  const handleServersideViewer = async (
    project_id: number,
    viewer: string,
    path: string
  ) => {
    const jsonBody = { viewer, path };
    const previewResult: any = await post<any>(
      `projects/${project_id}/preview_file/`,
      jsonBody
    );
  };

  const handleServersideViewerInCoot = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (project) {
        handleServersideViewer(
          project.id,
          "coot",
          menuNode?.path.slice(project.directory.length) || ""
        );
      }
      handleClose();
    },
    [project, menuNode]
  );

  const handleServersideViewerInViewHKL = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (project) {
        handleServersideViewer(
          project.id,
          "viewhkl",
          menuNode?.path.slice(project.directory.length) || ""
        );
      }
      handleClose();
    },
    [project, menuNode]
  );

  const handleServersideViewerInCcp4mg = useCallback(
    async (ev: SyntheticEvent) => {
      ev.stopPropagation();
      if (project) {
        handleServersideViewer(
          project.id,
          "ccp4mg",
          menuNode?.path.slice(project.directory.length) || ""
        );
      }
      handleClose();
    },
    [project, menuNode]
  );

  return (
    menuNode &&
    anchorEl && (
      <Menu
        anchorEl={anchorEl}
        open={isOpen}
        onClose={handleClose}
        anchorOrigin={{
          vertical: "top",
          horizontal: "left",
        }}
        transformOrigin={{
          vertical: "top",
          horizontal: "left",
        }}
        disablePortal={true}
        keepMounted={false}
        disableScrollLock={true}
        slotProps={{
          paper: {
            sx: {
              zIndex: 10000,
              minWidth: 150,
              maxWidth: "calc(100vw - 20px)",
              maxHeight: "calc(100vh - 20px)",
            },
          },
        }}
      >
        {menuNode?.type !== "directory" && (
          <MenuItem onClick={handleDownload}>
            <Download />
            Download
          </MenuItem>
        )}
        {menuNode?.name &&
          [
            "aln",
            "log",
            "xml",
            "json",
            "txt",
            "mmcif",
            "cif",
            "csv",
            "pdb",
            "dict",
            "mol",
            "mtz",
          ].includes(menuNode?.name?.split(".").at(-1) || "") && (
            <MenuItem onClick={handlePreview}>
              <Preview />
              Preview
            </MenuItem>
          )}
        {menuNode?.name &&
          ["pdb", "mmcif", "cif", "mtz"].includes(
            menuNode?.name?.split(".").at(-1) || ""
          ) && (
            <MenuItem onClick={handleServersideViewerInCoot}>
              <Preview />
              Coot
            </MenuItem>
          )}
        {menuNode?.name &&
          ["pdb", "mmcif", "cif", "mtz"].includes(
            menuNode?.name?.split(".").at(-1) || ""
          ) && (
            <MenuItem onClick={handleServersideViewerInCcp4mg}>
              <Preview />
              CCP4MG
            </MenuItem>
          )}
        {menuNode?.name &&
          ["mtz"].includes(menuNode?.name?.split(".").at(-1) || "") && (
            <MenuItem onClick={handleServersideViewerInViewHKL}>
              <Preview />
              ViewHKL
            </MenuItem>
          )}{" "}
      </Menu>
    )
  );
};
