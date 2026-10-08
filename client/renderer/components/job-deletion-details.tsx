"use client";
import { useEffect, useState } from "react";
import {
  Box,
  Checkbox,
  FormControlLabel,
  List,
  ListItem,
  Paper,
  Toolbar,
  Typography,
} from "@mui/material";
import { Job } from "../types/models";
import { CCP4i2JobAvatar } from "./job-avatar";
import { useDeleteDialog } from "../providers/delete-dialog";

/** A job as the delete preview names it alongside an imported file. */
export interface JobBrief {
  id: number;
  number: string;
  title: string;
  parent: number | null;
  status: number;
}

/** One file the jobs being deleted imported into the project. */
export interface ImportedFileAtStake {
  id: number;
  uuid: string;
  name: string;
  /** The file's name where it was imported from. */
  source_name: string;
  annotation: string;
  /** The job that imported it. */
  job: JobBrief;
  /** Jobs outside the selection that used it; they go if it goes. */
  used_by: JobBrief[];
}

/** What deleting removes, for one answer to "delete the imported files?". */
export interface DeletionPlan {
  additional_dependents: Job[];
  total_to_delete: number;
  has_active_dependents: boolean;
}

/** POST jobs/delete_preview/ */
export interface DeletePreview {
  selected_jobs: Job[];
  imported_files: ImportedFileAtStake[];
  keep_imported_files: DeletionPlan;
  delete_imported_files: DeletionPlan;
}

export function planFor(
  preview: DeletePreview,
  deleteImportedFiles: boolean
): DeletionPlan {
  return deleteImportedFiles
    ? preview.delete_imported_files
    : preview.keep_imported_files;
}

const plural = (n: number, one: string, many = `${one}s`) =>
  n === 1 ? one : many;

const jobList = (jobs: JobBrief[]) => jobs.map((j) => j.number).join(", ");

interface JobDeletionDetailsProps {
  preview: DeletePreview;
  /** Called with the current choice whenever it changes. Starts false. */
  onChange: (deleteImportedFiles: boolean) => void;
  /** Whether a plan may not go ahead (e.g. a dependent job is running). */
  isBlocked: (plan: DeletionPlan) => boolean;
}

/**
 * The body of the delete-job dialog: which other jobs go, and --- when the
 * jobs imported files --- whether those files go too.
 *
 * Qt-era CCP4i2 asked "delete imported files too?"; this asks the same and
 * defaults to keeping them (#587). Kept files stay in the project, the job
 * that imported them stays as a file holder, and jobs that used them survive.
 * Deleting them takes those jobs too, so the list of dependents follows the
 * choice, and so does whether Delete is allowed. Like the project dialog's
 * files option it keeps its own state, and the caller reads the choice back
 * through onChange when Delete is pressed.
 */
export function JobDeletionDetails({
  preview,
  onChange,
  isBlocked,
}: JobDeletionDetailsProps) {
  const [deleteImportedFiles, setDeleteImportedFiles] = useState(false);
  const deleteDialog = useDeleteDialog();
  const plan = planFor(preview, deleteImportedFiles);
  const blocked = isBlocked(plan);

  useEffect(() => {
    deleteDialog?.({ type: "update", deleteDisabled: blocked });
  }, [blocked, deleteDialog]);

  const dependents = plan.additional_dependents.filter(
    (job) => job.parent === null
  );
  const files = preview.imported_files;
  const importers = Array.from(
    new Map(files.map((f) => [f.job.id, f.job])).values()
  );
  const users = Array.from(
    new Map(
      files.flatMap((f) => f.used_by).map((j) => [j.id, j] as const)
    ).values()
  ).filter((j) => j.parent === null);

  return (
    <Box>
      {files.length > 0 && (
        <Box sx={{ mt: 2 }}>
          <FormControlLabel
            control={
              <Checkbox
                checked={deleteImportedFiles}
                onChange={(event) => {
                  setDeleteImportedFiles(event.target.checked);
                  onChange(event.target.checked);
                }}
              />
            }
            label={`Also delete the ${files.length} imported ${plural(
              files.length,
              "file"
            )}`}
          />
          <List dense disablePadding sx={{ maxHeight: "8rem", overflowY: "auto" }}>
            {files.map((file) => (
              <ListItem key={file.uuid} sx={{ py: 0 }}>
                <Typography variant="body2">
                  {file.source_name || file.name}
                  <Typography component="span" variant="body2" color="text.secondary">
                    {` (imported by job ${file.job.number}`}
                    {file.used_by.length > 0
                      ? `; used by ${plural(file.used_by.length, "job")} ${jobList(
                          file.used_by
                        )})`
                      : ")"}
                  </Typography>
                </Typography>
              </ListItem>
            ))}
          </List>
          <Typography variant="body2" color="text.secondary">
            {deleteImportedFiles
              ? users.length > 0
                ? `The files are removed from the project, and so ${plural(
                    users.length,
                    "is the job",
                    "are the jobs"
                  )} that used them (${jobList(users)}).`
                : "The files are removed from the project."
              : `The files stay in the project, for any job to use. ${plural(
                  importers.length,
                  "Job",
                  "Jobs"
                )} ${jobList(importers)} ${plural(
                  importers.length,
                  "stays",
                  "stay"
                )} in the job list as a file holder.`}
          </Typography>
        </Box>
      )}
      {dependents.length > 0 && (
        <Paper sx={{ mt: 2, maxHeight: "10rem", overflowY: "auto" }}>
          <Typography variant="body2" sx={{ m: 1 }}>
            The following {dependents.length} dependent{" "}
            {plural(dependents.length, "job")} would also be deleted:
          </Typography>
          <List dense>
            {dependents.map((job) => (
              <ListItem key={job.uuid}>
                <Toolbar>
                  <CCP4i2JobAvatar job={job} />
                  {`${job.number}: ${job.title}`}
                </Toolbar>
              </ListItem>
            ))}
          </List>
        </Paper>
      )}
      {blocked && (
        <Typography variant="body2" color="error" sx={{ mt: 1 }}>
          A job that would be deleted is still active; wait for it to finish
          {deleteImportedFiles && !isBlocked(preview.keep_imported_files)
            ? ", or keep the imported files."
            : "."}
        </Typography>
      )}
    </Box>
  );
}

type DeleteDialog = ReturnType<typeof useDeleteDialog>;

interface ConfirmJobDeletionArgs {
  api: { post: <T>(endpoint: string, body?: any) => Promise<T> };
  deleteDialog: DeleteDialog;
  jobIds: number[];
  /** Named in the title: "Delete <what>?" */
  what: string;
  isBlocked: (plan: DeletionPlan) => boolean;
  /** Called when Delete is pressed, with the choice and what it removes. */
  onDelete: (
    deleteImportedFiles: boolean,
    plan: DeletionPlan
  ) => void | Promise<void>;
  onCancel?: () => void;
}

/**
 * Ask the server what deleting ``jobIds`` would remove, then show the delete
 * dialog with the choice about imported files.
 */
export async function confirmJobDeletion({
  api,
  deleteDialog,
  jobIds,
  what,
  isBlocked,
  onDelete,
  onCancel,
}: ConfirmJobDeletionArgs) {
  if (!deleteDialog) return;
  const response = await api.post<{ data: DeletePreview }>(
    "jobs/delete_preview/",
    { job_ids: jobIds }
  );
  const preview = response.data;
  const choice = { deleteImportedFiles: false };
  deleteDialog({
    type: "show",
    what,
    children: [
      <JobDeletionDetails
        key="job-deletion-details"
        preview={preview}
        isBlocked={isBlocked}
        onChange={(value) => {
          choice.deleteImportedFiles = value;
        }}
      />,
    ],
    deleteDisabled: isBlocked(preview.keep_imported_files),
    onDelete: () =>
      onDelete(
        choice.deleteImportedFiles,
        planFor(preview, choice.deleteImportedFiles)
      ),
    onCancel,
  });
}
