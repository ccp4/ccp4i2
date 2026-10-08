/**
 * The delete-job dialog's "also delete the imported files" choice (#587):
 * off by default (Qt-era CCP4i2 asked, and keeping is the default here),
 * offered only when the jobs imported something, and the list of jobs that
 * go, and whether Delete is allowed, follow the choice.
 */
import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import {
  confirmJobDeletion,
  DeletePreview,
  JobDeletionDetails,
} from "@/components/job-deletion-details";
import { deleteDialogReducer } from "@/providers/delete-dialog";
import { Job } from "@/types/models";

vi.mock("@/components/job-avatar", () => ({ CCP4i2JobAvatar: () => null }));

const job = (id: number, number: string, title: string, status = 6) =>
  ({ id, uuid: `u${id}`, number, title, parent: null, status }) as unknown as Job;

const brief = (id: number, number: string) => ({
  id,
  number,
  title: `job ${number}`,
  parent: null,
  status: 6,
});

function preview(overrides: Partial<DeletePreview> = {}): DeletePreview {
  return {
    selected_jobs: [job(1, "1", "import job")],
    imported_files: [
      {
        id: 10,
        uuid: "f10",
        name: "X_imported.pdb",
        source_name: "model.pdb",
        annotation: "",
        job: brief(1, "1"),
        used_by: [brief(2, "2")],
      },
    ],
    keep_imported_files: {
      additional_dependents: [job(3, "3", "uses the output")],
      total_to_delete: 2,
      has_active_dependents: false,
    },
    delete_imported_files: {
      additional_dependents: [
        job(2, "2", "uses the import", 3),
        job(3, "3", "uses the output"),
      ],
      total_to_delete: 3,
      has_active_dependents: true,
    },
    ...overrides,
  };
}

const notBlocked = () => false;

describe("JobDeletionDetails", () => {
  it("offers the choice unticked, naming the files and who used them", () => {
    const onChange = vi.fn();
    render(
      <JobDeletionDetails preview={preview()} onChange={onChange} isBlocked={notBlocked} />
    );
    const box = screen.getByRole("checkbox") as HTMLInputElement;
    expect(box.checked).toBe(false);
    expect(screen.getByText("Also delete the 1 imported file")).toBeTruthy();
    expect(screen.getByText("model.pdb")).toBeTruthy();
    expect(screen.getByText(/imported by job 1; used by job 2/)).toBeTruthy();
    expect(screen.getByText(/stay in the project/)).toBeTruthy();
    expect(screen.getByText(/stays in the job list as a file holder/)).toBeTruthy();
    expect(screen.getByText("3: uses the output")).toBeTruthy();
    expect(screen.queryByText("2: uses the import")).toBeNull();
    expect(onChange).not.toHaveBeenCalled();
  });

  it("lists the jobs that used the files once they are to be deleted", () => {
    const onChange = vi.fn();
    render(
      <JobDeletionDetails preview={preview()} onChange={onChange} isBlocked={notBlocked} />
    );
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onChange).toHaveBeenLastCalledWith(true);
    expect(screen.getByText("2: uses the import")).toBeTruthy();
    expect(screen.getByText(/so is the job that used them \(2\)/)).toBeTruthy();
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onChange).toHaveBeenLastCalledWith(false);
  });

  it("offers no choice when nothing was imported", () => {
    render(
      <JobDeletionDetails
        preview={preview({ imported_files: [] })}
        onChange={() => {}}
        isBlocked={notBlocked}
      />
    );
    expect(screen.queryByRole("checkbox")).toBeNull();
  });

  it("says why Delete is blocked, and that keeping the files unblocks it", () => {
    render(
      <JobDeletionDetails
        preview={preview()}
        onChange={() => {}}
        isBlocked={(plan) => plan.has_active_dependents}
      />
    );
    expect(screen.queryByText(/still active/)).toBeNull();
    fireEvent.click(screen.getByRole("checkbox"));
    expect(screen.getByText(/or keep the imported files/)).toBeTruthy();
  });
});

describe("confirmJobDeletion", () => {
  it("asks the server, shows the dialog, and reports the choice on Delete", async () => {
    const post = vi.fn().mockResolvedValue({ data: preview() });
    const deleteDialog = vi.fn();
    const onDelete = vi.fn();
    await confirmJobDeletion({
      api: { post } as any,
      deleteDialog,
      jobIds: [1],
      what: "1: import job",
      isBlocked: (plan) => plan.has_active_dependents,
      onDelete,
    });
    expect(post).toHaveBeenCalledWith("jobs/delete_preview/", { job_ids: [1] });
    const shown = deleteDialog.mock.calls[0][0];
    expect(shown.type).toBe("show");
    expect(shown.deleteDisabled).toBe(false);
    shown.onDelete();
    expect(onDelete).toHaveBeenCalledWith(false, preview().keep_imported_files);
  });
});

describe("deleteDialogReducer update", () => {
  it("changes whether Delete is enabled only while the dialog is open", () => {
    const open = deleteDialogReducer(
      { open: false },
      { type: "show", what: "x", deleteDisabled: false }
    );
    expect(
      deleteDialogReducer(open, { type: "update", deleteDisabled: true })
        .deleteDisabled
    ).toBe(true);
    const closed = { open: false };
    expect(deleteDialogReducer(closed, { type: "update", deleteDisabled: true })).toBe(
      closed
    );
  });
});
