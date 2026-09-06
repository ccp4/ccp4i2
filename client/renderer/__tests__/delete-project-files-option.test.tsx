/**
 * The delete-project dialog's "also delete the files on disk" choice: off by
 * default (the Qt-era default of leaving files), reported to the caller when
 * it changes, and explained in the wording.
 */
import { describe, expect, it, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { DeleteProjectFilesOption } from "@/components/delete-project-files-option";

describe("DeleteProjectFilesOption", () => {
  it("starts unticked and says the folder is kept", () => {
    const onChange = vi.fn();
    render(<DeleteProjectFilesOption onChange={onChange} />);
    const box = screen.getByRole("checkbox") as HTMLInputElement;
    expect(box.checked).toBe(false);
    expect(screen.getByText(/left where it is/)).toBeTruthy();
    expect(onChange).not.toHaveBeenCalled();
  });

  it("reports the choice and changes the wording when ticked", () => {
    const onChange = vi.fn();
    render(<DeleteProjectFilesOption onChange={onChange} />);
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onChange).toHaveBeenLastCalledWith(true);
    expect(screen.getByText(/everything in it will be removed/)).toBeTruthy();
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onChange).toHaveBeenLastCalledWith(false);
  });

  it("speaks of several projects when asked about several", () => {
    render(<DeleteProjectFilesOption onChange={() => {}} count={3} />);
    expect(screen.getByText(/projects' files on disk/)).toBeTruthy();
  });
});
