/**
 * The bulk-tagging route: apply one tag to a whole selection of projects.
 *
 * The thing worth locking down is that it goes through the tag's bulk
 * `add_projects` action in a single request, rather than one request per
 * project — and that creating a new tag on the way through still ends up
 * applying it.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";

const SEP = "\u001f";

const post = vi.fn((endpoint: string) =>
  endpoint === "projecttags/"
    ? Promise.resolve({ id: 42, text: "fresh" })
    : Promise.resolve({ status: "success" })
);
const mutate = vi.fn(() => Promise.resolve());

const TAGS = [
  {
    id: 1,
    text: "Nutlin site",
    parent: null,
    path: "Nutlin site",
    display_path: "Nutlin site",
    projects: [],
  },
  {
    id: 2,
    text: "soaks",
    parent: 1,
    path: `Nutlin site${SEP}soaks`,
    display_path: "Nutlin site/soaks",
    projects: [],
  },
];

vi.mock("../api", () => ({
  useApi: () => ({ get: () => ({ data: TAGS, mutate }), post }),
}));

import TagSelectionDialog from "../components/tag-selection-dialog";

describe("TagSelectionDialog", () => {
  beforeEach(() => post.mockClear());

  it("shows nested tags by their full path, so two leaves sharing a name are distinguishable", async () => {
    render(
      <TagSelectionDialog
        open
        projectIds={[7, 8]}
        onClose={vi.fn()}
        onApplied={vi.fn()}
      />
    );
    fireEvent.click(screen.getByLabelText("Open"));
    expect(await screen.findByText("Nutlin site/soaks")).toBeInTheDocument();
  });

  it("applies an existing tag to every selected project in one request", async () => {
    const onApplied = vi.fn();
    render(
      <TagSelectionDialog
        open
        projectIds={[7, 8]}
        onClose={vi.fn()}
        onApplied={onApplied}
      />
    );

    fireEvent.click(screen.getByLabelText("Open"));
    fireEvent.click(await screen.findByText("Nutlin site/soaks"));
    fireEvent.click(screen.getByText("Apply"));

    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/2/add_projects/", {
        project_ids: [7, 8],
      })
    );
    await waitFor(() =>
      expect(onApplied).toHaveBeenCalledWith("Nutlin site/soaks", 2)
    );
  });

  it("creates a tag typed by the user, then applies it", async () => {
    render(
      <TagSelectionDialog
        open
        projectIds={[7]}
        onClose={vi.fn()}
        onApplied={vi.fn()}
      />
    );

    fireEvent.click(screen.getByLabelText("Open"));
    // The label matches the listbox too once it is open, so target the input
    // by its role.
    fireEvent.change(screen.getByRole("combobox"), {
      target: { value: "fresh" },
    });
    fireEvent.click(await screen.findByText('Create tag "fresh"'));
    fireEvent.click(screen.getByText("Apply"));

    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/", {
        text: "fresh",
        parent: null,
        projects: [],
      })
    );
    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/42/add_projects/", {
        project_ids: [7],
      })
    );
  });
});
