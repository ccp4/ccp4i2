/**
 * Dropping a tag from the tree onto a project in the list applies that tag to
 * that one project — the reverse of dragging a selection onto a tag node.
 *
 * The list is virtualised, and jsdom has no layout, so the virtualiser is
 * replaced with one that renders every row; the drag-and-drop wiring on the
 * rows is what this locks down, not the windowing.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";

const TAG_TYPE = "application/x-ccp4i2-tag";

const PROJECTS = [
  {
    id: 7,
    uuid: "u7",
    name: "Nutlin soak 12",
    description: "",
    tags: [],
    creation_time: "2026-01-01T00:00:00Z",
    last_access: "2026-01-02T00:00:00Z",
  },
];

const TREE = {
  tags: [
    {
      id: 3,
      text: "Mpro",
      parent: null,
      path: "Mpro",
      display_path: "Mpro",
      depth: 0,
      project_count: 0,
      total_project_count: 0,
      children: [],
    },
  ],
  untagged_project_count: 1,
};

const post = vi.fn(() => Promise.resolve({ status: "success" }));
const mutateProjects = vi.fn(() => Promise.resolve());
const mutateTree = vi.fn(() => Promise.resolve());
const setMessage = vi.fn();

vi.mock("../api", () => ({
  useApi: () => ({
    get: (endpoint: string | null) => {
      if (endpoint === "projects") return { data: PROJECTS, mutate: mutateProjects };
      if (endpoint === "projecttags/tree/") return { data: TREE, mutate: mutateTree };
      return { data: undefined, mutate: vi.fn() };
    },
    post,
    patch: vi.fn(),
    delete: vi.fn(),
  }),
}));

vi.mock("swr", () => ({ default: () => ({ data: {} }) }));
vi.mock("../api-fetch", () => ({ apiFetch: vi.fn() }));
vi.mock("../providers/popcorn-provider", () => ({
  usePopcorn: () => ({ setMessage }),
}));
vi.mock("../providers/delete-dialog", () => ({
  useDeleteDialog: () => vi.fn(),
}));
vi.mock("../utils/platform", () => ({ isElectron: () => false }));

// Render every row: there is no viewport to window against in jsdom.
vi.mock("@tanstack/react-virtual", () => ({
  useVirtualizer: ({ count }: { count: number }) => ({
    getVirtualItems: () =>
      Array.from({ length: count }, (_, index) => ({
        index,
        key: index,
        start: index * 60,
        end: (index + 1) * 60,
        size: 60,
      })),
    getTotalSize: () => count * 60,
    measureElement: () => {},
  }),
}));

import ProjectsTable from "../components/projects-table";

function tagDataTransfer(payload: unknown) {
  return {
    types: [TAG_TYPE],
    dropEffect: "none",
    getData: (type: string) =>
      type === TAG_TYPE ? JSON.stringify(payload) : "",
  };
}

describe("ProjectsTable: dropping a tag on a project", () => {
  beforeEach(() => {
    post.mockClear();
    mutateProjects.mockClear();
    mutateTree.mockClear();
    setMessage.mockClear();
  });

  it("applies the tag to that project and refreshes the list", async () => {
    render(<ProjectsTable />);

    const row = screen.getByText("Nutlin soak 12").closest("tr")!;
    const dataTransfer = tagDataTransfer({ id: 3, label: "Mpro" });
    fireEvent.dragOver(row, { dataTransfer });
    fireEvent.drop(row, { dataTransfer });

    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/3/add_projects/", {
        project_ids: [7],
      })
    );
    await waitFor(() => expect(mutateProjects).toHaveBeenCalled());
    expect(mutateTree).toHaveBeenCalled();
    expect(setMessage).toHaveBeenCalledWith(
      'Tagged "Nutlin soak 12" with "Mpro"',
      "success"
    );
  });

  it("accepts the drop on a card in the card view too", async () => {
    render(<ProjectsTable />);
    fireEvent.click(screen.getByLabelText("card view"));

    const card = (await screen.findByText("Nutlin soak 12")).closest(
      ".MuiCard-root"
    )!;
    expect(card).not.toBeNull();
    const dataTransfer = tagDataTransfer({ id: 3, label: "Mpro" });
    fireEvent.dragOver(card, { dataTransfer });
    fireEvent.drop(card, { dataTransfer });

    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/3/add_projects/", {
        project_ids: [7],
      })
    );
  });

  it("ignores a drop that is not a tag", async () => {
    render(<ProjectsTable />);

    const row = screen.getByText("Nutlin soak 12").closest("tr")!;
    const dataTransfer = { types: ["Files"], dropEffect: "none", getData: () => "" };
    fireEvent.dragOver(row, { dataTransfer });
    fireEvent.drop(row, { dataTransfer });

    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(post).not.toHaveBeenCalled();
  });

  it("reports a rejected drop instead of failing silently", async () => {
    post.mockImplementationOnce(() => Promise.reject(new Error("Tag not found")));
    render(<ProjectsTable />);

    const row = screen.getByText("Nutlin soak 12").closest("tr")!;
    const dataTransfer = tagDataTransfer({ id: 99, label: "gone" });
    fireEvent.dragOver(row, { dataTransfer });
    fireEvent.drop(row, { dataTransfer });

    await waitFor(() =>
      expect(setMessage).toHaveBeenCalledWith(
        'Failed to tag "Nutlin soak 12": Tag not found',
        "error"
      )
    );
    expect(mutateProjects).not.toHaveBeenCalled();
  });
});
