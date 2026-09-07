/**
 * Regression tests for the tag tree's editing affordances.
 *
 * These exist because of a specific bug: "Add tag inside" and "Rename" in a
 * node's context menu appeared to do nothing. The field they open focuses
 * itself, and MUI's focus trap stays active for the whole of the menu's exit
 * transition — it pulled focus back out of the new field, blurring it, which
 * cancelled it a frame after it appeared. The same actions reached by the
 * inline "+" button (no menu, so no trap) worked fine, which is what made the
 * fault look like a missing handler rather than a race.
 *
 * So the control case and the failing case are both locked down here.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen, waitFor } from "@testing-library/react";

const post = vi.fn(() => Promise.resolve({ id: 99 }));
const patch = vi.fn(() => Promise.resolve({}));
const del = vi.fn(() => Promise.resolve());
const mutate = vi.fn(() => Promise.resolve());

const SEP = "\u001f";

// A root with one nested child, plus a second root to move things onto:
//
//   Nutlin site
//     soaks
//   Mpro
const TREE = {
  tags: [
    {
      id: 1,
      text: "Nutlin site",
      parent: null,
      path: "Nutlin site",
      display_path: "Nutlin site",
      depth: 0,
      project_count: 2,
      total_project_count: 2,
      children: [
        {
          id: 2,
          text: "soaks",
          parent: 1,
          path: `Nutlin site${SEP}soaks`,
          display_path: "Nutlin site/soaks",
          depth: 1,
          project_count: 0,
          total_project_count: 0,
          children: [],
        },
      ],
    },
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
  untagged_project_count: 3,
};

const TAG_TYPE = "application/x-ccp4i2-tag";

/**
 * A stand-in for the browser's DataTransfer: `types` is what dragover can
 * read, and setData/getData carry the payload between dragstart and drop.
 */
function fakeDataTransfer(types: string[] = []) {
  const store: Record<string, string> = {};
  return {
    types,
    effectAllowed: "uninitialized",
    dropEffect: "none",
    setData: (type: string, value: string) => {
      store[type] = value;
      if (!types.includes(type)) types.push(type);
    },
    getData: (type: string) => store[type] ?? "",
  };
}

/** Drag the row labelled `from` and drop it on the row labelled `onto`. */
function dragTag(from: string, onto: string) {
  const dataTransfer = fakeDataTransfer();
  const source = screen.getByText(from);
  fireEvent.dragStart(source, { dataTransfer });
  const target = screen.getByText(onto);
  fireEvent.dragOver(target, { dataTransfer });
  fireEvent.drop(target, { dataTransfer });
  fireEvent.dragEnd(source, { dataTransfer });
  return dataTransfer;
}

vi.mock("../api", () => ({
  useApi: () => ({
    get: () => ({ data: TREE, mutate }),
    post,
    patch,
    delete: del,
  }),
}));

const setMessage = vi.fn();
vi.mock("../providers/popcorn-provider", () => ({
  usePopcorn: () => ({ setMessage }),
}));

vi.mock("../providers/delete-dialog", () => ({
  useDeleteDialog: () => vi.fn(),
}));

import ProjectTagTreePane from "../components/project-tag-tree";

function renderPane() {
  return render(
    <ProjectTagTreePane
      selected={{ kind: "all" }}
      onSelect={vi.fn()}
      totalCount={5}
    />
  );
}

describe("ProjectTagTreePane editing", () => {
  beforeEach(() => {
    post.mockClear();
    patch.mockClear();
    mutate.mockClear();
    setMessage.mockClear();
  });

  it("renders the forest with rolled-up counts", () => {
    renderPane();
    expect(screen.getByText("Nutlin site")).toBeInTheDocument();
    expect(screen.getByText("All projects")).toBeInTheDocument();
    expect(screen.getByText("Untagged")).toBeInTheDocument();
  });

  it("opens a name field from the inline + button (the control case)", async () => {
    renderPane();

    fireEvent.click(screen.getByLabelText("Add a tag inside Nutlin site"));

    expect(
      await screen.findByPlaceholderText("New tag name")
    ).toBeInTheDocument();
  });

  it("opens a name field from the context menu, and it survives the menu closing", async () => {
    renderPane();

    fireEvent.click(screen.getByLabelText("More actions for Nutlin site"));
    fireEvent.click(await screen.findByText("Add tag inside"));

    // The field must still be there once the menu has finished closing —
    // this is the assertion the focus-trap race used to fail.
    const field = await screen.findByPlaceholderText("New tag name");
    await waitFor(() =>
      expect(screen.queryByText("Add tag inside")).not.toBeInTheDocument()
    );
    expect(field).toBeInTheDocument();

    fireEvent.change(field, { target: { value: "soaks" } });
    fireEvent.keyDown(field, { key: "Enter" });
    await waitFor(() =>
      expect(post).toHaveBeenCalledWith("projecttags/", {
        text: "soaks",
        parent: 1,
        projects: [],
      })
    );
  });

  it("accepts a dropped project selection and reports the target tag", async () => {
    const onDropProjects = vi.fn();
    render(
      <ProjectTagTreePane
        selected={{ kind: "all" }}
        onSelect={vi.fn()}
        onDropProjects={onDropProjects}
        totalCount={5}
      />
    );

    const node = screen.getByText("Nutlin site");
    const dataTransfer = {
      types: ["application/x-ccp4i2-projects"],
      getData: () => "[1,2]",
    };
    fireEvent.dragOver(node, { dataTransfer });
    fireEvent.drop(node, { dataTransfer });

    expect(onDropProjects).toHaveBeenCalledWith(1, "Nutlin site");
  });

  it("ignores a drag that is not a project selection", async () => {
    const onDropProjects = vi.fn();
    render(
      <ProjectTagTreePane
        selected={{ kind: "all" }}
        onSelect={vi.fn()}
        onDropProjects={onDropProjects}
        totalCount={5}
      />
    );

    const node = screen.getByText("Nutlin site");
    const dataTransfer = { types: ["Files"], getData: () => "" };
    fireEvent.dragOver(node, { dataTransfer });
    fireEvent.drop(node, { dataTransfer });

    expect(onDropProjects).not.toHaveBeenCalled();
  });

  it("re-parents a tag dropped onto another tag", async () => {
    renderPane();

    const dataTransfer = dragTag("Nutlin site", "Mpro");

    expect(dataTransfer.types).toContain(TAG_TYPE);
    expect(JSON.parse(dataTransfer.getData(TAG_TYPE))).toEqual({
      id: 1,
      label: "Nutlin site",
    });
    await waitFor(() =>
      expect(patch).toHaveBeenCalledWith("projecttags/1/", { parent: 3 })
    );
    expect(mutate).toHaveBeenCalled();
  });

  it("makes a nested tag a root when it is dropped on All projects", async () => {
    renderPane();

    dragTag("soaks", "All projects");

    await waitFor(() =>
      expect(patch).toHaveBeenCalledWith("projecttags/2/", { parent: null })
    );
  });

  it("does nothing when a root tag is dropped on All projects (it already is one)", async () => {
    renderPane();
    dragTag("Nutlin site", "All projects");
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(patch).not.toHaveBeenCalled();
  });

  it("refuses to drop a tag on its own descendant", async () => {
    renderPane();
    dragTag("Nutlin site", "soaks");
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(patch).not.toHaveBeenCalled();
  });

  it("refuses to drop a tag on itself", async () => {
    renderPane();
    dragTag("Mpro", "Mpro");
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(patch).not.toHaveBeenCalled();
  });

  it("cleans up the drag ghost when the drag ends", () => {
    renderPane();
    const dataTransfer = { ...fakeDataTransfer(), setDragImage: vi.fn() };
    const source = screen.getByText("Mpro");

    fireEvent.dragStart(source, { dataTransfer });
    expect(dataTransfer.setDragImage).toHaveBeenCalledTimes(1);
    const ghost = dataTransfer.setDragImage.mock.calls[0][0] as HTMLElement;
    // The ghost is just the tag's name, not a snapshot of the list.
    expect(ghost.textContent).toBe("Mpro");
    expect(document.body.contains(ghost)).toBe(true);

    fireEvent.dragEnd(source, { dataTransfer });
    expect(document.body.contains(ghost)).toBe(false);
  });

  it("reports the server's reason when a move is rejected", async () => {
    patch.mockImplementationOnce(() =>
      Promise.reject(new Error("A tag with this text and parent already exists."))
    );
    renderPane();

    dragTag("soaks", "Mpro");

    await waitFor(() =>
      expect(setMessage).toHaveBeenCalledWith(
        "Could not move that tag: A tag with this text and parent already exists.",
        "error"
      )
    );
  });

  it("opens a rename field from the context menu, prefilled and surviving the close", async () => {
    renderPane();

    fireEvent.click(screen.getByLabelText("More actions for Nutlin site"));
    fireEvent.click(await screen.findByText("Rename"));

    const field = await screen.findByDisplayValue("Nutlin site");
    await waitFor(() =>
      expect(screen.queryByText("Rename")).not.toBeInTheDocument()
    );
    expect(field).toBeInTheDocument();

    fireEvent.change(field, { target: { value: "Nutlin" } });
    fireEvent.keyDown(field, { key: "Enter" });
    await waitFor(() =>
      expect(patch).toHaveBeenCalledWith("projecttags/1/", { text: "Nutlin" })
    );
  });
});
