/**
 * The job list's "select jobs to delete" mode.
 *
 * Only top-level jobs can be selected in that mode, so the per-row controls
 * that make sense outside it — the expand/collapse chevron and the burger
 * menu — are out of place inside it. Worse, the chevron's click used to fall
 * through to the tree's selection handler and tick the checkbox while also
 * expanding the sub-jobs (alpha tracker, Paul Bond). These tests lock down
 * the intended shape of a row in each mode.
 *
 * The heavy neighbours (API, context menus, avatars needing the theme
 * provider, dnd-kit) are stubbed; the tree view and the row component under
 * test are real.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";

const mocks = vi.hoisted(() => ({
  push: vi.fn(),
  setJobMenuAnchorEl: vi.fn(),
  setJob: vi.fn(),
  post: vi.fn(() => Promise.resolve({})),
  mutate: vi.fn(() => Promise.resolve()),
  params: {} as Record<string, string>,
  tree: undefined as unknown,
}));

const CHILD = {
  serverjob: 0,
  float_values: [],
  char_values: [],
  file_uses: [],
  xdatas: [],
  id: 11,
  uuid: "job-1-1",
  project: 7,
  parent: 1,
  number: "1.1",
  title: "Refmac step",
  status: 6,
  evaluation: 0,
  comments: "",
  creation_time: "",
  finish_time: "",
  task_name: "refmac",
  process_id: 0,
  files: [],
  kpis: { float_values: {}, char_values: {} },
  children: [],
};

const TOP = {
  ...CHILD,
  id: 1,
  uuid: "job-1",
  parent: null,
  number: "1",
  title: "Refine",
  task_name: "servalcat_pipe",
  children: [CHILD],
};

const TREE = { job_tree: [TOP], total_jobs: 2, total_files: 0 };

const topLevel = (id: number, title: string) => ({
  ...TOP,
  id,
  uuid: `job-${id}`,
  number: String(id),
  title,
  children: [],
});
const FOUR_JOBS = {
  job_tree: [topLevel(4, "Four"), topLevel(3, "Three"), topLevel(2, "Two"), topLevel(1, "One")],
  total_jobs: 4,
  total_files: 0,
};

vi.mock("next/navigation", () => ({
  useRouter: () => ({ push: mocks.push }),
  usePathname: () => "/",
  useSearchParams: () => new URLSearchParams(),
  useParams: () => mocks.params,
}));

vi.mock("../api", () => ({
  useApi: () => ({
    get: () => ({ data: undefined }),
    get_endpoint: () => ({ data: mocks.tree, isLoading: false, mutate: mocks.mutate }),
    post: mocks.post,
  }),
}));

vi.mock("../providers/job-context-menu", () => ({
  useJobMenu: () => ({
    setJobMenuAnchorEl: mocks.setJobMenuAnchorEl,
    setJob: mocks.setJob,
  }),
}));

vi.mock("../providers/file-context-menu", () => ({
  useFileMenu: () => ({ setFileMenuAnchorEl: vi.fn(), setFile: vi.fn() }),
}));

vi.mock("../providers/recently-started-jobs-context", () => ({
  useRecentlyStartedJobs: () => ({
    hasRecentlyStartedJobsInProject: () => false,
  }),
}));

vi.mock("../providers/delete-dialog", () => ({
  useDeleteDialog: () => vi.fn(),
}));

// The avatars need the app theme provider; they are decoration here.
vi.mock("../components/job-avatar", () => ({
  CCP4i2JobAvatar: React.forwardRef<HTMLSpanElement, any>(function Avatar(_p, ref) {
    return <span ref={ref} data-testid="job-avatar" />;
  }),
}));
vi.mock("../components/file-avatar", () => ({
  FileAvatar: () => <span data-testid="file-avatar" />,
}));

vi.mock("@dnd-kit/core", () => ({
  useDraggable: () => ({ attributes: {}, listeners: {}, setNodeRef: () => {} }),
}));

import { ClassicJobList, EXPAND_TOGGLE_ROLE } from "../components/classic-jobs-list";

const TOGGLE = `[data-role="${EXPAND_TOGGLE_ROLE}"]`;

function renderList() {
  return render(<ClassicJobList projectId={7} />);
}

function enterSelectMode() {
  fireEvent.click(screen.getByLabelText("Select jobs to delete"));
}

beforeEach(() => {
  mocks.push.mockClear();
  mocks.params = {};
  mocks.tree = TREE;
  mocks.setJobMenuAnchorEl.mockClear();
  mocks.setJob.mockClear();
});

describe("job list outside select mode", () => {
  it("shows a chevron and a burger on the row", () => {
    const { container } = renderList();
    expect(screen.getByText("1: Refine")).toBeInTheDocument();
    expect(container.querySelectorAll(TOGGLE)).toHaveLength(1);
    expect(screen.getAllByLabelText("Open job menu")).toHaveLength(1);
    expect(screen.queryByRole("checkbox")).toBeNull();
  });

  it("clicking the chevron expands the sub-jobs without selecting the row", () => {
    const { container } = renderList();
    expect(screen.queryByText("1.1: Refmac step")).toBeNull();

    fireEvent.click(container.querySelector(TOGGLE)!);

    expect(screen.getByText("1.1: Refmac step")).toBeInTheDocument();
    expect(mocks.push).not.toHaveBeenCalled();
  });

  it("clicking the row itself selects it (navigates)", () => {
    renderList();
    fireEvent.click(screen.getByText("1: Refine"));
    expect(mocks.push).toHaveBeenCalledWith("/ccp4i2/project/7/job/1");
  });
});

describe("job list in select mode", () => {
  it("hides the chevron and the burger, and offers a checkbox for top-level jobs only", () => {
    const { container } = renderList();
    // Expand first so a sub-job row is on screen when the mode is entered.
    fireEvent.click(container.querySelector(TOGGLE)!);
    expect(screen.getByText("1.1: Refmac step")).toBeInTheDocument();

    enterSelectMode();

    expect(container.querySelectorAll(TOGGLE)).toHaveLength(0);
    expect(screen.queryAllByLabelText("Open job menu")).toHaveLength(0);
    // One checkbox: the top-level job. The sub-job row has none.
    expect(screen.getAllByRole("checkbox")).toHaveLength(1);
  });

  it("a click anywhere on a top-level row toggles its checkbox and does not navigate", () => {
    renderList();
    enterSelectMode();
    const checkbox = screen.getByRole("checkbox") as HTMLInputElement;
    expect(checkbox.checked).toBe(false);

    fireEvent.click(screen.getByText("1: Refine"));
    expect(checkbox.checked).toBe(true);
    expect(screen.getByText("1 job selected")).toBeInTheDocument();

    fireEvent.click(screen.getByText("1: Refine"));
    expect(checkbox.checked).toBe(false);

    expect(mocks.push).not.toHaveBeenCalled();
  });

  it("clicking a sub-job row does nothing", () => {
    const { container } = renderList();
    fireEvent.click(container.querySelector(TOGGLE)!);
    enterSelectMode();

    fireEvent.click(screen.getByText("1.1: Refmac step"));

    expect((screen.getByRole("checkbox") as HTMLInputElement).checked).toBe(false);
    expect(screen.queryByText("1 job selected")).toBeNull();
    expect(mocks.push).not.toHaveBeenCalled();
  });

  it("right-click and double-click open nothing", () => {
    const open = vi
      .spyOn(window, "open")
      .mockImplementation(() => null as unknown as ReturnType<typeof window.open>);
    renderList();
    enterSelectMode();
    const row = screen.getByText("1: Refine");

    fireEvent.contextMenu(row);
    expect(mocks.setJobMenuAnchorEl).not.toHaveBeenCalled();
    expect(mocks.setJob).not.toHaveBeenCalled();

    fireEvent.doubleClick(row);
    expect(open).not.toHaveBeenCalled();
    open.mockRestore();
  });

  it("the chevron and burger come back on leaving select mode", () => {
    const { container } = renderList();
    enterSelectMode();
    expect(container.querySelectorAll(TOGGLE)).toHaveLength(0);

    fireEvent.click(screen.getByRole("button", { name: "Cancel" }));

    expect(container.querySelectorAll(TOGGLE)).toHaveLength(1);
    expect(screen.getAllByLabelText("Open job menu")).toHaveLength(1);
    expect(screen.queryByRole("checkbox")).toBeNull();
  });
});

describe("job list highlight", () => {
  const highlighted = () =>
    screen
      .getAllByRole("treeitem")
      .filter((item) => item.getAttribute("aria-selected") === "true")
      .map((item) => item.textContent);

  it("follows the job in the route", () => {
    mocks.tree = FOUR_JOBS;
    mocks.params = { id: "7", jobid: "2" };
    renderList();
    expect(highlighted()).toEqual([expect.stringContaining("2: Two")]);
  });

  it("is empty when the route names no job", () => {
    mocks.tree = FOUR_JOBS;
    mocks.params = { id: "7" };
    renderList();
    expect(highlighted()).toEqual([]);
  });
});

describe("selection shortcuts", () => {
  const checked = () =>
    screen
      .getAllByRole("checkbox")
      .filter((box) => (box as HTMLInputElement).checked)
      .map((box) => box.closest("li")!.textContent);

  beforeEach(() => {
    mocks.tree = FOUR_JOBS;
    mocks.params = { id: "7", jobid: "3" };
  });

  it("ctrl-click enters select mode with the open job and the clicked job ticked", () => {
    renderList();
    fireEvent.click(screen.getByText("1: One"), { ctrlKey: true });
    expect(screen.getByText("2 jobs selected")).toBeInTheDocument();
    expect(checked()).toEqual([
      expect.stringContaining("3: Three"),
      expect.stringContaining("1: One"),
    ]);
    expect(mocks.push).not.toHaveBeenCalled();
  });

  it("ctrl-click toggles within select mode", () => {
    renderList();
    fireEvent.click(screen.getByText("1: One"), { ctrlKey: true });
    fireEvent.click(screen.getByText("1: One"), { ctrlKey: true });
    expect(screen.getByText("1 job selected")).toBeInTheDocument();
  });

  it("shift-click ticks the range from the open job", () => {
    renderList();
    fireEvent.click(screen.getByText("1: One"), { shiftKey: true });
    expect(screen.getByText("3 jobs selected")).toBeInTheDocument();
    expect(checked()).not.toContainEqual(expect.stringContaining("4: Four"));
  });

  it("shift-click within select mode extends from the last toggled job", () => {
    renderList();
    enterSelectMode();
    fireEvent.click(screen.getByText("4: Four"));
    fireEvent.click(screen.getByText("2: Two"), { shiftKey: true });
    expect(checked()).toEqual([
      expect.stringContaining("4: Four"),
      expect.stringContaining("3: Three"),
      expect.stringContaining("2: Two"),
    ]);
  });

  it("Escape leaves select mode", () => {
    renderList();
    fireEvent.click(screen.getByText("1: One"), { ctrlKey: true });
    fireEvent.keyDown(screen.getByText("1: One"), { key: "Escape" });
    expect(screen.queryByRole("checkbox")).toBeNull();
    expect(screen.getByPlaceholderText("Search jobs…")).toBeInTheDocument();
  });
});
