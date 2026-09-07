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

vi.mock("next/navigation", () => ({
  useRouter: () => ({ push: mocks.push }),
  usePathname: () => "/",
  useSearchParams: () => new URLSearchParams(),
  useParams: () => ({}),
}));

vi.mock("../api", () => ({
  useApi: () => ({
    get: () => ({ data: undefined }),
    get_endpoint: () => ({ data: TREE, isLoading: false, mutate: mocks.mutate }),
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
