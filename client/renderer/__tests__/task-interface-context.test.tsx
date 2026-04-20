/**
 * Narrow smoke test for the task-interface context + hook plumbing.
 *
 * This does NOT mount real task interfaces — those pull in a deep module graph
 * (monaco-editor, rdkit, moorhen) that is expensive to stub and irrelevant to
 * the behaviour under test here.
 *
 * What it DOES verify, cheaply:
 *   1. <TaskInterfaceProvider> sets up the context with a job + useTaskItem.
 *   2. useTaskInterface() reads those values inside the provider.
 *   3. useTaskInterface() throws outside the provider (clear error for authors
 *      who forget to wrap).
 *   4. useTaskToggles batches CBoolean toggles keyed by item name.
 *
 * Coverage rationale: the wide regression risk in the context/rip refactor is
 * caught by `tsc --noEmit` (removing the `job` prop turns every stale
 * `job={job}` / `{...props}` into a compile error). What TS cannot see is
 * whether the provider is wrapping correctly — that is what this file locks
 * down. One focused test, ~instant to run.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";

import {
  TaskInterfaceProvider,
  useTaskInterface,
} from "../components/task/task-elements/task-interface-context";
import { useTaskToggles } from "../components/task/task-elements/shared-hooks";

const fakeJob: any = {
  id: 42,
  uuid: "test-uuid",
  task_name: "test_task",
  status: 1,
  project: 1,
};

// useJob is re-exported from ../utils; we don't need to mock it — the real
// implementation will kick in, but its SWR fetches will just idle (no backend
// in the test env). TaskInterfaceProvider itself only reads `useTaskItem` off
// the return value, which is a stable memoized function that tolerates
// missing container data.

describe("TaskInterfaceProvider + useTaskInterface", () => {
  it("exposes job through the context", () => {
    const ReadJob = () => {
      const { job } = useTaskInterface();
      return <span data-testid="job-id">{job.id}</span>;
    };
    render(
      <TaskInterfaceProvider job={fakeJob}>
        <ReadJob />
      </TaskInterfaceProvider>
    );
    expect(screen.getByTestId("job-id").textContent).toBe("42");
  });

  it("exposes useTaskItem through the context", () => {
    const ReadItem = () => {
      const { useTaskItem } = useTaskInterface();
      const result = useTaskItem("ANY_ITEM");
      return (
        <span data-testid="item-shape">
          {typeof result === "object" && "value" in result ? "ok" : "broken"}
        </span>
      );
    };
    render(
      <TaskInterfaceProvider job={fakeJob}>
        <ReadItem />
      </TaskInterfaceProvider>
    );
    expect(screen.getByTestId("item-shape").textContent).toBe("ok");
  });

  it("throws a helpful error when used outside a provider", () => {
    const Orphan = () => {
      useTaskInterface();
      return null;
    };
    const spy = vi.spyOn(console, "error").mockImplementation(() => {});
    expect(() => render(<Orphan />)).toThrow(/TaskInterfaceProvider/);
    spy.mockRestore();
  });
});

describe("useTaskToggles", () => {
  it("returns a keyed record with value+onChange for each name", () => {
    let captured: any = null;
    const Inspect = () => {
      captured = useTaskToggles(["FLAG_A", "FLAG_B", "FLAG_C"] as const);
      return null;
    };
    render(
      <TaskInterfaceProvider job={fakeJob}>
        <Inspect />
      </TaskInterfaceProvider>
    );
    expect(Object.keys(captured).sort()).toEqual(["FLAG_A", "FLAG_B", "FLAG_C"]);
    for (const key of ["FLAG_A", "FLAG_B", "FLAG_C"]) {
      expect(typeof captured[key].value).toBe("boolean");
      expect(typeof captured[key].onChange).toBe("function");
    }
  });
});
