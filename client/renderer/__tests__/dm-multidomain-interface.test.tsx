/**
 * dm_multidomain opens on something that works.
 *
 * The old interface opened on an empty assembly list and a blank rigid body,
 * and explained in prose that leaving the list empty meant "detect from the
 * model". That is magic by absence: the user is shown a blank box and told
 * what the blank means. Here the detection is written in, where it can be read
 * and corrected — so these tests are about the prefill and the two writes it
 * makes, which are the only part of the interface that touches the parameters
 * on its own.
 */
import React from "react";
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";

const mocks = vi.hoisted(() => ({
  assembly: [] as any[],
  domains: [] as any[],
  replaceAssembly: vi.fn(async () => undefined),
  replaceDomains: vi.fn(async () => undefined),
  preview: null as any,
  callPluginMethod: vi.fn(),
  xyzin: { dbFileId: "deadbeef" } as any,
}));

vi.mock("../utils", () => ({
  useJob: () => ({
    useTaskItem: (name: string) => ({
      value: name === "XYZIN" ? mocks.xyzin : name === "PHASE_SOURCE" ? "input" : null,
      item: null,
    }),
    callPluginMethod: mocks.callPluginMethod,
  }),
}));

vi.mock("../components/task/task-elements/hooks/useContainerList", () => ({
  useContainerList: ({ itemName }: any) =>
    itemName === "ASSEMBLY"
      ? { items: mocks.assembly, replaceArray: mocks.replaceAssembly }
      : { items: mocks.domains, replaceArray: mocks.replaceDomains },
}));

vi.mock("../components/task/task-elements/task-element", () => ({
  CCP4i2TaskElement: () => null,
}));
vi.mock("../components/task/task-elements/ccontainer", () => ({
  CCP4i2ContainerElement: ({ children }: any) => <div>{children}</div>,
}));
vi.mock("../components/task/task-elements/inline-field", () => ({
  InlineField: ({ children }: any) => <div>{children}</div>,
}));
vi.mock("../components/task/task-elements/tabs", () => ({
  // render every tab, so the Domains tab is reachable without clicking
  CCP4i2Tabs: ({ children }: any) => <div>{children}</div>,
  CCP4i2Tab: ({ children }: any) => <div>{children}</div>,
}));

import TaskInterface from "../components/task/task-interfaces/dm_multidomain";

const PREVIEW = {
  ok: true,
  model: {
    chains: [
      { id: "A", first: 1, last: 298, entity: 0 },
      { id: "B", first: 175, last: 432, entity: 1 },
      { id: "C", first: 1, last: 298, entity: 0 },
      { id: "D", first: 175, last: 432, entity: 1 },
    ],
    entities: [["A", "C"], ["B", "D"]],
    nCopiesDetected: 2,
  },
  suggestion: {
    assembly: ["A=A B=B", "A=C B=D"],
    segments: "A:1-298,B:175-432",
  },
  assembly: { source: "detected", rows: [], instances: [] },
  bodies: [],
  messages: [],
};

const job: any = { id: 7, status: 1, task_name: "dm_multidomain" };

const stringRow = (value: string) => ({ _value: value });
const bodyRow = (segments: string, mode = "average") => ({
  _value: { segments: { _value: segments }, mode: { _value: mode } },
});

beforeEach(() => {
  vi.clearAllMocks();
  mocks.assembly = [];
  mocks.domains = [bodyRow("")];          // a fresh job's one blank body
  mocks.xyzin = { dbFileId: "deadbeef" };
  mocks.callPluginMethod.mockResolvedValue(PREVIEW);
});

describe("opening a new job", () => {
  it("writes in the assembly the model implies, rather than leaving it blank", async () => {
    render(<TaskInterface job={job} />);
    await waitFor(() => expect(mocks.replaceAssembly).toHaveBeenCalled());
    expect(mocks.replaceAssembly).toHaveBeenCalledWith(["A=A B=B", "A=C B=D"]);
  });

  it("writes in a body covering the whole copy", async () => {
    render(<TaskInterface job={job} />);
    await waitFor(() => expect(mocks.replaceDomains).toHaveBeenCalled());
    expect(mocks.replaceDomains).toHaveBeenCalledWith([
      { segments: "A:1-298,B:175-432", mode: "average" },
    ]);
  });

  it("does not overwrite an assembly the user has already set", async () => {
    mocks.assembly = [stringRow("A=A B=B"), stringRow("A=C B=D")];
    mocks.domains = [bodyRow("A:1-100")];
    render(<TaskInterface job={job} />);
    await waitFor(() => expect(mocks.callPluginMethod).toHaveBeenCalled());
    expect(mocks.replaceAssembly).not.toHaveBeenCalled();
    expect(mocks.replaceDomains).not.toHaveBeenCalled();
  });

  it("asks for nothing and writes nothing until a model is chosen", async () => {
    mocks.xyzin = null;
    render(<TaskInterface job={job} />);
    expect(
      screen.getByText(/Choose a model on the Input Data tab/)
    ).toBeTruthy();
    await waitFor(() => expect(mocks.callPluginMethod).not.toHaveBeenCalled());
    expect(mocks.replaceAssembly).not.toHaveBeenCalled();
  });

  it("leaves the parameters alone on a job that has already run", async () => {
    render(<TaskInterface job={{ ...job, status: 3 }} />);
    await waitFor(() => expect(mocks.callPluginMethod).toHaveBeenCalled());
    expect(mocks.replaceAssembly).not.toHaveBeenCalled();
  });

  it("says what it found in the model", async () => {
    render(<TaskInterface job={job} />);
    await waitFor(() =>
      expect(
        screen.getByText(/4 protein chains in 2 entities; 2 copies/)
      ).toBeTruthy()
    );
  });

  it("says why it cannot help when the model cannot be read", async () => {
    mocks.callPluginMethod.mockResolvedValue({ ok: false, error: "no such file" });
    render(<TaskInterface job={job} />);
    await waitFor(() =>
      expect(screen.getByText(/Could not read the model: no such file/)).toBeTruthy()
    );
    expect(mocks.replaceAssembly).not.toHaveBeenCalled();
  });
});
