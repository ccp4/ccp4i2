/**
 * The dm_multidomain assembly grid, body editor and coverage strip.
 *
 * These three replace two free-text lists whose meaning could only be conveyed
 * in a paragraph above them. What is tested here is the part that made the
 * paragraph unnecessary: the reference row is labelled rather than implied, a
 * chain that is already spoken for says so rather than silently colliding, a
 * body reports its fit against each copy, and two bodies claiming the same
 * residues are drawn as a clash.
 *
 * All three are presentational — model in, next model out — so they render
 * with no job, no server and no mocks.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";

import {
  Model,
  parseAssemblyRows,
  parseBodies,
} from "../components/task/task-interfaces/dm_multidomain/dm-spec";
import { AssemblyGrid } from "../components/task/task-interfaces/dm_multidomain/assembly-grid";
import { BodiesEditor } from "../components/task/task-interfaces/dm_multidomain/bodies-editor";
import { CoverageStrip } from "../components/task/task-interfaces/dm_multidomain/coverage-strip";

const CHAINS = [
  { id: "A", first: 1, last: 298 },
  { id: "B", first: 175, last: 432 },
  { id: "C", first: 1, last: 298 },
  { id: "D", first: 175, last: 432 },
];
const ENTITIES = [
  ["A", "C"],
  ["B", "D"],
];

const heteroModel = (): Model => ({
  assembly: parseAssemblyRows(["CDK=A cyclin=B", "CDK=C cyclin=D"]),
  bodies: parseBodies([{ segments: "CDK:1-298", mode: "average" }]),
});

const homomerModel = (): Model => ({
  assembly: parseAssemblyRows(["A", "C"]),
  bodies: parseBodies([{ segments: "1-298", mode: "average" }]),
});

const grid = (model: Model, onChange = vi.fn()) => {
  render(
    <AssemblyGrid
      model={model}
      chains={CHAINS}
      entities={ENTITIES}
      onChange={onChange}
    />
  );
  return onChange;
};

describe("assembly grid", () => {
  it("labels the first row as the reference instead of implying it", () => {
    grid(heteroModel());
    expect(screen.getByText("reference")).toBeTruthy();
    expect(screen.getByText("copy 2")).toBeTruthy();
  });

  it("shows each entity as a column heading", () => {
    grid(heteroModel());
    expect(screen.getByText("CDK")).toBeTruthy();
    expect(screen.getByText("cyclin")).toBeTruthy();
  });

  it("names the unnamed entity after the reference chain", () => {
    grid(homomerModel());
    expect(screen.getByText("entity A")).toBeTruthy();
  });

  it("says why a chain cannot be picked rather than hiding it", () => {
    grid(heteroModel());
    // open the reference copy's CDK cell
    fireEvent.mouseDown(screen.getAllByRole("combobox")[0]);
    const option = screen.getByRole("option", { name: /chain C/ });
    expect(option.getAttribute("aria-disabled")).toBe("true");
    expect(within(option).getByText(/already in copy 2/)).toBeTruthy();
  });

  it("writes the chosen chain back into the model", () => {
    const model = heteroModel();
    const onChange = grid(model);
    fireEvent.mouseDown(screen.getAllByRole("combobox")[0]);
    fireEvent.click(screen.getByRole("option", { name: /not in this copy/ }));
    expect(onChange).toHaveBeenCalledTimes(1);
    const next: Model = onChange.mock.calls[0][0];
    expect(next.assembly.instances[0]).toEqual({ cyclin: "B" });
  });

  it("adds a copy as a visible empty row", () => {
    const onChange = grid(heteroModel());
    fireEvent.click(screen.getByRole("button", { name: /add copy/i }));
    const next: Model = onChange.mock.calls[0][0];
    expect(next.assembly.instances).toHaveLength(3);
    expect(next.assembly.instances[2]).toEqual({});
  });
});

describe("bodies editor", () => {
  const preview = [
    {
      index: 1,
      spec: "CDK:1-298",
      mode: "average",
      segments: [{ role: "CDK", lo: 1, hi: 298 }],
      copies: [{ label: "C+D", nCA: 298, rmsd: 0.84 }],
      error: null,
    },
  ];

  const editor = (model: Model, onChange = vi.fn(), withPreview = true) => {
    render(
      <BodiesEditor
        bodies={model.bodies}
        assembly={model.assembly}
        chains={CHAINS}
        preview={withPreview ? (preview as any) : undefined}
        onChange={onChange}
      />
    );
    return onChange;
  };

  it("reports the fit of the body against each copy", () => {
    editor(heteroModel());
    // the number that says whether these residues really move as one unit
    expect(screen.getByText("C+D: 298 CA, 0.84 Å")).toBeTruthy();
  });

  it("bounds a residue range by the reference chain's own numbering", () => {
    editor(heteroModel());
    const first = screen.getByLabelText(/segment 1 first residue/);
    expect(first.getAttribute("min")).toBe("1");
    expect(first.getAttribute("max")).toBe("298");
    expect(screen.getByText("of 1–298")).toBeTruthy();
  });

  it("clamps a residue typed outside the chain", () => {
    const onChange = editor({
      assembly: heteroModel().assembly,
      bodies: parseBodies([{ segments: "CDK:1-100", mode: "average" }]),
    });
    const last = screen.getByLabelText(/segment 1 last residue/);
    fireEvent.change(last, { target: { value: "9999" } });
    fireEvent.blur(last);
    const next = onChange.mock.calls[0][0];
    expect(next[0].segments[0].hi).toBe(298);
  });

  it("changes the averaging mode", () => {
    const onChange = editor(heteroModel());
    fireEvent.click(screen.getByRole("button", { name: "exclude" }));
    expect(onChange.mock.calls[0][0][0].mode).toBe("exclude");
  });

  it("offers a role per entity only when there is more than one", () => {
    const { unmount } = render(
      <BodiesEditor
        bodies={homomerModel().bodies}
        assembly={homomerModel().assembly}
        chains={CHAINS}
        onChange={vi.fn()}
      />
    );
    expect(screen.queryAllByRole("combobox")).toHaveLength(0);
    unmount();
    editor(heteroModel());
    expect(screen.getAllByRole("combobox").length).toBeGreaterThan(0);
  });

  it("shows an unreadable stored spec rather than an empty body", () => {
    const { container } = render(
      <BodiesEditor
        bodies={parseBodies([{ segments: "CDK:forty-two", mode: "average" }])}
        assembly={heteroModel().assembly}
        chains={CHAINS}
        onChange={vi.fn()}
      />
    );
    expect(container.querySelector("code")?.textContent).toBe("CDK:forty-two");
  });
});

describe("coverage strip", () => {
  it("draws a track per entity of the reference copy", () => {
    const model = heteroModel();
    const { container } = render(
      <CoverageStrip
        assembly={model.assembly}
        bodies={model.bodies}
        chains={CHAINS}
      />
    );
    expect(screen.getByText("CDK")).toBeTruthy();
    expect(screen.getByText("cyclin")).toBeTruthy();
    // the residue extent of the chain, so a gap is visible as a gap
    expect(screen.getByText("298")).toBeTruthy();
    expect(container.querySelectorAll("rect").length).toBeGreaterThan(2);
  });

  it("draws two bodies claiming the same residues as a clash", () => {
    const model = heteroModel();
    const { container } = render(
      <CoverageStrip
        assembly={model.assembly}
        bodies={parseBodies([
          { segments: "CDK:1-298", mode: "average" },
          { segments: "CDK:1-100", mode: "average" },
        ])}
        chains={CHAINS}
      />
    );
    const clash = [...container.querySelectorAll("rect")].filter(
      (r) => r.getAttribute("fill") === "url(#dm-clash)"
    );
    expect(clash).toHaveLength(1);
  });

  it("renders nothing when the model is not loaded yet", () => {
    const model = heteroModel();
    const { container } = render(
      <CoverageStrip assembly={model.assembly} bodies={model.bodies} chains={[]} />
    );
    expect(container.querySelector("svg")).toBeNull();
  });
});
