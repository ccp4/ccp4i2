/**
 * The tag chips on a project row or card show a nested tag by its full path,
 * the same string the tree filter label and the tag pickers use — so a tag
 * reads the same everywhere, and two nested tags that share a leaf name can
 * be told apart without hovering for a tooltip.
 */
import React from "react";
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";

import { ProjectTagChips } from "../components/project-tag-chips";

const SEP = "\u001f";

const TAGS = [
  {
    id: 1,
    text: "Nutlin site",
    parent: null as unknown as number,
    children: [],
    projects: [],
    path: "Nutlin site",
    display_path: "Nutlin site",
  },
  {
    id: 2,
    text: "soaks",
    parent: 1,
    children: [],
    projects: [],
    path: `Nutlin site${SEP}soaks`,
    display_path: "Nutlin site/soaks",
  },
];

describe("ProjectTagChips", () => {
  it("labels a nested tag with its full path", () => {
    render(<ProjectTagChips tags={TAGS} />);
    expect(screen.getByText("Nutlin site")).toBeInTheDocument();
    expect(screen.getByText("Nutlin site/soaks")).toBeInTheDocument();
    expect(screen.queryByText("soaks")).not.toBeInTheDocument();
  });

  it("falls back to the text when a tag carries no display path", () => {
    render(
      <ProjectTagChips
        tags={[{ id: 5, text: "plain", parent: 0, children: [], projects: [] }]}
      />
    );
    expect(screen.getByText("plain")).toBeInTheDocument();
  });

  it("folds tags beyond maxVisible into a +N chip", () => {
    render(<ProjectTagChips tags={TAGS} maxVisible={1} />);
    expect(screen.getByText("Nutlin site")).toBeInTheDocument();
    expect(screen.getByText("+1")).toBeInTheDocument();
    expect(screen.queryByText("Nutlin site/soaks")).not.toBeInTheDocument();
  });
});
