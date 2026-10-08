/**
 * The Program locations page offers only programs whose setting takes effect
 * (#402). Tasks that find their program themselves are named instead, so the
 * page does not imply a setting there would reach them.
 */
import React from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";

vi.mock("../api-fetch", () => ({ apiGet: vi.fn(), apiPatch: vi.fn() }));

import { SelfLocatedNote } from "../components/program-locations";

describe("SelfLocatedNote", () => {
  it("names each self-locating task and how it finds its program", () => {
    render(
      <SelfLocatedNote
        tasks={[
          { task: "arp_warp_classic", title: "ARP/wARP", how: "ARP/wARP's own setup script" },
          { task: "clustalw", title: "ClustalW", how: "CCP4 (libexec/clustalw2)" },
        ]}
      />
    );
    const note = screen.getByTestId("self-located-note");
    expect(note.textContent).toContain("Not configurable here");
    expect(note.textContent).toContain("ARP/wARP (ARP/wARP's own setup script)");
    expect(note.textContent).toContain("ClustalW (CCP4 (libexec/clustalw2))");
  });

  it("renders nothing when every listed program is relocatable", () => {
    const { container } = render(<SelfLocatedNote tasks={[]} />);
    expect(container.textContent).toBe("");
  });
});
