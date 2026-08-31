#!/usr/bin/env python3
"""Compose task icons for the Import* suite from their output datatype icons.

Each Import<X> task writes exactly one output datatype, so its icon is that
datatype's icon plus the house import glyph (the green arrow from
``import_arrow_new.svg``) as a badge in the bottom-left corner.

The base icon is embedded as a nested ``<svg>`` element, so this works whether
the base is a vector drawing (AsuDataFile) or a wrapped PNG (most of them).

Run from the repo root::

    python3 scripts/make_import_icons.py

Writes into client/assets/svgicons/ (the tracked copy). ``npm run copy-svgicons``
in client/ mirrors them into the gitignored renderer/public/svgicons/.
"""

import os
import re
import xml.etree.ElementTree as ET

SVGNS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVGNS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")

ICONS = os.path.join(os.path.dirname(__file__), "..", "client", "assets", "svgicons")

# Task -> the CDataFile class it writes to outputData, minus the leading "C".
TASKS = {
    "ImportCoordinate": "PdbDataFile",
    "ImportDictionary": "DictDataFile",
    "ImportSequence": "SeqDataFile",
    "ImportAsuContent": "AsuDataFile",
    "ImportMap": "MapDataFile",
    "ImportObs": "ObsDataFile",
    "ImportMapCoeffs": "MapCoeffsDataFile",
    "ImportFreeR": "FreeRDataFile",
    "ImportPhases": "PhsDataFile",
}

# The arrow from import_arrow_new.svg, normalised to a 560 x 360 box.
ARROW_D = "M 560,120 L 240,120 L 240,0 L 0,200 L 240,360 L 240,240 L 560,240 Z"
ARROW_W = 560.0


def nested_base(path, x, y, size):
    """The base icon as a nested <svg>, scaled to fit (x, y, size, size)."""
    root = ET.parse(path).getroot()
    viewbox = root.get("viewBox")
    if not viewbox:
        width = re.sub(r"[a-z]", "", root.get("width", "64"))
        height = re.sub(r"[a-z]", "", root.get("height", "64"))
        viewbox = f"0 0 {width} {height}"
    nested = ET.Element(
        f"{{{SVGNS}}}svg",
        {
            "x": str(x), "y": str(y), "width": str(size), "height": str(size),
            "viewBox": viewbox, "preserveAspectRatio": "xMidYMid meet",
            "overflow": "hidden",
        },
    )
    for child in root:
        # Inkscape bookkeeping; it renders nothing and only bloats the file.
        if isinstance(child.tag, str) and child.tag.endswith(("}metadata", "}namedview")):
            continue
        nested.append(child)
    return nested


def arrow(x, y, width, stroke_width):
    """The house import arrow, mirrored to point right (into the datatype)."""
    scale = width / ARROW_W
    return ET.Element(
        f"{{{SVGNS}}}path",
        {
            "d": ARROW_D,
            "transform": f"translate({x + width},{y}) scale({-scale},{scale})",
            "fill": "#22c55e", "stroke": "#14532d",
            "stroke-width": str(stroke_width),
            "stroke-linejoin": "round", "stroke-linecap": "round",
            "vector-effect": "non-scaling-stroke",
        },
    )


def build(base_name):
    root = ET.Element(
        f"{{{SVGNS}}}svg",
        {"xmlns:xlink": "http://www.w3.org/1999/xlink",
         "viewBox": "0 0 64 64", "width": "64", "height": "64"},
    )
    root.append(nested_base(os.path.join(ICONS, base_name + ".svg"), 0, 0, 64))
    # A white disc behind the arrow: several bases are dark or busy, and the
    # green-on-black ones are unreadable without it.
    root.append(ET.Element(f"{{{SVGNS}}}circle", {
        "cx": "16", "cy": "48", "r": "15.5",
        "fill": "#ffffff", "stroke": "#14532d", "stroke-width": "2",
    }))
    root.append(arrow(5, 41.5, 22, 2))
    return root


def main():
    for task, base_name in sorted(TASKS.items()):
        out = os.path.join(ICONS, task + ".svg")
        body = ET.tostring(build(base_name), encoding="unicode")
        with open(out, "w") as stream:
            stream.write('<?xml version="1.0" encoding="UTF-8"?>\n' + body + "\n")
        print(f"wrote {task}.svg from {base_name}.svg")


if __name__ == "__main__":
    main()
