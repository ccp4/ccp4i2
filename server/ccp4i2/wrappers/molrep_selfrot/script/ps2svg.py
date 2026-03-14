"""Convert molrep self-rotation function PostScript to SVG.

Parses the simple PS output from molrep (line segments, text labels,
line widths) and produces an SVG string suitable for embedding in reports.

The PS uses a small set of custom operators:
  /L  { moveto lineto stroke }          - draw a line segment
  /Print { ... findfont ... show }      - draw text (Times-Roman)
  /Center { ... }                       - centre text then draw
  moveto (string) <size> Print          - positioned text
  moveto (string) <size> Center (string) <size> Print - centred text
  <width> setlinewidth                  - change stroke width
"""

import re
from xml.sax.saxutils import escape as xml_escape


# Contour line widths used by molrep, mapped to visual stroke widths and colours.
# Thicker PS lines = higher contour level = darker colour.
CONTOUR_STYLES = {
    0.10: {'stroke': '#888888', 'width': 0.4},   # grid lines
    0.20: {'stroke': '#888888', 'width': 0.4},   # grid lines
    0.60: {'stroke': '#2166ac', 'width': 0.8},   # lowest contour
    0.80: {'stroke': '#4393c3', 'width': 1.0},
    1.00: {'stroke': '#d6604d', 'width': 1.2},
    1.25: {'stroke': '#b2182b', 'width': 1.5},
    1.50: {'stroke': '#67001f', 'width': 1.8},   # highest contour
    2.00: {'stroke': '#333333', 'width': 1.2},   # outer circle boundary
}

# Grid line widths (the polar grid, radial lines, outer circle)
GRID_WIDTHS = {0.10, 0.20, 0.80, 2.00}


def _parse_ps(ps_text):
    """Parse molrep PS into a list of drawing operations.

    Returns a list of dicts, each one of:
      {'op': 'line', 'x1':, 'y1':, 'x2':, 'y2':, 'linewidth':}
      {'op': 'text', 'x':, 'y':, 'text':, 'size':, 'centered': bool}
    """
    ops = []
    current_linewidth = 1.0
    # PS coordinate system: origin bottom-left, y increases upward.
    # We'll flip y during SVG generation.

    lines = ps_text.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Line segment: x1 y1 x2 y2 L
        m = re.match(
            r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+L$',
            line
        )
        if m:
            ops.append({
                'op': 'line',
                'x1': float(m.group(1)), 'y1': float(m.group(2)),
                'x2': float(m.group(3)), 'y2': float(m.group(4)),
                'linewidth': current_linewidth,
            })
            i += 1
            continue

        # setlinewidth
        m = re.match(r'(-?\d+\.?\d*)\s+setlinewidth$', line)
        if m:
            current_linewidth = float(m.group(1))
            i += 1
            continue

        # Centred text: moveto / (text) / size Center / (text) / size Print
        # Pattern spans 5 lines
        if line.endswith('moveto') and i + 4 < len(lines):
            line1 = lines[i + 1].strip()
            line2 = lines[i + 2].strip()
            line3 = lines[i + 3].strip()
            line4 = lines[i + 4].strip()
            m_pos = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+moveto$', line)
            m_text1 = re.match(r'\((.+)\)$', line1)
            m_center = re.match(r'(-?\d+\.?\d*)\s+Center$', line2)
            if m_pos and m_text1 and m_center:
                # This is a centred text block
                m_text2 = re.match(r'\((.+)\)$', line3)
                m_print = re.match(r'(-?\d+\.?\d*)\s+Print$', line4)
                if m_text2 and m_print:
                    ops.append({
                        'op': 'text',
                        'x': float(m_pos.group(1)),
                        'y': float(m_pos.group(2)),
                        'text': m_text2.group(1),
                        'size': float(m_print.group(1)),
                        'centered': True,
                    })
                    i += 5
                    continue

        # Regular text: moveto / (text) / size Print
        # Pattern spans 3 lines
        if line.endswith('moveto') and i + 2 < len(lines):
            m_pos = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+moveto$', line)
            line1 = lines[i + 1].strip()
            line2 = lines[i + 2].strip()
            m_text = re.match(r'\((.+)\)$', line1)
            m_print = re.match(r'(-?\d+\.?\d*)\s+Print$', line2)
            if m_pos and m_text and m_print:
                ops.append({
                    'op': 'text',
                    'x': float(m_pos.group(1)),
                    'y': float(m_pos.group(2)),
                    'text': m_text.group(1),
                    'size': float(m_print.group(1)),
                    'centered': False,
                })
                i += 3
                continue

        # Regular text with Bprint (bold): moveto / (text) / size Bprint
        if line.endswith('moveto') and i + 2 < len(lines):
            m_pos = re.match(r'(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+moveto$', line)
            line1 = lines[i + 1].strip()
            line2 = lines[i + 2].strip()
            m_text = re.match(r'\((.+)\)$', line1)
            m_bprint = re.match(r'(-?\d+\.?\d*)\s+Bprint$', line2)
            if m_pos and m_text and m_bprint:
                ops.append({
                    'op': 'text',
                    'x': float(m_pos.group(1)),
                    'y': float(m_pos.group(2)),
                    'text': m_text.group(1),
                    'size': float(m_bprint.group(1)),
                    'centered': False,
                    'bold': True,
                })
                i += 3
                continue

        i += 1

    return ops


def _ops_to_svg(ops, bbox):
    """Convert parsed operations to SVG string.

    bbox: (x_min, y_min, x_max, y_max) in PS coordinates.
    PS has y increasing upward; SVG has y increasing downward, so we flip.
    """
    x_min, y_min, x_max, y_max = bbox
    width = x_max - x_min
    height = y_max - y_min

    def flip_y(y):
        return y_max - (y - y_min)

    svg_parts = []
    svg_parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="{x_min} 0 {width} {height}" '
        f'width="{width}" height="{height}" '
        f'style="background: white; max-width: 100%;">'
    )

    # Style definitions
    svg_parts.append('<style>')
    svg_parts.append('  text { font-family: "Times New Roman", Times, serif; fill: #333; }')
    svg_parts.append('  text.bold { font-weight: bold; }')
    svg_parts.append('  line { stroke-linecap: round; }')
    svg_parts.append('</style>')

    # Group lines by linewidth for efficiency
    current_group_width = None
    current_group_color = None

    for op in ops:
        if op['op'] == 'line':
            lw = op['linewidth']
            style = CONTOUR_STYLES.get(lw, {'stroke': '#333333', 'width': lw})
            stroke = style['stroke']
            sw = style['width']

            if stroke != current_group_color or sw != current_group_width:
                if current_group_width is not None:
                    svg_parts.append('</g>')
                svg_parts.append(f'<g stroke="{stroke}" stroke-width="{sw}">')
                current_group_width = sw
                current_group_color = stroke

            x1, y1 = op['x1'], flip_y(op['y1'])
            x2, y2 = op['x2'], flip_y(op['y2'])
            # Skip zero-length lines (degenerate segments)
            if x1 == x2 and y1 == y2:
                continue
            svg_parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}"/>')

        elif op['op'] == 'text':
            if current_group_width is not None:
                svg_parts.append('</g>')
                current_group_width = None
                current_group_color = None

            x = op['x']
            y = flip_y(op['y'])
            size = op['size']
            text = xml_escape(op['text'])
            anchor = 'middle' if op.get('centered') else 'start'
            bold_class = ' class="bold"' if op.get('bold') else ''

            svg_parts.append(
                f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}"{bold_class} '
                f'text-anchor="{anchor}">{text}</text>'
            )

    if current_group_width is not None:
        svg_parts.append('</g>')

    svg_parts.append('</svg>')
    return '\n'.join(svg_parts)


def ps_to_svg(ps_text):
    """Convert molrep PostScript text to an SVG string.

    Args:
        ps_text: Contents of a molrep_rf.ps file.

    Returns:
        SVG string suitable for embedding in HTML/XML reports.
    """
    # Extract bounding box
    m = re.search(r'%%BoundingBox:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', ps_text)
    if m:
        bbox = (float(m.group(1)), float(m.group(2)),
                float(m.group(3)), float(m.group(4)))
    else:
        bbox = (0, 0, 600, 800)

    ops = _parse_ps(ps_text)
    return _ops_to_svg(ops, bbox)


def convert_file(ps_path, svg_path=None):
    """Convert a molrep PS file to SVG.

    Args:
        ps_path: Path to input .ps file.
        svg_path: Path to output .svg file. If None, returns SVG string.

    Returns:
        SVG string if svg_path is None, otherwise writes file and returns path.
    """
    with open(ps_path, 'r') as f:
        ps_text = f.read()

    svg = ps_to_svg(ps_text)

    if svg_path is not None:
        with open(svg_path, 'w') as f:
            f.write(svg)
        return svg_path

    return svg
