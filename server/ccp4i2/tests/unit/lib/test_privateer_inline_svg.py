"""Privateer's report embeds SVG files (its glycan drawings and the
Cremer-Pople sphere) inside report elements. Each file starts with an XML
declaration, which is legal only at the start of a document, so the whole
report failed to parse: "XML or text declaration not at start of entity".
"""
import xml.etree.ElementTree as ET

from ccp4i2.report.svg import fit_svg, inline_svg

SVG = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<!-- Generator: Privateer -->
<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><rect width="1" height="1"/></svg>
"""


def test_embedded_svg_parses():
    with_declaration = '<div style="float:left;">' + SVG + "</div>"
    try:
        ET.fromstring(with_declaration)
        raise AssertionError("the declaration was expected to break parsing")
    except ET.ParseError:
        pass
    div = ET.fromstring('<div style="float:left;">' + inline_svg(SVG) + "</div>")
    assert div.find("{http://www.w3.org/2000/svg}svg") is not None


def test_svg_without_declaration_unchanged():
    plain = '<svg xmlns="http://www.w3.org/2000/svg"/>'
    assert inline_svg(plain) == plain


def test_fixed_size_svg_scales_to_its_box():
    fitted = fit_svg('<svg xmlns="http://www.w3.org/2000/svg" width="430" height="430"><g/></svg>', 250)
    svg = ET.fromstring(fitted)
    assert (svg.get("width"), svg.get("height"), svg.get("viewBox")) == ("250", "250", "0 0 430 430")


def test_svg_with_viewbox_left_alone():
    already = '<svg viewBox="0 0 10 10" width="10" height="10"/>'
    assert fit_svg(already, 250) == already
