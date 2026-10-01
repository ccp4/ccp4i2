"""Embedding SVG files in a report."""
import re


def inline_svg(text):
    """An SVG file's text, fit to go inside a report element: without its XML
    declaration and DOCTYPE, which are legal only at the start of a document.
    Embedded as they came, they made the whole report fail to parse ("XML or
    text declaration not at start of entity")."""
    text = re.sub(r"^\s*<\?xml[^>]*\?>", "", text)
    return re.sub(r"<!DOCTYPE[^>\[]*(\[[^\]]*\])?\s*>", "", text).lstrip()


def fit_svg(text, size):
    """Scale an SVG drawn at a fixed size to fit a size-by-size box: give it
    a viewBox from its own width and height, then the box's dimensions. (The
    Cremer-Pople sphere is drawn at 430 px and overflowed its 250 px box
    across the plot beside it.)"""
    m = re.search(r'<svg\b[^>]*?\swidth="([\d.]+)"[^>]*?\sheight="([\d.]+)"', text)
    if not m or "viewBox" in text[:m.end()]:
        return text
    w, h = m.group(1), m.group(2)
    tag = m.group(0)
    fitted = (tag.replace(f'width="{w}"', f'width="{size}"', 1)
                 .replace(f'height="{h}"', f'height="{size}" viewBox="0 0 {w} {h}"', 1))
    return text.replace(tag, fitted, 1)
