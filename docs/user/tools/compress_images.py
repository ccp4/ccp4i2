"""Shrink the help's PNG figures in place: 256 colours, at most 2000 px wide.

    python3 compress_images.py [paths...]    (default: all of docs/user/source)

Screenshots are flat colour and text, which a 256-colour palette holds almost
exactly; the 2000 px cap is still twice the width the page shows. Across the
imported Qt set this took 45 MB to 12 MB. JPEG does worse on both counts: 26 MB,
and it blurs text. Needs Pillow. Files it cannot shrink are left alone.
"""
import io
import sys
from pathlib import Path

from PIL import Image

MAX_WIDTH = 2000


def shrink(path: Path) -> int:
    """Rewrite one PNG if that makes it smaller; return the bytes saved."""
    before = path.stat().st_size
    im = Image.open(path)
    im.load()
    if im.mode == "P":
        return 0  # already a palette image: compressed once
    if im.width > MAX_WIDTH:
        im = im.resize((MAX_WIDTH, round(im.height * MAX_WIDTH / im.width)),
                       Image.LANCZOS)
    rgba = im.convert("RGBA")
    if rgba.getextrema()[3][0] == 255:
        # Opaque: median cut on RGB keeps pale backgrounds their colour.
        out = rgba.convert("RGB").quantize(256, method=0)
    else:
        out = rgba.quantize(256, method=2)
    buf = io.BytesIO()
    out.save(buf, "PNG", optimize=True)
    if buf.tell() >= before:
        return 0
    path.write_bytes(buf.getvalue())
    return before - buf.tell()


def main():
    root = Path(__file__).resolve().parent.parent / "source"
    paths = [Path(p) for p in sys.argv[1:]] or [root]
    files = [f for p in paths
             for f in ([p] if p.is_file() else sorted(p.rglob("*.png")))]
    saved = sum(shrink(f) for f in files)
    print(f"{len(files)} PNGs, {saved / 1e6:.1f} MB saved")


if __name__ == "__main__":
    main()
