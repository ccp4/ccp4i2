from os import environ
from pathlib import Path


VERSION = "9.0"
if "CLIB" in environ:
    path = Path(environ["CLIB"], "ccp4", "MAJOR_MINOR")
    if path.exists():
        with path.open(encoding="utf-8") as text:
            VERSION = text.read().strip()


def get_revno():
    numeric = ''.join(VERSION.split('.'))
    return int(numeric) if numeric.isdigit() else 0


def get_ccp4_str():
    return f"CCP4-{VERSION}"
