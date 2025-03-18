from os import environ
from pathlib import Path


CCP4_VERSION = "9.0"
if "CLIB" in environ:
    path = Path(environ["CLIB"], "ccp4", "MAJOR_MINOR")
    if path.exists():
        with path.open(encoding="utf-8") as text:
            CCP4_VERSION = text.read().strip()
