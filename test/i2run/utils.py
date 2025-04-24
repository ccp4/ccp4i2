from contextlib import contextmanager
from multiprocessing import Process
from os.path import basename, join
from pathlib import Path
from random import choice
from shutil import rmtree
from string import ascii_letters, digits
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
from urllib.request import urlopen
from xml.etree import ElementTree as ET

import gemmi

from ccp4i2.core import CCP4I2Runner
from ccp4i2.core.CCP4Utils import getCCP4I2Dir


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a pathlib.Path object to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    urlName = unquote(basename(urlparse(url).path))
    with urlopen(url, timeout=30) as response:
        name = response.headers.get_filename() or urlName
        name = name.strip().replace(" ", "_")
        name = "".join(c for c in name if c.isalnum() or c in "-_.")
        with NamedTemporaryFile(suffix=f"_{name}", delete=False) as temp:
            while chunk := response.read(1_000_000):
                temp.write(chunk)
        path = Path(temp.name).resolve()
        try:
            yield str(path)
        finally:
            path.unlink(missing_ok=True)


@contextmanager
def i2run(args: list[str]):
    """
    Run a task through CCP4I2Runner with the given arguments,
    check the diagnostic.xml file does not contain any error reports
    and yield a Path object to the job directory (CCP4_JOBS/job_1).
    Use in a with statement and the project will be cleaned up
    as long as an error is not raised.
    """
    chars = ascii_letters + digits
    tmp_name = "tmp_" + "".join(choice(chars) for _ in range(10))
    args = ["i2run"] + args
    args += ["--projectName", tmp_name]
    args += ["--projectPath", tmp_name]
    args += ["--dbFile", f"{tmp_name}.sqlite"]
    # Must be a separate process because calling CCP4i2Runner.main
    # multiple times (even sequentially) is not supported
    process = Process(target=CCP4I2Runner.main, args=(args,))
    process.start()
    process.join()
    directory = Path(tmp_name, "CCP4_JOBS", "job_1")
    xml_path = directory / "diagnostic.xml"
    errors = ET.parse(xml_path).findall(".//errorReport")
    assert len(errors) == 0, "Error reports found in diagnostic.xml"
    # Below code not inside a try/finally block
    # So that the project is only removed if an error is not raised
    yield directory
    rmtree(tmp_name, ignore_errors=True)
    for extension in ("sqlite", "sqlite-shm", "sqlite-wal"):
        Path(f"{tmp_name}.{extension}").unlink(missing_ok=True)


def demoData(*paths):
    return join(getCCP4I2Dir(), "demo_data", *paths)


def hasLongLigandName(path):
    "Does the structure contains a residue with a name longer than 3 characters?"
    structure = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Mmcif)
    for model in structure:
        for chain in model:
            for residue in chain:
                if len(residue.name) > 3:
                    return True
    return False
