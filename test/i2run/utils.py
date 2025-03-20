from contextlib import contextmanager
from email.message import EmailMessage
from multiprocessing import Process
from os.path import basename
from pathlib import Path
from random import choice
from shutil import rmtree
from string import ascii_letters, digits
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
from xml.etree import ElementTree as ET

from requests import get, Response

from ccp4i2.core import CCP4I2Runner


def valid_filename_from_response(response: Response):
    """
    Extracts a valid filename from the response headers or URL.
    Ensures that the filename is safe to use by stripping whitespace,
    replacing spaces with underscores, and removing any characters
    that are not alphanumeric, dash, underscore or dot.
    """
    message = EmailMessage()
    for header, value in response.headers.items():
        message[header] = value
    url_name = unquote(basename(urlparse(response.url).path))
    name = message.get_filename() or url_name
    name = name.strip().replace(" ", "_")
    name = "".join(c for c in name if c.isalnum() or c in "-_.")
    return name


@contextmanager
def download(url: str):
    """
    Downloads a file from the given URL and saves it to a temporary file.
    Yields a pathlib.Path object to the temporary file.
    Use in a with statement to ensure the file is deleted afterwards.
    """
    response = get(url, allow_redirects=True, stream=True, timeout=30)
    response.raise_for_status()
    suffix = f"_{valid_filename_from_response(response)}"
    with NamedTemporaryFile(suffix=suffix, delete=False) as temp:
        for chunk in response.iter_content(chunk_size=1_000_000):
            temp.write(chunk)
    path = Path(temp.name).resolve()
    try:
        yield path
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
    assert len(list(ET.parse(xml_path).iter("errorReport"))) == 0
    # Below code not inside a try/finally block
    # So that the project is only removed if an error is not raised
    yield directory
    rmtree(tmp_name, ignore_errors=True)
    for extension in ("sqlite", "sqlite-shm", "sqlite-wal"):
        Path(f"{tmp_name}.{extension}").unlink(missing_ok=True)
