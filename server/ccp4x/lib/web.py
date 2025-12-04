from contextlib import contextmanager
from email.message import EmailMessage
from os.path import basename
from pathlib import Path
from tempfile import NamedTemporaryFile
from urllib.parse import urlparse, unquote
from django.utils.text import get_valid_filename
from requests import Response
from requests_cache import CachedSession
from ..config.settings import USER_DIR


_SESSION = CachedSession(USER_DIR / "http-cache")


@contextmanager
def download(url: str):
    response = _SESSION.get(url, allow_redirects=True, stream=True, timeout=30)
    response.raise_for_status()
    suffix = f"_{filename(response)}"
    with NamedTemporaryFile(suffix=suffix) as tmp_file:
        with open(tmp_file.name, "wb") as writer:
            for chunk in response.iter_content(chunk_size=1_000_000):
                writer.write(chunk)
        yield Path(tmp_file.name)


def filename(response: Response) -> str:
    url_name = unquote(basename(urlparse(response.url).path))
    message = EmailMessage()
    for header, value in response.headers.items():
        message[header] = value
    return get_valid_filename(message.get_filename() or url_name)
