import pathlib
import subprocess
from ..helpers.terminal import open_terminal_in_directory


def preview_job(viewer: str = None, file_path: pathlib.Path = None):
    """
    Preview a job using the appropriate viewer.

    Args:
        viewer (str, optional): The name of the viewer to use. Defaults to None.
        file_path (pathlib.Path, optional): The path to the job directory to be previewed. Defaults to None.

    Returns:
        None

    Raises:
        ValueError: If the viewer is not supported.
    """

    if viewer == "terminal":
        open_terminal_in_directory(str(pathlib.Path(file_path)))
        return {"status": "Success"}

    raise ValueError(f"Viewer {viewer} is not supported.")
