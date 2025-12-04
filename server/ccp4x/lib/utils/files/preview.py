import pathlib
import subprocess
from ..helpers.terminal import open_terminal_in_directory


def preview_file(viewer: str = None, file_path: pathlib.Path = None):
    """
    Preview a file using the appropriate viewer.

    Args:
        viewer (str, optional): The name of the viewer to use. Defaults to None.
        file_path (pathlib.Path, optional): The path to the file to be previewed. Defaults to None.

    Returns:
        None

    Raises:
        ValueError: If the viewer is not supported.
    """
    if viewer == "coot":
        # Here should dial into mechanism to use the user-specified path to coot as appropriate.
        # Also suspect will need diverse environment variables to be set for this to work.
        subprocess.Popen(
            ["coot", "--no-guano", str(file_path)],
            start_new_session=True,
        )
        return {"status": "Success"}
    if viewer == "ccp4mg":
        # Here should dial into mechanism to use the user-specified path to coot as appropriate.
        # Also suspect will need diverse environment variables to be set for this to work.
        subprocess.Popen(
            ["ccp4mg", str(file_path)],
            start_new_session=True,
        )
        return {"status": "Success"}
    if viewer == "viewhkl":
        # Here should dial into mechanism to use the user-specified path to coot as appropriate.
        # Also suspect will need diverse environment variables to be set for this to work.
        subprocess.Popen(
            ["viewhkl", str(file_path)],
            start_new_session=True,
        )
        return {"status": "Success"}
    if viewer == "terminal":
        open_terminal_in_directory(str(pathlib.Path(file_path).parent))
        return {"status": "Success"}

    raise ValueError(f"Viewer {viewer} is not supported.")
