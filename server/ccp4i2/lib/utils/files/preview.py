import pathlib
import platform
import subprocess
from ..helpers.terminal import open_terminal_in_directory


def get_postscript_viewer():
    """
    Get the appropriate PostScript/PDF viewer for the current platform.

    Returns:
        list: The command list for opening PostScript/PDF files, or None if unavailable.
    """
    system = platform.system()

    if system == "Darwin":  # macOS
        # Use "open -a Preview" which can handle both PS and PDF files
        return ["open", "-a", "Preview"]
    if system == "Linux":
        # Try common Linux PostScript/PDF viewers in order of preference
        viewers = [
            ["evince"],           # GNOME document viewer (supports PS)
            ["okular"],           # KDE document viewer
            ["gv"],               # Ghostview - dedicated PostScript viewer
            ["xdg-open"],         # System default handler
        ]
        for viewer in viewers:
            try:
                # Check if the viewer is available
                subprocess.run(
                    ["which", viewer[0]],
                    capture_output=True,
                    check=True
                )
                return viewer
            except subprocess.CalledProcessError:
                continue
        return ["xdg-open"]  # Fallback to system default
    if system == "Windows":
        # Windows - use start command to open with default handler
        return ["start", ""]
    return None


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
    if viewer == "postscript":
        ps_viewer = get_postscript_viewer()
        if ps_viewer is None:
            raise ValueError("No PostScript viewer found on this system.")

        subprocess.Popen(
            ps_viewer + [str(file_path)],
            start_new_session=True,
        )
        return {"status": "Success"}

    raise ValueError(f"Viewer {viewer} is not supported.")
