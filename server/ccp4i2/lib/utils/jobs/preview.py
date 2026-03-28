# Copyright (C) 2025 University of York
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
