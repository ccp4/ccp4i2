import os


def which(program):
    """Checks if `program` exists and finds its location. Analogy of the
    `which` GNU/Linux command.

    Args:
        program (str): Name of an executable

    Returns:
        str: Path of an executable location
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program) or is_exe(program + ".exe"):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
            exe_file = os.path.join(path, program + ".exe")
            if is_exe(exe_file):
                return exe_file
    return None