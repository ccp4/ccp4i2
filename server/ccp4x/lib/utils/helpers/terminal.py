import os
import subprocess
import sys


def open_terminal_in_directory(directory):
    # Make sure the directory exists
    if not os.path.isdir(directory):
        print(f"Error: The directory '{directory}' does not exist.")
        return

    # Normalize the path for different operating systems
    directory = os.path.abspath(directory)

    if sys.platform == "win32":  # For Windows
        subprocess.run(["start", "cmd", "/K", f"cd /d {directory}"], shell=True)

    elif sys.platform == "darwin":  # For macOS
        # AppleScript to tell Terminal to open and bring to the front
        applescript = f"""
        tell application "Terminal"
            activate
            do script "cd {directory} && clear"
        end tell
        """
        subprocess.run(["osascript", "-e", applescript])

    elif sys.platform == "linux" or sys.platform == "linux2":  # For Linux
        subprocess.run(["gnome-terminal", "--working-directory", directory])

    else:
        print("Unsupported OS")
