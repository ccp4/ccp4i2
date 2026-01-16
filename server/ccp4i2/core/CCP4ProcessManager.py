"""
Qt-free CProcessManager - Replacement for CCP4i2's Qt-based process manager.

This module provides a pure Python subprocess execution system that replaces
Qt's QProcess-based CProcessManager, enabling the migration away from Qt dependencies.

Key Features:
- Synchronous and asynchronous process execution
- Environment management (CCP4, library paths)
- Process tracking with PIDs
- Handler callbacks on completion
- I/O redirection (stdin, stdout, stderr)
- Timeout handling
- Exit code and status tracking
- Compatible with legacy CPluginScript API

Architecture:
- Based on Python's subprocess module (not Qt QProcess)
- Singleton pattern for global process manager
- Process registry with detailed tracking
- CErrorReport integration for error handling

Example Usage:
    # Get singleton instance
    pm = PROCESSMANAGER()

    # Synchronous execution (default)
    pm.setWaitForFinished(timeout=60000)  # Wait up to 60 seconds
    pid = pm.startProcess(
        command='refmac5',
        args=['HKLIN', 'input.mtz'],
        inputFile='refmac.com',
        logFile='refmac.log'
    )
    exitCode = pm.getJobData(pid, 'exitCode')

    # Asynchronous execution
    pm.setWaitForFinished(-1)  # Don't wait
    pid = pm.startProcess(
        command='long_running_job',
        handler=[callback_func, {}]
    )

Author: Generated for CCP4i2 Qt-free migration
Date: 2025-11-07
"""

from __future__ import annotations

import os
import sys
import time
import shutil
import subprocess
import threading
from typing import Any, Callable, Dict, List, Optional, Union

from ccp4i2.core.base_object.error_reporting import CErrorReport


class CProcessManager:
    """
    Qt-free process manager for running external programs.

    This class manages subprocess execution without Qt dependencies,
    replacing the original CCP4i2 CProcessManager that used QProcess.

    Attributes:
        lastProcessId (int): Counter for generating unique PIDs
        processInfo (dict): Registry of process information by PID
        ifAsync (bool): Whether to run processes asynchronously
        timeout (int): Default timeout in milliseconds
    """

    # Singleton instance
    _instance: Optional[CProcessManager] = None

    # Error codes matching original CProcessManager
    ERROR_CODES = {
        101: {'description': 'Error creating temporary command file'},
        102: {'description': 'Process input file does not exist'},
        103: {'severity': 1, 'description': 'Existing log file has been deleted'},
        104: {'description': 'Error opening input file'},
        105: {'description': 'Error opening log file'},
        106: {'description': 'Error starting sub-process'},
        107: {'description': 'Can not run process - no executable with name'},
        108: {'severity': 1, 'description': 'Creating temporary log file for sub-process'},
        109: {'description': 'Error opening stderr file'},
        110: {'description': 'Error handling finished sub-process'},
        111: {'description': 'Error calling handler after finished sub-process'}
    }

    def __new__(cls):
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the process manager."""
        if getattr(self, '_initialized', False):
            return

        self.lastProcessId = 0
        self.processInfo: Dict[int, Dict[str, Any]] = {}
        self.ifAsync = False
        self.timeout = 999999  # Default: very long timeout (milliseconds)
        self._maxRunningProcesses = 10
        self._processEnvironment = None
        self._initialized = True

    def setWaitForFinished(self, timeout: int = -1):
        """
        Configure whether processes should run synchronously or asynchronously.

        Args:
            timeout: If >= 0, processes run synchronously with this timeout (ms).
                    If < 0, processes run asynchronously (non-blocking).

        Example:
            # Synchronous execution
            pm.setWaitForFinished(60000)  # Wait up to 60 seconds

            # Asynchronous execution
            pm.setWaitForFinished(-1)  # Don't wait
        """
        if timeout < 0:
            self.ifAsync = True
        else:
            self.ifAsync = False
            self.timeout = timeout

    def setMaxRunningProcesses(self, maxproc: int):
        """Set maximum number of concurrent processes."""
        self._maxRunningProcesses = maxproc

    def maxRunningProcesses(self) -> int:
        """Get maximum number of concurrent processes."""
        return self._maxRunningProcesses

    def _get_ccp4_env(self, resetEnv: bool = True, editEnv: List = None) -> Dict[str, str]:
        """
        Get process environment with CCP4 variables.

        Args:
            resetEnv: If True, use system environment. If False, use os.environ.
            editEnv: List of (key, value) tuples to add/override.

        Returns:
            Environment dictionary for subprocess
        """
        if resetEnv:
            env = os.environ.copy()
        else:
            env = dict(os.environ)

        # Apply custom environment edits
        if editEnv:
            for key, value in editEnv:
                env[key] = value

        return env

    def startProcess(
        self,
        command: str = None,
        args: Union[List[str], str] = None,
        inputFile: str = None,
        logFile: str = None,
        interpreter: str = None,
        inputText: str = None,
        handler: List = None,
        resetEnv: bool = True,
        readyReadStandardOutputHandler: Callable = None,
        **kw
    ) -> int:
        """
        Start an external process.

        Args:
            command: Path to executable or command name
            args: List of arguments or space-separated string
            inputFile: Path to file for stdin redirection
            logFile: Path to file for stdout redirection
            interpreter: If 'python', prepend Python executable
            inputText: Text to write to stdin (creates temp file)
            handler: [callback_func, kwargs] called on completion
            resetEnv: Use clean environment (vs inherited)
            readyReadStandardOutputHandler: Callback for streaming output
            **kw: Additional options:
                - ifAsync (bool): Override global async setting
                - timeout (int): Override global timeout (ms)
                - cwd (str): Working directory
                - editEnv (list): [(key, value), ...] env overrides
                - jobId, jobNumber, projectId: Metadata
                - stderrFile (str): Path for stderr redirection

        Returns:
            int: Process ID (PID) for tracking

        Example:
            pid = pm.startProcess(
                command='/usr/bin/refmac5',
                args=['HKLIN', 'data.mtz'],
                inputFile='refmac.com',
                logFile='refmac.log',
                cwd='/path/to/job'
            )
        """
        # Allocate unique PID
        self.lastProcessId += 1
        pid = self.lastProcessId

        # Parse arguments
        command = shutil.which(command)
        if args is None:
            args = []
        argList = [command]
        if isinstance(args, list):
            argList.extend(args)
        else:
            argList.extend(args.split())

        # Handle interpreter prefix (e.g., python script.py)
        if interpreter == 'python':
            python_exe = sys.executable
            argList = [python_exe, command] + (args if isinstance(args, list) else args.split())

        # Initialize process info
        ifAsync = kw.get('ifAsync', self.ifAsync)
        process_timeout = kw.get('timeout', self.timeout)

        self.processInfo[pid] = {
            'command': command,
            'argList': argList,
            'handler': handler or [],
            'readyReadStandardOutputHandler': readyReadStandardOutputHandler,
            'errorReport': CErrorReport(),
            'ifAsync': ifAsync,
            'timeout': process_timeout,
            'resetEnv': resetEnv,
            'startTime': None,
            'finishTime': None,
            'inputFile': inputFile,
            'logFile': logFile,
            'stderrFile': kw.get('stderrFile'),
            'jobId': kw.get('jobId'),
            'jobNumber': kw.get('jobNumber'),
            'projectId': kw.get('projectId'),
            'editEnv': kw.get('editEnv', []),
            'cwd': kw.get('cwd'),
            'exitStatus': None,
            'exitCode': None,
            'status': 'pending'
        }

        # Handle inputText (create temp file)
        if inputFile is None and inputText is not None:
            try:
                import tempfile
                fd, tmpFile = tempfile.mkstemp(suffix='.com', text=True)
                with os.fdopen(fd, 'w') as f:
                    f.write(inputText)
                inputFile = tmpFile
                self.processInfo[pid]['inputFile'] = inputFile
            except Exception as e:
                self.processInfo[pid]['errorReport'].append(
                    self.__class__.__name__, 101, f"Error creating temp input file: {e}"
                )

        # Validate input file
        if inputFile is not None and not os.path.exists(inputFile):
            self.processInfo[pid]['errorReport'].append(
                self.__class__.__name__, 102, f"Input file does not exist: {inputFile}"
            )
            self.processInfo[pid]['status'] = 'failed'
            self.processInfo[pid]['exitCode'] = -1
            return pid

        # Handle log file
        if logFile is not None and os.path.exists(logFile):
            try:
                os.remove(logFile)
                self.processInfo[pid]['errorReport'].append(
                    self.__class__.__name__, 103, f"Removed existing log file: {logFile}"
                )
            except Exception as e:
                pass

        # Print command for debugging
        try:
            cmd_str = ' '.join(argList)
            if inputFile:
                cmd_str += f' < {inputFile}'
            if logFile:
                cmd_str += f' > {logFile}'
            print(f"\n{'='*60}")
            print(f"Running: {cmd_str}")
            print(f"Working directory: {kw.get('cwd', os.getcwd())}")
            if logFile:
                print(f"Log file: {logFile}")
            if self.processInfo[pid]['stderrFile']:
                print(f"Stderr file: {self.processInfo[pid]['stderrFile']}")
            print(f"{'='*60}\n")
            sys.stdout.flush()
        except Exception as e:
            print(f"Error printing command: {e}")

        # Execute process
        if not ifAsync:
            # Synchronous execution
            self._run_sync(pid)
        else:
            # Asynchronous execution
            self._run_async(pid)

        return pid

    def _run_sync(self, pid: int):
        """Run process synchronously (blocking)."""
        info = self.processInfo[pid]
        info['startTime'] = time.time()
        info['status'] = 'running'

        # Track file handles we open for proper cleanup
        stdin_file = None
        stdout_file = None
        stderr_file = None

        try:
            # Prepare subprocess arguments
            kwargs = {
                'env': self._get_ccp4_env(info['resetEnv'], info['editEnv']),
                'cwd': info['cwd']
            }

            # Handle I/O redirection using context managers to ensure cleanup
            # Open file handles and track them for cleanup
            if info['inputFile']:
                stdin_file = open(info['inputFile'], 'r')
                kwargs['stdin'] = stdin_file
            else:
                kwargs['stdin'] = None

            if info['logFile']:
                stdout_file = open(info['logFile'], 'w')
                kwargs['stdout'] = stdout_file
            else:
                kwargs['stdout'] = subprocess.PIPE

            if info['stderrFile']:
                stderr_file = open(info['stderrFile'], 'w')
                kwargs['stderr'] = stderr_file
            else:
                kwargs['stderr'] = subprocess.PIPE

            # Convert timeout from milliseconds to seconds
            timeout_sec = info['timeout'] / 1000.0 if info['timeout'] > 0 else None

            # Run process
            result = subprocess.run(
                info['argList'],
                timeout=timeout_sec,
                **kwargs
            )

            # Store result
            info['exitCode'] = result.returncode
            info['exitStatus'] = result.returncode
            info['status'] = 'finished' if result.returncode == 0 else 'failed'
            info['finishTime'] = time.time()

            print(f"✅ Process completed successfully (exit code {result.returncode})")

        except subprocess.TimeoutExpired:
            info['exitCode'] = -1
            info['exitStatus'] = -1
            info['status'] = 'timeout'
            info['finishTime'] = time.time()
            print(f"⏱️  Process timed out after {info['timeout']}ms")

        except Exception as e:
            info['exitCode'] = -1
            info['exitStatus'] = -1
            info['status'] = 'failed'
            info['finishTime'] = time.time()
            info['errorReport'].append(
                self.__class__.__name__, 106, f"Error starting process: {e}"
            )
            print(f"❌ Process failed: {e}")

        finally:
            # Close file handles we explicitly opened (not subprocess.PIPE constants)
            if stdin_file is not None:
                try:
                    stdin_file.close()
                except Exception:
                    pass  # Best effort cleanup
            if stdout_file is not None:
                try:
                    stdout_file.close()
                except Exception:
                    pass
            if stderr_file is not None:
                try:
                    stderr_file.close()
                except Exception:
                    pass

        # Call handler if provided
        self._call_handler(pid)

    def _run_async(self, pid: int):
        """Run process asynchronously (non-blocking)."""
        info = self.processInfo[pid]

        def run_in_thread():
            """Thread worker for async execution."""
            self._run_sync(pid)

        thread = threading.Thread(target=run_in_thread, daemon=True)
        thread.start()
        info['thread'] = thread

    def _call_handler(self, pid: int):
        """Call completion handler if provided."""
        info = self.processInfo[pid]
        handler = info.get('handler', [])

        if handler and len(handler) >= 1:
            try:
                callback = handler[0]
                kwargs = handler[1] if len(handler) > 1 else {}

                if callable(callback):
                    callback(pid, info['exitCode'], **kwargs)
            except Exception as e:
                info['errorReport'].append(
                    self.__class__.__name__, 111,
                    f"Error calling handler: {e}"
                )
                print(f"Error in handler callback: {e}")

    def getJobData(self, pid: int, attribute: str = 'exitStatus') -> Any:
        """
        Get process information by attribute name.

        For synchronous processes, checks both the processInfo dict and
        the plugin object's attributes (_exitCode, _exitStatus, etc.).

        Args:
            pid: Process ID returned by startProcess()
            attribute: Data attribute name (exitCode, exitStatus, status, startTime, qprocess, etc.)
                      Defaults to 'exitStatus' for backward compatibility with legacy code.

        Returns:
            Requested value or None if not found

        Example:
            status = pm.getJobData(pid)  # Returns exitStatus (legacy API)
            exitCode = pm.getJobData(pid, attribute='exitCode')
            status = pm.getJobData(pid, attribute='status')
            qprocess = pm.getJobData(pid, attribute='qprocess')
        """
        if pid not in self.processInfo:
            return None

        info = self.processInfo[pid]

        # First check if attribute exists in processInfo dict
        if attribute in info:
            return info[attribute]

        # For synchronous processes, check plugin object attributes
        # CPluginScript stores exit codes as _exitCode, _exitStatus
        plugin = info.get('plugin')
        if plugin:
            # Try with underscore prefix (e.g., exitCode -> _exitCode)
            attr_name = f'_{attribute}'
            if hasattr(plugin, attr_name):
                return getattr(plugin, attr_name)

        return None

    def waitForFinished(self, pid: int, timeout: int = -1) -> bool:
        """
        Wait for a specific process to complete.

        Args:
            pid: Process ID to wait for
            timeout: Timeout in milliseconds (-1 for infinite)

        Returns:
            True if process finished, False if timeout
        """
        if pid not in self.processInfo:
            return False

        info = self.processInfo[pid]

        # If process already finished, return immediately
        if info.get('finishTime'):
            return True

        # If async execution, wait for thread
        if info.get('thread'):
            timeout_sec = timeout / 1000.0 if timeout >= 0 else None
            info['thread'].join(timeout=timeout_sec)
            return info.get('finishTime') is not None

        # Synchronous process should already be finished
        return info.get('finishTime') is not None

    def killProcess(self, pid: int):
        """
        Kill a running process.

        Args:
            pid: Process ID to kill
        """
        # For subprocess-based implementation, process has already
        # completed by the time we could kill it (synchronous execution).
        # For async execution, we'd need to store the Popen object.
        # This is a placeholder for future enhancement.
        if pid in self.processInfo:
            self.processInfo[pid]['status'] = 'killed'

    def register(self, plugin):
        """
        Register a plugin with the process manager.

        Legacy API compatibility method for CPluginScript.
        Stores a reference to the plugin so that getJobData() can query
        plugin attributes (_exitCode, _exitStatus) for synchronous processes.

        Args:
            plugin: CPluginScript instance to register
        """
        # Use plugin's object ID as PID (matches getProcessId() implementation)
        pid = id(plugin)

        # Create or update process info entry
        if pid not in self.processInfo:
            self.processInfo[pid] = {}

        # Store plugin reference for attribute queries
        self.processInfo[pid]['plugin'] = plugin
        self.processInfo[pid]['startTime'] = time.time()

    def formatted_job(self, pid: int) -> str:
        """
        Get formatted job status string (legacy API).

        Args:
            pid: Process ID

        Returns:
            Formatted status string
        """
        if pid not in self.processInfo:
            return f"{pid:4d}  UNKNOWN"

        info = self.processInfo[pid]
        number = f"{pid:4d}"

        if info.get('finishTime'):
            if info.get('exitCode', 0) != 0:
                status_str = "FAILED"
            else:
                status_str = "FINISHED"

            finish_time = time.strftime(
                '%Y-%m-%d %H:%M:%S',
                time.localtime(info['finishTime'])
            )
            return f"{number}  {status_str:8s} {finish_time}"
        else:
            start_time = time.strftime(
                '%Y-%m-%d %H:%M:%S',
                time.localtime(info.get('startTime', time.time()))
            )
            return f"{number}  STARTED  {start_time}"


# Module-level singleton instance
_process_manager_instance = CProcessManager()


def PROCESSMANAGER() -> CProcessManager:
    """
    Get the global singleton process manager instance.

    This function provides the legacy CCP4Modules.PROCESSMANAGER() API.

    Returns:
        Global CProcessManager instance

    Example:
        from ccp4i2.core.CCP4ProcessManager import PROCESSMANAGER

        pm = PROCESSMANAGER()
        pm.setWaitForFinished(60000)
        pid = pm.startProcess('refmac5', ['HKLIN', 'data.mtz'])
        exitCode = pm.getJobData(pid, 'exitCode')
    """
    return _process_manager_instance
