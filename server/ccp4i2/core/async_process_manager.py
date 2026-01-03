"""
Async Process Manager - Qt-free replacement for CCP4ProcessManager.

This module provides async subprocess execution using Python's asyncio library,
replacing Qt's QProcess-based system with a pure Python implementation.

Key features:
- Async subprocess execution with asyncio
- Signal-based completion notification
- Process monitoring and timeout handling
- Compatible with existing CPluginScript API
- No Qt dependencies
"""

import asyncio
import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field


logger = logging.getLogger(__name__)


@dataclass
class ProcessInfo:
    """Information about a running or completed process."""

    pid: int
    command: str
    args: List[str] = field(default_factory=list)
    inputFile: Optional[str] = None
    logFile: Optional[str] = None
    cwd: Optional[str] = None
    env: Optional[Dict[str, str]] = None
    handler: Optional[Any] = None
    timeout: Optional[int] = None

    # Runtime info
    startTime: float = field(default_factory=time.time)
    finishTime: Optional[float] = None
    status: str = "pending"  # pending, running, finished, failed, timeout
    exitCode: Optional[int] = None
    exitStatus: Optional[int] = None
    error: Optional[str] = None

    # Process object
    process: Optional[asyncio.subprocess.Process] = None

    # File handles for log files (direct I/O, no buffering)
    stdout_file: Optional[Any] = None
    stderr_file: Optional[Any] = None


class AsyncProcessManager:
    """
    Manages async subprocess execution using asyncio.

    This replaces Qt's QProcess-based CProcessManager with a pure Python
    implementation using asyncio.

    Example:
        pm = AsyncProcessManager()

        # Start a process
        pid = await pm.startProcess(
            command='refmac5',
            args=['HKLIN', 'input.mtz', 'HKLOUT', 'output.mtz'],
            inputFile='refmac.com',
            logFile='refmac.log',
            handler=[callback_func, {}]
        )

        # Check status
        status = pm.getJobData(pid, 'status')
    """

    # Singleton instance
    _instance: Optional["AsyncProcessManager"] = None
    _loop: Optional[asyncio.AbstractEventLoop] = None
    _loop_thread: Optional[Any] = None

    def __new__(cls):
        """Ensure singleton instance."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the process manager."""
        if self._initialized:
            return

        self.processes: Dict[int, ProcessInfo] = {}
        self._next_pid = 1000
        self._lock = asyncio.Lock() if self._has_running_loop() else None
        self._executor = ThreadPoolExecutor(max_workers=4)
        self._max_concurrent = 10  # Maximum concurrent processes
        self._semaphore = None  # Created when loop available
        self._initialized = True

        # Ensure event loop is running
        self._ensure_event_loop()

        logger.info("AsyncProcessManager initialized")

    @staticmethod
    def _has_running_loop() -> bool:
        """Check if there's a running event loop."""
        try:
            asyncio.get_running_loop()
            return True
        except RuntimeError:
            return False

    def _ensure_event_loop(self):
        """Ensure there's an event loop running."""
        if self._has_running_loop():
            # Already have a loop
            return

        # Create a new event loop in a background thread
        import threading

        def run_event_loop():
            """Run event loop in background thread."""
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            AsyncProcessManager._loop = loop

            # Create lock and semaphore in this loop
            async def init_async():
                self._lock = asyncio.Lock()
                self._semaphore = asyncio.Semaphore(self._max_concurrent)

            loop.run_until_complete(init_async())

            logger.info("Event loop started in background thread")
            loop.run_forever()

        thread = threading.Thread(target=run_event_loop, daemon=True, name="AsyncProcessManager")
        thread.start()
        AsyncProcessManager._loop_thread = thread

        # Wait a bit for loop to start
        time.sleep(0.1)

    def _get_loop(self) -> asyncio.AbstractEventLoop:
        """Get the event loop."""
        if AsyncProcessManager._loop:
            return AsyncProcessManager._loop

        try:
            return asyncio.get_running_loop()
        except RuntimeError:
            # Create and return new loop
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            return loop

    async def _get_next_pid(self) -> int:
        """Get next available process ID."""
        if self._lock:
            async with self._lock:
                pid = self._next_pid
                self._next_pid += 1
                return pid
        else:
            pid = self._next_pid
            self._next_pid += 1
            return pid

    def startProcess(
        self,
        command: str = None,
        args: List[str] = None,
        inputFile: str = None,
        logFile: str = None,
        interpreter: str = None,
        cwd: str = None,
        env: Dict[str, str] = None,
        handler: Any = None,
        timeout: int = None,
        ifAsync: bool = True,
        **kwargs
    ) -> int:
        """
        Start a subprocess (sync interface that schedules async work).

        Args:
            command: Command to execute
            args: Command arguments
            inputFile: File to pipe to stdin
            logFile: File to write stdout/stderr
            interpreter: Optional interpreter (e.g., 'python')
            cwd: Working directory
            env: Environment variables
            handler: Completion handler [callback, kwargs]
            timeout: Timeout in milliseconds (-1 for no timeout)
            ifAsync: Whether to run asynchronously

        Returns:
            Process ID (int)

        Note: This is the sync interface compatible with existing code.
              It schedules async work in the event loop.
        """
        # Get event loop
        loop = self._get_loop()

        # Schedule the async start
        future = asyncio.run_coroutine_threadsafe(
            self._startProcess_async(
                command=command,
                args=args,
                inputFile=inputFile,
                logFile=logFile,
                interpreter=interpreter,
                cwd=cwd,
                env=env,
                handler=handler,
                timeout=timeout,
                ifAsync=ifAsync,
                **kwargs
            ),
            loop
        )

        # Wait for PID (blocking)
        pid = future.result(timeout=5.0)
        return pid

    async def _startProcess_async(
        self,
        command: str,
        args: List[str] = None,
        inputFile: str = None,
        logFile: str = None,
        interpreter: str = None,
        cwd: str = None,
        env: Dict[str, str] = None,
        handler: Any = None,
        timeout: int = None,
        ifAsync: bool = True,
        **kwargs
    ) -> int:
        """
        Start a subprocess asynchronously.

        Returns:
            Process ID
        """
        # Get next PID
        pid = await self._get_next_pid()

        # Build command with interpreter
        if interpreter:
            full_command = [interpreter, command] + (args or [])
        else:
            full_command = [command] + (args or [])

        # Create process info
        proc_info = ProcessInfo(
            pid=pid,
            command=command,
            args=args or [],
            inputFile=inputFile,
            logFile=logFile,
            cwd=cwd,
            env=env or os.environ.copy(),
            handler=handler,
            timeout=timeout
        )

        self.processes[pid] = proc_info

        # Start monitoring task
        if ifAsync:
            asyncio.create_task(self._monitor_process(pid, full_command))
        else:
            # Run synchronously
            await self._monitor_process(pid, full_command)

        logger.info(f"Started process {pid}: {command} {' '.join(args or [])}")
        return pid

    async def _monitor_process(self, pid: int, full_command: List[str]):
        """
        Monitor a subprocess until completion.

        Args:
            pid: Process ID
            full_command: Full command with args
        """
        proc_info = self.processes[pid]

        try:
            # Acquire semaphore to limit concurrent processes
            if self._semaphore:
                await self._semaphore.acquire()

            proc_info.status = "running"
            proc_info.startTime = time.time()

            # Prepare stdin
            stdin_dest = asyncio.subprocess.PIPE if proc_info.inputFile else None

            # Prepare stdout/stderr - use direct file descriptors for real-time logging
            stdout_dest = None
            stderr_dest = None
            stdout_file = None
            stderr_file = None

            if proc_info.logFile:
                # Create parent directory if needed
                Path(proc_info.logFile).parent.mkdir(parents=True, exist_ok=True)
                # Open log file for direct writing (real-time, no buffering)
                stdout_file = open(proc_info.logFile, 'w')
                stdout_dest = stdout_file
                # Stderr to separate file
                err_file_path = str(Path(proc_info.logFile).with_suffix('')) + '_err.txt'
                stderr_file = open(err_file_path, 'w')
                stderr_dest = stderr_file
            else:
                # No log file - still need to capture output
                stdout_dest = asyncio.subprocess.PIPE
                stderr_dest = asyncio.subprocess.PIPE

            # Launch subprocess
            try:
                process = await asyncio.create_subprocess_exec(
                    *full_command,
                    stdin=stdin_dest,
                    stdout=stdout_dest,
                    stderr=stderr_dest,
                    cwd=proc_info.cwd,
                    env=proc_info.env
                )
            except Exception as e:
                # Close file handles if subprocess launch failed
                if stdout_file:
                    stdout_file.close()
                if stderr_file:
                    stderr_file.close()
                raise

            proc_info.process = process
            # Store file handles so we can close them later
            proc_info.stdout_file = stdout_file
            proc_info.stderr_file = stderr_file

            # Send input file if needed
            if proc_info.inputFile and process.stdin:
                try:
                    with open(proc_info.inputFile, 'rb') as f:
                        input_data = f.read()
                    process.stdin.write(input_data)
                    await process.stdin.drain()
                    process.stdin.close()
                except Exception as e:
                    logger.error(f"Error reading input file {proc_info.inputFile}: {e}")

            # Wait for completion (with timeout)
            timeout_seconds = None
            if proc_info.timeout and proc_info.timeout > 0:
                timeout_seconds = proc_info.timeout / 1000.0

            try:
                # Wait for process to complete (no buffering - logs written in real-time)
                await asyncio.wait_for(
                    process.wait(),
                    timeout=timeout_seconds
                )
            except asyncio.TimeoutError:
                logger.warning(f"Process {pid} timed out, killing...")
                process.kill()
                await process.wait()
                proc_info.status = "timeout"
                proc_info.exitCode = -1
                proc_info.exitStatus = 1
                await self._handle_finish(pid, -1, 1)
                return
            finally:
                # Close file handles now that process is done
                if stdout_file:
                    stdout_file.close()
                if stderr_file:
                    stderr_file.close()

            # Update process info
            proc_info.exitCode = process.returncode
            proc_info.exitStatus = 0 if process.returncode == 0 else 1
            proc_info.finishTime = time.time()
            proc_info.status = "finished" if process.returncode == 0 else "failed"

            logger.info(
                f"Process {pid} finished with code {process.returncode} "
                f"({proc_info.finishTime - proc_info.startTime:.2f}s)"
            )

            # Call completion handler
            await self._handle_finish(pid, process.returncode, proc_info.exitStatus)

        except Exception as e:
            logger.error(f"Error in process {pid}: {e}", exc_info=True)
            proc_info.status = "failed"
            proc_info.error = str(e)
            proc_info.exitCode = -1
            proc_info.exitStatus = 1
            await self._handle_finish(pid, -1, 1)

        finally:
            # Release semaphore
            if self._semaphore:
                self._semaphore.release()

    async def _handle_finish(self, pid: int, exitCode: int, exitStatus: int):
        """
        Handle process completion.

        Args:
            pid: Process ID
            exitCode: Exit code from process
            exitStatus: Exit status (0=success, 1=failure)
        """
        proc_info = self.processes.get(pid)
        if not proc_info:
            logger.warning(f"Process {pid} not found in _handle_finish")
            return

        handler = proc_info.handler
        if not handler:
            return

        try:
            # Handler format: [callback, kwargs] or just callback
            if isinstance(handler, list) and len(handler) >= 1:
                callback = handler[0]
                kwargs = handler[1] if len(handler) > 1 else {}
            else:
                callback = handler
                kwargs = {}

            # Call handler with pid as first argument
            if asyncio.iscoroutinefunction(callback):
                await callback(pid, **kwargs)
            else:
                # Run sync function in executor to avoid blocking
                loop = asyncio.get_event_loop()
                await loop.run_in_executor(
                    self._executor,
                    lambda: callback(pid, **kwargs)
                )

        except Exception as e:
            logger.error(f"Error calling handler for process {pid}: {e}", exc_info=True)

    def getJobData(self, pid: int, attribute: str = 'exitStatus') -> Any:
        """
        Get data about a process.

        Args:
            pid: Process ID
            attribute: Attribute name to get

        Returns:
            Attribute value or None
        """
        proc_info = self.processes.get(pid)
        if not proc_info:
            return None

        # Map attribute names to ProcessInfo fields
        attr_map = {
            'exitStatus': 'exitStatus',
            'exitCode': 'exitCode',
            'status': 'status',
            'command': 'command',
            'logFile': 'logFile',
            'startTime': 'startTime',
            'finishTime': 'finishTime',
            'error': 'error'
        }

        field_name = attr_map.get(attribute, attribute)
        return getattr(proc_info, field_name, None)

    def deleteJob(self, pid: int):
        """
        Delete job data for a process.

        Args:
            pid: Process ID
        """
        if pid in self.processes:
            del self.processes[pid]

    def getRunningProcesses(self) -> List[int]:
        """Get list of currently running process IDs."""
        return [
            pid for pid, info in self.processes.items()
            if info.status == "running"
        ]

    def killProcess(self, pid: int):
        """
        Kill a running process.

        Args:
            pid: Process ID to kill
        """
        proc_info = self.processes.get(pid)
        if not proc_info or not proc_info.process:
            logger.warning(f"Process {pid} not found or not running")
            return

        try:
            proc_info.process.kill()
            logger.info(f"Killed process {pid}")
        except Exception as e:
            logger.error(f"Error killing process {pid}: {e}")


# Global singleton accessor
def ASYNC_PROCESSMANAGER() -> AsyncProcessManager:
    """Get the singleton AsyncProcessManager instance."""
    return AsyncProcessManager()
