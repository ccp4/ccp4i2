"""
Tests for AsyncProcessManager.

Tests the new asyncio-based process manager that replaces Qt's QProcess.
"""

import pytest
import asyncio
import time
from pathlib import Path
from ccp4i2.core.async_process_manager import AsyncProcessManager


class TestAsyncProcessManager:
    """Test AsyncProcessManager basic functionality."""

    def test_singleton(self):
        """Test that AsyncProcessManager is a singleton."""
        pm1 = AsyncProcessManager()
        pm2 = AsyncProcessManager()
        assert pm1 is pm2

    @pytest.mark.asyncio
    async def test_simple_command_async(self):
        """Test running a simple command asynchronously."""
        pm = AsyncProcessManager()

        # Track completion
        completed = asyncio.Event()
        result_data = {}

        def handler(pid):
            result_data['pid'] = pid
            result_data['exitCode'] = pm.getJobData(pid, 'exitCode')
            result_data['status'] = pm.getJobData(pid, 'status')
            completed.set()

        # Start process
        pid = await pm._startProcess_async(
            command='echo',
            args=['hello world'],
            handler=[handler, {}],
            ifAsync=True
        )

        assert pid > 0

        # Wait for completion
        await asyncio.wait_for(completed.wait(), timeout=5.0)

        # Verify
        assert result_data['exitCode'] == 0
        assert result_data['status'] == 'finished'

    def test_simple_command_sync_interface(self):
        """Test sync interface to async process manager."""
        pm = AsyncProcessManager()

        # Track completion
        completion_flag = {'completed': False, 'pid': None}

        def handler(pid):
            completion_flag['pid'] = pid
            completion_flag['completed'] = True

        # Start process using sync interface
        pid = pm.startProcess(
            command='echo',
            args=['test message'],
            handler=[handler, {}],
            ifAsync=True
        )

        assert pid > 0

        # Wait for completion (polling)
        for _ in range(50):  # Wait up to 5 seconds
            if completion_flag['completed']:
                break
            time.sleep(0.1)

        assert completion_flag['completed']
        assert completion_flag['pid'] == pid
        assert pm.getJobData(pid, 'exitCode') == 0

    @pytest.mark.asyncio
    async def test_command_with_output_file(self):
        """Test command that writes to log file."""
        import tempfile

        pm = AsyncProcessManager()

        # Create temp directory
        with tempfile.TemporaryDirectory() as tmpdir:
            log_file = Path(tmpdir) / 'test.log'

            # Track completion
            completed = asyncio.Event()

            def handler(pid):
                completed.set()

            # Start process
            pid = await pm._startProcess_async(
                command='echo',
                args=['log test message'],
                logFile=str(log_file),
                handler=[handler, {}],
                ifAsync=True
            )

            # Wait for completion
            await asyncio.wait_for(completed.wait(), timeout=5.0)

            # Verify log file was created
            assert log_file.exists()

            # Verify content
            content = log_file.read_text()
            assert 'log test message' in content

    @pytest.mark.asyncio
    async def test_command_failure(self):
        """Test handling of failing command."""
        pm = AsyncProcessManager()

        # Track completion
        completed = asyncio.Event()
        result_data = {}

        def handler(pid):
            result_data['exitCode'] = pm.getJobData(pid, 'exitCode')
            result_data['exitStatus'] = pm.getJobData(pid, 'exitStatus')
            completed.set()

        # Run command that will fail
        pid = await pm._startProcess_async(
            command='false',  # Always returns exit code 1
            handler=[handler, {}],
            ifAsync=True
        )

        # Wait for completion
        await asyncio.wait_for(completed.wait(), timeout=5.0)

        # Verify failure was captured
        assert result_data['exitCode'] == 1
        assert result_data['exitStatus'] == 1

    def test_get_job_data(self):
        """Test retrieving job data."""
        pm = AsyncProcessManager()

        # Track completion
        completion_flag = {'completed': False}

        def handler(pid):
            completion_flag['completed'] = True

        # Start process
        pid = pm.startProcess(
            command='echo',
            args=['data test'],
            handler=[handler, {}],
            ifAsync=True
        )

        # Wait for completion
        for _ in range(50):
            if completion_flag['completed']:
                break
            time.sleep(0.1)

        # Test getJobData
        assert pm.getJobData(pid, 'command') == 'echo'
        assert pm.getJobData(pid, 'exitCode') == 0
        assert pm.getJobData(pid, 'exitStatus') == 0
        assert pm.getJobData(pid, 'status') == 'finished'

    def test_multiple_processes(self):
        """Test running multiple processes concurrently."""
        pm = AsyncProcessManager()

        # Track completions
        completions = {}

        def make_handler(proc_id):
            def handler(pid):
                completions[proc_id] = pid
            return handler

        # Start multiple processes
        pids = []
        for i in range(3):
            pid = pm.startProcess(
                command='echo',
                args=[f'process {i}'],
                handler=[make_handler(i), {}],
                ifAsync=True
            )
            pids.append(pid)

        # Wait for all to complete
        for _ in range(100):  # Wait up to 10 seconds
            if len(completions) == 3:
                break
            time.sleep(0.1)

        # Verify all completed
        assert len(completions) == 3
        for i in range(3):
            assert i in completions

    def test_delete_job(self):
        """Test deleting job data."""
        pm = AsyncProcessManager()

        completion_flag = {'completed': False}

        def handler(pid):
            completion_flag['completed'] = True

        # Start and complete a process
        pid = pm.startProcess(
            command='echo',
            args=['delete test'],
            handler=[handler, {}],
            ifAsync=True
        )

        # Wait for completion
        for _ in range(50):
            if completion_flag['completed']:
                break
            time.sleep(0.1)

        # Verify data exists
        assert pm.getJobData(pid, 'command') == 'echo'

        # Delete
        pm.deleteJob(pid)

        # Verify data is gone
        assert pm.getJobData(pid, 'command') is None

    def test_get_running_processes(self):
        """Test getting list of running processes."""
        pm = AsyncProcessManager()

        completion_flags = {}

        def make_handler(i):
            def handler(pid):
                completion_flags[i] = True
            return handler

        # Start a slow process (sleep for 2 seconds)
        pid1 = pm.startProcess(
            command='sleep',
            args=['2'],
            handler=[make_handler(1), {}],
            ifAsync=True
        )

        # Give it a moment to start
        time.sleep(0.2)

        # Check running processes
        running = pm.getRunningProcesses()
        assert pid1 in running

        # Wait for completion
        time.sleep(2.5)

        # Check it's no longer running
        running = pm.getRunningProcesses()
        assert pid1 not in running


if __name__ == '__main__':
    pytest.main([__file__, '-v', '-s'])
