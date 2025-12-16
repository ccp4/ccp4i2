"""
Modern Python replacement for Qt's event loop and task execution system.

This module provides a clean, async-first event loop implementation that can replace
Qt's QEventLoop and QApplication functionality in the CCP4i2 Django backend. It supports:

- Async/await task execution
- Thread pool integration
- Signal-based inter-task communication
- Task lifecycle management
- Error handling and recovery
- Task monitoring and debugging
- Graceful shutdown handling
"""

import asyncio
import logging
import signal
import sys
import threading
import time
import traceback
import weakref
from concurrent.futures import ThreadPoolExecutor, Future
from contextlib import asynccontextmanager
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import (
    Any,
    Awaitable,
    Callable,
    Dict,
    List,
    Optional,
    Set,
    Union,
    TypeVar,
    Coroutine,
)
import uuid

# Import will be done relatively when used as module
try:
    from .signal_system import Signal, SignalManager
    from .hierarchy_system import HierarchicalObject
except ImportError:
    # For standalone testing
    from signal_system import Signal, SignalManager
    from hierarchy_system import HierarchicalObject

logger = logging.getLogger(__name__)

T = TypeVar("T")


class TaskState(Enum):
    """Task execution states."""

    PENDING = auto()
    RUNNING = auto()
    COMPLETED = auto()
    FAILED = auto()
    CANCELLED = auto()


@dataclass
class TaskInfo:
    """Information about a running task."""

    task_id: str
    name: str
    state: TaskState = TaskState.PENDING
    created_at: float = field(default_factory=time.time)
    started_at: Optional[float] = None
    completed_at: Optional[float] = None
    result: Any = None
    error: Optional[Exception] = None
    progress: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


class TaskRunner(HierarchicalObject):
    """
    Individual task runner that can execute async or sync functions.
    """

    def __init__(
        self, task_func: Callable, name: str = None, parent: "EventLoop" = None
    ):
        super().__init__(parent, name or f"Task_{uuid.uuid4().hex[:8]}")

        self.task_func = task_func
        self.task_info = TaskInfo(task_id=self.name, name=self.name)

        # Signals
        self.task_started = self.create_signal("task_started", TaskInfo)
        self.task_progress = self.create_signal("task_progress", dict)
        self.task_completed = self.create_signal("task_completed", TaskInfo)
        self.task_failed = self.create_signal("task_failed", TaskInfo)

        self._future: Optional[Future] = None
        self._cancelled = False

    async def execute(self, *args, **kwargs) -> Any:
        """Execute the task function."""
        if self._cancelled:
            self.task_info.state = TaskState.CANCELLED
            return None

        self.task_info.state = TaskState.RUNNING
        self.task_info.started_at = time.time()
        self.task_started.emit(self.task_info)

        try:
            if asyncio.iscoroutinefunction(self.task_func):
                result = await self.task_func(*args, **kwargs)
            else:
                # Run sync function in thread pool
                loop = asyncio.get_event_loop()
                result = await loop.run_in_executor(
                    None, self.task_func, *args, **kwargs
                )

            self.task_info.result = result
            self.task_info.state = TaskState.COMPLETED
            self.task_info.completed_at = time.time()
            self.task_completed.emit(self.task_info)

            return result

        except Exception as e:
            self.task_info.error = e
            self.task_info.state = TaskState.FAILED
            self.task_info.completed_at = time.time()
            self.task_failed.emit(self.task_info)
            raise

    def cancel(self) -> bool:
        """Cancel the task if possible."""
        self._cancelled = True
        if self._future and not self._future.done():
            return self._future.cancel()
        return True

    def update_progress(self, progress: float, message: str = None):
        """Update task progress (0.0 to 1.0)."""
        self.task_info.progress = max(0.0, min(1.0, progress))
        self.task_progress.emit(
            {
                "progress": self.task_info.progress,
                "message": message,
                "task_id": self.task_info.task_id,
            }
        )


class EventLoop(HierarchicalObject):
    """
    Modern replacement for Qt's QEventLoop.

    Provides async task execution, signal handling, and lifecycle management
    without Qt dependencies.
    """

    def __init__(self, name: str = "EventLoop"):
        super().__init__(name=name)

        self._loop: Optional[asyncio.AbstractEventLoop] = None
        self._thread_pool = ThreadPoolExecutor(max_workers=4)
        self._running = False
        self._shutdown_event = asyncio.Event()
        self._tasks: Dict[str, TaskRunner] = {}
        self._task_futures: Dict[str, Future] = {}

        # Signals
        self.loop_started = self.create_signal("loop_started")
        self.loop_stopped = self.create_signal("loop_stopped")
        self.task_scheduled = self.create_signal("task_scheduled", TaskRunner)
        self.error_occurred = self.create_signal("error_occurred", dict)

        # Setup signal handlers for graceful shutdown
        self._setup_signal_handlers()

    def _setup_signal_handlers(self):
        """Setup system signal handlers for graceful shutdown."""

        def signal_handler(signum, frame):
            logger.info("Received shutdown signal %d", signum)
            asyncio.create_task(self.shutdown())

        if sys.platform != "win32":
            signal.signal(signal.SIGINT, signal_handler)
            signal.signal(signal.SIGTERM, signal_handler)

    @property
    def is_running(self) -> bool:
        """Check if the event loop is running."""
        return self._running and self._loop is not None and self._loop.is_running()

    @property
    def task_count(self) -> int:
        """Number of active tasks."""
        return len(self._tasks)

    def get_task_info(
        self, task_id: str = None
    ) -> Union[TaskInfo, List[TaskInfo], None]:
        """Get information about tasks."""
        if task_id:
            task = self._tasks.get(task_id)
            return task.task_info if task else None
        else:
            return [task.task_info for task in self._tasks.values()]

    def schedule_task(
        self, task_func: Callable, *args, name: str = None, **kwargs
    ) -> str:
        """
        Schedule a task for execution.

        Args:
            task_func: Function to execute (sync or async)
            *args: Arguments for the function
            name: Optional task name
            **kwargs: Keyword arguments for the function

        Returns:
            Task ID string
        """
        task_runner = TaskRunner(task_func, name, parent=self)
        self._tasks[task_runner.name] = task_runner

        # Connect to task signals for monitoring
        task_runner.task_completed.connect(self._on_task_completed)
        task_runner.task_failed.connect(self._on_task_failed)

        # Schedule the task
        if self._loop:
            future = asyncio.ensure_future(
                task_runner.execute(*args, **kwargs), loop=self._loop
            )
            self._task_futures[task_runner.name] = future

        self.task_scheduled.emit(task_runner)
        logger.info("Scheduled task %s", task_runner.name)

        return task_runner.name

    def cancel_task(self, task_id: str) -> bool:
        """Cancel a specific task."""
        task = self._tasks.get(task_id)
        if task:
            success = task.cancel()
            if success:
                self._cleanup_task(task_id)
            return success
        return False

    def _on_task_completed(self, task_info: TaskInfo):
        """Handle task completion."""
        logger.info("Task %s completed successfully", task_info.task_id)
        self._cleanup_task(task_info.task_id)

    def _on_task_failed(self, task_info: TaskInfo):
        """Handle task failure."""
        logger.error("Task %s failed: %s", task_info.task_id, task_info.error)
        self.error_occurred.emit(
            {
                "task_id": task_info.task_id,
                "error": task_info.error,
                "traceback": traceback.format_exc() if task_info.error else None,
            }
        )
        self._cleanup_task(task_info.task_id)

    def _cleanup_task(self, task_id: str):
        """Clean up completed task."""
        self._tasks.pop(task_id, None)
        future = self._task_futures.pop(task_id, None)
        if future and not future.done():
            future.cancel()

    async def run_until_complete(self, coro: Awaitable[T]) -> T:
        """Run a coroutine until completion."""
        if not self.is_running:
            raise RuntimeError("Event loop is not running")

        return await coro

    def run_in_executor(self, func: Callable, *args) -> Future:
        """Run a function in the thread pool executor."""
        if not self._loop:
            raise RuntimeError("Event loop is not running")

        return self._loop.run_in_executor(self._thread_pool, func, *args)

    async def start(self):
        """Start the event loop."""
        if self._running:
            logger.warning("Event loop is already running")
            return

        self._loop = asyncio.get_event_loop()
        self._running = True
        self._shutdown_event.clear()

        logger.info("Starting event loop %s", self.name)
        self.loop_started.emit()

    async def stop(self):
        """Stop the event loop gracefully."""
        if not self._running:
            return

        logger.info("Stopping event loop %s", self.name)

        # Cancel all running tasks
        for task_id in list(self._tasks.keys()):
            self.cancel_task(task_id)

        # Wait a bit for tasks to finish
        await asyncio.sleep(0.1)

        self._running = False
        self._shutdown_event.set()
        self.loop_stopped.emit()

    async def shutdown(self):
        """Shutdown the event loop and cleanup resources."""
        await self.stop()

        # Shutdown thread pool
        self._thread_pool.shutdown(wait=True)

        # Cleanup
        self._tasks.clear()
        self._task_futures.clear()

        logger.info("Event loop %s shutdown complete", self.name)

    async def wait_for_shutdown(self):
        """Wait for shutdown signal."""
        await self._shutdown_event.wait()

    def run_forever(self):
        """Run the event loop forever (blocking)."""

        async def main():
            await self.start()
            await self.wait_for_shutdown()
            await self.shutdown()

        try:
            asyncio.run(main())
        except KeyboardInterrupt:
            logger.info("Received keyboard interrupt")
        except Exception as e:
            logger.error("Event loop error: %s", e)
            raise


class Application(EventLoop):
    """
    Top-level application class, similar to Qt's QApplication.

    Manages the main event loop and provides application-level services.
    """

    def __init__(self, name: str = "CCP4i2Application"):
        super().__init__(name)

        self._plugins: Dict[str, Any] = {}
        self._config: Dict[str, Any] = {}

        # Application signals
        self.plugin_loaded = self.create_signal("plugin_loaded", dict)
        self.config_changed = self.create_signal("config_changed", dict)

    def load_plugin(self, plugin_name: str, plugin_class: type, **config):
        """Load and initialize a plugin."""
        try:
            plugin_instance = plugin_class(parent=self, **config)
            self._plugins[plugin_name] = plugin_instance

            self.plugin_loaded.emit(
                {
                    "name": plugin_name,
                    "class": plugin_class.__name__,
                    "instance": plugin_instance,
                }
            )

            logger.info("Loaded plugin %s", plugin_name)
            return plugin_instance

        except Exception as e:
            logger.error("Failed to load plugin %s: %s", plugin_name, e)
            raise

    def get_plugin(self, plugin_name: str) -> Optional[Any]:
        """Get a loaded plugin by name."""
        return self._plugins.get(plugin_name)

    def unload_plugin(self, plugin_name: str) -> bool:
        """Unload a plugin."""
        plugin = self._plugins.pop(plugin_name, None)
        if plugin and hasattr(plugin, "destroy"):
            plugin.destroy()
        return plugin is not None

    def set_config(self, key: str, value: Any):
        """Set application configuration."""
        old_value = self._config.get(key)
        self._config[key] = value
        self.config_changed.emit(
            {"key": key, "old_value": old_value, "new_value": value}
        )

    def get_config(self, key: str, default: Any = None) -> Any:
        """Get application configuration."""
        return self._config.get(key, default)


# Example usage and testing
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Create application
    app = Application("TestApp")

    # Define some test tasks
    async def async_task(duration: float, name: str):
        """Test async task."""
        print(f"Starting async task {name}")
        for i in range(5):
            await asyncio.sleep(duration / 5)
            # Simulate progress updates if we had access to the task runner
            print(f"Task {name} progress: {(i + 1) * 20}%")
        print(f"Async task {name} completed")
        return f"Result from {name}"

    def sync_task(duration: float, name: str):
        """Test sync task."""
        print(f"Starting sync task {name}")
        time.sleep(duration)
        print(f"Sync task {name} completed")
        return f"Result from {name}"

    # Schedule some tasks
    async def test_application():
        await app.start()

        # Schedule tasks
        task1_id = app.schedule_task(async_task, 2.0, "AsyncTask1")
        task2_id = app.schedule_task(sync_task, 1.0, "SyncTask1")

        print(f"Scheduled tasks: {task1_id}, {task2_id}")

        # Wait a bit then check status
        await asyncio.sleep(0.5)

        task_infos = app.get_task_info()
        print(f"Active tasks: {len(task_infos)}")

        # Wait for tasks to complete
        await asyncio.sleep(3.0)

        await app.shutdown()

    # Run the test
    try:
        asyncio.run(test_application())
    except KeyboardInterrupt:
        print("Test interrupted")
