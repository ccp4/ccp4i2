"""
Qt-free QtCore replacement module.

Provides Qt-compatible APIs that bridge to the modern event_system.py implementation.
This allows legacy CCP4i2 plugins to run without Qt dependencies while using
our async-first event loop architecture.
"""

from typing import Any, Callable, TypeVar
import sys
import asyncio

FuncT = TypeVar('FuncT', bound=Callable[..., Any])

# Import our modern event system
try:
    from core.base_object.event_system import EventLoop as ModernEventLoop
except ImportError:
    # Fallback for module loading issues
    ModernEventLoop = None


def Slot(*arg_types, **kwargs) -> Callable[[FuncT], FuncT]:
    """
    Stub for Qt's Slot decorator.

    This is a no-op decorator that allows code using @Slot to be imported.
    It's compatible with our signal_system.Slot if needed, but for plugin
    discovery we just pass through the original function unchanged.

    Usage (from CCP4i2 plugins):
        @Slot()
        def my_slot(self):
            pass

        @Slot(str, int)
        def handle_data(self, text, value):
            pass
    """
    def decorator(func: FuncT) -> FuncT:
        # Just return the function unchanged - we only need imports to work
        return func
    return decorator


# Signal class with functional connect/emit support
class Signal:
    """
    Qt-compatible Signal class.

    Supports connecting slots and emitting signals for Qt-free operation.
    """
    def __init__(self, *args):
        """
        Initialize a signal.

        Args:
            *args: Type arguments (for Qt compatibility, ignored in our implementation)
        """
        self._slots = []
        self._arg_types = args

    def connect(self, slot):
        """
        Connect a slot (callback function) to this signal.

        Args:
            slot: Callable to invoke when signal is emitted
        """
        if slot not in self._slots:
            self._slots.append(slot)

    def disconnect(self, slot=None):
        """
        Disconnect a slot from this signal.

        Args:
            slot: Specific slot to disconnect, or None to disconnect all
        """
        if slot is None:
            self._slots.clear()
        elif slot in self._slots:
            self._slots.remove(slot)

    def emit(self, *args):
        """
        Emit the signal, calling all connected slots.

        Args:
            *args: Arguments to pass to connected slots
        """
        for slot in self._slots[:]:  # Copy list to avoid modification during iteration
            try:
                slot(*args)
            except Exception as e:
                import sys
                print(f"Error in signal slot: {e}", file=sys.stderr)
                import traceback
                traceback.print_exc()


# Stub QObject - base class for Qt objects
class QObject:
    """Stub QObject class."""
    def __init__(self, parent=None):
        self._parent = parent

    def parent(self):
        """Return the parent object (Qt API)."""
        return self._parent

    def setParent(self, parent):
        """Set the parent object (Qt API)."""
        self._parent = parent


# Other common QtCore items that might be imported
class QThread:
    """Stub QThread class."""
    pass


class QTimer:
    """Stub QTimer class."""
    pass


class Qt:
    """Stub Qt namespace class."""
    # Common Qt enum values
    AlignLeft = 0x0001
    AlignRight = 0x0002
    AlignCenter = 0x0004
    AlignTop = 0x0020
    AlignBottom = 0x0040


class QFileSystemWatcher:
    """
    Qt-free file system watcher implementation.

    Monitors directories and files for changes using platform-independent polling.
    Provides Qt-compatible API for legacy CCP4i2 plugins.
    """

    def __init__(self, paths=None, parent=None):
        """
        Initialize the file system watcher.

        Args:
            paths: Optional list of paths to watch
            parent: Optional parent object (Qt compatibility)
        """
        self._parent = parent
        self._watched_paths = set()
        self._watched_files = set()
        self._last_modified = {}
        self._monitoring = False
        self._monitor_task = None

        # Signals (Qt API)
        self.directoryChanged = Signal(str)
        self.fileChanged = Signal(str)

        if paths:
            for path in paths:
                self.addPath(path)

    def addPath(self, path):
        """
        Add a path to watch.

        Args:
            path: File or directory path to monitor

        Returns:
            bool: True if successfully added
        """
        import os

        if not os.path.exists(path):
            return False

        if os.path.isdir(path):
            self._watched_paths.add(os.path.abspath(path))
            self._last_modified[path] = self._get_dir_state(path)
        else:
            self._watched_files.add(os.path.abspath(path))
            self._last_modified[path] = os.path.getmtime(path)

        # Start monitoring if not already running
        if not self._monitoring:
            self._start_monitoring()

        return True

    def addPaths(self, paths):
        """
        Add multiple paths to watch.

        Args:
            paths: List of paths to monitor

        Returns:
            list: Paths that failed to be added
        """
        failed = []
        for path in paths:
            if not self.addPath(path):
                failed.append(path)
        return failed

    def removePath(self, path):
        """
        Remove a path from monitoring.

        Args:
            path: Path to stop watching

        Returns:
            bool: True if successfully removed
        """
        import os
        abs_path = os.path.abspath(path)
        removed = False

        if abs_path in self._watched_paths:
            self._watched_paths.discard(abs_path)
            removed = True

        if abs_path in self._watched_files:
            self._watched_files.discard(abs_path)
            removed = True

        if abs_path in self._last_modified:
            del self._last_modified[abs_path]

        # Stop monitoring if no paths left
        if not self._watched_paths and not self._watched_files:
            self._stop_monitoring()

        return removed

    def removePaths(self, paths):
        """
        Remove multiple paths from monitoring.

        Args:
            paths: List of paths to stop watching

        Returns:
            list: Paths that failed to be removed
        """
        failed = []
        for path in paths:
            if not self.removePath(path):
                failed.append(path)
        return failed

    def directories(self):
        """
        Get list of watched directories.

        Returns:
            list: Currently watched directory paths
        """
        return list(self._watched_paths)

    def files(self):
        """
        Get list of watched files.

        Returns:
            list: Currently watched file paths
        """
        return list(self._watched_files)

    def _get_dir_state(self, path):
        """
        Get the current state of a directory for change detection.

        Args:
            path: Directory path

        Returns:
            dict: Mapping of filenames to modification times
        """
        import os

        try:
            state = {}
            for entry in os.listdir(path):
                entry_path = os.path.join(path, entry)
                try:
                    state[entry] = os.path.getmtime(entry_path)
                except (OSError, IOError):
                    pass
            return state
        except (OSError, IOError):
            return {}

    def _check_changes(self):
        """Check for changes in watched paths and emit signals."""
        import os

        # Check directories
        for dir_path in list(self._watched_paths):
            if not os.path.exists(dir_path):
                continue

            current_state = self._get_dir_state(dir_path)
            last_state = self._last_modified.get(dir_path, {})

            if current_state != last_state:
                self._last_modified[dir_path] = current_state
                # Emit signal (in Qt-free mode, this is a no-op unless connected)
                if hasattr(self.directoryChanged, 'emit'):
                    self.directoryChanged.emit(dir_path)
                # Also try calling connected slots directly
                if hasattr(self, '_directory_changed_slots'):
                    for slot in self._directory_changed_slots:
                        try:
                            slot(dir_path)
                        except Exception as e:
                            import sys
                            print(f"Error in directory changed slot: {e}", file=sys.stderr)

        # Check files
        for file_path in list(self._watched_files):
            if not os.path.exists(file_path):
                continue

            try:
                current_mtime = os.path.getmtime(file_path)
                last_mtime = self._last_modified.get(file_path, 0)

                if current_mtime != last_mtime:
                    self._last_modified[file_path] = current_mtime
                    if hasattr(self.fileChanged, 'emit'):
                        self.fileChanged.emit(file_path)
                    if hasattr(self, '_file_changed_slots'):
                        for slot in self._file_changed_slots:
                            try:
                                slot(file_path)
                            except Exception as e:
                                import sys
                                print(f"Error in file changed slot: {e}", file=sys.stderr)
            except (OSError, IOError):
                pass

    def _start_monitoring(self):
        """Start the monitoring loop."""
        import threading
        import time

        if self._monitoring:
            return

        self._monitoring = True

        def monitor_loop():
            """Background thread that polls for changes."""
            while self._monitoring:
                self._check_changes()
                time.sleep(0.5)  # Poll every 500ms

        # Start monitoring in a background thread
        self._monitor_task = threading.Thread(target=monitor_loop, daemon=True)
        self._monitor_task.start()

    def _stop_monitoring(self):
        """Stop the monitoring loop."""
        self._monitoring = False
        if self._monitor_task:
            self._monitor_task = None


# Common property decorator
def Property(type, *args, **kwargs):
    """Stub Property decorator."""
    def decorator(func):
        return func
    return decorator


# Qt Event Loop - bridges to our modern event system
class QEventLoop:
    """
    Qt-compatible event loop that bridges to core.base_object.event_system.EventLoop.

    This provides the Qt QEventLoop API while using our modern async event system
    under the hood, allowing legacy CCP4i2 plugins to run without Qt.
    """

    # Exit codes matching Qt's QEventLoop
    AllEvents = 0x00
    ExcludeUserInputEvents = 0x01
    ExcludeSocketNotifiers = 0x02
    WaitForMoreEvents = 0x04

    def __init__(self, parent=None):
        self._parent = parent
        self._modern_loop = None
        self._exit_code = 0
        self._is_running = False

        # Create our modern event loop if available
        if ModernEventLoop:
            self._modern_loop = ModernEventLoop(name="QtCompat_EventLoop")

    def parent(self):
        """Return the parent object (Qt API)."""
        return self._parent

    def setParent(self, parent):
        """Set the parent object (Qt API)."""
        self._parent = parent

    def exec(self, flags=0):
        """
        Start the event loop (Qt API: exec_()).

        This is the main entry point for running async tasks in legacy code.
        Bridges to our modern EventLoop.run_forever().

        Returns:
            Exit code (0 for success)
        """
        return self.exec_(flags)

    def exec_(self, flags=0):
        """
        Start the event loop (Qt naming convention).

        Runs the modern event loop until exit() is called.
        """
        if self._modern_loop:
            self._is_running = True
            try:
                # Run the modern event loop
                self._modern_loop.run_forever()
            except KeyboardInterrupt:
                pass
            finally:
                self._is_running = False
        else:
            # Fallback: just run asyncio event loop
            try:
                self._is_running = True
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                loop.run_forever()
            except KeyboardInterrupt:
                pass
            finally:
                self._is_running = False
                loop.close()

        return self._exit_code

    def exit(self, returnCode=0):
        """
        Exit the event loop with the specified return code.

        Args:
            returnCode: Exit code to return from exec_()
        """
        self._exit_code = returnCode

        if self._modern_loop and hasattr(self._modern_loop, 'shutdown'):
            # Use modern loop's shutdown mechanism
            asyncio.create_task(self._modern_loop.shutdown())
        else:
            # Fallback: stop asyncio loop
            loop = asyncio.get_event_loop()
            loop.stop()

    def quit(self):
        """Exit the event loop with code 0 (Qt convenience method)."""
        self.exit(0)

    def isRunning(self):
        """Check if the event loop is running."""
        return self._is_running

    def processEvents(self, flags=AllEvents):
        """
        Process pending events (Qt API).

        In our async system, this yields control to allow other tasks to run.
        """
        if self._modern_loop:
            # Let other tasks run
            loop = asyncio.get_event_loop()
            if loop.is_running():
                # Schedule a brief pause to process events
                asyncio.create_task(asyncio.sleep(0))


class QCoreApplication:
    """
    Qt-compatible application object.

    Provides Qt's QCoreApplication API while using our event system.
    This is a singleton that manages the application lifecycle.
    """

    _instance = None

    def __init__(self, argv=None):
        if QCoreApplication._instance is not None:
            # Qt behavior: only one QCoreApplication instance allowed
            raise RuntimeError("QCoreApplication instance already exists")

        QCoreApplication._instance = self
        self.argv = argv or sys.argv
        self._event_loop = None

    @classmethod
    def instance(cls):
        """Get the singleton QCoreApplication instance (Qt API)."""
        return cls._instance

    def exec(self):
        """Start the application event loop."""
        return self.exec_()

    def exec_(self):
        """Start the application event loop (Qt naming)."""
        if self._event_loop is None:
            self._event_loop = QEventLoop()
        return self._event_loop.exec_()

    def exit(self, returnCode=0):
        """Exit the application with return code."""
        if self._event_loop:
            self._event_loop.exit(returnCode)

    def quit(self):
        """Exit the application with code 0."""
        self.exit(0)

    @staticmethod
    def processEvents(flags=QEventLoop.AllEvents):
        """Process pending events (static Qt API)."""
        # This is a class method in Qt that processes events for the app
        loop = asyncio.get_event_loop()
        if loop.is_running():
            asyncio.create_task(asyncio.sleep(0))


# Alias for compatibility
QApplication = QCoreApplication
