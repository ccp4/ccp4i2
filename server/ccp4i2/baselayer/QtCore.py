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
    from ccp4i2.core.base_object.event_system import EventLoop as ModernEventLoop
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

# Common property decorator
def Property(type, *args, **kwargs):
    """Stub Property decorator."""
    def decorator(func):
        return func
    return decorator


# Qt Event Loop - bridges to our modern event system
class _QEventLoop:
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
            self._event_loop = _QEventLoop()
        return self._event_loop.exec_()

    def exit(self, returnCode=0):
        """Exit the application with return code."""
        if self._event_loop:
            self._event_loop.exit(returnCode)

    def quit(self):
        """Exit the application with code 0."""
        self.exit(0)

    @staticmethod
    def processEvents(flags=_QEventLoop.AllEvents):
        """Process pending events (static Qt API)."""
        # This is a class method in Qt that processes events for the app
        loop = asyncio.get_event_loop()
        if loop.is_running():
            asyncio.create_task(asyncio.sleep(0))


# Alias for compatibility
QApplication = QCoreApplication
