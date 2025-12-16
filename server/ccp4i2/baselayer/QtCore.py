"""
Qt-free QtCore replacement module.

Provides Qt-compatible APIs that bridge to the modern event_system.py implementation.
This allows legacy CCP4i2 plugins to run without Qt dependencies while using
our async-first event loop architecture.
"""

from typing import Any, Callable, TypeVar

FuncT = TypeVar('FuncT', bound=Callable[..., Any])


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
