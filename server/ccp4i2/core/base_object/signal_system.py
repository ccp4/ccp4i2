"""
Modern Python replacement for Qt's Signal/Slot system.

This module provides a clean, type-safe signal/slot implementation that can replace
Qt's signal system in the CCP4i2 Django backend. It supports:

- Type-safe signal definitions
- Multiple slot connections
- Automatic disconnection on object destruction
- Weak references to prevent memory leaks
- Async/await support for modern Python patterns
- Thread-safe operations
- Signal chaining and filtering
"""

import asyncio
import logging
import threading
import weakref
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Optional,
    TypeVar,
    Generic,
    Union,
    Protocol,
    runtime_checkable,
)
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)

T = TypeVar("T")
CallableT = TypeVar("CallableT", bound=Callable[..., Any])


# Decorator system for signals and slots
class SlotInfo:
    """Information about a decorated slot method."""

    def __init__(
        self,
        result_type: type = None,
        arg_types: tuple = None,
        name: str = None,
        auto_connect: bool = True,
    ):
        self.result_type = result_type
        self.arg_types = arg_types or ()
        self.name = name
        self.auto_connect = auto_connect
        self.connections = []


def Slot(*arg_types, result: type = None, name: str = None, auto_connect: bool = True):
    """
    Decorator to mark methods as slots (equivalent to Qt's @Slot decorator).

    Args:
        *arg_types: Expected argument types for type checking
        result: Expected return type
        name: Custom slot name (defaults to function name)
        auto_connect: Whether to automatically connect when signals are created

    Usage:
        @Slot(str, int, result=bool)
        def my_slot(self, text: str, value: int) -> bool:
            return True

        @Slot()  # No arguments
        def simple_slot(self):
            pass
    """

    def decorator(func: CallableT) -> CallableT:
        slot_info = SlotInfo(
            result_type=result,
            arg_types=arg_types,
            name=name or func.__name__,
            auto_connect=auto_connect,
        )

        # Store slot metadata on the function
        func._slot_info = slot_info
        func._is_slot = True

        # Add type checking wrapper if types specified
        if arg_types:
            original_func = func

            def type_checked_wrapper(*args, **kwargs):
                # Skip 'self' argument for methods
                check_args = (
                    args[1:]
                    if len(args) > 0 and hasattr(args[0], func.__name__)
                    else args
                )

                if len(check_args) != len(arg_types):
                    raise TypeError(
                        f"Slot {func.__name__} expects {len(arg_types)} arguments, "
                        f"got {len(check_args)}"
                    )

                for i, (arg, expected_type) in enumerate(zip(check_args, arg_types)):
                    if not isinstance(arg, expected_type):
                        raise TypeError(
                            f"Slot {func.__name__} argument {i} expected {expected_type.__name__}, "
                            f"got {type(arg).__name__}"
                        )

                return original_func(*args, **kwargs)

            type_checked_wrapper._slot_info = slot_info
            type_checked_wrapper._is_slot = True
            type_checked_wrapper.__name__ = func.__name__
            type_checked_wrapper.__doc__ = func.__doc__

            return type_checked_wrapper

        return func

    return decorator


@runtime_checkable
class SlotCallable(Protocol):
    """Protocol for slot callables - functions that can be connected to signals."""

    def __call__(self, *args, **kwargs) -> Any: ...


@dataclass
class Connection:
    """Represents a connection between a signal and a slot."""

    slot: Union[Callable, weakref.ReferenceType]
    weak: bool = True
    once: bool = False
    thread_safe: bool = False
    priority: int = 0
    filter_func: Optional[Callable] = None
    connection_id: str = field(default_factory=lambda: f"conn_{id(object())}")

    def __post_init__(self):
        if self.weak and not isinstance(self.slot, weakref.ReferenceType):
            # Create weak reference if requested
            try:
                self.slot = weakref.ref(self.slot)
            except TypeError:
                # Some callables can't be weak-referenced (built-ins, lambdas)
                self.weak = False
                logger.warning(
                    f"Cannot create weak reference for {self.slot}, using strong reference"
                )

    def get_callable(self) -> Optional[Callable]:
        """Get the actual callable, handling weak references."""
        if self.weak and isinstance(self.slot, weakref.ReferenceType):
            return self.slot()
        return self.slot

    def is_valid(self) -> bool:
        """Check if the connection is still valid (not garbage collected)."""
        return self.get_callable() is not None


class Signal(Generic[T]):
    """
    Modern replacement for Qt's QSignal.

    Supports type-safe connections, weak references, async operations,
    and thread-safe emission.

    Example:
        # Define a signal
        data_changed = Signal[dict]()

        # Connect a slot
        def on_data_changed(data: dict):
            print(f"Data changed: {data}")

        data_changed.connect(on_data_changed)

        # Emit the signal
        data_changed.emit({"key": "value"})
    """

    def __init__(self, name: str = None):
        self._connections: List[Connection] = []
        self._name = name or f"Signal_{id(self)}"
        self._lock = threading.RLock()
        self._emission_count = 0

    @property
    def name(self) -> str:
        return self._name

    @property
    def connection_count(self) -> int:
        """Number of active connections."""
        with self._lock:
            return len([c for c in self._connections if c.is_valid()])

    def connect(
        self,
        slot: SlotCallable,
        weak: bool = False,
        once: bool = False,
        thread_safe: bool = False,
        priority: int = 0,
        filter_func: Optional[Callable] = None,
    ) -> str:
        """
        Connect a slot to this signal.

        Args:
            slot: The callable to connect
            weak: Use weak reference (default False - strong references prevent
                  garbage collection of signal handlers in async pipelines)
            once: Disconnect after first emission
            thread_safe: Execute slot in thread pool
            priority: Higher priority slots are called first
            filter_func: Optional filter function to control emission

        Returns:
            connection_id: Unique ID for this connection
        """
        with self._lock:
            connection = Connection(
                slot=slot,
                weak=weak,
                once=once,
                thread_safe=thread_safe,
                priority=priority,
                filter_func=filter_func,
            )

            # Insert in priority order (highest first)
            inserted = False
            for i, conn in enumerate(self._connections):
                if connection.priority > conn.priority:
                    self._connections.insert(i, connection)
                    inserted = True
                    break

            if not inserted:
                self._connections.append(connection)

            logger.debug(f"Connected {slot} to {self._name}")
            return connection.connection_id

    def disconnect(self, slot_or_id: Union[SlotCallable, str] = None) -> int:
        """
        Disconnect slot(s) from this signal.

        Args:
            slot_or_id: Specific slot callable or connection ID to disconnect.
                       If None, disconnects all slots.

        Returns:
            Number of connections removed
        """
        with self._lock:
            if slot_or_id is None:
                # Disconnect all
                count = len(self._connections)
                self._connections.clear()
                logger.debug(f"Disconnected all slots from {self._name}")
                return count

            removed = 0
            to_remove = []

            for i, connection in enumerate(self._connections):
                if isinstance(slot_or_id, str):
                    # Disconnect by connection ID
                    if connection.connection_id == slot_or_id:
                        to_remove.append(i)
                        removed += 1
                        break
                else:
                    # Disconnect by slot callable
                    slot = connection.get_callable()
                    if slot == slot_or_id:
                        to_remove.append(i)
                        removed += 1

            # Remove in reverse order to maintain indices
            for i in reversed(to_remove):
                del self._connections[i]

            logger.debug(f"Disconnected {removed} slots from {self._name}")
            return removed

    def emit(self, *args, **kwargs) -> List[Any]:
        """
        Emit the signal with given arguments.

        Returns:
            List of return values from all connected slots
        """
        with self._lock:
            self._emission_count += 1
            valid_connections = []
            results = []
            to_remove = []

            # Clean up invalid connections and collect valid ones
            for i, connection in enumerate(self._connections):
                if not connection.is_valid():
                    to_remove.append(i)
                    continue

                # Apply filter if present
                if connection.filter_func:
                    try:
                        if not connection.filter_func(*args, **kwargs):
                            continue
                    except Exception as e:
                        logger.error(f"Filter function error in {self._name}: {e}")
                        continue

                valid_connections.append(connection)

                # Mark for removal if it's a once-only connection
                if connection.once:
                    to_remove.append(i)

            # Remove invalid and once-only connections
            for i in reversed(to_remove):
                del self._connections[i]

            # Execute slots
            for connection in valid_connections:
                try:
                    slot = connection.get_callable()
                    if slot is None:
                        continue

                    if connection.thread_safe:
                        # Execute in thread pool (non-blocking)
                        future = ThreadPoolExecutor().submit(slot, *args, **kwargs)
                        results.append(future)
                    else:
                        # Execute synchronously
                        result = slot(*args, **kwargs)
                        results.append(result)

                except Exception as e:
                    import traceback
                    logger.error(f"Slot execution error in {self._name}: {e}")
                    logger.error(f"Traceback:\n{traceback.format_exc()}")
                    # Continue with other slots even if one fails

            logger.debug(f"Emitted {self._name} to {len(valid_connections)} slots")
            return results

    async def emit_async(self, *args, **kwargs) -> List[Any]:
        """
        Async version of emit that can handle async slots.
        """
        with self._lock:
            self._emission_count += 1
            valid_connections = [c for c in self._connections if c.is_valid()]
            results = []

        for connection in valid_connections:
            try:
                slot = connection.get_callable()
                if slot is None:
                    continue

                # Apply filter if present
                if connection.filter_func:
                    if not connection.filter_func(*args, **kwargs):
                        continue

                if asyncio.iscoroutinefunction(slot):
                    result = await slot(*args, **kwargs)
                else:
                    result = slot(*args, **kwargs)

                results.append(result)

            except Exception as e:
                logger.error(f"Async slot execution error in {self._name}: {e}")

        return results

    def __repr__(self) -> str:
        return f"Signal(name={self._name}, connections={self.connection_count})"


class SignalManager:
    """
    Manages signals for objects, providing automatic cleanup and organization.
    """

    def __init__(self):
        self._signals: Dict[str, Signal] = {}
        self._lock = threading.RLock()

    def create_signal(self, name: str, signal_type: type = None) -> Signal:
        """Create and register a new signal."""
        with self._lock:
            if name in self._signals:
                raise ValueError(f"Signal '{name}' already exists")

            signal = Signal[signal_type or Any](name=name)
            self._signals[name] = signal
            return signal

    def disconnect_all(self):
        """Disconnect all slots from all signals."""
        with self._lock:
            for signal in self._signals.values():
                signal.disconnect()

    def cleanup(self):
        """Full cleanup - disconnect and remove all signals."""
        with self._lock:
            self.disconnect_all()
            self._signals.clear()
