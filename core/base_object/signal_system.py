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


class SignalInfo:
    """Information about a decorated signal attribute."""

    def __init__(
        self, signal_type: type = None, arg_types: tuple = None, name: str = None
    ):
        self.signal_type = signal_type
        self.arg_types = arg_types or ()
        self.name = name


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


def SignalDecorator(signal_type: type = None, *arg_types, name: str = None):
    """
    Decorator/descriptor for creating signals (equivalent to Qt's Signal).

    Can be used as:
    1. Class attribute: my_signal = SignalDecorator(str)
    2. Method decorator: @SignalDecorator(str, int)

    Args:
        signal_type: Primary signal data type
        *arg_types: Additional argument types
        name: Custom signal name

    Usage:
        class MyClass:
            # As class attribute
            value_changed = SignalDecorator(int)
            data_ready = SignalDecorator(dict)

            # As property descriptor
            @SignalDecorator(str)
            def message_sent(self): pass
    """

    class SignalDescriptor:
        """Descriptor that creates Signal instances on first access."""

        def __init__(
            self, signal_type: type = None, arg_types: tuple = (), name: str = None
        ):
            self.signal_type = signal_type
            self.arg_types = arg_types
            self.name = name
            self._signals = weakref.WeakKeyDictionary()

        def __get__(self, obj, objtype=None):
            if obj is None:
                return self

            # Return placeholder - will be replaced when Signal class is defined
            return None

        def __set_name__(self, owner, name):
            if self.name is None:
                self.name = name

    # If used as decorator (@SignalDecorator(int))
    if callable(signal_type):
        func = signal_type
        signal_info = SignalInfo(name=name or func.__name__)

        # Mark function as signal creator
        func._signal_info = signal_info
        func._is_signal = True

        return SignalDescriptor(name=signal_info.name)

    # If used as descriptor (SignalDecorator(int))
    all_types = (signal_type,) + arg_types if signal_type else arg_types
    return SignalDescriptor(signal_type=signal_type, arg_types=all_types, name=name)


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
        self._blocked = False
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
        if self._blocked:
            return []

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
        if self._blocked:
            return []

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

    @contextmanager
    def blocked(self):
        """Context manager to temporarily block signal emission."""
        self._blocked = True
        try:
            yield
        finally:
            self._blocked = False

    def __repr__(self) -> str:
        return f"Signal(name={self._name}, connections={self.connection_count})"


# Utility functions for working with decorated slots and signals
def get_slots(obj) -> Dict[str, SlotInfo]:
    """Get all @Slot decorated methods from an object."""
    slots = {}
    for name in dir(obj):
        attr = getattr(obj, name)
        if callable(attr) and hasattr(attr, "_is_slot"):
            slots[name] = attr._slot_info
    return slots


def get_signals(obj) -> Dict[str, SignalInfo]:
    """Get all @SignalDecorator decorated attributes from an object."""
    signals = {}
    for name in dir(type(obj)):
        attr = getattr(type(obj), name)
        if hasattr(attr, "_is_signal"):
            signals[name] = attr._signal_info
    return signals


def auto_connect_slots(sender_obj, receiver_obj, signal_prefix: str = ""):
    """
    Automatically connect signals to slots based on naming convention.

    Connects signals like 'value_changed' to slots like 'on_value_changed'
    or 'handle_value_changed'.

    Args:
        sender_obj: Object with signals
        receiver_obj: Object with slots
        signal_prefix: Optional prefix for signal names
    """
    sender_signals = get_signals(sender_obj)
    receiver_slots = get_slots(receiver_obj)

    connections = []

    for signal_name, signal_info in sender_signals.items():
        # Try different slot naming patterns
        slot_patterns = [
            f"on_{signal_name}",
            f"handle_{signal_name}",
            f"{signal_name}_handler",
            signal_name,  # Exact match
        ]

        if signal_prefix:
            slot_patterns.extend(
                [
                    f"on_{signal_prefix}_{signal_name}",
                    f"handle_{signal_prefix}_{signal_name}",
                ]
            )

        for slot_pattern in slot_patterns:
            if slot_pattern in receiver_slots:
                # Get the actual signal and slot
                signal_obj = getattr(sender_obj, signal_name)
                slot_method = getattr(receiver_obj, slot_pattern)

                if isinstance(signal_obj, Signal):
                    conn = signal_obj.connect(slot_method)
                    connections.append(conn)
                    logger.info(f"Auto-connected {signal_name} -> {slot_pattern}")
                    break

    return connections


def create_signal_from_decorator(obj, attr_name: str, descriptor) -> Signal:
    """Create a Signal instance from a SignalDecorator descriptor."""
    if hasattr(descriptor, "signal_type"):
        signal = Signal(descriptor.signal_type, name=attr_name)
        # Store on object
        setattr(obj, f"_{attr_name}_signal", signal)
        return signal
    return None


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

    def get_signal(self, name: str) -> Optional[Signal]:
        """Get an existing signal by name."""
        return self._signals.get(name)

    def remove_signal(self, name: str) -> bool:
        """Remove and cleanup a signal."""
        with self._lock:
            if name in self._signals:
                signal = self._signals[name]
                signal.disconnect()  # Disconnect all slots
                del self._signals[name]
                return True
            return False

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

    @property
    def signal_names(self) -> List[str]:
        """Get list of all signal names."""
        return list(self._signals.keys())


# Decorator for automatic signal creation
def signal(signal_type: type = None, name: str = None):
    """
    Decorator to automatically create signals as class attributes.

    Example:
        class MyClass(SignalEmitter):
            @signal(dict)
            def data_changed(self): pass

            def some_method(self):
                self.data_changed.emit({"key": "value"})
    """

    def decorator(func: Callable) -> Signal:
        signal_name = name or func.__name__
        return Signal[signal_type or Any](name=signal_name)

    return decorator


# Example usage and testing
if __name__ == "__main__":
    # Demo of the signal system
    logging.basicConfig(level=logging.DEBUG)

    # Create a signal
    data_changed = Signal[dict](name="data_changed")

    # Define some slots
    def slot1(data):
        print(f"Slot 1 received: {data}")
        return "slot1_result"

    def slot2(data):
        print(f"Slot 2 received: {data}")
        return "slot2_result"

    async def async_slot(data):
        print(f"Async slot received: {data}")
        await asyncio.sleep(0.1)
        return "async_result"

    # Connect slots
    conn1 = data_changed.connect(slot1)
    conn2 = data_changed.connect(slot2, priority=10)  # Higher priority

    # Emit signal
    results = data_changed.emit({"test": "data"})
    print(f"Results: {results}")

    # Disconnect specific slot
    data_changed.disconnect(conn1)

    # Test async emission
    async def test_async():
        data_changed.connect(async_slot)
        results = await data_changed.emit_async({"async": "test"})
        print(f"Async results: {results}")

    # Run async test
    asyncio.run(test_async())
