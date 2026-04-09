"""
Modern Python replacement for Qt's parent/child object hierarchy.

This module provides a clean, memory-safe object hierarchy system that can replace
Qt's QObject parent/child relationships in the CCP4i2 Django backend. It supports:

- Automatic lifecycle management
- Hierarchical object organization
- Event propagation through the hierarchy
- Automatic cleanup on parent destruction
- Weak references to prevent cycles
- Type-safe parent/child relationships
"""

import logging
import threading
import weakref
from abc import ABC
from typing import Any, Dict, List, Optional, Type, TypeVar, Callable, Set
from dataclasses import dataclass, field
from enum import Enum, auto

# Import will be done relatively when used as module
try:
    from .signal_system import Signal, SignalManager
except ImportError:
    # For standalone testing
    from signal_system import Signal, SignalManager

logger = logging.getLogger(__name__)

T = TypeVar("T", bound="HierarchicalObject")


class ObjectState(Enum):
    """Object lifecycle states."""

    CREATED = auto()
    INITIALIZED = auto()
    ACTIVE = auto()
    DESTROYING = auto()
    DESTROYED = auto()


class HierarchicalObject(ABC):
    """
    Base class providing Qt-like parent/child relationships and lifecycle management.

    This replaces Qt's QObject functionality with:
    - Parent/child hierarchy with automatic cleanup
    - Signal system integration
    - Property system
    - Event handling and propagation
    - Thread-safe operations
    - Memory leak prevention through weak references
    """

    def __init__(
        self, parent: Optional["HierarchicalObject"] = None, name: str = None
    ):
        self._parent_ref: Optional[weakref.ReferenceType] = None
        self._children: Set[weakref.ReferenceType] = set()
        self._children_by_name: Dict[str, weakref.ReferenceType] = {}  # O(1) name lookup cache
        self._child_storage: Dict[str, Any] = {}  # Strong references to prevent GC of children
        self._name = name or f"{self.__class__.__name__}_{id(self)}"
        self._lock = threading.RLock()
        self._signal_manager = SignalManager()
        self._state = ObjectState.CREATED
        self._properties: Dict[str, Any] = {}
        self._event_handlers: Dict[str, List[Callable]] = {}

        # Core signals that all objects have
        self.destroyed = self._signal_manager.create_signal("destroyed", type(None))
        self.parent_changed = self._signal_manager.create_signal(
            "parent_changed", "HierarchicalObject"
        )
        self.child_added = self._signal_manager.create_signal(
            "child_added", "HierarchicalObject"
        )
        self.child_removed = self._signal_manager.create_signal(
            "child_removed", "HierarchicalObject"
        )

        # Set parent relationship
        if parent is not None:
            self.set_parent(parent)

        self._state = ObjectState.INITIALIZED
        logger.debug(f"Created {self._name}")

    # NOTE: No 'name' property to avoid collision with CData 'name' attributes
    # Use objectName() to get the hierarchical name, or _name directly in internal code

    @property
    def state(self) -> ObjectState:
        """Current object state."""
        return self._state

    def parent(self) -> Optional["HierarchicalObject"]:
        """
        Get the parent object (None if no parent or parent was destroyed).

        This is a METHOD for backward compatibility with legacy ccp4i2 code that calls parent().
        It can be used both as a method: obj.parent() or as a property-like: obj.parent

        Returns:
            Parent object or None if no parent or parent was destroyed
        """
        if self._parent_ref is None:
            return None
        parent_obj = self._parent_ref()
        if parent_obj is None:
            # Parent was garbage collected
            self._parent_ref = None
        return parent_obj

    def get_parent(self) -> Optional["HierarchicalObject"]:
        """
        Get the parent object (backward compatibility wrapper).

        Returns:
            Parent object or None if no parent or parent was destroyed

        Note:
            This is a wrapper around the `parent()` method for backward compatibility.
            Prefer using the `parent()` method directly in new code.
        """
        return self.parent()

    def set_parent(self, parent: Optional["HierarchicalObject"]) -> bool:
        """
        Set the parent object.

        Args:
            parent: New parent object or None to remove parent

        Returns:
            True if parent was successfully set
        """
        with self._lock:
            if self._state == ObjectState.DESTROYED:
                logger.warning(
                    f"Cannot set parent on destroyed object {self._name}"
                )
                return False

            old_parent = self.parent()

            # Remove from old parent
            if old_parent is not None:
                old_parent._remove_child(self)

            # Set new parent
            if parent is not None:
                # Check if parent is destroyed (only if it's a HierarchicalObject)
                if hasattr(parent, '_state') and parent._state == ObjectState.DESTROYED:
                    logger.warning(f"Cannot set destroyed object as parent")
                    return False

                self._parent_ref = weakref.ref(parent)

                # Only add child if parent supports hierarchy
                if hasattr(parent, '_add_child'):
                    parent._add_child(self)
            else:
                self._parent_ref = None

            # Emit signal (guard against GC ordering issues)
            if hasattr(self, 'parent_changed') and self.parent_changed is not None:
                self.parent_changed.emit(parent)

            # Log parent change (safely handle non-HierarchicalObject parents)
            parent_name = None
            if parent is not None:
                if hasattr(parent, '_name'):
                    parent_name = parent._name
                elif hasattr(parent, '__class__'):
                    parent_name = parent.__class__.__name__
                else:
                    parent_name = str(type(parent))
            logger.debug(f"Set parent of {self._name} to {parent_name}")
            return True

    def _add_child(self, child: "HierarchicalObject"):
        """Internal method to add a child (called by set_parent)."""
        with self._lock:
            # Clean up any dead references first
            self._cleanup_dead_children()

            child_ref = weakref.ref(child)
            self._children.add(child_ref)

            # Add to name lookup cache for O(1) access
            child_name = child._name
            if child_name:
                self._children_by_name[child_name] = child_ref
                # Store strong reference to prevent GC
                self._child_storage[child_name] = child

            # Guard against GC ordering issues
            if hasattr(self, 'child_added') and self.child_added is not None:
                self.child_added.emit(child)
            logging.debug(f"Added child {child._name} to {self._name}")

    def _remove_child(self, child: "HierarchicalObject"):
        """Internal method to remove a child (called by set_parent)."""
        with self._lock:
            # Find and remove the child reference
            # Use 'is' for identity comparison to avoid calling __eq__ during GC
            to_remove = None
            for child_ref in self._children:
                if child_ref() is child:
                    to_remove = child_ref
                    break

            if to_remove:
                self._children.remove(to_remove)

                # Remove from name lookup cache and strong storage
                child_name = child._name
                if child_name and child_name in self._children_by_name:
                    # Only remove if it's the same child (names can be reused)
                    cached_ref = self._children_by_name.get(child_name)
                    if cached_ref is not None and cached_ref() is child:
                        del self._children_by_name[child_name]
                        # Also remove from strong storage
                        self._child_storage.pop(child_name, None)

                # Guard against GC ordering issues - signal might be cleaned up already
                if hasattr(self, 'child_removed') and self.child_removed is not None:
                    self.child_removed.emit(child)
                logger.debug(
                    f"Removed child {child._name} from {self._name}"
                )

    def _cleanup_dead_children(self):
        """Remove weak references to destroyed children."""
        dead_refs = {ref for ref in self._children if ref() is None}
        self._children -= dead_refs

        # Also clean up dead entries in name cache and strong storage
        dead_names = [name for name, ref in self._children_by_name.items() if ref() is None]
        for name in dead_names:
            del self._children_by_name[name]
            self._child_storage.pop(name, None)

    def children(self) -> List["HierarchicalObject"]:
        """Get list of all child objects."""
        with self._lock:
            self._cleanup_dead_children()
            result = [ref() for ref in self._children if ref() is not None]
            return result

    def find_child(
        self, name: str, recursive: bool = False
    ) -> Optional["HierarchicalObject"]:
        """Find a child by name. Uses O(1) cache lookup for direct children."""
        with self._lock:
            # O(1) lookup in name cache for direct children
            child_ref = self._children_by_name.get(name)
            if child_ref is not None:
                child = child_ref()
                if child is not None:
                    return child
                else:
                    # Dead reference - clean it up
                    del self._children_by_name[name]

            # Recursive search if requested
            if recursive:
                for child in self.children():
                    found = child.find_child(name, recursive=True)
                    if found:
                        return found

            return None

    def find_children(
        self, object_type: Type[T] = None, recursive: bool = False
    ) -> List[T]:
        """Find children by type."""
        results = []
        children = self.children()

        for child in children:
            if object_type is None or isinstance(child, object_type):
                results.append(child)

            if recursive:
                results.extend(child.find_children(object_type, recursive=True))

        return results

    def root(self) -> "HierarchicalObject":
        """Get the root object (topmost ancestor)."""
        current = self
        while current.parent is not None:
            current = current.parent
        return current

    def path_from_root(self) -> List[str]:
        """Get the path from root to this object as list of names.

        Only includes objects that have a valid name (non-empty objectName()).
        Objects without names are skipped to produce clean paths.
        """
        path = []
        current = self
        while current is not None:
            # Use objectName() if available (CData objects), else fall back to _name
            if hasattr(current, 'objectName') and callable(current.objectName):
                name = current.objectName()
                if name:  # Only include non-empty names
                    path.insert(0, name)
            elif hasattr(current, '_name') and current._name:
                path.insert(0, current._name)
            # Skip objects without names - don't add <ClassName> placeholders

            # Navigate to parent - handle both property and method parent access
            if hasattr(current, 'parent'):
                if callable(current.parent):
                    current = current.parent()
                else:
                    current = current.parent
            else:
                current = None
        return path

    def object_path(self) -> str:
        """Get the dot-separated path from root to this object.

        Returns a clean path like "task_name.container.inputData.XYZIN".
        For list items, produces paths like "list_name[0].child" (no dot before index).
        Only includes objects with valid names (non-empty objectName()).
        """
        parts = self.path_from_root()
        if not parts:
            return ""

        # Build path, omitting dot before array indices (names starting with '[')
        result = parts[0]
        for part in parts[1:]:
            if part.startswith('['):
                # Array index - don't add dot
                result += part
            else:
                # Regular name - add dot separator
                result += '.' + part
        return result

    # Legacy camelCase alias for compatibility with CCP4i2 code
    def objectPath(self) -> str:
        """Legacy alias for object_path(). Use object_path() in new code."""
        return self.object_path()

    def find_by_path(self, path: str) -> Optional["HierarchicalObject"]:
        """
        Find an object by dot-separated path with array indexing support.

        Supports:
        - Simple paths: "child.grandchild"
        - Array indexing: "children[0].property"
        - Mixed: "database.tables[2].columns[0]"

        Args:
            path: Dot-separated path string (e.g., "child.grandchild[0]")

        Returns:
            Found object or None
        """
        if not path:
            return self

        return self._resolve_path_segment(path.split("."), 0)

    def _resolve_path_segment(
        self, segments: List[str], index: int
    ) -> Optional["HierarchicalObject"]:
        """Recursively resolve path segments with array indexing."""
        if index >= len(segments):
            return self

        segment = segments[index]

        # Handle array indexing: "name[index]"
        if "[" in segment and segment.endswith("]"):
            name_part, bracket_part = segment.split("[", 1)
            array_index = bracket_part[:-1]  # Remove closing ']'

            try:
                idx = int(array_index)

                # Find child by name
                target_child = None
                if name_part:  # "child[0]"
                    target_child = self.find_child(name_part)
                else:  # "[0]" - use this object's children
                    children = self.children()
                    if 0 <= idx < len(children):
                        target_child = children[idx]

                if target_child is None:
                    return None

                # Handle array indexing on the found child
                if hasattr(target_child, "_indexed_children"):
                    # Custom array-like access
                    indexed_items = target_child._indexed_children()
                    if 0 <= idx < len(indexed_items):
                        next_obj = indexed_items[idx]
                    else:
                        return None
                elif name_part:  # Regular child with index - use the child itself
                    next_obj = target_child
                else:
                    return None

            except (ValueError, IndexError):
                return None

        else:
            # Simple name lookup
            next_obj = self.find_child(segment)

        if next_obj is None:
            return None

        # Continue with next segment
        return next_obj._resolve_path_segment(segments, index + 1)

    # Signal system access
    def create_signal(self, name: str, signal_type: type = None) -> Signal:
        """Create a new signal for this object."""
        return self._signal_manager.create_signal(name, signal_type)

    # Lifecycle management
    def destroy(self):
        """Destroy this object and all its children."""
        if self._state == ObjectState.DESTROYED:
            return

        logger.debug(f"Destroying {self._name}")
        self._state = ObjectState.DESTROYING

        # Destroy all children first
        for child in self.children():
            child.destroy()

        # Remove from parent (only if parent supports hierarchy)
        parent = self.parent()
        if parent and hasattr(parent, '_remove_child'):
            parent._remove_child(self)
            self._parent_ref = None

        # NOTE: destroyed signal is not emitted because nothing in the codebase
        # connects to it, and emitting signals during garbage collection can
        # cause timing issues with attribute deletion
        # destroyed_signal = getattr(self, 'destroyed', None)
        # if destroyed_signal is not None:
        #     try:
        #         destroyed_signal.emit()
        #     except Exception:
        #         pass

        # Cleanup
        self._signal_manager.cleanup()
        self._children.clear()
        self._children_by_name.clear()
        self._child_storage.clear()  # Clear strong references to children
        self._properties.clear()
        self._event_handlers.clear()

        self._state = ObjectState.DESTROYED
        logger.debug(f"Destroyed {self._name}")

    def connectSignal(self, origin, signal_name: str, handler):
        """
        Connect a signal from an origin object to a handler.

        This provides Qt-compatible API for connecting signals, with automatic
        adaptation between modern (dict) and legacy (int) handler signatures.

        Supports both modern and legacy handler signatures:
        - Modern: handler(data: dict) receives full payload
        - Legacy: handler(status: int) receives just an int value

        Args:
            origin: Object that emits the signal
            signal_name: Name of the signal to connect to (e.g., 'finished')
            handler: Callable to invoke when signal is emitted

        Example:
            sub_plugin = self.makePluginObject('mtzdump')
            self.connectSignal(sub_plugin, 'finished', self.on_mtzdump_finished)

        Raises:
            ValueError: If signal_name doesn't exist on origin object
        """
        # Check if origin has the requested signal
        if not hasattr(origin, signal_name):
            raise ValueError(
                f"Object {origin.__class__.__name__} has no signal '{signal_name}'"
            )

        signal = getattr(origin, signal_name)

        # Verify it's actually a Signal object
        if not hasattr(signal, 'connect'):
            raise ValueError(
                f"Attribute '{signal_name}' on {origin.__class__.__name__} "
                f"is not a Signal object"
            )

        # For 'finished' signal specifically, check for int vs dict signature
        if signal_name == 'finished':
            import inspect
            sig = inspect.signature(handler)
            params = list(sig.parameters.values())

            # Skip 'self' parameter if it's a method
            if params and params[0].name == 'self':
                params = params[1:]

            # Check parameter type hint if available
            expects_int = False
            if params:
                param = params[0]
                if param.annotation != inspect.Parameter.empty:
                    # Has type hint - check if it's int
                    expects_int = param.annotation == int or param.annotation == 'int'

            # If handler expects int, wrap it to extract from dict
            if expects_int:
                def adapter(status_dict):
                    """Adapter for legacy handlers expecting int."""
                    # Handle both dict and int for flexibility
                    if isinstance(status_dict, dict):
                        status = status_dict.get('finishStatus', 0)
                    else:
                        status = status_dict
                    return handler(status)

                signal.connect(adapter, weak=False)
                return

        # Default: connect handler directly
        signal.connect(handler, weak=False)

    def find(self, path: str):
        """
        Find a child object by name or path in the hierarchy.

        Supports both simple names and dotted paths:
        - find("XYZIN") - Searches for first child named "XYZIN" (depth-first)
        - find("protein.XYZIN") - Finds child "protein", then its child "XYZIN"
        - find("container.inputData.HKLIN") - Multi-level path navigation

        Uses O(1) cache lookup for immediate children when possible.

        Args:
            path: Object name or dotted path (e.g., "protein.XYZIN")

        Returns:
            The found object, or None if not found

        Examples:
            >>> container = CContainer(name="root")
            >>> protein = CContainer(name="protein", parent=container)
            >>> xyzin = CPdbDataFile(name="XYZIN", parent=protein)
            >>> container.find("protein.XYZIN")  # Returns xyzin
            >>> container.find("XYZIN")  # Also returns xyzin (depth-first search)
        """
        if not path:
            return None

        # Split path into components
        path_parts = path.split('.')

        # If single name, check immediate children first using O(1) cache
        if len(path_parts) == 1:
            name = path_parts[0]

            # O(1) lookup in name cache
            child_ref = self._children_by_name.get(name)
            if child_ref is not None:
                child = child_ref()
                if child is not None:
                    # Verify not destroyed
                    if hasattr(child, '_state') and child.state == ObjectState.DESTROYED:
                        del self._children_by_name[name]
                    else:
                        return child
                else:
                    # Dead reference - clean it up
                    del self._children_by_name[name]

            # Not found in cache - search recursively in children
            for child_ref in list(self._children):
                child = child_ref() if isinstance(child_ref, weakref.ReferenceType) else child_ref
                if child is not None and hasattr(child, 'find'):
                    result = child.find(name)
                    if result is not None:
                        return result

            return None

        # Multi-part path: navigate step by step using O(1) lookup
        first_name = path_parts[0]
        remaining_path = '.'.join(path_parts[1:])

        # O(1) lookup for first component
        child_ref = self._children_by_name.get(first_name)
        if child_ref is not None:
            child = child_ref()
            if child is not None:
                if hasattr(child, '_state') and child.state == ObjectState.DESTROYED:
                    del self._children_by_name[first_name]
                elif hasattr(child, 'find'):
                    return child.find(remaining_path)
            else:
                del self._children_by_name[first_name]

        return None

    def __del__(self):
        """Ensure cleanup on garbage collection."""
        # Use getattr with default to avoid AttributeError during GC
        # when object was only partially initialized
        try:
            state = getattr(self, '_state', ObjectState.DESTROYED)
            if state != ObjectState.DESTROYED:
                self.destroy()
        except Exception:
            # Suppress all exceptions during GC - object state may be inconsistent
            pass

    def __repr__(self) -> str:
        # Defensive: parent might not be a HierarchicalObject (e.g., QEventLoop)
        parent = self.parent()
        if parent is None:
            parent_name = "None"
        elif hasattr(parent, '_name'):
            parent_name = parent._name
        elif hasattr(parent, '__class__'):
            parent_name = f"<{parent.__class__.__name__}>"
        else:
            parent_name = str(type(parent))

        child_count = len(self.children())
        return (
            f"{self.__class__.__name__}(name={self._name}, "
            f"parent={parent_name}, children={child_count}, state={self._state.name})"
        )
