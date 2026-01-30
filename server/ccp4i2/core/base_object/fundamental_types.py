
"""Fundamental CCP4i2 data types that form the base of the type system."""

from typing import List, Any, Optional, TYPE_CHECKING
from ..base_object.base_classes import CData, ValueState
from .class_metadata import cdata_class, attribute, AttributeType

if TYPE_CHECKING:
    import xml.etree.ElementTree as ET


@cdata_class(
    error_codes={
        "101": {"description": "below minimum"},
        "102": {"description": "above maximum"},
        "103": {"description": "not one of limited allowed values"}
    },
    qualifiers={
        "max": None,
        "min": None,
        "enumerators": [],
        "menuText": [],
        "onlyEnumerators": False
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'
    ],
    qualifiers_definition={
        "default": {"type": int},
        "max": {"type": int, "description": "The inclusive minimum allowed value"},
        "min": {"type": int, "description": "The inclusive maximum allowed value"},
        "enumerators": {"type": list, "listItemType": "<class 'int'>", "description": "A Python list of allowed or recommended values - see onlyEnumerators"},
        "menuText": {"type": list, "listItemType": "<class 'str'>", "description": "A Python list of strings, matching items in enumerators list, to appear on GUI menu"},
        "onlyEnumerators": {"type": bool, "description": "If this is true then the enumerators are obligatory - otherwise they are treated as recommended values"}
    },
    gui_label="CInt",
)
class CInt(CData):
    """Integer value type."""

    def __hash__(self):
        """Make CInt hashable by object identity for use in sets and as dict keys."""
        return hash(id(self))

    def __delattr__(self, name):
        """Prevent deletion of _value attribute during garbage collection."""
        if name == '_value':
            # Don't allow deletion of _value to prevent AttributeError in __str__
            return
        super().__delattr__(name)

    def __init__(self, value: int = None, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)

        # Handle value setting with proper state tracking
        if value is None:
            # Default initialization - set value but mark as NOT_SET
            super().__setattr__("_value", 0)
            if hasattr(self, "_value_states"):
                self._value_states["value"] = ValueState.NOT_SET
        else:
            # Explicit value provided - mark as EXPLICITLY_SET
            self.value = value

    def _validate_value(self, val):
        """Validate value against min/max qualifiers."""
        # Skip validation if flag is set (used during .def.xml parsing)
        if getattr(self, '_skip_validation', False):
            return val

        min_val = self.get_qualifier("min")
        max_val = self.get_qualifier("max")

        # Convert min/max to appropriate numeric type if they're strings (from XML deserialization)
        if min_val is not None:
            try:
                min_val = float(min_val) if isinstance(self, CFloat) else int(min_val)
            except (ValueError, TypeError):
                min_val = None  # Invalid min value, skip validation

        if max_val is not None:
            try:
                max_val = float(max_val) if isinstance(self, CFloat) else int(max_val)
            except (ValueError, TypeError):
                max_val = None  # Invalid max value, skip validation

        if min_val is not None and val < min_val:
            raise ValueError(f"Value {val} is below minimum {min_val}")
        if max_val is not None and val > max_val:
            raise ValueError(f"Value {val} is above maximum {max_val}")

        return val

    @property
    def value(self):
        """Get the integer value."""
        return getattr(self, "_value", 0)

    @value.setter
    def value(self, val):
        """Set the integer value with validation."""
        validated = self._validate_value(int(val))
        super().__setattr__("_value", validated)
        if hasattr(self, "_value_states"):
            self._value_states["value"] = ValueState.EXPLICITLY_SET

    def __str__(self):
        return str(self.value)

    def __int__(self):
        return int(self.value)

    def set(self, value: int):
        """Set the value directly using .set() method.

        Args:
            value: The integer value to set, or None to clear/unset
        """
        if value is None:
            # Clear the value and mark as unset
            super().__setattr__("_value", 0)
            self.unSet()
            return self
        self.value = value
        return self

    def get(self, child_name=None):
        """Get the primitive int value for legacy API compatibility.

        Args:
            child_name: Ignored for fundamental types (provided for API compatibility)

        Returns:
            int: The primitive integer value
        """
        return self.value

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = True, allSet: bool = True) -> bool:
        """Check if the value has been set.

        Args:
            field_name: Optional field name. If not provided, checks if 'value' is set.
            allowUndefined: If True, allow None/undefined values to be considered "set"
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types (unused for fundamental types)

        Returns:
            True if the value (or specified field) has been set, False otherwise.
        """
        if field_name is None:
            field_name = "value"
        return super().isSet(field_name, allowUndefined=allowUndefined,
                           allowDefault=allowDefault, allSet=allSet)

    def _is_value_type(self) -> bool:
        return True

    # Arithmetic operators
    def __add__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value + other.value)
        elif isinstance(other, CFloat):
            return CFloat(self.value + other.value)
        elif isinstance(other, int):
            return CInt(self.value + other)
        elif isinstance(other, float):
            return CFloat(self.value + other)
        return int(self.value) + other

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value - other.value)
        elif isinstance(other, CFloat):
            return CFloat(self.value - other.value)
        elif isinstance(other, int):
            return CInt(self.value - other)
        elif isinstance(other, float):
            return CFloat(self.value - other)
        return int(self.value) - other

    def __rsub__(self, other):
        if isinstance(other, CInt):
            return CInt(other.value - self.value)
        elif isinstance(other, CFloat):
            return CFloat(other.value - self.value)
        elif isinstance(other, int):
            return CInt(other - self.value)
        elif isinstance(other, float):
            return CFloat(other - self.value)
        return other - int(self.value)

    def __mul__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value * other.value)
        elif isinstance(other, CFloat):
            return CFloat(self.value * other.value)
        elif isinstance(other, int):
            return CInt(self.value * other)
        elif isinstance(other, float):
            return CFloat(self.value * other)
        return int(self.value) * other

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, (CInt, int)):
            return CFloat(
                self.value / (other.value if isinstance(other, CInt) else other)
            )
        elif isinstance(other, (CFloat, float)):
            return CFloat(
                self.value / (other.value if isinstance(other, CFloat) else other)
            )
        return int(self.value) / other

    def __rtruediv__(self, other):
        if isinstance(other, (CInt, int)):
            return CFloat(
                (other.value if isinstance(other, CInt) else other) / self.value
            )
        elif isinstance(other, (CFloat, float)):
            return CFloat(
                (other.value if isinstance(other, CFloat) else other) / self.value
            )
        return other / int(self.value)

    def __floordiv__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value // other.value)
        elif isinstance(other, int):
            return CInt(self.value // other)
        elif isinstance(other, CFloat):
            return CFloat(self.value // other.value)
        elif isinstance(other, float):
            return CFloat(self.value // other)
        return int(self.value) // other

    def __rfloordiv__(self, other):
        if isinstance(other, CInt):
            return CInt(other.value // self.value)
        elif isinstance(other, int):
            return CInt(other // self.value)
        elif isinstance(other, CFloat):
            return CFloat(other.value // self.value)
        elif isinstance(other, float):
            return CFloat(other // self.value)
        return other // int(self.value)

    def __mod__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value % other.value)
        elif isinstance(other, int):
            return CInt(self.value % other)
        elif isinstance(other, CFloat):
            return CFloat(self.value % other.value)
        elif isinstance(other, float):
            return CFloat(self.value % other)
        return int(self.value) % other

    def __rmod__(self, other):
        if isinstance(other, CInt):
            return CInt(other.value % self.value)
        elif isinstance(other, int):
            return CInt(other % self.value)
        elif isinstance(other, CFloat):
            return CFloat(other.value % self.value)
        elif isinstance(other, float):
            return CFloat(other % self.value)
        return other % int(self.value)

    def __pow__(self, other):
        if isinstance(other, CInt):
            return CInt(self.value**other.value)
        elif isinstance(other, int):
            return CInt(self.value**other)
        elif isinstance(other, CFloat):
            return CFloat(self.value**other.value)
        elif isinstance(other, float):
            return CFloat(self.value**other)
        return int(self.value) ** other

    def __rpow__(self, other):
        if isinstance(other, CInt):
            return CInt(other.value**self.value)
        elif isinstance(other, int):
            return CInt(other**self.value)
        elif isinstance(other, CFloat):
            return CFloat(other.value**self.value)
        elif isinstance(other, float):
            return CFloat(other**self.value)
        return other ** int(self.value)

    # Comparison operators
    def __eq__(self, other):
        if isinstance(other, CInt):
            return self.value == other.value
        elif isinstance(other, CFloat):
            return float(self.value) == other.value
        return int(self.value) == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, CInt):
            return self.value < other.value
        elif isinstance(other, CFloat):
            return float(self.value) < other.value
        return int(self.value) < other

    def __le__(self, other):
        if isinstance(other, CInt):
            return self.value <= other.value
        elif isinstance(other, CFloat):
            return float(self.value) <= other.value
        return int(self.value) <= other

    def __gt__(self, other):
        if isinstance(other, CInt):
            return self.value > other.value
        elif isinstance(other, CFloat):
            return float(self.value) > other.value
        return int(self.value) > other

    def __ge__(self, other):
        if isinstance(other, CInt):
            return self.value >= other.value
        elif isinstance(other, CFloat):
            return float(self.value) >= other.value
        return int(self.value) >= other

    def validity(self):
        """Validate the integer value against qualifiers.

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from .error_reporting import (CErrorReport, SEVERITY_ERROR,
                                      SEVERITY_WARNING)

        report = CErrorReport()

        # Check allowUndefined - if False and value not set, it's an error
        allow_undefined = self.get_qualifier('allowUndefined')
        if allow_undefined is False and not self.isSet():
            obj_path = self.object_path() if hasattr(self, 'object_path') else self.objectName()
            report.append(
                "CInt", 100, f"Required value not set: {obj_path}",
                obj_path, SEVERITY_ERROR
            )
            return report  # No point checking other constraints if not set

        # If optional (allowUndefined is True or None) and not explicitly set,
        # skip min/max validation. The script will check isSet() and not write
        # anything for this parameter if it's not set.
        is_optional = allow_undefined is None or allow_undefined is True
        if is_optional and not self.isSet():
            return report

        val = self.value

        # Check min/max constraints
        min_val = self.get_qualifier("min")
        if min_val is not None and val < min_val:
            report.append(
                "CInt", 101, f"Value {val} is below minimum {min_val}",
                self.objectName(), SEVERITY_ERROR
            )

        max_val = self.get_qualifier("max")
        if max_val is not None and val > max_val:
            report.append(
                "CInt", 102, f"Value {val} is above maximum {max_val}",
                self.objectName(), SEVERITY_ERROR
            )

        # Check enumerators if onlyEnumerators is True
        only_enumerators = self.get_qualifier("onlyEnumerators")
        enumerators = self.get_qualifier("enumerators")
        if only_enumerators and enumerators:
            # Convert string enumerators to int for proper comparison
            # (metadata often stores enumerators as strings like ['0', '1', '2'])
            numeric_enumerators = []
            for e in enumerators:
                try:
                    numeric_enumerators.append(int(e) if isinstance(e, str) else e)
                except (ValueError, TypeError):
                    numeric_enumerators.append(e)  # Keep original if not convertible
            if val not in numeric_enumerators:
                report.append(
                    "CInt", 103,
                    f"Value {val} not in allowed values {enumerators}",
                    self.objectName(), SEVERITY_ERROR
                )

        return report


@cdata_class(
    error_codes={
        "101": {"description": "below minimum"},
        "102": {"description": "above maximum"},
        "103": {"description": "not one of limited allowed values"}
    },
    qualifiers={
        "max": None,
        "min": None,
        "enumerators": [],
        "menuText": [],
        "onlyEnumerators": False
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'
    ],
    qualifiers_definition={
        "default": {"type": float},
        "max": {"description": "The inclusive maximum value"},
        "min": {"description": "The inclusive minimum value"},
        "enumerators": {"type": list, "description": "A Python list of allowed or recommended values - see onlyEnumerators"},
        "menuText": {"type": list, "listItemType": "<class 'str'>", "description": "A Python list of strings, matching items in enumerators list, to appear on GUI menu"},
        "onlyEnumerators": {"type": bool, "description": "If this is true then the enumerators are obligatory - otherwise they are treated as recommended values"}
    },
    gui_label="CFloat",
)
class CFloat(CData):
    """Float value type."""

    def __hash__(self):
        """Make CFloat hashable by object identity for use in sets and as dict keys."""
        return hash(id(self))

    def __delattr__(self, name):
        """Prevent deletion of _value attribute during garbage collection."""
        if name == '_value':
            # Don't allow deletion of _value to prevent AttributeError in __str__
            return
        super().__delattr__(name)

    def __init__(self, value: float = None, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)

        # Handle value setting with proper state tracking
        if value is None:
            # Default initialization - set value but mark as NOT_SET
            super().__setattr__("_value", 0.0)
            if hasattr(self, "_value_states"):
                self._value_states["value"] = ValueState.NOT_SET
        else:
            # Explicit value provided - mark as EXPLICITLY_SET
            self.value = value

    def _validate_value(self, val):
        """Validate value against min/max qualifiers."""
        # Skip validation if flag is set (used during .def.xml parsing)
        if getattr(self, '_skip_validation', False):
            return val

        min_val = self.get_qualifier("min")
        max_val = self.get_qualifier("max")

        # Convert min/max to appropriate numeric type if they're strings (from XML deserialization)
        if min_val is not None:
            try:
                min_val = float(min_val) if isinstance(self, CFloat) else int(min_val)
            except (ValueError, TypeError):
                min_val = None  # Invalid min value, skip validation

        if max_val is not None:
            try:
                max_val = float(max_val) if isinstance(self, CFloat) else int(max_val)
            except (ValueError, TypeError):
                max_val = None  # Invalid max value, skip validation

        if min_val is not None and val < min_val:
            raise ValueError(f"Value {val} is below minimum {min_val}")
        if max_val is not None and val > max_val:
            raise ValueError(f"Value {val} is above maximum {max_val}")

        return val

    @property
    def value(self):
        """Get the float value."""
        return getattr(self, "_value", 0.0)

    @value.setter
    def value(self, val):
        """Set the float value with validation."""
        validated = self._validate_value(float(val))
        old_value = getattr(self, "_value", None)
        super().__setattr__("_value", validated)
        if hasattr(self, "_value_states"):
            # Only mark as EXPLICITLY_SET if this is a real value change.
            # If setting to the same value while currently NOT_SET, keep it NOT_SET.
            # This prevents spurious state changes during .def.xml loading and merging.
            current_state = self._value_states.get("value", ValueState.NOT_SET)
            if current_state == ValueState.NOT_SET and old_value == validated:
                # Keep as NOT_SET - this is internal copying, not user assignment
                pass
            else:
                self._value_states["value"] = ValueState.EXPLICITLY_SET

    def __str__(self):
        return str(self.value)

    def __float__(self):
        return float(self.value)

    def __abs__(self):
        """Return absolute value of the float."""
        return abs(self.value)

    def __sub__(self, other):
        """Support subtraction for cell difference calculations."""
        if isinstance(other, (CFloat, CInt)):
            return self.value - other.value
        return self.value - other

    def __rsub__(self, other):
        """Support reverse subtraction."""
        if isinstance(other, (CFloat, CInt)):
            return other.value - self.value
        return other - self.value

    def set(self, value: float):
        """Set the value directly using .set() method.

        Args:
            value: The float value to set, or None to clear/unset
        """
        if value is None:
            # Clear the value and mark as unset
            super().__setattr__("_value", 0.0)
            self.unSet()
            return self
        self.value = value
        return self

    def get(self, child_name=None):
        """Get the primitive float value for legacy API compatibility.

        Args:
            child_name: Ignored for fundamental types (provided for API compatibility)

        Returns:
            float: The primitive float value
        """
        return self.value

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = True, allSet: bool = True) -> bool:
        """Check if the value has been set.

        Args:
            field_name: Optional field name. If not provided, checks if 'value' is set.
            allowUndefined: If True, allow None/undefined values to be considered "set"
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types (unused for fundamental types)

        Returns:
            True if the value (or specified field) has been set, False otherwise.
        """
        if field_name is None:
            field_name = "value"
        return super().isSet(field_name, allowUndefined=allowUndefined,
                           allowDefault=allowDefault, allSet=allSet)

    def _is_value_type(self) -> bool:
        return True


    # Arithmetic operators
    def __add__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value + other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value + other.value)
        elif isinstance(other, float):
            return CFloat(self.value + other)
        elif isinstance(other, int):
            return CFloat(self.value + other)
        return float(self.value) + other

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value - other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value - other.value)
        elif isinstance(other, float):
            return CFloat(self.value - other)
        elif isinstance(other, int):
            return CFloat(self.value - other)
        return float(self.value) - other

    def __rsub__(self, other):
        if isinstance(other, CFloat):
            return CFloat(other.value - self.value)
        elif isinstance(other, CInt):
            return CFloat(other.value - self.value)
        elif isinstance(other, float):
            return CFloat(other - self.value)
        elif isinstance(other, int):
            return CFloat(other - self.value)
        return other - float(self.value)

    def __mul__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value * other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value * other.value)
        elif isinstance(other, float):
            return CFloat(self.value * other)
        elif isinstance(other, int):
            return CFloat(self.value * other)
        return float(self.value) * other

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value / other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value / other.value)
        elif isinstance(other, float):
            return CFloat(self.value / other)
        elif isinstance(other, int):
            return CFloat(self.value / other)
        return float(self.value) / other

    def __rtruediv__(self, other):
        if isinstance(other, CFloat):
            return CFloat(other.value / self.value)
        elif isinstance(other, CInt):
            return CFloat(other.value / self.value)
        elif isinstance(other, float):
            return CFloat(other / self.value)
        elif isinstance(other, int):
            return CFloat(other / self.value)
        return other / float(self.value)

    def __floordiv__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value // other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value // other.value)
        elif isinstance(other, float):
            return CFloat(self.value // other)
        elif isinstance(other, int):
            return CFloat(self.value // other)
        return float(self.value) // other

    def __rfloordiv__(self, other):
        if isinstance(other, CFloat):
            return CFloat(other.value // self.value)
        elif isinstance(other, CInt):
            return CFloat(other.value // self.value)
        elif isinstance(other, float):
            return CFloat(other // self.value)
        elif isinstance(other, int):
            return CFloat(other // self.value)
        return other // float(self.value)

    def __mod__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value % other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value % other.value)
        elif isinstance(other, float):
            return CFloat(self.value % other)
        elif isinstance(other, int):
            return CFloat(self.value % other)
        return float(self.value) % other

    def __rmod__(self, other):
        if isinstance(other, CFloat):
            return CFloat(other.value % self.value)
        elif isinstance(other, CInt):
            return CFloat(other.value % self.value)
        elif isinstance(other, float):
            return CFloat(other % self.value)
        elif isinstance(other, int):
            return CFloat(other % self.value)
        return other % float(self.value)

    def __pow__(self, other):
        if isinstance(other, CFloat):
            return CFloat(self.value**other.value)
        elif isinstance(other, CInt):
            return CFloat(self.value**other.value)
        elif isinstance(other, float):
            return CFloat(self.value**other)
        elif isinstance(other, int):
            return CFloat(self.value**other)
        return float(self.value) ** other

    def __rpow__(self, other):
        if isinstance(other, CFloat):
            return CFloat(other.value**self.value)
        elif isinstance(other, CInt):
            return CFloat(other.value**self.value)
        elif isinstance(other, float):
            return CFloat(other**self.value)
        elif isinstance(other, int):
            return CFloat(other**self.value)
        return other ** float(self.value)

    # Comparison operators
    def __eq__(self, other):
        if isinstance(other, CFloat):
            return self.value == other.value
        elif isinstance(other, CInt):
            return self.value == float(other.value)
        return float(self.value) == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, CFloat):
            return self.value < other.value
        elif isinstance(other, CInt):
            return self.value < float(other.value)
        return float(self.value) < other

    def __le__(self, other):
        if isinstance(other, CFloat):
            return self.value <= other.value
        elif isinstance(other, CInt):
            return self.value <= float(other.value)
        return float(self.value) <= other

    def __gt__(self, other):
        if isinstance(other, CFloat):
            return self.value > other.value
        elif isinstance(other, CInt):
            return self.value > float(other.value)
        return float(self.value) > other

    def __ge__(self, other):
        if isinstance(other, CFloat):
            return self.value >= other.value
        elif isinstance(other, CInt):
            return self.value >= float(other.value)
        return float(self.value) >= other

    def validity(self):
        """Validate the float value against qualifiers.

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from .error_reporting import (CErrorReport, SEVERITY_ERROR,
                                      SEVERITY_WARNING)

        report = CErrorReport()

        # Check allowUndefined - if False and value not set, it's an error
        allow_undefined = self.get_qualifier('allowUndefined')
        if allow_undefined is False and not self.isSet():
            obj_path = self.object_path() if hasattr(self, 'object_path') else self.objectName()
            report.append(
                "CFloat", 100, f"Required value not set: {obj_path}",
                obj_path, SEVERITY_ERROR
            )
            return report  # No point checking other constraints if not set

        # If optional (allowUndefined is True or None) and not explicitly set,
        # skip min/max validation. The script will check isSet() and not write
        # anything for this parameter if it's not set.
        is_optional = allow_undefined is None or allow_undefined is True
        if is_optional and not self.isSet():
            return report

        val = self.value

        # Check min/max constraints
        min_val = self.get_qualifier("min")
        if min_val is not None and val < min_val:
            report.append(
                "CFloat", 101, f"Value {val} is below minimum {min_val}",
                self.objectName(), SEVERITY_ERROR
            )

        max_val = self.get_qualifier("max")
        if max_val is not None and val > max_val:
            report.append(
                "CFloat", 102, f"Value {val} is above maximum {max_val}",
                self.objectName(), SEVERITY_ERROR
            )

        # Check enumerators if onlyEnumerators is True
        only_enumerators = self.get_qualifier("onlyEnumerators")
        enumerators = self.get_qualifier("enumerators")
        if only_enumerators and enumerators:
            # Convert string enumerators to float for proper comparison
            # (metadata often stores enumerators as strings like ['0.0', '1.0', '2.0'])
            numeric_enumerators = []
            for e in enumerators:
                try:
                    numeric_enumerators.append(float(e) if isinstance(e, str) else e)
                except (ValueError, TypeError):
                    numeric_enumerators.append(e)  # Keep original if not convertible
            if val not in numeric_enumerators:
                report.append(
                    "CFloat", 103,
                    f"Value {val} not in allowed values {enumerators}",
                    self.objectName(), SEVERITY_ERROR
                )

        return report


@cdata_class(
    error_codes={
        "101": {"description": "String too short"},
        "102": {"description": "String too long"},
        "103": {"description": "not one of limited allowed values"},
        "104": {"description": "Contains disallowed characters"},
        "105": {"description": "Value does not match required pattern"}
    },
    qualifiers={
        "minLength": None,
        "maxLength": None,
        "enumerators": [],
        "menuText": [],
        "onlyEnumerators": False,
        "charWidth": -1,
        "allowedCharsCode": 0
    },
    qualifiers_order=[
        'minLength',
        'maxLength',
        'onlyEnumerators',
        'enumerators',
        'menuText',
        'allowedCharsCode'
    ],
    qualifiers_definition={
        "default": {"type": str},
        "maxLength": {"type": int, "description": "Maximum length of string"},
        "minLength": {"type": int, "description": "Minimum length of string"},
        "enumerators": {"type": list, "description": "A list of allowed or recommended values for string"},
        "menuText": {"type": list, "description": "A list of strings equivalent to the enumerators that will appear in the GUI"},
        "onlyEnumerators": {"type": bool, "description": "If this is true then the enumerators are obligatory - otherwise they are treated as recommended values"},
        "charWidth": {"type": int, "description": "Width of the string in characters (for GUI layout)"},
        "allowedCharsCode": {"type": int, "description": "Flag if the text is limited to set of allowed characters"},
        "patternRegex": {"type": str, "description": "Regular expression pattern that the value must match"},
        "patternErrorMessage": {"type": str, "description": "Custom error message when patternRegex validation fails"}
    },
    gui_label="CString",
)
class CString(CData):
    """String value type with Python string dunder methods."""
    def __init__(self, value: str = "", parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        self.value = value

    def __hash__(self):
        """Make CString hashable by object identity for use in sets and as dict keys.

        Identity-based hashing is required because:
        1. CString objects are mutable (value can change)
        2. Multiple CString objects can have the same value
        3. Children tracking uses sets which require stable hashing

        NOTE: CString objects cannot be used directly as dictionary keys to match
        plain string keys. Plugin code must convert to string using str() or .value:

        WRONG:  dict[cstring_obj]
        RIGHT:  dict[str(cstring_obj)]
        """
        return hash(id(self))

    def __delattr__(self, name):
        """Prevent deletion of _value attribute during garbage collection."""
        if name == '_value':
            # Don't allow deletion of _value to prevent AttributeError in __str__
            return
        super().__delattr__(name)

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return repr(self.value)

    def __eq__(self, other):
        if isinstance(other, CString):
            return self.value == other.value
        return self.value == other

    def __ne__(self, other):
        if isinstance(other, CString):
            return self.value != other.value
        return self.value != other

    def __lt__(self, other):
        if isinstance(other, CString):
            return self.value < other.value
        return self.value < other

    def __le__(self, other):
        if isinstance(other, CString):
            return self.value <= other.value
        return self.value <= other

    def __gt__(self, other):
        if isinstance(other, CString):
            return self.value > other.value
        return self.value > other

    def __ge__(self, other):
        if isinstance(other, CString):
            return self.value >= other.value
        return self.value >= other

    def __add__(self, other):
        if isinstance(other, CString):
            return CString(self.value + other.value)
        return CString(self.value + str(other))

    def __radd__(self, other):
        if isinstance(other, CString):
            return CString(other.value + self.value)
        return CString(str(other) + self.value)

    def __getitem__(self, key):
        return self.value[key]

    def __contains__(self, item):
        return item in self.value

    def __len__(self):
        if self.value is None:
            return 0
        return len(self.value)

    @property
    def value(self):
        """Get the string value."""
        return getattr(self, "_value", "")

    @value.setter
    def value(self, val):
        """Set the string value with state tracking."""
        from .base_classes import ValueState

        validated = str(val) if val is not None else ""
        super().__setattr__("_value", validated)
        if hasattr(self, "_value_states"):
            self._value_states["value"] = ValueState.EXPLICITLY_SET

    def set(self, value: str):
        """Set the value directly using .set() method.

        Args:
            value: The string value to set, or None to clear/unset
        """
        if value is None:
            # Clear the value and mark as unset
            super().__setattr__("_value", "")
            self.unSet()
            return self
        self.value = value
        return self

    def get(self) -> str:
        """Get the value as a plain string.

        This method provides Qt API compatibility for legacy code
        that expects property objects with .get() method.

        Returns:
            The string value
        """
        return self.value

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = True, allSet: bool = True) -> bool:
        """Check if the value has been set.

        For CString, empty strings are treated as "not set" to support legacy
        plugin code that uses empty string to mean "auto-detect" (e.g., x2mtz
        column selection).

        Args:
            field_name: Optional field name. If not provided, checks if 'value' is set.
            allowUndefined: If True, allow None/undefined values to be considered "set"
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types (unused for fundamental types)

        Returns:
            True if the value (or specified field) has been set to a non-empty string,
            False otherwise (including empty string case).
        """
        if field_name is None:
            field_name = "value"

        # First check parent's isSet logic (state tracking, etc.)
        parent_is_set = super().isSet(field_name, allowUndefined=allowUndefined,
                                     allowDefault=allowDefault, allSet=allSet)

        # If parent says it's not set, respect that
        if not parent_is_set:
            return False

        # If parent says it IS set, also check if the value is empty string
        # Empty strings are treated as "not set" for CString
        if field_name == "value" and self.value == "":
            return False

        return True

    def _is_value_type(self) -> bool:
        return True

    def validity(self):
        """Validate the string value against qualifiers.

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from .error_reporting import (CErrorReport, SEVERITY_ERROR,
                                      SEVERITY_WARNING)

        report = CErrorReport()

        # Check allowUndefined - if False and value not set, it's an error
        allow_undefined = self.get_qualifier('allowUndefined')
        if allow_undefined is False and not self.isSet():
            obj_path = self.object_path() if hasattr(self, 'object_path') else self.objectName()
            report.append(
                "CString", 100, f"Required value not set: {obj_path}",
                obj_path, SEVERITY_ERROR
            )
            return report  # No point checking other constraints if not set

        # If optional (allowUndefined is True or None) and not explicitly set,
        # skip minLength/maxLength validation. The script will check isSet()
        # and not write anything for this parameter if it's not set.
        is_optional = allow_undefined is None or allow_undefined is True
        if is_optional and not self.isSet():
            return report

        val = self.value

        # Check minLength constraint
        min_length = self.get_qualifier("minLength")
        if min_length is not None and len(val) < min_length:
            report.append(
                "CString", 101,
                f"String length {len(val)} below minimum {min_length}",
                self.objectName(), SEVERITY_ERROR
            )

        # Check maxLength constraint
        max_length = self.get_qualifier("maxLength")
        if max_length is not None and len(val) > max_length:
            report.append(
                "CString", 102,
                f"String length {len(val)} above maximum {max_length}",
                self.objectName(), SEVERITY_ERROR
            )

        # Check enumerators if onlyEnumerators is True
        only_enumerators = self.get_qualifier("onlyEnumerators")
        enumerators = self.get_qualifier("enumerators")
        if only_enumerators and enumerators:
            if val not in enumerators:
                report.append(
                    "CString", 103,
                    f"Value '{val}' not in allowed values {enumerators}",
                    self.objectName(), SEVERITY_ERROR
                )

        # Check patternRegex if set
        pattern_regex = self.get_qualifier("patternRegex")
        if pattern_regex and val:
            import re
            if not re.fullmatch(pattern_regex, val):
                # Use custom error message if provided, otherwise generic
                pattern_error = self.get_qualifier("patternErrorMessage")
                if not pattern_error:
                    pattern_error = f"Value does not match required pattern"
                report.append(
                    self.__class__.__name__, 105,
                    pattern_error,
                    self.objectName(), SEVERITY_ERROR
                )

        return report

    # Class variable for whitespace pattern (cached)
    RE_PATTERN_WHITESPACE = None

    def reWhiteSpacePattern(self):
        """
        Get compiled regex pattern for whitespace characters.

        Legacy API compatibility method for CRangeSelection and other string validators.
        Uses cached class variable to avoid recompiling the pattern.

        Returns:
            re.Pattern: Compiled regex pattern matching all whitespace characters
        """
        if CString.RE_PATTERN_WHITESPACE is None:
            import string
            import re
            pat = ''
            for item in string.whitespace:
                pat = pat + repr(item)[1:-1] + '|'
            CString.RE_PATTERN_WHITESPACE = re.compile(pat[0:-1])
        return CString.RE_PATTERN_WHITESPACE

    def removeWhiteSpace(self, arg):
        """
        Remove all whitespace characters from a string.

        Legacy API compatibility method for CRangeSelection and other validators.

        Args:
            arg: String to remove whitespace from

        Returns:
            str: String with all whitespace removed
        """
        p = self.reWhiteSpacePattern()
        arg = p.sub('', arg)
        return arg

    def split(self, sep=' ', maxsplit=-1):
        """
        Split string into list of substrings.

        Legacy API compatibility method that wraps Python's str.split().
        Used by CRangeSelection and other string-based classes.

        Args:
            sep: Separator to split on (default: space)
            maxsplit: Maximum number of splits (default: -1 = no limit)

        Returns:
            list: List of substrings
        """
        return self.value.split(sep, maxsplit)

    def splitlines(self, keepends=False):
        """
        Split string into list of lines.

        Legacy API compatibility method that wraps Python's str.splitlines().

        Args:
            keepends: If True, line breaks are included in results (default: False)

        Returns:
            list: List of lines
        """
        return self.value.splitlines(keepends)

    def startswith(self, prefix, start=None, end=None):
        """
        Check if string starts with prefix.

        Delegates to underlying string value's startswith method.

        Args:
            prefix: String prefix to check (or tuple of prefixes)
            start: Optional start index (default: None)
            end: Optional end index (default: None)

        Returns:
            bool: True if string starts with prefix, False otherwise
        """
        if self.value is None:
            return False
        return str(self.value).startswith(prefix, start, end)

    def endswith(self, suffix, start=None, end=None):
        """
        Check if string ends with suffix.

        Delegates to underlying string value's endswith method.

        Args:
            suffix: String suffix to check (or tuple of suffixes)
            start: Optional start index (default: None)
            end: Optional end index (default: None)

        Returns:
            bool: True if string ends with suffix, False otherwise
        """
        if self.value is None:
            return False
        return str(self.value).endswith(suffix, start, end)

    def upper(self):
        """Convert string to uppercase."""
        return CString(self.value.upper() if self.value else "")

    def lower(self):
        """Convert string to lowercase."""
        return CString(self.value.lower() if self.value else "")

    def strip(self, chars=None):
        """Remove leading and trailing characters."""
        return CString(self.value.strip(chars) if self.value else "")

    def lstrip(self, chars=None):
        """Remove leading characters."""
        return CString(self.value.lstrip(chars) if self.value else "")

    def rstrip(self, chars=None):
        """Remove trailing characters."""
        return CString(self.value.rstrip(chars) if self.value else "")

    def replace(self, old, new, count=-1):
        """Replace occurrences of substring."""
        return CString(self.value.replace(old, new, count) if self.value else "")

    def find(self, sub, start=None, end=None):
        """Find substring, return index or -1."""
        if self.value is None:
            return -1
        return self.value.find(sub, start, end)

    def rfind(self, sub, start=None, end=None):
        """Find substring from right, return index or -1."""
        if self.value is None:
            return -1
        return self.value.rfind(sub, start, end)

    def index(self, sub, start=None, end=None):
        """Find substring, raise ValueError if not found."""
        if self.value is None:
            raise ValueError("substring not found")
        return self.value.index(sub, start, end)

    def rindex(self, sub, start=None, end=None):
        """Find substring from right, raise ValueError if not found."""
        if self.value is None:
            raise ValueError("substring not found")
        return self.value.rindex(sub, start, end)


@cdata_class(
    error_codes={
        "101": {"description": "not allowed value"}
    },
    qualifiers={
        "menuText": ['NotImplemented', 'NotImplemented']
    },
    qualifiers_order=[
        'charWidth'
    ],
    qualifiers_definition={
        "default": {"type": bool},
        "menuText": {"type": list, "listItemType": "<class 'str'>", "description": "A list of two string descriptions for true and false"}
    },
    gui_label="CBoolean",
)
class CBoolean(CData):
    """Boolean value type."""

    def __init__(self, value: bool = None, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)

        # Handle value setting with proper state tracking
        if value is None:
            # Default initialization - set value but mark as NOT_SET
            super().__setattr__("_value", False)
            if hasattr(self, "_value_states"):
                self._value_states["value"] = ValueState.NOT_SET
        else:
            # Explicit value provided - mark as EXPLICITLY_SET
            self.value = value

    @property
    def value(self):
        """Get the boolean value."""
        return getattr(self, "_value", False)

    @value.setter
    def value(self, val):
        """Set the boolean value with state tracking."""
        # Convert string representations to bool
        if isinstance(val, str):
            val = val.lower() in ('true', '1', 'yes')
        validated = bool(val)
        old_value = getattr(self, "_value", None)
        super().__setattr__("_value", validated)
        if hasattr(self, "_value_states"):
            # Only mark as EXPLICITLY_SET if this is a real value change.
            # If setting to the same value while currently NOT_SET, keep it NOT_SET.
            # This prevents spurious state changes during .def.xml loading and merging.
            current_state = self._value_states.get("value", ValueState.NOT_SET)
            if current_state == ValueState.NOT_SET and old_value == validated:
                # Keep as NOT_SET - this is internal copying, not user assignment
                pass
            else:
                self._value_states["value"] = ValueState.EXPLICITLY_SET

    def __hash__(self):
        """Make CBoolean hashable for use in sets and as dict keys."""
        return hash(id(self))

    def __delattr__(self, name):
        """Prevent deletion of _value attribute during garbage collection."""
        if name == '_value':
            # Don't allow deletion of _value to prevent AttributeError in __str__
            return
        super().__delattr__(name)

    def __str__(self):
        return str(self.value)

    def __bool__(self):
        """Return True if this boolean parameter is set AND has a True value.

        For CBoolean, we check both:
        1. Is the parameter set? (explicitly or via default)
        2. Is the value True?

        This matches the common pattern where boolean flags should be both
        present and True to trigger behavior:
            if self.container.controlParameters.USE_FEATURE:
                # Feature is configured AND enabled

        Returns False if:
        - The parameter is not set (no value configured)
        - The parameter is set to False (explicitly disabled)

        Returns True if:
        - The parameter is set (explicitly or via default) to True

        To check only if set (regardless of value): use .isSet(allowDefault=True)
        To check only the value (regardless of set state): use bool(param.value)
        """
        return self.isSet(allowDefault=True) and bool(self.value)

    def set(self, value: bool):
        """Set the value directly using .set() method.

        Args:
            value: The boolean value to set, or None to clear/unset
        """
        if value is None:
            # Clear the value and mark as unset
            super().__setattr__("_value", False)
            self.unSet()
            return self
        self.value = value
        return self

    def get(self, child_name=None):
        """Get the primitive boolean value for legacy API compatibility.

        Legacy wrapper code expects .get() to return the primitive value, not a dict.
        Example: int(cpar.AWA_NCS_RESTRAINTS.get()) expects a boolean that can be converted to int.

        Args:
            child_name: Ignored for fundamental types (provided for API compatibility)

        Returns:
            bool: The primitive boolean value
        """
        return self.value

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = True, allSet: bool = True) -> bool:
        """Check if the value has been set.

        Args:
            field_name: Optional field name. If not provided, checks if 'value' is set.
            allowUndefined: If True, allow None/undefined values to be considered "set"
            allowDefault: If False, consider values that equal the default as "not set"
            allSet: For container types (unused for fundamental types)

        Returns:
            True if the value (or specified field) has been set, False otherwise.
        """
        if field_name is None:
            field_name = "value"
        return super().isSet(field_name, allowUndefined=allowUndefined,
                           allowDefault=allowDefault, allSet=allSet)

    def _is_value_type(self) -> bool:
        return True

    # Boolean and comparison operators
    def __eq__(self, other):
        return bool(self.value) == other

    def __ne__(self, other):
        return bool(self.value) != other

    def __and__(self, other):
        return bool(self.value) and other

    def __rand__(self, other):
        return other and bool(self.value)

    def __or__(self, other):
        return bool(self.value) or other

    def __ror__(self, other):
        return other or bool(self.value)

    def __invert__(self):
        return not bool(self.value)

    def validity(self):
        """Validate the boolean value.

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from .error_reporting import CErrorReport, SEVERITY_ERROR

        report = CErrorReport()

        # Check allowUndefined - if False and value not set, it's an error
        allow_undefined = self.get_qualifier('allowUndefined')
        if allow_undefined is False and not self.isSet():
            obj_path = self.object_path() if hasattr(self, 'object_path') else self.objectName()
            report.append(
                "CBoolean", 100, f"Required value not set: {obj_path}",
                obj_path, SEVERITY_ERROR
            )

        return report


@cdata_class(
    error_codes={
        "101": {"description": "List shorter than required minimum length"},
        "102": {"description": "List longer than required maximum length"},
        "103": {"description": "Consecutive values in list fail comparison test"},
        "104": {"description": "Attempting to add object of wrong type"},
        "105": {"description": "Attempting to add object of correct type but wrong qualifiers"},
        "106": {"description": "Attempting to add data which does not satisfy the qualifiers for a list item"},
        "107": {"description": "Deleting item will reduce list below minimum length"},
        "108": {"description": "Adding item will extend list beyond maximum length"},
        "109": {"description": "Invalid item class"},
        "110": {"description": "etree (XML) list item of wrong type"},
        "112": {"description": "No list item object set for list"}
    },
    qualifiers={
        "listMinLength": 0
    },
    qualifiers_order=[
        'listMinLength',
        'listMaxLength',
        'listCompare'
    ],
    qualifiers_definition={
        "default": {"type": list},
        "listMaxLength": {"type": int, "description": "Inclusive maximum length of list"},
        "listMinLength": {"type": int, "description": "Inclusive minimum length of list"},
        "listCompare": {"type": int, "description": "If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method."}
    },
    gui_label="CList",
)
class CList(CData):
    """List container type for collections of CData objects."""

    def __init__(
        self, items: Optional[List[Any]] = None, parent=None, name=None, **kwargs
    ):
        super().__init__(parent=parent, name=name, **kwargs)
        object.__setattr__(self, '_items', items or [])
        object.__setattr__(self, '_item_type', None)
        object.__setattr__(self, '_item_qualifiers', {})

        # Register existing items as children
        for i, item in enumerate(self._items):
            if isinstance(item, CData):
                item.set_parent(self)
                item._name = f"[{i}]"

    def append(self, item: Any) -> None:
        """Add an item to the list.

        Supports smart type conversion:
        - If appending a string to a list expecting CDataFile objects,
          creates a new file object and sets its path
        """
        from ..base_object.base_classes import ValueState, CDataFile
        index = len(self._items)

        # Smart type conversion for legacy code compatibility
        # If appending a string but subItem expects a CDataFile, create one
        sub_item_def = self.get_qualifier('subItem')
        if sub_item_def and isinstance(item, str):
            item_class = sub_item_def.get('class')
            if item_class and issubclass(item_class, CDataFile):
                # Create file object and set its path
                file_obj = item_class(parent=None, name=f"temp_item_{index}")
                file_obj.setFullPath(item)
                item = file_obj

        # If item is CData, register as child
        if isinstance(item, CData):
            item.set_parent(self)
            # Set hierarchical name to just the index (no parent name)
            # The parent's name will be included in path_from_root() traversal
            # This produces paths like "list_name[0]" not "list_name.list_name[0]"
            item._name = f"[{index}]"

        self._items.append(item)

        # Mark as explicitly set
        self._value_states["_items"] = ValueState.EXPLICITLY_SET

    def insert(self, index: int, item: Any) -> None:
        """Insert an item at specified index."""
        if isinstance(item, CData):
            item.set_parent(self)
            # Set hierarchical name to just the index
            item._name = f"[{index}]"

        self._items.insert(index, item)

        # Update names of subsequent items
        for i in range(index + 1, len(self._items)):
            if isinstance(self._items[i], CData):
                self._items[i]._name = f"[{i}]"

        # Mark as explicitly set
        self._value_states["_items"] = ValueState.EXPLICITLY_SET

    def remove(self, item: Any) -> None:
        """Remove an item from the list."""
        index = self._items.index(item)
        self._items.remove(item)

        # Update names of subsequent items
        for i in range(index, len(self._items)):
            if isinstance(self._items[i], CData):
                self._items[i]._name = f"[{i}]"

        # Mark as explicitly set
        self._value_states["_items"] = ValueState.EXPLICITLY_SET

    def pop(self, index: int = -1) -> Any:
        """Remove and return item at index."""
        item = self._items.pop(index)

        # Update names of subsequent items if needed
        if index >= 0:
            for i in range(index, len(self._items)):
                if isinstance(self._items[i], CData):
                    self._items[i]._name = f"[{i}]"

        # Mark as explicitly set
        self._value_states["_items"] = ValueState.EXPLICITLY_SET
        return item

    def clear(self) -> None:
        """Remove all items from the list."""
        self._items.clear()
        self._value_states["_items"] = ValueState.EXPLICITLY_SET

    def set(self, value=None, validate=False):
        """
        Replace list contents with items from value.

        Legacy API compatibility method for CList. Accepts a list or CList
        and replaces all current items with the new items.

        If items are dicts and this list has a subItem qualifier with a 'class',
        the dicts will be converted to proper typed objects.

        Args:
            value: List, CList, or single item to set. If None or empty list, clears the list.
            validate: If True, validate items before setting (default: False)

        Returns:
            self (for method chaining)
        """
        # Handle None or empty default
        if value is None:
            value = []

        # Convert single item to list
        if not isinstance(value, (list, CList)):
            value = [value]

        # Validate if requested
        if validate:
            from .error_reporting import SEVERITY_WARNING
            v = self.validity(value)
            if v.maxSeverity() > SEVERITY_WARNING:
                raise v

        # Get subItem class for dict conversion
        # First check instance qualifier, then fallback to SUBITEM class attribute
        sub_item_def = self.get_qualifier('subItem')
        if not sub_item_def:
            # Check for SUBITEM class attribute (legacy CCP4i2 pattern)
            sub_item_def = getattr(self.__class__, 'SUBITEM', None)
        item_class = sub_item_def.get('class') if isinstance(sub_item_def, dict) else None

        # Clear current items
        self._items.clear()

        # Add new items using append to ensure proper parent/name setup
        for item in value:
            # Convert dict to typed object if we have a subItem class
            if isinstance(item, dict) and item_class is not None:
                # Create a new instance of the item class
                new_item = item_class(parent=None, name=None)
                # Update the new item with dict values
                if hasattr(new_item, 'update'):
                    new_item.update(item)
                item = new_item
            elif isinstance(item, CData):
                # IMPORTANT: Deep-copy CData items to avoid re-parenting issues
                # When copying from another CList, we must create new instances
                # otherwise the original CList loses its items (they get re-parented)
                item_type = type(item)
                new_item = item_type(parent=None, name=None)
                # Copy data from source item
                if hasattr(item, 'get') and callable(item.get):
                    item_data = item.get()
                    if hasattr(new_item, 'set') and callable(new_item.set):
                        new_item.set(item_data)
                    elif hasattr(new_item, 'update') and callable(new_item.update):
                        new_item.update(item_data)
                item = new_item
            self.append(item)

        return self

    def __len__(self) -> int:
        return len(self._items)

    def __getitem__(self, index: int) -> Any:
        return self._items[index]

    def __setitem__(self, index: int, value: Any) -> None:
        if isinstance(value, CData):
            value.set_parent(self)
            value.name = f"{self.name}[{index}]"

        self._items[index] = value
        self._value_states["_items"] = ValueState.EXPLICITLY_SET

    def __iter__(self):
        return iter(self._items)

    def __contains__(self, item: Any) -> bool:
        return item in self._items

    def __delattr__(self, name):
        """Prevent deletion of _items attribute during garbage collection."""
        if name == '_items':
            # Don't allow deletion of _items to prevent AttributeError in __str__
            return
        super().__delattr__(name)

    def __str__(self) -> str:
        return f"CList({len(self._items)} items)"

    def makeItem(self) -> CData:
        """
        Create a new item for this list based on the subItem qualifier.

        This method provides backward compatibility with legacy CCP4i2 code.
        The subItem qualifier should contain:
            {
                'class': SomeClass,  # The class to instantiate
                'qualifiers': {...}  # Qualifiers for the new item
            }

        Returns:
            A new instance of the item class with qualifiers applied

        Raises:
            ValueError: If subItem qualifier is not properly configured

        Example:
            >>> my_list = CList(name="numbers")
            >>> my_list.set_qualifier("subItem", {
            ...     'class': CInt,
            ...     'qualifiers': {'min': 0, 'max': 100, 'default': 50}
            ... })
            >>> new_item = my_list.makeItem()
            >>> my_list.append(new_item)
        """
        # Get subItem qualifier or SUBITEM class attribute
        sub_item_def = self.get_qualifier('subItem')

        if not sub_item_def:
            # Check for SUBITEM class attribute (legacy CCP4i2 pattern)
            sub_item_def = getattr(self.__class__, 'SUBITEM', None)

        if not sub_item_def:
            # Default to CString for simple lists (common in legacy plugins)
            # This matches the behavior when lists hold string values like "CA"
            sub_item_def = {'class': CString, 'qualifiers': {}}

        if not isinstance(sub_item_def, dict):
            raise ValueError(
                f"CList '{self.objectName()}' subItem qualifier must be a dict, "
                f"got {type(sub_item_def).__name__}"
            )

        # Extract class and qualifiers
        item_class = sub_item_def.get('class')
        item_qualifiers = sub_item_def.get('qualifiers', {})

        if not item_class:
            raise ValueError(
                f"CList '{self.objectName()}' subItem qualifier missing 'class' key"
            )

        # Instantiate the item with qualifiers
        # Set parent to self and name to [?] placeholder for JSON encoder
        # Frontend will replace [?] with actual index when cloning for new items
        item = item_class(parent=self, name="[?]")

        # Apply qualifiers if provided
        if item_qualifiers and hasattr(item, '_qualifiers'):
            item._qualifiers.update(item_qualifiers)

        return item

    def dataOrder(self) -> list:
        """Return list of item names/indices for compatibility with legacy API.

        For CList objects, returns indices of items in the list.

        Returns:
            List of string indices representing items in the list
        """
        return [str(i) for i in range(len(self._items))]

    def getEtree(self, name: str = None, excludeUnset: bool = False, allSet: bool = False) -> 'ET.Element':
        """Override getEtree to serialize CList with proper item structure.

        For CList, we create a container element with child elements for each item.
        The child element tags are the class names of the items (CCP4i2 convention).

        IMPORTANT: For list items, presence in the list IS the explicit state.
        When a user adds an empty item to a list, that's an explicit action.
        We do NOT filter out items based on whether their fields are set -
        the list membership itself is the "set" state that must be preserved.
        The excludeUnset/allSet flags are passed to child getEtree() calls
        to filter the item's FIELDS, not the items themselves.

        Args:
            name: Optional element name (not used for CList, uses objectName instead)
            excludeUnset: If True, pass this flag to item's getEtree() calls (for field filtering)
            allSet: If True, pass this flag to item's getEtree() calls (for field filtering)

        Returns:
            ET.Element with list items as children

        Example XML structure:
            <UNMERGEDFILES>
                <CImportUnmerged>
                    <file>...</file>
                    <crystalName>...</crystalName>
                </CImportUnmerged>
                <CImportUnmerged>
                    ...
                </CImportUnmerged>
            </UNMERGEDFILES>
        """
        import xml.etree.ElementTree as ET

        # Create container element with the list's name
        # Use provided name if given, otherwise use objectName
        list_name = name or self.objectName() or "list"
        container_elem = ET.Element(list_name)

        # Add each item as a child element
        # NOTE: We always include ALL items in the list - do not filter based on
        # excludeUnset/allSet at the item level. The list membership itself is
        # the explicit user action. These flags only affect field-level filtering
        # within each item's getEtree() call.
        for item in self._items:
            # Use the item's class name as the XML tag (CCP4i2 convention)
            item_class_name = type(item).__name__

            if isinstance(item, CData):
                # If item is CData, recursively get its etree
                # Pass excludeUnset/allSet to filter FIELDS within the item
                if hasattr(item, 'getEtree'):
                    item_etree = item.getEtree(excludeUnset=excludeUnset, allSet=allSet)
                    # Change the tag to the class name
                    item_etree.tag = item_class_name
                    # Always append - the item exists in the list, so serialize it
                    # Even an empty <CImportUnmerged /> is valid - it represents
                    # a list slot the user created that they may fill in later
                    container_elem.append(item_etree)
                else:
                    # Fallback: create element with item's string representation
                    item_elem = ET.Element(item_class_name)
                    item_elem.text = str(item)
                    container_elem.append(item_elem)
            else:
                # For non-CData items (primitives), create simple element
                item_elem = ET.Element(item_class_name)
                item_elem.text = str(item)
                container_elem.append(item_elem)

        return container_elem

    def setEtree(self, element, ignore_missing: bool = False, preserve_state: bool = False):
        """Override setEtree to deserialize CList with proper item structure.

        For CList, child elements use class names as tags (CCP4i2 convention).
        We create new items using makeItem() and populate them from the XML.

        Args:
            element: xml.etree.ElementTree.Element containing list items
            ignore_missing: If True, ignore missing attributes (passed to item's setEtree)
            preserve_state: If True, don't mark values as EXPLICITLY_SET when deserializing

        Example XML structure:
            <UNMERGEDFILES>
                <CImportUnmerged>
                    <file>...</file>
                    <crystalName>...</crystalName>
                </CImportUnmerged>
                <CImportUnmerged>
                    ...
                </CImportUnmerged>
            </UNMERGEDFILES>
        """
        # Clear existing items
        self.clear()

        # Each child element represents a list item
        # The tag is the item's class name (e.g., <CImportUnmerged>)
        for child_elem in element:
            # Create a new item using makeItem()
            # makeItem() uses the subItem qualifier to determine the correct type
            new_item = self.makeItem()

            # Recursively deserialize the item's content
            if hasattr(new_item, 'setEtree'):
                new_item.setEtree(child_elem, ignore_missing=ignore_missing, preserve_state=preserve_state)

            # Add the item to the list
            self.append(new_item)

    def addItem(self):
        """
        Legacy API: Create and add a new item to the list.

        Returns:
            The newly created item

        Note: This is a legacy ccp4i2 API method. Modern code should use makeItem().
        """
        new_item = self.makeItem()
        self.append(new_item)
        return new_item

    def validity(self):
        """
        Validate the list against qualifiers and children.

        Checks:
        - List length >= listMinLength qualifier (if set)
        - List length <= listMaxLength qualifier (if set)
        - All children are valid (aggregates their validity errors)

        Returns:
            CErrorReport containing validation errors/warnings
        """
        from .error_reporting import CErrorReport, SEVERITY_ERROR

        report = CErrorReport()

        # Get full object path for error reporting (matches frontend _objectPath)
        obj_path = self.object_path() if hasattr(self, 'object_path') else ''

        # Check listMinLength
        min_length = self.get_qualifier('listMinLength')
        if min_length is not None and len(self) < min_length:
            report.append(
                klass=self.__class__.__name__,
                code=101,
                details=f'List must have at least {min_length} item(s), has {len(self)}',
                name=obj_path,
                severity=SEVERITY_ERROR
            )

        # Check listMaxLength
        max_length = self.get_qualifier('listMaxLength')
        if max_length is not None and len(self) > max_length:
            report.append(
                klass=self.__class__.__name__,
                code=102,
                details=f'List must have at most {max_length} item(s), has {len(self)}',
                name=obj_path,
                severity=SEVERITY_ERROR
            )

        # Validate each child and aggregate their errors
        for item in self._items:
            if hasattr(item, 'validity'):
                child_report = item.validity()
                report.extend(child_report)

        return report


# NOTE: Type aliases removed - all custom types now have proper stub classes
# in core/cdata_stubs/ and implementation classes in core/
# (COneWord, CUUID, CProjectId, CCellLength, CWavelength, etc. are all
# regular custom classes derived from fundamental types)
#CProjectName = CString
#CDatasetName = CString
#CFilePath = CString
#CHostName = CString
#CJobStatus = CInt