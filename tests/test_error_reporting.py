"""Test error reporting and validation for CData objects."""

import sys
import os

import pytest

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from core.base_object.base_classes import CContainer
from core.base_object.fundamental_types import CInt, CFloat, CString, CBoolean
from core.base_object.error_reporting import (
    CErrorReport, CException,
    SEVERITY_OK, SEVERITY_WARNING, SEVERITY_ERROR
)


def test_error_report_creation():
    """Test creating an error report."""
    report = CErrorReport()
    assert report.count() == 0
    assert report.maxSeverity() == SEVERITY_OK


def test_error_report_append():
    """Test adding errors to a report."""
    report = CErrorReport()
    report.append("CInt", 101, "Value below minimum", "myInt", SEVERITY_ERROR)

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR


def test_error_report_extend():
    """Test merging two error reports."""
    report1 = CErrorReport()
    report1.append("CInt", 101, "Error 1", "obj1", SEVERITY_ERROR)

    report2 = CErrorReport()
    report2.append("CFloat", 102, "Error 2", "obj2", SEVERITY_WARNING)

    report1.extend(report2)
    assert report1.count() == 2
    assert report1.maxSeverity() == SEVERITY_ERROR


def test_error_report_max_severity():
    """Test maxSeverity() returns the highest severity level."""
    report = CErrorReport()
    report.append("CInt", 101, "Warning", "obj1", SEVERITY_WARNING)
    assert report.maxSeverity() == SEVERITY_WARNING

    report.append("CFloat", 102, "Error", "obj2", SEVERITY_ERROR)
    assert report.maxSeverity() == SEVERITY_ERROR


def test_error_report_string_representation():
    """Test error report string formatting."""
    report = CErrorReport()
    report.append("CInt", 101, "Value too small", "myInt", SEVERITY_ERROR)

    output = str(report)
    assert "ERROR" in output
    assert "CInt" in output
    assert "myInt" in output
    assert "Value too small" in output
    assert "101" in output


def test_error_report_bool():
    """Test that error report evaluates to True if it has warnings/errors."""
    report = CErrorReport()
    assert not bool(report)  # No errors = False

    report.append("CInt", 101, "Warning", "obj", SEVERITY_WARNING)
    assert bool(report)  # Has warning = True


def test_cexception_creation():
    """Test CException creation."""
    exc = CException("CInt", 101, "Value error", "myInt", SEVERITY_ERROR)
    assert exc.count() == 1
    assert exc.maxSeverity() == SEVERITY_ERROR
    assert str(exc.args[0]) == "Value error"


def test_cint_validity_no_constraints():
    """Test CInt validation with no constraints."""
    num = CInt(value=42, name="myInt")
    report = num.validity()

    assert report.count() == 0
    assert report.maxSeverity() == SEVERITY_OK


def test_cint_validity_min_constraint_pass():
    """Test CInt validation passes min constraint."""
    num = CInt(value=10, name="myInt")
    num.set_qualifier("min", 5)
    report = num.validity()

    assert report.count() == 0


def test_cint_validity_min_constraint_fail():
    """Test CInt validation fails min constraint."""
    num = CInt(value=3, name="myInt")
    num.set_qualifier("min", 5)
    report = num.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "below minimum" in str(report)


def test_cint_validity_max_constraint_pass():
    """Test CInt validation passes max constraint."""
    num = CInt(value=10, name="myInt")
    num.set_qualifier("max", 20)
    report = num.validity()

    assert report.count() == 0


def test_cint_validity_max_constraint_fail():
    """Test CInt validation fails max constraint."""
    num = CInt(value=25, name="myInt")
    num.set_qualifier("max", 20)
    report = num.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "above maximum" in str(report)


def test_cint_validity_enumerators_pass():
    """Test CInt validation passes enumerators constraint."""
    num = CInt(value=5, name="myInt")
    num.set_qualifier("enumerators", [1, 5, 10, 20])
    num.set_qualifier("onlyEnumerators", True)
    report = num.validity()

    assert report.count() == 0


def test_cint_validity_enumerators_fail():
    """Test CInt validation fails enumerators constraint."""
    num = CInt(value=7, name="myInt")
    num.set_qualifier("enumerators", [1, 5, 10, 20])
    num.set_qualifier("onlyEnumerators", True)
    report = num.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "not in allowed values" in str(report)


def test_cfloat_validity_min_constraint():
    """Test CFloat validation with min constraint."""
    num = CFloat(value=2.5, name="myFloat")
    num.set_qualifier("min", 5.0)
    report = num.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR


def test_cfloat_validity_max_constraint():
    """Test CFloat validation with max constraint."""
    num = CFloat(value=15.7, name="myFloat")
    num.set_qualifier("max", 10.0)
    report = num.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR


def test_cstring_validity_no_constraints():
    """Test CString validation with no constraints."""
    text = CString(value="hello", name="myString")
    report = text.validity()

    assert report.count() == 0


def test_cstring_validity_min_length_pass():
    """Test CString validation passes minLength."""
    text = CString(value="hello", name="myString")
    text.set_qualifier("minLength", 3)
    report = text.validity()

    assert report.count() == 0


def test_cstring_validity_min_length_fail():
    """Test CString validation fails minLength."""
    text = CString(value="hi", name="myString")
    text.set_qualifier("minLength", 5)
    report = text.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "below minimum" in str(report)


def test_cstring_validity_max_length_pass():
    """Test CString validation passes maxLength."""
    text = CString(value="hello", name="myString")
    text.set_qualifier("maxLength", 10)
    report = text.validity()

    assert report.count() == 0


def test_cstring_validity_max_length_fail():
    """Test CString validation fails maxLength."""
    text = CString(value="very long string here", name="myString")
    text.set_qualifier("maxLength", 10)
    report = text.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "above maximum" in str(report)


def test_cstring_validity_enumerators_pass():
    """Test CString validation passes enumerators."""
    text = CString(value="option1", name="myString")
    text.set_qualifier("enumerators", ["option1", "option2", "option3"])
    text.set_qualifier("onlyEnumerators", True)
    report = text.validity()

    assert report.count() == 0


def test_cstring_validity_enumerators_fail():
    """Test CString validation fails enumerators."""
    text = CString(value="invalid", name="myString")
    text.set_qualifier("enumerators", ["option1", "option2", "option3"])
    text.set_qualifier("onlyEnumerators", True)
    report = text.validity()

    assert report.count() == 1
    assert report.maxSeverity() == SEVERITY_ERROR
    assert "not in allowed values" in str(report)


def test_cboolean_validity():
    """Test CBoolean validation (always valid)."""
    flag = CBoolean(value=True, name="myFlag")
    report = flag.validity()

    assert report.count() == 0
    assert report.maxSeverity() == SEVERITY_OK


def test_multiple_validation_errors():
    """Test object with multiple validation errors."""
    num = CInt(value=50, name="myInt")
    num.set_qualifier("min", 0)
    num.set_qualifier("max", 10)
    report = num.validity()

    # Should have 1 error (above maximum)
    assert report.count() == 1
    assert "above maximum" in str(report)


def test_cdata_validity_base():
    """Test CData base validity() method."""
    from core.base_object.base_classes import CData
    obj = CData(name="test")
    report = obj.validity()

    assert report.count() == 0
    assert report.maxSeverity() == SEVERITY_OK


def test_error_report_clear():
    """Test clearing an error report."""
    report = CErrorReport()
    report.append("CInt", 101, "Error", "obj", SEVERITY_ERROR)
    assert report.count() == 1

    report.clear()
    assert report.count() == 0
    assert report.maxSeverity() == SEVERITY_OK


def test_error_report_get_errors():
    """Test getting list of errors."""
    report = CErrorReport()
    report.append("CInt", 101, "Error 1", "obj1", SEVERITY_ERROR)
    report.append("CFloat", 102, "Error 2", "obj2", SEVERITY_WARNING)

    errors = report.getErrors()
    assert len(errors) == 2
    assert errors[0]['class'] == "CInt"
    assert errors[1]['class'] == "CFloat"
