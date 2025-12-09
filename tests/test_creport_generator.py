"""
Unit tests for the modernized CReportGenerator.

These tests verify that the CReportGenerator class works correctly after
being modernized to remove Qt dependencies and use Django for database access.
"""
import sys
import os
import pytest

# Ensure the main report module is used, not tests/report
# The main report/ directory should be in the path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
report_path = os.path.join(project_root, 'report')

# Remove tests/report from path if present
for p in sys.path[:]:
    if 'tests/report' in p or p.endswith('tests'):
        # Don't remove tests itself, but ensure report/ is earlier in path
        pass

# Add main report directory early in path
if report_path not in sys.path:
    sys.path.insert(0, project_root)

# Now import from the main report module
import importlib
if 'report.CCP4ReportGenerator' in sys.modules:
    # Reload to get the correct module
    del sys.modules['report.CCP4ReportGenerator']


class TestSignalClass:
    """Tests for the Qt-free Signal replacement class."""

    def test_signal_import(self):
        """Test that Signal can be imported without Qt."""
        from report.CCP4ReportGenerator import Signal
        assert Signal is not None

    def test_signal_instantiation(self):
        """Test that Signal can be instantiated."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()
        assert sig._callbacks == []

    def test_signal_connect(self):
        """Test connecting a callback to a signal."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()

        results = []
        def callback(value):
            results.append(value)

        sig.connect(callback)
        assert len(sig._callbacks) == 1

    def test_signal_emit(self):
        """Test emitting a signal calls connected callbacks."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()

        results = []
        def callback(value):
            results.append(value)

        sig.connect(callback)
        sig.emit("test_value")

        assert results == ["test_value"]

    def test_signal_emit_multiple_callbacks(self):
        """Test that emit calls all connected callbacks."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()

        results1 = []
        results2 = []

        sig.connect(lambda v: results1.append(v))
        sig.connect(lambda v: results2.append(v))

        sig.emit("data")

        assert results1 == ["data"]
        assert results2 == ["data"]

    def test_signal_disconnect(self):
        """Test disconnecting a callback from a signal."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()

        results = []
        def callback(value):
            results.append(value)

        sig.connect(callback)
        sig.disconnect(callback)

        sig.emit("test")
        assert results == []  # Should not have been called

    def test_signal_emit_with_exception_handling(self):
        """Test that exceptions in callbacks don't break emission."""
        from report.CCP4ReportGenerator import Signal
        sig = Signal()

        results = []

        def bad_callback(value):
            raise ValueError("Intentional error")

        def good_callback(value):
            results.append(value)

        sig.connect(bad_callback)
        sig.connect(good_callback)

        # Should not raise, and should still call the good callback
        sig.emit("test")
        assert results == ["test"]


class TestCReportGeneratorClass:
    """Tests for the CReportGenerator class itself."""

    def test_import_without_qt(self):
        """Test that CReportGenerator can be imported without Qt."""
        from report.CCP4ReportGenerator import CReportGenerator
        assert CReportGenerator is not None

    def test_instantiation(self):
        """Test basic instantiation of CReportGenerator."""
        from report.CCP4ReportGenerator import CReportGenerator

        # Should be able to create with just a job ID
        gen = CReportGenerator(jobId="test-job-id-123")

        assert gen.jobId == "test-job-id-123"
        assert gen.jobStatus == "Finished"
        assert gen.reportFile is None
        assert gen.jobInfo is None

    def test_instantiation_with_status(self):
        """Test instantiation with custom status."""
        from report.CCP4ReportGenerator import CReportGenerator

        gen = CReportGenerator(
            jobId="test-id",
            jobStatus="Running"
        )

        assert gen.jobStatus == "Running"

    def test_to_delete_status_becomes_finished(self):
        """Test that 'To delete' status becomes 'Finished'."""
        from report.CCP4ReportGenerator import CReportGenerator

        gen = CReportGenerator(
            jobId="test-id",
            jobStatus="To delete"
        )

        assert gen.jobStatus == "Finished"

    def test_job_number_kwarg(self):
        """Test that jobNumber can be passed as kwarg."""
        from report.CCP4ReportGenerator import CReportGenerator

        gen = CReportGenerator(
            jobId="test-id",
            jobNumber="1.2.3"
        )

        assert gen.jobNumber == "1.2.3"

    def test_finished_pictures_signal(self):
        """Test that each instance has its own FinishedPictures signal."""
        from report.CCP4ReportGenerator import CReportGenerator, Signal

        gen1 = CReportGenerator(jobId="job1")
        gen2 = CReportGenerator(jobId="job2")

        # Each should have its own signal instance
        assert gen1.FinishedPictures is not gen2.FinishedPictures
        assert isinstance(gen1.FinishedPictures, Signal)
        assert isinstance(gen2.FinishedPictures, Signal)

    def test_set_job_status(self):
        """Test setJobStatus method."""
        from report.CCP4ReportGenerator import CReportGenerator

        gen = CReportGenerator(jobId="test-id")
        gen.jobInfo = {"some": "data"}

        gen.setJobStatus("Running")
        assert gen.jobStatus == "Running"
        # Running is not a finished status, so jobInfo should remain
        assert gen.jobInfo is not None

    def test_set_job_status_finished_clears_job_info(self):
        """Test that setting a finished status clears jobInfo."""
        from report.CCP4ReportGenerator import CReportGenerator

        gen = CReportGenerator(jobId="test-id", jobStatus="Running")
        gen.jobInfo = {"some": "data"}

        gen.setJobStatus("Finished")
        assert gen.jobInfo is None

    def test_hierarchical_object_inheritance(self):
        """Test that CReportGenerator inherits from HierarchicalObject."""
        from report.CCP4ReportGenerator import CReportGenerator
        from ccp4i2.core.base_object.hierarchy_system import HierarchicalObject

        gen = CReportGenerator(jobId="test-id")

        # Should be an instance of HierarchicalObject
        assert isinstance(gen, HierarchicalObject)

        # Should have HierarchicalObject methods
        assert hasattr(gen, 'get_parent')
        assert hasattr(gen, 'children')
        assert hasattr(gen, 'object_path')

        # Name should be set
        assert gen._name is not None
        assert "test-id" in gen._name

    def test_hierarchical_parent_child(self):
        """Test parent-child relationship between generators."""
        from report.CCP4ReportGenerator import CReportGenerator

        parent_gen = CReportGenerator(jobId="parent-job")
        child_gen = CReportGenerator(jobId="child-job", parent=parent_gen)

        # Child should have parent set
        assert child_gen.get_parent() == parent_gen

        # Parent should have child in children list
        assert child_gen in parent_gen.children()


class TestGetReportJobInfoFunction:
    """Tests for the getReportJobInfo function."""

    def test_import(self):
        """Test that getReportJobInfo can be imported."""
        from report.CCP4ReportGenerator import getReportJobInfo
        assert callable(getReportJobInfo)


class TestModernizationCompleteness:
    """Tests verifying the modernization is complete."""

    def test_no_qt_imports(self):
        """Verify no Qt imports remain in the module (excluding comments/docstrings)."""
        import report.CCP4ReportGenerator as module
        import re

        source_file = module.__file__
        with open(source_file, 'r') as f:
            content = f.read()

        # Remove docstrings and comments for checking actual code
        # Remove triple-quoted strings (docstrings)
        code_only = re.sub(r'""".*?"""', '', content, flags=re.DOTALL)
        code_only = re.sub(r"'''.*?'''", '', code_only, flags=re.DOTALL)
        # Remove single-line comments
        code_only = re.sub(r'#.*$', '', code_only, flags=re.MULTILINE)

        # Should not have PySide2 or QtCore imports in actual code
        assert 'from PySide2' not in code_only, "PySide2 import still present"
        assert 'import PySide2' not in code_only, "PySide2 import still present"
        assert 'QtCore.QObject' not in code_only, "QtCore.QObject still present"
        assert 'QtCore.Signal' not in code_only, "QtCore.Signal still present"
        assert '@QtCore.Slot' not in code_only, "QtCore.Slot decorator still present"

    def test_no_legacy_dbapi_import(self):
        """Verify legacy dbapi import has been replaced."""
        import report.CCP4ReportGenerator as module

        source_file = module.__file__
        with open(source_file, 'r') as f:
            content = f.read()

        # Should not have legacy CCP4DbApi import
        assert 'from dbapi import CCP4DbApi' not in content, \
            "Legacy dbapi import still present"

    def test_uses_static_data_constants(self):
        """Verify the module uses ccp4x.db.ccp4i2_static_data constants."""
        import report.CCP4ReportGenerator as module

        source_file = module.__file__
        with open(source_file, 'r') as f:
            content = f.read()

        # Should import from ccp4i2_static_data
        assert 'ccp4i2_static_data' in content, \
            "Should import from ccp4x.db.ccp4i2_static_data"

    def test_error_codes_defined(self):
        """Verify ERROR_CODES are still defined."""
        from report.CCP4ReportGenerator import CReportGenerator

        assert hasattr(CReportGenerator, 'ERROR_CODES')
        assert isinstance(CReportGenerator.ERROR_CODES, dict)
        assert len(CReportGenerator.ERROR_CODES) >= 10
