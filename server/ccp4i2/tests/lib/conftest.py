"""Conftest for report fixture tests."""


def pytest_addoption(parser):
    """Add --rebaseline option for updating expected fixture output."""
    parser.addoption(
        "--rebaseline",
        action="store_true",
        default=False,
        help="Update expected_report.xml with actual output for failing fixtures.",
    )
