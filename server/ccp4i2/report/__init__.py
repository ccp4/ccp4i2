"""
CCP4i2 Report Module.

This module provides the report generation system for CCP4i2, converting
job output data into structured XML for rendering by the React frontend.

Core Classes:
    Report - Main report container (from CCP4ReportParser)
    ReportElement - Modern base class for report elements
    ReportContainer - Base class for container elements

Grid Layout:
    GridContainer - MUI Grid container for responsive layouts
    GridItem - Grid item with responsive column spans
    GridSpan - Column span configuration
    GridPresets - Common span configurations

Error Handling:
    DiagnosticCollector - Collects errors and warnings
    ReportDiagnostic - Individual diagnostic entry
    ReportErrorCodes - Standard error codes

Example:
    from ccp4i2.report import Report, GridContainer, GridItem, GridSpan

    class MyReport(Report):
        TASKNAME = "mytask"

        def __init__(self, xmlnode=None, jobInfo={}, **kw):
            super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)

            # Create a two-column layout
            grid = GridContainer(spacing=2)

            left = GridItem(span=GridSpan(xs=12, md=8))
            left.append(self.addTable(title="Results"))
            grid.append(left)

            right = GridItem(span=GridSpan(xs=12, md=4))
            right.append(self.addGraph(title="Plot"))
            grid.append(right)

            self.append(grid)
"""

# Re-export main classes from CCP4ReportParser for backward compatibility
from ccp4i2.report.CCP4ReportParser import (
    Report,
    Container,
    Fold,
    Results,
    Text,
    Pre,
    Table,
    Graph,
    FlotGraph,
    GraphGroup,
    FlotGraphGroup,
    Picture,
    PictureGroup,
    ObjectGallery,
    Title,
    JobDetails,
    Reference,
    ReferenceGroup,
    InputData,
    OutputData,
    ImportedFiles,
    Div,
    DrawnDiv,
    Progress,
    Generic,
    GenericElement,
    GenericReport,
    # Utilities
    PARSER,
    CCP4NS,
    XRTNS,
)

# New modern base classes
from ccp4i2.report.base import (
    ReportElement,
    ReportContainer,
    ElementRegistry,
)

# Grid layout system
from ccp4i2.report.grid import (
    GridContainer,
    GridItem,
    GridRow,
    GridSpan,
    GridPresets,
    GridDirection,
    GridJustify,
    GridAlign,
    # Convenience functions
    two_column_layout,
    equal_columns,
    stacked_layout,
)

# Error handling
from ccp4i2.report.errors import (
    DiagnosticCollector,
    ReportDiagnostic,
    DiagnosticLevel,
    ReportErrorCodes,
)

__all__ = [
    # Legacy classes (backward compatibility)
    'Report',
    'Container',
    'Fold',
    'Results',
    'Text',
    'Pre',
    'Table',
    'Graph',
    'FlotGraph',
    'GraphGroup',
    'FlotGraphGroup',
    'Picture',
    'PictureGroup',
    'ObjectGallery',
    'Title',
    'JobDetails',
    'Reference',
    'ReferenceGroup',
    'InputData',
    'OutputData',
    'ImportedFiles',
    'Div',
    'DrawnDiv',
    'Progress',
    'Generic',
    'GenericElement',
    'GenericReport',
    'PARSER',
    'CCP4NS',
    'XRTNS',
    # Modern base classes
    'ReportElement',
    'ReportContainer',
    'ElementRegistry',
    # Grid layout
    'GridContainer',
    'GridItem',
    'GridRow',
    'GridSpan',
    'GridPresets',
    'GridDirection',
    'GridJustify',
    'GridAlign',
    'two_column_layout',
    'equal_columns',
    'stacked_layout',
    # Error handling
    'DiagnosticCollector',
    'ReportDiagnostic',
    'DiagnosticLevel',
    'ReportErrorCodes',
]
