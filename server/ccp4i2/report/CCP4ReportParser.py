"""
CCP4i2 Report Parser — backward-compatibility shim.

This module previously contained all report element classes (~2950 lines).
It has been split into focused sub-modules:

    core.py      — ReportClass, Container, IfContainer, Loop, Report, utilities
    elements.py  — Text, Pre, Table, Fold, Div, Generic, Progress, etc.
    graphs.py    — Graph, FlotGraph, GraphGroup, PictureGroup, ObjectGallery, etc.
    pictures.py  — Picture
    io_data.py   — IODataList, InputData, OutputData, ImportedFiles
    metadata.py  — Title, JobDetails, JobLogFiles, Reference, ReferenceGroup, etc.
    actions.py   — Help, Launch, Download, CopyToClipboard, etc.

All public names are re-exported here so that existing imports like:
    from ccp4i2.report.CCP4ReportParser import Report
continue to work unchanged.
"""

# Core classes and utilities
from ccp4i2.report.core import (  # noqa: F401
    CURRENT_CSS_VERSION,
    CCP4NS,
    DOCTYPE,
    I2XmlParser,
    PARSER,
    XHTMLNS,
    XRTNS,
    Container,
    IfContainer,
    Loop,
    Report,
    ReportClass,
    applySelect,
    findChildren,
    findallEval,
    getChildObject,
    htmlBase,
    saveToFile,
    testPathExists,
    toBoolean,
)

# Element classes
from ccp4i2.report.elements import (  # noqa: F401
    BaseTable,
    Copy,
    Div,
    DrillDown,
    FetchPre,
    Fold,
    Generic,
    GenericElement,
    Plot,
    Pre,
    Progress,
    Results,
    Status,
    Table,
    Text,
    foldLinkLine,
    foldTitleLine,
    parse_from_unicode,
    _set_cell_content,
)

# Graph and chart classes
from ccp4i2.report.graphs import (  # noqa: F401
    DAGGraph,
    DrawnDiv,
    FlotGraph,
    FlotGraphGroup,
    Graph,
    GraphGroup,
    GraphLineChooser,
    ObjectGallery,
    PictureGroup,
)

# Picture classes
from ccp4i2.report.pictures import (  # noqa: F401
    Picture,
)

# I/O data classes
from ccp4i2.report.io_data import (  # noqa: F401
    IODataList,
    ImportedFiles,
    InputData,
    OutputData,
)

# Metadata classes
from ccp4i2.report.metadata import (  # noqa: F401
    GenericReport,
    JobDetails,
    JobLogFiles,
    Reference,
    ReferenceGroup,
    Title,
)

# Action classes
from ccp4i2.report.actions import (  # noqa: F401
    CopyToClipboard,
    CopyUrlToClipboard,
    Download,
    Help,
    Launch,
    LaunchTask,
)
