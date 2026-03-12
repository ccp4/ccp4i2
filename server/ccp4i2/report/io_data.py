"""
Input/Output data report elements.

IODataList, InputData, OutputData, ImportedFiles.

These elements list job input/output files in the report. The frontend
fetches full file metadata client-side; these elements only provide file UUIDs.
"""

from __future__ import annotations

import xml.etree.ElementTree as etree
from typing import Any

from ccp4i2.report.core import ReportClass


class IODataList(ReportClass):
    """Base class for input/output/imported file lists."""

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        super().__init__()
        if jobInfo is None:
            jobInfo = {}
        self.jobInfo: dict[str, Any] = {}
        self.jobInfo.update(jobInfo)

    def _file_ids_for_role(self, role: str) -> list[str]:
        """Return list of file IDs for the given role."""
        if role not in self.jobInfo:
            return []
        return [str(fi['fileId']) for fi in self.jobInfo[role]]

    def _build_data_etree(self, title: str, roles: list[str]) -> etree.Element:
        """Build minimal as_data_etree with just file IDs and title.

        The frontend (CCP4i2ReportInputOutputData) only needs:
        - <h5> elements for the accordion title
        - <div id="input_file_{UUID}"> elements to extract file UUIDs
        All file metadata is fetched client-side from the project's file list.
        """
        root = super().as_data_etree()
        inner = etree.Element('root')
        head = etree.Element('h5')
        head.text = title
        inner.append(head)
        for role in roles:
            for file_id in self._file_ids_for_role(role):
                div = etree.Element('div')
                div.set('id', 'input_file_' + file_id)
                inner.append(div)
        root.append(inner)
        return root


class InputData(IODataList):
    """Input files list for the report."""

    def as_data_etree(self) -> etree.Element:
        return self._build_data_etree('Input Data', ['inputfiles', 'importedfiles'])


class OutputData(IODataList):
    """Output files list for the report."""

    def as_data_etree(self) -> etree.Element:
        return self._build_data_etree('Output Data', ['outputfiles'])


class ImportedFiles(IODataList):
    """Imported files list for the report."""

    def as_data_etree(self) -> etree.Element:
        return self._build_data_etree('Imported Files', ['importedfiles'])
