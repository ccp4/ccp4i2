"""
Input/Output data report elements.

IODataList, InputData, OutputData, ImportedFiles.
"""

import xml.etree.ElementTree as etree

from ccp4i2.report.core import ReportClass


class IODataList(ReportClass):

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        super().__init__()
        self.jobInfo = {}
        self.jobInfo.update(jobInfo)

    def _file_ids_for_role(self, role):
        """Return list of file IDs for the given role."""
        if role not in self.jobInfo:
            return []
        return [str(fi['fileId']) for fi in self.jobInfo[role]]

    def _build_data_etree(self, title, roles):
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

    def as_data_etree(self):
        return self._build_data_etree('Input Data', ['inputfiles', 'importedfiles'])


class OutputData(IODataList):

    def as_data_etree(self):
        return self._build_data_etree('Output Data', ['outputfiles'])


class ImportedFiles(IODataList):

    def as_data_etree(self):
        return self._build_data_etree('Imported Files', ['importedfiles'])
