# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Picture and visual report elements.

Picture — wraps a CCP4mg/Moorhen scene for inline 3D visualization.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as etree
from typing import Any

from ccp4i2.core.CCP4ErrorHandling import CException


class Picture:
    """Scene visualization element for CCP4mg/Moorhen.

    When a ``sceneFile`` is provided, it is used directly without copying.
    When building a scene from scratch (via ``scene`` param),
    a new scene file is created with a deterministic name based on label.
    """

    ERROR_CODES: dict = {
        101: {
            'description': 'Error reading picture definition'}, 102: {
            'description': 'Error parsing xml from scene file'}, 103: {
                'description': 'Error parsing xml from scene description'}, 104: {
                    'description': 'No scene description provided'}}

    # Class-level counter for generating unique scene file names within a session
    _scene_counter: dict[str, int] = {}

    def __init__(
        self,
        xmlnode: etree.Element | None = None,
        jobInfo: dict[str, Any] | None = None,
        **kw: Any,
    ) -> None:
        if jobInfo is None:
            jobInfo = {}
        import copy
        from ccp4i2.report.actions import Launch
        self.id: str | None = kw.get('id', None)
        self.class_: str | None = kw.get('class_', None)
        self.launchList: list[Launch] = []

        # Track whether we're using an externally-provided scene file
        external_scene_file: str | None = None

        bodyEle = etree.Element('ccp4i2_body')
        sceneEle = etree.Element('scene')
        bodyEle.append(sceneEle)

        sceneRoot: etree.Element | None = None
        if kw.get('scene', None) is not None:
            try:
                sceneRoot = etree.fromstring(kw['scene'], etree.XMLParser(encoding='utf-8'))
            except BaseException:
                raise CException(self.__class__, 103, kw['scene'])
        elif kw.get('sceneFile', None) is not None:
            from ccp4i2.core import CCP4Utils
            fileName = kw.get('sceneFile')
            if fileName[0:8] == '$CCP4I2/':
                fileName = os.path.join(CCP4Utils.getCCP4I2Dir(), fileName[8:])
            if not os.path.exists(fileName):
                raise CException(self.__class__, 101, fileName)
            else:
                try:
                    sceneRoot = etree.fromstring(
                        open(fileName).read(), etree.XMLParser(encoding='utf-8'))
                    # Remember that we have an external scene file - no need to copy
                    external_scene_file = fileName
                except BaseException:
                    raise CException(self.__class__, 102, fileName)

        if sceneRoot is None:
            raise CException(self.__class__, 104, str(kw))

        for child in sceneRoot:
            sceneEle.append(copy.deepcopy(child))

        from ccp4i2.core import CCP4File

        # If an external scene file was provided, use it directly without copying
        if external_scene_file is not None:
            self.picDefFile = external_scene_file
        else:
            # Generate a deterministic scene file name based on job and label
            # Use label if available, otherwise use a counter per job
            job_id = jobInfo.get('jobid', 'unknown')
            label = kw.get('label', None)

            if label:
                # Create filename from label (sanitize for filesystem)
                safe_label = "".join(c if c.isalnum() or c in '_-' else '_' for c in label)
                scene_name = f'scene_{safe_label}.scene.xml'
            else:
                # Use counter for this job
                if job_id not in Picture._scene_counter:
                    Picture._scene_counter[job_id] = 0
                Picture._scene_counter[job_id] += 1
                scene_name = f'scene_{Picture._scene_counter[job_id]}.scene.xml'

            if 'fileroot' in jobInfo and jobInfo['fileroot'] is not None:
                scene_path = jobInfo['fileroot'] + scene_name
            else:
                scene_path = os.path.join(os.getcwd(), scene_name)

            self.picDefFile = CCP4File.CI2XmlDataFile(fullPath=scene_path)
            self.picDefFile.header.setCurrent()
            self.picDefFile.header.function.set('MGSCENE')
            self.picDefFile.header.projectId.set(jobInfo.get('projectid', None))
            self.picDefFile.header.projectName.set(
                jobInfo.get('projectname', None))
            self.picDefFile.header.jobId.set(jobInfo.get('jobid', None))
            self.picDefFile.header.jobNumber.set(jobInfo.get('jobnumber', None))
            self.picDefFile.saveFile(bodyEtree=bodyEle, useLXML=False)

        self.label: str | None = kw.get('label', None)

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'CCP4mg')
        launchNode.set('label', 'View in CCP4mg')
        if self.picDefFile is not None:
            launchNode.set('sceneFile', str(self.picDefFile))
        self.launchList.append(Launch(launchNode, jobInfo=jobInfo))

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'Coot')
        launchNode.set('label', 'View in Coot')
        self.launchList.append(Launch(launchNode, jobInfo=jobInfo))
