"""
Picture and visual report elements.

Picture, PictureGroup, ObjectGallery, DrawnDiv.
"""

import os
import xml.etree.ElementTree as etree

from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2.report.core import (
    ReportClass, Container,
    XRTNS,
    applySelect, PARSER,
)


class Picture:
    """
    Picture element for CCP4mg/Moorhen scene visualization.

    When a sceneFile is provided, it is used directly without copying.
    When building a scene from scratch (via xrtnode or scene param),
    a new scene file is created with a deterministic name based on label.
    """

    ERROR_CODES = {
        101: {
            'description': 'Error reading picture definition'}, 102: {
            'description': 'Error parsing xml from scene file'}, 103: {
                'description': 'Error parsing xml from scene description'}, 104: {
                    'description': 'No scene description provided'}}

    # Class-level counter for generating unique scene file names within a session
    _scene_counter = {}

    def __init__(self, xrtnode=None, xmlnode=None, jobInfo={}, **kw):
        import copy
        from ccp4i2.report.actions import Launch
        self.id = kw.get('id', None)
        self.class_ = kw.get('class_', None)
        self.launchList = []

        # Track whether we're using an externally-provided scene file
        external_scene_file = None

        bodyEle = etree.Element('ccp4i2_body')
        sceneEle = etree.Element('scene')
        bodyEle.append(sceneEle)

        sceneRoot = None
        if xrtnode is not None:
            sceneRoot = xrtnode
        elif kw.get('scene', None) is not None:
            try:
                sceneRoot = etree.fromstring(kw['scene'], PARSER())
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
                        open(fileName).read(), PARSER())
                    # Remember that we have an external scene file - no need to copy
                    external_scene_file = fileName
                except BaseException:
                    raise CException(self.__class__, 102, fileName)

        if sceneRoot is None:
            raise CException(self.__class__, 104, fileName)

        for child in sceneRoot:
            sceneEle.append(copy.deepcopy(child))
            bodyEle = applySelect(bodyEle, xmlnode, jobInfo)

        from ccp4i2.core import CCP4File

        # If an external scene file was provided, use it directly without copying
        if external_scene_file is not None:
            self.picDefFile = external_scene_file
        else:
            # Generate a deterministic scene file name based on job and label
            # Use label if available, otherwise use a counter per job
            job_id = jobInfo.get('jobid', 'unknown')
            label = kw.get('label', None)
            if xrtnode is not None:
                label = xrtnode.get('label', label)

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

        if xrtnode is not None:
            self.label = xrtnode.get('label', None)
        else:
            self.label = kw.get('label', None)

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'CCP4mg')
        launchNode.set('label', 'View in CCP4mg')
        if self.picDefFile is not None:
            launchNode.set('sceneFile', str(self.picDefFile))
        self.launchList.append(Launch(launchNode, jobInfo=jobInfo))

        launchNode = etree.Element('launch')
        launchNode.set('exe', 'Coot')
        launchNode.set('label', 'View in Coot')
        self.launchList.append(Launch(xrtnode=launchNode, jobInfo=jobInfo))
