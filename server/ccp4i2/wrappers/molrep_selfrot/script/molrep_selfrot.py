from __future__ import print_function

from ccp4i2.wrappers.molrep_mr.script import molrep_mr
from lxml import etree
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING

class molrep_selfrot(molrep_mr.molrep_mr):
    TASKNAME='molrep_selfrot'

    ERROR_CODES = {
        201: {'severity': SEVERITY_WARNING, 'description': 'Failed to set PERFORM parameter for self-rotation'},
        202: {'severity': SEVERITY_WARNING, 'description': 'Failed to call parent processInputFiles'},
        203: {'severity': SEVERITY_WARNING, 'description': 'Failed to call parent processOutputFiles'},
        204: {'severity': SEVERITY_WARNING, 'description': 'Failed to create MolrepResult XML element'},
        205: {'severity': SEVERITY_WARNING, 'description': 'Failed to scrape molrep.doc file'},
        206: {'severity': SEVERITY_WARNING, 'description': 'Failed to write program XML file'},
        207: {'severity': SEVERITY_WARNING, 'description': 'molrep.doc.txt file not found'},
        208: {'severity': SEVERITY_WARNING, 'description': 'Failed to parse Structure Factors block'},
        209: {'severity': SEVERITY_WARNING, 'description': 'Failed to parse Patterson block'},
        210: {'severity': SEVERITY_WARNING, 'description': 'Failed to parse Self-Rotation block'},
        211: {'severity': SEVERITY_WARNING, 'description': 'Failed to parse peak line'},
        212: {'severity': SEVERITY_WARNING, 'description': 'Failed to parse rotation line'},
    }

    def __init__(self, *args, **kws):
        super(molrep_selfrot, self).__init__(*args, **kws)

    def processInputFiles(self):
        # Set PERFORM to 'srf' (self-rotation function) before processing
        try:
            self.container.guiParameters.PERFORM.set('srf')
        except Exception as e:
            self.appendErrorReport(201, str(e))
            return CPluginScript.FAILED

        try:
            return super(molrep_selfrot, self).processInputFiles()
        except Exception as e:
            self.appendErrorReport(202, str(e))
            return CPluginScript.FAILED

    def processOutputFiles(self):
        try:
            result = super(molrep_selfrot, self).processOutputFiles()
            if result == CPluginScript.FAILED:
                return CPluginScript.FAILED
        except Exception as e:
            self.appendErrorReport(203, str(e))
            return CPluginScript.FAILED

        try:
            self.xmlnode = etree.Element('MolrepResult')
        except Exception as e:
            self.appendErrorReport(204, str(e))
            return CPluginScript.FAILED

        self.scrapeDocFile()

        try:
            with open(self.makeFileName('PROGRAMXML'), 'w') as xmlFile:
                CCP4Utils.writeXML(xmlFile, etree.tostring(self.xmlnode, pretty_print=True))
        except Exception as e:
            self.appendErrorReport(206, str(e))
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def scrapeDocFile(self):
        import os
        docFilepath = os.path.join(self.workDirectory, 'molrep.doc.txt')

        if not os.path.exists(docFilepath):
            self.appendErrorReport(207, docFilepath)
            return

        try:
            with open(docFilepath, 'r') as docFile:
                lines = docFile.readlines()
        except Exception as e:
            self.appendErrorReport(205, f'Failed to read file: {str(e)}')
            return

        nearStructureFactorBlock = False
        inStructureFactorBlock = False
        inPattersonBlock = False
        structureFactorNode = None
        pattersonNode = None

        for line in lines:
            strippedLine = line.strip()

            # Parse Structure Factors block
            if strippedLine.startswith('-- Structure Factors --'):
                try:
                    structureFactorNode = etree.SubElement(self.xmlnode, 'StructureFactors')
                    nearStructureFactorBlock = True
                except Exception as e:
                    self.appendErrorReport(208, f'Failed to create StructureFactors element: {str(e)}')
                    continue

            elif nearStructureFactorBlock and strippedLine.startswith('---------'):
                inStructureFactorBlock = True
                nearStructureFactorBlock = False

            elif inStructureFactorBlock:
                try:
                    if strippedLine.startswith('-------------'):
                        inStructureFactorBlock = False
                    elif strippedLine.startswith('Fobs: resolution'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'LowResProvided', tokens[3])
                        self.subElementOfTypeWithText(structureFactorNode, 'HighResProvided', tokens[4])
                    elif strippedLine.startswith('Completeness of Fobs :'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'Completeness', tokens[4])
                    elif strippedLine.startswith('B_overall of Fobs    :'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'BOverall', tokens[4])
                    elif strippedLine.startswith('Optical resolution   :'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'OpticalResolution', tokens[3])
                    elif strippedLine.startswith('Resmax (from Opt.res):'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'OpticalHighRes', tokens[3])
                    elif strippedLine.startswith('Ratio of Eigen values :'):
                        tokens = strippedLine.split()
                        self.subElementOfTypeWithText(structureFactorNode, 'EigenValueRatioH', tokens[5])
                        self.subElementOfTypeWithText(structureFactorNode, 'EigenValueRatioK', tokens[6])
                        self.subElementOfTypeWithText(structureFactorNode, 'EigenValueRatioL', tokens[7])
                    elif strippedLine.startswith('INFO'):
                        self.subElementOfTypeWithText(structureFactorNode, 'INFO', strippedLine[5:])
                except Exception as e:
                    self.appendErrorReport(208, f'Error parsing line "{strippedLine[:50]}": {str(e)}')

            # Parse Patterson block
            elif strippedLine.startswith('--- Patterson ---'):
                try:
                    inPattersonBlock = True
                    pattersonNode = etree.SubElement(self.xmlnode, 'Patterson')
                except Exception as e:
                    self.appendErrorReport(209, f'Failed to create Patterson element: {str(e)}')
                    continue

            elif inPattersonBlock:
                try:
                    splitLine = strippedLine.split()
                    if strippedLine.startswith('INFO'):
                        self.subElementOfTypeWithText(pattersonNode, 'INFO', strippedLine.split(':')[1])
                        inPattersonBlock = False
                    elif strippedLine.startswith('Number of peaks :'):
                        self.subElementOfTypeWithText(pattersonNode, 'NPeaks', strippedLine.split(':')[1])
                    elif len(splitLine) == 12 and splitLine[0].isdigit():
                        self.parsePeakLine(strippedLine, pattersonNode)
                except Exception as e:
                    self.appendErrorReport(209, f'Error parsing Patterson line "{strippedLine[:50]}": {str(e)}')

            # Parse Self-Rotation block
            elif strippedLine.startswith('Sol_Rf'):
                try:
                    existingNodes = self.xmlnode.xpath('//SelfRotation')
                    if existingNodes:
                        rotationsNode = existingNodes[0]
                    else:
                        rotationsNode = etree.SubElement(self.xmlnode, 'SelfRotation')
                    self.parseRotationLine(strippedLine, rotationsNode)
                except Exception as e:
                    self.appendErrorReport(210, f'Error parsing rotation line "{strippedLine[:50]}": {str(e)}')

    def parsePeakLine(self, line, pattersonNode):
        if line.strip().startswith('IX'):
            return
        headings = 'No IX IY IZ Xfrac Yfrac Zfrac Xort Yort Zort Dens Dens_sigma'.split()
        values = line.split()
        try:
            peakNode = etree.SubElement(pattersonNode, 'Peak')
            for iHeading in range(len(headings)):
                if iHeading < len(values):
                    self.subElementOfTypeWithText(peakNode, headings[iHeading], values[iHeading])
        except Exception as e:
            self.appendErrorReport(211, f'Line: "{line[:50]}": {str(e)}')

    def parseRotationLine(self, line, rotationNode):
        headings = 'No theta phi chi alpha beta gamma Rf Rf_sigma'.split()
        values = line.split()
        try:
            peakNode = etree.SubElement(rotationNode, 'Peak')
            for iHeading in range(len(headings)):
                if iHeading + 1 < len(values):
                    self.subElementOfTypeWithText(peakNode, headings[iHeading], values[iHeading + 1])
        except Exception as e:
            self.appendErrorReport(212, f'Line: "{line[:50]}": {str(e)}')

    def subElementOfTypeWithText(self, parent=None, key=None, text=None):
        newNode = etree.SubElement(parent, str(key))
        newNode.text = str(text)
        return newNode
