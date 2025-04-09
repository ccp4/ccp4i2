import os

from lxml import etree

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript
from ...molrep_mr.script import molrep_mr


class molrep_selfrot(molrep_mr.molrep_mr):
    TASKNAME='molrep_selfrot'

    def __init__(self, *args, **kws):
        super(molrep_selfrot, self).__init__(*args,**kws)

    def processOutputFiles(self):
        result = super(molrep_selfrot,self).processOutputFiles()
        if result == CPluginScript.FAILED: return CPluginScript.FAILED
        
        self.xmlnode = etree.Element('MolrepResult')
        self.scrapeDocFile()
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            CCP4Utils.writeXML(xmlFile,etree.tostring(self.xmlnode,pretty_print=True))
        return CPluginScript.SUCCEEDED

    def scrapeDocFile(self):
        docFilepath = os.path.join(self.workDirectory,'molrep.doc.txt')
        lines = []
        with open(docFilepath,'r') as docFile:
            lines = docFile.readlines()

        nearStructureFactorBlock = False
        inStructureFactorBlock = False
        inPattersonBlock = False
        
        for line in lines:
            strippedLine = line.strip()
            if strippedLine.startswith('-- Structure Factors --'):
                structureFactorNode = etree.SubElement(self.xmlnode,'StructureFactors')
                nearStructureFactorBlock = True
            elif nearStructureFactorBlock and strippedLine.startswith('---------'):
                inStructureFactorBlock = True
                nearStructureFactorBlock = False
            elif inStructureFactorBlock:
                if strippedLine.startswith('-------------'):
                    inStructureFactorBlock = False
                elif strippedLine.startswith('Fobs: resolution'):
                    tokens = strippedLine.split()
                    print(structureFactorNode)
                    self.subElementOfTypeWithText(structureFactorNode,'LowResProvided',tokens[3])
                    self.subElementOfTypeWithText(structureFactorNode,'HighResProvided',tokens[4])
                elif strippedLine.startswith('Completeness of Fobs :'):
                    tokens = strippedLine.split()
                    self.subElementOfTypeWithText(structureFactorNode,'Completeness',tokens[4])
                elif strippedLine.startswith('B_overall of Fobs    :'):
                    tokens = strippedLine.split()
                    self.subElementOfTypeWithText(structureFactorNode,'BOverall',tokens[4])
                elif strippedLine.startswith('Optical resolution   :'):
                    tokens = strippedLine.split()
                    self.subElementOfTypeWithText(structureFactorNode,'OpticalResolution',tokens[3])
                elif strippedLine.startswith('Resmax (from Opt.res):'):
                    tokens = strippedLine.split()
                    self.subElementOfTypeWithText(structureFactorNode,'OpticalHighRes',tokens[3])
                elif strippedLine.startswith('Ratio of Eigen values :'):
                    tokens = strippedLine.split()
                    self.subElementOfTypeWithText(structureFactorNode,'EigenValueRatioH',tokens[5])
                    self.subElementOfTypeWithText(structureFactorNode,'EigenValueRatioK',tokens[6])
                    self.subElementOfTypeWithText(structureFactorNode,'EigenValueRatioL',tokens[7])
                elif strippedLine.startswith('INFO'):
                    self.subElementOfTypeWithText(structureFactorNode,'INFO',strippedLine[5:])
            elif strippedLine.startswith('--- Patterson ---'):
                inPattersonBlock = True
                pattersonNode = etree.SubElement(self.xmlnode,'Patterson')
            elif inPattersonBlock:
                splitLine = strippedLine.split()
                if strippedLine.startswith('INFO'):
                    self.subElementOfTypeWithText(pattersonNode, 'INFO', strippedLine.split(':')[1])
                    inPattersonBlock = False
                elif strippedLine.startswith('Number of peaks :'):
                    self.subElementOfTypeWithText(pattersonNode,'NPeaks',strippedLine.split(':')[1])
                elif len(splitLine) == 12 and splitLine[0].isdigit():
                    self.parsePeakLine(strippedLine, pattersonNode)
            elif strippedLine.startswith('Sol_Rf'):
                
                try: rotationsNode = self.xmlnode.xpath('//SelfRotation')[0]
                except: rotationsNode = etree.SubElement(self.xmlnode,'SelfRotation')
                
                self.parseRotationLine(strippedLine, rotationsNode)
    
                        
    def parsePeakLine(self, line, pattersonNode):
        if not line.strip().startswith('IX'):
            headings = 'No IX  IY  IZ  Xfrac  Yfrac  Zfrac  Xort   Yort    Zort   Dens Dens_sigma'.split()
            values = line.split()
            try:
                peakNode = etree.SubElement(pattersonNode,'Peak')
                for iHeading in range(len(headings)):
                    try:
                        self.subElementOfTypeWithText(peakNode,headings[iHeading],values[iHeading])
                    except:
                        pass
                #pattersonNode.append(peakNode)
            except:
                print('failed to parse line',line)

    def parseRotationLine(self, line, rotationNode):
        headings=  'No theta    phi     chi    alpha    beta   gamma      Rf    Rf_sigma'.split()
        values = line.split()
        try:
            peakNode = etree.SubElement(rotationNode,'Peak')
            for iHeading in range(len(headings)):
                try:
                    self.subElementOfTypeWithText(peakNode,headings[iHeading],values[iHeading+1])
                except:
                    pass
            #pattersonNode.append(peakNode)
        except:
            print('failed to parse line',line)    

    def subElementOfTypeWithText(self, parent=None, key=None, text=None):
        print(parent, key, text)
        newNode = etree.SubElement(parent, str(key))
        newNode.text = str(text)
        return newNode
