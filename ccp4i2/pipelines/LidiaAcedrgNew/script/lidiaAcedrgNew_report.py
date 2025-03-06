from __future__ import print_function

from report.CCP4ReportParser import *
import sys
#from lxml import etree
import xml.etree.ElementTree as etree
import math
from wrappers.acedrgNew.script.acedrgNew_report import acedrgNew_report
class lidiaAcedrgNew_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'LidiaAcedrgNew'
    RUNNING = True

    def __init__(self, *args, **kws):
        '''
        kws['additionalJsFiles'] = ['jsme/jsme.nocache.js']
        kws['additionalScript'] = 'function jsmeOnLoad() { jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px");'
        '''
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        # parent.addDiv(style="float:left;")
        parent.addDiv(style='clear:both;')
        smilesNodes = self.xmlnode.findall('.//SMILES')
        for smilesNode in smilesNodes:
            parent.addText(text='SMILES String: '+smilesNode.text)
        parent.addDiv(style='clear:both;')
        structureFold = parent.addFold(label='2D Structures',initiallyOpen=True)
        structureGallery = structureFold.addObjectGallery(style='float:left;border:1px solid black;',height='360px', tableWidth='260px', contentWidth='360px')
        clearingDiv = parent.addDiv(style="clear:both;")
        svgNodes = self.xmlnode.findall('.//Lidia/SVGNode')
        for svgNode in svgNodes:
            #Append the pretty printed svg for each SVGNode
            structureDiv = structureGallery.addDiv(label='From Lidia', title='From Lidia',style='width:355px;height:355px;border:0px solid white;')
            structureDiv.append(etree.tostring(svgNode[0]).decode("utf-8"))
        svgNodes = self.xmlnode.findall('.//Acedrg/SVGNode')
        for svgNode in svgNodes:
            #Append the pretty printed svg for each SVGNode
            structureDiv = structureGallery.addDiv(label='Interpreted by Acedrg/Rdkit', title='Interpreted by Acedrg/Rdkit', style='width:355px;height:355px;border:0px solid white;')
            structureDiv.append(etree.tostring(svgNode[0]).decode("utf-8"))
        '''for node in self.xmlnode.findall('.//Acedrg/Warning'):
            parent.addText(text=node.text)
            parent.append('<br/>')
        '''
        parent.addDiv(style='clear:both;')
#FIXME - XML PICTURE
        #self.addPictures(parent)
        acedrgNodes = self.xmlnode.findall('.//Acedrg')
        for acedrgNode in acedrgNodes:
            acedrgReport = acedrgNew_report(xmlnode=acedrgNode, jobStatus='nooutput')
            acedrgReport.analyseGeometry(parent)

    def addPictures(self, parent=None):
        #Note...this cnnot be moved into the acedrg_report because it uses this jobs jobInfo
        if parent is None: parent = self
        from core import CCP4Utils
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        import os
        baseScenePath = os.path.join(ccp4i2_root,'wrappers','acedrg','script','acedrg.scene.xml')
        if not 'TLC' in self.jobInfo['filenames']: return
        tlc = self.jobInfo['filenames']['TLC'].upper()
        if 'filenames' in self.jobInfo and 'XYZOUT_LIST' in self.jobInfo['filenames'] and len(self.jobInfo['filenames']['XYZOUT_LIST'])>0:
            pictureFold = parent.addFold(label='Pictures',initiallyOpen=False)
            pictureGallery = pictureFold.addObjectGallery(style='float:left;border:1px solid black;',height='400px', tableWidth='170px', contentWidth='450px')
            clearingDiv = parent.addDiv(style="clear:both;")
            print(self.jobInfo['filenames']['XYZOUT_LIST'])
            for i in range (len(self.jobInfo['filenames']['XYZOUT_LIST'])):
                pdbPath = self.jobInfo['filenames']['XYZOUT_LIST'][i]
                scenePath = os.path.join(self.jobInfo['fileroot'],'my.scene_'+str(i)+'.xml')
                print(pdbPath)
                print(scenePath)
                with open(baseScenePath,'r') as baseScene:
                    baseSceneText = baseScene.read()
                    specializedText = baseSceneText.replace('SUBSTITUTEME',pdbPath)
                    rootNode = etree.fromstring(specializedText)
                    
                    
                    if 'filenames' in self.jobInfo and 'DICTOUT_LIST' in self.jobInfo['filenames'] and len(self.jobInfo['filenames']['DICTOUT_LIST'])>0:
                        dictPath =  self.jobInfo['filenames']['DICTOUT_LIST'][0]
                        tlc = self.jobInfo['filenames']['TLC'].upper()
        
                        molDataNode = rootNode.findall('/scene/data/MolData')[0]
                        with open(dictPath) as dictF:
                           dictData = dictF.read()
                        """
                        customResNode = etree.fromstring('''<customResCIFFiles> <cifmonomer> <name>'''+tlc+'''</name> <filedata>'''+dictData+'''</filedata> </cifmonomer> </customResCIFFiles>''')
                        molDataNode.append(customResNode)
                        """
                    
                    with open(scenePath,'w') as specializedScene:
                        if sys.version_info > (3,0):
                            specializedScene.write(etree.tostring(rootNode).decode())
                        else:
                            specializedScene.write(etree.tostring(rootNode))
                    pic = pictureGallery.addPicture(label=os.path.split(pdbPath.__str__())[1],title=os.path.split(pdbPath.__str__())[1],sceneFile=scenePath)
        return
    
