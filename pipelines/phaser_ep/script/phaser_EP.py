from __future__ import print_function

from core.CCP4PluginScript import CPluginScript
import sys, os
from core import CCP4ErrorHandling
from core import CCP4Modules
from lxml import etree
from core import CCP4Utils


class phaser_EP(CPluginScript):

    TASKNAME = 'phaser_EP'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    ERROR_CODES = {202:{'description':'Failed in harvest operation'},}

    def process(self):
        self.xmlroot = etree.Element('PhaserEP')
        if self.container.inputData.PARTIALMODELORMAP == 'SEARCH':
            rv = self.run_shelx()
            if rv == CPluginScript.SUCCEEDED:
                self.updateXml(self.shelxPlugin.makeFileName('PROGRAMXML'),'ShelxCD')
        elif self.container.inputData.PARTIALMODELORMAP in ['MODEL', 'MAP']:
            self.container.keywords.HAND.set('off')
        rv = self.run_phaser()
        if rv == CPluginScript.SUCCEEDED:
            self.updateXml(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserEpResults')
            pipelineOutputs = self.container.outputData
            pluginOutputs=self.phaserPlugin.container.outputData
            hands = [0, 1] if self.container.keywords.HAND == 'both' else [0]
            for i in hands:
                self.copyPluginOutput(pluginOutputs.XYZOUT[i], pipelineOutputs.XYZOUT)
                self.copyPluginOutput(pluginOutputs.ABCDOUT[i], pipelineOutputs.ABCDOUT)
                self.copyPluginOutput(pluginOutputs.MAPOUT[i], pipelineOutputs.MAPOUT)

        if self.container.controlParameters.RUNPARROT:
            if self.container.keywords.HAND in ['both', 'off']:
                self.parrotPluginOriginalHand = self.makeParrotPlugin(hand='original')
                parrot_original = self.parrotPluginOriginalHand.process()
                if parrot_original == CPluginScript.SUCCEEDED:
                    self.updateXml(self.parrotPluginOriginalHand.makeFileName('PROGRAMXML'),'ParrotResult', hand='original')
                    self.copyPluginOutput(self.parrotPluginOriginalHand.container.outputData.ABCDOUT, pipelineOutputs.ABCDOUT, annotation='Phases from density modification - original hand')
                    self.copyPluginOutput(self.parrotPluginOriginalHand.container.outputData.FPHIOUT, pipelineOutputs.FPHIOUT, annotation='Map coefficients from density modification - original hand')
            if self.container.keywords.HAND in ['both', 'on']:
                if self.container.keywords.HAND == 'both':
                    self.parrotPluginInvertedHand = self.makeParrotPlugin(hand='inverted')
                else:
                    self.parrotPluginInvertedHand = self.makeParrotPlugin(hand='inverted', inverted_only=True)
                parrot_inverted = self.parrotPluginInvertedHand.process()
                if parrot_inverted == CPluginScript.SUCCEEDED:
                    self.updateXml(self.parrotPluginInvertedHand.makeFileName('PROGRAMXML'),'ParrotResult', hand='inverted')
                    self.copyPluginOutput(self.parrotPluginInvertedHand.container.outputData.ABCDOUT, pipelineOutputs.ABCDOUT, annotation='Phases from density modification - reversed hand')
                    self.copyPluginOutput(self.parrotPluginInvertedHand.container.outputData.FPHIOUT, pipelineOutputs.FPHIOUT, annotation='Map coefficients from density modification - reversed hand')

        if self.container.controlParameters.RUNBUCCANEER:
            if self.container.keywords.HAND in ['both', 'off']:
                self.buccaneerPluginOriginalHand = self.makeBuccaneerPlugin(hand='original')
                buccaneer_original = self.buccaneerPluginOriginalHand.process()
                if buccaneer_original == CPluginScript.SUCCEEDED:
                    self.updateXml(self.buccaneerPluginOriginalHand.makeFileName('PROGRAMXML'),'BuccaneerBuildRefineResult', hand='original')
                    self.copyPluginOutput(self.buccaneerPluginOriginalHand.container.outputData.XYZOUT, pipelineOutputs.XYZOUT, annotation='Model built by Autobuild protein - original hand')
            if self.container.keywords.HAND in ['both', 'on']:
                self.buccaneerPluginInvertedHand = self.makeBuccaneerPlugin(hand='inverted')
                buccaneer_inverted = self.buccaneerPluginInvertedHand.process()
                if buccaneer_inverted == CPluginScript.SUCCEEDED:
                    self.updateXml(self.buccaneerPluginInvertedHand.makeFileName('PROGRAMXML'),'BuccaneerBuildRefineResult', hand='inverted')
                    self.copyPluginOutput(self.buccaneerPluginInvertedHand.container.outputData.XYZOUT, pipelineOutputs.XYZOUT, annotation='Model built by Autobuild protein - reversed hand')
        self.reportStatus(CPluginScript.SUCCEEDED)


    def run_shelx(self):
        self.shelxPlugin = self.makePluginObject('ShelxCD')
        self.shelxPlugin.container.inputData.SAD.set(self.container.inputData.F_SIGF)
        self.shelxPlugin.container.controlParameters.MODE == 'SAD'
        self.shelxPlugin.container.controlParameters.SFAC = self.container.inputData.SFAC
        self.shelxPlugin.container.controlParameters.NTRY = self.container.inputData.NTRY
        self.shelxPlugin.container.controlParameters.FIND = self.container.inputData.FIND
        self.shelxPlugin.doAsync = False
        rv = self.shelxPlugin.process()
        if rv == CPluginScript.FAILED:
            self.reportStatus(rv)
            return CPluginScript.FAILED
        else:
            self.container.inputData.XYZIN_HA.set(self.shelxPlugin.container.outputData.XYZOUT)
            self.container.inputData.PARTIALMODELORMAP.set('NONE')
            return rv

    def run_phaser(self):
        self.phaserPlugin = self.makePluginObject('phaser_EP_AUTO')
        # see comment in pipelines/phaser_pipeline/script/phaser_pipeline.py
        for attrName in self.phaserPlugin.container.keywords.dataOrder():
            if hasattr(self.container.keywords,attrName):
                attr = getattr(self.container.keywords,attrName)
                if hasattr(attr,'isSet') and attr.isSet(allSet=False):
                    setattr(self.phaserPlugin.container.keywords,attrName,attr)
        self.phaserPlugin.container.inputData=self.container.inputData
        self.phaserPlugin.doAsync = False
        rv = self.phaserPlugin.process()
        if rv == CPluginScript.FAILED:
            self.reportStatus(rv)
            return CPluginScript.FAILED
        else:
            return rv

    def makeParrotPlugin(self, hand, inverted_only=False):
        self.parrotPlugin =  self.makePluginObject('parrot')
        self.parrotPlugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
        self.parrotPlugin.container.inputData.ASUIN.set(self.container.inputData.ASUFILE)
        if hand == 'original':
            self.parrotPlugin.container.inputData.ABCD.set(self.phaserPlugin.container.outputData.ABCDOUT[0])
        elif hand == 'inverted':
            if not inverted_only:
                self.parrotPlugin.container.inputData.ABCD.set(self.phaserPlugin.container.outputData.ABCDOUT[1])
            else:
                self.parrotPlugin.container.inputData.ABCD.set(self.phaserPlugin.container.outputData.ABCDOUT[0])
        self.parrotPlugin.doAsync = False
        return self.parrotPlugin

    def makeBuccaneerPlugin(self, hand):
        self.buccaneerPlugin =  self.makePluginObject('buccaneer_build_refine_mr')
        self.buccaneerPlugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
        self.buccaneerPlugin.container.inputData.FREERFLAG.set(self.container.inputData.FREERFLAG)
        self.buccaneerPlugin.container.inputData.ASUIN.set(self.container.inputData.ASUFILE)
        self.buccaneerPlugin.container.controlParameters.BUCCANEER_PHSIN_TYPE.set('experimental')
        self.buccaneerPlugin.container.controlParameters.ITERATIONS.set(self.container.controlParameters.BUCCANEER_ITERATIONS)

        if hand == 'original':
            self.buccaneerPlugin.container.inputData.ABCD.set(self.parrotPluginOriginalHand.container.outputData.ABCDOUT)
        elif hand == 'inverted':
            self.buccaneerPlugin.container.inputData.ABCD.set(self.parrotPluginInvertedHand.container.outputData.ABCDOUT)
        self.buccaneerPlugin.doAsync = False
        return self.buccaneerPlugin

    def updateXml(self, xmlFilename, element, hand=None):
        if hand == None:
            for oldNode in self.xmlroot.xpath(element):
                self.xmlroot.remove(oldNode)
            newXML = CCP4Utils.openFileToEtree(xmlFilename)
            self.xmlroot.append(newXML)
        else:
            if len(self.xmlroot.xpath('//{}/{}'.format(hand,element))) > 0:
                self.xmlroot.remove(self.xmlroot.xpath('//{}/{}'.format(hand,element))[0])
            hand_node = etree.SubElement(self.xmlroot, hand)
            newXML = CCP4Utils.openFileToEtree(xmlFilename)
            hand_node.append(newXML)
        tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
        with open(tmpFilename,'w') as xmlfile:
            CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))
        self.renameFile(tmpFilename,self.makeFileName('PROGRAMXML'))

    def copyPluginOutput(self, pluginOutputItem, pipelineOutputList, annotation=None):
        pipelineOutputList.append(pipelineOutputList.makeItem())
        pipelineOutputItem = pipelineOutputList[-1]
        pipelineOutputItem.fullPath = os.path.join(self.workDirectory,os.path.basename(str(pluginOutputItem.fullPath)))
        try:
            import shutil
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            if annotation is None:
                pipelineOutputItem.annotation = pluginOutputItem.annotation
            else: 
                pipelineOutputItem.annotation = annotation
            pipelineOutputItem.contentFlag = pluginOutputItem.contentFlag
            pipelineOutputItem.subType = pluginOutputItem.subType
        except:
            self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
            self.reportStatus(CPluginScript.FAILED)
