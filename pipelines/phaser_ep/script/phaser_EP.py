import os
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils


class phaser_EP(CPluginScript):

    TASKNAME = 'phaser_EP'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    ERROR_CODES = {
        201: {'description': 'ShelxCD sub-job failed'},
        202: {'description': 'Failed in harvest operation'},
        203: {'description': 'phaser_EP_AUTO sub-job failed'},
        204: {'description': 'Exception in ShelxCD setup'},
        205: {'description': 'Exception in phaser_EP_AUTO setup'},
        206: {'description': 'Exception in Parrot setup'},
        207: {'description': 'Exception in ModelCraft setup'},
        208: {'description': 'Exception in XML update'},
        209: {'description': 'Parrot sub-job failed'},
        210: {'description': 'ModelCraft sub-job failed'},
    }

    def process(self):
        self.xmlroot = etree.Element('PhaserEP')
        if self.container.inputData.PARTIALMODELORMAP == 'SEARCH':
            rv = self.run_shelx()
            if rv == CPluginScript.SUCCEEDED:
                self.updateXmlFromFile(self.shelxPlugin.makeFileName('PROGRAMXML'),'ShelxCD')
        elif self.container.inputData.PARTIALMODELORMAP in ['MODEL', 'MAP']:
            self.container.keywords.HAND.set('off')
        rv = self.run_phaser()
        if rv == CPluginScript.SUCCEEDED:
            self.updateXmlFromFile(self.phaserPlugin.makeFileName('PROGRAMXML'),'PhaserEpResults')
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
                    self.updateXmlFromFile(self.parrotPluginOriginalHand.makeFileName('PROGRAMXML'),'ParrotResult', hand='original')
                    self.copyPluginOutput(self.parrotPluginOriginalHand.container.outputData.ABCDOUT, pipelineOutputs.ABCDOUT, annotation='Phases from density modification - original hand')
                    self.copyPluginOutput(self.parrotPluginOriginalHand.container.outputData.FPHIOUT, pipelineOutputs.FPHIOUT, annotation='Map coefficients from density modification - original hand')
            if self.container.keywords.HAND in ['both', 'on']:
                if self.container.keywords.HAND == 'both':
                    self.parrotPluginInvertedHand = self.makeParrotPlugin(hand='inverted')
                else:
                    self.parrotPluginInvertedHand = self.makeParrotPlugin(hand='inverted', inverted_only=True)
                parrot_inverted = self.parrotPluginInvertedHand.process()
                if parrot_inverted == CPluginScript.SUCCEEDED:
                    self.updateXmlFromFile(self.parrotPluginInvertedHand.makeFileName('PROGRAMXML'),'ParrotResult', hand='inverted')
                    self.copyPluginOutput(self.parrotPluginInvertedHand.container.outputData.ABCDOUT, pipelineOutputs.ABCDOUT, annotation='Phases from density modification - reversed hand')
                    self.copyPluginOutput(self.parrotPluginInvertedHand.container.outputData.FPHIOUT, pipelineOutputs.FPHIOUT, annotation='Map coefficients from density modification - reversed hand')

        if self.container.controlParameters.RUNBUCCANEER:
            if self.container.keywords.HAND in ['both', 'off']:
                modelCraftOriginal = self.makeModelCraftPlugin(hand='original')
                if modelCraftOriginal.process() == CPluginScript.SUCCEEDED:
                    element = etree.Element("ModelCraft")
                    element.text = modelCraftOriginal.getWorkDirectory()
                    self.updateXml(element, 'ModelCraft', hand='original')
                    self.copyPluginOutput(modelCraftOriginal.container.outputData.XYZOUT, pipelineOutputs.XYZOUT, annotation='Autobuilt model - original hand')
            if self.container.keywords.HAND in ['both', 'on']:
                modelCraftInverted = self.makeModelCraftPlugin(hand='inverted')
                if modelCraftInverted.process() == CPluginScript.SUCCEEDED:
                    element = etree.Element("ModelCraft")
                    element.text = modelCraftInverted.getWorkDirectory()
                    self.updateXml(element, 'ModelCraft', hand='inverted')
                    self.copyPluginOutput(modelCraftInverted.container.outputData.XYZOUT, pipelineOutputs.XYZOUT, annotation='Autobuilt model - reversed hand')
        self.reportStatus(CPluginScript.SUCCEEDED)

    def run_shelx(self):
        try:
            self.shelxPlugin = self.makePluginObject('ShelxCD')
            self.shelxPlugin.container.inputData.SAD.set(self.container.inputData.F_SIGF)
            self.shelxPlugin.container.controlParameters.MODE == 'SAD'
            self.shelxPlugin.container.controlParameters.SFAC = self.container.inputData.SFAC
            self.shelxPlugin.container.controlParameters.NTRY = self.container.inputData.NTRY
            self.shelxPlugin.container.controlParameters.FIND = self.container.inputData.FIND
            self.shelxPlugin.doAsync = False
        except Exception as e:
            self.appendErrorReport(204, 'Exception during ShelxCD setup: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        try:
            rv = self.shelxPlugin.process()
        except Exception as e:
            self.appendErrorReport(201, 'Exception during ShelxCD process: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        if rv == CPluginScript.FAILED:
            self.appendErrorReport(201, 'ShelxCD sub-job failed in phaser_EP pipeline')
            self.reportStatus(rv)
            return CPluginScript.FAILED
        else:
            try:
                self.container.inputData.XYZIN_HA.set(self.shelxPlugin.container.outputData.XYZOUT)
                self.container.inputData.PARTIALMODELORMAP.set('NONE')
            except Exception as e:
                self.appendErrorReport(204, 'Exception copying ShelxCD output: ' + str(e))
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            return rv

    def run_phaser(self):
        try:
            self.phaserPlugin = self.makePluginObject('phaser_EP_AUTO')
            # see comment in pipelines/phaser_pipeline/script/phaser_pipeline.py
            for attrName in self.phaserPlugin.container.keywords.dataOrder():
                if hasattr(self.container.keywords,attrName):
                    attr = getattr(self.container.keywords,attrName)
                    if hasattr(attr,'isSet') and attr.isSet(allSet=False):
                        setattr(self.phaserPlugin.container.keywords,attrName,attr)
            self.phaserPlugin.container.inputData=self.container.inputData
            self.phaserPlugin.doAsync = False
        except Exception as e:
            self.appendErrorReport(205, 'Exception during phaser_EP_AUTO setup: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        try:
            rv = self.phaserPlugin.process()
        except Exception as e:
            self.appendErrorReport(203, 'Exception during phaser_EP_AUTO process: ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        if rv == CPluginScript.FAILED:
            self.appendErrorReport(203, 'phaser_EP_AUTO sub-job failed in phaser_EP pipeline')
            self.reportStatus(rv)
            return CPluginScript.FAILED
        else:
            return rv

    def makeParrotPlugin(self, hand, inverted_only=False):
        try:
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
        except Exception as e:
            self.appendErrorReport(206, 'Exception during Parrot setup (hand=' + hand + '): ' + str(e))
            raise

    def makeModelCraftPlugin(self, hand):
        try:
            plugin = self.makePluginObject('modelcraft')
            plugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
            plugin.container.inputData.FREERFLAG.set(self.container.inputData.FREERFLAG)
            plugin.container.inputData.ASUIN.set(self.container.inputData.ASUFILE)
            if hand == 'original':
                plugin.container.inputData.PHASES.set(self.parrotPluginOriginalHand.container.outputData.ABCDOUT)
            elif hand == 'inverted':
                plugin.container.inputData.PHASES.set(self.parrotPluginInvertedHand.container.outputData.ABCDOUT)
            plugin.container.controlParameters.BASIC.set(True)
            plugin.container.controlParameters.USE_MODEL_PHASES.set(False)
            plugin.container.controlParameters.CYCLES.set(self.container.controlParameters.BUCCANEER_ITERATIONS)
            plugin.container.controlParameters.STOP_CYCLES.set(2)
            return plugin
        except Exception as e:
            self.appendErrorReport(207, 'Exception during ModelCraft setup (hand=' + hand + '): ' + str(e))
            raise

    def updateXmlFromFile(self, xmlFilename, element, hand=None):
        try:
            newXML = CCP4Utils.openFileToEtree(xmlFilename)
            self.updateXml(newXML, element, hand)
        except Exception as e:
            self.appendErrorReport(208, 'Exception in updateXmlFromFile (element=' + element + '): ' + str(e))
            # Don't fail pipeline for XML update issues - log and continue

    def updateXml(self, newXML, element, hand=None):
        try:
            if hand is None:
                for oldNode in self.xmlroot.xpath(element):
                    self.xmlroot.remove(oldNode)
                self.xmlroot.append(newXML)
            else:
                if len(self.xmlroot.xpath('//{}/{}'.format(hand,element))) > 0:
                    self.xmlroot.remove(self.xmlroot.xpath('//{}/{}'.format(hand,element))[0])
                hand_node = etree.SubElement(self.xmlroot, hand)
                hand_node.append(newXML)
            tmpFilename = self.makeFileName('PROGRAMXML')+'_tmp'
            with open(tmpFilename,'w') as xmlfile:
                CCP4Utils.writeXML(xmlfile,etree.tostring(self.xmlroot,pretty_print=True))
            self.renameFile(tmpFilename,self.makeFileName('PROGRAMXML'))
        except Exception as e:
            self.appendErrorReport(208, 'Exception in updateXml (element=' + element + '): ' + str(e))
            # Don't fail pipeline for XML update issues - log and continue

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
        except Exception as e:
            self.appendErrorReport(202, 'Failed to copy ' + str(pluginOutputItem.fullPath) + ' to ' + str(pipelineOutputItem.fullPath) + ': ' + str(e))
            self.reportStatus(CPluginScript.FAILED)
