"""
Copyright (C) 2015 STFC
"""

import os
import shutil
import subprocess as SP
import xml.etree.ElementTree as ET

from PySide2 import QtCore

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript


class molrep_pipe(CPluginScript):

    TASKNAME = 'molrep_pipe'
    MAINTAINER = 'liz.potterton@york.ac.uk'
    PERFORMANCECLASS = 'CRefinementPerformance'
    WHATNEXT = ['prosmart_refmac','modelcraft','coot_rebuild']
    ASYNCHRONOUS = True     # controls dynamic refmac table and graph only
    TIMEOUT_PERIOD = 240
    MAXNJOBS = 4
    ERROR_CODES = { 301 : { 'description' : 'Error reading program xml output from first molrep run' },
                    302 : { 'description' : 'No Laue results in program xml output from first molrep run' }
                    }
    PURGESEARCHLIST = [ [ 'molrep_mr%*/align.pdb' , 1],
                        [ 'molrep_mr%*/molrep_mtz.cif' , 1 ],
                        [ 'molrep_mr%*/molrep.doc.txt' , 5 ]
                        ]

    def process(self):
      self.dynrep = bool(self.container.inputData.DYNREP)
      if str(self.container.controlParameters.SG_OPTIONS) == 'specify':
        self.newspacegroup = str(self.container.controlParameters.SG)

      else:
        self.newspacegroup = str(self.container.inputData.F_SIGF.fileContent.spaceGroup)

      self.xmlroot = ET.Element('MolrepPipe')
      self.xmlroot.text = '\n'
      CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))

      if str(self.container.controlParameters.SG_OPTIONS) == 'laue':
        self.run1()

      else:
        self.run2()

      return CPluginScript.SUCCEEDED

    def run1(self):
      self.molrep1 = self.makePluginObject('molrep_mr')
      self.molrep1.container.inputData.copyData(self.container.inputData)
      self.molrep1.container.controlParameters.copyData(self.container.controlParameters)
      self.molrep1.container.guiParameters.copyData(self.container.guiParameters)
      self.connectSignal(self.molrep1,'finished', self.fin1)
      self.molrep1.process()

    def unique_tags(self, tree):
      cou = 0
      for e1 in tree.findall('RFpeaks'):
        e1.tag = 'RFpeaks' + str(cou)
        rf_list = list()
        for e2 in e1.findall('RFpeak'):
          rft = e2.findtext('RF')
          rf = float(rft)
          rf_list.append((rf, e2))

        f1 = ET.SubElement(tree, 'RFsorted' + str(cou))
        f1.text = e1.text
        f1.tail = e1.tail
        for rf, e2 in sorted(rf_list):
          f2 = ET.SubElement(f1, e2.tag)
          f2.text = e2.text
          f2.tail = e2.tail
          for e3 in e2:
            f3 = ET.SubElement(f2, e3.tag)
            f3.text = e3.text
            f3.tail = e3.tail

        cou += 1

    def appendxml(self, tree, tag):
      if len(self.xmlroot) and self.xmlroot[-1].tag == 'RefmacRunning':
        del self.xmlroot[-1]

      self.xmlroot.append(tree)
      tree.tail = '\n'
      tree.tag = tag
      CCP4Utils.writeXml(self.xmlroot, self.makeFileName('PROGRAMXML'))

    @QtCore.Slot(dict)
    def fin1(self, status):
      if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      xml = self.molrep1.makeFileName('PROGRAMXML')
      tree = ET.parse(xml).getroot()

      self.unique_tags(tree)

      laueE = tree.find('laue_group_alternatives')
      if laueE is None:
        self.appendErrorReport(302, xml)
        self.reportStatus(CPluginScript.FAILED)
        return

      laueData = []
      testList = laueE.findall('test')
      eleNames = ['space_group', 'score', 'contrast']
      for testE in testList:
        e1 = ET.SubElement(testE, 'selected')
        e1.text = '-'
        laueData.append({})
        for name in eleNames:
          e = testE.find(name)
          if name == 'space_group':
            laueData[-1][name] = str(e.text)
          else:
            laueData[-1][name] = float(e.text)

      best = 0
      for ind in range(1, len(laueData)):
        if (laueData[ind]['score'])>laueData[best]['score']: best = ind

      print('Proceeding for best space group', laueData[best]['space_group'])
      self.newspacegroup = laueData[best]['space_group']
      e1 = testList[best].find('selected')
      e1.text = 'yes'
      e0 = ET.SubElement(laueE, 'selected')
      e0.text = self.newspacegroup
      e0.tail = '\n'

      self.appendxml(tree, 'MolrepSpaceGroup')
      self.run2()

    def reindex(self, mtzin, mtzout):
      cmd = ('reindex', 'hklin', mtzin, 'hklout', mtzout)
      stdi = 'symm \'%s\'\nend\n' %self.newspacegroup
      sp = SP.Popen(cmd, stdin=SP.PIPE)
      sp.stdin.write(stdi.encode('ascii'))
      sp.stdin.close()
      return sp.wait()

    def run2(self):
      if self.newspacegroup == str(self.container.inputData.F_SIGF.fileContent.spaceGroup):
        self.fobs = str(self.container.inputData.F_SIGF)
        self.free = str(self.container.inputData.FREERFLAG)

      else:
        self.fobs = str(self.container.outputData.F_SIGF)
        self.free = str(self.container.outputData.FREERFLAG)
        ret1 = self.reindex(str(self.container.inputData.F_SIGF), self.fobs)
        ret2 = self.reindex(str(self.container.inputData.FREERFLAG), self.free)
        if ret1 or ret2:
          self.reportStatus(CPluginScript.FAILED)
          return
        else:
          self.container.outputData.F_SIGF.setContentFlag(reset=True)
          self.container.outputData.F_SIGF.annotation = 'Observed data reindexed to '+str(self.newspacegroup)
          self.container.outputData.FREERFLAG.annotation = 'FreeR reindexed to '+str(self.newspacegroup)

      self.molrep2 = self.makePluginObject('molrep_mr')
      self.molrep2.container.inputData.copyData(self.container.inputData)
      self.molrep2.container.inputData.F_SIGF.set(self.fobs)
      self.molrep2.container.controlParameters.copyData(self.container.controlParameters)
      self.molrep2.container.controlParameters.SG_OPTIONS = 'no'
      self.molrep2.container.guiParameters.copyData(self.container.guiParameters)
      self.connectSignal(self.molrep2, 'finished', self.fin2)
      self.molrep2.process()

    @QtCore.Slot(dict)
    def fin2(self, status):
      if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = ET.parse(self.molrep2.makeFileName('PROGRAMXML')).getroot()

      self.unique_tags(tree)
      self.appendxml(tree, 'MolrepSearch')

      self.xyz = str(self.molrep2.container.outputData.XYZOUT)
      if os.path.isfile(self.xyz):
        shutil.copyfile(str(self.molrep2.container.outputData.XYZOUT), str(self.container.outputData.XYZOUT_MOLREP))
        self.container.outputData.XYZOUT_MOLREP.annotation='Atomic model from molecular replacement'
        if self.container.inputData.RUNSHEETBEND:
          self.runSheetbend()

        else:
          self.run3()

      else:
        self.reportStatus(CPluginScript.SUCCEEDED)


    def runSheetbend(self):
      self.sheetbendPlugin = self.makePluginObject('sheetbend')
      self.sheetbendPlugin.container.inputData.XYZIN.set(self.xyz)
      self.sheetbendPlugin.container.inputData.F_SIGF.set(self.fobs)
      self.sheetbendPlugin.container.inputData.FREERFLAG.set(self.free)


      self.connectSignal(self.sheetbendPlugin, 'finished', self.finSheetbend)
      self.sheetbendPlugin.process()

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
      try:
        shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
        pipelineOutputItem.annotation = pluginOutputItem.annotation
        pipelineOutputItem.contentFlag = pluginOutputItem.contentFlag
        pipelineOutputItem.subType = pluginOutputItem.subType
      except:
        self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
        self.reportStatus(CPluginScript.FAILED)

    @QtCore.Slot(dict)
    def finSheetbend(self, status):
      if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = ET.parse(self.sheetbendPlugin.makeFileName('PROGRAMXML')).getroot()

      self.appendxml(tree, 'SheetbendResult')

      pluginOutputs=self.sheetbendPlugin.container.outputData
      pipelineOutputs = self.container.outputData
      self.harvestFile(pluginOutputs.XYZOUT, pipelineOutputs.XYZOUT_SHEETBEND)
      pipelineOutputs.XYZOUT_SHEETBEND.annotation='Atomic model after shift field refinement'

      self.xyz = str(self.sheetbendPlugin.container.outputData.XYZOUT)
      if os.path.isfile(self.xyz):
        self.run3()

      else:
        self.reportStatus(CPluginScript.SUCCEEDED)


    def run3(self):
      self.refmac = self.makePluginObject('refmac')
      self.refmac.container.inputData.XYZIN.set(self.xyz)
      self.refmac.container.inputData.F_SIGF.set(self.fobs)
      self.refmac.container.inputData.FREERFLAG.set(self.free)
      self.refmac.container.controlParameters.HYDROGENS = 'NO'
      self.refmac.container.controlParameters.NCYCLES = str(self.container.inputData.REFMAC_NCYC)
      self.refmac.container.controlParameters.PHOUT = False
      self.connectSignal(self.refmac, 'finished', self.fin3)

      if self.dynrep:
#     if self.ASYNCHRONOUS:
        self.refmac.doAsync = self.doAsync
        xml = self.refmac.makeFileName(format='PROGRAMXML')
        self.watchFile(xml, handler=self.handleXmlChanged, minDeltaSize=34, unwatchWhileHandling=True)

      self.refmac.process()

    @QtCore.Slot(dict)
    def handleXmlChanged(self, xmlFilename):
      tree = ET.parse(self.refmac.makeFileName('PROGRAMXML')).getroot()

      self.appendxml(tree, 'RefmacRunning')

    def fin3(self, status):
      if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = ET.parse(self.refmac.makeFileName('PROGRAMXML')).getroot()

      self.appendxml(tree, 'Refmac')

      try:
        shutil.copyfile(str(self.refmac.container.outputData.XYZOUT), str(self.container.outputData.XYZOUT))
        shutil.copyfile(str(self.refmac.container.outputData.DIFFPHIOUT), str(self.container.outputData.DIFFPHIOUT))
        shutil.copyfile(str(self.refmac.container.outputData.FPHIOUT), str(self.container.outputData.FPHIOUT))
        self.container.outputData.XYZOUT.annotation='Atomic model after restrained refinement'
        cycles = list()
        for e3 in tree.findall('./Overall_stats/stats_vs_cycle/new_cycle'):
          cycles.append((int(e3.findtext('cycle')), float(e3.findtext('r_factor')), float(e3.findtext('r_free'))))

        cycle_no, Rcrist, Rfree = sorted(cycles)[-1]
        self.container.outputData.PERFORMANCE.RFactor.set(Rcrist)
        self.container.outputData.PERFORMANCE.RFree.set(Rfree)
        self.container.outputData.FPHIOUT.annotation.set("Weighted map from restrained refinement")
        self.container.outputData.DIFFPHIOUT.annotation.set("Weighted difference map from restrained refinement")
        self.reportStatus(CPluginScript.SUCCEEDED)

      except:
        self.reportStatus(CPluginScript.FAILED)
