"""
TODO

* If 2 enantiomers don't bother with second one if first one gives R < 30%.

"""

import copy
import json
import os
import shutil
import sys

from lxml import etree

from ccp4i2.baselayer.QtCore import Slot
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4ErrorHandling import CException
from ccp4i2.core.CCP4PluginScript import CPluginScript


class dr_mr_modelbuild_pipeline(CPluginScript):

    TASKNAME = 'dr_mr_modelbuild_pipeline'
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    PERFORMANCECLASS = 'CRefinementPerformance'
    WHATNEXT = []
    ERROR_CODES = { 301 : { 'description' : 'Error reading program xml output from first molrep run' },
                    302 : { 'description' : 'No Laue results in program xml output from first molrep run' },
                    202 : { 'description' : 'Failed in harvesting file from lidia/acedrg' }
                    }
    PURGESEARCHLIST = [ [ 'molrep_mr%*/align.pdb' , 1],
                        [ 'molrep_mr%*/molrep_mtz.cif' , 1 ],
                        [ 'molrep_mr%*/molrep.doc.txt' , 5 ]
                        ]

    def validity(self):
        """Filter CSMILESString validation errors when mode is not SMILES."""
        from ccp4i2.core import CCP4ErrorHandling
        error = super(dr_mr_modelbuild_pipeline, self).validity()
        mode = str(self.container.controlParameters.LIGANDAS) if self.container.controlParameters.LIGANDAS.isSet() else ""
        if mode != 'SMILES':
            # Filter out SMILES validation errors when not using SMILES input
            filtered = CCP4ErrorHandling.CErrorReport()
            for err in error.getErrors():
                err_class = err.get('class', '')
                err_name = err.get('name', '')
                if err_class == 'CSMILESString' and 'SMILESIN' in err_name:
                    continue
                filtered.append(
                    err.get('class', ''),
                    err.get('code', 0),
                    err.get('details', ''),
                    err.get('name', ''),
                    err.get('severity', 0)
                )
            return filtered
        return error

    def process(self):
      self.runningJobs=[]
      self.newspacegroup = str(self.container.inputData.F_SIGF.fileContent.spaceGroup)

      self.xmlroot = etree.Element('CCP4i2DRMRMBPipe')
      self.xmlroot.text = '\n'
      with open(str(self.makeFileName('PROGRAMXML')), 'w') as ostream:
        CCP4Utils.writeXML(ostream,etree.tostring(self.xmlroot,pretty_print=True))

      if self.container.controlParameters.LIGANDAS.__str__() == 'DICT':
            self.dictToUse = self.container.inputData.DICTIN
            self.dictDone()
      elif self.container.controlParameters.LIGANDAS.__str__() == 'NONE':
            self.dictDone()
      elif self.container.controlParameters.LIGANDAS.__str__() != 'NONE':
            self.lidiaAcedrgPlugin = self.makePluginObject('LidiaAcedrgNew')
            self.lidiaAcedrgPlugin.container.inputData.MOLSMILESORSKETCH.set(self.container.controlParameters.LIGANDAS)
            if self.container.inputData.MOLIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.MOLIN.set(self.container.inputData.MOLIN)
            if self.container.inputData.SMILESIN.isSet():
                self.lidiaAcedrgPlugin.container.inputData.SMILESIN.set(self.container.inputData.SMILESIN)
            self.lidiaAcedrgPlugin.container.inputData.CONFORMERSFROM.set('RDKIT')
            self.lidiaAcedrgPlugin.container.inputData.TLC.set('DRG')
            self.connectSignal(self.lidiaAcedrgPlugin,'finished',self.lidiaAcedrg_finished)
            self.lidiaAcedrgPlugin.process()

    @Slot(dict)
    def lidiaAcedrg_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)
        pluginRoot = CCP4Utils.openFileToEtree(self.lidiaAcedrgPlugin.makeFileName('PROGRAMXML'))
        self.xmlroot.append(pluginRoot)
        self.flushXML()
        self.harvestFile(self.lidiaAcedrgPlugin.container.outputData.DICTOUT_LIST[0], self.container.outputData.DICTOUT)
        self.dictToUse = self.container.outputData.DICTOUT
        self.dictDone()

    def harvestFile(self, pluginOutputItem, pipelineOutputItem):
        try:
            shutil.copyfile(str(pluginOutputItem.fullPath), str(pipelineOutputItem.fullPath))
            pipelineOutputItem.annotation.set(pluginOutputItem.annotation)
            pipelineOutputItem.contentFlag.set(pluginOutputItem.contentFlag)
            #print '#harvestFile',pluginOutputItem.fullPath, pluginOutputItem.contentFlag
            pipelineOutputItem.subType.set(pluginOutputItem.subType)
        except:
            self.appendErrorReport(202,str(pluginOutputItem.fullPath)+' '+str(pipelineOutputItem.fullPath))
            self.finishWithStatus(CPluginScript.FAILED)
    
    def finishWithStatus(self, status=CPluginScript.SUCCEEDED):
        self.flushXML()
        self.reportStatus(status)

    def dictDone(self):
      if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
          self.aimlessPipe()
      elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
          self.importMergedPipe()
      else:
          self.aimlessPlugin_finished({"finishStatus":CPluginScript.SUCCEEDED})

      return CPluginScript.SUCCEEDED

    def importMergedPipe(self):
        self.aimlessPlugin = self.makePluginObject('import_merged')
        self.aimlessPlugin.container.inputData.HKLIN.set(self.container.inputData.HKLIN)
        self.connectSignal(self.aimlessPlugin,'finished',self.aimlessPlugin_finished)
        self.aimlessPlugin.process()

    def aimlessPipe(self):
        self.aimlessPlugin = self.makePluginObject('aimless_pipe')
        self.aimlessPlugin.container.controlParameters.AUTOCUTOFF.set(self.container.controlParameters.AUTOCUTOFF)
        self.aimlessPlugin.container.controlParameters.SCALING_PROTOCOL.set('DEFAULT')
        self.aimlessPlugin.container.controlParameters.ONLYMERGE.set(False)
        self.aimlessPlugin.container.inputData.copyData(self.container.inputData,['UNMERGEDFILES'])
        self.connectSignal(self.aimlessPlugin,'finished',self.aimlessPlugin_finished)
        self.aimlessPlugin.process()

    @Slot(dict)
    def aimlessPlugin_finished(self, status):
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(CPluginScript.FAILED)

        if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
            pluginRoot = CCP4Utils.openFileToEtree(self.aimlessPlugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(pluginRoot)
        self.flushXML()

        if self.container.inputData.XYZINORMRBUMP.__str__() == 'MRBUMP':
            """
            Run MrBUMP and then call self.aimlessCyclesFinished with MrBUMP output .
            """
            try:
                print("Create MrBUMP model prep plugin")
                self.mrbump = self.makePluginObject('mrbump_model_prep')
                print("Set MrBUMP model prep plugin option ASUIN")
                self.mrbump.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
                print("Set MrBUMP model prep plugin option F_SIGF")
                if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
                    self.mrbump.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
                    print("Set MrBUMP model prep plugin option FREERFLAG")
                    self.mrbump.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
                    print("Set MrBUMP model prep plugin option MRMAX")
                elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
                    self.mrbump.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
                    print("Set MrBUMP model prep plugin option FREERFLAG")
                    self.mrbump.container.inputData.FREERFLAG.set(self.aimlessPlugin.freerflag.container.outputData.FREEROUT)
                    print("Set MrBUMP model prep plugin option MRMAX")
                else:
                    self.mrbump.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
                    print("Set MrBUMP model prep plugin option FREERFLAG")
                    self.mrbump.container.inputData.FREERFLAG.set(self.container.inputData.FREER_IN)
                    print("Set MrBUMP model prep plugin option MRMAX")
                self.mrbump.container.inputData.MRMAX.set(self.container.inputData.MRMAX)
                print("Set MrBUMP model prep plugin option REDUNDANCYLEVEL")
                self.mrbump.container.inputData.REDUNDANCYLEVEL.set(self.container.inputData.REDUNDANCYLEVEL)
                print("Set MrBUMP model prep plugin option AFDBLEVEL")
                self.mrbump.container.inputData.AFDBLEVEL.set(self.container.inputData.AFDBLEVEL)
                print("Set MrBUMP model prep plugin option SEARCH_PDB")
                self.mrbump.container.inputData.SEARCH_PDB.set(self.container.inputData.SEARCH_PDB)
                print("Set MrBUMP model prep plugin option SEARCH_AFDB")
                self.mrbump.container.inputData.SEARCH_AFDB.set(self.container.inputData.SEARCH_AFDB)
                print("Run MrBUMP model prep plugin")
                rv = self.mrbump.process()
                if rv == CPluginScript.SUCCEEDED:
                    print("MrBUMP model prep plugin succeeded. Copying data ...")
                    self.XYZIN = str(self.mrbump.container.outputData.XYZOUT)
                    print("Move on to MR...")
                    pluginRoot = CCP4Utils.openFileToEtree(self.mrbump.makeFileName('PROGRAMXML'))
                    self.xmlroot.append(pluginRoot)
                    self.flushXML()
                    self.aimlessCyclesFinished()
                else:
                    self.reportStatus(CPluginScript.FAILED)
                    return CPluginScript.FAILED
            except:
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
        elif self.container.inputData.XYZINORMRBUMP.__str__() == 'MRPARSE':
            try:
                print("Create MrParse model prep plugin")
                self.mrparse = self.makePluginObject('mrparse_simple')
                print("Created MrParse model prep plugin")
                seqin = os.path.join(self.workDirectory,'SEQIN.fasta')
                self.container.inputData.ASUIN.writeFasta(seqin)
                self.mrparse.container.inputData.SEQIN.set(seqin)
                self.mrparse.container.options.DATABASE.unSet()
                self.mrparse.container.options.MAXHITS.unSet()
                self.mrparse.container.options.NPROC.unSet()
                rv = self.mrparse.process()
                print("MrParse rv",rv); sys.stdout.flush()
                if rv == CPluginScript.SUCCEEDED:
                    print("MrParse model prep plugin succeeded. Copying data ..."); sys.stdout.flush()
                    self.XYZIN = str(self.mrparse.container.outputData.XYZOUT[0])
                    print("Move on to MR...")
                    self.aimlessCyclesFinished()
                else:
                    print("MrParse model prep plugin FAILED"); sys.stdout.flush()
                    self.reportStatus(CPluginScript.FAILED)
                    return CPluginScript.FAILED
            except:
                self.reportStatus(CPluginScript.FAILED)
        else:
            self.XYZIN = self.container.inputData.XYZIN
            self.aimlessCyclesFinished()

    def aimlessCyclesFinished(self):
#At this point we check for enantiomorphs

        self.enantio = False
        self.enantiomorphs = []

        if len(self.xmlroot.xpath('/CCP4i2DRMRMBPipe/AIMLESS_PIPE/POINTLESS/SolutionWarning')) > 0:
            for el in self.xmlroot.xpath('/CCP4i2DRMRMBPipe/AIMLESS_PIPE/POINTLESS/SolutionWarning'):
                if "enantiomorph".lower() in el.text.lower():
                    self.enantio = True
                    break

        if self.enantio:
            best_sg = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/AIMLESS_PIPE/POINTLESS/BestSolution/GroupName')[0].text
            best_sg_prob = float(self.xmlroot.xpath('/CCP4i2DRMRMBPipe/AIMLESS_PIPE/POINTLESS/BestSolution/TotalProb')[0].text)
            sglist = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/AIMLESS_PIPE/POINTLESS/SpacegroupList/Spacegroup')
            for sg in sglist:
                sgname, total_prob = (sg.xpath('SpacegroupName')[0].text.strip(), float(sg.xpath('TotalProb')[0].text))
                if abs(total_prob-best_sg_prob)<1e-5:
                    self.enantiomorphs.append(sgname)
            #print("possible enantiomorphic space groups:",self.enantiomorphs)
            self.otherSG = list(set(self.enantiomorphs) - set([best_sg]))[0]
            #print("'Other' space group is",self.otherSG)


        self.useOtherSG = False

        if self.enantio and len(self.enantiomorphs)==2:
            self.molrep_job = self.makePluginObject('molrep_mr')
            self.molrep_job.container.inputData.XYZIN.set(self.XYZIN)
            self.molrep_job.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
            self.molrep_job.container.controlParameters.SG_OPTIONS.set('no')
            self.molrep_job.container.controlParameters.NMON.set(self.container.controlParameters.NMON)
            if self.container.inputData.XYZINORMRBUMP.__str__() == 'XYZINPUT':
                self.molrep_job.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
            rv = self.molrep_job.process()
            rvfin = self.fin_molrep(rv,self.molrep_job)

            print("Let's do it again with enantiomorph .....")

            self.reindexJob  = self.makePluginObject('pointless_reindexToMatch')
            self.reindexJob.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
            self.reindexJob.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
            self.reindexJob.container.controlParameters.USE_REINDEX.set(False)
            self.reindexJob.container.controlParameters.REFERENCE.set('SPECIFY')
            self.reindexJob.container.controlParameters.CHOOSE_SPACEGROUP.set(self.otherSG)
            rv = self.reindexJob.process()
            self.useOtherSG = True

            self.molrep_job2 = self.makePluginObject('molrep_mr')
            self.molrep_job2.container.inputData.XYZIN.set(self.XYZIN)
            self.molrep_job2.container.inputData.F_SIGF.set(self.reindexJob.container.outputData.F_SIGF_OUT)
            self.molrep_job2.container.controlParameters.SG_OPTIONS.set('specify')
            self.molrep_job2.container.controlParameters.SG.set(self.otherSG)

            self.molrep_job2.container.controlParameters.NMON.set(self.container.controlParameters.NMON)
            if self.container.inputData.XYZINORMRBUMP.__str__() == 'XYZINPUT':
                self.molrep_job2.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
            rv = self.molrep_job2.process()
            rvfin2 = self.fin_molrep(rv,self.molrep_job2)

            print("Now I have to work out what is best")

            try:
                # TODO: Will break without the old Buccaneer pipeline
                print("Get rfactors")
                rfactors = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/BuccaneerBuildRefineResult/FinalStatistics/r_factor')
                print("Get rfrees")
                rfrees = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/BuccaneerBuildRefineResult/FinalStatistics/r_free')
                r0 = rfactors[0].text
                r1 = rfactors[1].text
                rf0 = rfrees[0].text
                rf1 = rfrees[1].text

                if float(r0) < float(r1):

                    print("Copy FREER")
                    shutil.copyfile(str(self.aimlessPlugin.container.outputData.FREEROUT), str(self.container.outputData.FREEROUT))
                    print("Copy HKLOUT")
                    shutil.copyfile(str(self.aimlessPlugin.container.outputData.HKLOUT[0]), str(self.container.outputData.HKLOUT))
                    print("Copy FREEROUT_OTHERSG")
                    shutil.copyfile(str(self.reindexJob.container.outputData.FREERFLAG_OUT), str(self.container.outputData.FREEROUT_OTHERSG))
                    print("Copy HKLOUT_OTHERSG")
                    shutil.copyfile(str(self.reindexJob.container.outputData.F_SIGF_OUT), str(self.container.outputData.HKLOUT_OTHERSG))

                    print("Annotate HKLOUT_OTHERSG")
                    self.container.outputData.HKLOUT_OTHERSG.annotation.set("Observations in spacegroup "+self.otherSG)
                    print("Annotate HKLOUT")
                    self.container.outputData.HKLOUT.annotation.set("Observations in spacegroup "+best_sg+" BEST")
                    print("Annotate FREEROUT_OTHERSG")
                    self.container.outputData.FREEROUT_OTHERSG.annotation.set("FreeR in spacegroup "+self.otherSG)
                    print("Annotate FREEROUT")
                    self.container.outputData.FREEROUT.annotation.set("FreeR in spacegroup "+best_sg+" BEST")

                    modelcraft = self.modelcraft
                    print("Set KV RFactor")
                    self.container.outputData.PERFORMANCE.RFactor.set(r0)
                    print("Set KV RFree")
                    self.container.outputData.PERFORMANCE.RFree.set(rf0)

                    shutil.copyfile(str(modelcraft.container.outputData.XYZOUT), str(self.container.outputData.XYZOUT))
                    shutil.copyfile(str(modelcraft.container.outputData.DIFFPHIOUT), str(self.container.outputData.DIFFPHIOUT))
                    shutil.copyfile(str(modelcraft.container.outputData.FPHIOUT), str(self.container.outputData.FPHIOUT))

                    self.container.outputData.XYZOUT.annotation.set('Atomic model after modelcraft autobuild in ' + best_sg)
                    self.container.outputData.FPHIOUT.annotation.set("Weighted map from modelcraft autobuild in " + best_sg)
                    self.container.outputData.DIFFPHIOUT.annotation.set("Weighted difference map from modelcraft autobuild in " + best_sg)
                else:

                    print("Copy FREER_OTHERSG")
                    shutil.copyfile(str(self.aimlessPlugin.container.outputData.FREEROUT), str(self.container.outputData.FREEROUT_OTHERSG))
                    print("Copy HKLOUT_OTHERSG")
                    shutil.copyfile(str(self.aimlessPlugin.container.outputData.HKLOUT[0]), str(self.container.outputData.HKLOUT_OTHERSG))
                    print("Copy FREEROUT")
                    shutil.copyfile(str(self.reindexJob.container.outputData.FREERFLAG_OUT), str(self.container.outputData.FREEROUT))
                    print("Copy HKLOUT")
                    shutil.copyfile(str(self.reindexJob.container.outputData.F_SIGF_OUT), str(self.container.outputData.HKLOUT))

                    print("Annotate HKLOUT")
                    self.container.outputData.HKLOUT.annotation.set("Observations in spacegroup "+self.otherSG+" BEST")
                    print("Annotate HKLOUT_OTHERSG")
                    self.container.outputData.HKLOUT_OTHERSG.annotation.set("Observations in spacegroup "+best_sg)
                    print("Annotate FREEROUT")
                    self.container.outputData.FREEROUT.annotation.set("FreeR in spacegroup "+self.otherSG+" BEST")
                    print("Annotate FREEROUT_OTHERSG")
                    self.container.outputData.FREEROUT_OTHERSG.annotation.set("FreeR in spacegroup "+best_sg)

                    modelcraft = self.modelcraft2
                    print("Set KV RFactor")
                    self.container.outputData.PERFORMANCE.RFactor.set(r1)
                    print("Set KV RFree")
                    self.container.outputData.PERFORMANCE.RFree.set(rf1)

                    shutil.copyfile(str(modelcraft.container.outputData.XYZOUT), str(self.container.outputData.XYZOUT))
                    shutil.copyfile(str(modelcraft.container.outputData.DIFFPHIOUT), str(self.container.outputData.DIFFPHIOUT))
                    shutil.copyfile(str(modelcraft.container.outputData.FPHIOUT), str(self.container.outputData.FPHIOUT))

                    self.container.outputData.XYZOUT.annotation.set('Atomic model after modelcraft autobuild in '+self.otherSG)
                    self.container.outputData.FPHIOUT.annotation.set("Weighted map from modelcraft autobuild in "+self.otherSG)
                    self.container.outputData.DIFFPHIOUT.annotation.set("Weighted difference map from modelcraft autobuild in "+self.otherSG)

                if self.container.controlParameters.LIGANDAS.__str__() != 'NONE':
                    print("Add ligand in enantiomorphic case")
                    self.coordinatesForCoot = modelcraft.container.outputData.XYZOUT
                    self.mapToUse = modelcraft.container.outputData.FPHIOUT
                    return self.cootAddLigand()
                else:
                    print("Report status")
                    self.reportStatus(rvfin2)
            except:
                self.reportStatus(CPluginScript.FAILED)

        else:
            self.molrep_job = self.makePluginObject('molrep_mr')
            self.molrep_job.container.inputData.XYZIN.set(self.XYZIN)
            if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
                self.molrep_job.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
            elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
                self.molrep_job.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
            else:
                self.molrep_job.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
            self.molrep_job.container.controlParameters.SG_OPTIONS.set('no')
            self.molrep_job.container.controlParameters.NMON.set(self.container.controlParameters.NMON)
            if self.container.inputData.XYZINORMRBUMP.__str__() == 'XYZINPUT':
                self.molrep_job.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
            rv = self.molrep_job.process()
            rvfin = self.fin_molrep(rv,self.molrep_job)
            rfactors = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/BuccaneerBuildRefineResult/FinalStatistics/r_factor')
            rfrees = self.xmlroot.xpath('/CCP4i2DRMRMBPipe/BuccaneerBuildRefineResult/FinalStatistics/r_free')
            try:
                print("##################################################")
                print("##################################################")
                print("##################################################")
                print("FRACTORS FROM MODELCRAFT")
                print("Try to read json...")
                path = self.modelcraft.workDirectory / "modelcraft" / "modelcraft.json"
                print("...", path)
                with path.open() as stream:
                    result = json.load(stream)
                print("read ...")
                print(result)
                self.container.outputData.PERFORMANCE.RFactor.set(result["final"]["r_work"])
                self.container.outputData.PERFORMANCE.RFree.set(result["final"]["r_free"])
                print("##################################################")
                print("##################################################")
                print("##################################################")
            except:
                print("PARSING JSON FAILS!"); sys.stdout.flush()

#Do ligand fit now .... (TODO - add this to anantiomorphic version above)
            if self.container.controlParameters.LIGANDAS.__str__() != 'NONE':
                modelcraft = self.modelcraft
                self.coordinatesForCoot = modelcraft.container.outputData.XYZOUT
                self.mapToUse = modelcraft.container.outputData.FPHIOUT
                return self.cootAddLigand()
            else:
                self.reportStatus(rvfin)

    def fin_molrep(self, status, molrep_job):
      print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fin_molrep 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"); sys.stdout.flush()
      #if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
      if status is not None and status == CPluginScript.FAILED:
        self.reportStatus(status)
        return
      pluginRoot = CCP4Utils.openFileToEtree(molrep_job.makeFileName('PROGRAMXML'))
      self.xmlroot.append(pluginRoot)
      self.flushXML()

      if self.useOtherSG:
          self.sheetbendPlugin2 = self.makePluginObject('sheetbend')
          self.sheetbendPlugin2.container.inputData.XYZIN.set(molrep_job.container.outputData.XYZOUT)
          self.sheetbendPlugin2.container.inputData.F_SIGF.set(self.reindexJob.container.outputData.F_SIGF_OUT)
          print("Setting FREER of sheetbend plugin")
          self.sheetbendPlugin2.container.inputData.FREERFLAG.set(self.reindexJob.container.outputData.FREERFLAG_OUT)
          rv = self.sheetbendPlugin2.process()
          return self.finSheetbend(rv,self.sheetbendPlugin2)
      else:
          self.sheetbendPlugin = self.makePluginObject('sheetbend')
          self.sheetbendPlugin.container.inputData.XYZIN.set(molrep_job.container.outputData.XYZOUT)
          if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
              self.sheetbendPlugin.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
              self.sheetbendPlugin.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
          elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
              self.sheetbendPlugin.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
              self.sheetbendPlugin.container.inputData.FREERFLAG.set(self.aimlessPlugin.freerflag.container.outputData.FREEROUT)
          else:
              self.sheetbendPlugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
              self.sheetbendPlugin.container.inputData.FREERFLAG.set(self.container.inputData.FREER_IN)
          rv = self.sheetbendPlugin.process()
          return self.finSheetbend(rv,self.sheetbendPlugin)

    def finSheetbend(self, status, sheetbend_job):
      if status is not None and status == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = CCP4Utils.openFileToEtree(sheetbend_job.makeFileName('PROGRAMXML'))
      self.xmlroot.append(tree)
      self.flushXML()

      if self.useOtherSG:
          self.refmac2 = self.makePluginObject('refmac')
          self.refmac2.container.inputData.XYZIN.set(self.sheetbendPlugin2.container.outputData.XYZOUT)
          self.refmac2.container.inputData.F_SIGF.set(self.reindexJob.container.outputData.F_SIGF_OUT)
          print("Setting FREER of refmac job")
          self.refmac2.container.inputData.FREERFLAG.set(self.reindexJob.container.outputData.FREERFLAG_OUT)
          #self.refmac2.container.controlParameters.HYDROGENS.set('NO')
          self.refmac2.container.controlParameters.HYDR_USE.set(False)
          self.refmac2.container.controlParameters.NCYCLES.set(str(self.container.inputData.REFMAC_NCYC))
          self.refmac2.container.controlParameters.PHOUT.set(False)
          rv = self.refmac2.process()
          return self.fin_refmac(rv,self.refmac2)
      else:
          self.refmac = self.makePluginObject('refmac')
          self.refmac.container.inputData.XYZIN.set(self.sheetbendPlugin.container.outputData.XYZOUT)
          if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
              self.refmac.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
              self.refmac.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
          elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
              self.refmac.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
              self.refmac.container.inputData.FREERFLAG.set(self.aimlessPlugin.freerflag.container.outputData.FREEROUT)
          else:
              self.refmac.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
              self.refmac.container.inputData.FREERFLAG.set(self.container.inputData.FREER_IN)
          #self.refmac.container.controlParameters.HYDROGENS.set('NO')
          self.refmac.container.controlParameters.HYDR_USE.set(False)
          self.refmac.container.controlParameters.NCYCLES.set(str(self.container.inputData.REFMAC_NCYC))
          self.refmac.container.controlParameters.PHOUT.set(False)
          rv = self.refmac.process()
          return self.fin_refmac(rv,self.refmac)

    def fin_refmac(self, status, refmac_job):
      print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fin_refmac 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"); sys.stdout.flush()
      #if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
      if status is not None and status == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = CCP4Utils.openFileToEtree(refmac_job.makeFileName('PROGRAMXML'))
      self.xmlroot.append(tree)
      self.flushXML()

      try:
        self.container.outputData.XYZOUT.annotation='Atomic model after restrained refinement'
        cycles = list()
        for e3 in tree.findall('./Overall_stats/stats_vs_cycle/new_cycle'):
          cycles.append((int(e3.findtext('cycle')), float(e3.findtext('r_factor')), float(e3.findtext('r_free'))))

        cycle_no, Rcrist, Rfree = sorted(cycles)[-1]
        self.container.outputData.PERFORMANCE.RFactor.set(Rcrist)
        self.container.outputData.PERFORMANCE.RFree.set(Rfree)

        if self.container.inputData.RUNACORN:
            if self.useOtherSG:
                self.acorn2 = self.makePluginObject('acorn')
                self.acorn2.container.inputData.XYZIN.set(self.refmac2.container.outputData.XYZOUT)
                self.acorn2.container.inputData.F_SIGF.set(self.reindexJob.container.outputData.F_SIGF_OUT)
                self.acorn2.container.inputData.FREERFLAG.set(self.reindexJob.container.outputData.FREERFLAG_OUT)
                rv = self.acorn2.process()
                return self.fin_acorn(rv,self.acorn2)
            else:
                self.acorn = self.makePluginObject('acorn')
                self.acorn.container.inputData.XYZIN.set(self.refmac.container.outputData.XYZOUT)
                if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
                    self.acorn.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
                    self.acorn.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
                elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
                    self.acorn.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
                    self.acorn.container.inputData.FREERFLAG.set(self.aimlessPlugin.freerflag.container.outputData.FREEROUT)
                else:
                    self.acorn.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
                    self.acorn.container.inputData.FREERFLAG.set(self.container.inputData.FREER_IN)
                rv = self.acorn.process()
                return self.fin_acorn(rv,self.acorn)
        else:
            return self.run_modelcraft()
      except:
        self.reportStatus(CPluginScript.FAILED)

    def fin_acorn(self, status,acorn_job):
      #print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fin_acorn 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"); sys.stdout.flush()
      #if status is not None and status.get('finishStatus') == CPluginScript.FAILED:
      if status is not None and status == CPluginScript.FAILED:
        self.reportStatus(status)
        return

      tree = CCP4Utils.openFileToEtree(acorn_job.makeFileName('PROGRAMXML'))
      self.xmlroot.append(tree)
      self.flushXML()

      return self.run_modelcraft()

    def run_modelcraft(self):
      try:
        if self.useOtherSG:
            self.modelcraft2 = self.makePluginObject('modelcraft')
            self.modelcraft2.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
            self.modelcraft2.container.inputData.XYZIN.set(self.refmac2.container.outputData.XYZOUT)
            self.modelcraft2.container.inputData.F_SIGF.set(self.reindexJob.container.outputData.F_SIGF_OUT)
            self.modelcraft2.container.inputData.FREERFLAG.set(self.reindexJob.container.outputData.FREERFLAG_OUT)
            modelcraft = self.modelcraft2
            if self.container.inputData.RUNACORN:
                modelcraft.container.inputData.PHASES.set(self.acorn2.container.outputData.PHSOUT)
        else:
            self.modelcraft = self.makePluginObject('modelcraft')
            self.modelcraft.container.inputData.ASUIN.set(self.container.inputData.ASUIN)
            self.modelcraft.container.inputData.XYZIN.set(self.refmac.container.outputData.XYZOUT)
            if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
                self.modelcraft.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.HKLOUT[0])
                self.modelcraft.container.inputData.FREERFLAG.set(self.aimlessPlugin.container.outputData.FREEROUT)
            elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
                self.modelcraft.container.inputData.F_SIGF.set(self.aimlessPlugin.container.outputData.OBSOUT)
                self.modelcraft.container.inputData.FREERFLAG.set(self.aimlessPlugin.freerflag.container.outputData.FREEROUT)
            else:
                self.modelcraft.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF_IN)
                self.modelcraft.container.inputData.FREERFLAG.set(self.container.inputData.FREER_IN)
            modelcraft = self.modelcraft
            if self.container.inputData.RUNACORN:
                modelcraft.container.inputData.PHASES.set(self.acorn.container.outputData.PHSOUT)

        modelcraft.container.controlParameters.SHEETBEND.set(self.container.inputData.RUNSHEETBEND) 
        modelcraft.container.controlParameters.CYCLES.set(self.container.inputData.MODELCRAFT_NCYC)

        self.oldTree = copy.deepcopy(self.xmlroot)
        rv = self.processModelCraft(modelcraft)
        print("rv from modelcraft",rv)

        self.xmlroot = copy.deepcopy(self.oldTree)
        print("Do something else!"); sys.stdout.flush()
        tree = etree.Element("ModelCraft")
        tree.text = str(modelcraft.workDirectory / "modelcraft")
        self.xmlroot.append(tree)
        print("Done something else!"); sys.stdout.flush()
        self.flushXML()

        model_build_text = "ModelCraft"

        shutil.copyfile(str(modelcraft.container.outputData.XYZOUT), str(self.container.outputData.XYZOUT))
        shutil.copyfile(str(modelcraft.container.outputData.DIFFPHIOUT), str(self.container.outputData.DIFFPHIOUT))
        shutil.copyfile(str(modelcraft.container.outputData.FPHIOUT), str(self.container.outputData.FPHIOUT))
        self.container.outputData.XYZOUT.annotation.set('Atomic model after '+model_build_text)
        self.container.outputData.FPHIOUT.annotation.set('Weighted map from '+model_build_text)
        self.container.outputData.DIFFPHIOUT.annotation.set('Weighted difference map from '+model_build_text)
        if self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'UNMERGED':
            shutil.copyfile(str(self.aimlessPlugin.container.outputData.FREEROUT), str(self.container.outputData.FREEROUT))
            shutil.copyfile(str(self.aimlessPlugin.container.outputData.HKLOUT[0]), str(self.container.outputData.HKLOUT))
        elif self.container.controlParameters.MERGED_OR_UNMERGED.__str__() == 'MERGED_F':
            shutil.copyfile(str(self.aimlessPlugin.freerflag.container.outputData.FREEROUT), str(self.container.outputData.FREEROUT))
            shutil.copyfile(str(self.aimlessPlugin.container.outputData.OBSOUT), str(self.container.outputData.HKLOUT))
        else:
            #This might be unnecessary
            shutil.copyfile(str(self.container.inputData.FREER_IN), str(self.container.outputData.FREEROUT))
            shutil.copyfile(str(self.container.inputData.F_SIGF_IN), str(self.container.outputData.HKLOUT))

        return rv

      except:
        raise
        self.reportStatus(CPluginScript.FAILED)

    @Slot(dict)
    def cootPlugin_finished(self, status):
        print("\n\n1", status)
        if status.get('finishStatus') == CPluginScript.FAILED:
            self.reportStatus(status)
        print("\n\n1","beyond")

        self.harvestFile(self.cootPlugin.container.outputData.XYZOUT[0], self.container.outputData.XYZOUT)
        self.container.outputData.XYZOUT.annotation.set('Atomic model after modelcraft and coot add ligand')
        self.reportStatus(status)

    def cootAddLigand(self):
        self.cootPlugin = self.makePluginObject('coot_script_lines')
        xyzinList = self.cootPlugin.container.inputData.XYZIN
        xyzinList.append(xyzinList.makeItem())
        xyzinList[-1].set(self.coordinatesForCoot)
        fphiinList = self.cootPlugin.container.inputData.FPHIIN
        fphiinList.append(fphiinList.makeItem())
        fphiinList[-1].set(self.mapToUse)
        self.cootPlugin.container.inputData.DICT.set(self.dictToUse)
        self.cootPlugin.container.controlParameters.SCRIPT.set('''#Script to fit ligand into density
monomerMolNo = get_monomer('DRG')
add_ligand_clear_ligands()
set_ligand_search_protein_molecule(MolHandle_1)
set_ligand_search_map_molecule(MapHandle_1)
add_ligand_search_wiggly_ligand_molecule(monomerMolNo)
#Execute search
nToCopy = 0
ligandsFound=execute_ligand_search()
if ligandsFound is not False:
    nToCopy = len(ligandsFound)
#Check on ncs
equivs = ncs_chain_ids(0)
if equivs is not False and len(equivs)>0:
    nToCopy = min(len(ligandsFound),len(equivs[0]))
if nToCopy > 0:
    ligandsToCopy = ligandsFound[0:nToCopy]
    merge_molecules(ligandsToCopy,0)

write_pdb_file(MolHandle_1,os.path.join(dropDir,"output.pdb"))''')
        self.connectSignal(self.cootPlugin,'finished',self.cootPlugin_finished)
        self.cootPlugin.process()

    def processModelCraft(self, plugin, **kw):
        #I am reimplementing this because I want to be able to reproduce the top part of ModelCraft pipeline so that I can get its XML.
        ''' Check input data is set, create program command script (by calling makeCommandAndScript
        which should be implemented in sub-class and call startProcess '''
        #print 'CPluginScript.process',plugin.objectName()
        #plugin.loadProjectDefaults()
        try:
            unsetData = plugin.checkInputData()
        except:
            plugin.appendErrorReport(41)
            return plugin.reportStatus(CPluginScript.FAILED)
        #print 'CPluginScript.process unsetData',unsetData
        if len(unsetData) > 0:
            return plugin.reportStatus(CPluginScript.FAILED)
        try:
            rv = plugin.checkOutputData(plugin.container)
            #print 'CPluginScript.process unsetOutputData',e
        except Exception as e:
            plugin.appendErrorReport(42, exc_info=sys.exc_info())
        else:
            if len(rv) > 0:
                plugin.extendErrorReport(rv)
        try:
            status = plugin.processInputFiles()
        except CException as e:
            return plugin.reportStatus(CPluginScript.FAILED)
        except Exception as e:
            plugin.appendErrorReport(43, exc_info=sys.exc_info())
            plugin.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        else:
            #print 'CPluginScript.process processInputFiles',status
            if status == CPluginScript.FAILED:
                return plugin.reportStatus(CPluginScript.FAILED)

        if plugin.editComFile:
            plugin.displayEditor()
            return
        try:
            rv = self.startModelCraftProcess(plugin)
        except:
            plugin.appendErrorReport(48, exc_info=sys.exc_info())
            return plugin.reportStatus(CPluginScript.FAILED)
        else:
            if rv == CPluginScript.FAILED:
                return plugin.reportStatus(rv)
        if not plugin._ifAsync:
            return plugin.postProcess(processId=plugin._runningProcessId)
        else:
            return CPluginScript.SUCCEEDED

    def startModelCraftProcess(self, plugin):
        print("##################################################")
        print("startModelCraftProcess",plugin)
        print("##################################################")

        self.oldTree = copy.deepcopy(self.xmlroot)
        self.xmlroot = copy.deepcopy(self.oldTree)
        tree = etree.Element("ModelCraft")
        tree.text = str(plugin.workDirectory / "modelcraft")
        self.xmlroot.append(tree)
        self.flushXML()
        rv = plugin.process()
        if rv == CPluginScript.FAILED:
            return plugin.reportStatus(rv)
        else:
            plugin.reportStatus(CPluginScript.SUCCEEDED)
            return CPluginScript.SUCCEEDED

        plugin.initialiseProperties()
        plugin.initialiseXML()
        if plugin.restarted: plugin.handleRestart()

        if plugin.shouldCalculateModelPhases():
            plugin.refmacABCD()

        for plugin.cycle in range(plugin.cycle, plugin.ncycles):
            plugin.processedPlugins[plugin.cycle] = []
            plugin.cleanUpIntermediateFiles()
            plugin.addXMLCycle()

            print('Starting cycle', plugin.cycle + 1, 'of', plugin.ncycles)
            try:
                plugin.buildRefineCycle()
            except:
                print('buildRefineCycle failed, reporting error'); sys.stdout.flush()
                return self.reportStatus(CPluginScript.FAILED)
            print('Finished cycle', plugin.cycle + 1, 'of', plugin.ncycles)

            plugin.checkForImprovement()
            if plugin.auto_stop and plugin.cycles_without_improvement == 4: break

            if plugin.testForInterrupt() and plugin.cycle + 1 < plugin.ncycles:
                plugin.alignment()
                plugin.container.interruptStatus.LASTCYCLE.set(plugin.cycle)
                plugin.reportStatus(CPluginScript.INTERRUPTED)

            self.xmlroot = copy.deepcopy(self.oldTree)
            tree = CCP4Utils.openFileToEtree(plugin.makeFileName('PROGRAMXML'))
            self.xmlroot.append(tree)
            self.flushXML()

        print("plugin.writeFinalXML()")
        plugin.writeFinalXML()
        print("plugin.alignment()")
        plugin.alignment()
        print("cleanup")
        print("report status")
        plugin.reportStatus(CPluginScript.SUCCEEDED)

    def flushXML(self):
        with open(self.makeFileName('PROGRAMXML'),'w') as programXML:
            CCP4Utils.writeXML(programXML,etree.tostring(self.xmlroot,pretty_print=True))
