from __future__ import print_function

import sys
import os
import zipfile
import shutil
import json

from lxml import etree

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4XtalData import CMapCoeffsDataFile, CObsDataFile, CPhsDataFile
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core import CCP4Modules

from . import test_api

class pdb_redo_api(CPluginScript):

    TASKTITLE='PDB-REDO Web services'
    TASKNAME = 'pdb_redo_api'
    TASKMODULE= 'refinement'
    TASKCOMMAND = ''
    TASKVERSION= 0.0
    COMLINETEMPLATE = None
    COMTEMPLATE = None
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS=False
    PERFORMANCECLASS = 'CRefinementPerformance'

    def startProcess(self,command,**kw):
        print("pdb_redo_api.startProcess")
        inp = self.container.inputData
        xyzin = str( inp.XYZIN.fullPath )
        print(self.hklin, xyzin)

        token_id = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID)
        token_secret = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET)

        sequence=None
        restraints=None

        if inp.SEQIN.isSet():
            sequence = str(inp.SEQIN.fullPath)

        if inp.DICT.isSet():
            restraints = str(inp.DICT.fullPath)

        params = {
            'paired': int(self.container.controlParameters.PAIRED),
            'noloops': int(self.container.controlParameters.NOLOOPS),
            'nopepflip': int(self.container.controlParameters.NOPEPFLIP),
            'noscbuild': int(self.container.controlParameters.NOSCBUILD),
            'nocentrifuge': int(self.container.controlParameters.NOCENTRIFUGE),
            'nosugarbuild': int(self.container.controlParameters.NOSUGARBUILD),
            'norebuild': int(self.container.controlParameters.NOREBUILD),
            'newmodel': int(self.container.controlParameters.NEWMODEL),
            'isotropic': int(self.container.controlParameters.ISOTROPIC),
            'anisotropic': int(self.container.controlParameters.ANISOTROPIC),
            'notls': int(self.container.controlParameters.NOTLS),
            'tighter': int(self.container.controlParameters.TIGHTER),
            'looser': int(self.container.controlParameters.LOOSER),
        }

        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))

        print("Submit PDB-REDO job")
        self.pdb_redo_job_id = test_api.submit(xyzin,self.hklin,token_id,token_secret,sequence=sequence,restraints=restraints,params=params)
        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            pdbRedoId = etree.SubElement(xmlStructure,'PDB_REDO_JOB_ID')
            pdbRedoId.text = str(self.pdb_redo_job_id)
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))

        print("Wait for PDB-REDO job")
        try:
            test_api.monitor(self.pdb_redo_job_id,token_id,token_secret)
        except:
            try:
                print("Fetching/extracting from zip after job failure"); sys.stdout.flush()
                output_zip = os.path.join(self.getWorkDirectory(),"pdb_redo_results.zip")
                test_api.do_fetch(self.pdb_redo_job_id,token_id,token_secret,output_zip)
                print("Got",output_zip); sys.stdout.flush()
                with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
                    xmlStructure = etree.Element("pdb_redo_api")
                    pdbRedoId = etree.SubElement(xmlStructure,'PDB_REDO_JOB_ID')
                    pdbRedoId.text = str(self.pdb_redo_job_id)
                    with zipfile.ZipFile(output_zip) as myzip:
                        for f in myzip.namelist():
                            myzip.extract(f,self.getWorkDirectory())
                            pdbRedoDir = etree.SubElement(xmlStructure,'PDB_REDO_RESULTS_DIR')
                            redoDir = os.path.dirname(f)
                            pdbRedoDir.text = str(redoDir)
                    CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure))
            except:
                print("Failed to get pdb_redo_results.zip after job failure")
            
            return CPluginScript.FAILED
        
        print("PDB-REDO job has finished, leave rest to processOutputFiles")

        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        miniMtzs = [
            ["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN],
            ["FREERFLAG", None],
        ]
        self.hklin, self.columns, error = self.makeHklin0(miniMtzs)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        token_id = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_ID)
        token_secret = str(CCP4Modules.PREFERENCES().PDB_REDO_TOKEN_SECRET)

        print("Extracting from zip"); sys.stdout.flush()
        output_zip = os.path.join(self.getWorkDirectory(),"pdb_redo_results.zip")
        test_api.do_fetch(self.pdb_redo_job_id,token_id,token_secret,output_zip)
        print("Got",output_zip); sys.stdout.flush()
        
        redoDir = None
        finalRefmacLog = None
        pdbRedoLog = None

        with zipfile.ZipFile(output_zip) as myzip:
            outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']
            infolist = myzip.infolist()
            for finfo in infolist:
                if finfo.filename.endswith("_besttls.mtz"):
                    outputFiles = ['FPHIOUT_BESTTLS','DIFFPHIOUT_BESTTLS']
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    hkloutFile= os.path.join(self.getWorkDirectory(), finfo.filename)
                    print("Set hkloutFile",hkloutFile); sys.stdout.flush()
                    print("Splitting..."); sys.stdout.flush()
                    self.splitHklout(outputFiles,outputColumns,infile=hkloutFile)
                    print("Split",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith("_final.mtz"):
                    outputFiles = ['FPHIOUT_FINAL','DIFFPHIOUT_FINAL']
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    hkloutFile=os.path.join(self.getWorkDirectory(), finfo.filename)
                    print("Set hkloutFile",hkloutFile); sys.stdout.flush()
                    print("Splitting..."); sys.stdout.flush()
                    self.splitHklout(outputFiles,outputColumns,infile=hkloutFile)
                    print("Split",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith("_besttls.pdb"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    outputPDB      = os.path.normpath(os.path.join(self.getWorkDirectory(),finfo.filename))
                    outputFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),os.path.basename(finfo.filename)))
                    shutil.copyfile(outputPDB, outputFilePath)
                    self.container.outputData.XYZOUT_BESTTLS.setFullPath(outputFilePath)
                if finfo.filename.endswith("_final.pdb"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    outputPDB      = os.path.normpath(os.path.join(self.getWorkDirectory(),finfo.filename))
                    outputFilePath = os.path.normpath(os.path.join(self.getWorkDirectory(),os.path.basename(finfo.filename)))
                    shutil.copyfile(outputPDB, outputFilePath)
                    self.container.outputData.XYZOUT_FINAL.setFullPath(outputFilePath)
                if finfo.filename.endswith("_final.log"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    finalRefmacLog = finfo.filename
                if finfo.filename.endswith(".log") and not finfo.filename.endswith("_final.log") and os.path.basename(finfo.filename) != "process.log":
                    print("Extracting log",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted log",finfo.filename); sys.stdout.flush()
                    pdbRedoLog = finfo.filename
                if finfo.filename.endswith(".html"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                if finfo.filename.endswith(".json"):
                    print("Extracting",finfo.filename); sys.stdout.flush()
                    myzip.extract(finfo,self.getWorkDirectory())
                    print("Extracted",finfo.filename); sys.stdout.flush()
                    if os.path.basename(finfo.filename) == "data.json":
                        print("Reading",os.path.join(self.getWorkDirectory(),finfo.filename))
                        with open(os.path.join(self.getWorkDirectory(),finfo.filename)) as f:
                            print("Loading as JSON")
                            redoDir = os.path.dirname(finfo.filename)
                            j = json.load(f)
                            if "properties" in j and "RFIN" in j["properties"]:
                                print("Set performance R")
                                self.container.outputData.PERFORMANCE.RFactor.set(j["properties"]["RFIN"])
                            if "properties" in j and "RFFIN" in j["properties"]:
                                print("Set performance RFree")
                                self.container.outputData.PERFORMANCE.RFree.set(j["properties"]["RFFIN"])

        with open(self.makeFileName("PROGRAMXML"),"w") as programXMLFile:
            xmlStructure = etree.Element("pdb_redo_api")
            if "properties" in j:
                for k,v in j["properties"].items():
                    ele = etree.SubElement(xmlStructure,k)
                    ele.text = str(v)
            pdbRedoDir = etree.SubElement(xmlStructure,'PDB_REDO_RESULTS_DIR')
            pdbRedoDir.text = str(redoDir)
            pdbRedoId = etree.SubElement(xmlStructure,'PDB_REDO_JOB_ID')
            pdbRedoId.text = str(self.pdb_redo_job_id)
            if finalRefmacLog:
                finalRefmacLogEle = etree.SubElement(xmlStructure,'PDB_REDO_FINAL_REFMAC_LOG_FILE')
                finalRefmacLogEle.text = str(finalRefmacLog)
            if pdbRedoLog:
                pdbRedoLogEle = etree.SubElement(xmlStructure,'PDB_REDO_LOG_FILE')
                pdbRedoLogEle.text = str(pdbRedoLog)
            CCP4Utils.writeXML(programXMLFile,etree.tostring(xmlStructure,encoding='utf-8', pretty_print=True))

        return CPluginScript.SUCCEEDED
