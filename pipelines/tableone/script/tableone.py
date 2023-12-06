from __future__ import print_function

import os
import re
import sys
import math
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from core import CCP4ErrorHandling
from core import CCP4Utils
import clipper
from wrappers.validate_protein.script import validate_protein
from iotbx import mtz
import ccp4mg
import mmdb2
from mmtbx.command_line import molprobity

class tableone(CPluginScript):

    TASKMODULE = 'test'                 # Gui menu location
    TASKTITLE = 'Generate Table One'          # Short title for Gui
    TASKNAME = 'tableone'                     # Task name - same as class name
    TASKCOMMAND = 'xia2.merging_statistics'   # The command to run the executable - there is more than one here.
    TASKVERSION = 1.0                         # plugin version
    COMTEMPLATE = None
    COMTEMPLATEFILE = None
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = True
    MAINTAINER = 'Kyle.Stevenson@stfc.ac.uk'


    def __init__(self, *args, **kwargs):
        self.xml_root = etree.Element("tableone")
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        CPluginScript.process(self)
        valpro_pth1 = os.path.join(CCP4Utils.getCCP4Dir(), 'share', 'ccp4i2', 'wrappers', 'validate_protein', 'script')
        sys.path.append(valpro_pth1)
        # Be careful with this. validate_protein may well change, best to co-ordinate this.
        vprotein = validate_protein.validate_protein()
        l1, x1 = vprotein.b_averages(str(self.container.inputData.XYZIN))
        l2, x2 = vprotein.ramachandran_maps(str(self.container.inputData.XYZIN))
        self.run_molprobity()
        self.xml_root.append(x1)
        self.xml_root.append(x2)
        f = clipper.MMDBfile()
        f.read_file(str(self.container.inputData.XYZIN))
        mmol = clipper.MiniMol()
        f.import_minimol(mmol)
        spg = mmol.spacegroup()
        cel = mmol.cell()
        spgt = spg.symbol_hm()
        clenw = ["a", "b", "c"]
        cangw = ["alpha", "beta", "gamma"]
        clen = []
        cang = []
        clen.append(str(cel.a()))
        clen.append(str(cel.a()))
        clen.append(str(cel.c()))
        cang.append(str(math.degrees(cel.alpha())))
        cang.append(str(math.degrees(cel.beta())))
        cang.append(str(math.degrees(cel.gamma())))
        xmlpdb = etree.SubElement(self.xml_root, "PdbInfo")
        etree.SubElement(xmlpdb, "SpaceGroup").text = spgt
        for t1, v1, t2, v2 in zip(clenw, clen, cangw, cang):
            etree.SubElement(xmlpdb, t1).text = v1
            etree.SubElement(xmlpdb, t2).text = v2
        return CPluginScript.SUCCEEDED

    def PreparePdbForMolprobity(self, pdbin, workdir):
        # Again had to split off & replicate.
        file_name = os.path.basename(pdbin)
        wrkdirnm = os.path.join(workdir, file_name)
        sanitized_pdbin, pdbin_ext = os.path.splitext(wrkdirnm)
        sanitized_pdbin += "_as_PDB.pdb"
        file_root, file_ext = os.path.splitext(wrkdirnm)
        #Use mmdb to do some sanitization
        mmdb2.InitMatType()
        m = mmdb2.Manager()
        #This line ought to just do it, but seems not to...
        m.SetFlag(mmdb2.MMDBF_IgnoreSegID)
        m.ReadCoorFile(pdbin)
        #Remove any SEGIDs (which can confuse molprobity if heterogenous)
        sel = m.NewSelection()
        m.SelectAtoms(sel, 0, "*", mmdb2.ANY_RES, "*", mmdb2.ANY_RES, "*", "*", "*", "*", "*", mmdb2.SKEY_OR)
        selindexp = mmdb2.intp()
        selAtoms = mmdb2.GetAtomSelIndex(m, sel, selindexp)
        nSelAtoms = selindexp.value()
        # ... but this certainly does.
        for i in range(nSelAtoms):
            at = mmdb2.getPCAtom(selAtoms, i)
            at.segID = "    "
        m.FinishStructEdit()
        m.WritePDBASCII(sanitized_pdbin)
        return sanitized_pdbin, file_root

    def run_molprobity(self):
        # Unfortunate, but need to replicate part of the validate_protein code for time being.
        sanitized_pdbin, file_root = self.PreparePdbForMolprobity(str(self.container.inputData.XYZIN), self.workDirectory)
        self.mlog_filename = os.path.join(self.workDirectory, "molprobity.log")
        with open(self.mlog_filename, "w") as molprobity_logfile :
            molprobity.run(["input.pdb.file_name={}".format(sanitized_pdbin),
                            "output.prefix={}".format(file_root)], out=molprobity_logfile )
        return

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.parseLogfile()
        with open(self.makeFileName('PROGRAMXML'), 'w') as xml_file:
            CCP4Utils.writeXML(xml_file,etree.tostring(self.xml_root, pretty_print=True))
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        self.appendCommandLine(str(self.container.inputData.UNMERGED.fullPath))
        self.appendCommandLine("labels=I,SIGI")
        # Get the resolution bounds from the merged mtz file
        mobj = mtz.object(str(self.container.inputData.F_SIGF.fullPath))
        resol = mobj.max_min_resolution()
        maxr = resol[0]
        minr = resol[1]
        self.appendCommandLine("high_resolution=%s"%str(minr))
        self.appendCommandLine("low_resolution=%s"%str(maxr))
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

    def parseLogfile(self):
        print("Parsing the logfiles and creating xml")
        logf = open(self.makeFileName('LOG'), 'r')
        lrin = logf.read()
        logf.close()
        val = {}
        val.update({"Observ" : int(re.search(r'Observations:\s*(.*)', lrin).group(1))})
        val.update({"Unique" : int(re.search(r'Unique\sreflections:\s*(.*)\s*(.*)', lrin).group(1))})
        val.update({"Multi" : float(re.search(r'Redundancy:\s(\d.\d)', lrin).group(1))})
        val.update({"Compl" : float(re.search(r'Completeness:\s*(\d*.\d*)', lrin).group(1))})
        val.update({"SigN" : float(re.search(r'Mean\sI.sigma.I.:\s*(\d*.\d*)', lrin).group(1))})
        val.update({"Rmer" : float(re.search(r'R-merge:\s*(\d*.\d*)', lrin).group(1))})
        val.update({"Rmeas" : float(re.search(r'R-meas:\s*(\d*.\d*)', lrin).group(1))})
        val.update({"Rpim" : float(re.search(r'R-pim:\s*(\d*.\d*)', lrin).group(1))})
        # Get the table info for correlation info. The last tab line has what you want.
        lrinl = lrin.splitlines()
        ftab = False
        tabnums = []
        standtablen = 13  # This is the length of the logfile table. If it changes this will stop working.
        for lin in lrinl:
            if ftab == True:
                tabobj = re.findall(r"[-+]?\d+[\.]?\d*", lin)
                if len(tabobj) == standtablen:
                    tabnums.append(tabobj)
                else:
                    break
            if ftab == False:
                ftab = bool(re.search(r"cc1/2\s+cc_ano$",lin))
        cchalf = tabnums[-1][-2]
        ccanom = tabnums[-1][-1]
        val.update({"cchalf" : cchalf})
        val.update({"ccanom" : ccanom})
        # Now extract the molprobity gubbins
        logmolf = open(self.mlog_filename, 'r')
        lrinmol = logmolf.read()
        logmolf.close()
        val.update({"Rotamers" : float(re.search(r"Rotamer outliers\s+=\s+(\d+.\d+)", lrinmol).group(1))})
        val.update({"Clashscore" : float(re.search(r"Clashscore\s+=\s+(\d+.\d+)", lrinmol).group(1))})
        xmlsec = etree.SubElement(self.xml_root, "DataInfo")
        for x,y in val.items():
            etree.SubElement(xmlsec, x).text = str(y)
        return
