from report.CCP4ReportParser import *
from core import CCP4Modules
import sys, os
from xml.etree import ElementTree as ET


class metalCoord_report(Report):
    TASKNAME='metalCoord'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        # projectid = self.jobInfo.get("projectid", None)
        # jobNumber = self.jobInfo.get("jobnumber", None)
        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
        self.jobLog = os.path.join(jobDirectory, "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:  # elif jobStatus in ["Finished"]:
            self.defaultReport(parent=self)


    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        if os.path.isfile(self.jobLog):
            jobLogFold = parent.addFold(label="MetalCoord log", initiallyOpen=True)
            jobLogFold.addPre("MetalCoord is running...")


    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title

        def makeAddress(entry, node=""):
            if node:  # for ligand, ligand1, ligand2
                nodePrefix = node + "/"
            else:
                nodePrefix = ""
            atomAddress = entry.findall(nodePrefix + "chain")[0].text + "/" + \
                entry.findall(nodePrefix + "residue")[0].text + " " + \
                entry.findall(nodePrefix + "sequence")[0].text
            if entry.findall(nodePrefix + "icode")[0].text == ".":
                icode = ""
            if icode: atomAddress += "." + icode
            if entry.findall(nodePrefix + "name"):
                name_or_metal = "name"   # for ligand, ligand1, ligand2
            else:  # if entry.findall(nodePrefix + "metal"):
                name_or_metal = "metal"  # for metal
            atomAddress += "/" + entry.findall(nodePrefix + name_or_metal)[0].text
            if entry.findall(nodePrefix + "altloc")[0].text:
                atomAddress += "." + entry.findall(nodePrefix + "altloc")[0].text
            return atomAddress
            
        for i, site in enumerate(self.xmlnode.findall(".//site")):
            metalAddress = makeAddress(site, node="")
            for j, symmClass in enumerate(site.findall(".//ligands")):
                table = parent.addTable(xmlnode=symmClass)
                table.addData(title="Symmetry class", select='class')
                table.addData(title="Procrustes", select='procrustes')
                table.addData(title="Coordination", select='coordination')
                table.addData(title="Count", select='count')
                table.addData(title="Description", select='description')
                indentDiv = parent.addDiv(style="margin-left:3em;")
                for k, entry in enumerate(symmClass.findall(".//base")):
                    n_options = len(entry.findall("std"))
                    #if n_options >= 2:
                    #    for l in range(n_options - 1):
                    #        entry.append(copy.deepcopy(entry.findall("ligand")[0]))
                    for l in range(n_options):
                        ligandAddress = makeAddress(entry, node="ligand")
                        ligandAddressElement = ET.SubElement(entry, "ligandAtomAddress")
                        ligandAddressElement.text = ligandAddress
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                if symmClass.findall(".//base"):
                    table = indentDiv.addTable(xmlnode=symmClass)
                    table.addData(title="Metal site", select='base/metalAtomAddress')
                    table.addData(title="Atom", select='base/ligandAtomAddress')
                    table.addData(title="Distance (&Aring;)", select='base/distance')
                    table.addData(title="St. dev. (&Aring;)", select='base/std')

                for k, entry in enumerate(symmClass.findall(".//pdb")):
                    n_options = len(entry.findall("std"))
                    for l in range(n_options):
                        ligandAddress = makeAddress(entry, node="ligand")
                        ligandAddressElement = ET.SubElement(entry, "ligandAtomAddress")
                        ligandAddressElement.text = ligandAddress
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                if symmClass.findall(".//pdb"):
                    table = indentDiv.addTable(xmlnode=symmClass)
                    table.addData(title="Metal site", select='pdb/metalAtomAddress')
                    table.addData(title="Atom", select='pdb/ligandAtomAddress')
                    table.addData(title="Distance (&Aring;)", select='pdb/distance')
                    table.addData(title="St. dev. (&Aring;)", select='pdb/std')

                for k, entry in enumerate(symmClass.findall(".//angles")):
                    n_options = len(entry.findall("std"))
                    for l in range(n_options):
                        ligand1Address = makeAddress(entry, node="ligand1")
                        ligand1AddressElement = ET.SubElement(entry, "ligand1AtomAddress")
                        ligand1AddressElement.text = ligand1Address
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                        ligand2Address = makeAddress(entry, node="ligand2")
                        ligand2AddressElement = ET.SubElement(entry, "ligand2AtomAddress")
                        ligand2AddressElement.text = ligand2Address
                if symmClass.findall(".//angles"):
                    table = indentDiv.addTable(xmlnode=symmClass)
                    table.addData(title="Atom", select='angles/ligand1AtomAddress')
                    table.addData(title="Metal site", select='angles/metalAtomAddress')
                    table.addData(title="Atom", select='angles/ligand2AtomAddress')
                    table.addData(title="Angle (&deg;)", select='angles/angle')
                    table.addData(title="St. dev. (&deg;)", select='angles/std')

        self.addDiv(style="clear:both;")
