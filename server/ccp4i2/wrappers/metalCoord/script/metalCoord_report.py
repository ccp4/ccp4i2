from ccp4i2.report.CCP4ReportParser import *
from ccp4i2.core import CCP4Modules
import os
from xml.etree import ElementTree as ET


class metalCoord_report(Report):
    TASKNAME='metalCoord'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
        self.jobLog = os.path.join(jobDirectory, "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:
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
            else:
                name_or_metal = "metal"  # for metal
            atomAddress += "/" + entry.findall(nodePrefix + name_or_metal)[0].text
            if entry.findall(nodePrefix + "altloc")[0].text:
                atomAddress += "." + entry.findall(nodePrefix + "altloc")[0].text
            return atomAddress

        def makeSymmetry(entry, node="ligand"):
            atomSymmetry = entry.findall(node + "/symmetry")[0].text
            if int(atomSymmetry) == 0:
                return "-"
            else:
                return str(atomSymmetry)

        for site in self.xmlnode.findall(".//site"):
            metalAddress = makeAddress(site, node="")
            siteFold = parent.addFold(label="Metal site " + metalAddress, initiallyOpen=True)
            n_classes = len(site.findall(".//ligands"))
            noteDiv = siteFold.addDiv(style='font-size:110%;')
            if n_classes == 1:
                noteDiv.append("Only one coordination class is reported for this metal site.")
            elif n_classes > 1:
                noteDiv.append(f"{n_classes} possible coordination classes are reported. The class with the lowest procrustes distance is the most favourable, has been used for restraints generation and is listed as the first.")
            for j, symmClass in enumerate(site.findall(".//ligands")):
                indentDiv = siteFold.addDiv(style="margin-left:3em;")
                if j == 0: initiallyOpen=True
                else: initiallyOpen=False
                label = ""
                if len(symmClass.findall(".//class")) > 0:
                    if symmClass.findall(".//class")[0].text:
                        label = symmClass.findall(".//class")[0].text + " coordination class"
                classFold = indentDiv.addFold(
                    label=label,
                    initiallyOpen=initiallyOpen)
                table = classFold.addTable(xmlnode=symmClass)
                table.addData(title="Procrustes distance", select='procrustes')
                table.addData(title="Coordination number", select='coordination')
                table.addData(title="No. reference structures", select='count')
                table.addData(title="Note", select='description')

                headerDiv = classFold.addDiv(style='font-size:110%;font-weight:bold;')
                headerDiv.append("Ideal distances")
                for entry in symmClass.findall(".//base"):
                    n_options = len(entry.findall("std"))
                    for l in range(n_options):
                        ligandAddress = makeAddress(entry, node="ligand")
                        ligandAddressElement = ET.SubElement(entry, "ligandAtomAddress")
                        ligandAddressElement.text = ligandAddress
                        ligandSymmetry = makeSymmetry(entry, node="ligand")
                        ligandSymmetryElement = ET.SubElement(entry, "ligandAtomSymmetry")
                        ligandSymmetryElement.text = ligandSymmetry
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                if symmClass.findall(".//base"):
                    table = classFold.addTable(xmlnode=symmClass)
                    table.addData(title="Metal site", select='base/metalAtomAddress')
                    table.addData(title="Atom", select='base/ligandAtomAddress')
                    table.addData(title="Symmetry?", select='base/ligandAtomSymmetry')
                    table.addData(title="Distance (&Aring;)", select='base/distance')
                    table.addData(title="St. dev. (&Aring;)", select='base/std')

                for entry in symmClass.findall(".//pdb"):
                    n_options = len(entry.findall("std"))
                    for l in range(n_options):
                        ligandAddress = makeAddress(entry, node="ligand")
                        ligandAddressElement = ET.SubElement(entry, "ligandAtomAddress")
                        ligandAddressElement.text = ligandAddress
                        ligandSymmetry = makeSymmetry(entry, node="ligand")
                        ligandSymmetryElement = ET.SubElement(entry, "ligandAtomSymmetry")
                        ligandSymmetryElement.text = ligandSymmetry
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                if symmClass.findall(".//pdb"):
                    table = classFold.addTable(xmlnode=symmClass)
                    table.addData(title="Metal site", select='pdb/metalAtomAddress')
                    table.addData(title="Atom", select='pdb/ligandAtomAddress')
                    table.addData(title="Symmetry?", select='pdb/ligandAtomSymmetry')
                    table.addData(title="Distance (&Aring;)", select='pdb/distance')
                    table.addData(title="St. dev. (&Aring;)", select='pdb/std')

                headerDiv = classFold.addDiv(style='font-size:110%;font-weight:bold;')
                headerDiv.append("Ideal angles")
                for entry in symmClass.findall(".//angles"):
                    n_options = len(entry.findall("std"))
                    for l in range(n_options):
                        ligand1Address = makeAddress(entry, node="ligand1")
                        ligand1AddressElement = ET.SubElement(entry, "ligand1AtomAddress")
                        ligand1AddressElement.text = ligand1Address
                        ligand1Symmetry = makeSymmetry(entry, node="ligand1")
                        ligand1SymmetryElement = ET.SubElement(entry, "ligand1AtomSymmetry")
                        ligand1SymmetryElement.text = ligand1Symmetry
                        metalAddressElement = ET.SubElement(entry, "metalAtomAddress")
                        metalAddressElement.text = metalAddress
                        ligand2Address = makeAddress(entry, node="ligand2")
                        ligand2AddressElement = ET.SubElement(entry, "ligand2AtomAddress")
                        ligand2AddressElement.text = ligand2Address
                        ligand2Symmetry = makeSymmetry(entry, node="ligand2")
                        ligand2SymmetryElement = ET.SubElement(entry, "ligand2AtomSymmetry")
                        ligand2SymmetryElement.text = ligand2Symmetry
                if symmClass.findall(".//angles"):
                    table = classFold.addTable(xmlnode=symmClass)
                    table.addData(title="Atom", select='angles/ligand1AtomAddress')
                    table.addData(title="Symmetry?", select='angles/ligand1AtomSymmetry')
                    table.addData(title="Metal site", select='angles/metalAtomAddress')
                    table.addData(title="Atom", select='angles/ligand2AtomAddress')
                    table.addData(title="Symmetry?", select='angles/ligand2AtomSymmetry')
                    table.addData(title="Angle (&deg;)", select='angles/angle')
                    table.addData(title="St. dev. (&deg;)", select='angles/std')

        self.addDiv(style="clear:both;")
