from __future__ import print_function

from report.CCP4ReportParser import *
from core import CCP4Utils
import sys
from xml.etree import ElementTree as ET

class qtpisa_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'qtpisa'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        results = self.addResults()
        self.addText(text='The QtPISA session is finished.')
        pdbsWrittenPath = './/number_output_files'
        pdbsWrittenStringS = xmlnode.findall(pdbsWrittenPath)
        if len(pdbsWrittenStringS) > 0:
          pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
          self.append('<br/>')
          self.addText(text=pdbsWrittenString+' PDB files were written to the CCP4i2 database during the session.')
        else:
          self.append('<br/>')
          self.addText(text='No PDB files were written to the CCP4i2 database during the session.')

        print("SELF.XMLNODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",self.xmlnode)

        t = ET.tostring(self.xmlnode)

        ser_no_interface = []
        sym_id_interface = []
        symop_interface = []
        area_interface = []
        delta_g_interface = []
        nhbonds_interface = []
        nsbridges_interface = []
        ndsbonds_interface = []

        ser_no_assembly = []
        delta_g_assembly = []
        composition_assembly = []
        olig_state_assembly = []
        asa_assembly = []
        bsa_assembly = []

        for pisa_file in self.xmlnode.findall(".//qtpisa_report/pisa_file"):
          if pisa_file.findall("type")[0].text == "assembly":
            ser_no_assembly.append(pisa_file.findall("ser_no")[0].text)
            delta_g_assembly.append(pisa_file.findall("delta_g")[0].text)
            composition_assembly.append(pisa_file.findall("composition")[0].text)
            olig_state_assembly.append(pisa_file.findall("oligomeric_state")[0].text)
            asa_assembly.append(pisa_file.findall("asa")[0].text)
            bsa_assembly.append(pisa_file.findall("bsa")[0].text)
          elif pisa_file.findall("type")[0].text == "interface":
            ser_no_interface.append(pisa_file.findall("ser_no")[0].text)
            sym_id_interface.append(pisa_file.findall("sym_id")[0].text)
            symop_interface.append(pisa_file.findall("symop")[0].text)
            area_interface.append(pisa_file.findall("area")[0].text)
            delta_g_interface.append(pisa_file.findall("delta_g")[0].text)
            nhbonds_interface.append(pisa_file.findall("nhbonds")[0].text)
            nsbridges_interface.append(pisa_file.findall("nsbridges")[0].text)
            ndsbonds_interface.append(pisa_file.findall("ndsbonds")[0].text)
             
        if len(ser_no_interface)>0:
          interfaceFold = results.addFold ( label='Saved Interfaces', initiallyOpen=True )
          interfaceTable = interfaceFold.addTable(transpose=False)
          interfaceTable.addData ( title = "Ser. No.", data = ser_no_interface )
          interfaceTable.addData ( title = "Sym. Id.", data = sym_id_interface )
          interfaceTable.addData ( title = "Sym. Op.", data = symop_interface)
          interfaceTable.addData ( title = "Area", data = area_interface)
          interfaceTable.addData ( title = "Delta G", data = delta_g_interface)
          interfaceTable.addData ( title = "H-Bonds", data = nhbonds_interface)
          interfaceTable.addData ( title = "Salt Bridges", data = nsbridges_interface)
          interfaceTable.addData ( title = "DS Bonds", data = ndsbonds_interface)

        if len(ser_no_assembly)>0:
          assemblyFold = results.addFold ( label='Saved Assemblies', initiallyOpen=True )
          assemblyTable = assemblyFold.addTable(transpose=False)
          assemblyTable.addData ( title = "Ser. No.", data = ser_no_assembly )
          assemblyTable.addData ( title = "Comp.", data = composition_assembly )
          assemblyTable.addData ( title = "Olig. State", data = olig_state_assembly )
          assemblyTable.addData ( title = "Delta G", data = delta_g_assembly )
          assemblyTable.addData ( title = "ASA", data = asa_assembly )
          assemblyTable.addData ( title = "BSA", data = bsa_assembly )

        #self.drawPictures()
        #self.drawWarnings()

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.findall('.//Warnings/Warning')
        warningsFolder = parent.addFold(label='Warnings', initiallyOpen=True)
        if len(warnings)>0:
            for warning in warnings:
                warningsFolder.addText(text=warning.text,style='color:red;')
                warningsFolder.append('<br/>')
        else:
            warningsFolder.addText(text='No warnings from qtpisa')
        return

    def drawPictures(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=True):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        baseScenePath = os.path.join(ccp4i2_root,'wrappers','qtpisa','script','qtpisa.scene.xml')

        pictureFold = parent.addFold(label='Pictures', initiallyOpen=initiallyOpen,brief='Picture')
        pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='450px', tableWidth='260px', contentWidth='450px')
        jobDirectory = jobInfo['fileroot']

        if self.jobInfo and "filenames" in self.jobInfo and "XYZOUT" in self.jobInfo["filenames"]:
          i = 0
          for fname in self.jobInfo['filenames']["XYZOUT"]:
             baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath)
             et = ET.ElementTree(baseSceneXML)
             filename_element = et.findall(".//scene/data/MolData/filename")[0]
             del filename_element.attrib["database"]
             filename_element.text = fname
             ET.indent(et)
             #print(ET.tostring(et))
             sceneFilePath = os.path.join(jobDirectory,'qtpisa_scene'+str(i)+'.scene.xml')
             ET.indent(et)
             et.write(sceneFilePath)
             pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of structure '+str(i+1))
             i = i + 1
