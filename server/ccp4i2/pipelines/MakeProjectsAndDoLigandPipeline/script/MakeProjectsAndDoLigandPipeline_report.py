import xml.etree.ElementTree as etree

from ccp4i2.core import CCP4Utils
from ccp4i2.report.CCP4ReportParser import Report


class MakeProjectsAndDoLigandPipeline_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'MakeProjectsAndDoLigandPipeline'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        if len(self.xmlnode.findall(".//Message")) > 0:
            textDiv = parent.addDiv()
            textDiv.addText(text=self.xmlnode.findall(".//Message/text()")[-1],style="font-size:150%;")
        clearingDiv = self.addDiv(style="clear:both;")
        if len(self.xmlnode.findall(".//ProjectName")) > 0:
            import os

            from rdkit import Chem
            from rdkit.Chem import AllChem, Draw
            from rdkit.Chem.Draw import MolDraw2DSVG
            newFold = parent.addFold(label="Result summaries", initiallyOpen=True)
            
            #Create a list of svg file names
            useSVG=True
            try:
                for datasetNode in self.xmlnode.findall(".//Dataset"):
                    projectName = datasetNode.findall("ProjectName/text()")[-1]
                    smilesString = datasetNode.findall("SMILES/text()")[-1]
                    svgFilename = os.path.join(self.jobInfo['fileroot'], projectName+".svg")
                    from ccp4i2.wrappers.acedrg.script import acedrg
                    mol = Chem.MolFromSmiles(smilesString)
                    confId2D = AllChem.Compute2DCoords(mol)
                    svgStructure = acedrg.svgFromMol(mol)
                    svgText = etree.tostring(svgStructure,pretty_print=True)
                    svgText = svgText.replace("<svg>","<svg:svg>").replace("</svg>","</svg:svg>").replace("<line>","<svg:line>").replace("</line>","</svg:line>").replace("<text>","<svg:text>").replace("</text>","</svg:text>").replace("<polygon>","<svg:polygon>").replace("</polygon>","</svg:polygon>").replace('xmlns="http://www.w3.org/2000/svg"','xmlns:svg="http://www.w3.org/2000/svg"')
                    with open(svgFilename,"w") as svgFile:
                        CCP4Utils.writeXML(svgFile,svgText)
            except:
                useSVG=False
            
            #Add a summary table
            tableDiv = newFold.addDiv(style="float:left;width:300px;");
            tableDiv.addText(text="Project summaries",style="font-size:125%;");
            table = tableDiv.addTable(style="float:left;")
            table.addData(title="Name",data=self.xmlnode.findall(".//ProjectName/text()"))
            table.addData(title="Res",data=self.xmlnode.findall(".//MaximumResolution/text()"))
            table.addData(title="Rfac",data=self.xmlnode.findall(".//Rfactor/text()"))
            table.addData(title="Rfree",data=self.xmlnode.findall(".//Rfree/text()"))
            if useSVG:
                table.addData(title="SMILES",data=['<img src="./'+node.text+'.svg" width="150" height="150" style="width:100;height:100;"/>' for node in self.xmlnode.findall(".//ProjectName")])
            else:
                table.addData(title="SMILES",data=self.xmlnode.findall(".//SMILES/text()"))

            clearingDiv = self.addDiv(style="clear:both;")
            newFold2 = parent.addFold(label="Warnings", initiallyOpen=True)
            for warningsNode in self.xmlnode.findall(".//Warnings"):
                newFold2.addText(text=warningsNode.text, style="font-size:125%;font-color:orange;");
