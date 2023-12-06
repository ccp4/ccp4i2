from report.CCP4ReportParser import *
import sys

class coordinate_selector_report(Report):
    TASKNAME = 'coordinate_selector'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus, **kw)

        self.addDiv(style="clear:both;")
        summaryFold = self.addFold(label="Summary",brief="Summary", initiallyOpen=True)

        chainList = []
        polymerList = []
        polymerTypeList = []
        nonPolymerList = []
        ligandList = []
        metalList = []
        WaterList = []
        UnknownList = []
        haveThis = False
        if len(xmlnode.findall('.//ModelComposition'))>0:
            compNode = xmlnode.findall('.//ModelComposition')[0]
            for node in compNode.findall("Chain"):
                 chainList.append(node.get('id'))
                 haveThis = True
                 polymerList.append(str(node.findall('PolymerLength')[0].text))
                 polyType = node.findall('PolymerTypes')[0].text
                 if polyType:
                     polymerTypeList.append(str(polyType))
                 else:
                     polymerTypeList.append("-")
                 nonPolymerList.append(str(node.findall('NonPolymerLength')[0].text))
                 ligandList.append(str(node.findall('LigandLength')[0].text))
                 metalList.append(str(node.findall('MetalLength')[0].text))
                 WaterList.append(str(node.findall('WaterLength')[0].text))
                 UnknownList.append(str(node.findall('UnknownLength')[0].text))

        if haveThis:
              t = summaryFold.addTable(transpose=True)
              t.addData(title="",data=chainList)
              t.addData(title="Polymers",data=polymerList)
              t.addData(title="Polymer Type(s)",data=polymerTypeList)
              t.addData(title="Non-polymers",data=nonPolymerList)
              t.addData(title="(Ligands)",data=ligandList)
              t.addData(title="(Metals)",data=metalList)
              t.addData(title="Water",data=WaterList)
              t.addData(title="Unknown",data=UnknownList)

        numModels = len(xmlnode.findall('.//ModelComposition/Model'))
        summaryFold.addDiv(style="clear:both;")
        if numModels > 1:
            summaryFold.addText(text="This file contains "+str(numModels)+" models")

        self.addDiv(style="clear:both;")

        numLigands = len(xmlnode.findall('.//ModelComposition/Ligand'))
        if numLigands > 0:
            ligandList = [str(x.attrib['id']) for x in xmlnode.findall('.//ModelComposition/Ligand')]
            ligandsFold = self.addFold(label="Ligands",brief="Ligands", initiallyOpen=False)
            t = ligandsFold.addTable(transpose=False)
            t.addData(title="Ligands",data=ligandList)

        self.addDiv(style="clear:both;")

        numMetals = len(xmlnode.findall('.//ModelComposition/Metal'))
        if numMetals > 0:
            ligandList = [str(x.attrib['id']) for x in xmlnode.findall('.//ModelComposition/Metal')]
            ligandsFold = self.addFold(label="Metals",brief="Metals", initiallyOpen=False)
            t = ligandsFold.addTable(transpose=False)
            t.addData(title="Metals",data=ligandList)
