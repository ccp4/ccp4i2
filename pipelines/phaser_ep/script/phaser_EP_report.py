import copy
import json
import os
import xml.etree.ElementTree as etree
from report.CCP4ReportParser import Report
from wrappers.modelcraft.script.modelcraft_report import modelcraft_report


class phaser_EP_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_EP'
    RUNNING = True

    def remove_element(self,name, value, root):
      """
      Iterates through the @root element and removes elements
      where the @name != @value.
      """
      for element in root:
          if element.attrib.get(name) != value:
              root.remove(element)

    def remove_siblings_of(self,name, value, root):
      """
      Recursively removes from the @root element all elements which (1) do
      not have @name == @value but (2) do have a sibling where @name == @value.
      """
      for element in root:
          if element.attrib.get(name) == value:
              self.remove_element(name, value, root)  # need to reiterate through element now to remove previous siblings
          if len(element):
              self.remove_siblings_of(name, value, element)
      return root

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        if jobStatus is None or jobStatus.lower() =='nooutput':
            return
        self.drawContent(jobStatus, self)

    def drawContent(self, jobStatus=None, parent=None):
        if parent is None: parent = self
        if len(self.xmlnode.findall('ShelxCD'))>0:
            shelxNode = self.xmlnode.findall('ShelxCD')[0]
            from wrappers.ShelxCDE.script.ShelxCD_report import ShelxCD_report
            shelx_report = ShelxCD_report (xmlnode=shelxNode, jobStatus='nooutput')
            self.addDiv(style='clear:both;')
            if len(shelxNode.findall('.//Shelxd'))>0:
                shelx_report.shelXCReport(self, initiallyOpen=True )
            else:
                shelx_report.shelXCReport(self, initiallyOpen=False )
                shelx_report.shelXDReport(self, initiallyOpen=False)

        treePhaser = copy.deepcopy(self.xmlnode)
        treePhaser.findall(".//PhaserEpResults")[0].attrib["myid"] = "phaser_EP_res_tree"
        self.remove_siblings_of("myid","phaser_EP_res_tree",treePhaser)
        phaserNode = treePhaser.findall('.//PhaserEpResults')[0]
        phaserDiv = self.addDiv( xmlnode=phaserNode,style="width:100%;border-width: 0px; border-color: black; clear:both; margin:1px; padding:1px;")

        if phaserNode is not None:
            from pipelines.phaser_pipeline.wrappers.phaser_EP_AUTO.script.phaser_EP_AUTO_report import phaser_EP_AUTO_report
            phaser_report = phaser_EP_AUTO_report(xmlnode=phaserNode, jobStatus='nooutput')
            phaser_report.drawContent(jobStatus=jobStatus, parent=phaserDiv)

        if len(self.xmlnode.findall('.//original/ParrotResult'))>0:
            parrotOriginalHandNode = self.xmlnode.findall('.//original/ParrotResult')[0]
            from wrappers.parrot.script.parrot_report import parrot_report
            parrotOriginalNode = etree.fromstring(etree.tostring(parrotOriginalHandNode))
            parrot_original_report = parrot_report(xmlnode=parrotOriginalNode, jobStatus='nooutput')
            parrot_original_hand = parent.addFold(label='Density modification: Original hand', initiallyOpen=False)
            parrot_original_report.defaultReport(parent=parrot_original_hand)
            parrot_original_hand.addDiv(style="clear:both;")
        if len(self.xmlnode.findall('.//inverted/ParrotResult'))>0:
            parrotInvertedHandNode = self.xmlnode.findall('.//inverted/ParrotResult')[0]
            from wrappers.parrot.script.parrot_report import parrot_report
            parrotInvertedNode = etree.fromstring(etree.tostring(parrotInvertedHandNode))
            parrot_inverted_report = parrot_report(xmlnode=parrotInvertedNode, jobStatus='nooutput')
            parrot_inverted_hand = parent.addFold(label='Density modification: Inverted hand', initiallyOpen=False)
            parrot_inverted_report.defaultReport(parent=parrot_inverted_hand)
            parrot_inverted_hand.addDiv(style="clear:both;")

        buccxml = self.xmlnode.findall('.//BuccaneerBuildRefineResult')

        buccsg1 = None
        buccsg2 = None

        if len(buccxml)>1:
            buccsg1 = etree.fromstring(etree.tostring(buccxml[0]))
            buccsg2 = etree.fromstring(etree.tostring(buccxml[1]))
        elif len(buccxml)>0:
            buccsg1 = buccxml[0]

        if buccsg1 is not None:
            buccDiv1 = self.addDiv( xmlnode=buccsg1,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            from pipelines.buccaneer_build_refine_mr.script.buccaneer_build_refine_mr_report import bucref_report
            buccaneer_original_report = bucref_report(xmlnode=buccsg1, jobStatus='nooutput')
            buccaneer_original_hand = buccDiv1.addFold(label='ModelBuilding: Original hand', initiallyOpen=False)
            buccaneer_original_report.finishedText(parent=buccaneer_original_hand)
            buccaneer_original_report.table1(parent=buccaneer_original_hand)
            buccaneer_original_report.completenessGraph(buccaneer_original_hand, 450, 300)
            buccaneer_original_hand.addDiv(style="clear:both;")

        if buccsg2 is not None:
            buccDiv2 = self.addDiv( xmlnode=buccsg2,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            from pipelines.buccaneer_build_refine_mr.script.buccaneer_build_refine_mr_report import bucref_report
            buccaneer_inverted_report = bucref_report(xmlnode=buccsg2, jobStatus='nooutput')
            buccaneer_inverted_hand = buccDiv2.addFold(label='ModelBuilding: Inverted hand', initiallyOpen=False)
            buccaneer_inverted_report.finishedText(parent=buccaneer_inverted_hand)
            buccaneer_inverted_report.table1(parent=buccaneer_inverted_hand)
            buccaneer_inverted_report.completenessGraph(buccaneer_inverted_hand, 450, 300)
            buccaneer_inverted_hand.addDiv(style="clear:both;")

        for hand in ["original", "inverted"]:
            if (node := self.xmlnode.find(f".//{hand}/ModelCraft")) is not None:
                div = self.addDiv(style="width: 100%; border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
                fold = div.addFold(label=f"Model building: {hand} hand", initiallyOpen=False)
                report = modelcraft_report(jobInfo={"fileroot": ""})  # Empty so it doesn't load yet
                if jobStatus in ["Running", "Running remotely"]:
                    fold.append("<p><b>The job is currently running.</b></p>")
                json_path = os.path.join(node.text, "modelcraft", "modelcraft.json")
                if os.path.exists(json_path):
                    with open(json_path, encoding="utf-8") as stream:
                        report.json = json.load(stream)
                    report.add_running_job(parent=fold)
                    report.add_table(parent=fold)
                    report.add_message(jobStatus, parent=fold)
                    report.add_graph(parent=fold)
                fold.addDiv(style="clear:both;")
