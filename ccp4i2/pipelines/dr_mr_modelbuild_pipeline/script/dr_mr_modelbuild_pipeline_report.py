from __future__ import print_function
import sys
import copy
import json
from report.CCP4ReportParser import *
from ccp4i2.wrappers.refmac_i2.script.refmac_report import refmac_report
#from lxml import etree
import xml.etree.ElementTree as etree
from ccp4i2.wrappers.sheetbend.script.sheetbend_report import sheetbend_report
from ccp4i2.wrappers.pointless.script import pointless_report
from ccp4i2.wrappers.aimless.script import aimless_report
from ccp4i2.wrappers.ctruncate.script import ctruncate_report
from ccp4i2.wrappers.phaser_analysis.script import phaser_analysis_report
from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import *
from ccp4i2.pipelines.buccaneer_build_refine_mr.script import buccaneer_build_refine_mr_report
from ccp4i2.pipelines.aimless_pipe.script import aimless_pipe_report
from ccp4i2.wrappers.modelcraft.script import modelcraft_report

class MyRefmacReport(refmac_report):
    def addSummary(self, xmlnode=None, parent=None, withTables=True):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        
        summaryFold = parent.addFold(label='Summary of refinement', brief='Summary', initiallyOpen=True)
        import uuid
        uuid._uuid_generate_time = None
        uuid._uuid_generate_random = None
        if sys.version_info >= (3,0):
            uuid_str = uuid.uuid4().hex
        else:
            uuid_str = uuid.uuid4().get_hex()
        
        self.addScrollableDownloadableTable1(parent=summaryFold,internalId=uuid_str)
        self.addProgressGraph(parent=summaryFold)
        if withTables: self.addTables(parent=summaryFold)

class dr_mr_modelbuild_pipeline_report(Report):
  TASKNAME = 'dr_mr_modelbuild_pipeline'
  RUNNING = True
  SEPARATEDATA=True

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

  def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
    self.xmlnode = xmlnode
    Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

    #import traceback
    #traceback.print_stack()

    aimlessdone = len(xmlnode.findall('.//AIMLESS_PIPE'))>0
    molrepdone = len(xmlnode.findall('.//MolrepResult'))>0
    sheetbenddone = len(xmlnode.findall('.//SheetbendResult'))>0
    refmacdone = len(xmlnode.findall('.//REFMAC'))>0
    acorndone = len(xmlnode.findall('.//acorn'))>0
    buccdone = len(xmlnode.findall('.//BuccaneerBuildRefineResult'))>0
    modelcraftdone = len(xmlnode.findall('.//ModelCraft'))>0

    molrepdoneNodes = xmlnode.findall('.//MolrepResult')
    sheetbendNodes = xmlnode.findall('.//SheetbendResult')
    refxml = xmlnode.findall(".//REFMAC")
    buccxml = xmlnode.findall(".//BuccaneerBuildRefineResult")
    mcxml = xmlnode.findall(".//ModelCraft")
    acornNodes = xmlnode.findall(".//acorn")
    mrbumpNodes = xmlnode.findall(".//mrbump_model_prep")

    self.enantio = False
    self.enantiomorphs = []

    if len(self.xmlnode.findall('.//AIMLESS_PIPE/POINTLESS/SolutionWarning')) > 0:
        for el in self.xmlnode.findall('.//AIMLESS_PIPE/POINTLESS/SolutionWarning'):
            if "enantiomorph".lower() in el.text.lower():
                self.enantio = True
                break

    if self.enantio:
        best_sg = self.xmlnode.findall('.//AIMLESS_PIPE/POINTLESS/BestSolution/GroupName')[0].text
        best_sg_prob = float(self.xmlnode.findall('.//AIMLESS_PIPE/POINTLESS/BestSolution/TotalProb')[0].text)
        sglist = self.xmlnode.findall('.//AIMLESS_PIPE/POINTLESS/SpacegroupList/Spacegroup')
        for sg in sglist:
            sgname, total_prob = (sg.findall('SpacegroupName')[0].text.strip(), float(sg.findall('TotalProb')[0].text))
            if abs(total_prob-best_sg_prob)<1e-5:
                self.enantiomorphs.append(sgname)
        print("possible enantiomorphic space groups:",self.enantiomorphs)
        self.otherSG = list(set(self.enantiomorphs) - set([best_sg]))[0]
        print("'Other' space group is",self.otherSG)

    if len(sheetbendNodes)>1:
        treeSG1 = copy.deepcopy(xmlnode)
        treeSG1.findall(".//SheetbendResult")[0].attrib["myid"] = "sg1"
        treeSG1.findall(".//SheetbendResult")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg1",treeSG1)
        sheetbendsg1 = treeSG1.findall(".//SheetbendResult")[0]
        treeSG2 = copy.deepcopy(xmlnode)
        treeSG2.findall(".//SheetbendResult")[0].attrib["myid"] = "sg1"
        treeSG2.findall(".//SheetbendResult")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg2",treeSG2)
        sheetbendsg2 = treeSG2.findall(".//SheetbendResult")[0]
    elif len(sheetbendNodes)>0:
        sheetbendsg1 = sheetbendNodes[0]

    if len(acornNodes)>1:
        treeSG1 = copy.deepcopy(xmlnode)
        treeSG1.findall(".//acorn")[0].attrib["myid"] = "sg1"
        treeSG1.findall(".//acorn")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg1",treeSG1)
        acornsg1 = treeSG1.findall(".//acorn")[0]
        treeSG2 = copy.deepcopy(xmlnode)
        treeSG2.findall(".//acorn")[0].attrib["myid"] = "sg1"
        treeSG2.findall(".//acorn")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg2",treeSG2)
        acornsg2 = treeSG2.findall(".//acorn")[0]
    elif len(acornNodes)>0:
        acornsg1 = acornNodes[0]

    if len(refxml)>1:
        treeSG1 = copy.deepcopy(xmlnode)
        treeSG1.findall(".//REFMAC")[0].attrib["myid"] = "sg1"
        treeSG1.findall(".//REFMAC")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg1",treeSG1)
        refsg1 = treeSG1.findall(".//REFMAC")[0]
        treeSG2 = copy.deepcopy(xmlnode)
        treeSG2.findall(".//REFMAC")[0].attrib["myid"] = "sg1"
        treeSG2.findall(".//REFMAC")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg2",treeSG2)
        refsg2 = treeSG2.findall(".//REFMAC")[0]
    elif len(refxml)>0:
        refsg1 = refxml[0]

    if len(buccxml)>1:
        treeSG1 = copy.deepcopy(xmlnode)
        treeSG1.findall(".//BuccaneerBuildRefineResult")[0].attrib["myid"] = "sg1"
        treeSG1.findall(".//BuccaneerBuildRefineResult")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg1",treeSG1)
        buccsg1 = treeSG1.findall(".//BuccaneerBuildRefineResult")[0]
        treeSG2 = copy.deepcopy(xmlnode)
        treeSG2.findall(".//BuccaneerBuildRefineResult")[0].attrib["myid"] = "sg1"
        treeSG2.findall(".//BuccaneerBuildRefineResult")[1].attrib["myid"] = "sg2"
        self.remove_siblings_of("myid","sg2",treeSG2)
        buccsg2 = treeSG2.findall(".//BuccaneerBuildRefineResult")[0]
    elif len(buccxml)>0:
        buccsg1 = buccxml[0]

    summaryDiv = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

    originalText = ""
    otherText = ""

    if self.enantio:
        summaryDiv.addText(text='The input data have been processed in an entiomorphic space group.')
        summaryDiv.addDiv(style="width:100%;border-width: 0px; clear:both; margin:0px; padding:0px;")
        if jobStatus in ['Running','Running remotely']:
            summaryDiv.addText(text='Because of this the pipeline is being run for both enantiomers: '+best_sg+" and "+self.otherSG)
        else:
            summaryDiv.addText(text='Because of this the pipeline has been run for both enantiomers: '+best_sg+" and "+self.otherSG)
            summaryDiv.addDiv(style="width:100%;border-width: 0px; clear:both; margin:0px; padding:0px;")
            rfactors = xmlnode.findall('.//BuccaneerBuildRefineResult/FinalStatistics/r_factor')
            r0 = rfactors[0].text
            r1 = rfactors[1].text
            summaryDiv.addText(text='The final R-factor from model building in '+best_sg+' is '+r0)
            clearingDiv = summaryDiv.addDiv(style="clear:both;")
            summaryDiv.addText(text='The final R-factor from model building in '+self.otherSG+' is '+r1)
            clearingDiv = summaryDiv.addDiv(style="clear:both;")
            if float(r0) < float(r1):
                summaryDiv.addText(text='Based on this the space group '+best_sg+' appears to have given the best solution')
                originalText = " (space group "+best_sg+" BEST)"
                otherText = " (space group "+self.otherSG+")"
            else:
                summaryDiv.addText(text='Based on this the space group '+self.otherSG+' appears to have given the best solution')
                originalText = " (space group "+best_sg+")"
                otherText = " (space group "+self.otherSG+" BEST)"
            clearingDiv = summaryDiv.addDiv(style="clear:both;")

    if aimlessdone and jobStatus in ['Running','Running remotely']:
        pointlessxml = xmlnode.findall(".//AIMLESS_PIPE/POINTLESS")[0]
        havePointlessReport = (pointlessxml != None) and (len(pointlessxml) > 0)
        if not molrepdone:
            summaryfold = summaryDiv.addFold(label='Aimless Pipeline Key summary', brief='Aimless', initiallyOpen=True)
        else:
            summaryfold = summaryDiv.addFold(label='Aimless Pipeline Key summary', brief='Aimless', initiallyOpen=False)

        if havePointlessReport:
            self.pointlessreport = pointless_report.pointless_report(xmlnode=pointlessxml, jobStatus='nooutput')
            fail = self.pointlessreport.Errors(summaryfold)   #  report any fatal errors
            if fail:
                return
            self.pointlessreport.keyText(summaryfold)

        aimlessxml = None
        if len(self.xmlnode.findall(".//AIMLESS_PIPE/AIMLESS"))>0:
            aimlessxml = self.xmlnode.findall(".//AIMLESS_PIPE/AIMLESS")[0]
        self.twoaimless = False
        if len(self.xmlnode.findall(".//AIMLESS_PIPE/AIMLESS1"))>0:
            self.twoaimless = True
            # XML from 1st Aimless, if it was run twice. Not used yet
            aimlessFirstxml = self.xmlnode.findall(".//AIMLESS_PIPE/AIMLESS1")[0]
        haveAimlessReport = aimlessxml is not None and (len(aimlessxml) > 0)
        if haveAimlessReport:
            self.aimlessreport = aimless_report.aimless_report(xmlnode=aimlessxml, jobStatus='nooutput',jobNumber=self.jobNumber)
            self.aimlessreport.setup()
            self.aimlessreport.keyText(None, summaryfold)

        self.phaseranalysisxml = None
        if len(self.xmlnode.findall(".//AIMLESS_PIPE/PHASER_ANALYSES"))>0:
          self.phaseranalysisxml = self.xmlnode.findall(".//AIMLESS_PIPE/PHASER_ANALYSES")[0]
        self.phaseranalysis1xml = None
        if len(self.xmlnode.findall(".//AIMLESS_PIPE/PHASER_ANALYSES1"))>0:
          self.phaseranalysis1xml = self.xmlnode.findall(".//AIMLESS_PIPE/PHASER_ANALYSES1")[0]
        havePhaserReport = self.phaseranalysisxml is not None

        ctruncatexmlsnode = self.xmlnode.findall(".//AIMLESS_PIPE/CTRUNCATES")[0]
        haveCtruncateReport = (ctruncatexmlsnode != None) and (len(ctruncatexmlsnode) > 0)
        if haveCtruncateReport:
            self.ctruncatexmlnodelist = ctruncatexmlsnode.findall("CTRUNCATE")

            self.phaserxmlnodelist = None
            self.phaserreports = None
            self.phasertNCS = None
            if havePhaserReport:
              self.phaserxmlnodelist = self.phaseranalysisxml.findall('PHASER_ANALYSIS')
              self.phaserreports = []
              for phaserxml in self.phaserxmlnodelist:
                self.phaserreports.append(phaser_analysis_report.phaser_analysis_report(
                  xmlnode = phaserxml, jobStatus='nooutput'))
              try:
                  self.phasertNCS = self.phaserreports[0].getdatavalue('tNCS')
              except:
                  self.phasertNCS = None

            self.ctruncatereports = None
            if (self.ctruncatexmlnodelist != None):
              self.ctruncatereports = []
              for ctruncatexmlnode in self.ctruncatexmlnodelist:
                self.ctruncatereports.append(
                  ctruncate_report.ctruncate_report(xmlnode=ctruncatexmlnode, jobStatus='nooutput') )
              #print("!! ctruncatereports !! ", self.ctruncatereports)

            if self.ctruncatereports != None:
              self.ctruncatereports[0].keyText(summaryfold, phasertncs=self.phasertNCS)

    elif aimlessdone:
        pointlessxml = xmlnode.findall(".//AIMLESS_PIPE/POINTLESS")[0]
        summaryfold = summaryDiv.addFold(label='Aimless Pipeline', brief='Aimless', initiallyOpen=False)
        self.aimless_report = aimless_pipe_report.aimless_pipe_report(xmlnode=self.xmlnode.findall(".//AIMLESS_PIPE")[0], jobStatus=self.jobStatus,jobNumber=self.jobInfo["jobnumber"]+".1")
        self.aimless_report.defaultReport(parent=summaryfold)
        clearingDiv = summaryfold.addDiv(style="clear:both;")

    if len(mrbumpNodes)>0:
        if not molrepdone:
            summaryfold = summaryDiv.addFold(label='MrBUMP model search and preparation'+originalText, brief='MrBUMP', initiallyOpen=True)
        else:
            summaryfold = summaryDiv.addFold(label='MrBUMP model search and preparation'+originalText, brief='MrBUMP', initiallyOpen=False)
        selectString = ".//mrbump_model_prep/model"
        summaryfold.addText(text='Summary of the MrBUMP MR model search. The highest ranked file will be used as input to the following molecular replacement (MOLREP) job.')
        summaryfold.addDiv(style="clear:both;")
        mrbumpModel  = xmlnode.findall('.//mrbump_model_prep/bestModel')[0].text
        summaryfold.addText(text='The chosen model is '+mrbumpModel)
        summaryfold.addDiv(style="width:100%;border-width: 0px; clear:both; margin:0px; padding:0px;")
        progressTable = summaryfold.addTable(select=selectString, style="height:250px; width:260px;float:left;")
        progressTable.addData(title="Rank", select="rank")
        progressTable.addData(title="Chain source", select="chainSource")
        progressTable.addData(title="Coverage", select="coverage",expr='"%.2f" % float(x)')
        progressTable.addData(title="eLLG", select="eLLG")
        progressTable.addData(title="Resolution", select="resolution")
        progressTable.addData(title="Score", select="score")
        progressTable.addData(title="Sequence identity", select="seqID")
        progressTable.addData(title="Search source", select="source")
        progressTable.addData(title="Model type", select="type")

    if molrepdone:
        if not sheetbenddone:
            summaryfold = summaryDiv.addFold(label='Molecular replacement'+originalText, brief='MR', initiallyOpen=True)
        else:
            summaryfold = summaryDiv.addFold(label='Molecular replacement'+originalText, brief='MR', initiallyOpen=False)

        graph = summaryfold.addFlotGraph( xmlnode=molrepdoneNodes[0], title="Best TF peak vs RF peak No" , select="RFpeaks/RFpeak" )
        graph.addData (title="RF_peak_No"  , select="RF" )
        graph.addData (title="Score" , select="Score" )
        graph.addData (title="Tf_sig" , select="TF_sig" )
        graph.addPlot ( plot = '''<plot>
<title>Best TF peak score vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="2">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>blue</colour>
</plotline>
</plot> ''' )
        graph.addPlot ( plot = '''<plot>
<title>TF/sig(TF) vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="3">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>red</colour>
</plotline>
</plot> ''' )

        
    if sheetbenddone:
        summaryDivSheet = self.addDiv(xmlnode=sheetbendsg1,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        sheetbendFold = summaryDivSheet.addFold(label='Shift field refinement'+originalText,initiallyOpen=not refmacdone,brief='Shift field')
        sheetbendReport = sheetbend_report(xmlnode=sheetbendsg1, jobStatus='nooutput')
        sheetbendReport.defaultReport(parent=sheetbendFold)
        clearingDiv = sheetbendFold.addDiv(style="clear:both;")

    if refmacdone:
      summaryDivRef = self.addDiv(xmlnode=refsg1,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      if not acorndone and not buccdone and not modelcraftdone:
          summaryfold = summaryDivRef.addFold(label='Refinement'+originalText, brief='Refinement', initiallyOpen=True)
      else:
          summaryfold = summaryDivRef.addFold(label='Refinement'+originalText, brief='Refinement', initiallyOpen=False)

      self.ref_report = MyRefmacReport(xmlnode=refsg1, jobStatus=self.jobStatus)
      self.ref_report.addSummary(xmlnode=refsg1,parent=summaryfold)
      clearingDiv = summaryfold.addDiv(style="clear:both;")

#FIXME - XML PICTURE
      """
      pictureFold = summaryfold.addFold(label='Picture', initiallyOpen=False)
      pictureFold.addText(text='View of the best model')
      if jobStatus in ['Finished']:
         pic = pictureFold.addPicture(label='', sceneFile="$CCP4I2/wrappers/molrep_mr/script/molrep_mr_1.scene.xml")
      """

    if acorndone:
        summaryDivAcorn = self.addDiv( xmlnode=acornsg1,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        if not buccdone and not modelcraftdone:
            summaryfold = summaryDivAcorn.addFold(label='Phase refinement'+originalText, brief='Phase refinement', initiallyOpen=True)
        else:
            summaryfold = summaryDivAcorn.addFold(label='Phase refinement'+originalText, brief='Phase refinement', initiallyOpen=False)
        summaryfold.append("<p>Results for Acorn Run</p>")
        graph_height = 300
        graph_width = 500
        
        graph = summaryfold.addFlotGraph( xmlnode=acornsg1, title="Results by Cycle", select="RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width) )
        graph.addData (title="Cycle",  select="NCycle" )
        graph.addData (title="Cycle",  select="CorrelationCoef" )
        
        p = graph.addPlotObject()
        p.append('title', 'Correlation Coefficient by Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','Corr.Coef.')
        
        l = p.append('plotline',xcol=1,ycol=2)
        l.append('label','Correlation Coefficient')
        l.append('colour','teal')
        summaryDiv.addDiv(style="width:100%;border-width: 0px; border-color: black; clear:both; margin:0px; padding:0px;")

    if modelcraftdone:
        summaryDivMC = self.addDiv(xmlnode=mcxml[0],style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        summaryfold = summaryDivMC.addFold(label='ModelCraft autobuild'+originalText, brief='ModelCraft', initiallyOpen=len(mcxml)<2)
        directory = mcxml[0].text
        self.mcreport = modelcraft_report.modelcraft_report(jobInfo={"fileroot":directory},jobStatus=self.jobStatus)

        if self.jobStatus in ["Running", "Running remotely"]:
            self.mcreport.append("<p><b>The job is currently running.</b></p>")
        json_path = os.path.join(directory, "modelcraft.json")
        if os.path.exists(json_path):
            with open(json_path) as stream:
                self.mcreport.json = json.load(stream)
            self.mcreport.add_running_job(parent=summaryfold)
            self.mcreport.add_table(parent=summaryfold)
            self.mcreport.add_message(self.jobStatus,parent=summaryfold)
            self.mcreport.add_graph(parent=summaryfold)
        clearingDiv = summaryfold.addDiv(style="clear:both;")

    if buccdone:
        summaryDivBucc = self.addDiv(xmlnode=buccsg1,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        summaryfold = summaryDivBucc.addFold(label='Buccaneer autobuild'+originalText, brief='Buccaneer', initiallyOpen=len(buccxml)<2)
        summaryfold.append("<p>Results for Automated Model Building with Buccaneer</p>")
        self.buccreport = buccaneer_build_refine_mr_report.bucref_report(xmlnode=buccsg1, jobStatus='nooutput')
        if jobStatus in ['Running','Running remotely']:
          self.buccreport.runningText(summaryfold)
          self.buccreport.summaryTable(summaryfold)
          self.buccreport.completenessGraph(summaryfold, 400, 265)
          clearingDiv = summaryfold.addDiv(style="clear:both;")
        else:
          results = self.buccreport.addResults()
          self.buccreport.finishedText(summaryfold)
          print("Call table1(1)")
          self.buccreport.table1(summaryfold)
          self.buccreport.completenessGraph(summaryfold, 450, 300)
          clearingDiv = summaryfold.addDiv(style="clear:both;")
          self.buccreport.details(summaryfold)
          clearingDiv = summaryfold.addDiv(style="clear:both;")
          self.buccreport.alignment()
          clearingDiv = summaryfold.addDiv(style="clear:both;")
#FIXME - XML PICTURE
          #self.buccreport.picture(summaryfold)
          clearingDiv = summaryfold.addDiv(style="clear:both;")

    if len(molrepdoneNodes)>1:
            summaryDivMolrep = self.addDiv(style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            if len(sheetbendNodes)<2:
                summaryfold = summaryDivMolrep.addFold(label='Molecular replacement'+otherText, brief='MR', initiallyOpen=True)
            else:
                summaryfold = summaryDivMolrep.addFold(label='Molecular replacement'+otherText, brief='MR', initiallyOpen=False)
            graph = summaryfold.addFlotGraph( xmlnode=molrepdoneNodes[1], title="Best TF peak vs RF peak No" , select="RFpeaks/RFpeak" )
            graph.addData (title="RF_peak_No"  , select="RF" )
            graph.addData (title="Score" , select="Score" )
            graph.addData (title="Tf_sig" , select="TF_sig" )
            graph.addPlot ( plot = '''<plot>
<title>Best TF peak score vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="2">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>blue</colour>
</plotline>
</plot> ''' )
            graph.addPlot ( plot = '''<plot>
<title>TF/sig(TF) vs RF peak No</title>
<plottype>xy</plottype>
<plotline xcol="1" ycol="3">
<linestyle>.</linestyle>
<markeredgewidth>0</markeredgewidth>
<colour>red</colour>
</plotline>
</plot> ''' )

    if len(sheetbendNodes)>1:
            summaryDivSheet2 = self.addDiv(xmlnode=sheetbendsg2,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            sheetbendFold = summaryDivSheet2.addFold(label='Shift field refinement'+otherText,initiallyOpen=len(refxml)<2,brief='Shift field')
            sheetbendReport = sheetbend_report(xmlnode=sheetbendsg2, jobStatus='nooutput')
            sheetbendReport.defaultReport(parent=sheetbendFold)
            clearingDiv = sheetbendFold.addDiv(style="clear:both;")

    if len(refxml)>1:
          summaryDivRef2 = self.addDiv(xmlnode=refsg2,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
          if len(acornNodes)<2 and len(buccxml)<2:
              summaryfold = summaryDivRef2.addFold(label='Refinement'+otherText, brief='Refinement', initiallyOpen=True)
          else:
              summaryfold = summaryDivRef2.addFold(label='Refinement'+otherText, brief='Refinement', initiallyOpen=False)
          self.ref_report2 = MyRefmacReport(xmlnode=refsg2, jobStatus=self.jobStatus)
          self.ref_report2.addSummary(xmlnode=refsg2,parent=summaryfold)
          clearingDiv = summaryfold.addDiv(style="clear:both;")

    if len(acornNodes)>1:
            summaryDivAcorn2 = self.addDiv( xmlnode=acornsg2,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            if len(buccxml)<2:
                summaryfold = summaryDivAcorn2.addFold(label='Phase refinement'+otherText, brief='Phase refinement', initiallyOpen=True)
            else:
                summaryfold = summaryDivAcorn2.addFold(label='Phase refinement'+otherText, brief='Phase refinement', initiallyOpen=False)
            summaryfold.append("<p>Results for Acorn Run</p>")
            graph_height = 300
            graph_width = 500
        
            graph = summaryfold.addFlotGraph( xmlnode=acornsg2, title="Results by Cycle", select="RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width) )
            graph.addData (title="Cycle",  select="NCycle" )
            graph.addData (title="Cycle",  select="CorrelationCoef" )
        
            p = graph.addPlotObject()
            p.append('title', 'Correlation Coefficient by Cycle')
            p.append('plottype','xy')
            p.append('xintegral','true')
            p.append('xlabel','Cycle')
            p.append('ylabel','Corr.Coef.')
        
            l = p.append('plotline',xcol=1,ycol=2)
            l.append('label','Correlation Coefficient')
            l.append('colour','teal')
            summaryDivAcorn.addDiv(style="width:100%;border-width: 0px; border-color: black; clear:both; margin:0px; padding:0px;")

    if len(buccxml)>1:
            summaryDivBucc2 = self.addDiv(xmlnode=buccsg2,style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
            summaryfold2 = summaryDivBucc2.addFold(label='Buccaneer autobuild'+otherText, brief='Buccaneer', initiallyOpen=jobStatus in ['Running','Running remotely'])
            summaryfold2.append("<p>Results for Automated Model Building with Buccaneer</p>")
            self.buccreport2 = buccaneer_build_refine_mr_report.bucref_report(xmlnode=buccsg2, jobStatus='nooutput')
            if jobStatus in ['Running','Running remotely']:
              self.buccreport2.runningText(summaryfold2)
              self.buccreport2.summaryTable(summaryfold2)
              self.buccreport2.completenessGraph(summaryfold2, 400, 265)
              clearingDiv = summaryfold2.addDiv(style="clear:both;")
            else:
              results = self.buccreport2.addResults()
              self.buccreport2.finishedText(summaryfold2)
              print("Call table1(2)")
              self.buccreport2.table1(summaryfold2)
              self.buccreport2.completenessGraph(summaryfold2, 450, 300)
              clearingDiv = summaryfold2.addDiv(style="clear:both;")
              self.buccreport2.details(summaryfold2)
              clearingDiv = summaryfold2.addDiv(style="clear:both;")
              self.buccreport2.alignment()
              clearingDiv = summaryfold2.addDiv(style="clear:both;")
#FIXME - XML PICTURE
              #self.buccreport2.picture(summaryfold2)
              clearingDiv = summaryfold2.addDiv(style="clear:both;")

if __name__ == "__main__":
  import sys
  dr_mr_modelbuild_pipeline_report(xmlFile=sys.argv[1], jobId=sys.argv[2])


