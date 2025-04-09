import os
import sys

from ccp4i2.core.CCP4ErrorHandling import Severity
from ccp4i2.report.CCP4ReportParser import Report


class bucref_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'buccaneer_build_refine_mr'
  # Flag that a 'Running' mode is supported
  RUNNING = True
  # Indicate css version that this code was written against
  CSS_VERSION = '0.1.0'
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,taskVersion=None,**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)
    if self.errorReport().maxSeverity()>Severity.WARNING:
      print('FAILED instantiating buccaneer_build_refine_mr report generator')
      self.errorReport().report()
      return
    
    # 'nooutput' mode would be used by another report class that wanted
    # to use some method(s) from this class for its own report
    print('bucref_report',jobStatus, jobStatus in ['Running','Running remotely'])
    if jobStatus is not None and jobStatus.lower() == 'nooutput':
      return
    elif jobStatus in ['Running','Running remotely']:
      self.runningText(self)
      self.summaryTable(self)
      self.completenessGraph(self, 400, 265)
      clearingDiv = self.addDiv(style="clear:both;")

    else:
      results = self.addResults()
      self.finishedText()
      self.table1(self)
      self.completenessGraph(self, 450, 300)
      clearingDiv = self.addDiv(style="clear:both;")
      self.details(self)
      self.alignment()
#FIXME - XML PICTURE
      #self.picture(self)
  
  def summaryTable(self, parent=None):
    if parent is None:
        parent=self
    if len(self.xmlnode.findall("BuildRefineCycle"))>0:
        xmlPath = './/BuildRefineCycle'
        xmlNodes = self.xmlnode.findall(xmlPath)
        if len(xmlNodes)>0:
            selectString = ".//BuildRefineCycle[1]"
            if len(xmlNodes)>2:
                selectString += " | .//BuildRefineCycle[%d]" % (len(xmlNodes)-1)
            if len(xmlNodes)>1:
                selectString += " | .//BuildRefineCycle[%d]" % len(xmlNodes)
            tableDiv = parent.addDiv(style="height:20em; width:20em; float:left; border:0px;")
            progressTable = tableDiv.addTable(select=selectString)
            progressTable.addData(title="Cycle", select="Number")
            progressTable.addData(title="Completeness by residue", select="BuccaneerResult/Final/CompletenessByResiduesBuilt",expr='round(x,2)')
            #progressTable.addData(title="Completeness by chains", select="BuccaneerResult/Final/CompletenessByChainsBuilt",expr='round(x,2)')
            progressTable.addData(title="R<sub>Work</sub>", select="RefmacResult/r_factor")
            progressTable.addData(title="R<sub>Free</sub>",   select="RefmacResult/r_free",   expr="x if float(x)>=0.0 else '-'")

  def completenessGraph(self, parent=None, graph_width=450, graph_height=300):

    graph = parent.addFlotGraph( title="Progress by build-refine iteration", select=".//BuildRefineCycle",style="height:%dpx; width:%dpx; float:right; border:0px;" % (graph_height, graph_width) )
    graph.addData (title="Cycle",  select="Number" )
    graph.addData (title="Compl<sub>res</sub>", select="BuccaneerResult/Final/CompletenessByResiduesBuilt")
    graph.addData (title="Compl<sub>chn</sub>", select="BuccaneerResult/Final/CompletenessByChainsBuilt")
    graph.addData (title="Residues Built", select="BuccaneerResult/Final/ResiduesBuilt")
    graph.addData (title="Residues Sequenced", select="BuccaneerResult/Final/ResiduesSequenced")
    graph.addData (title="Longest Fragment", select="BuccaneerResult/Final/ResiduesLongestFragment")
    graph.addData (title="R<sub>Work</sub>", select="RefmacResult/r_factor")
    graph.addData (title="R<sub>Free</sub>",  select="RefmacResult/r_free" )
    graph.addData (title="RMS<sub>Bonds</sub>", select="RefmacResult/rmsBONDx100", expr="x/100.0" )
    graph.addData (title="RMS<sub>Angles</sub>",  select="RefmacResult/rmsANGLE" ) 

    p = graph.addPlotObject()
    p.append('title','Completeness and R-factors VS cycle')
    p.append('plottype','xy')
    p.append('xintegral','true')
    p.append('xlabel','Cycle')
    p.append('ylabel','Completeness')
    p.append('yrange',rightaxis='false',min='0.0',max='1.0')
    l = p.append('plotline',xcol=1,ycol=2)
    l.append('label','Compl<sub>res</sub>')
    l.append('colour','blue')
    l = p.append('plotline',xcol=1,ycol=8,rightaxis='true')
    l.append('label','R<sub>Free</sub>')
    l.append('colour','gold')
    l = p.append('plotline',xcol=1,ycol=7,rightaxis='true')
    l.append('label','R<sub>Work</sub>')
    l.append('colour','lightblue')
    if len(p.validate())>0: print(p.validate())

    p = graph.addPlotObject()
    p.append('title','R-factors after each iteration')
    p.append('plottype','xy')
    p.append('xintegral','true')
    p.append('xlabel','Cycle')
    p.append('ylabel', 'R-factor')
    l = p.append('plotline',xcol=1,ycol=8)
    l.append('label','R<sub>Free</sub>')
    l.append('colour','gold')
    l = p.append('plotline',xcol=1,ycol=7)
    l.append('label','R<sub>Work</sub>')
    l.append('colour','lightblue')
    if len(p.validate())>0: print(p.validate())

    p = graph.addPlotObject()
    p.append('title','Residues built after each iteration')
    p.append('plottype','xy')
    p.append('xlabel','Cycle')
    p.append('xintegral','true')
    p.append('ylabel', 'Number of Residues')
    p.append('yintegral','true')
    l = p.append('plotline',xcol=1,ycol=4)
    l.append('label','Built')
    l.append('colour','blue')
    l = p.append('plotline',xcol=1,ycol=5)
    l.append('label','Sequenced')
    l.append('colour','lightblue')
    l = p.append('plotline',xcol=1,ycol=6)
    l.append('label','Longest')
    l.append('colour','green')
    if len(p.validate())>0: print(p.validate())

    p = graph.addPlotObject()
    p.append('title','Geometry after each iteration')
    p.append('plottype','xy')
    p.append('xintegral','true')
    p.append('xlabel','Cycle')
    l = p.append('plotline',xcol=1,ycol=9, rightaxis='true')
    l.append('label','RMS<sub>Bonds</sub>')
    l.append('colour','firebrick')
    l = p.append('plotline',xcol=1,ycol=10)
    l.append('label','RMS<sub>Angles</sub>')
    l.append('colour','green')
    if len(p.validate())>0: print(p.validate())

    p = graph.addPlotObject()
    p.append('title','Completeness after each iteration')
    p.append('plottype','xy')
    p.append('xintegral','true')
    p.append('xlabel','Cycle')
    p.append('ylabel', 'Completeness')
    p.append('yrange', min='0.0',max='1.0')
    l = p.append('plotline',xcol=1,ycol=2)
    l.append('label','Compl<sub>res</sub>')
    l.append('colour','blue')
    l = p.append('plotline',xcol=1,ycol=3)
    l.append('label','Compl<sub>chn</sub>')
    l.append('colour','lightblue')
    if len(p.validate())>0: print(p.validate())

  def details(self,parent=None):
    if parent is None:
        parent = self
    fold = parent.addFold(label="Detailed progress by iteration")
    table = fold.addTable(select=".//BuildRefineCycle", transpose=True, downloadable=True,id='details')
    
    table.addData ( title='Iteration', select='Number', expr='int(x)' )

    for title,select in [ [ "Completeness by residue" ,"BuccaneerResult/Final/CompletenessByResiduesBuilt" ],
                          [ "Completeness by chains" , "BuccaneerResult/Final/CompletenessByChainsBuilt"  ]]:
        table.addData(title=title,select=select,expr='round(x,2)')
    for title,select in  [[ "Number of chains"  , "BuccaneerResult/Final/ChainsBuilt"],
                          [ "Residues built" , "BuccaneerResult/Final/ResiduesBuilt" ],
                          [ "Residues sequenced" , "BuccaneerResult/Final/ResiduesSequenced"],
                          [ "Longest fragment" , "BuccaneerResult/Final/ResiduesLongestFragment" ],
                          [ "Number of fragments" , "BuccaneerResult/Final/FragmentsBuilt" ] ]:
        table.addData(title=title,select=select)
    for title,select,expr in  [[ "R<sub>Work</sub>"  , "RefmacResult/r_factor", "x"],
                          [ "R<sub>Free</sub>" , "RefmacResult/r_free","x if float(x)>0.0 else '-' " ],
                          [ "RMS<sub>Bonds</sub>" , "RefmacResult/rmsBONDx100", "round(x/100,3)"],
                          [ "RMS<sub>Angles</sub>", "RefmacResult/rmsANGLE","x" ] ]:
        table.addData(title=title,select=select,expr=expr)
            
  def table1(self,parent=None):
    if parent is None:
        parent = self
    tableDiv = parent.addDiv(style="height:30em;width:20em;float:left;border:0px;")
    table = tableDiv.addTable(select=".", transpose=True, id='table_1') 
    for title,select in  [[ "Completeness by residue" ,"FinalStatistics/CompletenessByResiduesBuilt" ],
                          [ "Completeness by chains" , "FinalStatistics/CompletenessByChainsBuilt"  ]]:
        table.addData(title=title,select=select,expr='round(x,2)')
    for title,select in  [[ "Number of chains"  , "FinalStatistics/ChainsBuilt"],
                          [ "Residues built" , "FinalStatistics/ResiduesBuilt" ],
                          [ "Residues sequenced" , "FinalStatistics/ResiduesSequenced"],
                          [ "Longest fragment" , "FinalStatistics/ResiduesLongestFragment" ],
                          [ "Number of fragments" , "FinalStatistics/FragmentsBuilt" ] ]:
        table.addData(title=title,select=select)
    for title,select,expr in  [[ "R<sub>Work</sub>"  , "FinalStatistics/r_factor", "x"],
                          [ "R<sub>Free</sub>" , "FinalStatistics/r_free","x if float(x)>0.0 else '-' " ],
                          [ "RMS<sub>Bonds</sub>" , "FinalStatistics/rmsBONDx100", "round(x/100,3)"],
                          [ "RMS<sub>Angles</sub>", "FinalStatistics/rmsANGLE","x" ] ]:
        table.addData(title=title,select=select,expr=expr)

  def picture(self,parent=None) :
    pic = parent.addPicture(label="Autobuilt structure",sceneFile="$CCP4I2/pipelines/buccaneer_build_refine_mr/script/buccaneer_1.scene.xml",id='autobuild_1')

  def runningText(self,parent=None) :
    if parent is None: parent = self
    parent.append( "<p><b>The job is currently running. Updates will be shown here for each iteration of model building and refinement.</b></p>" )
    try:
      cres = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/CompletenessByResiduesBuilt')[-1].text)
      cchn = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/CompletenessByChainsBuilt')[-1].text)
      frgb = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/FragmentsBuilt')[-1].text)
      chnb = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/ChainsBuilt')[-1].text)
      resb = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/ResiduesBuilt')[-1].text)
      ress = float(self.xmlnode.findall('BuildRefineCycle/BuccaneerResult/Final/ResiduesSequenced')[-1].text)
      parent.append( self.buccaneer_text(cres, cchn, frgb, chnb, resb, ress) )
    except Exception as e:
      parent.append( "<p>Model building results are not yet available.</p>" )
      print("ERROR buccaneer_build_refine_mr_report ",e, file=sys.stderr)
    try:
      rwrk = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/r_factor')[-1].text)
      rfre = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/r_free')[-1].text)
      rbnd = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/rmsBONDx100')[-1].text)*0.01
      rang = float(self.xmlnode.findall('BuildRefineCycle/RefmacResult/rmsANGLE')[-1].text)
      parent.append( self.refmac_text(rwrk, rfre, rbnd, rang) )
    except Exception as e:
      parent.append( "<p>Refinement results are not yet available.</p>" )
       
  def finishedText(self,parent=None) :
    if parent is None: parent = self
    try:
      best_cycle = self.xmlnode.findall('FinalStatistics/BestCycle')[-1].text
      parent.append( "<p>The final model is taken from cycle %s as this had the lowest free-R factor.</p>" % best_cycle )
      cres = float(self.xmlnode.findall('FinalStatistics/CompletenessByResiduesBuilt')[-1].text)
      cchn = float(self.xmlnode.findall('FinalStatistics/CompletenessByChainsBuilt')[-1].text)
      frgb = float(self.xmlnode.findall('FinalStatistics/FragmentsBuilt')[-1].text)
      chnb = float(self.xmlnode.findall('FinalStatistics/ChainsBuilt')[-1].text)
      resb = float(self.xmlnode.findall('FinalStatistics/ResiduesBuilt')[-1].text)
      ress = float(self.xmlnode.findall('FinalStatistics/ResiduesSequenced')[-1].text)
      parent.append( self.buccaneer_text(cres, cchn, frgb, chnb, resb, ress) )
      rwrk = float(self.xmlnode.findall('FinalStatistics/r_factor')[-1].text)
      rfre = float(self.xmlnode.findall('FinalStatistics/r_free')[-1].text)
      rbnd = float(self.xmlnode.findall('FinalStatistics/rmsBONDx100')[-1].text)*0.01
      rang = float(self.xmlnode.findall('FinalStatistics/rmsANGLE')[-1].text)
      parent.append( self.refmac_text(rwrk, rfre, rbnd, rang) )
    except Exception as e:
      print("ERROR buccaneer_build_refine_mr_report ",e,file=sys.stderr)

  def buccaneer_text(self, cres, cchn, frgb, chnb, resb, ress):
    return "<p>%d residues were built in %d fragments. Of these, %d residues were assigned to the sequence.</p><p>The number of chains is estimated to be %d. Of these chains, %5.1f%% of the residues have been built.<br/>Of the residues that were built, %5.1f%% were assigned to a chain.</p>"%(resb,frgb,ress,chnb,100*cchn,100*cres)

  def refmac_text(self, rwrk, rfre, rbnd, rang):
    if rwrk > 0.5:    s = "the model is very incomplete or wrong"
    elif rwrk > 0.4:  s = "the model is substantially incomplete and may contain incorrect regions"
    elif rwrk > 0.35: s = "the model is likely to contain correct regions but requires further work"
    else:             s = "the model is approaching completion"
    return "<p>The refinement R-factor is %5.2f, and the free-R factor is %5.2f. The RMS bond deviation is %5.3f A.<br/>On the basis of the refinement statistics, %s.</p>"%(rwrk,rfre,rbnd,s)

  def alignment(self):
    alignChainNodes = self.xmlnode.findall('.//AlignChain')
    fold = self.addFold(label="Alignments for model and AU content sequences",initiallyOpen=False)
    fold.append("<p>Showing alignment of each chain in model v. best match to input sequence chain</p>")
    for alignChainNode in alignChainNodes:
      chainIdNode = alignChainNode.findall('ChainId')[0]
      fold.append("<p>Alignment for chain "+chainIdNode.text+'</p>')
      fold.addPre(text = alignChainNode.findall('Alignment')[0].text)

def test(xmlFile=None,jobId=None,reportFile=None):
  if reportFile is None:
    if xmlFile is not None:
      reportFile = os.path.join(os.path.split(xmlFile)[0],'report.html')
    else:
      reportFile = os.path.join(os.getcwd(),'report.html')
  r = bucref_report(xmlFile=xmlFile,jobId=jobId)
  r.as_html_file(reportFile)
  if len(r.errorReport())>0: print('ERRORS:',r.errorReport())
  # Temporary hard-wire to save typing..
  #r = bucref_report(xmlFile='/Users/lizp/Desktop/test_projects/t10/CCP4_JOBS/job_2/program.xml',jobId='1b86b099978d11e2822c3c0754185dfb')
  #r.as_html_file('/Users/lizp/Desktop/test_projects/t10/CCP4_JOBS/job_2/report.html')
