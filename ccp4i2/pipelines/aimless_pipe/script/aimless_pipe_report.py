import sys

from ....report.CCP4ReportParser import Report
from ....wrappers.aimless.script import aimless_report
from ....wrappers.ctruncate.script import ctruncate_report
from ....wrappers.phaser_analysis.script import phaser_analysis_report
from ....wrappers.pointless.script import pointless_report
from .aimless_pipe_utils import colourText, displayFileList


class aimless_pipe_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'aimless_pipe'
  RUNNING = True
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,style="overflow:auto;",**kw)

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    if kw.get('jobStatus',None) is not None and kw.get('jobStatus').lower() == 'nooutput':
      return
    elif kw.get('jobStatus',None) is not None and kw.get('jobStatus').lower().count('running'):
      self.pointlessReport()
    else:
      self.defaultReport()

  # - - - - - - - - - - - - - - - - - - - - - - - 
  def defaultReport(self, parent=None):
    if parent is None: parent = self
##    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,style="width:1280px;overflow:auto;",**kw)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Instantiate aimless and pointless and ctruncate reports from which to cherry pick
    # The nooutput flag makes sure thay don't do anything

    # First check if any parts of the report exist

    errorspresent = False
    pointlessFatalErrors = False
    aimlessFatalErrors = False
    
    #  0) PIPELINE_ERROR  fatal error somewhere
    if len(self.xmlnode.findall('PIPELINE_ERROR'))>0:
      errorspresent = True
      errormessage =  self.xmlnode.findall('PIPELINE_ERROR')[0].text
      #print("errormessage", errormessage)
      errorDiv = parent.addDiv(
        style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
      errorDiv.addText(text='FATAL ERROR',
                       style='font-weight:bold; font-size:150%; color:red;')
      errorDiv.append(errormessage)

      if len(self.xmlnode.findall('AIMLESS/DisasterMessage'))>0:
        dmessage = self.xmlnode.findall('AIMLESS/DisasterMessage')[0].text
        errorDiv.append(dmessage)
        errorDiv.append("The pipeline after scaling with Aimless has been abandoned")
        errorDiv.append("You should examine carefully the statistics below, your images and the integration process")
      
    #  1) IMPORT_LOG for mmCIF files
    importlogxml = None
    if len(self.xmlnode.findall('IMPORT_LOG'))>0:
      importlogxml = self.xmlnode.findall('IMPORT_LOG')[0]

    #  2) POINTLESS
    pointlessxml = self.xmlnode.findall("POINTLESS")[0]
    havePointlessReport = (pointlessxml != None) and (len(pointlessxml) > 0)

    #  3) AIMLESS
    #     Note: if Aimless has been run twice, for automatic resolution cutoff,
    #     then there will be two AIMLESS elements:
    #     the first in the XML file is the final (post-cutoff) one,
    #     the second is from the first run
    aimlessxml = None
    if len(self.xmlnode.findall("AIMLESS"))>0:
      aimlessxml = self.xmlnode.findall("AIMLESS")[0]
    self.twoaimless = False
    if len(self.xmlnode.findall("AIMLESS1"))>0:
      self.twoaimless = True
      # XML from 1st Aimless, if it was run twice. Not used yet
      aimlessFirstxml = self.xmlnode.findall("AIMLESS1")[0]
    haveAimlessReport = aimlessxml is not None and (len(aimlessxml) > 0)

    #  4) phaser_analysis
    self.phaseranalysisxml = None
    if len(self.xmlnode.findall("PHASER_ANALYSES"))>0:
      self.phaseranalysisxml = self.xmlnode.findall("PHASER_ANALYSES")[0]
    self.phaseranalysis1xml = None
    if len(self.xmlnode.findall("PHASER_ANALYSES1"))>0:
      self.phaseranalysis1xml = self.xmlnode.findall("PHASER_ANALYSES1")[0]
    havePhaserReport = self.phaseranalysisxml is not None
    
    #  5) CTRUNCATE
    ctruncatexmlsnode = self.xmlnode.findall("CTRUNCATES")
    if len(ctruncatexmlsnode) == 0:
      haveCtruncateReport = False
      ctruncatexmlsnode = None
    else:
      ctruncatexmlsnode = self.xmlnode.findall("CTRUNCATES")[0]
      haveCtruncateReport = True

    # Empty XML file
    if not havePointlessReport and not haveAimlessReport and not haveCtruncateReport:
      parent.addText(text='No Report data available, probably due to failure of job',
                     style='color:red;')
      return

    self.merged = False   # files are usually unmerged
    if havePointlessReport:
      self.pointlessreport = \
          pointless_report.pointless_report(xmlnode=pointlessxml, jobStatus='nooutput')
      self.pointlessreport.setFileRoot(self.fileroot) # pass on fileroot
      #print("self.pointlessreport", self.pointlessreport)
      self.merged = self.pointlessreport.isMerged()  # in case all files are merged

    if (haveAimlessReport):
      self.aimlessreport = \
               aimless_report.aimless_report(xmlnode=aimlessxml, jobStatus='nooutput',jobNumber=self.jobNumber,projectId=self.projectId)
      self.aimlessreport.setFileRoot(self.fileroot) # pass on fileroot
    
    self.phaserxmlnodelist = None
    self.phaserreports = None
    self.phasertNCS = None
    self.phaserTwin = None
    self.phaserFail = False
    self.phaserErrorMessage = None
    if havePhaserReport:
      self.phaserxmlnodelist = self.phaseranalysisxml.findall('PHASER_ANALYSIS')
      self.phaserreports = []
      for phaserxml in self.phaserxmlnodelist:
        self.phaserreports.append(phaser_analysis_report.phaser_analysis_report(
          xmlnode = phaserxml, jobStatus='nooutput'))

      # check for Phaser fail
      for phaserreport in self.phaserreports:
        if phaserreport.hasfailed():
          self.phaserFail = True
          self.phaserErrorMessage = phaserreport.errorMessage
          havePhaserReport = False
          self.phaserxmlnodelist = None
        
    if havePhaserReport:
      self.phasertNCS = None
      self.phaserTwin = 0
      for phaserreport in self.phaserreports:
        # tNCS "True" or "False", True trumps False
        phtncs = phaserreport.getdatavalue('tNCS')
        if self.phasertNCS is None:
          self.phasertNCS = phtncs
        elif self.phasertNCS != 'True':
          self.phasertNCS = phtncs

        # Twin = 0 untwinned, +2 probably twinned, +1 possible twin
        #  Keep maximum value
        self.phaserTwin = max(self.phaserTwin, phaserreport.istwin())
      
    self.ctruncatexmlnodelist = None
    self.ctruncatereports = None
    if haveCtruncateReport:
      self.ctruncatexmlnodelist = ctruncatexmlsnode.findall("CTRUNCATE")

      if (self.ctruncatexmlnodelist != None):
        self.ctruncatereports = []
        for ctruncatexmlnode in self.ctruncatexmlnodelist:
          self.ctruncatereports.append(
            ctruncate_report.ctruncate_report(xmlnode=ctruncatexmlnode, jobStatus='nooutput') )
        #print("!! ctruncatereports !! ", self.ctruncatereports)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Process fatal error messages, if any
    
    if havePointlessReport and self.pointlessreport.isFatalError():
        errorspresent = True
        pointlessFatalErrors = True

    if haveAimlessReport and self.aimlessreport.isFatalError():
        errorspresent = True
        aimlessFatalErrors = True

    self.numberofdatasets = 1
    if (haveAimlessReport):
      self.numberofdatasets = self.aimlessreport.numberOfDatasets()

    if errorspresent:
      if pointlessFatalErrors:
        errorDiv = parent.addDiv(
          style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
        errorDiv.append('<br/>')
        errorDiv.addText(text='Errors from Pointless', style='font-size:120%;color:red;')
        self.pointlessreport.Errors(errorDiv)   #  report any fatal errors
        havePointlessReport = False
      if aimlessFatalErrors:
        errorDiv = parent.addDiv(
          style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
        errorDiv.append('<br/>')
        errorDiv.addText(text='Errors from Aimless', style='font-size:120%;color:red;')
        self.aimlessreport.Errors(errorDiv)   #  report any fatal errors
        self.displayLogFiles(parent)
        haveAimlessReport = False
        
    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    #   Overall summaries, layout different if > 1 dataset
    if (haveAimlessReport):

      self.aimlessreport.setup()

      self.datasetIndex = datasetIndex(aimlessxml)
      if self.numberofdatasets > 1:
        # For multiple datasets, setup indexing of datasets for phaser_analysis
        # and ctruncate, which may be in a different order to those in Aimless
        self.datasetIndex.indexCtruncate(self.ctruncatexmlnodelist)
        if self.phaserxmlnodelist is not None and not self.phaserFail:
          self.datasetIndex.indexPhaser(self.phaserxmlnodelist)

      # Autocutoff block, if present
      self.autocutoffnode = None
      if len(self.xmlnode.findall('Autocutoff'))>0:
        print("** Autocutoff node found")
        self.autocutoffnode = self.xmlnode.findall('Autocutoff')[0]

      # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
      # Extract resolution information into list for each dataset
      # Store back in aimless_report
      self.estimatesofresolution = self.aimlessreport.getEstimatesofResolution()
      if self.estimatesofresolution is not None:
        self.addtoEstimatesofResolution()  # add autocutoff and some Phaser things
        self.aimlessreport.storeEstimatesofResolution(self.estimatesofresolution)
    
    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Some headers if present
    if havePointlessReport:
      fail = self.pointlessreport.Errors(parent)   #  report any fatal errors
      if fail:
        return

    if importlogxml is not None:
      topDiv = parent.addDiv(style="font-size:90%;border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;")
      mmcifDiv = topDiv.addDiv(style="padding:5px;background-color:#CCEEFF;")
      self.mmcifreport(importlogxml, mmcifDiv)

    # Summary of summaries
    summaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    summaryfold = summaryDiv.addFold(label='Key summary', brief='Headline',
                                     initiallyOpen=True)

    if havePointlessReport:
      if self.merged:
        self.pointlessreport.keyTextMerged(summaryfold)
      else:
        self.pointlessreport.keyText(summaryfold)

    if haveAimlessReport:
      if self.merged:
        phaserinfo = None  # maybe add this
        self.aimlessreport.keyTextMerged(summaryfold, phaserinfo)
      else:
        # Pass Autocutoff stuff to Aimless report
        self.aimlessreport.keyText(self.autocutoffnode, summaryfold)

      # Report from optional automatic cutoff
      self.cutoffText(self.xmlnode, summaryfold)
      
    # ctruncate
    if self.ctruncatereports != None:
      nctruncates = len(self.ctruncatereports)
      for ctruncatereport in self.ctruncatereports:
        if nctruncates > 1:
          text = 'Key warnings for dataset'+ctruncatereport.getDatasetid()
          summaryfold.addText(text=text,
                              style='font-weight:bold; font-size:120%;')
          summaryfold.append('<br/>')
        ctruncatekey = summaryfold.addDiv(\
          style="margin-left: 40px;")

        ctruncatereport.keyText(ctruncatekey,
                                phasertncs=self.phasertNCS,
                                phasertwin=self.phaserTwin)

    freerxml = self.xmlnode.findall('FREERFLAG')
    if len(freerxml) > 0:
      freerxml = freerxml[0]
      self.addFreerReports(freerxml, summaryfold)

    if self.phaserErrorMessage is not None:
      warnDiv = summaryfold.addDiv(style="color:darkorange;font-size:120%;")
      warnDiv.append("WARNING  Phaser was not run")
      warnDiv.append(self.phaserErrorMessage)

    if not errorspresent:  # or already done
      self.displayLogFiles(summaryfold)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Main summary
    overallsummaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    overallfold = overallsummaryDiv.addFold(label='Overall summary',
                                            brief='Summary',initiallyOpen=True)
    # Put the key messages about spacegroup and resolution at the very top
    headlineDiv = overallfold.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:5px;")
    leftDiv = headlineDiv.addDiv(style="width:49%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;border:0px;")
    if not self.merged:  
      # for merged files, not usually much intereseting in Pointless output
      if havePointlessReport:
        self.addPointlessSummary(leftDiv)

        # Pointless element table in:
        #    (1) if 1 dataset, in leftDiv with Pointless summary
        #    (2) if >1 dataset in new Div on right of Pointless summary
        if self.numberofdatasets == 1:
          nextDiv = leftDiv.addDiv(style="border: 1px solid black; margin: 1px; padding:1px;")
        else:
          nextDiv = headlineDiv.addDiv(style="width:49%;float:left;margin:1px;text-align:center; \
          padding:3px;border:1px solid black;")

          self.pointlessreport.ElementScoresTable(nextDiv)

    if haveAimlessReport:
      # Aimless summary
      if self.numberofdatasets == 1:
        statsDiv = headlineDiv.addDiv(style="width:46%;float:left;text-align:center;margin:0px; padding:10px;border:1px solid black;")
      else:
        statsDiv = overallfold.addDiv(style="border:1px solid black; clear:both; margin:5px; padding:10px;")

      table_info = None
      if havePhaserReport:
        # getinformation content for the Table 1, in resolution bins
        # List for each dataset, list of [Overall, Inner, Outer]
        table_info = [self.infocontent()]
        #print("table_info", table_info)

      if self.merged:  
        self.addAimlessSummaryMerged(statsDiv, table_info)
      else:
        self.aimlessreport.addAimlessSummary(statsDiv, table_info)

      if self.numberofdatasets > 1:
        interdatasetDiv = overallfold.addDiv(style="border-width: 1px; border-color: black; \
            clear:both; margin:0px; padding:0px;")
        aHeaderDiv = interdatasetDiv.addDiv(style=\
                         "clear:both;font-weight:bold; font-size:130%;text-align:center;")
        aHeaderDiv.append('<br/>Comparison of differences between datasets')
        self.aimlessreport.addInterDatasetGraphs(interdatasetDiv)

        self.aimlessreport.addScaleBfacReport(interdatasetDiv)

      # Main graphs as function of resolution and batch
      self.aimlessreport.importantGraphs(headlineDiv, self.merged)

      nextDiv = overallfold.addDiv(style="clear:both;")
      nextleftDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
      nextrightDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")

      if self.ctruncatereports != None:
        self.addTruncateWilson(nextleftDiv)

      # Phaser information plot
      if self.phaserxmlnodelist is not None:
        self.addPhaserInfoGraphs(nextrightDiv)

      # Phaser analysis
      if self.phaserxmlnodelist is not None:
        phaserDiv = parent.addDiv(style='width:100%; clear:both; margin:0px; padding:0px;')
        phaserfold = phaserDiv.addFold(label='Analyses of twinning, tNCS and anisotropy from Phaser',
                                            brief='Phaser',initiallyOpen=True)
        headerDiv = phaserfold.addDiv(style='width:100%;text-align:center;border-width: 1px; border-color: solid black;')
        self.phaseranalysis(headerDiv)

      if self.ctruncatereports != None:
        toptruncateDiv = parent.addDiv(\
          style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        truncatefold = toptruncateDiv.addFold(label='Analysis of twinning from ctruncate',
                                            brief='Truncate',initiallyOpen=True)

        ctruncateDiv = truncatefold.addDiv(style='width:100%;text-align:center;')
        headerDiv = ctruncateDiv.addDiv(style='width:100%;text-align:center;font-weight: bold; font-size:130%;')
##        headerDiv = ctruncateDiv.addDiv(style='width:100%;text-align:center;font-weight: bold; font-size:130%;')
        headerDiv.addText(text='Analysis for twinning from ctruncate, more details in Istats section')


        nctruncates = len(self.ctruncatereports)
        for ctruncatereport in self.ctruncatereports:
          twinned = ctruncatereport.addWarningTwin(headerDiv,check=True)
          ctruncateDiv.append('<br/>')
          if nctruncates > 1:
            text = 'Statistics for dataset'+ctruncatereport.getDatasetid()
            ctruncateDiv.addText(text=text)
            ctruncateDiv.append('<br/>')
          if twinned:
            ctruncateDiv.addText(text='This dataset is probably twinned',
                          style='color:red;')
          else:
            ctruncateDiv.addText(text='This dataset is probably NOT twinned',
                              style='font-style:italic; color:green;')

          twinDiv = ctruncateDiv.addDiv()
          twinleftDiv = twinDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
          twinleftDiv.addText(text='L-test for twinning',style="font-weight:bold;x font-size:130%;")
          twinleftDiv.append('<br/>')
          ctruncatereport.CtruncateLtest(twinleftDiv)

          twinrightDiv = twinDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
          ctruncatereport.acentricMoments(twinrightDiv)


    if havePointlessReport or haveAimlessReport:
      details1Div = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      details1Div.append('<br/>')
      if havePointlessReport:
        self.pointlessreport.Details(details1Div)

    if haveAimlessReport:
      nextgraphsDiv = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      nextfold = nextgraphsDiv.addFold(label='Other merging graphs', brief='Merging', initiallyOpen=True)
      self.aimlessreport.moreGraphs(nextfold)

      if not self.merged:
        sdDiv = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        sdfold =sdDiv.addFold(label='SD analysis', brief='SDanalysis', initiallyOpen=True)
        self.aimlessreport.sdAnalysis(sdfold)
      
      details2Div = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      fold = details2Div.addFold(label="Details of merging", brief='Details')
      #  Resolution estimates
      leftDiv = fold.addDiv(style="width:50%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;")
      aHeaderDiv = leftDiv.addDiv(style="clear:both;font-weight:bold; font-size:130%;")
      aHeaderDiv.append('Resolution estimates')      
      self.aimlessreport.resolutionEstimateTable(leftDiv)

      # Other stuff
      rightDiv = fold.addDiv(style="width:50%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;")
      self.aimlessreport.otherStatistics(rightDiv)
      
      fold.append('<br/>')
      crossDiv = fold.addDiv(style="width:100%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;")
      crossHDiv = crossDiv.addDiv(style="clear:both;font-weight:bold; font-size:130%;text-align:left")

      self.aimlessreport.interRunTable(crossDiv)
      
      truncateDiv = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      #truncateDiv.append('<br/>')
      fold = truncateDiv.addFold(label="Intensity statistics: twinning tNCS etc",brief='Istats')

      if self.ctruncatereports != None:
        self.addCtruncateReports(fold)

      #The following is to return the html to doing one item per line
      parent.addDiv(style="width:95%; clear:both;")

  # - - - - - - - - - - - - - - - - - - - - - - - 
  def addAimlessSummaryMerged(self, parent=None, extraitems=None):
    aHeaderDiv = parent.addDiv(
      style="clear:both;font-weight:bold; font-size:130%;margin:0px;padding:0px;")
    aHeaderDiv.append('Data internal consistency statistics')

    #print("addAimlessSummaryMerged", extraitems)

    self.aimlessreport.ResultTableMerged(parent, extraitems)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  def infocontent(self):
    # getinformation content for the Table 1, in resolution bins
    # Overall, Inner, Outer
    table_info = ['Information content bits/reflection',
                  '']
    # first we need the resolution bins used by Aimless
    #  [low,high] for each bin for each dataset
    limits = self.aimlessreport.getResolutionBins()
    #print("limits",limits)

    for dsetidx in range(self.numberofdatasets):  # Loop datasets
      phaseridx = self.datasetIndex.phaserIndex(dsetidx)  # index into phaser lists
      dsetlimits = limits[dsetidx]
      dsetinfo = []  # Overall, Inner, Outer
      
      for binidx in range(3):
        resrange = [float(r) for r in dsetlimits[binidx]]
        value = self.phaserreports[phaseridx].getInformationValues(resrange)
        dsetinfo.append("{:5.2f}".format(value))

      #print("dsetinfo",dsetinfo)

      table_info.append(dsetinfo)
        
    return table_info

  # - - - - - - - - - - - - - - - - - - - - - - - 
  def cutoffText(self, xmlroot, parent):
    # Display message about AUTOCUTOFF option
    # Autocutoff block is in self.autocutoffnode

    cutoffDiv = parent.addDiv(style='color:Blue;font-size:105%;margin-left: 20px;')

    # Extract CChalf threshold
    threshold = -1.0
    estimates = xmlroot.findall('Result/Dataset/ResolutionLimitEstimate[@type="CChalf"]')
    for estimate in estimates:
      direction = estimate.findall('Direction')[0].text
      if direction =='Overall':
        threshold = float(estimate.findall('Threshold')[0].text)

    s = ""
    if self.autocutoffnode is not None:
      acdone = self.autocutoffnode.findall('Cutoff')[0].text

      cutoffdatasets = self.autocutoffnode.findall('Dataset')
      ndatasets = len(cutoffdatasets)

      inputresolution = 10000.
      for cdataset in cutoffdatasets:
        name = cdataset.attrib['name']
        res = float(cdataset.findall('InputResolution')[0].text)
        inputresolution = min(inputresolution, res)
      if inputresolution > 1000.:
        inputresolution = None
        
      scut = " CC(1/2) > "+str(threshold)  # CC(1/2)
      if len(self.autocutoffnode.findall('InformationThreshold'))>0:
        ICthreshold = self.autocutoffnode.findall('InformationThreshold')[0].text
        scut = " information content > "+ICthreshold

      achires = ''
      s += \
          "Automatic resolution cutoff was requested"
      if acdone == "Yes":
        # cutoff was applied
        if len(self.autocutoffnode.findall('ResolutionCutoff'))>0:
          achires = self.autocutoffnode.findall('ResolutionCutoff')[0].text
        s += " and applied, based on 1st Aimless run "
        s += "<br/>Resolution was cut automatically to "+achires+"&#197; from"+\
             scut
        if ndatasets > 1:
          s += ", from highest input resolution "+str(inputresolution)+"&#197;"
        else:
          s += ", from input resolution "+str(inputresolution)+"&#197;"
      else:
        #  not done
        if ndatasets > 1:
          s += ", <br/> but was not applied as at least one dataset goes to the maximum input resolution"
        else:
          s += ", <br/> but was not applied as the data go to the maximum input resolution"

        if inputresolution is not None:
          s += " of "+str(inputresolution)+"&#197;"

      cutoffDiv.append(s)

      elimit = cutoffDiv.addDiv(style='color:Blue;font-size:95%;margin-left: 40px;')
      if ndatasets == 1:
        name = cutoffdatasets[0].attrib['name']
        dtsreslimit = cutoffdatasets[0].findall('ResolutionEstimate')[0].text
        elimit.append("Resolution limit estimate for dataset ["+name+
                      "] is "+dtsreslimit+"&#197;, "+scut)
        
      elif ndatasets > 1:
        elimit.append("Resolution limit estimates for each dataset, "+
                      "with "+scut)
        table = elimit.addTable(transpose=True, class_="center")
        labels=["Resolution limit &#197;", ""]
        table.addData(title="Dataset", data=labels)
      
        for cdataset in cutoffdatasets:
          name = cdataset.attrib['name']
          dtsreslimit = cdataset.findall('ResolutionEstimate')[0].text

          message = " "
          if len(cdataset.findall('Message'))>0:
            message = cdataset.findall('Message')[0].text

          table.addData(title=name, data=[dtsreslimit, message])

    else:
      s += \
          "No automatic resolution cutoff was requested"
      cutoffDiv.append(s)
    
  # - - - - - - - - - - - - - - - - - - - - - - - 
  def displayLogFiles(self, parent):
    numphaser = self.numberofdatasets   # Phaser run for each dataset
    if self.twoaimless:
      # Aimless was run twice
      naimless = "{:1d}".format(3+numphaser)
      nctruncate = "{:1d}".format(4+2*numphaser)
      aimlesslog = 'job_'+naimless+'/log.txt'
      ctruncatelog = 'job_'+nctruncate+'/log.txt'
    else:
      # Aimless was run once
      aimlesslog = 'job_2/log.txt'
      ctruncatelog = 'job_4/log.txt'
      
    displayFileList(self.fileroot, parent,
                    [['job_1/log.txt', ' Show Pointless logfile'],
                     [aimlesslog, ' Show Aimless logfile'],
                     [ctruncatelog, ' Show Ctruncate logfile']], projectid=self.projectId, jobNumber=self.jobNumber)

  # - - - - - - - - - - - - - - - - - - - - - - - 
  def mmcifreport(self, importlogxml, mmcifDiv):
    #  information about the mmCIF file
    text = []

    if len(importlogxml.findall("mmcifblockinfo"))>0:
      info = importlogxml.findall("mmcifblockinfo")[0].text
      info = info.splitlines()
      bname = importlogxml.findall("mmcifblock")[0].text
      text.append("Data types in cifblock "+bname+": "+info[0])
      text.append(info[1]+": mmCIF data are assumed to be already scaled")
        
    if len(importlogxml.findall("mmcifblockdetails"))>0:
      details = importlogxml.findall("mmcifblockdetails")[0].text
      text.append(details)

    if len(importlogxml.findall("mmcifblockcolumns"))>0:
      cifcolumns = importlogxml.findall("mmcifblockcolumns")[0].text
      text.append( "mmCIF columns used: "+cifcolumns)

    mmcifDiv.append("<br/>")
    mmcifDiv.addText(text='Information about the imported  mmCIF file :',
                        style='font-weight:bold; font-size:110%; color:blue;')

    #mmcifDiv.append("<br/>")
    for line in text:
      mmcifDiv.addText(text=line)
      #mmcifDiv.append("<br/>")

  # - - - - - - - - - - - - - - - - - - - - - - - 
  def pointlessReport(self, parent=None):
    # just for Pointless

    if parent is None: parent = self
##    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,style="width:1280px;overflow:auto;",**kw)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Instantiate pointless report from which to cherry pick
    # The nooutput flag makes sure thay don't do anything
    # First check if any parts of the report exist

    errorspresent = False
    pointlessFatalErrors = False
    
    #  0) PIPELINE_ERROR  fatal error somewhere
    if len(self.xmlnode.findall('PIPELINE_ERROR'))>0:
      errorspresent = True
      errormessage =  self.xmlnode.findall('PIPELINE_ERROR')[0].text
      print("errormessage", errormessage)
      errorDiv = parent.addDiv(
        style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
      errorDiv.addText(text='FATAL ERROR',
                       style='font-weight:bold; font-size:150%; color:red;')
      errorDiv.append(errormessage)
      
    #  1) POINTLESS
    pointlessxml = self.xmlnode.findall("POINTLESS")[0]
    havePointlessReport = (pointlessxml != None) and (len(pointlessxml) > 0)

    # Empty XML file
    if not havePointlessReport:
      parent.addText(text='No Report data available, probably due to failure of job',
                     style='color:red;')
      return

    # Pointless report exists
    self.pointlessreport = \
           pointless_report.pointless_report(xmlnode=pointlessxml, jobStatus='nooutput')
    self.pointlessreport.setFileRoot(self.fileroot) # pass on fileroot
    #print("self.pointlessreport", self.pointlessreport)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Process fatal error messages, if any
    if havePointlessReport and self.pointlessreport.isFatalError():
      errorspresent = True
      pointlessFatalErrors = True

    if errorspresent:
      if pointlessFatalErrors:
        errorDiv = parent.addDiv(
          style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
        errorDiv.append('<br/>')
        errorDiv.addText(text='Errors from Pointless', style='font-size:120%;color:red;')
        self.pointlessreport.Errors(errorDiv)   #  report any fatal errors
        
    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Summary of summaries
    summaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    if havePointlessReport:
      fail = self.pointlessreport.Errors(summaryDiv)   #  report any fatal errors
      if fail:
        return

      summaryDiv.addText(text="Report from Pointless while Aimless is running",
                         style='font-weight:bold; font-size:150%; color:red;')

      self.pointlessreport.keyText(summaryDiv)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . . 
    # Main summary
    overallsummaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

    # Put the key messages about spacegroup and resolution at the very top
    headlineDiv = overallsummaryDiv.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:5px;")
    leftDiv = headlineDiv.addDiv(style="width:49%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;border:0px;")
    self.addPointlessSummary(leftDiv)
    rightDiv = headlineDiv.addDiv(style="width:50%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;border:1px solid black;")

    if havePointlessReport:
      self.pointlessreport.ElementScoresTable(rightDiv)
    
      details1Div = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      details1Div.append('<br/>')
      if havePointlessReport:
        self.pointlessreport.Details(details1Div, usefold=False, elementscores=False,
                                     all=True, open1=True)
      
      #The following is to return the html to doing one item per line
      parent.addDiv(style="width:95%; clear:both;")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addPointlessSummary(self, parent=None):
    spacegroupDiv = parent.addDiv(
      style="text-align:center;margin:0px auto; padding:3px; \
            line-height:100%;border:1px solid black;")
    spacegroupDiv.addText(text='Space group determination',style="font-size:130%;")
    self.pointlessreport.MainMessages(spacegroupDiv)
    self.pointlessreport.Errors(spacegroupDiv)   #  report any fatal errors

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addCtruncateReports(self, parent=None):
    ndatasets = len(self.ctruncatereports)
    for ctruncatereport in self.ctruncatereports:
      ctruncate_report.ctruncate_report(xmlnode=ctruncatereport, jobStatus='nooutput')
      ctruncatereport.addCtruncateReport(parent, ndatasets)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def addTruncateWilson(self, parent=None):
    ndts = len(self.ctruncatexmlnodelist)
    title = "Wilson plot B: "
    for ctruncatenode in self.ctruncatexmlnodelist:
      if ndts > 1:  # multiple datasets
        dts = ctruncatenode.findall('ReflectionFile/CrystalDatasetId')[0].text
        title += dts.split('/')[2]+" "
      title += ctruncatenode.findall('DataStatistics/WilsonB')[0].text+" "
    WilsonDiv = parent.addDiv(
      style="text-align:center;margin:0px auto; padding:1px;")
    WilsonDiv.addText(text=title,style="font-weight:bold; font-size:125%;")
    WilsonDiv.append('<br/>')
    #print("addTruncateWilson",self.ctruncatereports[0].xmlnode)
    for ctreport in self.ctruncatereports:
      ctreport.CtruncateWilson(WilsonDiv)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addFreerReports(self, freerxml, parent=None):
    """
    Cases:
    Status 0, Validity True:  all OK
    Status 1, Validity True:  failed due to symmetry mismatch
    Status 1, Validity False: failed due to cell mismatch
    """
    if freerxml != None:
      freeStatus = int(freerxml.findall('status')[0].text)
      validity = 'True'
      freeRmessage = []
      OK = True
      if len(freerxml.findall('ObsFreeCellComparison'))>0:
        validity = freerxml.findall('ObsFreeCellComparison/validity')[0].text

      printDetails = False

      if freeStatus== 0:
        # All OK, success, but possibly due to overriding the test
        freefraction = '0.05'   # default fraction in freerflag program
        if len(freerxml.findall('ObsFreeCellComparison'))>0:
          # Extending FreeR set
          freeRmessage.append(colourText('<p>A free-R set has been copied (and extended)</p>', 'green'))
        else:
          if len(freerxml.findall('FreeRfraction'))>0:
            freefraction = freerxml.findall('FreeRfraction')[0].text
          if len(freerxml.findall('FREERFLAGINFO/Fraction'))>0:
            freefraction = freerxml.findall('FREERFLAGINFO/Fraction')[0].text
          freeRmessage.append(colourText('<p>A free-R set has been created, fraction of the data = '+\
                                         freefraction+'</p>', 'green'))
        if validity == 'False':
          # cell check was overridden
          printDetails = True
          freeRmessage.append(colourText(
            '<p><b>WARNING: the discrepancy between the unit cell of the extended FreeR set and the observed data has been overridden</b></p>',
            'red'))
          freeRmessage.append(colourText(
            '<p><b>Did you mean to accept this?</b></p>',
            'red'))
          OK = False
      else:
        # Not OK
        OK = False
        # basic message
        message = '<p>'+\
           'WARNING: the FreeR set has not been created, because the input FreeR set is incompatible with the new data'+\
           '<br/>An input FreeR set for copying or extending must match the current data in cell and Laue group'+\
           '<br/>You should create an appropriate FreeR set in a separate task, or if you are confident it is OK, activate the cell-test override in "Additional options"'
        freeRmessage.append(colourText(message+'</p>', 'red'))
        printDetails = True

      if printDetails and len(freerxml.findall('ObsFreeCellComparison'))>0:
        # Only if there is a problem
        #print('ObsFreeCellComparison')

        # some XML names were changed in Nov 2022, allow both PRE
        if len(freerxml.findall('ObsFreeCellComparison/cellObserved'))>0:
          # old style
          sg1 = 'SGnameObserved'
          cl1 = 'cellObserved'
          sg2 = 'SGnameFreeR'
          cl2 = 'cellFreeR'
        else:
          # New style
          sg1 = 'sgname1'
          cl1 = 'cell1'
          sg2 = 'sgname2'
          cl2 = 'cell2'

        sgname1 = freerxml.findall('ObsFreeCellComparison/'+sg1)[0].text
        sgname2 = freerxml.findall('ObsFreeCellComparison/'+sg2)[0].text
        cell1 = freerxml.findall('ObsFreeCellComparison/'+cl1)[0].text
        cell2 = freerxml.findall('ObsFreeCellComparison/'+cl2)[0].text
        
        obs = "Observed data: SG "+sgname1+" Cell: "+cell1
        free = "FreeR data: SG "+sgname2+" Cell: "+cell2
        diff = "Cell difference: "+\
               freerxml.findall('ObsFreeCellComparison/CellDifference')[0].text+"&#197;"+\
               ", maximum acceptable resolution for free-R extension: "+\
               freerxml.findall('ObsFreeCellComparison/MaxAcceptableResolution')[0].text+"&#197;"
        if validity == 'True':
          message = colourText('The datasets belong to incompatible Laue groups',
                                     'red', fontstyle='italic', fontsize="110%")
        else:
          message = colourText('The datasets have incompatible unit cells',
                                     'red', fontstyle='italic', fontsize="110%")
          
        freeRmessage.append(message)
        freeRmessage.append(colourText('<br/>'+obs+'<br/>'+free+'<br/>'+diff, 'red'))

      if OK:
        if len(freerxml.findall('FREERFLAGINFO/FreerCutResolution'))>0:
          cutres = freerxml.findall('FREERFLAGINFO/FreerCutResolution')[0].text
          message = \
             "The resolution of the FreeR set was cut to match the data, "+cutres+" A"
          freeRmessage.append(colourText(message, 'green'))

          
      if len(freeRmessage) > 0:
        place = parent
        if not OK:
          place = parent.addDiv(
            style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
  
        for message in freeRmessage:
          #print('line: ', message)
          place.append(message)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addtoEstimatesofResolution(self):
    # add autocutoff and some Phaser things
    #  to self.estimatesofresolution (list for each dataset)
    # Autocutoff XML in self.autocutoffnode
    # Phaser XML in self.phaseranalysisxml
    print("addtoEstimatesofResolution")

    ndatasets = len(self.estimatesofresolution)
    # Index self.estimatesofresolution entries by dataset
    idx = 0
    eorindex = {}
    for idx in range(ndatasets):
      dsetname = self.estimatesofresolution[idx].datasetname
      eorindex[dsetname] = idx
    
    # Extract autocutoff stuff, if present
    originalresolution = None
    ICthreshold = None
    resolutioncutoff = None
    cutoffdone = False
    cutoffinfo = {}
    # this block will only be present if autocutoff was requested
    if self.autocutoffnode is not None:
      # Autocutoff stuff from phaser_analysis
      if len(self.autocutoffnode.findall('InformationThreshold'))>0:
        cutoffinfo['InformationThreshold'] = \
                      self.autocutoffnode.findall('InformationThreshold')[0].text
        cutoffinfo['ResolutionCutoff'] = \
                          self.autocutoffnode.findall('ResolutionCutoff')[0].text
      cutoffinfo['Cutoff'] = self.autocutoffnode.findall('Cutoff')[0].text  #  Yes/No

      cutoffdatasetnodes = self.autocutoffnode.findall('Dataset')
      if len(cutoffdatasetnodes) > 0:
        if len(cutoffdatasetnodes) != ndatasets:
          print("!!Help, cutoffdatasetnodes")
          return
        for cdataset in cutoffdatasetnodes:
          pxdname = cdataset.attrib['name']
          cutoffinfo['pxdname'] = pxdname
          cutoffinfo['InputResolution'] = cdataset.findall('InputResolution')[0].text
          cutoffinfo['ResolutionEstimate'] = \
                 cdataset.findall('ResolutionEstimate')[0].text
          cutoffinfo['Message'] = None
          if len(cdataset.findall('Message'))>0:
            cutoffinfo['Message'] = cdataset.findall('Message')[0].text

          # Add autocutoff stuff into EstimatesofResolution
          idx = eorindex[pxdname]
          self.estimatesofresolution[idx].addAutocutoff(cutoffinfo)

    #  Phaser block, just the resolution bits, for each dataset
    if self.phaseranalysisxml is not None:
      phaserdatasetnodes = self.phaseranalysisxml.findall('PHASER_ANALYSIS')
      if len(phaserdatasetnodes) != ndatasets:
        print("!!Help, phaserdatasetnodes")
        return
      
      for phaserdatasetnode in phaserdatasetnodes:
        pxdname = phaserdatasetnode.attrib['name']
        phaserinfo = {}
        phaserinfo['pxdname'] = pxdname
        estimatestuff = phaserdatasetnode.findall('ResolutionEstimate')[0]
        phaserinfo['type'] = estimatestuff.attrib['type']
        phaserinfo['Threshold'] = estimatestuff.findall('Threshold')[0].text
        phaserinfo['ResolutionLimitEstimate'] = \
                            estimatestuff.findall('ResolutionLimitEstimate')[0].text
        phaserinfo['Message'] = None
        phaserinfo['Maxres'] = False
        if len(estimatestuff.findall('Message'))>0:
          phaserinfo['Message'] = estimatestuff.findall('Message')[0].text
        if phaserinfo['Message'] == "== maximum resolution":
          phaserinfo['Maxres'] = True
          # data go to edge, replace ResolutionLimitEstimate by actual maximum
          phaserreso = \
            phaserdatasetnode.findall('Analysis/Resolution/InputResolution')[0].text
          phaserinfo['ResolutionLimitEstimate'] = phaserreso
        # Add autocutoff stuff into EstimatesofResolution
        idx = eorindex[pxdname]
        self.estimatesofresolution[idx].addPhaserInfo(phaserinfo)
  #     
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addPhaserInfoGraphs(self, parent=None):
    # <CCP4Table groupID="Graph" id="informationcontent" title=" Information content of data">

    icDiv = parent.addDiv(
      style="text-align:center;margin:3px auto; padding:1px; ;")
    icDiv.addText(text="Information content by resolution",
                  style="font-weight:bold; font-size:125%;")
    icDiv.append('<br/>')
    graphgroup = icDiv.addFlotGraphGroup(style="width:300px;  height:270px;border:1px solid black:padding:3px")
    for dsetidx in range(self.numberofdatasets):  # Loop datasets
      phaseridx = self.datasetIndex.phaserIndex(dsetidx)  # index into phaser lists
      thisgraph = self.phaserxmlnodelist[phaseridx].findall('CCP4Table[@id="informationcontent"]')[0]
      title = thisgraph.get("title")
      if self.numberofdatasets > 1:
        title += " "+self.datasetIndex.getDname(dsetidx)
        titlenode = thisgraph.find("plot").find("title")
        titlenode.text = title

      graph = graphgroup.addFlotGraph(xmlnode=thisgraph, title=title)
      graph = graph.addPimpleData(xmlnode=thisgraph)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def phaseranalysis(self, parent):
    # Analyses from Phaser

    for dsetidx in range(self.numberofdatasets):  # Loop datasets
      phaseridx = self.datasetIndex.phaserIndex(dsetidx)  # index into phaser lists
      dname = self.datasetIndex.getDname(dsetidx)
      phaserdtsDiv = parent.addDiv(style='width:100%;text-align:center;clear:both;')
      self.phaserreports[phaseridx].phaserAnalysisReport(phaserdtsDiv, False)


############################################################################
class datasetIndex:
  # Store dataset indices for Ctruncate and Phaser_analysis XML based on
  # the order in Aimless XML
  #
  # Multiple datasets may be in a different order when processed separately by
  #  ctruncate and phaser_analysis
  def __init__(self, aimlessXML):
    self.dnames = []
    self.dnameindex = {}
    self.ctruncateidx = [0]
    self.phaseridx = [0]
    
    #  Initialise from AIMLESS XML element
    datasets = aimlessXML.findall('Result/Dataset')
    # Index Aimless datasets by dname, extracted from pxdname
    i = 0
    for dataset in datasets:
      pxdname = dataset.attrib['name']
      dname = pxdname.split('/')[2]
      #print(dname)
      self.dnames.append(dname)
      self.dnameindex[dname] = i
      i += 1

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def getDname(self, dsetidx):
    return self.dnames[dsetidx]

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def indexCtruncate(self, ctruncatexmlnodelist):
    # index ctruncate XML blocks from nodelist
    # ctruncatexmlnodelist should contain list of CTRUNCATE blocks
    # print("indexCtruncate", ctruncatexmlnodelist)
    if ctruncatexmlnodelist is None: return
    
    self.ctruncateidx = [None]*len(self.dnames)

    i = 0
    for ctrblock in ctruncatexmlnodelist:
      xdname = ctrblock.findall('ReflectionFile/CrystalDatasetId')[0].text
      dname = xdname.split('/')[2]
      self.ctruncateidx[self.dnameindex[dname]] = i
      i += 1
    #print("ctruncateidx", self.ctruncateidx)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def indexPhaser(self, phaserxmlnodelist):
    # index phaser_analysis XML blocks from nodelist
    #  phaserxmlnodelist should contain list of PHASER_ANALYSIS blocks
    self.phaseridx = [None]*len(self.dnames)
    #print("phaserxmlnodelist", phaserxmlnodelist)
    i = 0
    for phsblock in phaserxmlnodelist:
      pxdname = phsblock.attrib['name']
      dname = pxdname.split('/')[2]
      self.phaseridx[self.dnameindex[dname]] = i
      i += 1
    #print("self.phaseridx", self.phaseridx)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def phaserIndex(self, aimlessidx):
    # Returns index for Phaser_analysis given the Aimless index
    return self.phaseridx[aimlessidx]


############################################################################
if __name__ == "__main__":
  report = aimless_pipe_report(xmlFile = sys.argv[1],jobStatus="Finished" )
  tree= report.as_etree()
  #print(etree.tostring(tree,pretty_print=True))
  report.as_html_file(fileName='./test-pipeline.html')
  if len(report.errorReport())>0: print('ERRORS:',r.errorReport())
  
