import sys

from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import (
    SDcorrectionData,
    colourText,
    displayFile,
    html_linebreak,
    median,
    selectGraphs,
)
from ccp4i2.report import Report


class EstimatesofResolution:
  """ Storage of resolution estimates for one dataset """
  def __init__(self, message, datasetname,
               thresholdCC, thresholdIovsd, thresholdIovsd2=None):
    self.message = message
    self.datasetname = datasetname
    self.thresholdCC = thresholdCC
    self.thresholdIovsd = thresholdIovsd
    if thresholdIovsd2 is not None:
      self.thresholdIovsd2 = thresholdIovsd2
    self.directionlimits = []   # list of lists for each direction
    self.cutoffinfo = None
    self.phaserinfo = None
    
  # - - - - - - - - - - - - - - - - -
  def addLimits(self, direction, reslims,
                nodata, maxres):
    """ store:
    direction message
    resolution limit list: (CC(1/2), I//sigI)
    nodata: T/F if  insufficient data
    maxres: T/F if at maximum resolution
    """
    self.directionlimits.append([direction, reslims,
                nodata, maxres])

  # - - - - - - - - - - - - - - - - -
  def addAutocutoff(self, cutoffinfo):
    # Add dictionary of data from Autocutoff (default None)
    self.cutoffinfo = cutoffinfo
    
  # - - - - - - - - - - - - - - - - -
  def addPhaserInfo(self, phaserinfo):
    # Add dictionary of data from Phaser (default None)
    # just resolution things
    self.phaserinfo = phaserinfo
    
  # - - - - - - - - - - - - - - - - -
  def formatBrief(self, inputresolution):
    """ return (a) overall limits sr (b) list of anisotropic limits """

    # Cases:
    #  Autocutoff was requested in aimless_pipe
    #    cutoff was done and aimless/phaser_analysis were run again
    #    no cutoff was done as the data go to the "edge" anyway
    #  No autocutoff was requested
    #    data go to edge  (beyond estimate)
    #    estimate < input
    cutoff = False
    sr = ""
    if self.cutoffinfo is not None:
      # Autocutoff was requested
      originalresolution = self.cutoffinfo['InputResolution']
      thresholdIC = None
      if 'InformationThreshold' in self.cutoffinfo:
        thresholdIC = self.cutoffinfo['InformationThreshold']
      if self.cutoffinfo['Cutoff'] == "Yes":
        # cutoff was done
        originalresolution = self.cutoffinfo['InputResolution']
        sr += "Automatic resolution cutoff was requested and applied,"+\
              " original "+originalresolution+\
              "&#197;, cut to "+inputresolution+"&#197;"
        if thresholdIC is not None:
          sr += ", from information content > "+thresholdIC+"\n - "
        cutoff = True
      else:
        # cutoff was requested but not done
        sr += "Automatic resolution cutoff was requested but not applied, "+\
              "from "+originalresolution+"&#197;"
        if thresholdIC is not None:
          sr += " as all data had information content > "+\
                thresholdIC+"\n - "
    else:
      # No autocutoff
      sr += "Estimates of resolution: input resolution "+\
            inputresolution+"&#197;"

    # Estimates from final phaser_analysis, if done
    if self.phaserinfo is not None:
      thresholdIC = self.phaserinfo['Threshold']
      reslimitIC = self.phaserinfo['ResolutionLimitEstimate']
      maxresIC = self.phaserinfo['Maxres']

      # Analysis was done with phaser_analysis using information content
      # criterion, so report this
      if self.cutoffinfo is None:
        sr += "\n - "
      if maxresIC:
        sr += " beyond "
      sr += reslimitIC+"&#197; from information content > "+thresholdIC

    # Estimates from CC(1/2) and I/sigI
    # Overall resolution estimate, first in list
    if self.directionlimits[0][0] == "Overall":
      directionlimit = self.directionlimits[0]
      reslimit = directionlimit[1] # overall limits from CC(1/2), I/sigI
      nodata = directionlimit[2]   # T/F if insufficient data CC(1/2),(I/sigI_
      maxres = directionlimit[3]   # True if maximum resolutionfrom CC(1/2)
      if maxres[0]:
        sr += ", beyond "+reslimit[0]+"&#197; from CC(1/2) > "+\
              str(self.thresholdCC)
      elif nodata[0]:
        sr += ", estimate undetermined from CC(1/2), insufficient data"
      else:
        sr += reslimit[0]+"&#197;"+" from CC(1/2) > "+\
              str(self.thresholdCC)
      sr +=  ", "+reslimit[1]+"&#197; from I/&#963; >"+\
             self.thresholdIovsd
      sr = html_linebreak(sr, False)

    # Anisotropy
    if len(self.directionlimits) <= 1:
      # no aniosotropy (cubic)
      return sr, None
      
    saniso = "Anisotropic limits: "
    for directionlimit in self.directionlimits:
      label = directionlimit[0]
      if label != "Overall":  # skip overall value
        reslimit = directionlimit[1]  # limits from CC(1/2), I/sigI
        nodata = directionlimit[2]
        maxres = directionlimit[3]
        saniso += label
        if maxres[0]:
          saniso += reslimit[0]+"&#197;"+" CC(1/2),"+\
                    reslimit[1]+"&#197; I/&#963;"
        elif nodata[0]:
          saniso += "undetermined"
        else:
          saniso += reslimit[0]+"&#197; CC(1/2),"+\
                    reslimit[1]+"&#197; I/&#963;"
    

    return sr, saniso

# - - - - - - - - - - - - - - - - -
class aimless_report(Report):
  TASKNAME='aimless'

  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=None,**kw)

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    self.datasetsProcessed = False
    self.estimatesofresolution = None
    self.referencedata = None
    
    # 'nooutput' mode would be used by another report class that wanted
    # to use some method(s) from this class for its own report
    if jobStatus is not None and jobStatus.lower() == 'nooutput':
      return
  
    fail = self.Errors(self)

    if fail == False:
      self.justAimless(self)

  # - - - - - - - - - - - - - - - - -
  def setup(self):
    # just extract datasetinformation
    self.processDatasets()
    self.getResolutionEstimates(self.datasetresultnodes)

    # Get info about reference data file
    # If HKLREF defined, sets self.referencedata
   #   name: filename, scaletoreference: True (or absent)
    reflectionfiles = self.xmlnode.findall("ReflectionFile")
    self.referencedata = None
    
    reffile = False
    for rf in reflectionfiles:
      # print("ReflectionFile", rf.items())
      reffile = False
      fname = ''
      for itm in rf.items():
        if 'stream' in itm:
          if 'HKLREF' in itm[1]:
            reffile = True
        if 'name' in itm:
          fname = itm[1]
      if reffile:  # HKLREF found
        self.referencedata = {'name': fname}

    if reffile and len(self.xmlnode.findall("ScalingType"))>0:
        # Scaling relative to reference
        self.referencedata['scaletoreference'] = True

    # print("Reffile info:", self.referencedata)

  # - - - - - - - - - - - - - - - - -
  def getEstimatesofResolution(self):
    return self.estimatesofresolution
  # - - - - - - - - - - - - - - - - -
  def storeEstimatesofResolution(self, estimatesofresolution):
    self.estimatesofresolution = estimatesofresolution
  # - - - - - - - - - - - - - - - - -
  def justAimless(self, parent=None):
    ''' report for single Aimless step '''
    
    parent.addText(text='AIMLESS', style='font-size: 150%;')

    displayFile(self.fileroot, parent,
                ['job_2/log.txt', './log.txt'], 'Show log file', self.projectId, self.jobNumber)
    summaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    self.keyText(None, summaryDiv)
    
    mainDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; margin:0px; padding:10px 0px 10px 0px;")

    if self.numberofdatasets == 1:
      # Use grid layout for two-column display
      leftDiv, rightDiv = mainDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
      self.addAimlessSummary(rightDiv)
      leftDiv.append('Resolution estimates')
      self.resolutionEstimateTable(leftDiv)
      self.otherStatistics(leftDiv)

    else:  #  > 1 datasets
      tableDiv = mainDiv.addDiv(style="width:100%;text-align:center;margin:0px; padding:0px; border:1px solid black;")
      self.addAimlessSummary(tableDiv)
      # Use grid layout for resolution/other stats
      leftDiv, rightDiv = mainDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
      leftDiv.append('Resolution estimates')
      self.resolutionEstimateTable(leftDiv)
      self.otherStatistics(rightDiv)

      self.addInterDatasetGraphs(mainDiv)

    nextDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:10px 0px 10px 0px;")
    self.importantGraphs(nextDiv)

    gDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:10px 0px 10px 0px;")

    self.moreGraphs(gDiv)
    rDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:10px 0px 10px 0px;")

    self.ReflectionDataDetails(rDiv, infold=False)

    
    self.sdAnalysis(parent)

  # - - - - - - - - - - - - - - - - -
  def setFileRoot(self, fileroot):
    ''' set file root from superior script '''
    self.fileroot = fileroot
    
  # - - - - - - - - - - - - - - - - -
  def keyText(self,autocutoffnode, parent=None):

    if not self.datasetsProcessed:
      self.setup()
      
    if self.numberofdatasets == 1:
      parent.addText(text="Key statistics for "+self.datasetheader,
                     style='font-weight:bold; font-size:120%;')

    cells = self.xmlnode.findall("Result/Dataset/cell")

    s = ""
    jdts = 0
    innerCompleteness = []
    overallCompleteness = []

    if len(self.datasetresultnodes) == 0:
      return


    for datasetresultnode in self.datasetresultnodes:  # loop datasets
      if self.numberofdatasets > 1:
        s += "<span style='font-weight:bold; font-size:120%;'>Key statistics for "+datasetresultnode.get("name")+"</span>"
      
      s += '<div style="margin-left: 40px;">'

      if len(cells) > jdts:
        celltext = self.formatCell(cells[jdts], astext=True)
        s += "Unit cell: "+celltext
        wavelengths = self.getWavelength()   # may be None
        if wavelengths is not None:
          s += ", wavelength "+wavelengths[jdts]+"&#197;"
        s += "<br/>"
    
      datareso = datasetresultnode.findall("ResolutionHigh/Overall")[0].text

      if self.estimatesofresolution is not None:
        sr, saniso = \
            self.estimatesofresolution[jdts].formatBrief(datareso)

        s += "<span style='color:Blue;font-size:100%;'>"+sr+"</span><br/>"
        if saniso is not None:
          s += saniso+"<br/>"

      s += " Rmeas: overall "+datasetresultnode.findall("Rmeas/Overall")[0].text
      s += ", inner bin "+datasetresultnode.findall("Rmeas/Inner")[0].text
      s += "<br/>"
      s += "In outer bin: Mean(I/&#963;I)"+\
           datasetresultnode.findall("MeanIoverSD/Outer")[0].text
      s += "    CC(1/2) "+\
           datasetresultnode.findall("CChalf/Outer")[0].text
      s += "<br/>"

      if len(datasetresultnode.findall("MeanChiSq"))>0:
        s += "Overall filtered Mean(chi^2): "+\
             datasetresultnode.findall("MeanChiSq/Overall")[0].text+"<br/>"

      # Completeness in inner shell and overall
      innerCompleteness.append(datasetresultnode.findall("Completeness/Inner")[0].text)
      overallCompleteness.append(datasetresultnode.findall("Completeness/Overall")[0].text)
      
      if self.hasAnomalous:
        s += "Anomalous CC(1/2) in inner bin  "+\
           datasetresultnode.findall("AnomalousCChalf/Inner")[0].text
        if len(datasetresultnode.findall('AnomalousLimitEstimate'))>0:
          anomresolution = datasetresultnode.findall("AnomalousLimitEstimate/MaximumResolution")[0].text
          ##print "anomresolution", anomresolution
          anomthreshold = datasetresultnode.findall("AnomalousLimitEstimate/Threshold")[0].text
          if float(anomresolution) == 0.0:
            s += "<br/>No significant anomalous signal detected<br/>"
          else:
            s += "<br/>Significant anomalous signal extends to a resolution of"+\
                 anomresolution+"A (above CCanom threshold "+anomthreshold+")<br/>"
      else:
        s += "No anomalous data present<br/>"
        
      s += "</div>"
      jdts += 1

    parent.append(s)

    self.getOutlierControlData()

    # Outliers
    rejectnumberunique = '0'
    rejectnumberfriedel = '0'
    rejectnumberemax = '0'
    if len(self.xmlnode.findall("Outliers/RejectNumberUnique"))>0:
      rejectnumberunique = self.xmlnode.findall("Outliers/RejectNumberUnique")[0].text
    if len(self.xmlnode.findall("Outliers/RejectNumberFriedel"))>0:
      rejectnumberfriedel = self.xmlnode.findall("Outliers/RejectNumberFriedel")[0].text
    if len(self.xmlnode.findall("Outliers/RejectNumberEmax"))>0:
      rejectnumberemax = self.xmlnode.findall("Outliers/RejectNumberEmax")[0].text

    # may be different SDrej limits for multiple datasets
    anomrejsd = ''
    for oca in self.outlierControlAnom:
      anomrejsd += oca[0] + ","
    anomrejsd = anomrejsd[:-1]   # remove trailing comma

    #print("SOCEMX",self.outlierControlEmax, self.outlierControlMain)
    s = "Number of rejected outliers: "+rejectnumberunique
    if self.outlierControlMain is not None:
        s += " (> "+self.outlierControlMain[0]+"&#963;)"+\
        "; between Friedel pairs: "+rejectnumberfriedel+\
        " (> "+anomrejsd+"&#963;)"+\
        "; too large: "+rejectnumberemax+\
        " (E >"+self.outlierControlEmax[0]+")"
    s = '<span style="margin-left: 20px;">'+s+'</span>'
    parent.append(s)

    if int(rejectnumberemax) > 0:
      # Warn if there are any Emax rejects
      s = "Warning: "+rejectnumberemax+\
          " observations were rejected as too large "
      if self.outlierControlEmax is not None:
        s += "(E >"+\
             self.outlierControlEmax[0]+"),"+\
             " check details below (Merging)"
      outlierdiv = parent.addDiv(style='color: red; margin-left: 20px;font-size:115%;')
      outlierdiv.append(s)
    if len(self.xmlnode.findall("OnlyMerge"))>0:
      s = "\nNOTE: no scaling was done, just merging"
      s = html_linebreak(s, False)
      s = '<span style="color: blue">'+s+'</span>'
      parent.append(s)
      
    # Check overloads
    noverloads = 0
    acceptedoverloads = 0
    if len(self.xmlnode.findall("ObservationFlags/ProfileFittedOverloads/NumberFlagged"))>0:
      noverloads = self.xmlnode.findall("ObservationFlags/ProfileFittedOverloads/NumberFlagged")[0].text
      acceptedoverloads = self.xmlnode.findall("ObservationFlags/ProfileFittedOverloads/NumberAccepted")[0].text
    #print "Overloads", noverloads,acceptedoverloads,type(noverloads)
    minInnerCompleteness = float(min(innerCompleteness))
    minOverallCompleteness = float(min(overallCompleteness))
    #print " Min completeness", minInnerCompleteness, type(minInnerCompleteness)
    COMPLETENESS_WARNING_THRESHOLD = 0.95; # fraction of overall completeness, for warning
    completeness_warning = (minInnerCompleteness < minOverallCompleteness*COMPLETENESS_WARNING_THRESHOLD)

    if (int(noverloads) > 0) or completeness_warning:
      s1 = ""
      if self.numberofdatasets > 1:
        s1 = " for "+str(self.numberofdatasets)+" datasets"
      s2 = ""
      if int(noverloads) > 0:
        s2 = noverloads+" observations were overloaded, "

      s = "Note: "+s2+"inner shell completeness: "+\
           ",".join(innerCompleteness)+"%"+s1

      if int(acceptedoverloads) > 0:
        s += "<br/>"+acceptedoverloads+" profile-fitted overloads were accepted"
        if completeness_warning:
          s += "<br/>Warning: low inner shell completeness persists after accepting overloaded strong observations"
      else:
        if completeness_warning:
          s += "<br/>Warning: low inner shell completeness may indicate that overloaded strong observations have been lost"

      overloaddiv = parent.addDiv(style='color: red;margin-left: 40px;')
      overloaddiv.append(s)

    # Refinement against reference
    if self.referencedata is not None:
      if 'scaletoreference' in self.referencedata:
        reffilename = self.referencedata['name']
        s = 'Scaling was done relative to a reference dataset from file:\n'+\
            reffilename
        s = html_linebreak(s, False)
        s = '<span style="color: blue">'+s+'</span>'
        parent.append(s)

    #  SD correction things
    sdDiv = parent.addDiv(style="margin-left: 20px;")
    s = colourText("SD correction information:\n","Black",
                   style='font-size:120%')
    self.sdcpar = self.SDcorrectionParameters(self.xmlnode)
    sdcstatus, sdcmessage = self.sdcpar.status()
    if sdcstatus > 0:
      if sdcstatus == +1:
        colour = 'darkorange'
      elif sdcstatus == +2:
        colour = 'red'
      s += colourText(sdcmessage + ", see SD analysis panel below for more details\n",
                        colour)
      s += self.getISa()
    elif  sdcstatus == 0:
      s += self.getISa()
    elif  sdcstatus < 0:
      # No SD correction refinement
      s += "SD correction parameters were not refined\n"

    #self.sdAnalysis(parent)   ### TESTING, move to fold

    s = html_linebreak(s, False)
    sdDiv.append(s)

  # - - - - - - - - - - - - - - - - -
  def keyTextMerged(self,parent=None, phaserinfo=None):
    # cut down version for merged files, one dataset
    print("keyTextMerged")
    if not self.datasetsProcessed:
      self.setup()
    self.numberofdatasets = self.numberOfDatasets()

    parent.append("Key statistics for "+self.datasetheader)
    self.sdcpar = None  # ignore SD corrections (should be absent)

    # usually no CC(1/2) information for merged files
    s = ""
    jdts = 0
    for datasetresultnode in self.datasetresultnodes:  # loop datasets
      s += '<div style="margin-left: 40px;">'
      datareso = datasetresultnode.findall("ResolutionHigh/Overall")[0].text

      s += "Resolution of input data: "+ datareso + "&#197;"
      if phaserinfo is not None:
        s += ", resolution estimate "+\
             phaserinfo['ResolutionLimitEstimate']+"&#197;"+\
             " for infomation content > "+\
             phaserinfo['Threshold']+" bits/reflection"
      
      reslimitnodes = \
       datasetresultnode.findall("ResolutionLimitEstimate")
      for reslimitnode in reslimitnodes:
        if reslimitnode.findall("[@type]")[0].attrib["type"] == "I/sd":
          thresholdIsd = reslimitnode.findall("Threshold")[0].text
          break
      s += "<br/>Estimates from I/sigI > "+thresholdIsd+"<br/>"
      for reslimitnode in reslimitnodes:
        if reslimitnode.findall("[@type]")[0].attrib["type"] == "I/sd":
          reslimit = reslimitnode.findall("MaximumResolution")[0].text
          direction = reslimitnode.findall("Direction")[0].text
          if direction == 'Overall':
            s += " - "
          else:
            s += " - Along "
          s += direction+", resolution estimate "+reslimit+"&#197;"
          s += "<br/>"
          
      s += "In outer bin: Mean(I/sdI)"+\
           datasetresultnode.findall("MeanIoverSD/Outer")[0].text
      s += "<br/>"

      if self.hasAnomalous:
        anomslope = datasetresultnode.findall("AnomalousNPslope")[0].text
        s += "Anomalous Q-Q plot slope = " + anomslope
        if (float(anomslope) < 1.1):
          s += "<br/>No significant anomalous signal detected<br/>"
        else:
          s += "<br/>Possibly significant anomalous signal<br/>"
      else:
        s += "No anomalous data present<br/>"
        
      s += "</div>"
      jdts += 1
            
    parent.append(s)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addAimlessSummary(self, parent=None, extraitems=None):
      aHeaderDiv = parent.addDiv(
        style="clear:both;font-weight:bold; font-size:130%;margin:0px;padding:0px;")
      aHeaderDiv.append('Data internal consistency statistics')
      self.ResultTable(parent, extraitems)

  # - - - - - - - - - - - - - - - - -
  def importantGraphs(self,parent=None, merged=False):
    """ Main graphs by resolution and batch """

    # print("Aimless importantGraphs", merged)

    summaryGraphsDiv = parent.addDiv(
      style="width:100%; border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    summaryGraphsDiv.append('<br/>')
    if self.numberofdatasets > 1:
      titleDiv = summaryGraphsDiv.addDiv(
        style="width:100%; border-width: 1px; border-color: black; clear:both; text-align:center; margin:0px; padding:0px;")
      titleDiv.addText(text='Note that there are separate graphs for each dataset',
                              style='font-size:150%; color: orange;')

    # Put the key graphs that report on data quality and image-to-image quality at the top
    if not merged:
      # Use grid layout for resolution and batch graphs side by side
      byResolutionDiv, byBatchDiv = summaryGraphsDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
    else:
      # For merged data, just resolution graph
      byResolutionDiv = summaryGraphsDiv.addDiv(style="width:100%;text-align:center;margin:6px; padding:0px; ")

    byResolutionDiv.addText(text='Analysis as a function of resolution',style='font-size:130%;font-weight:bold;')
    byResolutionDiv.append('<br/>')

    if merged:
      self.ByResolutionGraphsMerged(byResolutionDiv)
    else:
      resmessage = "Plot of CC(1/2) vs. resolution may indicate a suitable resolution cutoff, and indicate presence of an anomalous signal"
      byResolutionDiv.addText(text=resmessage,style='font-size:100%;font-style:italic;')
      byResolutionDiv.append('<br/>')
      resmessage = "(but check anisotropy)"
      byResolutionDiv.addText(text=resmessage,style='font-size:100%;font-style:italic;')
      byResolutionDiv.append('<br/>')
      self.ByResolutionGraphs(byResolutionDiv)

    if not merged:
      byBatchDiv.addText(text='Analysis as a function of batch',style='font-size:130%;font-weight:bold;')
      byBatchDiv.append('<br/>')
      batchmessage = "Analyses against Batch may show radiation damage, and which parts of the data should be removed"
      byBatchDiv.addText(text=batchmessage,style='font-size:100%;font-style:italic;')
      byBatchDiv.append('<br/>')
      batchmessage = " (but consider completeness)"
      byBatchDiv.addText(text=batchmessage,style='font-size:100%;font-style:italic;')
      byBatchDiv.append('<br/>')
      self.ByBatchGraphs(byBatchDiv)

    if self.referencedata is not None:
      if not merged:
        # We may have a graph of statistics against reference vs batch
        referenceDiv = summaryGraphsDiv.addDiv(style="width:100%;text-align:center;margin:6px; padding:0px; ")
        referenceDiv.addText(text="Analysis of agreement with reference by batch",
                             style="font-weight:bold; font-size:125%;")
        self.referenceGraph(referenceDiv)
      else:
        referenceDiv = summaryGraphsDiv.addDiv(style="width:100%;text-align:center;margin:6px; padding:0px; ")

      # and reference vs resolution
      referenceDiv.addText(text="Analysis of agreement with reference by resolution",
                           style="font-weight:bold; font-size:125%;")
      refresolist = self.xmlnode.findall("CCP4Table[@id='Graph-DatasetRefStatsVsReso']")
      if len(refresolist) > 0:
        self.referenceGraphReso(referenceDiv)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def moreGraphs(self, parent=None):
    ''' all graphs and scatterplots etc '''
    # Use grid layout for graphs and scatter plots side by side
    othergraphDiv, scatterplotDiv = parent.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
    othergraphDiv.addText(text='Other merging statistics graphs',style='font-weight:bold; font-size:130%;')
    othergraphDiv.append('<br/>')
    self.Graphs(othergraphDiv, select="All")

    scatterplotDiv.append('Scatter plots etc')
    scatterplotDiv.append('<br/>')
    self.aimlessScatterPlots(scatterplotDiv)

    displayFile(self.fileroot, scatterplotDiv,
                ['job_2/ROGUES.log', './ROGUES.log', './ROGUES'],
                'Show list of outliers', self.projectId, self.jobNumber)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addInterDatasetGraphs(self,parent=None):
    tableDiv = parent.addDiv(style="border-width: 1px; border-color: black; margin:0px; padding:0px;")

    # Intensity correlations in grid layout
    leftDiv, rightDiv = tableDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
    self.interdatasetIntensityTable(leftDiv)
    self.interdatasetIntensityGraph(rightDiv)

    # Anomalous/Dispersion section in grid layout
    leftDiv2, rightDiv2 = tableDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
    self.interdatasetAnomalousTable(leftDiv2)

    if self.numberofdatasets > 2:
      self.interdatasetDispersionTable(rightDiv2)
      # Anomalous and dispersive graphs in another row
      leftDiv3, rightDiv3 = tableDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
      self.interdatasetAnomalousGraph(leftDiv3)
      self.interdatasetDispersionGraph(rightDiv3)
    else:
      # if no dispersive differences, put anomalous graph on right
      self.interdatasetAnomalousGraph(rightDiv2)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def aimlessScatterPlots(self, parent=None):
    #    parent.append('Normal probability plot')
    self.NormalProbabilityPlot(parent)
    #    parent.append('Normal probability of anomalous differences')
    self.AnomalousNormalProbabilityPlot(parent)
    #    parent.append('DelAnom/RMS scatter plot')
    self.AnomalousScatterPlot(parent)
    #    parent.append('Outliers on detector (horizontal rotation axis)')
    self.RoguePlot(parent)
    
  # - - - - - - - - - - - - - - - - -
  def extractOutlierControlData(self, datablock):
    data = []
    for tag in ['SDrej','SDrej2','Reject2policy']:
      data.append(datablock.findall(tag)[0].text)
    return data
    
  # - - - - - - - - - - - - - - - - -
  def getOutlierControlData(self):
    # Get data from OutlierControl block
    self.outlierControlMain = None
    self.outlierControlAnom = []  # for each data set
    self.outlierControlEmax = None
    if len(self.xmlnode.findall("OutlierControl"))>0:
      outliercontrol = self.xmlnode.findall("OutlierControl")[0]
      mainblock = outliercontrol.findall("Main")[0]
      self.outlierControlMain = self.extractOutlierControlData(mainblock)
      anomblocks = outliercontrol.findall("Anom")  # for each dataset
      for anomblock in anomblocks:
        adata = self.extractOutlierControlData(anomblock)
        dataset = anomblock.get('dataset')
        adata.append(dataset)
        self.outlierControlAnom.append(adata)
      # Emax
      emaxblock = outliercontrol.findall("EmaxTest")[0]
      self.outlierControlEmax = []
      for tag in ['EmaxAcentric','EmaxCentric']:
        self.outlierControlEmax.append(emaxblock.findall(tag)[0].text)

  # - - - - - - - - - - - - - - - - -
  def Details(self,parent=None):
    """ Folded details """
    fold = parent.addFold(label="Details")
    self.ReflectionDataDetails(fold)

  # - - - - - - - - - - - - - - - - -
  def isFatalError(self):
    # True if there is a FatalErrorMessage
    if len(self.xmlnode.findall('FatalErrorMessage'))>0:
      return True
    if len(self.xmlnode.findall('ScaleModelError'))>0:
      return True
    return False

  # - - - - - - - - - - - - - - - - -
  def Errors(self,parent=None):
    fail = False
    if len(self.xmlnode.findall('FatalErrorMessage'))>0:
      message = self.xmlnode.findall('FatalErrorMessage')[0].text
      message = html_linebreak(message)
      parent.append(message)
      fail = True
    if len(self.xmlnode.findall('FatalErrorMessage2'))>0:
      message = self.xmlnode.findall('FatalErrorMessage2')[0].text
      message = html_linebreak(message)
      parent.append(message)
      fail = True

    # also add ScaleModelError
    if len(self.xmlnode.findall('ScaleModelError'))>0:
      message = self.xmlnode.findall('ScaleModelError')[0].text
      if len(self.xmlnode.findall('ScaleModelFail'))>0:
        message += self.xmlnode.findall('ScaleModelFail')[0].text
      message = html_linebreak(message)
      parent.append(message)
      fail = True

    return fail

  # - - - - - - - - - - - - - - - - -
  def Get3numbers(self, parent=None, tag=""):
    """ get a set of 3 numbers for Overall, Inner, Outer from tag """
    list = []
    if len(parent.findall(tag+"/Overall"))>0:
      list = [parent.findall(tag+"/Overall")[0].text]
      list.append(parent.findall(tag+"/Inner")[0].text)
      list.append(parent.findall(tag+"/Outer")[0].text)
    elif len(parent.findall(tag))>0:  # just one entry
      list = [parent.findall(tag)[0].text, "",""]
    return list

  # - - - - - - - - - - - - - - - - -
  def AddToTable(self, datasetresultnodes="",\
                 table="", label="", tag="", tip=""):
    data = []
    ndts = len(datasetresultnodes)
    count = 0
    found = False
    for datasetresultnode in datasetresultnodes:  # loop datasets
      if len(datasetresultnode.findall(tag))>0:
        num3 = self.Get3numbers(datasetresultnode, tag)
        data += num3
        count += 1
        if (count < ndts): data += " "
        found = True

    if found:
      table.addData(title=label, data=data, tip=tip)
    return found

  # - - - - - - - - - - - - - - - - -
  def numberOfDatasets(self):
    if not len(self.xmlnode.findall("ReflectionData/NumberDatasets"))>0:
      return 0
    return int(self.xmlnode.findall("ReflectionData/NumberDatasets")[0].text)

  # - - - - - - - - - - - - - - - - -
  def processDatasets(self):
    if self.datasetsProcessed: return
    self.datasetsProcessed = True
    self.onebatchperdataset = False

    numberofdatasets = self.numberOfDatasets()
    if numberofdatasets == 0: return
    # List of Results for each dataset
    datasetresultnodes = self.xmlnode.findall("Result/Dataset")
    self.hasAnomalous = False

    self.numberofdatasets = numberofdatasets
    self.datasetresultnodes = datasetresultnodes
    if len(datasetresultnodes) == 0:
      self.numberofdatasets = 0
      return
    
    self.datasetnames = []
    if numberofdatasets == 1:
      datasetname = datasetresultnodes[0].get("name")
      self.datasetnames.append(datasetname)
      self.datasetheader = "Dataset: "+datasetname
      if len(datasetresultnodes[0].findall("AnomalousCompleteness"))>0:
        anomcompleteness = datasetresultnodes[0].findall("AnomalousCompleteness/Overall")[0]
        if float(anomcompleteness.text) > 0.0:
          self.hasAnomalous = True
    else:
      #  Multiple datasets
      self.datasetheader = "Datasets: "
      self.datasetlabellist = []
      notfirst = False
      for datasetresultnode in datasetresultnodes:
        if len(datasetresultnode.findall("AnomalousCompleteness"))>0:
          anomcompleteness = datasetresultnode.findall("AnomalousCompleteness/Overall")[0]
          if float(anomcompleteness.text) > 0.0:
            self.hasAnomalous = True
        datasetname = datasetresultnode.get("name")
        self.datasetnames.append(datasetname)
        if notfirst:
          self.datasetheader += ", "
        notfirst = True
        self.datasetheader += datasetname
        self.datasetlabellist.append(datasetname)

        # Is the number of batches == number of datasets?
        nbatches = int(self.xmlnode.findall("ReflectionData/NumberBatches")[0].text)
        if nbatches == numberofdatasets:
          self.onebatchperdataset = True

  # - - - - - - - - - - - - - - - - -
  def getWavelength(self):   # may be None
    """ list of wavelengths for each dataset """
    if len(self.xmlnode.findall("ReflectionData/Dataset/Wavelength"))>0:
      datasets = self.xmlnode.findall("ReflectionData/Dataset")
      wavelengths = []
      for dataset in datasets:
        wavelengths.append(dataset.findall("Wavelength")[0].text)
      return wavelengths

    return None
  # - - - - - - - - - - - - - - - - -
  def ResultTable(self,parent=None, extraitems=None):
    # extraitems if present, for each dataset, a list for each item
    # [title, tip, [for each dataset, [data, 3 numbers Overall, Inner, Outer]]
    self.processDatasets()
    datasettext = ""
    datasetnumberlist = []
    if self.numberofdatasets == 0:
      return
    
    if self.numberofdatasets == 1:
      datasetname = self.datasetresultnodes[0].get("name")
      parent.append("Summary of merging statistics for dataset "+"<br/>"+datasetname)
    else:
      parent.append("Summary of merging statistics for multiple datasets")
      #print("self.datasetlabellist", self.datasetlabellist)
      datasettext = "Datasets: "
      jset = int(0)
      notfirst = False
      for label in self.datasetlabellist:
        jset += 1
        numbers = []
        if notfirst:
          datasettext += ", "
          numbers = [""]
        notfirst = True
        datasettext += "("+str(jset)+") "+label
        numbers.extend([str(jset), str(jset), str(jset)])
        datasetnumberlist.extend(numbers)
      
      #print("datasettext", datasettext)
      parent.append(datasettext)

    table = parent.addTable(select="Result", transpose=True,
                            style="line-height:100%; font-size:100%;",downloadable=True, class_="center")

    if self.numberofdatasets > 1:
      table.addData(title="Dataset number", data=datasetnumberlist)

    headers = []
    for i in range(self.numberofdatasets):
      headers += ["Overall", "Inner", "Outer"]
      if (i+1 < self.numberofdatasets): headers += " "

    table.addData(title="", data=headers)

    # List of header, tag, tool-tip
    taglist = \
        [["Low resolution limit", "ResolutionLow",""],
         ["High resolution limit", "ResolutionHigh",""],
#         ["Rmerge(within I+/I-)*", "Rmerge",
#          '\u2211 \u2211 | Ihl - <Ih> |/ \u2211  <Ih>'],
#         ["Rmerge(all I+ and I-)*", "RmergeOverall",
#          '\u2211 \u2211 | Ihl - <Ih> |/ \u2211  <Ih>'],
#         ["Rmeas (within I+/I-)*", "Rmeas",
#          '\u2211 \u2211 \u221A(n/n-1) | Ihl - <Ih> |/ \u2211  <Ih>'],
#         ["Rmeas (all I+ & I-)*", "RmeasOverall",
#          '\u2211 \u2211 \u221A(n/n-1) | Ihl - <Ih> |/ \u2211  <Ih>'],
#         ["Rpim (within I+/I-)", "Rpim",
#          '\u2211 \u2211 \u221A(1/n-1) | Ihl - <Ih> |/ \u2211  <Ih>'],
#         ["Rpim (all I+ & I-)", "RpimOverall",
#          '\u2211 \u2211 \u221A(1/n-1) | Ihl - <Ih> |/ \u2211  <Ih>'],
         ["Rmerge(within I+/I-)*", "Rmerge", ''],
         ["Rmerge(all I+ and I-)*", "RmergeOverall", ''],
         ["Rmeas (within I+/I-)*", "Rmeas", ''],
         ["Rmeas (all I+ & I-)*", "RmeasOverall", ''],
         ["Rpim (within I+/I-)", "Rpim", ''],
         ["Rpim (all I+ & I-)", "RpimOverall",  ''],
         ["Rmerge in top intensity bin*", "RmergeTopI",''],
         ["Number of observations", "NumberObservations",''],
         ["Number unique", "NumberReflections",''],
         ["Mean((I)/sd(I))", "MeanIoverSD",''],
         ["Half-set correlation CC(1/2)", "CChalf",''],
         ["Completeness %", "Completeness",''],
         ["Multiplicity", "Multiplicity",''],
         ["Filtered Mean(chi^2)", "MeanChiSq",'']]

    for label, tag, tip in taglist:
      self.AddToTable(datasetresultnodes=self.datasetresultnodes, \
                      table=table, label=label, tag=tag, tip=tip)

    # extra items if present
    #print("extraitems",extraitems)
    if extraitems is not None:
      for extra in extraitems:
        label = extra[0]
        tip = extra[1]
        data = []
        for i in range(self.numberofdatasets):
          num3 = extra[i+2]
          data += num3
          if i < self.numberofdatasets-1:
            data += ' '
      table.addData(title=label, data=data, tip=tip)
      
    if self.hasAnomalous:
      taglist = \
         [["Anomalous completeness %", "AnomalousCompleteness",''],
         ["Anomalous multiplicity", "AnomalousMultiplicity",''],
         ["DelAnom CC(1/2)", "AnomalousCChalf",''],
         ["Mid-Slope of Anom Probability", "AnomalousNPslope",'']]

    for label, tag, tip in taglist:
      self.AddToTable(datasetresultnodes=self.datasetresultnodes, \
                      table=table, label=label, tag=tag, tip=tip)

    if len(self.xmlnode.findall("AnomalousStatus"))>0:
      parent.append(self.xmlnode.findall("AnomalousStatus")[0].text)

    parent.append("Note: statistics marked '*' apply to unaveraged data, others to the merged data")

  # - - - - - - - - - - - - - - - - -
  def getResolutionBins(self):
    # Extract resolution limits for Overall, Inner, Outer, for each dataset

    limits = []  # for each dataset
    for datasetresultnode in self.datasetresultnodes:  # loop datasets
      tag = "ResolutionLow"
      low = self.Get3numbers(datasetresultnode, tag)
      tag = "ResolutionHigh"
      high = self.Get3numbers(datasetresultnode, tag)
      lh = []
      for i in range(3):
        lh.append([low[i], high[i]])
      limits.append(lh)
    
    return limits

  # - - - - - - - - - - - - - - - - -
  def ResultTableMerged(self,parent=None, extraitems=None):
    # extraitems if present, a list for each item
    # [title,[for each dataset, [data, 3 numbers Overall, Inner, Outer], tip]
    # print("ResultTableMerged")
    
    self.processDatasets()
    datasettext = ""
    datasetnumberlist = []
    ndts = self.numberofdatasets

    datasetnames = []
    dtsnames = ""
    for node in self.datasetresultnodes:
      datasetnames.append(node.get("name").split("/")[2])
      dtsnames += datasetnames[-1] + ", "
    dtsnames = dtsnames[:-2]

    if ndts == 1:
      parent.append("Summary of statistics for merged dataset "+\
                    datasetnames[0])
    else:
      parent.append("Summary of statistics for merged datasets "+\
                    dtsnames)
      

    table = parent.addTable(select="Result", transpose=True,
                style="line-height:100%; font-size:100%;",downloadable=True, class_="center")

    count = 0
    headers = []
    dtsheaders = []
    for i in range(ndts):
      headers += ["Overall", "Inner", "Outer"]
      dtsheaders += [datasetnames[i], datasetnames[i], datasetnames[i]]
      count += 1
      if count < ndts:
        headers += [' ']
        dtsheaders += [' ']
    table.addData(title="Dataset", data=dtsheaders)
    table.addData(title="Range", data=headers)

    # List of header, tag, tool-tip
    taglist = \
        [["Low resolution limit", "ResolutionLow",""],
         ["High resolution limit", "ResolutionHigh",""],
         ["Number unique reflections", "NumberReflections",''],
         ["Mean((I)/sd(I))", "MeanIoverSD",''],
         ["Half-set correlation CC(1/2)", "CChalf",''],
         ["Completeness %", "Completeness",'']]

    for label, tag, tip in taglist:
      self.AddToTable(datasetresultnodes=self.datasetresultnodes, \
                      table=table, label=label, tag=tag, tip=tip)

    # extra items if present
    if extraitems is not None:
      for extra in extraitems:
        label = extra[0]
        tip = extra[1]
        data = []
        for i in range(ndts):
          num3 = extra[i+2]
          data += num3
          if i < ndts-1:
            data += ' '
      table.addData(title=label, data=data, tip=tip)

    if self.hasAnomalous:
      taglist = \
         [["Anomalous completeness %", "AnomalousCompleteness",''],
         ["Mid-Slope of Anom Probability", "AnomalousNPslope",'']]

      for label, tag, tip in taglist:
        self.AddToTable(datasetresultnodes=self.datasetresultnodes, \
                        table=table, label=label, tag=tag, tip=tip)

    parent.append('Note that CC(1/2) is derived from correlating I+ and I-')


  # - - - - - - - - - - - - - - - - -
  def getResolutionEstimates(self, datasetresultnodes=""):
    """ Pick up and store resolution estimates into
    EstimatesofResolution class, instance self.estimatesofresolution
    No output from here
    """

    numberofdatasets = len(datasetresultnodes)
    if numberofdatasets == 0:
      return
    
    # 2 entry types for each direction, CC1/2 and I/sd (or 3 for overall)
    numberoftypes = 2

    #  datasets may have different numbers of directions if some are missing!
    nreslimits = len(datasetresultnodes[0].findall("ResolutionLimitEstimate"))

    # get thresholds
    thresholdCC = 0.0
    thresholdIovsd = 0.0

    self.estimatesofresolution = []  # for each dataset

    for datasetresultnode in datasetresultnodes:  # loop datasets
      reslimitnodes = datasetresultnode.findall("ResolutionLimitEstimate")
      # First 2 or 3 will be Overall
      noverall = 2
      # Thresholds will be the same for each dataset, but are stored separately
      if len(reslimitnodes[0].findall(".//*[@type='CChalf']"))>0:
        thresholdCC = reslimitnodes[0].findall("Threshold")[0].text
      try:
          thresholdIovsd = reslimitnodes[1].findall("Threshold")[0].text
      except:
          pass
      # May be 3rd Overall one for different threshold
      thresholdIovsd2 = None
      try:
        if reslimitnodes[2].findall("Direction")[0].text == "Overall" and len(reslimitnodes[2].findall(".//*[@type='I/sd"))>0:
          noverall = 3
          thresholdIovsd2 = reslimitnodes[2].findall("Threshold")[0].text
      except:
        pass
      
      message = "Estimates of limits from CC(1/2) (threshold"+\
              str(thresholdCC)+") <br/>and Mn(I/sd) (threshold"+str(thresholdIovsd)+")"
      datasetname = datasetresultnode.findall("[@name]")[0].attrib["name"]

      self.estimatesofresolution.append(EstimatesofResolution(
        message, datasetname,
        thresholdCC, thresholdIovsd, thresholdIovsd2))

      # Overall limits, for each dataset
      nodata = 0
      maxres = 0

      # Overall + anisotropy (0, 2 or 3)
      ndirections = (len(reslimitnodes)-noverall)//numberoftypes + 1
    
      for j in range(ndirections):  # loop directions and types
        # for this direction, for each type and dataset
        if j == 0:
          k = 0
        else:
          k = (j-1)*numberoftypes + noverall  #  3, 5, 7 eg
        # assume numberoftypes = 2
        reslimitnodes12 = [reslimitnodes[k], \
                           reslimitnodes[k+1]]    #  CC, I/sd 
        direction = reslimitnodes12[0].findall("Direction")[0].text
        direction2 = reslimitnodes12[1].findall("Direction")[0].text

        if direction == "Overall":
          label = "Overall"
        else:
          if "plane" in direction:
            label = " - In  "+direction
          else:
            label = " - Along "+direction

        nodata = [False, False]
        maxres = [False, False]
        reslims = []
        for i in range(2):
          reslim = reslimitnodes12[i].findall("MaximumResolution")[0].text
          if float(reslim) == 0.0:
            reslim = "-"

          message = None
          if len(reslimitnodes12[i].findall("Message"))>0:
            message = reslimitnodes12[i].findall("Message")[0].text
            if message == "insufficient data":
              nodata[i] = True
            elif message == " == maximum resolution":
              maxres[i] = True

          # limit from CC(1/2), IovSd, +'*' if maximum
          reslims.append(reslim)

        # store direction, message
        self.estimatesofresolution[-1].addLimits(label, reslims,
                                                 nodata, maxres)
        
  # - - - - - - - - - - - - - - - - -
  def resolutionEstimateTable(self, parent=None):
    """ Generate table from self.estimatesofresolution """

    ndatasets = len(self.estimatesofresolution)
    if ndatasets == 0: return

    # Get header stuff from 1st dataset
    estimatedreso = self.estimatesofresolution[0]
    message = estimatedreso.message
    parent.append(message)

    table = parent.addTable(transpose=True, class_="center")
    pad = ["",""]   # 2 columns of padding

    # Overall limits, for each dataset
    data = ["CC(1/2)","Mn(I/sd)"]
    table.addData(title="", data=data)

    for estimatedreso in self.estimatesofresolution:  # loop datasets
      limits = estimatedreso.directionlimits
      datasetname = estimatedreso.datasetname

      if ndatasets > 1:
        table.addData(title=">>>> Dataset: "+datasetname, data=pad)

      isnodata = False
      ismaxres = False
      for resolimits in estimatedreso.directionlimits:
        # resolimits is list of:
        #  label, [reslimlist], [nodata], [maxres]
        #   lists have 2 items
        label = resolimits[0]
        reslim = resolimits[1]
        nodata = resolimits[2]
        maxres = resolimits[3]
        for i in range(2):
          if nodata[i]:
            isnodata = True
            reslim[i] = "-"
          if maxres[i]:
            ismaxres = True
            reslim[i] += "*"

        table.addData(title=label, data=reslim)

    s = ""
    if isnodata:
      s += "entries with insufficient data for an estimate are marked '-'"
    if ismaxres:
      if s != "":
        s += "<br/>"
      s += "estimates extending to the maximum resolution marked '*'"

    if s != "":
      parent.append(s)

  # - - - - - - - - - - - - - - - - -
  def getDataList(self, nodes, tag, npad):
    # return a list of parameters for one or more data sets
    data = []
    for node in nodes:
      data.append(node.findall(tag)[0].text)

    # Pad to npad items
    if len(data) < npad:
      for i in range(npad-len(data)):
        data.append("")

    return data

  # - - - - - - - - - - - - - - - - -
  def formatCell(self, cellxml, astext=False):
    # astext = False, return list[[a,b,c,][alpha,beta,gamma]
    # astext = True,  return plain text
    a = cellxml.findall('a')[0].text
    b = cellxml.findall('b')[0].text
    c = cellxml.findall('c')[0].text
    alpha = cellxml.findall('alpha')[0].text
    beta = cellxml.findall('beta')[0].text
    gamma = cellxml.findall('gamma')[0].text
    if astext:
      return a+b+c+alpha+beta+gamma
    else:
      return [[a,b,c],[alpha,beta,gamma]]
    
  # - - - - - - - - - - - - - - - - -
  def otherStatisticsMultidataset(self, cells, parent=None):
    # More than one dataset
    parent.append('Unit cells for each dataset:')
    dnames = []
    for dataset in self.datasetresultnodes:
      dnames.append(dataset.get('name'))   # dataset names

    dnamelines = []
    for dname in dnames:
      dnamelines.append(dname)
      dnamelines.append("   (angles)")

    celllines = []
    for cell in cells:
      celldata = self.formatCell(cell)  # [[a,b,c],[alpha,beta,gamma]]
      celllines.append(celldata[0])
      celllines.append(celldata[1])

    if len(dnamelines) == len(celllines):
      table1 = parent.addTable(transpose=True, class_="center")
      for i in range(len(dnamelines)):
        table1.addData(title=dnamelines[i], data=celllines[i])


  # - - - - - - - - - - - - - - - - -
  def otherStatistics(self,parent=None):
    # Other miscellaneous things
    if self.numberofdatasets == 0:
      return

    cells = self.xmlnode.findall("Result/Dataset/cell")

    if self.numberofdatasets > 1:
      self.otherStatisticsMultidataset(cells, parent)

    else:
      # one dataset, usual case
      dname = (self.datasetresultnodes[0].get('name'))   # dataset name
      parent.append('Dataset:'+dname)
      table = parent.addTable(transpose='true', class_="center")
      cell = self.formatCell(cells[0])
      table.addData(title='Unit cell: a,b,c', data=cell[0])
      table.addData(title='  alpha,beta,gamma', data=cell[1])

      data = self.getDataList(self.datasetresultnodes,
                              "Mosaicity", 3)
      table.addData(title='Mosaicity', data=data)
      data = self.getDataList(self.datasetresultnodes,
                              'NumberReflections/Overall', 3)
      table.addData(title='Number of reflections', data=data)
      data = self.getDataList(self.datasetresultnodes,
                              'NumberObservations/Overall', 3)
      table.addData(title='Number of observations', data=data)

      outliers = self.xmlnode.findall('Outliers')

      data = self.getDataList(outliers, 'RejectNumberUnique', 3)
      table.addData(title='Number of rejected outliers', data=data)
      data = self.getDataList(outliers, 'RejectNumberEmax', 3)
      table.addData(title='Number of Emax rejects', data=data)


  # - - - - - - - - - - - - - - - - -
  def ReflectionDataDetails(self,parent=None, infold=True):
    if infold:
      fold = parent.addFold(label="Details of reflection data",brief='InputData')
    else:
      fold = parent
    
    fold.append("Summary of input reflection data")

    table = fold.addTable( select = "ReflectionData", transpose=False)
    for title,select in [ [ "Max resolution", "ResolutionHigh" ],
                          [ "Nreflections", "NumberReflections" ],
                          [ "NObservations", "NumberObservations" ],
                          [ "Nparts", "NumberParts" ],
                          [ "Nbatches", "NumberBatches" ],
                          [ "Ndatasets", "NumberDatasets" ] ]:
      table.addData( title=title , select = select )

    # list of run data
    datasetlist = self.xmlnode.findall("ReflectionData/Dataset")
    table2 = fold.addTable()

    # List of data wanted
    datasetnames = []
    runnumbers = []
    batchranges = []
    excludedbatches = [] 
    for dataset in datasetlist:  # loop datasets
      runlist = dataset.findall("Run")
      for run in runlist:  # loop runs
        datasetnames.append(dataset.findall("[@name]")[0].attrib["name"])
        runnumbers.append(run.findall("number")[0].text)
        batchranges.append(run.findall("BatchRange")[0].text)
        if len(run.findall("ExcludedBatches"))>0:
            excludedbatches.append(run.findall("ExcludedBatches")[0].text)
        
    table2.addData(title="DatasetName", data=datasetnames)
    table2.addData(title="RunNumber", data=runnumbers)
    table2.addData(title="Batch range", data=batchranges)
    table2.addData(title="Excluded batches", data=excludedbatches)

  # - - - - - - - - - - - - - - - - -
  def Graphs(self,parent=None, select="All", prioritise="CompletenessVsResolution"):
    'Plot All or selected graphs:'
    '  selections = "All", "Batch", "notBatch", "Completeness"'
    ' prioritise CompletenessVsResolution if set (default)'
    # A GraphGroup is a group of graphs displayed in the same graph viewer widget
    graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;")
    # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph

    # Loop over all Graph tables in the program output and add to the GraphGroup
    # The plotting instructions are provided as xml text
    graphlist = self.xmlnode.findall("CCP4Table[@groupID='Graph']")

    # Find prioritise graph if wanted
    priorityID = None
    if prioritise is not None:
      for thisgraph in graphlist:
        graphID = thisgraph.get("id")
        if prioritise in graphID:
          priorityID = graphID
          graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
          graph = graph.addPimpleData(xmlnode=thisgraph)

    if select == "Completeness" and "Completeness" in priorityID:
      return  # Skip the rest

    for thisgraph in graphlist:
      graphID = thisgraph.get("id")
      # Batch graph IDs contain either "RotationRange" or "Batch"
      isBatch = False
      if ("RotationRange" in graphID) or ("Batch" in graphID):
        isBatch = True
      plotit = False
      if prioritise is not None:
        if prioritise in graphID:
          continue   # skip this graph in main list

      if select == "All": plotit = True
      if select == "Batch" and isBatch: plotit = True
      if select == "notBatch" and not isBatch: plotit = True

      if plotit:
        # print("Graph ", graphID, isBatch)
        #print "Graph type", type(thisgraph)
        graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
        graph = graph.addPimpleData(xmlnode=thisgraph)
      
  # - - - - - - - - - - - - - - - - -
  def ByResolutionGraphs(self,parent=None):
      # A GraphGroup is a group of graphs displayed in the same graph viewer widget
      #      graphgroup = parent.addGraphGroup(style="width:370px;  height:300px;")
      graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;")
      # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        
      # Loop over all Graph tables in the program output and add to the GraphGroup
      # The plotting instructions are provided as xml text

      graphlist = self.xmlnode.findall("CCP4Table[@id='Graph-CChalf']")
      graphlist.extend(self.xmlnode.findall("CCP4Table[@id='Graph-StatsVsResolution']"))
      graphlist.extend(self.xmlnode.findall("CCP4Table[@id='Graph-Anisotropy']"))
      
      for thisgraph in graphlist:
        #  print "Aimless ByResolutionGraphs thisgraph", thisgraph.get("title")
        graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
        graph = graph.addPimpleData(xmlnode=thisgraph)

        # - - - - - - - - - - - - - - - - -
  def ByResolutionGraphsMerged(self,parent=None):
      # A GraphGroup is a group of graphs displayed in the same graph viewer widget
      #      graphgroup = parent.addGraphGroup(style="width:370px;  height:300px;")
      # print("ByResolutionGraphsMerged NDTS", self.numberofdatasets)
      byresoDiv = parent.addDiv(
        style="text-align:center;margin:0px auto; padding:3px;")
      graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;")
      # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        
      # Loop over all Graph tables in the program output and add to the GraphGroup
      # The plotting instructions are provided as xml text
      graphlist = self.xmlnode.findall("CCP4Table[@id='Graph-StatsVsResolution']")
      if len(self.xmlnode.findall("CCP4Table[@id='Graph-Anisotropy']"))>0:
        graphlist.extend(self.xmlnode.findall("CCP4Table[@id='Graph-Anisotropy']"))
      
      for thisgraph in graphlist:
        title = thisgraph.get("title")
        graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.attrib["title"] )
        graph = graph.addPimpleData(xmlnode=thisgraph)
    
  # - - - - - - - - - - - - - - - - -
  def ByBatchGraphs(self,parent=None):
      #print("ByBatchGraphs")
      # A GraphGroup is a group of graphs displayed in the same graph viewer widget
##      graphgroup = parent.addGraphGroup(style="width:370px; height:300px;")
      graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;")
      # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
      
      # Loop over all Graph tables in the program output and add to the GraphGroup
      # The plotting instructions are provided as xml text
      # Make Rmerge vs. batch the top graph
      selection = "CCP4Table[@id='Graph-StatsVsBatch']"
      graphlist = self.xmlnode.findall(selection)
      if not len(self.xmlnode.findall("OnlyMerge"))>0:
        selection = "CCP4Table[@id='Graph-ScalesVsRotationRange']"
        graphlist.append(self.xmlnode.findall(selection)[0])

      for thisgraph in graphlist:
        #print("Aimless ByBatchGraphs thisgraph", thisgraph.get("title"))
        graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
        graph = graph.addPimpleData(xmlnode=thisgraph)
  # - - - - - - - - - - - - - - - - -
  def referenceGraph(self,parent=None):
    # graph vs batch
    # A GraphGroup is a group of graphs displayed in the same graph viewer widget

    style = "width:700px;  height:270px;"
    graphgroup = parent.addFlotGraphGroup(style=style)
    # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
      
    # Loop over all Graph tables in the program output and add to the GraphGroup
    # The plotting instructions are provided as xml text

    # for Batch make Rmerge vs. batch the top graph
    selection = "CCP4Table[@id='Graph-RefStatsVsBatch']"
    graphs = self.xmlnode.findall(selection)
    refdatalist = []
    for thisgraph in graphs:
      refdatalist.append(self.getRefGraphInfo(thisgraph))
      graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
      graph = graph.addPimpleData(xmlnode=thisgraph)

      s = 'Overall comparison with reference data: correlation coefficients and R-factors\n'
      for refdata in refdatalist:
        s += "Dataset: {};  Mn(CCref) {:7.3f}; Mn(Rref) {:7.3f}; number {:9d}\n".\
            format(refdata[0], refdata[1], refdata[2], refdata[3])

      s = html_linebreak(s, False)
      parent.append(s)

    # - - - - - - - - - - - - - - - - -
  def referenceGraphReso(self,parent=None):
    # resolution
    style = "width:370px; height:300px;"
      
    graphgroup = parent.addFlotGraphGroup(style=style)
    # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
      
    # Should be just one graph
    # The plotting instructions are provided as xml text

    selection =  "CCP4Table[@id='Graph-DatasetRefStatsVsReso']"
    graphs = self.xmlnode.findall(selection)
    thisgraph = graphs[0]
    refdata = self.getRefGraphResoInfo(thisgraph)
    graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
    graph = graph.addPimpleData(xmlnode=thisgraph)

    s = 'Overall comparison with reference data: correlation coefficients and R-factors\n'
    s = html_linebreak(s, False)
    parent.append(s)

    table = parent.addTable(transpose=True,
                            style="line-height:100%; font-size:100%;"
                            ,downloadable=True, class_="center")
    ndts = len(refdata[0])
    headers = ['CCref', 'Rref', 'Num']
    table.addData(title="Dataset", data=headers)

    for i in range(ndts):
      cc = "{:.3f}".format(refdata[0][i])
      r  = "{:.3f}".format(refdata[1][i])
      num = "{}".format(refdata[2][i])
      data = [cc, r, num]

      dname = self.datasetnames[i].split('/')[-1]
      table.addData(title=dname, data=data)


  # - - - - - - - - - - - - - - - - -
  def getRefGraphInfo(self, refgraph):
    # Extract information from graph of agreement with reference by batch
    title = refgraph.get('title')
    # print(">> getRefGraphInfo", title)
    dname = title.split()[-1]
    # print('dname', dname)
    headers = refgraph.findall('headers')[0].text.split()
    # print(headers)
    fail = False
    avgrref = 0.0
    avgccref = 0.0
    ndata = 0
    try:
      irref = headers.index('Rref')
      iccref = headers.index('CCref')
      inum = headers.index('Number')
    except:
      fail = True
    if not fail:
      data = refgraph.findall('data')[0].text
      # get weighted average of Rref and CCref
      #  there are clever Pythonic ways of doing averages, but this will do
      lines = data.splitlines()
      sumwrrr = 0.0
      sumwccr = 0.0
      sumw = 0.0
      for line in lines:
        fields = line.split()
        if len(fields) > 0:
          num = float(fields[inum])
          if num > 0.0:
            w = 1.0/num
            sumwrrr += w * float(fields[irref])
            sumwccr += w * float(fields[iccref])
            sumw += w
            ndata += num

      if sumw > 0.0:
        avgrref = sumwrrr/sumw
        avgccref = sumwccr/sumw

    return [dname, avgccref, avgrref, int(ndata)]
      
  # - - - - - - - - - - - - - - - - -
  def getRefGraphResoInfo(self, refgraph):
    # Extract information from graph of agreement with reference by resolution
    title = refgraph.get('title')
    # print(">> getRefGraphResoInfo", title)
    headers = refgraph.findall('headers')[0].text.split()
    fail = False

    ndts = 0
    irref = []
    iccref = []
    inumref = []
    # eg N  1/d^2   Dmid     R-1     CC-1     N-1    R-2     CC-2     N-2
    for idx, label in enumerate(headers):
      if label[0] == 'R':
        irref.append(idx)
      elif label[0] == 'C':
        iccref.append(idx)
      elif label[0:2] == 'N-':
        inumref.append(idx)
    # lists of indices to the columns we want
    fail = False
    if len(irref) == 0 or len(iccref) == 0: fail = True

    if not fail:
      ndts = len(irref)
      data = refgraph.findall('data')[0].text
      # get weighted average of R-x and CC-x
      #  there are clever Pythonic ways of doing averages, but this will do
      lines = data.splitlines()
      sumwrrr = [0.0] * ndts
      sumwccr = [0.0] * ndts
      sumw = [0.0] * ndts
      ndata = [0] * ndts
      for line in lines:
        fields = line.split()
        if len(fields) > 0:
          for i in range(ndts):
            num = self.numconvert(fields[inumref[i]])
            if num is not None and num > 0.0:
              w = 1.0/num
              r = self.numconvert(fields[irref[i]])
              cc = self.numconvert(fields[iccref[i]])
              if r is not None and cc is not None:
                sumwrrr[i] += w * r
                sumwccr[i] += w * cc
                sumw[i] += w
                ndata[i] += int(num)

        avgrref = []
        avgccref = []
        for i in range(ndts):
          if sumw[i] > 0.0:
            avgrref.append(sumwrrr[i]/sumw[i])
            avgccref.append(sumwccr[i]/sumw[i])

    return [avgccref, avgrref, ndata]
      
  # - - - - - - - - - - - - - - - - -
  def numconvert(self, s):
    # return float if s is a number, else None
    try:
      f = float(s)
    except:
      f = None
    return f
  # - - - - - - - - - - - - - - - - -
  def formatTable(self, tablein, title, celltitle, datalabel,
                  includenumbers, parent=None):
    # format a table such as a cross-correlation table
    # tablein is list of elements
    # title is a header for the table
    # celltitle is a label for columns & rows
    # datalabel is the tag for the value
    # includenumbers True to add the numbers in brackets
    # <columnsheaders> contains <label>
    # <rows> contains <row> contains <datalabel>[, <Number>]

    columns = tablein.findall('columnheaders')[0]
    columnlabels = columns.findall('label')

    rowlabs = []
    for lab in columnlabels:
      rowlabs.append(lab.text)
    collabs = rowlabs
    
    # truncate columns if needed
    MAXNUMINLINE = 18  #  maximum number / line (all rows)
    ninline = min(len(columnlabels), MAXNUMINLINE)
    if ninline < len(columnlabels):
      collabs = collabs[0:ninline]
      collabs.append(' ...')

    rows = tablein.findall('row')
    #  true if the table has Number elements in rows
    hasNumber = False
    #print "rows", type(rows)

    rowlines = []
    rowlabels = []
    for row in rows:
      rowdata = row.findall(datalabel)
      rownumbers = row.findall('Number')
      #print "rl", row.findall('label').text
      rowlabels.append(row.findall('label')[0].text)
      line = []
      for i in range(ninline):
        d = rowdata[i]
        lab = ""
        if (d.text != " "):
          lab = d.text
          if includenumbers and (len(rownumbers) > 0):
            n = rownumbers[i]
            lab += " (" + n.text + ")"
            hasNumber = True
        line.append(lab)
      
      if ninline < len(rowdata):
        line.append(' ...')
    
      rowlines.append(line)

    header = title
    if hasNumber:
      header += "<br/>(numbers in brackets)"

    parent.append(header)

    table = parent.addTable(transpose="true", class_="center")
    table.addData(title=celltitle, data=collabs)

    for rowlabel, line in zip(rowlabels, rowlines):
      #print rowlabel
      table.addData(title=rowlabel, data=line)

  # - - - - - - - - - - - - - - - - -
  def interdatasetAnomalousTable(self,parent=None):
    #  Graph and table of interdataset correlations
    intertable = self.xmlnode.findall('crosstable[@id="DatasetAnomalousCorrelation"]')
    if intertable != None and len(intertable) != 0:
      self.formatTable(intertable[0],
           "Correlation of DelAnom between datasets", "Dataset", "CC",
                       True, parent)

  # - - - - - - - - - - - - - - - - -
  def interdatasetDispersionTable(self,parent=None):
    #  Table of interdataset correlations
    intertable = self.xmlnode.findall('crosstable[@id="DatasetDispersiveCorrelation"]')
    if intertable != None and len(intertable) != 0:
      self.formatTable(intertable[0],
          "Correlation of dispersive difference between datasets",
                       "Dataset", "CC", True, parent)
  # - - - - - - - - - - - - - - - - -
  def interdatasetIntensityTable(self,parent=None):
    #  Table of interdataset correlations
    intertable = self.xmlnode.findall('crosstable[@id="DatasetIntensityCorrelation"]')
    if intertable != None and len(intertable) != 0:
      self.formatTable(intertable[0],
          "Correlation of intensity difference between datasets",
                       "Dataset", "CC", True, parent)

  # - - - - - - - - - - - - - - - - -
  def interRunTable(self,parent=None):
    #  Table of run-run correlations
    intertable = self.xmlnode.findall('crosstable[@id="RunCorrelation"]')
    if (intertable is not None) and (len(intertable) != 0):
      parent.append('Correlation between runs')
      self.formatTable(intertable[0],
          "", "Run", "CC", False, parent)

  # - - - - - - - - - - - - - - - - -
  def interdatasetIntensityGraph(self,parent=None):
    #  Graph of interdataset correlations
    #
    #gt = "=== Correlation of intensities between datasets"
    #gpath = "CCP4Table[@title='" + gt + "']"
    #graph = self.xmlnode.findall(gpath)[0]
    #  change to this for Aimless 0.7.13+
    graph = self.xmlnode.findall("CCP4Table[@id='Graph-Intensities']")
    if len(graph) > 0:
      graph = graph[0]
      title =  graph.attrib['title']
      agraph = parent.addFlotGraph(xmlnode=graph, title=title)
##                             style="width:370px;" )
      agraph.addPimpleData(xmlnode=graph)

  # - - - - - - - - - - - - - - - - -
  def interdatasetAnomalousGraph(self,parent=None):
    #  Graph of interdataset correlations
    graph = self.xmlnode.findall("CCP4Table[@id='Graph-AnomalousDifferences']")
    #graph = self.xmlnode.findall("CCP4Table[@id='Graph-AnomalousDifferences']")[0].text
    if len(graph) > 0:
      graph = graph[0]
      title =  graph.attrib['title']
      agraph = parent.addFlotGraph(xmlnode=graph, title=title)
##                             style="width:370px;" )
      agraph.addPimpleData(xmlnode=graph)

  # - - - - - - - - - - - - - - - - -
  def interdatasetDispersionGraph(self,parent=None):
    #  Graph of interdataset correlations
    graph = self.xmlnode.findall("CCP4Table[@id='Graph-DispersiveDifferences']")
    if len(graph) > 0:
      graph = graph[0]
      title =  graph.attrib['title']
      agraph = parent.addFlotGraph(xmlnode=graph, title=title)
##                             style="width:370px;" )
      agraph.addPimpleData(xmlnode=graph)

  # - - - - - - - - - - - - - - - - -
  def addScaleBfacReport(self, parent=None):
    #  Report on relative scales and B-factors for multiple datasets
    #  but only if there is a single batch/rotationrange per dataset
    # print("addScaleBfacReport", self.onebatchperdataset)
    if not self.onebatchperdataset: return
    graphs = self.xmlnode.findall("CCP4Table[@id='Graph-ScalesVsRotationRange']")
    if len(graphs) < 2: return  # Only if there are multiple datasets

    headerlist = []
    datalist = []
    datasetlist = []
    for graph in graphs:
      headerlist.append(graph.findall("headers")[0].text)
      datalist.append(graph.findall("data")[0].text)
      title = graph.attrib['title']
      datasetlist.append(title.split()[-1])

    scales = []
    bfacs = []
    for dts, hdr, data in zip(datasetlist, headerlist, datalist):
      dl = data.split()
      scales.append(dl[hdr.split().index('0k')])
      bfacs.append(dl[hdr.split().index('Bfactor')])

    kbDiv = parent.addDiv(style="width:100%;text-align:center;margin:0px; padding:0px; border:1px solid black;")

    # Use grid layout for scale/bfactor display
    leftDiv, rightDiv = kbDiv.addTwoColumnLayout(left_span=7, right_span=5, spacing=2)

    leftDiv.append("Relative scales and B-factors between merged datasets")
    table = rightDiv.addTable(transpose=True,
                            style="line-height:100%; font-size:100%;",
                            class_="center")
    headers = ['Scale', 'RelativeBfactor']
    table.addData(title="Dataset", data=headers)
    for dname, scale, bfac in zip(datasetlist, scales, bfacs):
      data = [scale, bfac]
      table.addData(title=dname, data=data)

  # - - - - - - - - - - - - - - - - -
  def NormalProbabilityPlot(self,parent=None):
    if len(self.xmlnode.findall(".//CCP4Table[@title='Normal probability plot']"))==0: return
    nmplot = self.xmlnode.findall(\
      ".//CCP4Table[@title='Normal probability plot']")[0]
    if nmplot is not None:
      graph = parent.addFlotGraph(launcher='Normal probability')
      graph.addPimpleData(xmlnode=nmplot)

  # - - - - - - - - - - - - - - - - -
  def AnomalousNormalProbabilityPlot(self,parent=None):
    if len(self.xmlnode.findall(".//CCP4Table[@title='Anomalous differences']"))==0: return
    nmplot = self.xmlnode.findall(\
      ".//CCP4Table[@title='Anomalous differences']")[0]
    if nmplot is not None:
      graph = parent.addFlotGraph(launcher='Anomalous Q-Q plot')
      graph.addPimpleData(xmlnode=nmplot)

  # - - - - - - - - - - - - - - - - -
  def AnomalousScatterPlot(self,parent=None):
    if len(self.xmlnode.findall(".//CCP4Table[@title='DelAnom/RMS scatter plot']"))==0: return
    nmplot = self.xmlnode.findall(\
      ".//CCP4Table[@title='DelAnom/RMS scatter plot']")[0]
    if nmplot is not None:
      graph = parent.addFlotGraph(launcher='DelAnom scatterplot')
      graph.addPimpleData(xmlnode=nmplot)

  # - - - - - - - - - - - - - - - - -
  def RoguePlot(self,parent=None):
    if len(self.xmlnode.findall(".//CCP4Table[@title='Outliers on detector (horizontal rotation axis)']"))==0: return
    nmplot = self.xmlnode.findall(\
      ".//CCP4Table[@title='Outliers on detector (horizontal rotation axis)']")[0]
    if nmplot is not None:
      graph = parent.addFlotGraph(launcher='Outlier positions')
      graph.addPimpleData(xmlnode=nmplot)

  # - - - - - - - - - - - - - - - - -
  def getISa(self):
    """ list of ISa values (often just one) """
    s = "Asymptotic I/sigI, ISa: "
    s += self.sdcpar.displayISa()
    return s
    
  # - - - - - - - - - - - - - - - - -
  def sdAnalysis(self, parent=None):
    ''' Plots and tables for analysis of sig(I) '''

    if self.sdcpar is None: return
    #print("SDANALYSIS")
    sdanaldiv = parent.addDiv()
    sdanaldiv.addText(text="Analysis of sd(I)", style='font-size: 150%;')
    s = "Error estimates sd(I) are analysed as functions of intensity (and resolution and batch),\n"+\
        "as sigma(delta(I)/sd(I)), where delta(I) = Ihl-<Ih>,"+\
        " and as reduced chi^2 (goodness of fit), ([delta(I)/sd(I)]^2)/(n-m)"
    s = html_linebreak(s)
    sdanaldiv.append(s)

    sdcstatus, sdcmessage = self.sdcpar.status()
    if sdcstatus > 0:
      if sdcstatus == +1:
        colour = 'darkorange'
      elif sdcstatus == +2:
        colour = 'red'
      sdanaldiv.addText(text=sdcmessage+", see flags in table", style='color:'+colour+';')

    # Table of SD correction parameters
    self.sdcpar.sdcorrectionParameterTable(sdanaldiv)

    # With multiple Runs, there may be a graph for each run, then one for "All runs"
    #  do we have a multiple runs?
    sdanals = self.xmlnode.findall("CCP4Table[@id='Graph-SDanalysis']")
    ngraphs = len(sdanals)

    if ngraphs > 0:
      div1 = sdanaldiv.addDiv(style="width:45%; margin:0px; padding:0px;display:inline-block;")
      if ngraphs > 1:
        div1a = sdanaldiv.addDiv(style="width:45%; margin:0px; padding:0px;display:inline-block;")

        div1.addText(text="Analysis against intensity", style='font-size: 110%;')

      graphs = self.getSelectedGraphs(graphID="Graph-SDanalysis", graphTitle=" All runs",
                                      plotTitleList=None)
      if graphs is not None:
        div1.append("For all runs")
        self.addSelectedGraphs(graphs, div1)

      graphs = self.getSelectedGraphs(graphID="Graph-SDanalysis", graphTitle=" Run ",
                                    plotTitleList=None)
      if ngraphs > 1:
        # An all runs graph has been plotted
        div1a.append("Plots for individual runs")
        self.addSelectedGraphs(graphs, div1a)
      elif ngraphs > 0:
        self.addSelectedGraphs(graphs, div1)

    resographs = self.getSelectedGraphs(graphID="Graph-StatsVsResolution",
                                        plotTitleList=["Filtered Mean(Chi^2)",
                                                       "Mean(Chi^2), Mean(Chi^2)"])
    batchgraphs = self.getSelectedGraphs(graphID="Graph-StatsVsBatch",
                                        plotTitleList=["Filtered Mean(Chi^2)",
                                                       "Mean(Chi^2), Mean(Chi^2)"])

    if resographs is not None:
      div2a = sdanaldiv.addDiv(style="width:48%; margin:0px; padding:0px;display:inline-block;")
      div2a.append("Analysis against resolution")
      self.addSelectedGraphs(resographs, div2a)
    if batchgraphs is not None:
      div2b = sdanaldiv.addDiv(style="width:48%; margin:0px; padding:0px;display:inline-block;")
      div2b.append("Analysis against batch")
      self.addSelectedGraphs(batchgraphs, div2b)

  # - - - - - - - - - - - - - - - - -
  def getSelectedGraphs(self, graphID=None, graphTitle=None, plotTitleList=None):
    ''' Select graphs (plots in titleList) from CCP4Table with id=graphID'''
    if graphID is None: return None
    #print "getSelectedGraphs", graphID, graphTitle, plotTitleList
    return selectGraphs(self.xmlnode, graphID=graphID, graphTitle=graphTitle,
                        plotTitleList=plotTitleList)
  
  # - - - - - - - - - - - - - - - - -
  def addSelectedGraphs(self, graphs, parent=None):
    ''' Select graphs (plots in titleList) from CCP4Table with id=graphID and add to parent '''
    ngraphs = 0
    if graphs is not None:
      # if there is more than one graph, put them into a group
      where = parent
      if len(graphs) > 1:
        where = parent.addFlotGraphGroup(style="width:300px;  height:270px;")

      for thisgraph in graphs:
        graph = where.addFlotGraph(xmlnode=thisgraph,style="width:300px;  height:270px;")
        graph.addPimpleData(xmlnode=thisgraph)
        ngraphs += 1

    return ngraphs  # number plotted
  # - - - - - - - - - - - - - - - - -
  # - - - - - - - - - - - - - - - - -
  class SDcorrectionParameters:

    def __init__(self, xmlnode):
      ''' Extract SD parameters, allowing for unfortunate change of format '''

      self.sdc_error = -1
      # list of tagged SD correction parameters
      #   SDcorrection object: SDfac, SDb, SDadd, ISa
      self.sdcorrparams = []
      if not len(xmlnode.findall('SDcorrection'))>0: return

      sdcp = xmlnode.findall('SDcorrection') # all SDcorrections
      # Old versions of Aimless write multiple blocks, and we want the last one
      #   but old style is not dealt with
      oldstyle = False
      if len(sdcp) > 1:
        oldstyle = True
      if oldstyle: return  # Old style is not dealt with

      sdcp = sdcp[-1]
      allruns = False
      #print "sdcp", sdcp, type(sdcp), len(sdcp)

      # Note: "Allruns" is present when there are multiple runs all treated as same
      # with just one run, the block is flagged as "Run 1"
      if len(sdcp.findall('AllRuns'))>0:
        #print "allruns"
        allruns = True

      if allruns:
        run = sdcp.findall('AllRuns')[0]
        for tag in ['Fulls', 'Partials']:
          px = run.findall(tag)
          if len(px) > 0:
            # Fulls or partials found
            for pp in px:
              label = "AllRuns"
              sdcr = self.getSDset(pp, label, tag)
              #print "sdcr allruns",sdcr.as_string()
              self.sdcorrparams.append(sdcr)
      else:
        runs = sdcp.findall('Run')
        for run in runs:
          runnum = int(run.findall('number')[0].text)
          ntag = len(run.findall('Fulls'))
          ntag += len(run.findall('Partials'))
          for tag in ['Fulls', 'Partials']:
            px = run.findall(tag)
            if len(px) > 0:
              # Fulls or partials found
              for pp in px:  # Are there ever more than one?
                label = "Run {}".format(runnum)
                # suppress Full|Partial if only one
                fptag = tag
                if ntag == 1: fptag = ""
                sdcr = self.getSDset(pp, label, fptag)
                self.sdcorrparams.append(sdcr)
                
      # We now have a list self.sdcorrparams of SDcorrectionData objects,
      # so analyse them (class SDcorrectionData is in aimless_pipe_utils)
      #print "SDC", len(self.sdcorrparams)

      nparsets = len(self.sdcorrparams)  # excluding All runs

      # Get median values of SDfac and SDadd
      if nparsets > 1:
        sdfacs = []
        sdadds = []
        for sdc in self.sdcorrparams:
          sdfacs.append(sdc.SDfac)
          sdadds.append(sdc.SDadd)
        mediansdfac = median(sdfacs)
        mediansdadd = median(sdadds)

      self.sdc_error = 0  # flag for any errors
      self.error_flag = []   # this will be the same length as sdcorrparams

      for sdc in self.sdcorrparams:
        #print "SDC contents", sdc.as_string()
        #print "SDC contents", sdc.as_list()
        if nparsets > 1:
          sdc.validityCompare(mediansdfac, mediansdadd)

        # = 0 if no error, +1 minor, +2 major
        self.sdc_error = max(self.sdc_error, sdc.valid())
        self.error_flag.append(sdc.status)

      #print "sdcerrors", self.error_flag, self.sdc_error
    
    # - - - - - - - - - - - - - - - - -
    def status(self):
      ''' return status and warning message, = '' if none '''
      #  = 0 if no error, +1 minor, +2 major '''
      s = ''
      if self.sdc_error <= 0:
        return self.sdc_error, s

      if self.sdc_error == +1:
        s = "WARNING: SD correction parameters are somewhat outside expected limits"
      elif self.sdc_error == +2:
        s = "SEVERE WARNING: SD correction parameters are outside expected limits"
      return self.sdc_error, s
      
    # - - - - - - - - - - - - - - - - -
    def getSDset(self, element, label, tag):
      ''' make SDcorrection object: SDfac, SDb, SDadd, ISa '''
      sdfac = float(element.findall('SDfac')[0].text)
      sdb   = float(element.findall('SDb')[0].text)
      sdadd = float(element.findall('SDadd')[0].text)
      isa = float(element.findall('ISa')[0].text)

      sdcd = SDcorrectionData(label, tag, sdfac, sdb, sdadd, isa)
      return sdcd

    # - - - - - - - - - - - - - - - - -
    def displaySDcorrectionParameters(self,parent=None):
      ''' call SDcorrectionParameters first to extract values '''
      if len(self.sdcorrparams) == 0: return
      sdheader = "Parameters for improvement of sd(I) estimates: "
      sdheader += "sd'(I) = SdFac * Sqrt[sd(I)^2 + SdB I + (SdAdd I)^2]\n"
      sdheader += SDcorrectionData.header((len(self.sdcorrparams)>1))
      parent.append(html_linebreak(sdheader))
      
      for sdc in self.sdcorrparams:
        parent.append(html_linebreak(sdc.as_string()))

    # - - - - - - - - - - - - - - - - -
    def displayISa(self):
      """ ISa values """
      s = ""
      for sdc in self.sdcorrparams:
        label, isa = sdc.labelISa()
        if isa <= 0.0:
          s += ' {}: cannot determine ISa'.format(label)
        else:
          s += ' {}: ISa = {},'.format(label,isa)
          s = s[:-1]   # strip final comma
      return s
      
    # - - - - - - - - - - - - - - - - -
    def sdcorrectionParameterTable(self,parent=None):
      ''' call SDcorrectionParameters first to extract values '''
      if len(self.sdcorrparams) == 0:
        parent.append("SD correction parameters were not refined")
        return
      
      sdheader = "Parameters for improvement of sd(I) estimates: "
      sdheader += "sd'(I) = SdFac * Sqrt[sd(I)^2 + SdB I + (SdAdd I)^2]\n"

      sdheader += SDcorrectionData.header((len(self.sdcorrparams)>1))
      parent.append(html_linebreak(sdheader))
      
      table = parent.addTable(transpose=True, class_="center")
      table.addData(title='',
                    data=['SdFac','flag','SdB','SdAdd','flag','ISa'])

      for sdc in self.sdcorrparams:
        label, data = sdc.as_list()
        table.addData(title=label, data=data)
