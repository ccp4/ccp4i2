from __future__ import print_function

import os,sys
try:
  from ccp4i2.report.CCP4ReportParser import *
except:
  exec(compile(open(os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc')).read(), os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc'), 'exec'))
  from ccp4i2.report.CCP4ReportParser import *

try:
  import pointless_report
  import aimless_report
  import ctruncate_report
  from aimless_pipe_utils import *
  import phaser_analysis_report
except:
  from ccp4i2.wrappers.pointless.script import pointless_report
  from ccp4i2.wrappers.aimless.script import aimless_report
  from ccp4i2.wrappers.ctruncate.script import ctruncate_report
  from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import *
  from ccp4i2.wrappers.phaser_analysis.script import phaser_analysis_report

class import_merged_report(Report):
  TASKNAME = 'import_merged'
  RUNNING = True

  #def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
    ##self.aimlessreport = aimless_pipe_report.aimless_pipe_report(xmlnode=xmlnode,jobStatus='nooutput')

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    if kw.get('jobStatus',None) is not None and kw.get('jobStatus').lower() == 'nooutput':
      return
    elif kw.get('jobStatus',None) is not None and kw.get('jobStatus').lower() == 'running':
      self.runningReport()
    else:
      self.finalReport()


  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def runningReport(self, parent=None):
    if parent is None: parent = self

    topDiv = parent.addDiv(style="font-size:90%;")
    self.importReport(topDiv)

    # FreeR fail?
    self.omitmessage = None
    self.freeRmessage(parent)

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def reportstatus(self):
    '''
    Returns status, messagelist:
    status is
      "EMPTY"     nothing there
      "IMPORTING" running the import step
      "DRPIPE_RUNNING" running DR pipeline, aimless_pipe
      "FINISHED"  finished
    '''
    #print("reportstatus")
    status = "EMPTY"
    if not len(self.xmlnode.findall("IMPORT_LOG"))>0: return status

    text = []
    self.importxml = self.xmlnode.findall("IMPORT_LOG")[0] # all information from import_merged task

    # what is it doing now?
    if len(self.xmlnode.findall("AIMLESS_PIPE"))>0:
      status = "FINISHED"
    elif len(self.importxml.findall("DRPIPE_RUNNING"))>0:
      # this may be there after finishing, but is overriden by AIMLESS_PIPE above
      status = "DRPIPE_RUNNING"
    else:
      status = "IMPORTING"

    return status

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def importReport(self, parent):
    # Format stuff from import step IMPORT_LOG
    print("importReport")

    status = self.reportstatus()

    message1 = 'file read'

    # what is it doing now?
    if status == "EMPTY":
      return
    elif status == "IMPORTING":
      parent.addText(text="File import is running",
                     style='font-weight:bold; font-size:150%; color:red;')
      message1 = 'file reading'

    elif status == "DRPIPE_RUNNING":
      parent.addText(text="Data reduction pipeline is running to generate statistics",
                     style='font-weight:bold; font-size:150%; color:red;')
    elif status == "FINISHED":
      pass
      #parent.addText(text="Data reduction pipeline has finished",
      #               style='font-weight:bold; font-size:150%; color:red;')


    self.importxml = self.xmlnode.findall("IMPORT_LOG")
    if len(self.importxml) == 0: return
    self.importxml = self.importxml[0] # all information from import_merged task

    self.IorFtype = "Unknown"

    merged = 'Unmerged' # shouldn't ever be this, always merged!
    if self.importxml.findall('merged')[0].text == 'True': merged = 'Merged'

    # possible formats:
    #   'mtz', 'mmcif', 'shelx', 'sca'
    fformat = self.importxml.findall('fileformat')[0].text
    message = ''
    if fformat == 'mtz':
      columnlabels = self.importxml.findall('columnlabels')[0].text
      message = 'Merged MTZ ' + message1
    elif fformat == 'mmcif':
      message = 'mmCIF ' + message1
    elif fformat == 'shelx':
      message = 'ShelX ' + message1
    elif fformat == 'sca':
      message = 'Merged Scalepack ' + message1

    importDiv = parent.addDiv(style="border: 1px solid black; margin: 20px 1px 1px 1px; padding:5px;font-size:95%;line-height:1.7")

    importDiv.addText(text='Information about the import file :',
                       style='font-weight:bold; font-size:110%; color:blue;')

    #importDiv.append("<br/>")
    #importDiv.append(message)
    importDiv.addText(text=message)
    if len(self.importxml.findall('filename'))>0:
      filename = self.importxml.findall('filename')[0].text
      message = "Input file name: " + os.path.basename(filename)
      #importDiv.append("<br/>")
      importDiv.addText(text=message)
    
    #importDiv.append("<br/>")

    if fformat == 'mtz':
      if columnlabels is None:
        # try to get labels from X2MTZ block
        if len(self.importxml.findall("X2MTZ/inputcolumnames"))>0:
          columnlabels = self.importxml.findall("X2MTZ/inputcolumnames")[0].text
      if columnlabels is not None:
        message = 'Column labels selected: '+columnlabels
        importDiv.addText(text=message)
        #importDiv.append("<br/>")

    self.IorFtype = "Unknown"
    if len(self.importxml.findall("IorFtype"))>0:
      self.IorFtype = self.importxml.findall("IorFtype")[0].text
      if self.IorFtype == "Unknown":
        message = "Imported values are of unknown type, maybe intensities or amplitudes"
      elif self.IorFtype == "Amplitude":
        message = \
                "Imported values are amplitudes, and will be squared for statistical analysis only"
      elif self.IorFtype == "Intensity":
        message = "Imported values are intensities"
      importDiv.addText(text=message)
      #importDiv.append(message)

    if fformat == 'sca':
      s = self.reportResolution(self.importxml)
      if s != '':
        importDiv.append('Selected r'+s[1:])

    if fformat == 'mtz':
      message = self.mtzreport()  # for X2MTZ block from Gemmi import
      #importDiv.append("<br/>")
      for line in message:
        importDiv.addText(text=line)
        #importDiv.append("<br/>")
    elif fformat == 'mmcif':
      mmcifDiv = importDiv.addDiv(style="padding:5px;background-color:#CCEEFF;")
      mmcifDiv.addText(text='Information about the mmCIFfile :',
                        style='font-weight:bold; font-size:110%; color:blue;')

      message = self.mmcifreport()   #  information about the mmCIF file
      #mmcifDiv.append("<br/>")
      for line in message:
        mmcifDiv.addText(text=line)
        #mmcifDiv.append("<br/>")
      #mmcifDiv.append(message)

    parent.addText(text="NB data from the original file are kept for use in downstream processing",
                   style='font-weight:bold; font-size:120%; color:blue;')


    return

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  # FreeR fail?
  def freeRmessage(self, parent=None):
    '''
    Check whether a FreeR set was created
    Two XML blocks, both may occur twice:
      a) newFreeR, == True new set, False completing or completed
      b) FreeRfailed == True if failed and tried again
         == Again if failed finally
    Also
      c) freeRsource: if present may be "Explicit" for explicit freeR source,
         or "Input" if from main input data,
         or "None" if FreeR generation was skipped entirely

    Cases:
     1) Completed existing FreeR succesfully on 1st attempt
        newFreeR == "False", no FreeRfailed block
        but maybe ObsFreeCellComparison block if incompatbility was overridden
     2) New FreeR set on first attempt, maybe because of mismatch
        newFreeR == "True", no FreeRfailed block
     3) Complete failed on 1st attempt, running 2nd attempt to make new set
         newFreeR == "False", FreeRfailed == "True"
     4) Complete failed on 1st attempt, 2nd attempt to make new set succeeded
        two newFreeR == "False", "True", FreeRfailed == "True"
     5) Complete failed on 1st attempt, 2nd attempt to make new set failes
        two newFreeR == "False", "True", two FreeRfailed == "True", Again
        Probably should not occur
     6) No generation of FreeR (SKIP_FREER), probably for StarAniso data
        freeRsource = 'None'
    '''
    print("freeRmessage")
    headlines = []
    text = []

    frfs = [] # FreeRfailed values
    nfrs = [] # newFreeR values

    freeRsource = None
    OK = False
    if len(self.importxml.findall('freeRsource'))>0:
      # 'Explicit', 'Input', or 'None'
      freeRsource =  self.importxml.findall('freeRsource')[0].text

    freeRcolumnLabel = None
    if len(self.importxml.findall('freeRcolumnLabel'))>0:
      freeRcolumnLabel = self.importxml.findall('freeRcolumnLabel')[0].text

    if len(self.importxml.findall('FreeRfailed'))>0:
      freerfailed = self.importxml.findall('FreeRfailed')
      for frf in freerfailed:
        frfs.append(frf.text)

    newfreermade = False
    if len(self.importxml.findall('newFreeR'))>0:
      newfreer = self.importxml.findall('newFreeR')
      for nfr in newfreer:
        nfrs.append(nfr.text)
        if nfr.text == 'True':
          newfreermade = True

    obsfreecellcomparison = None
    if len(self.importxml.findall('ObsFreeCellComparison'))>0:
      obsfreecellcomparison = self.importxml.findall('ObsFreeCellComparison')[0]

    if len(frfs) == 2:
      # Case 5, shouldn't happen
      #print "Case 5"
      if (frfs[0] == 'True') and (frfs[0] == 'Again'):
        headlines.append("NOTE WARNING: completely failed to generate a FreeR set")
      else:
        headlines.append("SCRIPT ERROR: two values of FreeRfailed "+
                         frfs[0]+", "+frfs[1])
    elif len(frfs) == 1:
      # Case 3 or 4
      if len(nfrs) == 1:
        # Case 3, should be "False", during rerun
        if nfrs[0] == "False":
          #print "Case 3"
          freeRmessage = self.makeFreeRwarningMessage(obsfreecellcomparison)
          ##          headlines.append("WARNING: the imported FreeR set could not be extended, because all input values are the same")
          headlines.append("Generating a new FreeR set instead")
        else:
          headlines.append("SCRIPT ERROR: wrong value of newFreeR: "+nfrs[0])
      else:
        # Case 4, newFreeR should be "False", "True"
        if (nfrs[0] == "False") and (nfrs[1] == "True"):
          #print "Case 4"
          freeRmessage = self.makeFreeRwarningMessage(obsfreecellcomparison)
          for message in freeRmessage:
            headlines.append(message)
          headlines.append("A new FreeR set has been generated instead")
        else:
          headlines.append("SCRIPT ERROR: wrong values of newFreeR: "+
                           nfrs[0]+","+nfrs[1])
    else:
      # Success, case 1 or 2
      OK = True
      if len(nfrs) > 0:
        if nfrs[0] == 'False':
          print("Case 1", len(nfrs))
          message = ''
          if freeRcolumnLabel is None:
            message = "The imported FreeR set has been copied and completed"
          else:
            if freeRsource is not None and freeRsource != 'Explicit':
              message = "The imported FreeR set has been copied and completed from column " + freeRcolumnLabel
          text.append(message)
          # Was the freer set cut in resolution?
          if len(self.importxml.findall('FREERFLAGINFO/FreerCutResolution'))>0:
            cutres = self.importxml.findall('FREERFLAGINFO/FreerCutResolution')[0].text
            message = \
             "The resolution of the FreeR set was cut to match the data, "+cutres+" A"
            text.append(message)
          # But was there an incompatibility warning?
          freeRmessage = self.freerwarning(newfreermade, freeRsource, obsfreecellcomparison)
          if freeRmessage is not None:
            for message in freeRmessage:
              headlines.append(message)
            if newfreermade: headlines.append("A new FreeR set has been generated instead")
        else:
          #print("Case 2")
          #  but maybe an incompatible explicit FreeR set was rejected
          freeRmessage = self.freerwarning(newfreermade, freeRsource, obsfreecellcomparison)
          if freeRmessage is not None:
            for message in freeRmessage:
              headlines.append(message)
              OK = False
            headlines.append("A new FreeR set has been generated instead")
          else:
            text.append("A new FreeR set has been constructed")

    if freeRsource is not None:
      staraniso = False
      if len(self.importxml.findall('StarAniso'))>0:
        if str(self.importxml.findall('StarAniso')[0].text) == 'True':
          staraniso = True
      failed = ""
      if not OK:
        failed = "failed "
      # 'Explicit' or 'Input'
      if freeRsource == 'Explicit':
        text.append("The "+failed+"FreeR set was taken from explicit data object or file")
        if staraniso:
          text[-1] += ", overriding FreeR set from StarAniso"
      elif freeRsource == 'Input':
        text.append("The "+failed+"FreeR set was taken from main input data")
        if staraniso:
          text[-1] += ", and completed, even though it came from StarAniso"
      elif freeRsource == 'None':
        text.append("The "+failed+"FreeR set was taken from main input data unchanged")
        if staraniso:
          text[-1] += ", since it came from StarAniso"
      else:
        headlines.append("SCRIPT ERROR: wrong value of freeRsource: "+
                           freeRsource)  # shouldn't happen

    if len(headlines) > 0:
      # errors
      place = parent.addDiv(
        style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
      for line in headlines:
        place.append(line)

    if len(text) > 0:
      for line in text:
        #parent.append(line)
        #parent.append("<br/>")
        parent.addText(text=line)
      #parent.append("<br/>")
      if self.omitmessage is not None:
        msg = html_linebreak(self.omitmessage)
        parent.append(msg)
      #parent.append("<br/>")

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def freerwarning(self, newfreermade, freeRsource, obsfreecellcomparison):
    freeRmessage = None
    if freeRsource == 'Explicit':
      if obsfreecellcomparison is not None:
        if obsfreecellcomparison.findall('validity')[0].text == 'False':
          freeRmessage = self.makeFreeRwarningMessage(newfreermade, obsfreecellcomparison)
          return freeRmessage
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def makeFreeRwarningMessage(self, newfreermade, obsfreecellcomparison):
    ''' called if failed to extend FreeR set, returns message, cf aimless_pipe_report '''
    if obsfreecellcomparison is None:
      return ["WARNING: the imported FreeR set could not be extended, because all input values are the same"]
    else:
      # we have cell etc comparison
      freeRmessage = []

      if newfreermade:
        message = '<p>'+\
                  'WARNING: the FreeR set has not been extended,'+\
                  '  because the input FreeR set is incompatible with the observed data'+\
                  '<br/>An input FreeR set for copying or extending must match the current data in cell and Laue group'
      else:
        # cell discrepancy overridden
        message = '<p>'+\
                  'WARNING: the input FreeR set is incompatible with the observed data'+\
                  ' but acceptance was explicitly allowed<br/>'+\
                  '  BEWARE check that this is OK'+\
                  '<br/>An input FreeR set for copying or extending should match the current data in cell and Laue group'
        
      freeRmessage.append(colourText(message+'</p>', 'red'))
      #      print "** frm", type(obsfreecellcomparison)
      # some XML names were changed in Nov 2022, allow both PRE
      if len(obsfreecellcomparison.findall('cellObserved'))>0:
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

      sgname1 = obsfreecellcomparison.findall(sg1)[0].text
      sgname2 = obsfreecellcomparison.findall(sg2)[0].text
      cell1 = obsfreecellcomparison.findall(cl1)[0].text
      cell2 = obsfreecellcomparison.findall(cl2)[0].text

      obs = "Observed data: SG "+sgname1+" Cell: "+cell1
      free = "FreeR data: SG "+sg2+" Cell: "+cell2
      diff = "Cell difference: "+\
      obsfreecellcomparison.findall('CellDifference')[0].text+"&#197;"+\
      ", maximum acceptable resolution for free-R extension: "+\
      obsfreecellcomparison.findall('MaxAcceptableResolution')[0].text+"&#197;"

      validity = obsfreecellcomparison.findall('validity')[0].text
      if validity == 'True':
        message = colourText('The datasets belong to incompatible Laue groups',
        'red', fontstyle='italic', fontsize="110%")
      else:
        message = colourText('The datasets have incompatible unit cells',
                                     'red', fontstyle='italic', fontsize="110%")
      freeRmessage.append(message)
      freeRmessage.append(colourText('<br/>'+obs+'<br/>'+free+'<br/>'+diff, 'red'))

      return freeRmessage

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def finalReport(self, parent=None):
    print("finalReport")
    # Adapted from aimless_pipe_report

    if parent is None: parent = self

    self.omitmessage = None

    topDiv = parent.addDiv(style="font-size:100%;border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;")

    status = self.reportstatus()

    # what is it doing now?
    if not status == "FINISHED": return

    topDiv.addText(text="Import of merged data",
                     style='font-weight:bold; font-size:150%; color:blue;')
    #topDiv.append(' <br/>')
    self.importReport(topDiv)

    text = []
    text.append("The report below from the data reduction pipeline is purely for analysis")
    text.append("No files from this will be used in further jobs")
    for line in text:
      #topDiv.append("<br/>")
      topDiv.addText(text=line)

    #topDiv.append("<br/>")
    self.freeRmessage(topDiv)

    topDiv.append('<br/>')
    topDiv.addText(text="Report from data reduction, for analysis only",
                     style='font-weight:bold; font-size:140%; color:blue;')

    # Instantiate aimless and pointless and ctruncate reports from which to cherry pick
    # The nooutput flag makes sure thay don't do anything

    drpipelinexml = self.xmlnode.findall("AIMLESS_PIPE")[0]

    #  1) POINTLESS
    pointlessxml = drpipelinexml.findall("POINTLESS")[0]
    if pointlessxml != None:
      self.pointlessreport = \
          pointless_report.pointless_report(xmlnode=pointlessxml, jobStatus='nooutput')
      if self.fileroot is not None:
        self.pointlessreport.setFileRoot(self.fileroot) # pass on fileroot

    #  2) AIMLESS
    aimlessxml = drpipelinexml.findall("AIMLESS")
    if len(aimlessxml) == 0:
      aimlessxml = None
    else:
      aimlessxml = aimlessxml[0]
    if (aimlessxml != None):
      self.aimlessreport = \
       aimless_report.aimless_report(xmlnode=aimlessxml, jobStatus='nooutput')

    #  3) PHASER_ANALYSIS
    self.phaserreport = None
    phaserfailmessage = None
    phaserxml = drpipelinexml.findall("PHASER_ANALYSES/PHASER_ANALYSIS")
    if len(phaserxml) == 0:
      phaserxml = None
    else:
      phaserxml = phaserxml[0]
    phaserinfo = None
    table_info = None  # from Phaser
    phaserOK = False
    self.phasertNCS = None
    if phaserxml is not None:
      if len(phaserxml.findall("PhaserFailMessage"))>0:
        phaserfailmessage =  phaserxml.findall("PhaserFailMessage")[0].text
        #print(phaserfailmessage)
        phaserxml = None
        phaserOK = False  # Phaser failed
      else:
        if len(phaserxml) == 0:
          phaserxml = None
          phaserOK = False  # Phaser failed
        else:
          phaserOK = True
          self.phaserreport = phaser_analysis_report.phaser_analysis_report(
            xmlnode = phaserxml, jobStatus='nooutput')
          phaserinfo = self.getPhaserthings(phaserxml)
          self.phasertNCS = self.phaserreport.getdatavalue('tNCS')


    #  4) CTRUNCATE
    ctruncatexmlsnode = drpipelinexml.findall("CTRUNCATES")
    if len(ctruncatexmlsnode) == 0:
      ctruncatexmlsnode = None
    else:
      ctruncatexmlsnode = drpipelinexml.findall("CTRUNCATES")[0]
    ctruncatexmlnodelist = None
    self.ctruncatereports = None
    if ctruncatexmlsnode != None:
      ctruncatexmlnodelist = ctruncatexmlsnode.findall("CTRUNCATE")

    if (ctruncatexmlnodelist != None):
      self.ctruncatereports = []
      for ctruncatexmlnode in ctruncatexmlnodelist:
        self.ctruncatereports.append(
          ctruncate_report.ctruncate_report(xmlnode=ctruncatexmlnode, jobStatus='nooutput') )
        if len(self.ctruncatereports) == 0:
          self.ctruncatereports = None

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .

    # Error messages, if any
    errorspresent = ''
    if pointlessxml != None:
      if self.pointlessreport.isFatalError():
        errorspresent = errorspresent+'Fatal'

      warning = False
      if self.pointlessreport.isWarningMessage():
        errorspresent = errorspresent+'Warning'

    if len(drpipelinexml.findall('PIPELINE_ERROR'))>0:
      errorspresent += ' Pipeline Error'

      if errorspresent != '':
        errorDiv = parent.addDiv(style="border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;color:red")
        if 'Pipeline' in errorspresent:
          message = 'WARNING: Pipeline error: '+\
                    drpipelinexml.findall('PIPELINE_ERROR')[0].text
          errorDiv.append(message)

        # This takes priority over later Fatal errors
        self.pointlessreport.checkforlatticeWarning(errorDiv)

        if 'Fatal' in errorspresent:
          errorDiv.append('ERRORS')
          self.pointlessreport.Errors(errorDiv)   #  report any fatal errors

        if 'Warning' in errorspresent:
          warning = True
          self.pointlessreport.dataMissing(errorDiv)
          warningmessages = self.pointlessreport.WarningMessages()
          self.pointlessreport.addWarningMessages(errorDiv)

      else:
          self.pointlessreport.dataMissing(topDiv, colour="darkorange")

      if 'Fatal' in errorspresent:
        return
      
    warningmessage = None
    # Did Phaser run?
    if phaserOK:
      # Phaser ran successfully
      if self.IorFtype == "Amplitude":
        warningmessage = \
           "Phaser analysis was run on squared amplitudes, so may be biased"
      elif self.IorFtype == "Unknown":
        warningmessage = \
           "Phaser analysis was run on unknown data which could be amplitudes, if so may be biased"
    else:
      warn = "<br/>NOTE: analysis with Phaser was not done"+\
             " as it detected biased squared amplitudes<br/>"
      if self.IorFtype == 'Amplitude':
        warningmessage = warn
        #print("phaserfailmessage",phaserfailmessage)
        if phaserfailmessage is not None:
          if "French-Wilson" in phaserfailmessage:
            warningmessage = warn
          else:
            warn = phaserfailmessage

    if warningmessage is not None:
      warnDiv = parent.addDiv(style="color:darkorange;font-size:120%;")
      #print(warningmessage)
      warnDiv.append(warningmessage)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
    #   Overall summaries
    self.numberofdatasets = 1   # always = 1 here, I hope
    # if (aimlessxml != None):
    #    self.numberofdatasets = self.aimlessreport.numberOfDatasets()

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
    # Summary of summaries
    summaryDiv = parent.addDiv(\
           style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    summaryfold = summaryDiv.addFold(label='Key summary', brief='Headline',
                                     initiallyOpen=True)
    if pointlessxml != None:
      fail = self.pointlessreport.Errors(summaryfold)   #  report any fatal errors

    self.pointlessreport.keyTextMerged(summaryfold)

    if (aimlessxml != None):
      self.aimlessreport.keyTextMerged(summaryfold, phaserinfo)

    if phaserxml is not None:
      # getinformation content for the Table 1, in resolution bins
      # List for one dataset, list of [Overall, Inner, Outer]
      table_info = [self.infocontent(phaserxml)]
      #print("table_info",table_info)

    # ctruncate
    if self.ctruncatereports != None:
      self.ctruncatereports[0].keyText(summaryfold, phasertncs=self.phasertNCS)

    # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
    # Main summary
    overallsummaryDiv = parent.addDiv(\
         style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    overallfold = overallsummaryDiv.addFold(label='Overall summary',
                                                brief='Summary',initiallyOpen=True)

    # Put the key messages about spacegroup and resolution at the very top
    #headlineDiv = overallfold.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:5px;")
    headlineDiv = overallfold.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:5px;")
    leftDiv = headlineDiv.addDiv(style="width:49%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;border:0px;")
    statsDiv = headlineDiv.addDiv(style="width:46%;float:left;text-align:center;margin:0px; padding:10px;border:1px solid black;")
    # not usually anything useful from Pointless
    #if pointlessxml != None:
    #  self.addPointlessSummaryMerged(leftDiv)

    if (aimlessxml != None):
      # Aimless summary
      self.addAimlessSummaryMerged(statsDiv, table_info)

    #  Aimless resolution graphs in left div
    if (aimlessxml != None):
      # Put the key graphs that report on data quality at the top
      leftDiv.addText(text='Analysis as a function of resolution',style='font-size:130%;font-weight:bold;')
      leftDiv.append('<br/>')
      self.aimlessreport.ByResolutionGraphsMerged(leftDiv)
      self.aimlessreport.Graphs(leftDiv, select='Completeness')

    nextdiv = overallfold.addDiv()
    nextleftDiv = nextdiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
    nextrightDiv = nextdiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")

    if self.ctruncatereports != None:
      wilsonDiv = nextleftDiv
      wilsonDiv.append("Wilson plot")
      self.addTruncateWilson(wilsonDiv)

    # Phaser information plot
    if phaserxml is not None:
      phaserDiv = nextrightDiv
      phaserDiv.append("Information content of data")
      self.addPhaserInfoGraph(phaserxml, phaserDiv)

    if phaserfailmessage is not None:
      parent.addText(text=phaserfailmessage)
    if self.phaserreport != None:
      phaserDiv = parent.addDiv(style='width:100%; clear:both; margin:0px; padding:0px;')
      phaserfold = phaserDiv.addFold(label='Analyses of twinning, tNCS and anisotropy from Phaser',
                                            brief='Phaser',initiallyOpen=True)
      phaserDiv = phaserfold.addDiv(style="width:100%;margin:0px; padding:1px;")
      self.phaserreport.phaserAnalysisReport(phaserDiv, False)

    if self.pointlessreport is not None:
      pointlessDiv = parent.addDiv(style='width:100%; clear:both; margin:0px; padding:0px;')
      pointlessfold = pointlessDiv.addFold(\
        label='Details of space group determination',
        brief='Symmetry analysis',
        initiallyOpen=False)
      self.pointlessreport.shortreport(pointlessfold)

    # Truncate L-test next
    if self.ctruncatereports != None:
      toptruncateDiv = parent.addDiv(
        style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

      truncatefold = toptruncateDiv.addFold(label='Analysis of twinning from ctruncate',
                                    brief='Truncate',initiallyOpen=True)
      headerDiv = truncatefold.addDiv(style="width:100%; border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;text-align:center")
      headerDiv.append('<br/>')

      headerDiv.addText(text='Graphs for detecting twinning etc, more details below',
                        style='text-align:center;font-weight: bold; font-size:130%;')

      twinned = self.ctruncatereports[0].addWarningTwin(headerDiv,check=True)
      headerDiv.append('<br/>')
      if twinned:
        headerDiv.addText(text='This dataset is probably twinned',
                          style='text-align:center;font-weight: bold; font-style:italic; font-size:130%;color:red;')
      else:
        headerDiv.addText(text='This dataset is probably NOT twinned',
                            style='text-align:center;font-weight: bold; font-style:italic; font-size:130%;color:green;')

      truncateDiv = headerDiv.addDiv(style="clear:both;")
      leftDiv = truncateDiv.addDiv(style="width:48%;float:left;text-align:center;margin:0px; padding:0px;")
      rightDiv = truncateDiv.addDiv(style="width:48%;float:right;text-align:center;margin:0px; padding:0px;")

      self.addTruncateLtest(leftDiv)

      momentDiv = rightDiv.addDiv(style="width:100%;")
      for ctruncatereport in self.ctruncatereports:
        ctruncate_report.ctruncate_report(xmlnode=ctruncatereport, jobStatus='nooutput')
        ctruncatereport.acentricMoments(momentDiv)


      truncateDiv = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
      fold = truncateDiv.addFold(label="Intensity statistics: twinning tNCS etc",brief='Istats')
      self.addCtruncateReports(fold)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def mmcifreport(self):
    #  information about the mmCIF file
    text = []

    if len(self.importxml.findall("mmcifblockinfo"))>0:
      info = self.importxml.findall("mmcifblockinfo")[0].text
      info = info.splitlines()
      bname = self.importxml.findall("mmcifblock")[0].text
      text.append("Data types in cifblock "+bname+": "+info[0])
      text.append(info[1])
      for line in info:
        if 'WARNING' in line:
          text.append(line)
        
    if len(self.importxml.findall("mmcifblockdetails"))>0:
      details = self.importxml.findall("mmcifblockdetails")[0].text
      text.append(details)

    outputcolumnnames = None
    # X2MTZ block for mmcif is now replaced by MMCIF_CONVERT,
    #  but may be in old files
    if len(self.importxml.findall("X2MTZ"))>0:
      outputcolumnnames = self.importxml.findall("X2MTZ/outputcolumnnames")[0].text

    if len(self.importxml.findall('MMCIF_CONVERT'))>0:
      outputcolumnnames = self.importxml.findall("MMCIF_CONVERT/outputcolumnnames")[0].text
    if outputcolumnnames is not None:
      text.append( "Columns written to output file: "+outputcolumnnames)

    if len(self.importxml.findall("mmcifblockcolumns"))>0:
      cifcolumns = self.importxml.findall("mmcifblockcolumns")[0].text
      text.append( "mmCIF columns used: "+cifcolumns)
      
    if len(self.importxml.findall('MMCIF_CONVERT'))>0:
      outputdatatype = self.importxml.findall("MMCIF_CONVERT/datatype")[0].text
      nreflections = self.importxml.findall("MMCIF_CONVERT/nrefoutput")[0].text
      text.append("Number of reflections output: " + nreflections)
      # Resolution things
      if len(self.importxml.findall("MMCIF_CONVERT/ResolutionRange"))>0:
        rfindall = self.importxml.findall('MMCIF_CONVERT')[0]
        s = self.reportResolution(rfindall)
        text.append(s)
      cifmessage = self.importxml.findall("MMCIF_CONVERT/message")[0].text
      if cifmessage == "combine anomalous lines":
        text.append(outputdatatype+" generated from Bijvoet equivalents on different lines")
      if len(self.importxml.findall('MMCIF_CONVERT/reducemessage'))>0:
        text.append(self.importxml.findall("MMCIF_CONVERT/reducemessage")[0].text)

      self.omitmessage = None
      if len(self.importxml.findall('MMCIF_CONVERT/OmittedMissing'))>0:
        nomitted = self.importxml.findall('MMCIF_CONVERT/OmittedMissing')[0].text
        self.omitmessage = str(nomitted) + \
            ' reflections were omitted since they contained no data (all missing)'
        if len(self.importxml.findall('MMCIF_CONVERT/OmittedMissingValidFreeR'))>0:
          omitedValidFreeR = \
              self.importxml.findall('MMCIF_CONVERT/OmittedMissingValidFreeR')[0].text
          self.omitmessage += '\n  of these, '+str(omitedValidFreeR)+' had a valid FreeR flag'
        
      #text.append("Number of reflections output: "+str(nreflections))
      # FreeR things
      if len(self.importxml.findall('MMCIF_CONVERT/FreeRinFile'))>0:
        filefreer = self.importxml.findall("MMCIF_CONVERT/FreeRinFile")[0].text
        freerused = self.importxml.findall("MMCIF_CONVERT/FreeRused")[0].text
        text.append("Output FreeR generated from mmCIF column: "+freerused)
        frlist = filefreer.split()
        if len(frlist) > 1:
          s = "Note that the mmCIF file contains both "
          for frc in frlist:
            s += frc + " and "
          text.append(s[:-5])

    return text
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def reportResolution(self, rfindall):
    # Return text describing resolution range and maybe cutoff
    rrpaths = rfindall.findall('ResolutionRange')
    if len(rrpaths) == 0:
      return ''
    if len(rrpaths) > 1:
      # file resolution and cutoff resolution
      resorange = self.formatresorange(rrpaths[0])
      s = "Resolution range in file (A): " + resorange
      # Second one should be labelled "cutresolution"
      rrp = rrpaths[1].attrib
      attr = ''
      if 'id' in rrp:
        attr = rrp['id']
        if attr == "cutresolution":  # just checking
          resorange = self.formatresorange(rrpaths[1])
          s += ", trimmed to " + resorange
    else:
      # Just file resolution
      resorange = self.formatresorange(rrpaths[0])
      s = "Resolution range (A): " + resorange
      
    return s
  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def formatresorange(self, rrpath):
    dmax = None
    dmin = None
    n = 0
    s = ''
    if len(rrpath.findall('min'))>0:
      dmax = rrpath.findall('min')[0].text
      n += 1
    if len(rrpath.findall('max'))>0:
      dmin = rrpath.findall('max')[0].text
      n += 1
    if n == 2:
      # both rsent
      s = dmax + " to " + dmin
    elif dmax is None:
      s = 'high resolution '+dmin
    elif dmin is None:
      s = 'low resolution '+dmax
    return s

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def mtzreport(self):
    # for X2MTZ block from Gemmi import
    text = []
    if not len(self.importxml.findall("X2MTZ"))>0:
      return text
    print("mtzreport")
    if len(self.importxml.findall("X2MTZ/ResolutionRange"))>0:
      # X2MTZ block from Gemmi python mtzimport
      rfindall = self.importxml.findall('X2MTZ')[0]
      s = self.reportResolution(rfindall)
      text.append(s)
      nreflections = self.importxml.findall("X2MTZ/nrefoutput")[0].text
      text.append("Number of reflections output: " + nreflections)
      if len(self.importxml.findall("X2MTZ/freercolumnname"))>0:
        freercolumnname = self.importxml.findall("X2MTZ/freercolumnname")[0].text
        text.append("FreeR imported from column " + freercolumnname)
      else:
        text.append("No FreeR column in input file")
        
    return text

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def addTruncateLtest(self, parent=None):
    LtestDiv = parent.addDiv(
      style="text-align:center;margin:0px auto; padding:3px;")
    LtestDiv.addText(text='L-test for twinning',style="font-weight:bold; font-size:130%;")
    self.ctruncatereports[0].CtruncateLtest(LtestDiv)

  # . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  def addTruncateWilson(self, parent=None):
    icDiv = parent.addDiv(
      style="text-align:center;margin:0px auto; padding:3px;")
    self.ctruncatereports[0].CtruncateWilson(icDiv)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addAimlessSummaryMerged(self, parent=None, extraitems=None):
    aHeaderDiv = parent.addDiv(
      style="clear:both;font-weight:bold; font-size:130%;margin:0px;padding:0px;")
    aHeaderDiv.append('Data internal consistency statistics')

    #print("addAimlessSummaryMerged", extraitems)

    self.aimlessreport.ResultTableMerged(parent, extraitems)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addCtruncateReports(self, parent=None):
    for ctruncatereport in self.ctruncatereports:
      ctruncate_report.ctruncate_report(xmlnode=ctruncatereport, jobStatus='nooutput')
      ctruncatereport.addCtruncateReport(parent)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def getPhaserthings(self, phaserdatasetnode):
    phaserinfo = {}
    pxdname = phaserdatasetnode.attrib['name']
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
        #print("phaserreso", phaserreso)
        phaserinfo['ResolutionLimitEstimate'] = phaserreso

    return phaserinfo

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def addPhaserInfoGraph(self, phaserxml, parent=None):
    # <CCP4Table groupID="Graph" id="informationcontent" title=" Information content of data">

    icDiv = parent.addDiv(
      style="text-align:center;margin:0px auto; padding:3px;")
    graphgroup = icDiv.addFlotGraphGroup(style="width:300px;  height:270px;border:1px solid black;")
    thisgraph = phaserxml.findall('CCP4Table[@id="informationcontent"]')[0]
    graph = graphgroup.addFlotGraph(xmlnode=thisgraph, title=thisgraph.attrib["title"] )
    graph = graph.addPimpleData(xmlnode=thisgraph)

  # - - - - - - - - - - - - - - - - - - - - - - -
  def infocontent(self, phaserxml):
    # getinformation content for the Table 1, in resolution bins
    # Overall, Inner, Outer
    table_info = ['Information content bits/reflection',
                  '']
    # first we need the resolution bins used by Aimless
    #  [low,high] for each bin for each dataset
    limits = self.aimlessreport.getResolutionBins()[0]
    #print("limits",limits)

    dsetinfo = []
    for binidx in range(3):
      resrange = [float(r) for r in limits[binidx]]
      value = self.phaserreport.getInformationValues(resrange)
      dsetinfo.append("{:5.2f}".format(value))
    table_info.append(dsetinfo)

    return table_info


############################################################################
if __name__ == "__main__":

  report = import_merged_report(xmlFile = sys.argv[1],jobStatus="Finished" )
  tree= report.as_etree()
  #print etree.tostring(tree,pretty_print=True)
  report.as_html_file(fileName='./test-import.html')
  if len(report.errorReport())>0: print('ERRORS:',r.errorReport())

