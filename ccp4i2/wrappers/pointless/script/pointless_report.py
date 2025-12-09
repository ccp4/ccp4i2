from __future__ import print_function

# Normally would use findall (http://www.w3schools.com/findall/) to access the program output but
# I've added some convenience functions (haspath(),ifselect(),select() etc) to tidy up the Python code
# So here if the program output has a TwinWarning append some text to the report.  append() parses
# the text as xml - if that fails it tries adding <p> tag around text and creates a Generic() python
# object automatically and appends that.
#


import os,sys
try:
  from ccp4i2.report.CCP4ReportParser import *
except:
  exec(compile(open(os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc')).read(), os.path.join(os.environ['CCP4I2_TOP'],'bin/ccp4i2.pythonrc'), 'exec'))
  from ccp4i2.report.CCP4ReportParser import *

from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import *

# - - - - - - - - - - - - - - - - -
class pointless_report(Report):
  TASKNAME='pointless'

  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    if self.errorReport().maxSeverity()>SEVERITY_WARNING:
      print('FAILED instantiating Pointless report generator')
      self.errorReport().report()
      return

    # Testing
    #    if len(self.xmlnode.findall('POINTLESS'))>0:
    #      print "has POINTLESS"
    #    if len(self.xmlnode.findall('POINTLESS/BestSolution'))>0:
    #      print "has POINTLESS/BestSolution"
    #    if len(self.xmlnode.findall('BestSolution')>0):
    #      print "has BestSolution"

    # 'nooutput' mode would be used by another report class that wanted
    # to use some method(s) from this class for its own report
    if jobStatus is not None and jobStatus.lower() == 'nooutput':
      return

    # Comes here for report on Pointless subtask
    fail = self.Errors(self)

    projectid = jobInfo.get('projectid',None)
    jobNumber = jobInfo.get('jobnumber',None)

    if fail == False:
      self.justPointless(self,projectid=projectid,jobNumber=jobNumber)

  # - - - - - - - - - - - - - - - - -
  def justPointless(self, parent=None, extratext=None, projectid=None, jobNumber=None):
    #print "justPointless"

    summaryDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

    summaryDiv.addText(text='POINTLESS', style='font-size: 150%;text-align:center')
    if self.fileroot is not None:
      displayFile(self.fileroot, summaryDiv,
                  ['./log.txt'], 'Show log file',projectid=projectid,jobNumber=jobNumber)

    if extratext is not None:
      summaryDiv.addText(text=extratext,
                 style='font-weight:bold;font-size: 130%;text-align:center')

    self.keyText(summaryDiv)

    self.Warnings(summaryDiv)

    # do we have element scores?
    haveElements = len(self.xmlnode.findall(".//ElementScores"))>0
    
    nextDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    if haveElements:
      # we need 2 divs, main and right
      mainDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
      rightDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
    else:
      # otherwise just one
      mainDiv = nextDiv.addDiv(style="width:100%;float:left;text-align:center;margin:6px; padding:0px; ")

    #mainDiv.append("Left")
    #rightDiv.append("Right")

    self.twinWarning(mainDiv)
    self.BestSolutionType(mainDiv)
    self.BestReindex(mainDiv)
    #    self.CopyMessage(mainDiv)

    if haveElements:    
      self.ElementScoresTable(rightDiv)

    detailDiv = nextDiv.addDiv(style="width:100%;float:left;text-align:center;margin:6px; padding:0px; ")
    detailText = ""   # default use default
    if len(self.xmlnode.findall('CopyMessage'))>0:
      detailText = "Details"

    self.Details(detailDiv, usefold=False, elementscores=False, all=True,
                 text=detailText)

  # - - - - - - - - - - - - - - - - -
  def MainMessages(self,parent=None):
    """ Main messages, for different types of Pointless runs """

    #print "MainMessages.parent", parent
#    print parent.xmlnode
#    print parent.xmlnode.tag

    # not all calls here will do anything
    self.twinWarning(parent)
#    if len(self.xmlnode.findall('BestSolution'))>0:
#      print "MainMessages parent has BestSolution"
    self.BestSolutionType(parent)
    self.BestReindex(parent)
    self.CopyMessage(parent)

  # - - - - - - - - - - - - - - - - -
  def solutionWarning(self,parent=None):
    if len(self.xmlnode.findall('SolutionWarning'))>0:
      solutionwarning = html_linebreak(self.xmlnode.findall('SolutionWarning')[0].text)
      parent.append('<div style="color:red">'+solutionwarning+'</div>')

  # - - - - - - - - - - - - - - - - -
  def extraLatticeWarning(self, parent=None, addTrailer=True):
    if len(self.xmlnode.findall('ExtraLatticeCentering'))>0:
      warning = self.xmlnode.findall('ExtraLatticeCentering/Message')[0].text

      warning = colourText(warning, 'red', '110%')
      if addTrailer:
        warning += colourText(
          "(See details of possible lattice centering under 'Details of space group determination' below)",
          'red')
      warning = html_linebreak(warning,mutate=False)
      parent.append(warning)

  # - - - - - - - - - - - - - - - - -
  def extraLatticeHeadline(self, parent=None, addTrailer=True):
    # Headline version
    if len(self.xmlnode.findall('LatticeMessage'))>0:
      warning = self.xmlnode.findall('LatticeMessage')[0].text
      nl = warning[-1]
      if nl != '\n':
        warning += '\n'   #  this should be fixed in Pointless!
      lattype = self.getExtraLatticeType()
      if lattype != ' ':
        warning += 'Lattice type found: ' + lattype
      
      warning = colourText(warning, 'red', '110%')
      if addTrailer:
        warning += colourText(
          "(See details of possible lattice centering under 'Details of space group determination' below)",
          'red', '90%')

      if len(self.xmlnode.findall('StripLattice'))>0:
        warning += colourText(
          '\nCheck your indexing: if lattice is truly centred then you should re-integrate the images in the centred lattice',
          'Purple', '90%', style='; font-weight: bold')
        
      warning = html_linebreak(warning,mutate=False)
      parent.append(warning)

  
  # - - - - - - - - - - - - - - - - -
  def getExtraLatticeType(self):
    # extract found lattice type from message, if present
    lattype = ' '
    if len(self.xmlnode.findall('ExtraLatticeCentering'))>0:
      warning = self.xmlnode.findall('ExtraLatticeCentering/Message')[0].text
      fields = warning.split(' ')
      lattype = fields[fields.index('type')+1]
      
    return lattype
      
  # - - - - - - - - - - - - - - - - -
  def setFileRoot(self, fileroot):
    ''' set file root from superior script '''
    self.fileroot = fileroot
    
  # - - - - - - - - - - - - - - - - -
  def getValue(self, tag):
    ''' return value of tag if present, else None '''
    if len(self.xmlnode.findall(tag))>0:
      return self.xmlnode.findall(tag)[0].text
    return None
  # - - - - - - - - - - - - - - - - -
  def keyText(self,parent=None):
    #print "keyText"
    if len(self.xmlnode.findall('StripLattice'))>0:
      #print("***Strip lattice")
      oldLattype = self.xmlnode.findall('StripLattice/OriginalLatticeType')[0].text
      newLattype = self.xmlnode.findall('StripLattice/NewLatticeType')[0].text
      newSG = self.xmlnode.findall('StripLattice/NewSpaceGroup')[0].text
      nremoved = self.xmlnode.findall('StripLattice/Nremoved')[0].text
      nbefore = self.xmlnode.findall('StripLattice/Nbefore')[0].text
      #print(oldLattype, newLattype, newSG, nremoved, nbefore)
      s = "\nCentred lattice absences have been REMOVED for lattice type "+\
          newLattype+", original lattice type "+oldLattype
      s += "\nNumber of observation parts removed "+nremoved+" from "+nbefore
      s += "\nAfter lattice removal, space group first changed to "+newSG
      s = html_linebreak(s)
      parent.append("<br/>")
      parent.append('<span style="color:blue">'+s+'</span>')

    self.extraLatticeHeadline(parent)
      
    if len(self.xmlnode.findall("ReflectionFile[@stream='HKLIN']/MergedData"))>0:
      merged = self.xmlnode.findall("ReflectionFile[@stream='HKLIN']/MergedData")[0].text
      if merged == "True":
        s = "WARNING: Input data are merged, so space group determination is "+\
            "limited to testing for under-merging"
        parent.addText(text=s, style="color:darkorange;")
        parent.append("<br/>")
        datasetlist = self.xmlnode.findall("ReflectionData/Dataset")
        if len(datasetlist) < 2:
          s = "No scaling will be done, just analysis"
          parent.addText(text=s, style="color:darkorange;")

    if len(self.xmlnode.findall('BestSolutionType'))>0:
      solutionmessage = html_linebreak(self.xmlnode.findall("SolutionMessage")[0].text)
      parent.append(solutionmessage)
      self.solutionWarning(parent)

      stotalprob = self.xmlnode.findall("BestSolution/TotalProb")[0].text
      sconfidence = self.xmlnode.findall("BestSolution/Confidence")[0].text
      solstring = "Solution probability: %s,   Confidence %s" % (stotalprob, sconfidence)
      if len(self.xmlnode.findall('ResolutionUsed'))>0:
        testresolution = self.xmlnode.findall("ResolutionUsed/ResolutionHigh")[0].text
        solstring += "   (high resolution limit for symmetry testing %s&#197;)" % testresolution

      parent.append(solstring)

    
    if len(self.xmlnode.findall('BestReindex'))>0:
      source = ""
      if len(self.xmlnode.findall("BestReindex/HKLREF"))>0:
        source = self.xmlnode.findall("BestReindex/HKLREF")[0].text
        #print "has BestReindex/HKLREF", source
      elif len(self.xmlnode.findall("BestReindex/XYZIN"))>0:
        source = self.xmlnode.findall("BestReindex/XYZIN")[0].text

      if source is not None:
        parent.append("Reindex operator to match reference data from "+source+
                      ": "+self.xmlnode.findall("BestReindex/ReindexOperator")[0].text+
                      ", probability "+self.xmlnode.findall("BestReindex/Likelihood")[0].text)
      else:
        parent.append("<em>Alternative indexing comparisons to first file</em>")
        filenames   =  self.xmlnode.findall("BestReindex/Filename")
        reindexops  =  self.xmlnode.findall("BestReindex/ReindexOperator")
        likelihoods =  self.xmlnode.findall("BestReindex/Likelihood")
        
        table = parent.addTable(class_="center")

        numbers = []
        fnames = []
        reindexs = []
        lklihoods = []
        for i in range(len(filenames)):
          numbers.append(i+2)
          fnames.append(filenames[i].text)
          reindexs.append(reindexops[i].text)
          lklihoods.append( likelihoods[i].text)

        table.addData(title="File", data=numbers)
        table.addData(title="Filename", data=fnames)
        table.addData(title="Likelihood", data=lklihoods)
        table.addData(title="Reindex operator", data=reindexs)

      if len(self.xmlnode.findall('ResolutionUsed'))>0:
        testresolution = self.xmlnode.findall("ResolutionUsed/ResolutionHigh")[0].text
        parent.append("High resolution limit for symmetry testing %sA" % testresolution)

    else:
      # No BestReindex element, probably no alternative index to test
      hklreffile = None
      if len(self.xmlnode.findall("ReflectionFile[@stream='HKLREF']"))>0:
        hklreffile = self.xmlnode.findall("ReflectionFile[@stream='HKLREF']")[0]
      elif len(self.xmlnode.findall("ReflectionFile[@stream='XYZIN']"))>0:
        hklreffile = self.xmlnode.findall("ReflectionFile[@stream='XYZIN']")[0]
      #print "hklreffile", hklreffile
      if hklreffile is not None:
        aindx = self.xmlnode.findall("AlternativeIndexing")[0]
        # blank if no alternative indexing
        # print "ax",aindx
        s = ""
        if aindx is None:
          s = 'No alternative indexing to test relative to reference file '+\
              hklreffile.get('name')
        else:
          aindx = aindx.text
          s = 'Sole alternative indexing '+aindx +' relative to reference file\n'+\
              hklreffile.get('name')
        if s != "":
          parent.append(html_linebreak(s))

    self.CopyMessage(parent)

    self.AlternativeIndexWarning(parent)

  # - - - - - - - - - - - - - - - - -
  def keyTextMerged(self,parent=None):
    # for analysis of merged data

    if len(self.xmlnode.findall('SquaredFs'))>0:
      if self.xmlnode.findall('SquaredFs')[0].text == 'true':
        parent.append('NOTE: Input amplitudes (F) were squared to intensities (I) for analysis')

    #print 'cellpath', self.xmlnode.findall('ReflectionFile/cell')[0]
    #print 'cell', self.formatCell(self.xmlnode.findall('ReflectionFile/cell'), astext=True)[0]
    #print 'SGpath', self.xmlnode.findall('ReflectionFile/SpacegroupName')
    cell = self.xmlnode.findall('ReflectionFile/cell')[0]
    #print 'cell', type(cell), cell
    parent.append("Unit cell: "+self.formatCell(cell, astext=True))
    parent.append("Space group: "+self.xmlnode.findall('ReflectionFile/SpacegroupName')[0].text)
    
    self.extraLatticeWarning(parent)

    self.CopyMessage(parent)

    self.AlternativeIndexWarning(parent)

  # - - - - - - - - - - - - - - - - -
  def Details(self,parent=None, usefold=True, elementscores=True,
              all=False, open1=False, text="", fullversion=True):
    """    Folded details
    If fullversion False, some details omitted
    """

    if text == "":
      text = "Details of space group determination"
    brieftext = "SG details"
    if usefold:
      fold = parent.addFold(label=text,
                            brief=brieftext)
    else:
      fold = parent.addDiv(
        style="padding:10px 0px 0px 0px;")
      fold.addText(text=text , style="font-size:120%;")

    self.IndexScores(fold, open1=open1)
    self.ExtralatticeScores(fold)
    if elementscores: self.ElementScores(fold, open1=open1)
    self.LaueGroupScores(fold, open1=open1)
    self.ZoneScores(fold, open1=open1)
    self.SpaceGroupScores(fold, open1=open1)
    if fullversion: self.TestDataDetails(fold, open1=open1)
    if all:
      self.AdditionalGraphs(fold)

  # - - - - - - - - - - - - - - - - -
  def shortreport(self, parent=None, fromImportMerged=True):
    #print("shortreport", self.fileroot)

    projectid = self.jobInfo.get('projectid',None)
    jobNumber = self.jobInfo.get('jobnumber',None)

    if self.fileroot is not None and projectid is not None and jobNumber is not None:
      displayFile(self.fileroot, parent,
    [ 'job_3/job_1/log.txt', 'job_1/log.txt', './log.txt'], 'Show log file', projectid=projectid, jobNumber=jobNumber)

    # do we have element scores?
    haveElements = len(self.xmlnode.findall(".//ElementScores"))>0
    
    nextDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
    if haveElements:
      # we need 2 divs, main and right
      mainDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
      rightDiv = nextDiv.addDiv(style="width:48%;float:left;text-align:center;margin:6px; padding:0px; ")
    else:
      # otherwise just one
      mainDiv = nextDiv.addDiv(style="width:100%;float:left;text-align:center;margin:6px; padding:0px; ")

    self.BestSolutionType(mainDiv)
    self.BestReindex(mainDiv)

    if haveElements:    
      self.ElementScoresTable(rightDiv)

    detailDiv = nextDiv.addDiv(style="width:100%;float:left;text-align:center;margin:6px; padding:0px; ")
    detailText = ""   # default use default
    if len(self.xmlnode.findall('CopyMessage'))>0:
      detailText = "Details"

    self.Details(detailDiv, usefold=False, elementscores=False, all=False,
                 open1=True, text=detailText, fullversion=False)

  # - - - - - - - - - - - - - - - - -
  def isFatalError(self):
    # True if there is a FatalErrorMessage
    if len(self.xmlnode.findall('FatalErrorMessage'))>0:
      return True
    return False

  # - - - - - - - - - - - - - - - - -
  def isWarningMessage(self):
    # True if there is a warningmessage
    if len(self.xmlnode.findall(".//*[@class='warningmessage']"))>0:
      #      print ("isWarningMessage:")
      return True
    return False
  # - - - - - - - - - - - - - - - - -
  def WarningMessages(self):
    # return list of warning messages
    return self.xmlnode.findall(".//*[@class='warningmessage']")
  # - - - - - - - - - - - - - - - - -
  def isMerged(self):
    # Returns true if all HKLIN files are merged
    merged = False
    mergeddata = self.xmlnode.findall("ReflectionFile[@stream='HKLIN']/MergedData")
    if len(mergeddata) > 0:
      merged = True
      for md in mergeddata:
        if md.text != "True":
          merged = False
    self.merged = merged
    return merged
  # - - - - - - - - - - - - - - - - -
  def dataMissing(self,parent=None, colour='red'):
    # Process DataMissing block
    if len(self.xmlnode.findall('DataMissing'))>0:
      # DataMissing block found
      datamissing = self.xmlnode.findall('DataMissing')[0]
      att = datamissing.attrib
      if (len(att) > 0 and 'warningmessage' in att['class']):
        # severe warning message
        message = "<br/>SEVERE WARNING: "+datamissing.text
        if colour != '':
          message = '<div style="color:'+colour+';font-size:130%">'+message+'</div>'
        parent.append(message)
        parent.append("<br/>")
      else:
        # just some missing data
        message = "<br/>NB: "+datamissing.text+"<br/>"
        # I+, I- counts
        nIplus = self.xmlnode.findall('N_Iplus')[0].text
        nIminus = self.xmlnode.findall('N_Iminus')[0].text
        if nIplus != '' and nIminus != '':
          difference = abs(int(nIplus) - int(nIminus))
          message = message + \
            "Counts:  number I+ = "+nIplus+", number I- = "+nIminus+\
            ", difference = "+str(difference)+"<br/>"
        if colour != '':
          message = '<div style="color:'+colour+';font-size:130%">'+message+'</div>'
        parent.append(message)
        parent.append("<br/>")
      
  # - - - - - - - - - - - - - - - - -
  def isTwinWarning(self):
    if len(self.xmlnode.findall('TwinWarning'))>0:
      return True
    return False
      
  # - - - - - - - - - - - - - - - - -
  def Errors(self,parent=None, colour=True):
    fail = False
    if len(self.xmlnode.findall('FatalErrorMessage'))>0:
      message = "FATAL ERROR: "+\
                self.xmlnode.findall('FatalErrorMessage')[0].text
      message = html_linebreak(message)
      if colour:
        message = '<div style="color:red">'+message+'</div>'
      parent.append(message)
      fail = True
    if len(self.xmlnode.findall('FatalErrorMessage2'))>0:
      message = self.xmlnode.findall('FatalErrorMessage2')[0].text
      message = html_linebreak(message)
      if colour:
        message = '<div style="color:red">'+message+'</div>'
      parent.append(message)
      fail = True
    return fail
  # - - - - - - - - - - - - - - - - -
  def addWarningMessages(self,parent=None, colour=True):
    # Add warning messages to output
    if len(self.xmlnode.findall(".//*[@class='warningmessage']"))>0:
      warnings =  self.xmlnode.findall(".//*[@class='warningmessage']")
      for warning in warnings:
        # print("*!* WARNING:", warning.text)
        message = warning.text
        message = html_linebreak(message)
        if colour:
          message = '<div style="color:darkorange">'+message+'</div>'
        parent.append(message)
  # - - - - - - - - - - - - - - - - -
  def checkforlatticeWarning(self, parent):
    # Check for lattice discrepancy warning
    # eg WARNING: Lattice type P deduced from reflection data is different
    # from that in the input symmetry H
    # Special option added to deal with a weird file from EJD

    pattern = "deduced from reflection data is different from that in the input symmetry"
    if len(self.xmlnode.findall("*[@class='warningmessage']"))>0:
      warnings =  self.xmlnode.findall("*[@class='warningmessage']")
      for warning in warnings:
        if pattern in warning.text:
          # Make a better error message
          words = warning.text.split()
          lat_deduced = words[3]
          lat_input = words[16]
          # print( lat_deduced, lat_input)
          message = "SERIOUS PROBLEM: apparent lattice type "+lat_deduced+\
                    " does not match type specified in the input file "+\
                    lat_input
          message += "\nYou need to resolve this discrepancy in the input file\n"
          message = html_linebreak(message)
          message = '<div style="color:red">'+message+'</div>'
          parent.append(message)


  # - - - - - - - - - - - - - - - - -
  def twinWarning(self,parent=None):
    if len(self.xmlnode.findall('TwinWarning'))>0:
      parent.append('<div style="color:darkorange">'+
                    self.xmlnode.findall('TwinWarning')[0].text
                    + '<br/>\n' + 
                   'Rough estimated twin fraction:'+
                    self.xmlnode.findall("TwinFraction")[0].text+'</div>' )

  # - - - - - - - - - - - - - - - - -
  def AlternativeIndexWarning(self, parent=None):
    # add warning message if the output data has an alternative indexing and there if no reference
    if len(self.xmlnode.findall('NumberPossibleReindexing'))>0:
      nreindex = self.xmlnode.findall('NumberPossibleReindexing')[0].text
      if not (len(self.xmlnode.findall("ReflectionFile[@stream='HKLREF']"))>0) or \
             (len(self.xmlnode.findall("ReflectionFile[@stream='XYZIN']"))>0):

        alternativelist = self.xmlnode.findall('PossibleReindexing')
        allcellssame = True
        reindexingoperators = ""
        reindexingoperatorsanddiffs = ""
        notfirst = False
        
        for alternative in alternativelist:
          operator = alternative.findall('ReindexOperator')[0].text
          celldiff = alternative.findall('CellDiff')[0].text
          differentcell = (alternative.findall('DifferentCell')[0].text == 'true')
          #print 'alternatives ', operator, celldiff, differentcell, alternative.findall('DifferentCell')
          if celldiff != '0.00':
            allcellssame = False
          if notfirst:
            reindexingoperators += ", "
            reindexingoperatorsanddiffs += ", "
          notfirst = True
          reindexingoperators += operator
          reindexingoperatorsanddiffs += operator+" ("+celldiff+")"
          
        if allcellssame:
          message1 = "NOTE: the final selected symmetry has alternative indexing schemes, but no reference data has been given"
          reindexingoperators = "Possible alternative indexing operators: " + reindexingoperators
        else:
          message1 = "NOTE: the final selected symmetry and cell have alternative indexing schemes, but no reference data has been given"
          reindexingoperators = "Possible alternative indexing operators (with cell differences in A): "+reindexingoperatorsanddiffs

        message2 = "If you already have a matching dataset, you should choose it as a reference set to get consistent indexing"
        parent.append('<div style="color:blue"><br/>'+message1+'<br/>'+reindexingoperators+'<br/>'+message2+'</div>')
    
  # - - - - - - - - - - - - - - - - -
  def BestSolutionType(self,parent=None):
    # Dependent on a program output parameter create a Table object and call addData() to load the table
    # Here the append() method is passed the table object.
    
    if len(self.xmlnode.findall('BestSolutionType'))>0:
      solutionmessage = html_linebreak(self.xmlnode.findall("SolutionMessage")[0].text)
      parent.append(solutionmessage)
      self.solutionWarning(parent)

      parent.append('<br/> Solution type: '+
                    self.xmlnode.findall("BestSolutionType")[0].text )
      table = parent.addTable( select = ".//BestSolution", transpose=True, style="margin: 0 auto;border:1px solid orange;", class_="center" )
      for title,select in [ [ "Group name" , "GroupName" ],
                            [ "Reindex" , "ReindexOperator" ],
                            [ "Space group confidence" , "Confidence" ],
                            [ "Laue group confidence" , "LGconfidence" ],
                            [ "Laue group probability" , "LGProb" ],
                               [ "Systematic absence probability" , "SysAbsProb" ] ] :
        table.addData( title=title , select = select )
  # - - - - - - - - - - - - - - - - -
  def BestReindex(self,parent=None):
    # -- Indexing relative to reference file --
    if len(self.xmlnode.findall('.//BestReindex'))>0:
      #  XYZIN reference
      if len(self.xmlnode.findall(".//XYZIN"))>0:
        reffile, refSG = self.refSpaceGroup('XYZIN')
        if refSG is not None:
          parent.append('Reference reflection list generated from coordinate file'+ reffile+'<br/>' + \
                    'Space group '+ refSG )
      else:
        # HKLREF reference
        reffile, refSG = self.refSpaceGroup('HKLREF')
        if refSG is not None:
          parent.append('Determining best alternative indexing relative to reference file: ' + 
                        '<br/> ' + reffile + '<br/>' + 
                        'Space group '+ refSG )
        else:
            parent.append('Determining best alternative indexing relative to first input file: ')

      tablenode = None
      if len(self.xmlnode.findall("BestReindex"))>0:
          tablenode = self.xmlnode.findall("BestReindex")[0]
      if tablenode is not None:
        table = parent.addTable(select="BestReindex", class_="center")
        for title,select in [ [ "File number", "Filenumber" ],
                              [ "Reindex operator" , "ReindexOperator" ],
                              [ "Confidence" , "Confidence" ],
                              [ "Likelihood" , "Likelihood" ],
                              [ "CC" , "CC" ]] :
          table.addData(tablenode , title=title , select = select )

  # - - - - - - - - - - - - - - - - -
  def AxialGraphs(self,parent=None):
    # Select only axial graphs
    # Add a Graph - add a table of data and plot instructions to the graph

    # Loop over all Graph tables in the program output and add to the GraphGroup
    # The plotting instructions are provided as xml text

    self.Graphs(parent, select="AxialReflections")

  # - - - - - - - - - - - - - - - - -
  def Graphs(self,parent=None, select=None):
    # Select only specified graphs
    # Add a Graph - add a table of data and plot instructions to the graph
    # Loop over all Graph tables in the program output and add to the GraphGroup
    # The plotting instructions are provided as xml text
    if select is None: return

    graphXmlnodeList = self.xmlnode.findall("CCP4Table[@groupID='Graph']")

    for graphXmlnode in graphXmlnodeList:
      idt = graphXmlnode.findall("[@id]")[0].attrib["id"]
      if idt.find(select) >= 0:
        graph = parent.addFlotGraph( xmlnode=graphXmlnode, title=graphXmlnode.findall("[@title]")[0].attrib["title"], class_="center" )
        graph.addTable( select="data", headers = 'headers' )
        graph.addPlot(  select= "plot" )

  # - - - - - - - - - - - - - - - - -
  def CopyMessage(self,parent=None):
    if len(self.xmlnode.findall('CopyMessage'))>0:
      s = self.xmlnode.findall('CopyMessage')[0].text
      if len(self.xmlnode.findall('ReindexOperator'))>0:
        s += "  with reindex operator "+self.xmlnode.findall('ReindexOperator')[0].text
      parent.append(html_linebreak(s + '\n'))

  # - - - - - - - - - - - - - - - - -
  def Warnings(self, parent=None):
    colour = "chocolate"
    if len(self.xmlnode.findall('ReindexWarning'))>0:
      s ="WARNING: " + self.xmlnode.findall('ReindexWarning')[0].text
      s = html_linebreak(s)
      parent.append('<div style="color:'+colour+'">'+s+'</div>')

    printSGchange = True
    if self.getValue('SpacegroupChanged') == "False":
      printSGchange = False
    
    if printSGchange:
      sgn = self.getValue('NewSpacegroupName')
      if sgn is not None:
        sg1 = self.getValue('ReflectionFile/SpacegroupName')
        if sg1 is not None:
          s = "\nWARNING: Space group changed from input "+sg1+" to new group "+sgn
          s = html_linebreak(s)
          parent.append('<div style="color:'+colour+'">'+s+'</div>')

    if len(self.xmlnode.findall('ChangedSpaceGroupMessage'))>0:
      s = "\nWARNING: "+self.getValue('ChangedSpaceGroupMessage')
      s = html_linebreak(s)
      parent.append('<div style="color:'+colour+'">'+s+'</div>')

    if len(self.xmlnode.findall('ReindexDeterminantMessage'))>0:
      s = "\nWARNING: "+self.getValue('ReindexDeterminantMessage')
      #print("::",self.getValue('ReindexDeterminantMessage'))

      ndiscarded = self.getValue('FractionalIndicesDiscarded')
      if ndiscarded is not None:
        ntot = self.getValue('ReflectionData/NumberReflections')
        if ntot is not None:
          s += "\nNumber of reflections discarded = "+ndiscarded+\
               ", from the total number "+ntot
      s = html_linebreak(s)
      parent.append('<div style="color:'+colour+'">'+s+'</div>')


  # - - - - - - - - - - - - - - - - -
  def TestDataDetails(self,parent=None, open1=False):

    fold = parent.addFold(label="Details of test data (in Pointless)", brief='TestData',
                          initiallyOpen=open1)
    merged = "False"
    # data is merged if "ReflectionData/NumberObservations" is absent,
    #  or if "ReflectionData/Merged" = True
    if len(self.xmlnode.findall(".//MergedData"))>0:
      merged = self.xmlnode.findall(".//MergedData")

    if not (len(self.xmlnode.findall("ReflectionData/NumberObservations"))>0):
      merged = "True"

    if merged == "False":
      fold.append("Summary of unmerged test reflection data")

      table = fold.addTable( select = ".//POINTLESS/ReflectionData", transpose=False, class_="center" )
      for title,select in [ [ "Max resolution", "ResolutionHigh" ],
                            [ "Nreflections", "NumberReflections" ],
                            [ "NObservations", "NumberObservations" ],
                            [ "Nparts", "NumberParts" ],
                            [ "Nbatches", "NumberBatches" ],
                            [ "Ndatasets", "NumberDatasets" ] ]:
        table.addData( title=title , select = select )

      # Dataset list
      table2 = fold.addTable(class_="center")
      datasetlist = self.xmlnode.findall("ReflectionData/Dataset")

      #  Data lists for each title
      alldatalist = {}
      titlelist = []
      for title in ["DatasetName", "RunNumber", "Batch range",
                    "Phi range", "Excluded batches", "Wavelength"]:
        titlelist.append(title)
        alldatalist[title] = []   # dictionary of data lists
        
      for dataset in datasetlist:
        wavelength = dataset.findall("Wavelength")[0].text
        # list of run data
        runlist = dataset.findall("Run")
        # List of data wanted
        for title,select in [["DatasetName", "Datasetname"],
                             ["RunNumber", "number"],
                             ["Batch range", "BatchRange"],
                             ["Phi range", "PhiRange"],
                             ["Excluded batches", "ExcludedBatches"]]:
          for run in runlist:  # loop runs
            if select == 'BatchRange' or select == 'PhiRange':
              item = formatRange(run.findall(select)[0].text)
              #print "Range", item
            else:
              item = run.findall(select)[0].text
            alldatalist[title].append(item)

        # Add in Wavelength for each run
        alldatalist["Wavelength"] += [wavelength] * len(runlist)

      for title in titlelist:
        table2.addData( title=title , data=alldatalist[title])
        

    if (merged == "True") or len(self.xmlnode.findall("ReflectionData/MergedData"))>0:
      fold.append("Summary of merged test reflection data")
      table3 = fold.addTable(select="ReflectionData",class_="center")
      for title,select in [ [ "Max resolution", "ResolutionHigh" ],
                            [ "Nreflections", "NumberReflections" ]]:
        table3.addData( title=title , select = select )

  # - - - - - - - - - - - - - - - - -
  def ExtralatticeScores(self, parent=None):
    if len(self.xmlnode.findall('ExtraLatticeCentering'))>0:
      fold = parent.addFold(label="Warning of possible lattice centering",
                            brief='ExtraLatticeCentering',
                            initiallyOpen=True)
      self.extraLatticeWarning(fold, addTrailer=False)
      meanE2absent = self.xmlnode.findall("ExtraLatticeCentering/MeanE2absent")[0].text
      possiblelattice = self.xmlnode.findall("ExtraLatticeCentering/LatticeType")[0].text
      #print "possiblelattice, meanE2absent", possiblelattice, meanE2absent
      text = "Mean(E^2) for reflections which would be missing in a lattice of type "+\
             possiblelattice
      text += " is "+ meanE2absent + " instead of the ideal value 1.0"
      fold.append(text)
      self.Graphs(fold, select="LatticeCentering")

  # - - - - - - - - - - - - - - - - -
  def ElementScores(self,parent=None, open1=False):
    if len(self.xmlnode.findall(".//BestSolutionType"))>0:
      fold = parent.addFold(label="Scores for each symmetry element",
                            brief='SymmetryElements',
                            initiallyOpen=True)
      fold.append("Lattice group name "+\
                  self.xmlnode.findall("LatticeSymmetry/LatticegroupName")[0].text)
      if len(self.xmlnode.findall("LatticeSymmetry/LatticeReindex"))>0:
        fold.append("Reindex operator from input to lattice: "+\
                    self.xmlnode.findall("LatticeSymmetry/LatticeReindex")[0].text)

      table = fold.addTable(select=".//ElementScores",class_="center")
      for title,select in [["Likelihood", "Element/Likelihood"],
                           ["CC", "Element/CC"],
                           ["Z-CC", "Element/ZCC"],
                           ["N", "Element/NCC"],
                           ["R", "Element/R"],
                           ["", "Element/ElementScoreStars"],
                           ["Symmetry", "Element/SymmetryElementString"]]:
        table.addData( title=title , select = select )

  # - - - - - - - - - - - - - - - - -
  def shortsym(self,symstring):
    if symstring.count('{') < 2:
        # only one operator, leave unchanged
        return symstring
    # otherwise remove the second opreator
    f = symstring.rsplit('{')
    return f[0]+'{'+f[1]
  # - - - - - - - - - - - - - - - - -
  def ElementScoresTable(self,parent=None):
    # short form for embedding in other div

    parent.addText(text="Scores for each symmetry element", style="font-size:120%;")
    if len(self.xmlnode.findall(".//BestSolutionType"))>0:
      parent.append("Lattice group name "+\
                  self.xmlnode.findall("LatticeSymmetry/LatticegroupName")[0].text)
      if len(self.xmlnode.findall("LatticeSymmetry/LatticeReindex"))>0:
        parent.append("Reindex operator from input to lattice: "+\
                    self.xmlnode.findall("LatticeSymmetry/LatticeReindex")[0].text)

    xlikelihoods = self.xmlnode.findall(".//ElementScores/Element/Likelihood")
    xCCs = self.xmlnode.findall(".//ElementScores/Element/CC")
    xRs = self.xmlnode.findall(".//ElementScores/Element/R")
    xstars = self.xmlnode.findall(".//ElementScores/Element/ElementScoreStars")
    xsymmetries = self.xmlnode.findall(".//ElementScores/Element/SymmetryElementString")

    likelihoods = []
    CCs = []
    Rs = []
    stars = []
    symmetries = []
    for i in range(len(xCCs)):
      likelihoods.append(xlikelihoods[i].text)
      CCs.append(xCCs[i].text)
      Rs.append(xRs[i].text)
      stars.append(xstars[i].text)
      # Abbreviate symmetry element to fit
      symmetries.append(self.shortsym(xsymmetries[i].text))

    #print "syms", symmetries
    table = parent.addTable(class_="center")
    table.addData(title="Likelihood", data=likelihoods)
    table.addData(title="CC", data=CCs)
    table.addData(title="R", data=Rs)
    table.addData(title="", data=stars)
    table.addData(title="Symmetry", data=symmetries)

  # - - - - - - - - - - - - - - - - -
  def LaueGroupScores(self,parent=None, open1=False):
    if len(self.xmlnode.findall(".//LaueGroupScoreList"))>0:
      fold = parent.addFold(label="Scores for each Laue group", brief='Laue group',
                            initiallyOpen=open1)
      note = 'Accept flags: ">", accepted; "=" original accepted; "-" original rejected<br/>'+\
             'Dcell: maximum cell difference in degrees'
      fold.append(note)
      table = fold.addTable(select=".//LaueGroupScoreList",class_="center")
      for title,select in [["Accept", "LaueGroupScore/LaueGroupScoreString"],
                           ["Group", "LaueGroupScore/LaueGroupName"],
                           ["", "LaueGroupScore/LaueGroupScoreStars"],
                           ["Lklhd", "LaueGroupScore/Likelihood"],
                           ["NetZ CC", "LaueGroupScore/NetZCC"],
                           ["ZCC+", "LaueGroupScore/ZCC_plus"],
                           ["ZCC-", "LaueGroupScore/ZCC_minus"],
                           ["CC", "LaueGroupScore/CC"],
                           ["R", "LaueGroupScore/R"],
                           ["Dcell", "LaueGroupScore/CellDelta"],
                           ["Reindex", "LaueGroupScore/ReindexOperator"]]:
        table.addData( title=title , select = select )

  # - - - - - - - - - - - - - - - - -
  def ZoneScores(self,parent=None, open1=False):
    if len(self.xmlnode.findall(".//ZoneScoreList"))>0:
      fold = parent.addFold(label="Systematic absence analysis for each zone",
                            brief='Absences', initiallyOpen=open1)

      zonelauegroups = self.xmlnode.findall("ZoneScoreList/ZoneLaueGroup")

      for zlg in zonelauegroups:  # loop Laue groups (usually only 1)
        fold.append("Systematic absence analysis for each zone, Laue group "+\
                    zlg.findall("[@name]")[0].attrib["name"])

        # list of zones in group
        zonenodelist = zlg.findall("Zone")

        stars = []
        missing = False
        empty = False

        for znode in zonenodelist:          # loop zones
          nobs = znode.findall("Nobs")[0].text

          if len(znode.findall("MissingData"))>0:
            # Case 1, systematic missing data
            stars.append("Missing")
            missing = True
          elif (int(nobs) == 0):
            # Case 2, no observations
            stars.append("Empty")
            empty = True
          else:
            # Case 3, data found and used
            stars.append(znode.findall("ZoneScoreStars")[0].text)

          # end loop zones

        message = ""
        if empty:
          message += "'Empty' means no observations for this zone<br/>"
        if missing:
          message += \
            "'Missing' means data for this zone are systematically absent,"+\
            " may have been removed previously"
        if message != "":
          fold.append(message)

        # table for zones in this group
        table = fold.addTable(xmlnode=zlg,class_="center")
        table.addData(title="Zone", select="Zone/ZoneType")
        table.addData(title="PeakHeight", select="Zone/PeakHeight")
        table.addData(title="SD", select="Zone/SDPkHt")
        table.addData(title="", data=stars)
        table.addData(title="Probability", select="Zone/Prob")
        table.addData(title="Nobs", select="Zone/Nobs")
        table.addData(title="ReflectionCondition", select="Zone/Condition")

      # end loop Laue groups

      self.AxialGraphs(fold)

  # - - - - - - - - - - - - - - - - -
  def SpaceGroupScores(self,parent=None, open1=False):
    if len(self.xmlnode.findall(".//SolutionMessage"))>0:
      fold = parent.addFold(label="Scores for possible space groups",brief='SpaceGroupScores',
                            initiallyOpen=open1)

      fold.append("<p>"+self.xmlnode.findall("SolutionMessage")[0].text+"</p>")

      msg = " 'Reindex' is the operator to convert from the input hklin frame to the standard spacegroup frame.<br/>"+\
            " 'TotalProb' is a total probability estimate (unnormalised)<br/>"+\
            " 'SysAbsProb' is an estimate of the probability of the space group from observed systematic absences.<br/>"+\
            " 'Conditions' are the reflection conditions (absences)"
      fold.append(msg)


      table = fold.addTable(select=".//SpacegroupList",class_="center")
      for title,select in [["Spacegroup", "Spacegroup/SpacegroupName"],
                           ["ITnumber", "Spacegroup/SGnumber"],
                           ["TotalProb", "Spacegroup/TotalProb"],
                           ["SysAbsProb", "Spacegroup/SysAbsProb"],
                           ["Reindex", "Spacegroup/ReindexOperator"],
                           ["Conditions", "Spacegroup/Condition"]]:
        table.addData( title=title , select = select )

  # - - - - - - - - - - - - - - - - -
  def refSpaceGroup(self, source):
    #  source  ==  XYZIN or HKLREF
    reffile = self.xmlnode.findall(".//BestReindex/"+source)[0].text
    refSG = self.xmlnode.findall(".//BestReindex/RefSpaceGroup")[0].text
    return reffile, refSG    # may be None
  # - - - - - - - - - - - - - - - - -
  def IndexScores(self,parent=None, open1=False):
    if len(self.xmlnode.findall(".//IndexScores"))>0:
      fold = parent.addFold(label="Alternative index scores",brief='IndexScores',
                            initiallyOpen=open1)
      if len(self.xmlnode.findall("XYZREF"))>0:
        fold.append('Reference reflection list generated from coordinate file'+ \
                    self.xmlnode.findall("XYZREF")[0].text+'<br/>' + \
                    'Space group '+self.xmlnode.findall("XYZREFspacegroup")[0].text )

      else:
        reffile, refSG = self.refSpaceGroup('HKLREF')
        if refSG is not None:
            fold.append('Determining best alternative indexing relative to reference file: ' + \
                        reffile + '<br/>' + \
                        'Space group '+ refSG  )
        else:
            fold.append('Determining best alternative indexing relative to first input file: ')            

      if len(self.xmlnode.findall(".//AlternativeIndexing"))>0:
          fold.append('Possible reindex operators: '+ self.xmlnode.findall(".//AlternativeIndexing")[0].text )

      # there may be multiple IndexScores blocks
      tablenodes = self.xmlnode.findall("IndexScores")  # list of blocks
      if tablenodes is not None:
        for tablenode in tablenodes:
          table = fold.addTable(select="IndexScores",class_="center")
          for title,select in [ [ "Reindex operator" , "Index/ReindexOperator" ],
                                [  "Likelihood" , "Index/Likelihood" ],
                                [ "CC" , "Index/CC" ],
                                [ "CellDeviation", "Index/CellDeviation"]] :
            table.addData(tablenode , title=title , select = select )


  # - - - - - - - - - - - - - - - - -
  def AdditionalGraphs(self,parent=None):
    """ Resolution and twinning graphs"""
    id1 = "Resolution estimate from Pointless"
    id2 = "L-test"
    if (len(self.xmlnode.findall(".//CCP4Table[@id='"+id1+"']"))>0) or \
       (len(self.xmlnode.findall(".//CCP4Table[@id='"+id2+"']"))>0):
      fold = parent.addFold(label="Resolution analysis for score cutoff and twinning",
                            brief='Cutoff/Twin',
                            initiallyOpen=True)

      leftDiv = fold.addDiv(style="width:47%;float:left;text-align:left;margin:0px; padding:0px; line-height:100%; font-size:100%;")
      rightDiv = fold.addDiv(style="width:48%;float:left;text-align:left;margin:0px; padding:10px;")
      self.ResolutionGraphs(leftDiv, idt=id1)
      self.TwinningGraphs(rightDiv, idt=id2)
    
  # - - - - - - - - - - - - - - - - -
  def ResolutionGraphs(self,parent=None, idt=None):
    if idt is None: return
    if len(self.xmlnode.findall(".//CCP4Table[@id='"+idt+"']"))>0:
      self.Graphs(parent, select=idt)

    if len(self.xmlnode.findall('CutoffMessage'))>0:
      s = self.xmlnode.findall('CutoffMessage')[0].text + "\n"
      s += "Highest resolution used : "+\
           self.xmlnode.findall('ResolutionUsed/ResolutionHigh')[0].text+"A\n"
      s += " Estimate from I/sigI  > " + self.xmlnode.findall('IovSigCutoff')[0].text+\
           " : " + self.xmlnode.findall('ResolutionHighIovsigI')[0].text + "A\n"
      if len(self.xmlnode.findall('ResolutionHighCCHalf'))>0 and len(self.xmlnode.findall('CCHalfCutoff'))>0:
          s += " Estimate from CC(1/2) > " + self.xmlnode.findall('CCHalfCutoff')[0].text+\
               " : " + self.xmlnode.findall('ResolutionHighCCHalf')[0].text + "A\n"
      parent.append(html_linebreak(s))
  # - - - - - - - - - - - - - - - - -
  def TwinningGraphs(self,parent=None, idt=None):
    if idt is None: return
    self.Graphs(parent, select=idt)

    if len(self.xmlnode.findall('TwinWarning'))>0:
      s = self.xmlnode.findall('TwinWarning')[0].text+"\n"+\
          " Estimated fraction: "+self.xmlnode.findall('TwinFraction')[0].text
      parent.append(html_linebreak(s))
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


############################################################################
if __name__ == "__main__":
#  report = PointlessReport(xmlFile = os.path.join(os.environ['CCP4I2_TOP'],'test','report_test','gam_1.xml' ))

  report = pointless_report(xmlFile = sys.argv[1] )
  tree= report.as_etree()
  #  print etree.tostring(tree,pretty_print=True)
  report.as_html_file(fileName='./test.html')
