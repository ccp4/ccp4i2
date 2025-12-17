
import os
from ccp4i2.report.CCP4ReportParser import *

# - - - - - - - - - - - - - - - - -
def displayFile(fileroot, parent, filenames, text, projectid=None, jobNumber=None):
  ''' display message with link to open file
  fileroot   root of file names, None if not set
  parent     where to put it
  filenames  list of possible names
  text       to display
  '''
  p = GenericElement('p')
  filefound = False
  fnames = filenames

  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print(fileroot, parent, filenames, text, projectid, jobNumber)
  print("%%%%%%%")
  
  if fileroot == None:
      fileroot = ""
  
  if fileroot != None:
    # prpend fileroot on list
    #print 'fileroot', fileroot
    fnames = []
    for filename in filenames:
      os.path.isfile(fileroot+filename)
      filefound = True
      break

  if filefound:
    if len(filename.split('/')) > 1:
        subJobNumber = filename.split('/')[0][4:]
        try:
          if "." in jobNumber:
            jobNumber = jobNumber.split(".")[0]
          href = "/database/?getProjectJobFile?projectId="+projectid+"?fileName="+filename.split('/')[1]+"?jobNumber="+jobNumber+"?subJobNumber="+subJobNumber
        except Exception as err:
          href = "about:blank"
    else:
        href = "/database/?getProjectJobFile?projectId="+projectid+"?fileName="+filename+"?jobNumber="+jobnumber
    p.append(GenericElement('a',text,href=href))
    parent.append(p)

# - - - - - - - - - - - - - - - - -
def displayFileList(fileroot, parent, items, box=True, projectid=None, jobNumber=None):
  ''' display message with link to open file
  fileroot   root of file names, None if not set
  parent     where to put it
  items      list of [filename, text]
  '''

  fileDiv = parent.addDiv(style="width:80%;display:flex")

  if fileroot == None:
      fileroot = ""

  for item in items:
    filename = item[0]
    text = item[1]

    if os.path.isfile(fileroot+filename):
      nextDiv = fileDiv.addDiv(
        style="width:25%;border: 2px solid blue; text-align:center; margin:3px; padding:6px;display:inline-block;")
      p = GenericElement('p')
      if len(filename.split('/')) > 1:
          subJobNumber = filename.split('/')[0][4:]
          href = "/database/?getProjectJobFile?projectId="+projectid+"?fileName="+filename.split('/')[1]+"?jobNumber="+jobNumber+"?subJobNumber="+subJobNumber
      else:
          href = "/database/?getProjectJobFile?projectId="+projectid+"?fileName="+filename+"?jobNumber="+jobnumber
      p.append(GenericElement('a',text,href=href))
      nextDiv.append(p)

# - - - - - - - - - - - - - - - - -
def colourText(text, colour, fontsize=None, fontstyle=None, style=None):
  fsz = ""
  if fontsize is not None:
    fsz += '; font-size: '+fontsize
  if fontstyle is not None:
    fsz += '; font-style: '+fontstyle
  if style is not None:
    fsz += "; "+style
  return '<span style="color:'+colour+fsz+'">'+text+'</span>'

# - - - - - - - - - - - - - - - - -
def html_linebreak(line, mutate=True):
  # delineate lines by "<br/>"
  # mutate characters '&', '<', '>' into safe equivalents

  pline = ''
  for c in line:
    if mutate:
      if c == '&': pline += '&amp;'
      elif c == '<': pline += '&lt;'
      elif c == '>': pline += '&gt;'
      else: pline += c
    else: pline += c
    
  lines = pline.splitlines()
  if len(lines) <= 1:
    return pline
  linenew = ''
  for i, line in enumerate(lines):
    if i > 0: linenew += '<br/>'
    linenew += line
  return linenew
# - - - - - - - - - - - - - - - - -
def formatRange(range):
  # format two numbers as "n1 - n2"
  values = range.split()
  if len(values) == 2:
    return values[0] + ' - ' + values[1]
  return range

# - - - - - - - - - - - - - - - - -

def median(lst):
    quotient, remainder = divmod(len(lst), 2)
    if remainder:
        return sorted(lst)[quotient]
    return sum(sorted(lst)[quotient - 1:quotient + 1]) / 2.

# - - - - - - - - - - - - - - - - -
class SDcorrectionData:
    ''' Storage of SDcorrection data '''

    #  Warning levels
    MAXSDFAC1 = 2.0
    MAXSDFAC2 = 3.0  # serious
    MAXSDADD1 = 0.06
    MAXSDADD2 = 0.10
    MAXSDFACDIFF = 0.25  # from median
    MAXSDADDDIFF = 0.015

    # Bit flags
    SDFAC_FLAG_MINOR = 1
    SDFAC_FLAG_MAJOR = 2
    SDADD_FLAG_MINOR = 4
    SDADD_FLAG_MAJOR = 8
    SDADD_FLAG_NEGATIVE = 16
    SDFAC_FLAG_DIFF = 32
    SDADD_FLAG_DIFF = 64
    
    # - - - - - - - - - - - - - - - - -
    def __init__(self, label, tag, SDfac, SDB, SDadd, ISa):
      """ label is Run, tag is Full or Partial """
      self.label = label
      self.tag = tag
      self.SDfac = SDfac
      self.SDB = SDB
      self.SDadd = SDadd
      self.ISa = ISa
      self.status = 0
      
      self.validity()  # set status
      #print "Data", self.tag, self.SDfac, self.SDB, self.SDadd, self.ISa

    # - - - - - - - - - - - - - - - - -
    def validity(self):
        ''' Warning levels suggested by Andrew Leslie
        1. Moderate warning: SDfac > 2.0 and/or SDadd>0.06
        2. Serious warning: SDfac> 3.0 and/or SDadd>0.10 or SDadd < 0.0
        '''
        if self.SDfac > SDcorrectionData.MAXSDFAC2:
          self.status |= SDcorrectionData.SDFAC_FLAG_MAJOR
        elif self.SDfac > SDcorrectionData.MAXSDFAC1:
          self.status |= SDcorrectionData.SDFAC_FLAG_MINOR

        if self.SDadd > SDcorrectionData.MAXSDADD2:
          self.status |= SDcorrectionData.SDADD_FLAG_MAJOR
        elif self.SDadd > SDcorrectionData.MAXSDADD1:
          self.status |= SDcorrectionData.SDADD_FLAG_MINOR

        if self.SDadd < 0.0:
          self.status |= SDcorrectionData.SDADD_FLAG_NEGATIVE

    # - - - - - - - - - - - - - - - - -
    def valid(self):
      # return 0 OK, +1 minor error, +2 major error
      if self.status == 0:
        return 0
      if self.status & (SDcorrectionData.SDFAC_FLAG_MAJOR |
                        SDcorrectionData.SDADD_FLAG_MAJOR):
        return +2
      if self.status & SDcorrectionData.SDADD_FLAG_NEGATIVE:
        return +2
      if self.status & (SDcorrectionData.SDFAC_FLAG_MINOR |
                        SDcorrectionData.SDADD_FLAG_MINOR):
        return +1
      if self.status & (SDcorrectionData.SDFAC_FLAG_DIFF |
                        SDcorrectionData.SDADD_FLAG_DIFF):
        return +1
      return 0  # shouldn't get here
      
    # - - - - - - - - - - - - - - - - -
    def validityCompare(self, mediansdfac, mediansdadd):
      #  Andrew Leslie's suggestion
      #  SDFAC differs from median by > 0.25, SDADD differs by > 0.015
      if abs(self.SDfac - mediansdfac) > SDcorrectionData.MAXSDFACDIFF:
        self.status |= SDcorrectionData.SDFAC_FLAG_DIFF
      if abs(self.SDadd - mediansdadd) > SDcorrectionData.MAXSDADDDIFF:
        self.status |= SDcorrectionData.SDADD_FLAG_DIFF

    # - - - - - - - - - - - - - - - - -
    def statusFlags(self):
      # status characters for SDfac, SDadd
      sdfacflag = '  '
      sdaddflag = '  '
      if self.status & SDcorrectionData.SDFAC_FLAG_MINOR:
        sdfacflag = '*'
      if self.status & SDcorrectionData.SDFAC_FLAG_MAJOR:
        sdfacflag = 'x'
      if self.status & SDcorrectionData.SDFAC_FLAG_DIFF:
        sdfacflag = sdfacflag+'#'
      else:
        sdfacflag = sdfacflag+' '

      if self.status & SDcorrectionData.SDADD_FLAG_MINOR:
        sdaddflag = '*'
      if self.status & SDcorrectionData.SDADD_FLAG_MAJOR:
        sdaddflag = 'x'
      if self.status & SDcorrectionData.SDADD_FLAG_NEGATIVE:
        sdaddflag = 'x'
      if self.status & SDcorrectionData.SDADD_FLAG_DIFF:
        sdaddflag = sdaddflag+'#'
      else:
        sdaddflag = sdaddflag+' '

      return sdfacflag, sdaddflag
      
    # - - - - - - - - - - - - - - - - -
    def as_string(self):
      sdfacflag, sdaddflag = self.statusFlags()
      text = self.label + " " + self.tag
      s = '{}: SDfac: {}{}, SdB: {}, SDadd: {}{}, ISa: {}'.format(
        text, self.SDfac, sdfacflag, self.SDB,
        self.SDadd, sdaddflag, self.ISa)
      return s

    # - - - - - - - - - - - - - - - - -
    def as_list(self):
      ''' return label, list for table'''
      sdfacflag, sdaddflag = self.statusFlags()
      text = self.label + " " + self.tag
      sdclist = [
        self.SDfac, sdfacflag, self.SDB,
        self.SDadd, sdaddflag, self.ISa]
      return text, sdclist

    # - - - - - - - - - - - - - - - - -
    def labelISa(self):
      """ return label, ISa """
      return self.label, self.ISa
    # - - - - - - - - - - - - - - - - -
    @staticmethod
    def header(multiple=False):
      s = "SDcorrection parameters are flagged with moderate (*) or "+\
          "severe (x) warning, if SDfac > "+str(SDcorrectionData.MAXSDFAC1)+\
          " or "+str(SDcorrectionData.MAXSDFAC2)+",\n"
      s += " or if SDadd is negative or  > "+str(SDcorrectionData.MAXSDADD1)+\
          " or "+str(SDcorrectionData.MAXSDADD2)+" respectively\n"
      if multiple:
        s += "Parameters are flagged as (#) if they differ from median values"+\
             " of SdFac by > "+str(SDcorrectionData.MAXSDFACDIFF)+\
             " or SdAdd by > "+str(SDcorrectionData.MAXSDADDDIFF)+"\n"
      s += "ISa is the asymptotic maximum I/sig(I) = 1/(SdFac*SdAdd)"
      return s

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
def selectGraphs(xmlnode, baseElement="CCP4Table",
                 graphID=None,
                 graphTitle=None,
                 plotTitleList=None):
  import xml.etree.ElementTree as etree
  '''
  Select from xmlnode (usually = <CCP4Table> element) graphs
  with id= graphID
  Then make subset of "plot"s with (partially) matching plot title strings

   - baseElement     main tag, usually (always?) CCP4Table
   - graphID         partial id for selecting graphs (from attributes)
   - graphTitle      partial title for selecting graphs (from attributes)
   - plotTitleList   list of plots to accept
   
  graphTitle and plotTitleList may be None, to select all graphs

  Return list of graphs with selected subset of plots,
  with other attributes unchanged
  graphID may match multiple graphs, so return list
  Return None if nothing found
  '''
  if graphID is None: return None

  newCCP4TableList = []

  graphid = baseElement+"[@id='"+graphID+"']"

  # list of graphs matching graphID
  graphlist = xmlnode.findall(graphid)
  if len(graphlist) == 0: return None  # no graphs found

  for graph in graphlist:
    plotlist = graph.findall('plot')
    attrib = graph.attrib

    accept = True
    if graphTitle is not None:
      accept = False
      if graphTitle in attrib['title']:
        accept = True
    if accept:
      newplot = []
      for plot in plotlist:
        plottitle = plot.findall('title')[0].text
        include = True
        if plotTitleList is not None:
          include = False
          for title in plotTitleList:
            if title in plottitle:
              # this is the graph we want
              include = True
        if include:
          newplot.append(plot)

      if len(newplot) != 0:
        newCCP4Table = etree.Element("CCP4Table", attrib=attrib)
        for plot in newplot:
          newCCP4Table.append(plot)

        # copy headers and data
        headers = graph.findall('headers')
        newCCP4Table.append(headers[0])
        newCCP4Table.append((graph.findall('data'))[0])
        newCCP4TableList.append(newCCP4Table)  # add to list
  # end graph loop

  if len(newCCP4TableList) == 0: return None
  return newCCP4TableList

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
class CellCheck:
  ''' Check cell compatibility '''

  # - - - - - - - - - - - - - - - - -
  def __init__(self, mtzContent1 ,mtzContent2, tolerance=None):
    '''  Compare cells
    mtzContent1   observed data, mtzContent2   freeR data (usually)
    These are from CObsDataFile.fileContent
    '''
    # print('CellCheck')
    self.obscontents = mtzContent1
    self.obscontents2 = mtzContent2
    self.tolerance = tolerance

  # - - - - - - - - - - - - - - - - -
  def checks(self):
    # Mostly for comparing observed data with freer data,
    #  but also for observed v. observed
    cellsAreTheSame = self.obscontents.clipperSameCell(self.obscontents2,self.tolerance)
    # print('cellsAreTheSame', cellsAreTheSame)
    # Make XML block to report cells and their compatibility
    reportXML = self.cellCompatibilityXML(self.obscontents, self.obscontents2,
                                               cellsAreTheSame)
    return cellsAreTheSame, reportXML
  # - - - - - - - - - - - - - - - - -
  def isValid(self):
    cellsAreTheSame = self.obscontents.clipperSameCell(self.obscontents2,self.tolerance)
    return cellsAreTheSame['validity']

  # - - - - - - - - - - - - - - - - -
  def addElement(self, containerXML, elementname, elementtext):
    from lxml import etree
    e2 = etree.Element(elementname)
    e2.text = elementtext
    containerXML.append(e2)

  # - - - - - - - - - - - - - - - - -
  def cellCompatibilityXML (self,mtzContent1 ,mtzContent2, cellsAreTheSame):
    from lxml import etree
    " Make XML report of cells and their compatibility"
    # Mostly for comparing observed data with freer data,
    #  but also for observed v. observed
    sgname1 = 'Unk'
    sgname2 = 'Unk'
    if mtzContent1.spaceGroup.isSet(): sgname1 = mtzContent1.spaceGroup.__str__()
    if mtzContent2.spaceGroup.isSet(): sgname2 = mtzContent2.spaceGroup.__str__()
    """
    cell1, cell2    cells
    cellsAreTheSame        result dictionary from SameCell
    'validity'           True if cells are simlar within resolution of tolerance
    'maximumResolution1' maximum allowed resolution in cell1
    'maximumResolution2' maximum allowed resolution in cell2
    'difference'  average cell difference in A        
    'tolerance'   in A
    """

    cellReportXML= etree.Element('ObsFreeCellComparison')
    
    cellformat = CellFormat()
    self.addElement(cellReportXML, 'cell1',cellformat.shortformatCell(mtzContent1.cell))
    self.addElement(cellReportXML, 'cell2', cellformat.shortformatCell(mtzContent2.cell))
    self.addElement(cellReportXML,'tolerance', ("%7.2f" % cellsAreTheSame['tolerance']).strip())
    # Validity is true if cells are within tolerance
    self.addElement(cellReportXML, 'validity', str(cellsAreTheSame['validity']))
    # Cell difference
    self.addElement(cellReportXML, 'CellDifference', ("%7.2f" % cellsAreTheSame['difference']).strip())
    # Maximum resolution that extension from FreeR set would be valid
    #  calculated from FreeR cell to data cell
    self.addElement(cellReportXML, 'MaxAcceptableResolution',
                    ("%7.2f" % cellsAreTheSame['maximumResolution2']).strip())
    self.addElement(cellReportXML, 'sgname1', sgname1)
    self.addElement(cellReportXML, 'sgname2', sgname2)
    

    #print "cellReportXML",etree.tostring(cellReportXML,pretty_print=True)
    return cellReportXML
# end class CellCheck

class CellFormat:
  # - - - - - - - - - - - - - - - - -
  def formatCellLength(self, p):
    return "%7.1f" % float(p)

  # - - - - - - - - - - - - - - - - -
  def formatCellAngle(self, p):
    if float(p) < 10.0:
      return "%7.1f" % (float(p) * 57.29577951308233)
    else:
      return "%7.1f" % float(p)

  # - - - - - - - - - - - - - - - - -
  def shortformatCell(self, cell):
    # print 'shortformatCell',cell
    s = ""
    s += self.formatCellLength(cell.a).strip()+', '
    s += self.formatCellLength(cell.b).strip()+', '
    s += self.formatCellLength(cell.c).strip()+', '
    s += self.formatCellAngle(cell.alpha).strip()+', '
    s += self.formatCellAngle(cell.beta).strip()+', '
    s += self.formatCellAngle(cell.gamma).strip()
    return s

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#if __name__ == "__main__":
#    sdcd = SDcorrectionData("label", 1.1, 5.3, 0.01, 20.0)
#    print sdcd.tag, sdcd.SDfac
#    print sdcd.format()
