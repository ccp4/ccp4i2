import sys
import xml.etree.ElementTree as etree

from ....report.CCP4ReportParser import Report


class scaleit_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'scaleit'
  RUNNING = True
  
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,style="overflow:auto;",**kw)

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    self.finalReport()


  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def finalReport(self, parent=None):
      if parent is None: parent = self

      self.filetypes = None
      if len(self.xmlnode.findall('InputFileTypes'))>0:
          self.filetypes = self.xmlnode.findall('InputFileTypes')[0]

      if len(self.xmlnode.findall('SCALEITLOG'))>0:
          self.scaleitxml = self.xmlnode.findall('SCALEITLOG')[0]
      else:
          self.scaleitxml = self.xmlnode

      if len(self.xmlnode.findall('FATAL_ERROR'))>0:
        text = self.xmlnode.findall('FATAL_ERROR')[0].text
        errorDiv = parent.addDiv(
          style="width:90%;border: 2px solid red; clear:both; margin:3px; padding:6px;")
        errorDiv.addText(text='FATAL ERROR',
                         style='font-weight:bold; font-size:150%; color:red;')
        errorDiv.append(str(message))
        return

      topDiv = parent.addDiv(style="font-size:100%;border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;")

      topDiv.addText(text="Comparison of datasets using SCALEIT",
                     style='font-weight:bold; font-size:150%; color:blue;')
      topDiv.append(' <br/>')

      self.nderivatives = self.scaleitxml.findall('Nderivatives')[0].text
      #  Derivative data blocks
      self.derivatives  = self.scaleitxml.findall('Derivative')

      # Get column names for each derivative, index graphs from their titlea
      self.dcolnames = []
      for deriv in self.derivatives:
          self.dcolnames.append(deriv.findall('ColName')[0].text)

      # Graph list
      gpath = 'SCALEITGRAPHS/CCP4ApplicationOutput/CCP4Table'
      self.graphlist = self.scaleitxml.findall(gpath)

      # Dataset names
      self.nativename = self.scaleitxml.findall('NativeDname')[0].text
      self.dnames = self.dataList('Name', self.derivatives)

      message = 'Compare each "derivative" dataset against the first "native" dataset'
      topDiv.append('<br/>')
      topDiv.addText(text=message,
                     style='font-weight:bold; font-size:110%; color:blue;')

      message = '  Dataset ' + self.nativename + ' is treated as "native"'
      if len(self.dnames) == 1:
          message += '\n  Dataset ' + self.dnames[0] + ' is treated as "derivative"'
      else:
          message += '\n  Datasets ' + ' & '.join(self.dnames) + ' are treated as "derivatives"'
      message = html_linebreak(message)
      topDiv.append(message)
      
      # Input file content, maybe converted
      if self.filetypes is not None:
        self.inputFileTypes(topDiv)
        topDiv.append('<br/>')

      if len(self.scaleitxml.findall('ResolutionMax'))>0:
        maxres = self.scaleitxml.findall('ResolutionMax')[0].text
        topDiv.append('<br/>')
        message = 'Note: Resolution was cut to ' + maxres + ' \xc5'
        topDiv.addText(text=message,style='color: blue')

      scaleDiv = parent.addDiv(style="font-size:100%;border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;")
      text = "Results from scaling"
      scaleDiv.append('<br/>')
      scaleDiv.addText(text=text,
                     style='font-weight:bold; font-size:120%; color:blue;')
      self.scaleResults(scaleDiv)

      resoDiv = parent.addDiv(style="font-size:100%;border-width: 1px; border-color: black;clear:both; margin:0px; padding:0px;")
      text = "Analysis by resolution:"
      resoDiv.append('<br/>')

      resoDiv.addText(text=text,
                     style='font-weight:bold; font-size:120%; color:blue;')
      self.resolutionAnalysis(resoDiv)

      if len(self.scaleitxml.findall('NormalProbability'))>0:
          self.drawNormalProbabilityGraphs(parent)
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def inputFileTypes(self, parent):
    # Describe the content of each input file
    #   self.filetypes is XML block
    nativecontent = str(self.filetypes.findall('Native')[0].text)
    parent.addText(text='Input "native" '+self.nativename+\
                   ' data were '+
                   self.fileContentComment(nativecontent))
    derivs = self.filetypes.findall('Derivative')
    for i in range(len(derivs)):
      parent.append(' <br/>')
      content = derivs[i].text
      message = 'Input "derivative" '+self.dnames[i]+\
                ' data were '+\
                self.fileContentComment(content)
      parent.addText(text=message)
    
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def fileContentComment(self, filecontent):
    s = ''
    if filecontent == 'Fmean':
      s = 'Fmean, no need to convert for SCALEIT'
    else:
      s = filecontent+', converted to Fmean for SCALEIT'
    return s
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def scaleResults(self, scaleDiv):
      for i in range(len(self.dnames)):
          scale = self.derivatives[i].findall('ScaleFactor')[0].text
          tag = 'AnisoB'
          aniso = True
          if len(self.derivatives[i].findall('IsoB'))>0:
              tag = 'IsoB'
              aniso = False
          bfactor = self.derivatives[i].findall(tag)[0].text
          scaleDiv.append('<br/>')
          message = 'Derivative ' + self.dnames[i] +\
                    ', scale factor relative to '+\
                    self.nativename + ' = ' + scale
          scaleDiv.addText(text=message)
          scaleDiv.append('<br/>')
          if aniso:
              message = 'Anisotropic relative B-factors = '
          else:
              message = 'Isotropic relative B-factor = '
          scaleDiv.addText(text=message+bfactor)
          scaleDiv.append('<br/>')

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def resolutionAnalysis(self, resoDiv):

      # Get longest headers
      tableheaders = []
      hasanom = False
      for derivblock in self.derivatives:
          resoElement = derivblock.findall('ResolutionAnalysis')[0]
          headers = self.validStringList(resoElement.findall('Headers')[0].text)
          # remove first 2
          headers = headers[2:]
          if len(headers) > len(tableheaders):
              tableheaders = headers

          # Is there anomalous?
          if (derivblock.findall('Anomalous')[0].text == 'true'):
            hasanom = True  # Any with anomalous

      if not hasanom:
        # Prune out anomalous headers
        tableheaders = tableheaders[:-5]

      resoDiv.append('<br/>')
      resoDiv.addText(text="Totals over all resolution ranges")
      resoDiv.append('<br/>')

      totaltable = resoDiv.addTable(transpose=True,
          style="line-height:80%%; font-size:75%;",
                                    downloadable=False)
      totaltable.addData(title="Dataset", data=tableheaders)

      for i in range(len(self.dnames)):  # Loop derivatives
          dname = self.derivatives[i].findall('Name')[0].text
          data = resoElement.findall('Totals')[0].text.split()
          # remove first 2
          data = data[2:]
          totaltable.addData(title=dname, data=data)

          s = 'Resolution analysis for derivative ' + self.dnames[i] +\
              ' relative to ' + self.nativename
          resoDiv.append('<br/>')
          resoDiv.addText(text=s,
                     style='font-weight:bold; font-size:120%; color:blue;')
          resoDiv.append('<br/>')
          resoElement = self.derivatives[i].findall('ResolutionAnalysis')[0]

          # Graphs
          graphgroup = resoDiv.addFlotGraphGroup\
                       (style="width:500px;  height:300px;")
          thisgraph = self.selectResolutionGraph(i)  # for i'th derivative
          plottitles = [\
            'Rfac/Wted_R v resolution',
            '<diso> and Max(diso) v resolution',
            'Kraut_sc and RMS(FP/FHP)  v resolution']
          if hasanom:
            plottitles.append('<dano> and Max(dano) v resolution')
          thisgraph = self.editPlot(thisgraph, plottitles=plottitles)

          graph = graphgroup.addFlotGraph(xmlnode=thisgraph,
                                          title=thisgraph.get("title") )
          graph = graph.addPimpleData(xmlnode=thisgraph)

          
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def selectResolutionGraph(self, ideriv):
          # Pick out graphs for this dataset
          #print("selectResolutionGraph", self.graphlist)
          dcolname = self.dcolnames[ideriv]
          for thisgraph in self.graphlist:
              title = thisgraph.get('title')
              if "Analysis v resolution" in title and dcolname in title:
                  return thisgraph
          return None
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def drawNormalProbabilityGraphs(self, npDiv):

      npDiv.append('<br/>')
      npDiv.addText(\
          text="Normal probability analysis of isomorphous differences",
          style='font-weight:bold; font-size:120%; color:blue;')

      text = "\nGradient should be 1 for two equivalent isomorphous data sets\n"+\
      "Gradient should be >> 1 for different data sets\n\n"+\
      "Gradient2 etc are for filtered data excluding extremes of distribution, cutoff 0.9" 
      npDiv.append(html_linebreak(text))

      normalprobability  = self.scaleitxml.findall('NormalProbability')[0]
      totaltable = npDiv.addTable(transpose=True,
          style="line-height:80%%; font-size:75%;",
                                   downloadable=False)
      tableheaders = normalprobability.findall('Headers')[0].text.split()
      totaltable.addData(title="", data=tableheaders)
      data = normalprobability.findall('Centric')[0].text.split()
      totaltable.addData(title='Centric', data=data)
      data = normalprobability.findall('Acentric')[0].text.split()
      totaltable.addData(title='Acentric', data=data)

      # Graphs
      divstyle = "width:47%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;"
      #divstyle = "width:47%;float:left;text-align:center;margin:0px; padding:0px; line-height:100%; font-size:100%;border:1px solid black;"
      leftDiv = npDiv.addDiv(style=divstyle)
      rightDiv = npDiv.addDiv(style=divstyle)
      
      graphstyle = "width:300px;  height:300px;"
      thisgraph = self.getGraph("Centric Normal probability")
      thisgraph = self.editPlot(thisgraph, addresol=True)
      graphgroup = leftDiv.addFlotGraphGroup(style=graphstyle)
      graph = graphgroup.addFlotGraph(xmlnode=thisgraph,
                                      title=thisgraph.get("title") )
      graph = graph.addPimpleData(xmlnode=thisgraph)
      thisgraph = self.getGraph("Acentric Normal probability")
      thisgraph = self.editPlot(thisgraph, addresol=True)
      graphgroup = rightDiv.addFlotGraphGroup(style=graphstyle)
      graph = graphgroup.addFlotGraph(xmlnode=thisgraph,
                                      title=thisgraph.get("title") )
      graph = graph.addPimpleData(xmlnode=thisgraph)

      npDiv.append('<br/>')

  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def getGraph(self, tag):
      for gr in self.graphlist:
          title = gr.get('title')
          if tag in title:
              return gr
      return None
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def editPlot(self, thisgraph, plottitles=None, addresol=True):
    # Edit the plot elements in thisgraph
    # If plottitles is given, reorder the plot elements into that order
    # If addresol True, add <xscale>oneoversqrt</xscale>
    
    newplotelements = []
    plotelements = thisgraph.findall('plot') # plot elements
    if plottitles is not None:
      for pt in plottitles:
        for plotelement in plotelements:
          if pt in plotelement.findall('title')[0].text:
            newplotelements.append(plotelement)
            break
    else:
      newplotelements = plotelements

    if addresol:
      temp = []
      for plotelement in newplotelements:
        addElement(plotelement, 'xscale', 'oneoversqrt')
        temp.append(plotelement)
      newplotelements = temp

      newblock = etree.Element(thisgraph.tag)
      newblock.attrib['title'] = thisgraph.get('title')
      for elem in thisgraph:
          if elem.tag != 'plot':
              newblock.append(elem)
      for elem in newplotelements:
          newblock.append(elem)

      return newblock
      
  # . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .  . . . .
  def dataList(self, tag, xmlelements):
      #  xmlelements  list of xml elements
      data = []
      for element in xmlelements:
          text = element.findall(tag)[0].text
          data.append(text)
      return data

  # - - - - - - - - -  - - - - - - - - -  - - - - - - - - -
  def validStringList(self, labels):
      # Split labels into list and make XML-valid
      slist = labels.split()
      vlist = []
      for label in slist:
          vlist.append(self.validTag(label))
      return vlist
  # - - - - - - - - -  - - - - - - - - -  - - - - - - - - -
  def validTag(self, tag):
      # Strip XML invalid characters from tag
      invalid = ['<', '>', '(', ')', '|']
      t = ''
      for c in tag:
          if c not in invalid:
              t +=c
      # Check for '|', if so prepend 'Mod'
      if '|' in tag:
          t = 'Mod'+t
      if '<' in tag:
          t = 'Mean'+t
      return t
############################################################################
# - - - - - - - - -  - - - - - - - - -  - - - - - - - - - 
def addElement(containerXML, elementname, elementtext):
    #print 'addElement', elementname, type(elementtext), elementtext 
    e2 = etree.Element(elementname)
    e2.text = elementtext
    containerXML.append(e2)

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



############################################################################
if __name__ == "__main__":

  report = scaleit_report(xmlFile = sys.argv[1],jobStatus="Finished" )
  tree= report.as_etree()
  #print etree.tostring(tree,pretty_print=True)
  report.as_html_file(fileName='./test-scaleit.html')
  if len(report.errorReport())>0: print('ERRORS:',r.errorReport())
