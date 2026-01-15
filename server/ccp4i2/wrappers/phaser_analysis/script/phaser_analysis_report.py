from ccp4i2.report import Report
from ccp4i2.wrappers.phaser_analysis.script.phaser_analysis_utils import Tabledata


class phaser_analysis_report(Report):
  TASKNAME='phaser_analysis'

  # - - - - - - - - - - - - - - - - -
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
    Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=None,**kw)

    # For now (testing), isolate PHASER_ANALYSES:1 block 
    select1 = False   # True to use PHASER_ANALYSES1 if present
    if self.xmlnode.tag == 'IMPORT_MERGED':
      self.xmlnode = self.xmlnode.findall('AIMLESS_PIPE/PHASER_ANALYSES/PHASER_ANALYSIS')[0]
    elif self.xmlnode.tag == 'AIMLESS_PIPE':
      if select1 and len(self.xmlnode.findall('PHASER_ANALYSES1')) > 0:
        self.xmlnode = self.xmlnode.findall('PHASER_ANALYSES1/PHASER_ANALYSIS')[0]
      else:
        self.xmlnode = self.xmlnode.findall('PHASER_ANALYSES/PHASER_ANALYSIS')[0]

    try:
      self.fileroot = self.jobInfo['fileroot']
    except:
      self.fileroot = None

    self.datasetsProcessed = False

    # Check for fail message
    self.fail = False
    if len(self.xmlnode.findall('PhaserFailMessage'))>0:
      self.fail = True
      self.errorMessage = self.xmlnode.findall('PhaserFailMessage')[0]

    if not self.fail: self.harvest()  # Harvest all data into internal storage
    ##self.dump()  # for testing
    
    
    # 'nooutput' mode would be used by another report class that wanted
    # to use some method(s) from this class for its own report
    if jobStatus is not None and jobStatus.lower() == 'nooutput':
      return
  
    #fail = self.Errors(self)
    #if fail == False:

    self.phaserAnalysisReport(self, False)

  # - - - - - - - - - - - - - - - - -
  def phaserAnalysisReport(self, parent=None, drawresol=True):
    ''' report for single Phaser_analysis step '''
    # Data for one dataset has been harvested into self.data and self.graphs

    if parent is None: parent = self

    heading = 'Analysis of dataset ['+self.pxdname+'] using Phaser'
    parent.addText(text=heading, style='text-align:center: border:1px solid black;')

    mainDiv = parent.addDiv(style="width:100%;margin:0px; padding:0px;")

    # Use grid layout for twin analysis and resolution side by side
    leftDiv, rightDiv = mainDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)

    if drawresol:
      self.drawresolution(rightDiv)

    self.drawtwin(leftDiv)

    # Use grid layout for NCS and anisotropy analysis side by side
    ncsDiv, anisoDiv = mainDiv.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
    self.drawtNCS(ncsDiv)
    self.drawanisotropy(anisoDiv)

  # - - - - - - - - - - - - - - - - -
  def getInformationValues(self, resrange):
    # Get information values for resolution range for table 1
    mn = self.meanInformation(resrange,
                              'informationcontent', 'Selected')
    return mn

  # - - - - - - - - - - - - - - - - -
  def getdatavalue(self, tag):
    # return data[tag] if exists, else None
    if self.fail: return
    if tag in self.data:
      return self.data[tag]
    return None
  # - - - - - - - - - - - - - - - - -
  def drawresolution(self,parent=None):
    # Report from phaser_analysis on merged data

    angstrom = "&#197;"
    resall = self.data['resolutionAll']
    ressel = self.data['selectedresolution']
    resest = self.data['resolutionestimate']
    s  = "Input resolution:   "+resall+angstrom
    if float(ressel)-float(resall) > 0.05:
      s += "  (no significant data beyond "+ressel+angstrom+")"
    s += "\nResolution estimate: "+self.data['resolutionestimate']+angstrom
    s += "\nthreshold on information content "+\
         self.data['threshold']+" bits/reflection"

    s = html_linebreak(s, True)
    parent.append(s) 
    graphDiv = parent.addDiv(\
      style="width:100%;border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

    self.drawgraph('informationcontent', graphDiv)

  # - - - - - - - - - - - - - - - - -
  def drawgraph(self, graphid, parent=None):

    # Loop over all Graph tables in the program output and add to the GraphGroup
    # The plotting instructions are provided as xml text
    # Valid graphid = '2ndmoment', 'informationcontent'

    graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;")
    if graphid not in self.graphs: return
    thisgraph = self.graphs[graphid]
    graph = graphgroup.addFlotGraph(xmlnode=thisgraph,
                                    title=thisgraph.get("title") )
    graph = graph.addPimpleData(xmlnode=thisgraph)
      
  # - - - - - - - - - - - - - - - - -
  def drawtNCS(self, parent):
    if 'tNCS' not in self.data: return
    if self.data['tNCS'] == 'False':
      parent.append('No translational NCS detected')
      return

    s  = "Translational NCS detected:\n"
    s += "  Native Patterson peak height: "+self.data['PattPeakHeight']+\
         "% of origin peak\n"
    s += "  Patterson vector: "+self.data['PattVector']+"\n"
    s += "  Range of D-values due to tNCS: "+self.data['DvalueRange']
    s = html_linebreak(s, True)
    parent.append(s) 

  # - - - - - - - - - - - - - - - - -
  def drawtwin(self, parent):
    pvalues = self.data['Pvalues']
    pv = pvalues.split()

    s = ''
    jtwin = self.istwin()
    # return = 0 untwinned, +2 probably twinned, +1 possible twin
    
    if jtwin == +2:
      s = 'Data may be twinned:\n P(untwinned) = '+pv[0]+\
          '\n P(twin fraction < 5%) = '+pv[1]
    elif jtwin == +1:
      s = 'Some possibility that the data are twinned:\n P(untwinned) = '+pv[0]+\
          '\n P(twin fraction < 5%) = '+pv[1]
    else:
      s  = 'Data appear to be untwinned'
      s += '\nP(untwinned) = '+pv[0]+\
          '\n P(twin fraction < 5%) = '+pv[1]
          
    s = html_linebreak(s, True)
    parent.append(s)
    self.momentTable(parent)

    self.drawgraph('2ndmoment', parent)
    self.drawgraph('intensitydistribution', parent)

  # - - - - - - - - - - - - - - - - -
  def istwin(self):
    # return = 0 untwinned, +2 probably twinned, +1 possible twin
    pvalues = self.data['Pvalues']
    pv = pvalues.split()

    flag = 0
    if float(pv[1]) < 0.01:
      flag = +2
    elif float(pv[1]) < 0.4:
      flag = +1
    else:
      flag = 0
    return flag

  # - - - - - - - - - - - - - - - - -
  def momentTable(self, parent):
    table = parent.addTable(class_="center")
    rows = ['Theoretical untwinned', 'Theoretical twinned','Observed','SDobserved']
    table.addData(title='2nd moments',data=rows)
    col = [self.data['CentricTheoretical'],
           self.data['CentricTheoreticalTwin'],
           self.data['CentricObserved'], '-']
    table.addData(data=col, title='Centric')
    col = [self.data['AcentricTheoretical'],
           self.data['AcentricTheoreticalTwin'],
           self.data['AcentricObserved'],
           self.data['AcentricSDobserved']]
    table.addData(data=col, title='Acentric')

  # - - - - - - - - - - - - - - - - -
  def drawanisotropy(self, parent):
    
    s  = "Anisotropic deltaB (i.e. range of principal components): "+\
        self.data['AnisotropicDeltaB']+"\n"
    s += "Eigenvectors are direction cosines in orthogonal coordainates"
    s = html_linebreak(s, True)
    parent.append(s)
    self.anisotropyTable(parent)

  # - - - - - - - - - - - - - - - - -
  def anisotropyTable(self, parent):
    table = parent.addTable(class_="center")
    eigenvalues = self.data['AnisotropicEigenvalues']
    evls = eigenvalues.split()
    colvl = []
    eigenvectors = self.data['AnisotropicEigenvectors']
    evcs = eigenvectors.split()
    colvc  = []
    for i in range(3):
      colvl.append(evls[i])
      i1 = i*3
      colvc.append('{:8s} {:8s} {:8s}'.format(*evcs[i1:i1+3]))

    table.addData(title='Eigenvalues', data=colvl)
    table.addData(title='Eigenvectors', data=colvc)
    
  # - - - - - - - - - - - - - - - - -
  def harvest(self):
    # Harvest all XML data into internal storage
    # Creates two dictionaries:
    #   self.graphs    graphs
    #   self.data      everything else

    self.pxdname = self.xmlnode.attrib['name']

    # Graphs
    self.graphs = {}
    graphlist = self.xmlnode.findall("CCP4Table")
    for thisgraph in graphlist:
      graphid = thisgraph.attrib['id']
      self.graphs[graphid] = thisgraph

    self.data = {}

    # Resolution stuff
    analysisnode = self.xmlnode.findall('Analysis')[0]
    self.data['resolutionAll'] = self.getitem(analysisnode, 'Resolution/ResolutionAll')
    # "selected resolution" is the reslution to which Phaser thinks there
    # may be some information (maybe > 0) and is selected for analysis
    # It is distinct from resolutionestimate which comes from a
    # looking for a specified threshold
    self.data['selectedresolution'] = self.getitem(analysisnode, 'Resolution/SelectedResolution')
    self.data['totalbits'] = self.getitem(analysisnode,'InformationContent/TotalBits')
    self.data['nreflections'] = self.getitem(analysisnode,'InformationContent/Nreflections')
    self.data['averagebits'] = self.getitem(analysisnode,'InformationContent/Averagebits')

    self.data['resolutionestimate'] = self.getitem(self.xmlnode,
                       'ResolutionEstimate/ResolutionLimitEstimate')
    self.data['threshold'] = self.getitem(self.xmlnode,
                                      'ResolutionEstimate/Threshold')
    self.tNCS(analysisnode)
    self.twinning(analysisnode)
    self.anisotropy(analysisnode)
    

  # - - - - - - - - - - - - - - - - -
  def tNCS(self, analysisnode):
    tncsnode = analysisnode.findall('tNCS')[0]
    if tncsnode is None: return
    self.data['tNCS'] = tncsnode.attrib['tNCS']
    if self.data['tNCS'] == 'False':  #  NB str type!
      return

    self.data['PattPeakHeight'] = self.getitem(tncsnode,
                                               'NonOriginPatterson/PeakHeight')
    self.data['PattVector'] = self.getitem(tncsnode,
                                               'NonOriginPatterson/Vector')
    self.data['PattAngle'] = self.getitem(tncsnode,
                                               'NonOriginPatterson/Angle')
    self.data['DvalueRange'] = self.getitem(tncsnode,'DvalueRange')

  # - - - - - - - - - - - - - - - - -
  def twinning(self, analysisnode):
    twinningnode = analysisnode.findall('Twinning')[0]

    # Two Pvalues, probability of !1) twinned (2) twin < 5%
    self.data['Pvalues'] = self.getitem(twinningnode, 'Pvalues')

    for cen in ['Centric', 'Acentric']:
      moments = cen+'Moments'
      self.data[cen+'Theoretical'] = self.getitem(twinningnode,
                                                  moments+'/Theoretical')
      self.data[cen+'TheoreticalTwin'] = self.getitem(twinningnode,
                                                      moments+'/TheoreticalTwin')
      self.data[cen+'Observed'] = self.getitem(twinningnode,
                                               moments+'/Observed')
      self.data[cen+'SDobserved'] = ''
      if cen == 'Acentric':
        self.data[cen+'SDobserved'] = self.getitem(twinningnode,
                                                   moments+'/SDobserved')
  # - - - - - - - - - - - - - - - - -
  def anisotropy(self, analysisnode):
    anisotropynode = analysisnode.findall('Anisotropy')[0]

    self.data['AnisotropicDeltaB'] = self.getitem(anisotropynode,
                                                  'deltaB')
    self.data['AnisotropicEigenvalues'] = self.getitem(anisotropynode,
                                                  'Eigenvalues')
    self.data['AnisotropicEigenvectors'] = self.getitem(anisotropynode,
                                                  'Eigenvectors')

  # - - - - - - - - - - - - - - - - -
  def getitem(self, node, tag):
    # return item from node
    return node.findall(tag)[0].text

  # - - - - - - - - - - - - - - - - -
  def dump(self):
    # for debugging
    print("\n****phaser_analysis_report.dump")
    for key in self.graphs:
      print("Graph: ",key)

    for key, content in self.data.items():
      print("Data: ", key,"; ",content)

  # - - - - - - - - - - - - - - - - -
  def hasfailed(self):
    return self.fail

  # - - - - - - - - - - - - - - - - -
  def meanInformation(self, resrange, graphname, column):
    #print("meanInformation", resrange, graphname)

    if graphname not in self.graphs:
      print("Unrecognised graph name", graphname)
      return None
    
    # convert range from d to 1/d^2
    dstarrange = [1.0/(resrange[0]*resrange[0]),
                  1.0/(resrange[1]*resrange[1])]

    table = Tabledata(self.graphs[graphname])
    columns = table.getcollabels()
    if column not in columns:
      print("Unrecognised column", column)
      return None

    meanvalue = table.meanValue(dstarrange, column, xcol='1/d^2')

    return meanvalue
    

######################################################################
    
def html_linebreak(line, mutate=True):
  # delineate lines by "<br/>"
  # Optionally mutate characters '&', '<', '>' into safe equivalents
  #    but accept '&#' combination

  pline = ''

  if mutate:
      ll = len(line)
      for i, c in enumerate(line):
          if c == '&':
              if line[min(i+1, ll-1)] == '#':
                  pline += c
              else:
                  pline += '&amp;'
          elif c == '<': pline += '&lt;'
          elif c == '>': pline += '&gt;'
          else: pline += c
  else:
      pline = line

  # check for "\n"
  lines = pline.splitlines()
  if len(lines) <= 1:
    return pline
  linenew = ''
  for line in lines:
    linenew += line + '<br/>'

  return linenew
