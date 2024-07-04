from __future__ import print_function
import sys
from report.CCP4ReportParser import *
from wrappers.refmac_i2.script.refmac_report import refmac_report
#from lxml import etree
from xml.etree import ElementTree as ET
from wrappers.sheetbend.script.sheetbend_report import sheetbend_report

class molrep_pipe_report(refmac_report):
  jobStatus = None
  TASKNAME = 'molrep_pipe'
  RUNNING = True
  SEPARATEDATA=True
  COLOUR = (
    ('ice blue',     '#9cb0fe', '156,176,254,255'),
    ('gold',         '#b3b13d', '179,177,61,255'),
    ('coral',        '#ff7f50', '255,127,80,255'),
    ('grey',         '#808080', '128,128,128,255'),
    ('pink',         '#ff92ff', '255,146,255,255'),
    ('sea green',    '#7fbbb5', '127,187,181,255'),
    ('pale brown',   '#a97d5e', '169,125,94,255'),
    ('lilac',        '#ae87b9', '174,135,185,255'),
    ('lemon',        '#ffff80', '255,255,128,255'),
    ('lawn green',   '#459c4f', '69,156,79,255'),
    ('pale crimson', '#d23c3e', '210,60,62,255'),
    ('light blue',   '#419ae1', '65,154,225,255'),
    ('tan',          '#780000', '120,0,0,255'),
    ('light green',  '#9aff9a', '154,255,154,255'),
    ('yellow',       '#ffff00', '255,255,0,255'),
    ('white',        '#ffffff', '255,255,255,255'),
    ('blue',         '#0000ff', '0,0,255,255'),
    ('red',          '#ff0000', '255,0,0,255'),
    ('green',        '#00ff00', '0,255,0,255'),
    ('magenta',      '#ff00ff', '255,0,255,255'),
    ('cyan',         '#00ffe1', '0,255,225,255'),
    ('purple',       '#9400ff', '148,0,255,255'),
    ('dark purple',  '#922057', '146,32,87,255'),
    ('dark cyan',    '#1095a6', '16,149,166,255'),
    ('black',        '#000000', '0,0,0,255'),
  )

  def molrep_plot_xml(self, title, nmon, yincr, label_list):
    e0 = ET.Element('plot')
    e0.text = '\n'
    e1 = ET.SubElement(e0, 'title')
    e1.text = title
    e1.tail = '\n'
    e1 = ET.SubElement(e0, 'plottype')
    e1.text = 'xy'
    e1.tail = '\n'
    for imon in range(nmon):
      e1 = ET.SubElement(e0, 'plotline')
      e1.attrib['xcol'] = str(3* imon + 1)
      e1.attrib['ycol'] = str(3* imon + yincr)
      e1.text = '\n'
      e1.tail = '\n'
      e2 = ET.SubElement(e1, 'colour')
      e2.text = self.COLOUR[imon][1]
      e2.tail = '\n'
      e2 = ET.SubElement(e1, 'label')
      e2.text = label_list[imon]
      e2.tail = '\n'

    return ET.tostring(e0)

  def molrep_report(self, parent, prefix, sgtest=False):
    label_list = list()
    nmon = 0
    while len(self.xmlnode.findall(prefix + '/RFsorted' + str(nmon) + '/RFpeak/RF'))>0:
      nmon += 1
      label_list.append('Copy %d' %nmon)

    if sgtest:
       table = parent.addTable(select=prefix + '/laue_group_alternatives')
       for title, select in [['Space group', 'test/space_group'],
                             ['Score', 'test/score'],
                             ['Contrast', 'test/contrast'],
                             ['Selected', 'test/selected']]:
        table.addData(title=title, select=select)
        node_list = self.xmlnode.findall(prefix + '/laue_group_alternatives/test/space_group')
        if len(node_list) == nmon:
          del label_list[:]
          for node in node_list:
            label_list.append(node.text.strip())

    title = 'Best TF peak vs RF peak No'
    style = 'width:450px;'
    graph = parent.addFlotGraph(title=title, select=prefix, style=style)
    for imon in range(nmon):
      graph.addData(title='RF_peak_No', select='RFsorted%d/RFpeak/RF' %imon)
      graph.addData(title='Score', select='RFsorted%d/RFpeak/Score' %imon)
      graph.addData(title='Tf_sig', select='RFsorted%d/RFpeak/TF_sig' %imon)

    title = 'Best TF peak score vs RF peak No'
    graph.addPlot(plot=self.molrep_plot_xml(title, nmon, 2, label_list))
    title = 'TF/sig(TF) vs RF peak No'
    graph.addPlot(plot=self.molrep_plot_xml(title, nmon, 3, label_list))

    fold = parent.addFold(label='Show details')
    for imon in range(nmon):
      fold.addText(text=label_list[imon])
      table = fold.addTable(select=prefix + "/RFsorted%d/RFpeak" %imon)
      for title, select in [["RF", "RF"],
                            ["TF", "TF"],
                            ["Tf_sig", "TF_sig"],
                            ["TFcntrst", "TFcntrst"],
                            ["PFind", "PFind"],
                            ["PF", "PF"],
                            ["PFmin", "PFmin"],
                            ["wRfac", "wRfac"],
                            ["Score", "Score"],
                            ["Cntrst", "Cntrst"],
                            ["For", "for"]]:
        table.addData(title=title, select=select)

  def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
    self.outputXml = False if self.jobStatus is None else self.jobStatus.lower().count('running')
    print("MOLREP REPORT"); sys.stdout.flush()
    self.xmlnode = xmlnode
    Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
    molrep1done = len(xmlnode.findall('.//MolrepSpaceGroup'))>0
    molrep2done = len(xmlnode.findall('.//MolrepSearch'))>0
    sheetbendNode = None
    if len(self.xmlnode.findall('.//SheetbendResult'))>0:
        sheetbendNode = self.xmlnode.findall('.//SheetbendResult')[0]
    refmacrunning = len(xmlnode.findall('.//RefmacRunning'))>0
    refmacdone = len(xmlnode.findall('.//Refmac'))>0
    sheetbenddone = sheetbendNode is not None or refmacdone or refmacrunning

    if not refmacdone:
      self.addText(text='Job is running; see subjob folders for details')

    if molrep1done or molrep2done:
      results = self.addResults()

    if molrep1done:

      try:
          space_groups = self.xmlnode.findall(".//MolrepSpaceGroup/laue_group_alternatives/test/space_group")
          SGselecteds = self.xmlnode.findall(".//MolrepSpaceGroup/laue_group_alternatives/test/selected")
          sgChanged = False
          idxs = [i for i, x in enumerate(SGselecteds) if x.text == "yes"]
          if len(idxs) == 1 and idxs[0] != 0:
              sgChanged = True
          if sgChanged:
              results.append("<p><b><font color='orange'>Warning: Space group changed from input data's "+space_groups[0].text+" to "+space_groups[idxs[0]].text+" </font></b></p>")
      except:
          pass #This is just information. If for whatever reason the above fails it is probably not important.

      molrep1 = results.addFold(label='Space group selection', initiallyOpen=not molrep2done, brief='Space group')
      self.molrep_report(molrep1, './/MolrepSpaceGroup', True)

    if molrep2done:
      molrep2 = results.addFold(label='Molecular replacement', initiallyOpen=not sheetbenddone, brief='Molrep')
      self.molrep_report(molrep2, './/MolrepSearch')

    if sheetbendNode is not None:
      sheetbendReport = sheetbend_report(xmlnode=sheetbendNode, jobStatus='nooutput')
      opened = not refmacdone and not refmacrunning
      sheetbendFold = self.addFold(label='Shift field refinement',initiallyOpen=opened,brief='Shift field')
      sheetbendReport.defaultReport(parent=sheetbendFold,select=".//SheetbendResult")

    if refmacdone:
      self.addSummary()
      pictureFold = self.addFold(label='Picture', initiallyOpen=False)
      pictureFold.addText(text='View of the best model')
      pic = pictureFold.addPicture(label='', sceneFile="$CCP4I2/wrappers/molrep_mr/script/molrep_mr_1.scene.xml")

    elif refmacrunning:
      refmacFold = self.addFold(label='Refmac', initiallyOpen=True)
      self.addProgressTable2(refmacFold)
      self.addProgressGraph2(refmacFold)

  def addProgressGraph2(self, parent, internalId="SummaryGraph"):
      progressGraph = parent.addFlotGraph(title="Running refmac",select=".//RefmacRunning/Cycle",style="height:250px; width:400px;float:left;border:0px;",outputXml=self.outputXml,internalId=internalId)
      progressGraph.addData(title="Cycle",    select="number")
      progressGraph.addData(title="R_Factor", select="r_factor")
      progressGraph.addData(title="R_Free",  select="r_free")
      plot = progressGraph.addPlotObject()
      plot.append('title','Running refmac R-factors')
      plot.append('plottype','xy')
      plot.append('yrange', rightaxis='false')
      plot.append( 'xlabel', 'Cycle' )
      plot.append( 'xintegral', 'true' )
      plot.append( 'ylabel', 'R-factor' )
      plot.append( 'rylabel', 'Geometry' )
      for coordinate, colour in [(2,'blue'),(3,'green')]:
          plotLine = plot.append('plotline',xcol=1,ycol=coordinate,rightaxis='false',colour=colour)

      rmsBonds = self.xmlnode.findall('.//RefmacRunning/Cycle/rmsBonds')
      if len(rmsBonds)> 0:
          plot.append('yrange', rightaxis='true')
          cycleNodes = self.xmlnode.findall('.//RefmacRunning/Cycle')
          data = []
          for cycleNode in cycleNodes:
              try: data.append(cycleNode.findall('rmsBonds')[0].text)
              except: data.append(None)
          progressGraph.addData(title="rmsBonds",  data=data)
          plotLine = plot.append('plotline',xcol=1,ycol=4,rightaxis='true',colour='red')

  def addProgressTable2(self, parent, internalId='SummaryTable'):
      progressTableDiv = parent.addDiv(style='border:0px solid black; height:250px; width:260px; float:left; margin-top:1px; margin-right:0px;overflow:auto;')
      xmlNodes = self.xmlnode.findall('.//RefmacRunning/Cycle')
      nCycles = len(xmlNodes)
#FIXME - XPATH LOGIC
      if len(xmlNodes)>0:
          selectString = ".//RefmacRunning/Cycle[1] "
          if len(xmlNodes)>2:
              selectString += " | .//RefmacRunning/Cycle[%d]" % (nCycles-1)
          if len(xmlNodes)>1:
              selectString += " | .//RefmacRunning/Cycle[%d]" % nCycles
          progressTable = progressTableDiv.addTable(select=selectString, style="height:250px; width:260px;float:left;", outputXml=self.outputXml, internalId=internalId)
          progressTable.addData(title="Cycle", select="number")
          progressTable.addData(title="R-factor", select="r_factor")
          progressTable.addData(title="R-free",   select="r_free",   expr="x if float(x)>=0.0 else '-'")

          rmsBonds = self.xmlnode.findall('.//RefmacRunning/Cycle/rmsBonds')
          if len(rmsBonds)> 0:
              cycleNodes = self.xmlnode.findall('.//RefmacRunning/Cycle')
              data = []
              for iCycleNode, cycleNode in enumerate(cycleNodes):
                  if iCycleNode == 0 or iCycleNode == nCycles-2 or iCycleNode == nCycles-1:
                      try: data.append(cycleNode.findall('rmsBonds')[0].text)
                      except: data.append('-')
              progressTable.addData(title="RMS Deviation", subtitle="Bond", data=data)

def test(xmlFile=None,jobId=None,reportFile=None):
    import sys,os
    try:
        text = open( xmlFile ).read()
        xmlnode = ET.fromstring( text )
    except:
        print('FAILED loading XML file:', kw['xmlFile'])
    if reportFile is None and xmlFile is not None:
        reportFile = os.path.join(os.path.split(xmlFile)[0],'report.html')
    r = molrep_pipe_report(xmlFile=xmlFile,jobId=jobId, xmlnode=xmlnode)
    r.as_html_file(reportFile)

if __name__ == "__main__":
  import sys
  test(xmlFile=sys.argv[1], jobId=sys.argv[2], reportFile=sys.argv[3])

