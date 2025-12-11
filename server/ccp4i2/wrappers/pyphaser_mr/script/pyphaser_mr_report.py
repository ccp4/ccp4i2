from __future__ import print_function
from ccp4i2.report.CCP4ReportParser import *

class pyphaser_mr_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'pyphaser_mr'
  # Flag that a 'Running' mode is supported
  RUNNING = True
  def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
    if self.errorReport().maxSeverity()>SEVERITY_WARNING:
      print('FAILED instantiating pyphaser_mr report generator')
      self.errorReport().report()
      return

    # 'nooutput' mode would be used by another report class that wanted
    # to use some method(s) from this class for its own report
    if jobStatus is not None and jobStatus.lower() == 'nooutput':
      return
    elif jobStatus.count('Running'):
      results = self.addResults()
      self.inputData(xmlnode,results)
    else:
      results = self.addResults()
      self.finalResults(xmlnode,results)

  def inputData(self, xmlnode, parent):

    textDiv = parent.addDiv()
    textDiv.addText( text="Phaser has analysed your input data.")

    textDiv = parent.addDiv()
    textDiv.addText( text="Input spacegroup is "+xmlnode.findall(".//SPG")[0].text )
    textDiv = parent.addDiv()
    textDiv.addText( text="Input unit cell is "+xmlnode.findall(".//CELL_A")[0].text+" "+xmlnode.findall(".//CELL_B")[0].text+" "+xmlnode.findall(".//CELL_C")[0].text+" "+xmlnode.findall(".//CELL_ALPHA")[0].text+" "+xmlnode.findall(".//CELL_BETA")[0].text+" "+xmlnode.findall(".//CELL_GAMMA")[0].text+" " )

  def finalResults(self, xmlnode, parent):

    textDiv = parent.addDiv()
    textDiv.addText( text='There were '+xmlnode.findall('.//numSolutions')[0].text )
    textDiv.addText( text=' solutions found in spacegroup: '+xmlnode.findall('.//Solution/SPG')[0].text )

    graph = parent.addGraph( title="Solutions from Phaser" , select=".//Solutions/Solution", style="width:300px; height:300px;" )
    graph.addData (title="Solution" , select="ISOL" )
    graph.addData (title="Components_found" , select="NCOMPONENTS" )
    graph.addData (title="Overall_LLG" , select="overallLLG" )
    graph.addData (title="Overall_TFZ" , select="overallTFZ" )

    p = graph.addPlotObject ()
    p.append('title','Overall LLG')
    p.append('plottype','xy')
    p.append('xintegral','true')
    l = p.append('plotline',xcol=1,ycol=3)
    l.append('label','Overall LLG')
    l.append('colour','blue')

    p = graph.addPlotObject ()
    p.append('title','Overall TFZ-equivalent')
    p.append('plottype','xy')
    p.append('xintegral','true')
    l = p.append('plotline',xcol=1,ycol=4)
    l.append('label','Overall TFZ')
    l.append('colour','red')

    p = graph.addPlotObject ()
    p.append('title','Number of components found')
    p.append('plottype','xy')
    p.append('xintegral','true')
    l = p.append('plotline',xcol=1,ycol=2)
    l.append('label','Components found')
    l.append('colour','green')

    pic = parent.addPicture(label="Final structure",sceneFile="$CCP4I2/wrappers/pyphaser_mr/script/pyphaser_mr_1.scene.xml")

        
