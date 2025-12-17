from ccp4i2.report.CCP4ReportParser import *


class parrot_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'parrot'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None: parent = self

        results = parent.addResults()

        try:
          solc = float(self.xmlnode.findall('SolventContent/SolventContent')[0].text)
          v = 1.0-float(self.xmlnode.findall('SolventContent/NMols/SolventContent')[0].text)
          nseq = int(self.xmlnode.findall('SolventContent/NMolsFromSequence')[0].text)
          pseq = 1.0
          if nseq > 0: pseq = float(self.xmlnode.findall('SolventContent/NMols/Probability')[nseq-1].text)
          n = (1.0-solc)/v
          s = ""
          if pseq > 0.98:
            parent.append( "<p>The solvent content is set to %4.1f%% corresponding to %5.1f copies of the sequence in the asymmetric unit.</p>"%(100*solc,n) )
          else:
            parent.append( "<p>The solvent content is set to %4.1f%% corresponding to %5.1f copies of the sequence in the asymmetric unit. <i>This is only an estimate. Please check.</i></p>"%(100*solc,n) )
        except Exception as e:
          parent.append( "<p>THERE WAS A PROBLEM REPORTING THE SOLVENT CONTENT %s</p>"%(e,) )

        try:
          nncs = int(self.xmlnode.findall('NncsMR')[0].text)
          nmol = int(0.5*(1.0+(1.0+4.0*nncs)**0.5))
          if nncs > 0: parent.append( "<p>%d NCS operators were identified from the input MR or partial model, suggesting the presence of at least %d molecules.</p>"%(nncs,nmol) )
          else:        parent.append( "<p><b>Warning:</b> No NCS operators were determined from the input MR or partial model.</p>"%(nncs) )
        except Exception as e:
          print(e)

        try:
          nncs = int(self.xmlnode.findall('NncsHA')[0].text)
          nmol = int(0.5*(1.0+(1.0+4.0*nncs)**0.5))
          if nncs > 0: parent.append( "<p>%d NCS operators were identified from the input heavy atom model, suggesting the presence of at least %d molecules.</p>"%(nncs,nmol) )
          else:        parent.append( "<p><b>Warning:</b> No NCS operators were determined from the input heavy atom model.</p>"%(nncs) )
        except Exception as e:
          print(e)

        try:
          vol = float(self.xmlnode.findall('Final/NCSvolmean')[0].text)
          opc = self.xmlnode.findall('Final/Operators/NCScorrel')
          opv = self.xmlnode.findall('Final/Operators/NCSvolume')
          nncs = min(len(opc),len(opv))
          print(nncs)
          ngood = 0
          for i in range(nncs):
            print(i, opc[i].text, opv[i].text)
            if float(opc[i].text) > 0.4 and float(opv[i].text) > vol/2.0:
              ngood += 1
          if nncs > 0:
            if ngood == nncs:
              parent.append( "<p>Of the %d NCS operators, all had reasonable correlations and mask volumes at the end of the calculation.</p>"%(nncs) )
            elif ngood > 0:
              parent.append( "<p><b>Of the %d NCS operators, only %d had reasonable correlations and mask volumes at the end of the calculation.</b></p>"%(nncs,ngood) )
            else:
              parent.append( "<p><b>Of the %d NCS operators, none had reasonable correlations and mask volumes at the end of the calculation.</b></p>"%(nncs) )
        except Exception as e:
          print(e)

        try:
          final_fom = float(self.xmlnode.findall('Final/MeanFOM')[0].text)
          if final_fom > 0.600:
            parent.append( "<p>The final figure-of-merit is <b>%4.2f</b>, which suggests that the map is good enough for model building. However the figure-of-merit from density modification can be seriously overestimated.</p>"%(final_fom,) )
          elif final_fom > 0.400:
            parent.append( "<p>The final figure-of-merit is <b>%4.2f</b>, which suggests that the map is marginal for model building.</p>"%(final_fom,) )
          else:
            parent.append( "<p>The final figure-of-merit is <b>%4.2f</b>, which suggests that the map is too poor for model building.</p>"%(final_fom,) )
        except Exception as e:
          parent.append( "<p>THERE WAS A PROBLEM REPORTING THE FOM %s</p>"%(e,) )


        if len( self.xmlnode.findall('Cycles/Cycle') ) > 0:
          tableDiv = parent.addDiv(style="height:250px;width:20em;float:left;border:0px;")
          table = tableDiv.addTable(transpose=True)
          table.addData(title="Solvent content",select='SolventContent/SolventContent')
          table.addData(title="Initial FOM",select='Cycles/Cycle[1]/MeanFOM')
          table.addData(title="Final FOM",select='Cycles/Cycle[last()]/MeanFOM')

          graph = parent.addFlotGraph( title="Progress by cycle", xmlnode=self.xmlnode, select=".//Cycles/Cycle", style="width:450px;height:300px;margin:0 auto;float:left;border:0px;" )
          graph.addData(title="Cycle", select="Number" )
          graph.addData(title="Mean_FOM",               select="MeanFOM")
          graph.addData(title="Fcorrel<sub>work</sub>", select="Fcorrel")
          graph.addData(title="Fcorrel<sub>free</sub>", select="FreeFcorrel" )
          graph.addData(title="NCS_correlation_(mean)",     select="NCScormean" )
          graph.addData(title="NCS_volume_(mean)", select="NCSvolmean" )
          graph.addData(title="NCS_volume_(max)",  select="NCSvolmax" )
          graph.addData(title="NCS_volume_(min)",  select="NCSvolmin" )
          p = graph.addPlotObject()
          p.append( 'title', 'Reflection statistics by cycle' )
          p.append( 'plottype', 'xy' )
          p.append( 'xintegral', 'true' )
          p.append( 'xlabel', 'Cycle' )
          l = p.append('plotline',xcol=1,ycol=2)
          l = p.append('plotline',xcol=1,ycol=3)
          l = p.append('plotline',xcol=1,ycol=4)
          p = graph.addPlotObject()
          p.append( 'title', 'NCS correlation statistics by cycle' )
          p.append( 'plottype', 'xy' )
          p.append( 'xintegral', 'true' )
          p.append( 'xlabel', 'Cycle' )
          p.append( 'xrange', min='1' )
          p.append( 'yrange', min='0.0' )
          l = p.append('plotline',xcol=1,ycol=5)
          p = graph.addPlotObject()
          p.append( 'title', 'NCS volume statistics by cycle' )
          p.append( 'plottype', 'xy' )
          p.append( 'xintegral', 'true' )
          p.append( 'xlabel', 'Cycle' )
          p.append( 'xrange', min='1' )
          p.append( 'yrange', min='0.0' )
          l = p.append('plotline',xcol=1,ycol=6)
          l = p.append('plotline',xcol=1,ycol=7)
          l = p.append('plotline',xcol=1,ycol=8)
          clearingDiv = parent.addDiv(style="clear:both;")

        # add the solvent content analysis
        solContentFold = parent.addFold(label='Solvent content analysis')
        if len( self.xmlnode.findall('SolventContent/NMols') ) > 0:
          graph = solContentFold.addFlotGraph( title="Cell content analysis", xmlnode=self.xmlnode, select=".//SolventContent/NMols", style="width:450px; height:300px;margin: 0 auto;border:0px;" )
          graph.addData(title="Number_of_molecules", select="N")
          graph.addData(title="Solvent_fraction", select="SolventContent")
          graph.addData(title="Probability", select="Probability")
          p = graph.addPlotObject()
          p.append( 'title', 'Cell content analysis' )
          p.append( 'plottype', 'xy' )
          p.append( 'xintegral', 'true' )
          p.append( 'xlabel', 'Number of molecules' )
          p.append( 'yrange', min='0.0', max='1.0' )
          l = p.append('plotline',xcol=1,ycol=2)
          l = p.append('plotline',xcol=1,ycol=3)
          #l = p.append('histogram',col=3)
          #l.append('binwidth','.5')
        clearingDiv = parent.addDiv(style="clear:both;")

        # Add the rest of the graphs
        otherGraphsFold = parent.addFold(label='Other graphs from log file')
        # A GraphGroup is a group of graphs displayed in the same graph viewer widget
        graphgroup = otherGraphsFold.addFlotGraphGroup(style="width:450px; height:300px;margin: 0 auto;border:0px;")
        # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        
        # Loop over all Graph tables in the program output and add to the GraphGroup
        graphlist = self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table")
        for thisgraph in graphlist:
            graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)
        clearingDiv = parent.addDiv(style="clear:both;")
