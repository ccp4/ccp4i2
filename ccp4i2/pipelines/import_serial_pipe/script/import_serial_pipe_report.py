from ....report.CCP4ReportParser import Report


class import_serial_pipe_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'import_serial_pipe'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus is None or jobStatus.lower() is 'nooutput': return
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None: parent = self

        try:
            parent.append("<br />")
            #cc = float(self.xmlnode.findall('overall/cc')[0].text)
            #overallStats.append("<p>Overall statistics</p>")
            overallStats = parent.addFold(label='Overall statistics', initiallyOpen=True)

            #tableDiv = overallStats.addDiv(style="height:250px;width:20em;float:left;border:0px;")
            tableDiv = overallStats.addDiv(style="width:20em;float:left;border:0px;")
            table = tableDiv.addTable(transpose=True)
            table.addData(title="Low resolution limit",select='overall/d_max')
            table.addData(title="High resolution limit",select='overall/d_min')
            table.addData(title="Number observed reflections",select='overall/n_obs')
            table.addData(title="Number unique reflections",select='overall/n_unique')
            table.addData(title="Completeness %",select='overall/completeness')
            table.addData(title="Multiplicity",select='overall/multiplicity')
            table.addData(title="Mean(I)",select='overall/I')
            table.addData(title="Mean((I)/sd(I))",select='overall/IsigI')
            table.addData(title="Half-set correlation CC(1/2)",select='overall/cc')
            table.addData(title="CC*",select='overall/CCstar')
            table.addData(title="Rsplit",select='overall/rsplit')
            clearingDiv = overallStats.addDiv(style="clear:both;")
            parent.append("<br />")


            binnedStats = parent.addFold(label='Statistics vs. resolution', initiallyOpen=True)

            graph = binnedStats.addFlotGraph(
                title="Multiplicity and completeness v resolution",
                xmlnode=self.xmlnode, select="binned/bin",
                style="width:450px;height:300px;margin:0 auto;float:left;border:0px;" )
            # graph.addData(title="Bin", select="n_bin" )
            # graph.addData(title="Resolution(A)", select="d_max")
            # graph.addData(title="Resolution(A)", select="d_min")
            graph.addData(title="Resolution(A)", select="one_over_d_min_sq") # 1 (x)
            graph.addData(title="#ObservedReflections", select="n_obs")      # 2
            graph.addData(title="#UniqueReflections", select="n_unique")     # 3
            graph.addData(title="Completeness(%)", select="completeness")    # 4
            graph.addData(title="Multiplicity", select="multiplicity")       # 5
            p = graph.addPlotObject()
            p.append( 'title', 'Multiplicity' )
            p.append( 'plottype', 'xy' )
            #p.append( 'xintegral', 'true' )
            p.append( 'xscale', 'oneoversqrt' )
            p.append( 'xlabel', 'Resolution (A)' )
            l = p.append('plotline',xcol=1,ycol=5)
            p = graph.addPlotObject()
            p.append( 'title', 'Completeness (%)' )
            p.append( 'plottype', 'xy' )
            p.append( 'xscale', 'oneoversqrt' )
            p.append( 'xlabel', 'Resolution (A)' )
            l = p.append('plotline',xcol=1,ycol=4)
            p = graph.addPlotObject()
            p.append( 'title', '#ObservedReflections' )
            p.append( 'plottype', 'xy' )
            p.append( 'xscale', 'oneoversqrt' )
            p.append( 'xlabel', 'Resolution (A)' )
            l = p.append('plotline',xcol=1,ycol=2)
            p = graph.addPlotObject()
            p.append( 'title', '#UniqueReflections' )
            p.append( 'plottype', 'xy' )
            p.append( 'xscale', 'oneoversqrt' )
            p.append( 'xlabel', 'Resolution (A)' )
            l = p.append('plotline',xcol=1,ycol=3)
            clearingDiv = binnedStats.addDiv(style="clear:both;")
            parent.append("<br />")

            graph = binnedStats.addFlotGraph(title="Data quality indicators v resolution", xmlnode=self.xmlnode, select="binned/bin", style="width:450px;height:300px;margin:0 auto;float:left;border:0px;" )
            graph.addData(title="Resolution(A)", select="one_over_d_min_sq") # 1 (x)
            # graph.addData(title="MeanI", select="I" )
            graph.addData(title="Mn(I/sdI)", select="IsigI" )           # 2
            if bool(self.xmlnode.findall('overall/cc')):
                graph.addData(title="CC1/2", select="cc" )                       # 3
                graph.addData(title="CC*", select="CCstar" )                     # 4
                graph.addData(title="Rsplit", select="rsplit" )                  # 5

                p = graph.addPlotObject()
                p.append( 'title', 'Half-set correlation CC(1/2)' )
                p.append( 'plottype', 'xy' )
                p.append( 'xscale', 'oneoversqrt' )
                p.append( 'xlabel', 'Resolution (A)' )
                l = p.append('plotline',xcol=1,ycol=3)
                l = p.append('plotline',xcol=1,ycol=4)
            p = graph.addPlotObject()
            p.append( 'title', 'Mean I/sigma(I)' )
            p.append( 'plottype', 'xy' )
            p.append( 'xscale', 'oneoversqrt' )
            p.append( 'xlabel', 'Resolution (A)' )
            l = p.append('plotline',xcol=1,ycol=2)
            if bool(self.xmlnode.findall('overall/cc')):
                p = graph.addPlotObject()
                p.append( 'title', 'Rsplit' )
                p.append( 'plottype', 'xy' )
                p.append( 'xscale', 'oneoversqrt' )
                p.append( 'xlabel', 'Resolution (A)' )
                l = p.append('plotline',xcol=1,ycol=5)
            clearingDiv = binnedStats.addDiv(style="clear:both;")
            parent.append("<br />")

        except Exception as e:
            parent.append("<p>THERE WAS A PROBLEM REPORTING THE STATISTICS: %s</p>"%(e,))
