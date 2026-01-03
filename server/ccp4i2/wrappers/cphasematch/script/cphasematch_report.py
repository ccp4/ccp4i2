from ccp4i2.report import Report


class cphasematch_report(Report):
    TASKNAME = 'cphasematch'
    
    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        results = self.addResults()
        results.append( 'The mean phase error is ' + self.xmlnode.findall ( ".//phaseError" )[0].text + " degrees" )
        results.append( 'However this may include phases which are weighted to zero. Weighted mean phase errors and map correlations may be more informative.' )

        # Summary table
        tableDiv = results.addDiv(style="height:30em;width:20em;float:left;border:0px;")
        table = tableDiv.addTable(select=".", transpose=True, id='table_1')
        for title,select in [["Mean phase error","phaseError"],[".. weighted by FOM1","weightedPhaseError1"],[".. weighted by FOM2","weightedPhaseError2"],["F-map correlation","reflectionFCorrelation"],["E-map correlation","reflectionECorrelation"]]:
          table.addData(title=title,select=select)

        # A GraphGroup is a group of graphs displayed in the same graph viewer widget
        graphgroup = results.addFlotGraphGroup(style="height:300px; width:450px; border:0px;")
        # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        
        # Loop over all Graph tables in the program output and add to the GraphGroup
        graphlist = self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table")
        for thisgraph in graphlist:
            graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)
