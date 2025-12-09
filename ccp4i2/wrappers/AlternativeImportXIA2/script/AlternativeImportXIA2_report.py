from report.CCP4ReportParser import *

class AlternativeImportXIA2_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'AlternativeImportXIA2'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, **kw)
        
        runSummaryFold = self.addFold(label="Xia2 runs imported", initiallyOpen=True)
        runTable = runSummaryFold.addTable(title="Import table", transpose=True)
        
        
        names = []
        spaceGroups = []
        cells = []
        statistics = {'meanIOverSigI':'I/sigI','completeness':'Completeness','rPimAllIPlusIMinus':'rPim(pooled Friedels)','multiplicity':'Multiplicity','rPimWithinIPlusIMinus':'rPim\n(separate Friedels)','rMeasAllIPlusIMinus':'rMeas(pooled Friedels)','resolutionLimitLow':'Lower res.','resolutionLimitHigh':'High res.','nTotalObservations':'nObs','anomalousMultiplicity':'Anom. Mult.','rMerge':'RMerge','anomalousCompleteness':'Anom Compl','nTotalUniqueObservations':'nUnique'}
        orderedKeys = ['resolutionLimitLow','resolutionLimitHigh','meanIOverSigI','completeness','anomalousCompleteness','rMerge','rPimAllIPlusIMinus','rPimWithinIPlusIMinus','rMeasAllIPlusIMinus','multiplicity','anomalousMultiplicity','nTotalObservations','nTotalUniqueObservations']
        rounder = {'resolutionLimitLow':'{0:.2f}','resolutionLimitHigh':'{0:.2f}','meanIOverSigI':'{0:.2f}',
                   'completeness':'{0:.2f}','anomalousCompleteness':'{0:.2f}','rMerge':'{0:.3f}',
                   'rPimAllIPlusIMinus':'{0:.3f}','rPimWithinIPlusIMinus':'{0:.3f}','rMeasAllIPlusIMinus':'{0:.3f}',
                   'multiplicity':'{0:.2f}','anomalousMultiplicity':'{0:.2f}','nTotalObservations':'{0:.0f}',
                   'nTotalUniqueObservations':'{0:.0f}'}
        statisticDict = {}
        for statistic in statistics:
            statisticDict[statistic] = []

        for runNode in self.xmlnode.findall('.//XIA2Run'):
            names.append(runNode.get('name'))
            try:
                spaceGroups.append(runNode.findall('AutoProcContainer/AutoProc/spaceGroup')[0].text)
            except:
                spaceGroups.append ('Not found')
            try:
                cell = []
                for dimension in ['a','b','c','alpha','beta','gamma']:
                    value = float(runNode.findall('AutoProcContainer/AutoProc/refinedCell_'+dimension)[0].text)
                    cell += ['{0:.2f}'.format(value)]
                cells.append(', '.join(cell))
            except:
                cells.append ('Not found')
            for statistic in statistics:
                try:
                    ovVa = float(runNode.findall('AutoProcContainer/AutoProcScalingContainer/AutoProcScalingStatistics/'+statistic)[0].text)
                    overallValue = rounder[statistic].format(ovVa)
                    ouVa = float(runNode.findall('AutoProcContainer/AutoProcScalingContainer/AutoProcScalingStatistics/'+statistic)[2].text)
                    outerValue = rounder[statistic].format(ouVa)
                    statisticDict[statistic].append(overallValue+'  ('+outerValue+')')
                except:
                    statisticDict[statistic].append('Not known')

        runTable.addData(title='Run name', data=names)
        runTable.addData(title='Space group', data=spaceGroups)
        runTable.addData(title='Unit cell', data=cells)
        for statistic in orderedKeys:
            runTable.addData(title=statistics[statistic], data=statisticDict[statistic])

        for runNode in self.xmlnode.findall('.//XIA2Run'):
            runSummaryFold = self.addFold(label="Details of run "+runNode.get('name'), initiallyOpen=False)
            for program in ['POINTLESS','AIMLESS','TRUNCATE']:
                programNodes = runNode.findall(program)
                for programNode in programNodes:
                    summaryNodes = programNode.findall('CCP4Summary')
                    if len(summaryNodes) > 0:
                        runSummaryFold.addText(text='Summaries from  '+ program.lower(), style="font-size:125%;clear:both;")
                        runSummaryFold.append('<br/>')
                    for summaryNode in summaryNodes:
                        if len(summaryNode.text) < 2000:
                            font_color = "black"
                            if "ARNING" in summaryNode.text.upper(): font_color="orange"
                            runSummaryFold.addPre(text = summaryNode.text, font_color=font_color)

                    tableNodes = programNode.findall('CCP4ApplicationOutput/CCP4Table')
                    if len(tableNodes) > 0:
                        runSummaryFold.addText(text='Graphs from  '+ program.lower(), style="font-size:125%;clear:both;")
                        gallery = runSummaryFold.addObjectGallery(title="A gallery", tableWidth="32em", contentWidth="350px", height="300px", style="border:1px solid black;")
                    for tableNode in tableNodes:
                        graph = gallery.addFlotGraph( xmlnode=tableNode, title=tableNode.get("title"), style="width:340px; height:290px; border:0px solid black;" )
                        graph = graph.addPimpleData(xmlnode=tableNode)
                    runSummaryFold.addDiv(style="clear:both;")
                    runSummaryFold.append("<br/>")


