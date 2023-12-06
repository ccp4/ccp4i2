from report.CCP4ReportParser import *
import sys

#Use the RUN_TITLES from the script rather than the run code name
from . import xia2_run

class xia2_run_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'xia2_run'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, **kw)
        
        runSummaryFold = self.addFold(label="Xia2 runs", initiallyOpen=True)
        runTable = runSummaryFold.addTable(title="Data reduction summary", transpose=True)
        
        
        names = []
        titles = []
        spaceGroups = []
        cells = []
        statistics = {'meanIOverSigI':'I/sigI','completeness':'Completeness','rPimAllIPlusIMinus':'rPim(pooled Friedels)','multiplicity':'Multiplicity','rPimWithinIPlusIMinus':'rPim\n(separate Friedels)','rMeasAllIPlusIMinus':'rMeas(pooled Friedels)','resolutionLimitLow':'Lower res.','resolutionLimitHigh':'High res.','nTotalObservations':'nObs','anomalousMultiplicity':'Anom. Mult.','rMerge':'RMerge','anomalousCompleteness':'Anom Compl','nTotalUniqueObservations':'nUnique'}
        orderedKeys = ['resolutionLimitLow','resolutionLimitHigh','meanIOverSigI','completeness','anomalousCompleteness','rMerge','rPimAllIPlusIMinus','rPimWithinIPlusIMinus','rMeasAllIPlusIMinus','multiplicity','anomalousMultiplicity','nTotalObservations','nTotalUniqueObservations']
        statisticDict = {}
        for statistic in statistics:
            statisticDict[statistic] = []

        for runNode in self.xmlnode.findall('.//XIA2Run'):
            names.append(runNode.get('name'))
            titles.append(xia2_run.RUN_TITLES.get(names[-1],'Unknown'))
            try:
                spaceGroups.append (runNode.findall('AutoProcContainer/AutoProc/spaceGroup')[0].text)
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
                    overallValue = runNode.findall('AutoProcContainer/AutoProcScalingContainer/AutoProcScalingStatistics/'+statistic)[0].text
                    outerValue = runNode.findall('AutoProcContainer/AutoProcScalingContainer/AutoProcScalingStatistics/'+statistic)[2].text
                    statisticDict[statistic].append(overallValue+'('+outerValue+')')
                except:
                    statisticDict[statistic].append('Not known')

        runTable.addData(title='Run name', data=names)
        runTable.addData(title='Run mode', data=titles)
        runTable.addData(title='Space group', data=spaceGroups)
        runTable.addData(title='Unit cell', data=cells)
        for statistic in orderedKeys:
            runTable.addData(title=statistics[statistic], data=statisticDict[statistic])

        for runNode in self.xmlnode.findall('.//XIA2Run'):
            title = runNode.get('name')
            if runNode.get('name') in xia2_run.RUN_TITLES: title += ': '+xia2_run.RUN_TITLES[runNode.get('name')]
            runSummaryFold = self.addFold(label="Details of run "+title, initiallyOpen=False)
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




