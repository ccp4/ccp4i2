from __future__ import print_function

from report.CCP4ReportParser import *
import sys
from xml.etree import ElementTree as ET
from numpy import sign


def isnumber(n):
    is_number = True
    try:
        num = float(n)
        # check for "nan" floats
        is_number = num == num   # or use `math.isnan(num)`
    except ValueError:
        is_number = False
    return is_number


class servalcat_xtal_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'servalcat_xtal'
    TASKTITLE = 'Servalcat - Macromolecular refinement'
    RUNNING = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        # 'nooutput' mode would be used by another report class that wanted
        # to use some method(s) from this class for its own report
        self.outputXml = jobStatus is not None and jobStatus.lower().count('running')
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        
        self.addDiv(style='clear:both;')

        if jobStatus.lower().count('running'):
            self.GraphPerCycle(parent=self)
        else:
            # Raid the servalcat report of finished jobs
            self.addSummary()

    def addSummary(self, xmlnode=None, parent=None, withTables=True):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        summaryFold = parent.addFold(label='Summary of refinement', brief='Summary', initiallyOpen=True)
        cycle_data = self.getCycleData()
        # self.addScrollableDownloadableTable1(cycle_data, parent=summaryFold)
        self.addTablePerCycle(cycle_data, parent=summaryFold, initialFinalOnly=True)
        self.addGraphPerCycle(parent=summaryFold)
        if withTables:
            clearingDiv = parent.addDiv(style="clear:both;")
            perCycleFold = parent.addFold(label='Per cycle statistics', brief='Per cycle', initiallyOpen=False)
            self.addTablePerCycle(cycle_data, parent=perCycleFold, initialFinalOnly=False)
            self.addGraphsVsResolution()
            self.addOutlierAnalysis()
    
    def addRunningProgressGraph_ToDelete(self, parent=None, xmlnode=None):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        # I *do not know* why this is needed
        if len(xmlnode.findall("cycle")) > 0:
            progressGraph = parent.addFlotGraph(title="Running servalcat", style="height:250px; width:400px;float:left;", outputXml=self.outputXml,internalId="SummaryGraph") # , select="Cycle")
            progressGraph.addData(title="Cycle", select=".//cycle/Ncyc")
            # progressGraph.addData(title="R_Factor", select="r_factor")
            # progressGraph.addData(title="R_Free",   select="r_free")
            if len(xmlnode.findall('.//cycle[last()]/data/summary/Rwork')) > 0:
                progressGraph.addData(title="R-work", select=".//cycle/data/summary/Rwork", expr="x if float(x)>=0.0 else ''")
                if len(xmlnode.findall('.//cycle[last()]/data/summary/Rfree')) > 0:
                    progressGraph.addData(title="R-free", select=".//cycle/data/summary/Rfree", expr="x if float(x)>=0.0 else '-'")
            elif len(xmlnode.findall('.//cycle[last()]/data/summary/R2work')) > 0:
                progressGraph.addData(title="R2work", select=".//cycle/data/summary/R2work", expr="x if float(x)>=0.0 else ''")
                if len(xmlnode.findall('.//cycle[last()]/data/summary/R2free')) > 0:
                    progressGraph.addData(title="R2free", select=".//cycle/data/summary/R2free", expr="x if float(x)>=0.0 else '-'")
            plot = progressGraph.addPlotObject()
            plot.append('title', 'Running servalcat R-values')
            plot.append('plottype', 'xy')
            plot.append('yrange', rightaxis='false')
            plot.append('xlabel', 'Cycle')
            plot.append('xintegral', 'true')
            plot.append('ylabel', 'R-factor')
            plot.append('rylabel', 'Geometry')
            for coordinate, colour in [(2,'blue'),(3,'green')]:
                plotLine = plot.append('plotline', xcol=1, ycol=coordinate, rightaxis='false', colour=colour)
            
            rmsBonds = self.xmlnode.findall('.//cycle/geom/summary/rmsd/Bond_distances_non_H')
            if len(rmsBonds) > 0:
                plot.append('yrange', rightaxis='true')
                # cycleNodes = self.xmlnode.findall('.//cycle/Ncyc')
                # data = []
                # for cycleNode in cycleNodes:
                #     try: data.append(cycleNode.findall('rmsBonds')[0].text)
                #     except: data.append(None)
                # progressGraph.addData(title="rmsBonds",  data=data)
                progressGraph.addData(title="RMSDbondx100", select=".//cycle/geom/summary/rmsd/Bond_distances_non_H")
                plotLine = plot.append('plotline', xcol=1, ycol=4, rightaxis='true', colour='red')

    def addGraphPerCycle(self, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        # Note that when I add the progressgraph, I have to ensure that the select is rooted in my own xmlnode
        progressGraph = parent.addFlotGraph(
            title="Refinement results",
            xmlnode=self.xmlnode,
            style="height:250px;width:400px;float:left;") # select = ".//Overall_stats/stats_vs_cycle/new_cycle", 
        progressGraph.addData(title="Cycle", select=".//cycle/Ncyc") # ycol=1
        progressGraph.addData(title="-LL", select=".//cycle/data/summary/minusLL", expr="x if float(x)>=0.0 else ''") # ycol=2
        if len(xmlnode.findall('.//cycle[last()]/data/summary/Rwork')) > 0:
            progressGraph.addData(title="Rwork", select=".//cycle/data/summary/Rwork", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="CCFwork_avg", select=".//cycle/data/summary/CCFworkavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
            if len(xmlnode.findall('.//cycle[last()]/data/summary/Rfree')) > 0:
                progressGraph.addData(title="Rfree", select=".//cycle/data/summary/Rfree", expr="x if float(x)>=0.0 else '-'")  # ycol=5
                progressGraph.addData(title="CCFfree_avg", select=".//cycle/data/summary/CCFfreeavg", expr="x if float(x)>=-1.0 else '-'")  # ycol=6
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/R2work')) > 0:
            progressGraph.addData(title="R2work", select=".//cycle/data/summary/R2work", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="CCIwork_avg", select=".//cycle/data/summary/CCIworkavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
            if len(xmlnode.findall('.//cycle[last()]/data/summary/R2free')) > 0:
                progressGraph.addData(title="R2free", select=".//cycle/data/summary/R2free", expr="x if float(x)>=0.0 else '-'")  # ycol=5
                progressGraph.addData(title="CCIfree_avg", select=".//cycle/data/summary/CCIfreeavg", expr="x if float(x)>=-1.0 else '-'") # ycol=6
        # For lines that don''t have a value for each point, the trick is to replace missing values with '-'.
        # Out of refmac, they are flagged with a value of -999.
        # progressGraph.addData(title="RMSDbondx100", select=".//cycle/geom/summary/rmsd/Bond_distances_non_H", expr="str(100.*float(x)) if float(x)>=0.0 else '-'")
        # progressGraph.addData(title="RMSDangle", select=".//cycle/geom/summary/rmsd/Bond_angles_non_H", expr="x if float(x)>=0.0 else '-'")

        #progressGraph = parent.addFlotGraph( title="Refinement results", xmlnode=self.xmlnode, select = ".//Overall_stats/stats_vs_cycle/new_cycle",style="height:250px; width:400px;float:left;")
        #progressGraph.addData(title="Cycle",   select=".//cycle")
        #progressGraph.addData(title="R-free",   select=".//r_free", expr="x if float(x)>=0.0 else '-'")
        #progressGraph.addData(title="R-factor", select=".//r_factor", expr="x if float(x)>=0.0 else ''")
        ##For  lines that dont have a value for each point, the trick is to replace missing values with '-'.
        ##Out of refmac, they are flagged with a value of -999.
        #progressGraph.addData(title="rmsBONDx100", select=".//rmsBOND", expr="str(100.*float(x)) if float(x)>=0.0 else '-'")
        #progressGraph.addData(title="rmsANGLE", select=".//rmsANGLE", expr="x if float(x)>=0.0 else '-'")
        
        plotR = progressGraph.addPlotObject()
        plotR.append('title', 'R-values')
        plotR.append('plottype', 'xy')
        # plotR.append('yrange', rightaxis='false')
        plotR.append('xlabel', 'Cycle')
        plotR.append('ylabel', 'R-value')
        # plotR.append('rylabel', 'Geometry')  ### NOT VISIBLE !?
        plotR.append('xintegral', 'true')
        # for coordinate, colour in [(2,'blue'),(3,'green')]:
        #     plotLine = plotR.append('plotline', xcol=1, ycol=coordinate, rightaxis='false', colour=colour)
        plotLine = plotR.append('plotline', xcol=1, ycol=3)
        plotLine.append('colour', 'orange')
        plotLine.append('symbolsize', '0')
        plotLine = plotR.append('plotline', xcol=1, ycol=5)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')
        # plot.append('yrange', rightaxis='true')
        # for coordinate, colour in [(4,'red'),(5,'purple')]:
        #     plotLine = plot.append('plotline', xcol=1, ycol=coordinate, rightaxis='true', colour=colour)
        plotCC = progressGraph.addPlotObject()
        plotCC.append('title', 'Correlations')
        plotCC.append('plottype', 'xy')
        plotCC.append('xlabel', 'Cycle')
        plotCC.append('ylabel', 'Correlation')
        plotCC.append('xintegral', 'true')
        plotLine = plotCC.append('plotline', xcol=1, ycol=4)
        plotLine.append('colour', 'orange')
        plotLine.append('symbolsize', '0')
        plotLine = plotCC.append('plotline', xcol=1, ycol=6)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')

        plotLL = progressGraph.addPlotObject()
        plotLL.append('title', '-LL')
        plotLL.append('plottype', 'xy')
        plotLL.append('xlabel', 'Cycle')
        plotLL.append('ylabel', '-LL')
        plotLL.append('xintegral', 'true')
        plotLine.append('colour', 'orange')
        plotLine = plotLL.append('plotline', xcol=1, ycol=2)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')

        clearingDiv = parent.addDiv(style="float:left;width=100px;")
        clearingDiv.append('&#160;&#160;')

        if len(xmlnode.findall('.//cycle/geom/summary/rmsd/Bond_distances_non_H')) > 0:
            progressGraph2 = parent.addFlotGraph(
                title="Deviation from standard geometry",
                xmlnode=self.xmlnode,
                style="height:250px;width:400px;float:left;")
            # TO DO Could be improved - colour of RMSDangle and right y-axis label
            progressGraph2.addData(title="Cycle", select=".//cycle/Ncyc")
            progressGraph2.addData(title="RMSD_bond", select=".//cycle/geom/summary/rmsd/Bond_distances_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSD_angle", select=".//cycle/geom/summary/rmsd/Bond_angles_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSZ_bond", select=".//cycle/geom/summary/rmsZ/Bond_distances_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSZ_angle", select=".//cycle/geom/summary/rmsZ/Bond_angles_non_H", expr="x if float(x)>=0.0 else '-'")
            plotRmsd = progressGraph2.addPlotObject()
            plotRmsd.append('title', 'RMS Deviations')
            plotRmsd.append('plottype', 'xy')
            plotRmsd.append('yrange', rightaxis='false')
            plotRmsd.append('xlabel', 'Cycle')
            plotRmsd.append('ylabel', ' ')
            plotRmsd.append('xintegral', 'true')
            plotLine = plotRmsd.append('plotline', xcol=1, ycol=2, rightaxis='false')
            plotRmsd.append('yrange', rightaxis='true')
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')
            plotLine = plotRmsd.append('plotline', xcol=1, ycol=3, rightaxis='true')
            plotLine.append('colour', 'red')
            plotLine.append('symbolsize', '0')

            plotRmsz = progressGraph2.addPlotObject()
            plotRmsz.append('title', 'RMS Z-scores')
            plotRmsz.append('plottype', 'xy')
            plotRmsz.append('xlabel', '')
            plotRmsz.append('ylabel', '')
            plotRmsz.append('xintegral', 'true')
            plotLine = plotRmsz.append('plotline', xcol=1, ycol=4, rightaxis='false')
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')
            plotLine = plotRmsz.append('plotline', xcol=1, ycol=5, rightaxis='false')
            plotLine.append('colour', 'red')
            plotLine.append('symbolsize', '0')

    def getCycleData(self, xmlnode=None):
        if xmlnode is None: xmlnode = self.xmlnode
        R2WorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/R2work')
        #R2FreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/R2free')
        #CCIWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCIworkavg')
        #CCIFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCIfreeavg')
        #RWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/Rwork')
        #RFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/Rfree')
        #CCFWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCFworkavg')
        #CCFFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCFfreeavg')
        #RMSBondsNodes = xmlnode.findall('.//cycle[last()]/geom/summary/rmsd/Bond_distances_non_H')
        #RMSAnglesNodes = xmlnode.findall('.//cycle[last()]/geom/summary/rmsd/Bond_angles_non_H')
        all_cycles = xmlnode.findall('.//cycle')
        ncyc = len(all_cycles)
        cycle_data = {'mode':['-']*ncyc,
                      'cycle':['-']*ncyc,
                      '-LL':['-']*ncyc,
                      'Rwork':['-']*ncyc,
                      'Rfree':['-']*ncyc,
                      'CCFworkavg':['-']*ncyc,
                      'CCFfreeavg':['-']*ncyc,
                      'R2work':['-']*ncyc,
                      'R2free':['-']*ncyc,
                      'CCIworkavg':['-']*ncyc,
                      'CCIfreeavg':['-']*ncyc,
                      'rmsBOND':['-']*ncyc,
                      'rmsANGLE':['-']*ncyc,
                      'rmsCHIRAL':['-']*ncyc,
                      'zBOND':['-']*ncyc,
                      'zANGLE':['-']*ncyc,
                      'zCHIRAL':['-']*ncyc}
        idx = 0
        for cycle in all_cycles:
            #try:  # TO DO
            #   if len(xmlnode.findall("RigidMode"))>0:
            #      cycle_data['mode'][idx] = 'Rigid'
            #   elif float(cycle.findall('rmsBOND')[0].text) >= 0.0:
            #      cycle_data['mode'][idx] = 'Restr'
            #   elif str(float(cycle.findall('rmsBOND')[0].text)) == '-999.0' and len(xmlnode.findall("TLSMode"))>0:
            #      cycle_data['mode'][idx] = 'TLS'
            #   else: raise
            #except:
            #   print("*** ERROR - UNABLE TO INTERPRET REFMAC5 LOG FILE TO CONSTRUCT REFINEMENT TABLES - PLEASE CONTACT CCP4 WITH THIS ERROR (ccp4@ccp4.ac.uk) ***")
            #   return None
            cycle_data['mode'][idx] = 'Restr'
            try: cycle_data['cycle'][idx] = str(int(cycle.findall('Ncyc')[0].text))
            except: pass
            try: cycle_data['-LL'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/-LL')[0].text))
            except: pass
            if len(R2WorkNodes) > 0:  # Refinement against intensities
                try: cycle_data['R2work'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R2work')[0].text))
                except: pass
                try: cycle_data['R2free'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R2free')[0].text))
                except: pass
                try: cycle_data['CCIworkavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCIworkavg')[0].text))
                except: pass
                try: cycle_data['CCIfreeavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCIfreeavg')[0].text))
                except: pass
            else:                     # Refinement against amplitudes
                try: cycle_data['Rwork'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/Rwork')[0].text))
                except: pass
                try: cycle_data['Rfree'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/Rfree')[0].text))
                except: pass
                try: cycle_data['CCFworkavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCFworkavg')[0].text))
                except: pass
                try: cycle_data['CCFfreeavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCFfreeavg')[0].text))
                except: pass
            if cycle_data['mode'][idx] == 'Restr':
                try: cycle_data['rmsBOND'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsd/Bond_distances_non_H')[0].text))
                except: pass
                try: cycle_data['rmsANGLE'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsd/Bond_angles_non_H')[0].text))
                except: pass
                try: cycle_data['rmsCHIRAL'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsd/Chiral_centres')[0].text))
                except: pass
                try: cycle_data['zBOND'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsZ/Bond_distances_non_H')[0].text))
                except: pass
                try: cycle_data['zANGLE'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsZ/Bond_angles_non_H')[0].text))
                except: pass
                try: cycle_data['zCHIRAL'][idx] = "{:.4f}".format(float(cycle.findall('geom/summary/rmsZ/Chiral_centres')[0].text))
                except: pass
            idx += 1
        return cycle_data

    def addTablePerCycle(self, cycle_data, parent=None, initialFinalOnly=False):
        if parent is None: parent = self
        # if xmlnode is None: xmlnode = self.xmlnode
        clearingDiv = parent.addDiv(style="clear:both;")
        if initialFinalOnly:
            # Pick only data for the first and last cycle
            cycle_data_sel = {}
            for key, values in cycle_data.items():
                cycle_data_sel[key] = [values[0], values[-1]]
            cycle_data_sel['cycle'] = ["Initial", "Final"]
        else:
            cycle_data_sel = cycle_data

        RigidIdx = [idx for idx, val in enumerate([mode == 'Rigid' for mode in cycle_data_sel['mode']]) if val]
        TLSIdx = [idx for idx, val in enumerate([mode == 'TLS' for mode in cycle_data_sel['mode']]) if val]
        RestrIdx = [idx for idx, val in enumerate([mode == 'Restr' for mode in cycle_data_sel['mode']]) if val]

        """
        # Only display first and last in the summary tables
        if len(RigidIdx)>2: RigidIdx = [RigidIdx[0],RigidIdx[-1]]
        if len(TLSIdx)>2: TLSIdx = [TLSIdx[0],TLSIdx[-1]]
        if len(RestrIdx)>2: RestrIdx = [RestrIdx[0],RestrIdx[-1]]
        start_end = ['Initial', 'Final']

        rigid_data = {}
        tls_data = {}
        restr_data = {}
        for key,val in cycle_data_local.items():
           if len(RigidIdx)>0: rigid_data[key] = [val[i] for i in RigidIdx]
           if len(TLSIdx)>0: tls_data[key] = [val[i] for i in TLSIdx]
           if len(RestrIdx)>0: restr_data[key] = [val[i] for i in RestrIdx]

        if len(RigidIdx)>0:
           RigidDiv = perCycleFold.addDiv(style="float:left;")
           RigidDiv.addText(text="Rigid body refinement")
           RigidTable = RigidDiv.addTable()
           RigidTable.addData(title="",data=start_end)
           RigidTable.addData(title="Cycle", data=rigid_data['cycle'])
           RigidTable.addData(title="R-factor", data=rigid_data['r_factor'])
           RigidTable.addData(title="R-free", data=rigid_data['r_free'])
           #RigidTable.addData(title="FSCAv", data=rigid_data['fscAver'])
           #RigidTable.addData(title="FSCAv-free", data=rigid_data['fscAverFree'])

        if len(TLSIdx)>0:
           TLSDiv = perCycleFold.addDiv(style="float:left; width:250px;")
           TLSDiv.addText(text="TLS refinement")
           TLSTable = TLSDiv.addTable()
           TLSTable.addData(title="",data=start_end)
           TLSTable.addData(title="Cycle", data=tls_data['cycle'])
           TLSTable.addData(title="R-factor", data=tls_data['r_factor'])
           TLSTable.addData(title="R-free", data=tls_data['r_free'])
           #TLSTable.addData(title="FSCAv", data=tls_data['fscAver'])
           #TLSTable.addData(title="FSCAv-free", data=tls_data['fscAverFree'])

        if len(RestrIdx)>0:
           RestrDiv = perCycleFold.addDiv(style="float:left;")
           RestrDiv.addText(text="Full atom refinement - XYZ, B (Occ)")
           RestrTable = RestrDiv.addTable()
           RestrTable.addData(title="",data=start_end)
           RestrTable.addData(title="Cycle", data=restr_data['cycle'])
           RestrTable.addData(title="R-factor", data=restr_data['r_factor'])
           RestrTable.addData(title="R-free", data=restr_data['r_free'])
           #RestrTable.addData(title="FSCAv", data=restr_data['fscAver'])
           #RestrTable.addData(title="FSCAv-free", data=restr_data['fscAverFree'])
           RestrTable.addData(title="RMSD (bond/angle/chiral)", subtitle="Bond", data=restr_data['rmsBOND'])
           RestrTable.addData(subtitle="Angle", data=restr_data['rmsANGLE'])
           RestrTable.addData(subtitle="Chiral", data=restr_data['rmsCHIRAL'])
        """

        clearingDiv = parent.addDiv(style="clear:both;")
        # detailFold = parent.addFold(label='All cycles', initiallyOpen=False, brief='All cycles')
        # fullTable = detailFold.addTable()
        fullTable = None
        fullTable = parent.addTable()
        if len(TLSIdx) > 0 and len(RestrIdx) > 0:
           mode = ['TLS']*len(TLSIdx) + ['Full Atom']*len(RestrIdx)
        fullTable.addData(title="Cycle", data=cycle_data_sel['cycle'])
        if isnumber(cycle_data_sel['R2work'][-1]):  # Refinement against intensities
            fullTable.addData(title="R2work", data=cycle_data_sel['R2work'])
            if isnumber(cycle_data_sel['R2free'][-1]):
                fullTable.addData(title="R2free", data=cycle_data_sel['R2free'])
            fullTable.addData(title="CCIwork_avg", data=cycle_data_sel['CCIworkavg'])
            if isnumber(cycle_data_sel['CCIfreeavg'][-1]):
                fullTable.addData(title="CCIfree_avg", data=cycle_data_sel['CCIfreeavg'])
        else:                                # Refinement against amplitudes
            fullTable.addData(title="Rwork", data=cycle_data_sel['Rwork'])
            if isnumber(cycle_data_sel['Rfree'][-1]):
                fullTable.addData(title="Rfree", data=cycle_data_sel['Rfree'])
            fullTable.addData(title="CCFwork_avg", data=cycle_data_sel['CCFworkavg'])
            if isnumber(cycle_data_sel['CCFfreeavg'][-1]):
                fullTable.addData(title="CCFfree_avg", data=cycle_data_sel['CCFfreeavg'])
        if isnumber(cycle_data_sel['rmsANGLE'][-1]):
            fullTable.addData(title="RMSD (bond/angle/chiral)", subtitle="Bond", data=cycle_data_sel['rmsBOND'])
            fullTable.addData(subtitle="Angle", data=cycle_data_sel['rmsANGLE'])
            fullTable.addData(subtitle="Chiral", data=cycle_data_sel['rmsCHIRAL'])
            fullTable.addData(title="RMSZ (bond/angle/chiral)", subtitle="Bond", data=cycle_data_sel['zBOND'])
            fullTable.addData(subtitle="Angle", data=cycle_data_sel['zANGLE'])
            fullTable.addData(subtitle="Chiral", data=cycle_data_sel['zCHIRAL'])

    def addSymmetryAnalysis(self, parent=None):
        if parent is None:
            parent=self
        equivalenceNodes = self.xmlnode.findall('.//NCS/equivalence')
        if len(equivalenceNodes) == 0: return
        symmetryFold = parent.addFold(label='Non-crystallographic symmetry identified by Refmac',brief='Non-X symm')
        table = symmetryFold.addTable(select='.//NCS/equivalence')
            
        table.addData(title='',data=['Selection 1','Selection 2','Residues aligned','Score','RMS','Average(RMS loc)'])
        for equivalenceNode in equivalenceNodes:
            table.addData(title='',data=[equivalenceNode.get(property) for property in ['selection1','selection2','nEquivalent','score','rms','ave_rmsLoc']])
        return

    def addTwinningAnalysis(self, parent=None):
        if parent is None:
            parent=self
        try: twinningNode = self.xmlnode.findall('.//Twinning')[-1]
        except: return
        
        twinningFold = parent.addFold(label='Refmac analysis of twinning',brief='Twinning')
        criteriaText = ""
        try:
            rMergeNode = twinningNode.findall('RmergeFilter')[-1]
            criteriaText += ('Twins included if Rmerge < %s.'%rMergeNode.text)
        except: pass
        try:
            fractionNode = twinningNode.findall('.//Twinning/FractionFilter')[-1]
            criteriaText += ('Twins included if Fraction > %s.'%fractionNode.text)
        except: pass
        if len(criteriaText)>0:
            textDiv = twinningFold.addDiv()
            newText = textDiv.addText(text=criteriaText)
        
        try:
           rMergeNodes = twinningNode.findall('TwinningSymmetryOperator')
           nrmerge = len(rMergeNodes)
           rmerge_data = {'operator':['-']*nrmerge,'rmerge':['-']*nrmerge}
           
           idx = 0
           for operator in rMergeNodes:
              try: rmerge_data['operator'][idx] = str(operator.findall('SymmetryOperator')[0].text)
              except: pass
              try: rmerge_data['rmerge'][idx] = str(operator.findall('Rmerge')[0].text)
              except: pass
              idx = idx+1
                 
           RMergeTable = twinningFold.addTable()
           RMergeTable.addData(title="Symmetry Operator",data=rmerge_data['operator'])
           RMergeTable.addData(title="RMerge",data=rmerge_data['rmerge'])
        except: pass
        
        try:
           fractionNodes = twinningNode.findall('TwinOperator')
           nfraction = len(fractionNodes)
           fraction_data = {'operator':['-']*nfraction,'fraction':['-']*nfraction}
           
           idx = 0
           for operator in fractionNodes:
              try: fraction_data['operator'][idx] = str(operator.findall('SymmetryOperator')[0].text)
              except: pass
              try: fraction_data['fraction'][idx] = str(operator.findall('Fraction')[0].text)
              except: pass
              idx = idx+1
           
           FractionTable = twinningFold.addTable()
           FractionTable.addData(title="Twinning Operator",data=fraction_data['operator'])
           FractionTable.addData(title="Fraction",data=fraction_data['fraction'])
        except: pass
        
#        table = twinningFold.addTable(select='//Twinning[last()]/TwinningSymmetryOperator',style="width:250px;float:left;",internalId="TwinningOperatorTable",outputXml=self.outputXml)
#        for title,select in  [[ "Symmetry operator" ,"SymmetryOperator" ], [ "Rmerge"  , "Rmerge" ]]:
#            table.addData(title=title,select=select)
#        
#        table = twinningFold.addTable(select='//Twinning[last()]/TwinOperator',style="width:250px;float:left;",internalId="TwinningFractionTable",outputXml=self.outputXml)
#        for title,select in  [[ "Twinning operator" ,"SymmetryOperator" ], [ "Fraction"  , "Fraction" ]]:
#            table.addData(title=title,select=select)

        clearingDiv = parent.addDiv(style="clear:both;")
        return

    def addGraphsVsResolution(self, parent=None, xmlnode=None, internalIdPrefix=''):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        reportFold = parent.addFold(label='Statistics vs. resolution', brief='Other')

        gallery = reportFold.addObjectGallery(
            height='310px', contentWidth='420px', tableWidth='360px', style='float:left;width:800px;')
        galleryGraphStyle = "width:410px;height:290px;"

        graphRtitle = "R-values"
        graphR = gallery.addFlotGraph(
            xmlnode=xmlnode,
            title=graphRtitle,
            internalId=graphRtitle,
            outputXml=self.outputXml,
            label=graphRtitle,
            style=galleryGraphStyle,
            initiallyDrawn=True)
        graphR.addData(title="Resolution(A)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
        if len(xmlnode.findall('.//cycle[last()]/data/binned/Rwork')) > 0:
            graphR.addData(title="Rwork", select=".//cycle[last()]/data/binned/./Rwork")
            if len(xmlnode.findall('.//cycle[last()]/data/binned/Rfree')) > 0:
                graphR.addData(title="Rfree", select=".//cycle[last()]/data/binned/./Rfree")
        else:
            graphR.addData(title="R2work", select=".//cycle[last()]/data/binned/./R2work")
            if len(xmlnode.findall('.//cycle[last()]/data/binned/R2free')) > 0:
                graphR.addData(title="R2free", select=".//cycle[last()]/data/binned/./R2free")
        plotR = graphR.addPlotObject()
        plotR.append('title', graphRtitle)
        plotR.append('plottype', 'xy')
        plotR.append('xlabel', 'Resolution (A)')
        plotR.append('ylabel', 'R-value')
        plotR.append('xscale', 'oneoversqrt')
        plotR.append('legendposition', x=0, y=0)
        plotR.append('line', x1=0, x2=1, y1=0.42, y2=0.42, linecolour="red", linestyle="--")
        plotR.append('line', x1=0, x2=1, y1=0.63, y2=0.58, linecolour="red", linestyle="--")
        plotLine = plotR.append('plotline', xcol=1, ycol=2)
        plotLine.append('colour', 'orange')
        plotLine.append('symbolsize', '0')
        plotLine = plotR.append('plotline', xcol=1, ycol=3)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')

        graphCCtitle = "Correlations"
        graphCC = gallery.addFlotGraph(
            xmlnode=xmlnode,
            title=graphCCtitle,
            internalId=graphCCtitle,
            outputXml=self.outputXml,
            label=graphCCtitle,
            style=galleryGraphStyle)
        graphCC.addData(title="Resolution(A)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
        if len(xmlnode.findall('.//cycle[last()]/data/binned/CCFwork')) > 0:
            graphCC.addData(title="CCFwork", select=".//cycle[last()]/data/binned/./CCFwork")
            if len(xmlnode.findall('.//cycle[last()]/data/binned/CCFfree')) > 0:
                graphCC.addData(title="CCFfree", select=".//cycle[last()]/data/binned/./CCFfree")
        else:
            graphCC.addData(title="CCIwork", select=".//cycle[last()]/data/binned/./CCIwork")
            if len(xmlnode.findall('.//cycle[last()]/data/binned/CCIfree')) > 0:
                graphCC.addData(title="CCIfree", select=".//cycle[last()]/data/binned/./CCIfree")
        plotCC = graphCC.addPlotObject()
        plotCC.append('title', graphCCtitle)
        plotCC.append('plottype', 'xy')
        plotCC.append('xlabel', 'Resolution (A)')
        plotCC.append('ylabel', 'Correlation')
        plotCC.append('xscale', 'oneoversqrt')
        plotLine = plotCC.append('plotline', xcol=1, ycol=2)
        plotLine.append('colour', 'orange')
        plotLine.append('symbolsize', '0')
        plotLine = plotCC.append('plotline', xcol=1, ycol=3)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')

        graphNtitle = "Number of reflections"
        graphN = gallery.addFlotGraph(
            xmlnode=xmlnode,
            title=graphNtitle,
            internalId=graphNtitle,
            outputXml=self.outputXml,
            label=graphNtitle,
            style=galleryGraphStyle)
        graphN.addData(title="Resolution(A)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
        graphN.addData(title="Nobs", select=".//cycle[last()]/data/binned/./n_obs")
        graphN.addData(title="Nwork", select=".//cycle[last()]/data/binned/./n_work")
        graphN.addData(title="Nfree", select=".//cycle[last()]/data/binned/./n_free")
        plotN = graphN.addPlotObject()
        plotN.append('title', graphNtitle)
        plotN.append('plottype', 'xy')
        plotN.append('xlabel', 'Resolution (A)')
        # plotN.append('ylabel', '')
        plotN.append('xscale', 'oneoversqrt')
        plotLine = plotN.append('plotline', xcol=1, ycol=2)
        plotLine.append('colour', 'orange')
        plotLine.append('symbolsize', '0')
        plotLine = plotN.append('plotline', xcol=1, ycol=3)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')
        plotN.append('yrange', rightaxis='true')
        plotLine = plotN.append('plotline', xcol=1, ycol=4) # , rightaxis='true')
        plotLine.append('colour', 'red')
        plotLine.append('symbolsize', '0')

        graphDtitle = "Mean |D0*FC0| and |D1*FCbulk|"
        graphD = gallery.addFlotGraph(
            xmlnode=xmlnode,
            title=graphDtitle,
            internalId=graphDtitle,
            outputXml=self.outputXml,
            label=graphDtitle,
            style=galleryGraphStyle)
        graphD.addData(title="Resolution(A)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
        graphD.addData(title="Mean|D0*FC0|", select=".//cycle[last()]/data/binned/./MnD0FC0")
        graphD.addData(title="Mean|D1*FCbulk|", select=".//cycle[last()]/data/binned/./MnD1FCbulk")
        plotD = graphD.addPlotObject()
        plotD.append('title', graphDtitle)
        plotD.append('plottype', 'xy')
        plotD.append('xlabel', 'Resolution (A)')
        # plotD.append('ylabel', '')
        plotD.append('xscale', 'oneoversqrt')
        plotLine = plotD.append('plotline', xcol=1, ycol=2)
        plotLine.append('colour', 'blue')
        plotLine.append('symbolsize', '0')
        plotD.append('yrange', rightaxis='true')
        plotLine = plotD.append('plotline', xcol=1, ycol=3, rightaxis='true')
        plotLine.append('colour', 'red')
        plotLine.append('symbolsize', '0')

        clearingDiv = parent.addDiv(style="clear:both;")

        #reportFold = parent.addFold(label='Picture', brief='Other')
        #reportFold = parent.addFold(label='Outliers identified by servalcat', brief='Other')
        #reportFold = parent.addFold(label='Validation', brief='Other')
        #reportFold = parent.addFold(label='Verdict', brief='Other')

    def addSmartieGraphs_old(self, parent=None, internalIdPrefix=''):
        if parent is None:
            parent=self
        reportFold = parent.addFold(label='Other plots from log file',brief='Other')
        
        """
        #MN Here explicitly correct the absence of annotation of xscale as oneoversqrt in plots versus resln
        for fixmeNode in self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table[contains(@title,'resln')]/plot"):
            fixmeNode.append(ET.fromstring('<xscale>oneoversqrt</xscale>'))
        """
        #SJM - A python xml library version of above.    
        for tableNode in self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table"):
            if "title" in tableNode.attrib and "resln" in tableNode.attrib["title"]:
                plots = tableNode.findall("./plot")
                for plot in plots:
                    plot.append(ET.fromstring('<xscale>oneoversqrt</xscale>'))
        graphTableList = self.xmlnode.findall('SmartieGraphs/CCP4ApplicationOutput/CCP4Table')
        gallery = reportFold.addObjectGallery(height='300px',contentWidth='420px',tableWidth='360px',style='float:left;width:800px;')
        
        isFirstGraph = True
        for iGraph, graphTableNode in enumerate(graphTableList):
            graph = gallery.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title"), internalId=internalIdPrefix+'SmartiePlot'+str(iGraph), outputXml=self.outputXml, label=graphTableNode.get("title"),style="width:410px;height:290px;", initiallyDrawn=isFirstGraph )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
            isFirstGraph = False

    def addOutlierAnalysis(self, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        # All possible keys are for now: 'bond', 'angle', 'torsion', 'chir', 'plane', 'staca', 'stacd', 'vdw'
        # ADP will be added
        # 'per_atom' not to be used
        outlierFold = parent.addFold(label="Outliers identified by Servalcat", brief='Outliers')
        #outlierNodes = self.xmlnode.findall('.//cycle[last()]/geom/outliers/outliers')
        #outBond = outlierNodes.findall('bond')
        outBond = xmlnode.findall('.//cycle[last()]/geom/outliers/bond')
        outAngle = xmlnode.findall('.//cycle[last()]/geom/outliers/angle')
        outTorsion = xmlnode.findall('.//cycle[last()]/geom/outliers/torsion')
        outChir = xmlnode.findall('.//cycle[last()]/geom/outliers/chir')
        outPlane = xmlnode.findall('.//cycle[last()]/geom/outliers/plane')
        outStaca = xmlnode.findall('.//cycle[last()]/geom/outliers/staca')
        outStacd = xmlnode.findall('.//cycle[last()]/geom/outliers/stacd')
        outVdw = xmlnode.findall('.//cycle[last()]/geom/outliers/vdw')
        # outPerAtom = xmlnode.findall('.//cycle[last()]/geom/outliers/per_atom')

        if len(outBond) > 0:
            n_outliers = len(outBond)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Bond length outliers:")
            outData = {'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers,
                       'percent': ["-"]*n_outliers,
                       'type': ["-"]*n_outliers,
                       'note': ["-"]*n_outliers}
            for i, outlier in enumerate(outBond):
                try: outData['atom1'][i] = str(outlier.findall('atom1')[0].text)
                except: outData['atom1'][i] = '-'
                try: outData['atom2'][i] = str(outlier.findall('atom2')[0].text)
                except: outData['atom2'][i] = '-'
                try:
                    outType = int(outlier.findall('type')[0].text)
                    #if outType == 0:
                    #    outData['note'][i] = "Bond type restraint"
                    #elif outType == 1:
                    #    outData['note'][i] = "Additional restraint"
                    if outType == 2:
                        outData['note'][i] = "External restraint"
                    outData['type'][i] = -outType
                except:
                    outData['type'][i] = '-'
                    outData['note'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    # difference = | value - ideal |
                    difference = abs(value - ideal)
                    outData['difference'][i] = "{:.2f}".format(difference)
                    # sigma = | value - ideal | / z
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                    # relative deviation in per cent = 100 * difference / ideal
                    percent = 100 * abs(value - ideal) / ideal
                    outData['percent'][i] = "{:.2f}".format(percent)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['sigma'][i] = '-'
                    outData['percent'][i] = '-'
            outDataZip = list(zip(outData['type'], outData['percent'], outData['sigma'], outData['atom1'], outData['atom2'],
                                  outData['value'], outData['ideal'], outData['difference'], outData['z']))
            outDataZip.sort(reverse=True)
            outData['type'], outData['percent'], outData['sigma'], outData['atom1'], outData['atom2'], \
                outData['value'], outData['ideal'], outData['difference'], outData['z'] = zip(*outDataZip)
    
            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Deviation<br>(in %)", data=outData['percent'])
            fullTable.addData(title="Bond<br>length (A)", data=outData['value'])
            fullTable.addData(title="Ideal<br>length (A)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (A)", data=outData['difference'])
            fullTable.addData(title="Sigma (A)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
            fullTable.addData(title="Note", data=outData['note'])

        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No bond length outliers observed.")
            
        if len(outAngle) > 0:
            n_outliers = len(outAngle)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Bond angle outliers:")
            outData = {'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'atom3': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers}
            for i, outlier in enumerate(outAngle):
                try: outData['atom1'][i] = str(outlier.findall('atom1')[0].text)
                except: outData['atom1'][i] = '-'
                try: outData['atom2'][i] = str(outlier.findall('atom2')[0].text)
                except: outData['atom2'][i] = '-'
                try: outData['atom3'][i] = str(outlier.findall('atom3')[0].text)
                except: outData['atom3'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['difference_float'], outData['sigma'], outData['atom1'], outData['atom2'], outData['atom3'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['sigma'], outData['atom1'], outData['atom2'], outData['atom3'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Atom 3", data=outData['atom3'])
            fullTable.addData(title="Bond<br>angle (°)", data=outData['value'])
            fullTable.addData(title="Ideal<br>angle (°)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (°)", data=outData['difference'])
            fullTable.addData(title="Sigma (°)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No bond angle outliers observed.")

        if len(outTorsion) > 0:
            n_outliers = len(outTorsion)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Torsion angle outliers:")
            outData = {'label': ["-"]*n_outliers,
                       'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'atom3': ["-"]*n_outliers,
                       'atom4': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'per': ["-"]*n_outliers,
                       'ideal_per': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers}
            for i, outlier in enumerate(outTorsion):
                try: outData['label'][i] = str(outlier.findall('label')[0].text)
                except: outData['label'][i] = '-'
                try: outData['atom1'][i] = str(outlier.findall('atom1')[0].text)
                except: outData['atom1'][i] = '-'
                try: outData['atom2'][i] = str(outlier.findall('atom2')[0].text)
                except: outData['atom2'][i] = '-'
                try: outData['atom3'][i] = str(outlier.findall('atom3')[0].text)
                except: outData['atom3'][i] = '-'
                try: outData['atom4'][i] = str(outlier.findall('atom4')[0].text)
                except: outData['atom4'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    periodicity = int(outlier.findall('per')[0].text)
                    outData['per'][i] = periodicity
                    # difference = | value - ideal |
                    differences = []
                    if periodicity > 0:
                        differences = []
                        for j in range(2 * periodicity):
                            ideal_with_per = ideal + j * 360 / periodicity
                            while ideal_with_per >= 360:
                                ideal_with_per = ideal_with_per - 360
                            differences.append(abs(ideal_with_per - value))
                        for j in range(2 * periodicity):
                            ideal_with_per = ideal - j * 360 / periodicity
                            while ideal_with_per <= -360:
                                ideal_with_per = ideal_with_per + 360
                            differences.append(abs(ideal_with_per - value))
                        periodicity_angle = int(360 / periodicity)
                        outData['ideal_per'][i] = "{:.2f}".format(ideal) + " &#177 n*" + str(periodicity_angle)
                    else:
                        differences.append(abs(ideal - value))
                        differences.append(abs(ideal + 360 - value))
                        differences.append(abs(ideal - 360 - value))
                        outData['ideal_per'][i] = "{:.2f}".format(ideal)
                    difference = min(differences)
                    sigma = abs(difference / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                    outData['difference'][i] = "{:.2f}".format(difference)
                    outData['difference_float'][i] = difference
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['per'][i] = '-'
                    outData['z'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'

            outDataZip = list(zip(outData['difference_float'], outData['label'], outData['atom1'], outData['atom2'], outData['atom3'],
                                  outData['atom4'], outData['value'], outData['ideal'], outData['per'], outData['ideal_per'], outData['z'], outData['difference'], outData['sigma']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['label'], outData['atom1'], outData['atom2'], outData['atom3'], \
                outData['atom4'], outData['value'], outData['ideal'], outData['per'], outData['ideal_per'], outData['z'], outData['difference'], outData['sigma'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Label", data=outData['label'])
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Atom 3", data=outData['atom3'])
            fullTable.addData(title="Atom 4", data=outData['atom4'])
            fullTable.addData(title="Torsion<br>angle (°)", data=outData['value'])
            fullTable.addData(title="Ideal angle (°)", data=outData['ideal_per'])
            fullTable.addData(title="Difference<br>from ideal (°)", data=outData['difference'])
            fullTable.addData(title="Sigma (°)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No torsion angle outliers observed.")

        if len(outChir) > 0:
            n_outliers = len(outChir)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Chiral volume outliers:")
            outData = {'atomc': ["-"]*n_outliers,
                       'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'atom3': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'both': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers,
                       # 'percent': ["-"]*n_outliers,
                       'signum': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers}
            for i, outlier in enumerate(outChir):
                try: outData['atomc'][i] = str(outlier.findall('atomc')[0].text)
                except: outData['atomc'][i] = '-'
                try: outData['atom1'][i] = str(outlier.findall('atom1')[0].text)
                except: outData['atom1'][i] = '-'
                try: outData['atom2'][i] = str(outlier.findall('atom2')[0].text)
                except: outData['atom2'][i] = '-'
                try: outData['atom3'][i] = str(outlier.findall('atom3')[0].text)
                except: outData['atom3'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    if bool(outlier.findall('both')[0].text) == "True":
                        both = True
                    else:
                        both = False
                    outData['both'][i] = str(both)
                    # difference = | value - ideal |
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    # sigma = | value - ideal | / z
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                    # percent = 100 * abs((value - ideal) / ideal) # be careful if ideal == 0
                    if int(sign(value) * sign(ideal)) != 1 and not both:
                        signum = "No"  # "Ideal value has an opposite sign"
                    else:
                        signum = 'Yes'
                    # outData['percent'][i] = "{:.2f}".format(percent)
                    outData['signum'][i] = signum
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['both'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
                    # outData['percent'][i] = '-'
                    outData['signum'][i] = '-'
            outDataZip = list(zip(outData['difference_float'], outData['atomc'], outData['atom1'], outData['atom2'], outData['atom3'],
                                  outData['value'], outData['ideal'], outData['z'], outData['both'], outData['difference'], outData['sigma'], outData['signum']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['atomc'], outData['atom1'], outData['atom2'], outData['atom3'], outData['value'], outData['ideal'], outData['z'], outData['both'], outData['difference'], outData['sigma'], outData['signum'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Chiral atom", data=outData['atomc'])
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Atom 3", data=outData['atom3'])
            # fullTable.addData(title="Deviation<br>(in %)", data=outData['percent'])
            fullTable.addData(title="Chiral<br>volume (A<sup>3</sup>)", data=outData['value'])
            fullTable.addData(title="Ideal<br>value (A<sup>3</sup>)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (A<sup>3</sup>)", data=outData['difference'])
            fullTable.addData(title="Z", data=outData['z'])
            fullTable.addData(title="Correct sign?", data=outData['signum'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No chirality outliers observed.")

        if len(outPlane) > 0:
            n_outliers = len(outPlane)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Planarity outliers:")
            outData = {'atom': ["-"]*n_outliers,
                       'label': ["-"]*n_outliers,
                       'dev': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers}
            for i, outlier in enumerate(outPlane):
                try: outData['label'][i] = str(outlier.findall('label')[0].text)
                except: outData['label'][i] = '-'
                try: outData['atom'][i] = str(outlier.findall('atom')[0].text)
                except: outData['atom'][i] = '-'
                try: outData['dev'][i] = "{:.2f}".format(float(outlier.findall('dev')[0].text))
                except: outData['dev'][i] = '-'
                try: outData['z'][i] = "{:.2f}".format(float(outlier.findall('z')[0].text))
                except: outData['z'][i] = '-'
            # Does not need to be sorted
            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Label", data=outData['label'])
            fullTable.addData(title="Atom", data=outData['atom'])
            fullTable.addData(title="Deviation (A)", data=outData['dev'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No planarity outliers observed.")

        if len(outVdw) > 0:
            n_outliers = len(outVdw)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Van der Waals outliers indicating close contacts between non-bonding atoms:")
            outData = {'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'type': ["-"]*n_outliers,
                       'note': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers}
            for i, outlier in enumerate(outVdw):
                try: outData['atom1'][i] = str(outlier.findall('atom1')[0].text)
                except: outData['atom1'][i] = '-'
                try: outData['atom2'][i] = str(outlier.findall('atom2')[0].text)
                except: outData['atom2'][i] = '-'
                try: outData['value'][i] = "{:.2f}".format(float(outlier.findall('value')[0].text))
                except: outData['value'][i] = '-'
                try: outData['ideal'][i] = "{:.2f}".format(float(outlier.findall('ideal')[0].text))
                except: outData['ideal'][i] = '-'
                try: outData['z'][i] = "{:.2f}".format(float(outlier.findall('z')[0].text))
                except: outData['z'][i] = '-'
                try:
                    outType = int(outlier.findall('type')[0].text)
                    if outType == 1:
                        outData['note'][i] = "Van der Waals"
                    elif outType == 2:
                        outData['note'][i] = "Torsion"
                    elif outType == 3:
                        outData['note'][i] = "Hydrogen bond"
                    elif outType == 4:
                        outData['note'][i] = "Metal"
                    elif outType == 5:
                        outData['note'][i] = "Dummy-nondummy"
                    elif outType == 6:
                        outData['note'][i] = "Dummy-nondummy"
                    elif outType > 6:
                        outData['note'][i] = "Symmetry related"
                    outData['type'][i] = -outType
                except:
                    outData['type'][i] = '-'
                    outData['note'][i] = '-'
                try:
                    # difference = | value - ideal |
                    difference = abs(float(outlier.findall('value')[0].text) - float(outlier.findall('ideal')[0].text))
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                except:
                    outData['difference'] = '-'
                    outData['difference_float'] = '-'

            outDataZip = list(zip(outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Distance (A)", data=outData['value'])
            fullTable.addData(title="Critical<br>distance (A)", data=outData['ideal'])
            fullTable.addData(title="Difference from<br>critical (A)", data=outData['difference'])
            fullTable.addData(title="Z", data=outData['z'])
            fullTable.addData(title="Type", data=outData['note'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No Van der Waals outliers observed.")

        if len(outStacd) > 0:
            n_outliers = len(outStacd)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Stacking distance outliers:")
            outData = {'plane1': ["-"]*n_outliers,
                       'plane2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers}
            for i, outlier in enumerate(outStacd):
                try: outData['plane1'][i] = str(outlier.findall('plane1')[0].text)
                except: outData['plane1'][i] = '-'
                try: outData['plane2'][i] = str(outlier.findall('plane2')[0].text)
                except: outData['plane2'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    # difference = | value - ideal |
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    # sigma = | value - ideal | / z
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Plane 1", data=outData['plane1'])
            fullTable.addData(title="Plane 2", data=outData['plane2'])
            fullTable.addData(title="Stacking<br>distance (A)", data=outData['value'])
            fullTable.addData(title="Ideal<br>distance (A)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (A)", data=outData['difference'])
            fullTable.addData(title="Sigma (A)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No stacking distance outliers.")

        if len(outStaca) > 0:
            n_outliers = len(outStaca)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Stacking angle outliers:")
            outData = {'plane1': ["-"]*n_outliers,
                       'plane2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'difference': ["-"]*n_outliers,
                       'difference_float': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'sigma': ["-"]*n_outliers}
            for i, outlier in enumerate(outStaca):
                try: outData['plane1'][i] = str(outlier.findall('plane1')[0].text)
                except: outData['plane1'][i] = '-'
                try: outData['plane2'][i] = str(outlier.findall('plane2')[0].text)
                except: outData['plane2'][i] = '-'
                try:
                    value = float(outlier.findall('value')[0].text)
                    outData['value'][i] = "{:.2f}".format(value)
                    ideal = float(outlier.findall('ideal')[0].text)
                    outData['ideal'][i] = "{:.2f}".format(ideal)
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    # difference = | value - ideal |
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    # sigma = | value - ideal | / z
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            fullTable = None
            fullTable = outlierFold.addTable()
            fullTable.addData(title="Plane 1", data=outData['plane1'])
            fullTable.addData(title="Plane 2", data=outData['plane2'])
            fullTable.addData(title="Stacking<br>angle (°)", data=outData['value'])
            fullTable.addData(title="Ideal<br>angle (°)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (°)", data=outData['difference'])
            fullTable.addData(title="Sigma (°)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No stacking angle outliers.")

    def addOutlierAnalysisOldToDelete(self, parent=None):
        if parent is None: parent=self
        
        outlierFold = parent.addFold(label="Outliers identified by Refmac",brief='Outliers')
        #Check to see if outliers are captured in the program XML
        outliersByCriteriaNodes = self.xmlnode.findall('.//OutliersByCriteria')
        if len(outliersByCriteriaNodes) == 0 or len(outliersByCriteriaNodes[0])==0:
            outlierFold.append('<span style="font-size:110%">No outliers observed </span>')
        else:
            outlierFold.append('<span style="font-size:110%">Residues failing one or more outlier test are flagged below with corresponding Z-score </span>')
            #identify a unique list of amino acids flagged by Refmac, and their associated set of deviations, flagged by the most deviant score observed
            naughtyBits = {}
            criteriaThatFail = set()
            for outlierDictNode in outliersByCriteriaNodes[-1
                                                           ]:
                outliers = outlierDictNode.findall('Outlier')
                criteriaThatFail.add(outlierDictNode.tag)
                for outlier in outliers:
                  try:
                    identifier1 = outlier.get('chainId1') + outlier.get('resId1') + outlier.get('ins1')
                    identifier2 = outlier.get('chainId2') + outlier.get('resId2') + outlier.get('ins2')
                    deviationAsSigma = (float(outlier.get('mod'))-float(outlier.get('ideal'))) / float(outlier.get('sigma'))
                    for identifier in [identifier1, identifier2]:
                        if identifier not in naughtyBits:
                            naughtyBits[identifier] = {}
                        if outlierDictNode.tag not in naughtyBits[identifier]:
                            naughtyBits[identifier][outlierDictNode.tag] = deviationAsSigma
                        elif abs(naughtyBits[identifier][outlierDictNode.tag]) < abs(deviationAsSigma):
                            naughtyBits[identifier][outlierDictNode.tag] = deviationAsSigma
                  except:
                    pass
            if len(criteriaThatFail) > 0:
                offendersTable = outlierFold.addTable(title='Ooh',label='Aah')
                identifiers = []
                for naughtyBit in naughtyBits:
                    identifiers.append(naughtyBit)
                offendersTable.addData(title='Residue',data=identifiers)
                for criterion in criteriaThatFail:
                    transgression = []
                    for naughtyBit in naughtyBits:
                        if criterion in naughtyBits[naughtyBit]:
                            transgression.append(str(naughtyBits[naughtyBit][criterion]))
                        else:
                            transgression.append('-')
                    offendersTable.addData(title=criterion,data=transgression)
                
        #Allow user to open a sub-folder in which criteria used in screening are supplied
        if len(outliersByCriteriaNodes) > 0 and len(outliersByCriteriaNodes[0]) > 0:
            criteriaFold = outlierFold.addFold(label='Criteria used to spot outliers',brief='Criteria')
            criteriaTable = criteriaFold.addTable(title='Ooh',label='Aah',select='.//OutliersByCriteria')
            
            types = []
            criteria = []
            for outlierTypeNode in outliersByCriteriaNodes[-1]:
                types.append(str(outlierTypeNode.tag))
                criteria.append(str(outlierTypeNode.findall('Criteria')[0].text))
            criteriaTable.addData(title='Interaction type',data=types)
            criteriaTable.addData(title='Criterion',data=criteria)

    def appendMolDisp(self, molDataNode=None, selectionText='all', carbonColour='yellow', othersByElement=True,style='CYLINDER',bondOrder=False):
        if molDataNode is None: return
        molDispNode = ET.SubElement(molDataNode,'MolDisp')
        if bondOrder:
            drawingStyleNode = ET.fromstring('<drawing_style><show_multiple_bonds>1</show_multiple_bonds><deloc_ring>1</deloc_ring></drawing_style>')
            molDispNode.append(drawingStyleNode)
        
        selectNode = ET.SubElement(molDispNode,'select')
        selectNode.text = selectionText
        '''
            selectionParametersNode = ET.SubElement(molDispNode,'selection_parameters')
            selectNode = ET.SubElement(selectionParametersNode,'select')
            selectNode.text = 'cid'
            
            cidNode = ET.SubElement(selectionParametersNode,'cid')
            cidNode.text = selectionText
            '''
        colourParametersNode = ET.SubElement(molDispNode,'colour_parameters')
        
        colourModeNode = ET.SubElement(colourParametersNode,'colour_mode')
        colourModeNode.text = 'one_colour'
        if othersByElement:
            nonCNode = ET.SubElement(colourParametersNode,'non_C_atomtype')
            nonCNode.text='1'
        oneColourNode = ET.SubElement(colourParametersNode,'one_colour')
        oneColourNode.text=carbonColour
        #
        styleNode = ET.SubElement(molDispNode,'style')
        styleNode.text=style

    def addMap(self, mapDataNode = None, fPhiObjectName='FPHIOUT', fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=1.0, isDifmap=False):
        if mapDataNode is None: return
        filenameNode = ET.SubElement(mapDataNode,'filename',database='filenames/'+fPhiObjectName)
        columnsNode = None
        if isDifmap:
            columnsNode = ET.SubElement(mapDataNode,'difColumns')
            isDifferenceMapNode = ET.SubElement(mapDataNode,'isDifferenceMap')
            isDifferenceMapNode.text='True'
        else:
            columnsNode = ET.SubElement(mapDataNode,'columns')
            isDifferenceMapNode = ET.SubElement(mapDataNode,'isDifferenceMap')
            isDifferenceMapNode.text='False'
        fColNode = ET.SubElement(columnsNode,'F')
        fColNode.text = fCol
        phiColNode = ET.SubElement(columnsNode,'PHI')
        phiColNode.text = phiCol
        
        modelNode = ET.SubElement(mapDataNode,'model')
        modelNode.text = model
        gridSizeNode = ET.SubElement(mapDataNode,'gridSize')
        gridSizeNode.text=str(gridSize)
        contourUnitsNode = ET.SubElement(mapDataNode,'contourUnits')
        contourUnitsNode.text = contourUnits
        mapDispNode = ET.SubElement(mapDataNode,'MapDisp')
        contourLevelNode = ET.SubElement(mapDispNode,'contourLevel')
        contourLevelNode.text = str(contourLevel)
        differenceNode = ET.SubElement(mapDispNode,'difference')
        colourNode = ET.SubElement(mapDispNode,'colour')
        if isDifmap:
            differenceNode.text='1'
            colourNode.text='green'
            colourNode2 = ET.SubElement(mapDispNode,'second_colour')
            colourNode2.text='red'
        else:
            differenceNode.text='0'
            colourNode.text='ice blue'

    def addRefinementPictures(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=False):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        
        #Provide default filenames to identify database object
        coordObjectName = objectNameMap.get('XYZ','XYZOUT')
        mapObjectName = objectNameMap.get('MAP','FPHIOUT')
        difmapObjectName = objectNameMap.get('DIFMAP','DIFFPHIOUT')
        dictObjectName = objectNameMap.get('DICT','DICTOUT')
        
        #I *do not know* why This is needed
        clearingDiv = parent.addDiv(style="clear:both;")
        
        pictureFold = parent.addFold(label='Picture', initiallyOpen=initiallyOpen,brief='Picture')
        pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='550px', tableWidth='260px', contentWidth='450px')
        clearingDiv = parent.addDiv(style="clear:both;")
        jobDirectory = jobInfo['fileroot']
        from core import CCP4Utils
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        import os
        baseScenePath = os.path.join(ccp4i2_root,'pipelines','prosmart_refmac','script','prosmart_refmac_1.scene.xml')
        monomerNodes = xmlnode.findall('.//ModelComposition[last()]/Monomer')
        
        #Subsets for snapshotting are each observed monomer and the whole molecule
        interestingBitsSet = set([monomerNode.get('id') for monomerNode in monomerNodes])
        for filterString in ['(SO4)','(SUL)','(GOL)','(EDO)','(CA)','(MG)']:
            interestingBitsSet -= set([monomerNode.get('id') for monomerNode in monomerNodes if filterString in monomerNode.get('id')])
        interestingBits=list(interestingBitsSet)+['all']
        
        for iMonomer, interestingBit in enumerate(interestingBits):
            baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath,useLXML=False) #This is lxml, not xml ...
            sceneNode = baseSceneXML#.findall('.//scene')[0]
            
            #Define data and associated display objects
            dataNode = ET.SubElement(sceneNode,'data')
            
            #Load model and draw representations
            molDataNode = ET.SubElement(dataNode,'MolData',id='id1')
            fileNode = ET.SubElement(molDataNode,'filename',database='filenames/'+coordObjectName)
            
            #Handle custom monomers
            try:
                dictPath = self.jobInfo['filenames'][dictObjectName]
                if dictPath is not None and dictPath != "":
                    from wrappers.acedrgNew.script.MyCIFDigestor import MyCIFFile
                    myCIFFile = MyCIFFile(filePath=dictPath)
                    #Find and act on list of residues in the dictionary provided
                    for compListBlock  in [block for block in myCIFFile.blocks if block.category == 'data_comp_list']:
                        for loop in compListBlock.loops():
                            for loopline in loop:
                                tlc  = loopline['three_letter_code']
                                customResNode = ET.fromstring('''<customResCIFFiles>
                                    <cifmonomer>
                                    <name>'''+tlc+'''</name>
                                    <filename>'''+dictPath+'''</filename>
                                    </cifmonomer>
                                    </customResCIFFiles>''')
                                molDataNode.append(customResNode)
            except BaseException as err:
                print('Failed analysing DICT'+"Base error: {0}".format(err))
            
            #Define view (possibly some advantage in doing this after loading the moelcule from which it is inferred)
            viewNode = ET.SubElement(sceneNode,'View')
            autoScaleNode = self.subElementWithText(viewNode,'scale_auto',1)
            autoScaleNode = self.subElementWithText(viewNode,'slab_enabled',1)
            centreMolDataNode = self.subElementWithText(viewNode,'centre_MolData','id1')
            centreSelectionNode = self.subElementWithText(viewNode,'centre_selection',interestingBit)
            
            orientationAutoNode = ET.SubElement(viewNode,'orientation_auto')
            orientationAutoMolDataNode = self.subElementWithText(orientationAutoNode,'molData','id1')
            orientationAutoSelectionNode = self.subElementWithText(orientationAutoNode,'selection', interestingBit)
            
            #Draw stuff around the interesting bit
            selectionText='neighb cid="'+interestingBit+'" maxd=10.0 group=all excl=central'
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='green', othersByElement=True,style='CYLINDERS')
            
            #Draw surface of surroundings
            surfDispObj = self.subElementWithText(molDataNode,"SurfaceDispobj","")
            self.subElementWithText(surfDispObj,"transparent",1)
            self.subElementWithText(surfDispObj,"visible",1)
            self.subElementWithText(surfDispObj,"opacity",0.7)
            selectionText='neighb cid="'+interestingBit+'" maxd=8.0 group=all excl=central,solvent'
            self.subElementWithText(surfDispObj,"select",selectionText)
            styleParams = self.subElementWithText(surfDispObj,"style_parameters","")
            self.subElementWithText(styleParams,"select","cid")
            cidNode = self.subElementWithText(styleParams,"cid","not {"+interestingBit+" or solvent}")
            
            #Draw the interesting bit
            selectionText=interestingBit
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='yellow', othersByElement=True,style='BALLSTICK',bondOrder=True)
            
            #Load and draw 2Fo-Fc
            mapDataNode = ET.SubElement(dataNode,'MapData',id='id3')
            self.addMap(mapDataNode, fPhiObjectName=mapObjectName, fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=1.0, isDifmap=False)
            
            #Load and draw Fo-Fc
            mapDataNode = ET.SubElement(dataNode,'MapData',id='id4')
            self.addMap(mapDataNode, fPhiObjectName=difmapObjectName, fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=3.0, isDifmap=True)
            
            # Dump out the XML
            et = ET.ElementTree(baseSceneXML)
            sceneFilePath = os.path.join(jobDirectory,'monomer'+str(iMonomer)+'.scene.xml')
            et.write(sceneFilePath)
            # And add the picture
            pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of selection "'+interestingBit+'"')

        #And finally the full monty picture :-)
        baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath,useLXML=False)
        sceneNode = baseSceneXML#.findall('scene')[0]
        
        #Define data and associated display objects
        dataNode = ET.SubElement(sceneNode,'data')
        
        #Load model and draw representations
        molDataNode = ET.SubElement(dataNode,'MolData',id='id1')
        fileNode = ET.SubElement(molDataNode,'filename',database='filenames/'+coordObjectName)
        
        #Handle custom monomers
        try:
            dictPath = self.jobInfo['filenames'][dictObjectName]
            if dictPath is not None and dictPath != "":
                from wrappers.acedrgNew.script.MyCIFDigestor import MyCIFFile
                myCIFFile = MyCIFFile(filePath=dictPath)
                #Find and act on list of residues in the dictionary provided
                for compListBlock  in [block for block in myCIFFile.blocks if block.category == 'data_comp_list']:
                    for loop in compListBlock.loops():
                        for loopline in loop:
                            tlc  = loopline['three_letter_code']
                            customResNode = ET.fromstring('''<customResCIFFiles>
                                <cifmonomer>
                                <name>'''+tlc+'''</name>
                                <filename>'''+dictPath+'''</filename>
                                </cifmonomer>
                                </customResCIFFiles>''')
                            molDataNode.append(customResNode)
        except BaseException as err:
            print('Failed analysing DICT'+"Base error: {0}".format(err))
        
        #Define view (possibly some advantage in doing this after loading the moelcule from which it is inferred)
        viewNode = ET.SubElement(sceneNode,'View')
        autoScaleNode = self.subElementWithText(viewNode,'scale_auto',1)
        autoScaleNode = self.subElementWithText(viewNode,'slab_enabled',1)
        centreMolDataNode = self.subElementWithText(viewNode,'centre_MolData','id1')
        centreSelectionNode = self.subElementWithText(viewNode,'centre_selection',interestingBit)
        
        orientationAutoNode = ET.SubElement(viewNode,'orientation_auto')
        orientationAutoMolDataNode = self.subElementWithText(orientationAutoNode,'molData','id1')
        orientationAutoSelectionNode = self.subElementWithText(orientationAutoNode,'selection', interestingBit)
        
        #Draw all atoms as sticks
        selectionText='all'
        self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='green', othersByElement=True,style='FATBONDS')
        
        for monomer in monomerNodes:
            interestingBit = monomer.get("id")
            #Draw stuff around the interesting bit
            selectionText='neighb cid="'+interestingBit+'" maxd=10.0 group=all excl=central'
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='green', othersByElement=True,style='CYLINDERS')
            
            #Draw surface of surroundings
            surfDispObj = self.subElementWithText(molDataNode,"SurfaceDispobj","")
            self.subElementWithText(surfDispObj,"transparent",1)
            self.subElementWithText(surfDispObj,"visible",1)
            self.subElementWithText(surfDispObj,"opacity",0.7)
            selectionText='neighb cid="'+interestingBit+'" maxd=8.0 group=all excl=central,solvent'
            self.subElementWithText(surfDispObj,"select",selectionText)
            styleParams = self.subElementWithText(surfDispObj,"style_parameters","")
            self.subElementWithText(styleParams,"select","cid")
            cidNode = self.subElementWithText(styleParams,"cid","not {"+interestingBit+" or solvent}")
            
            #Draw the interesting bit
            selectionText=interestingBit
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='yellow', othersByElement=True,style='BALLSTICK',bondOrder=True)
        
        #Load and draw 2Fo-Fc
        mapDataNode = ET.SubElement(dataNode,'MapData',id='id3')
        self.addMap(mapDataNode, fPhiObjectName=mapObjectName, fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=1.0, isDifmap=False)
        
        #Load and draw Fo-Fc
        mapDataNode = ET.SubElement(dataNode,'MapData',id='id4')
        self.addMap(mapDataNode, fPhiObjectName=difmapObjectName, fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=3.0, isDifmap=True)
        
        # Dump out the XML
        et = ET.ElementTree(baseSceneXML)
        sceneFilePath = os.path.join(jobDirectory,'fullMonty.scene.xml')
        et.write(sceneFilePath)
        # And add the picture
        pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of Full Monty')


    def subElementWithText(self, parent, tag, text):
        result = ET.SubElement(parent, tag)
        result.text = str(text)
        return result

    def addScrollableDownloadableTable1(self, cycle_data, xmlnode=None, parent=None, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        #create a "shell" div to contain the scrollable table and the hyperlink
        scrollableDownloadableTableDiv = parent.addDiv(style="height:250px; width:315px;float:left;margin-top:2px;")
        #place a scrollable div into the shell: the table will be inserted into this div
        scrollableTableDiv = scrollableDownloadableTableDiv.addDiv(style="height:225px; width:300px;clear:both;overflow:auto;")
        #Put table1 into this (autoscrolling) div
        table1 = self.addTable1(cycle_data, parent=scrollableTableDiv,internalId=internalId)
        #scrollableDownloadableTableDiv.addDiv(style="height:10px;width:15px; float:right;")
        #Add a hyperlink to this table in the "tables_as_csv_files" directory
        #print '\n\nJobInfo',self.jobInfo
        download = scrollableDownloadableTableDiv.addDownload(jobInfo=self.jobInfo,dataName=table1.id)

    def addTable1(self, cycle_data, xmlnode=None, parent=None, downloadable=False, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        table1 = parent.addTable(xmlnode=xmlnode, style="width:240px;float:left;", downloadable=downloadable, outputXml=self.outputXml, internalId=internalId)
        
        statisticNames = []
        statisticInitial = []
        statisticFinal = []

        if isnumber(cycle_data["R2work"][-1]):  # refinement against intensities
            statisticNames.append('R2work')
            statisticInitial.append(cycle_data['R2work'][0])
            statisticFinal.append(cycle_data['R2work'][-1])
            if isnumber(cycle_data["R2free"][-1]):
                statisticNames.append('R2free')
                statisticInitial.append(cycle_data['R2free'][0])
                statisticFinal.append(cycle_data['R2free'][-1])
            statisticNames.append('CCIwork_avg')
            statisticInitial.append(cycle_data['CCIworkavg'][0])
            statisticFinal.append(cycle_data['CCIworkavg'][-1])
            if isnumber(cycle_data["CCIfreeavg"][-1]):
                statisticNames.append('CCIfree_avg')
                statisticInitial.append(cycle_data['CCIfreeavg'][0])
                statisticFinal.append(cycle_data['CCIfreeavg'][-1])
        else:                               # refinement against amplitudes
            statisticNames.append('Rwork')
            statisticInitial.append(cycle_data['Rwork'][0])
            statisticFinal.append(cycle_data['Rwork'][-1])
            if isnumber(cycle_data["Rfree"][-1]):
                statisticNames.append('Rfree')
                statisticInitial.append(cycle_data['Rfree'][0])
                statisticFinal.append(cycle_data['Rfree'][-1])
            statisticNames.append('CCFwork_avg')
            statisticInitial.append(cycle_data['CCFworkavg'][0])
            statisticFinal.append(cycle_data['CCFworkavg'][-1])
            if isnumber(cycle_data["CCFfreeavg"][-1]):
                statisticNames.append('CCFfree_avg')
                statisticInitial.append(cycle_data['CCFfreeavg'][0])
                statisticFinal.append(cycle_data['CCFfreeavg'][-1])
        if isnumber(cycle_data['rmsBOND'][-1]):
            statisticNames.append('<i>RMS Deviations</i>')
            statisticInitial.append(' ')
            statisticFinal.append(' ')
            statisticNames.append('&nbsp;&nbsp;&nbsp;bond')
            statisticInitial.append(cycle_data['rmsBOND'][0])
            statisticFinal.append(cycle_data['rmsBOND'][-1])
            statisticNames.append('&nbsp;&nbsp;&nbsp;angle')
            statisticInitial.append(cycle_data['rmsANGLE'][0])
            statisticFinal.append(cycle_data['rmsANGLE'][-1])
            statisticNames.append('&nbsp;&nbsp;&nbsp;chiral')
            statisticInitial.append(cycle_data['rmsCHIRAL'][0])
            statisticFinal.append(cycle_data['rmsCHIRAL'][-1])

            statisticNames.append('<i>RMS Z-scores</i>')
            statisticInitial.append(' ')
            statisticFinal.append(' ')
            statisticNames.append('&nbsp;&nbsp;&nbsp;bond')
            statisticInitial.append(cycle_data['zBOND'][0])
            statisticFinal.append(cycle_data['zBOND'][-1])
            statisticNames.append('&nbsp;&nbsp;&nbsp;angle')
            statisticInitial.append(cycle_data['zANGLE'][0])
            statisticFinal.append(cycle_data['zANGLE'][-1])
            statisticNames.append('&nbsp;&nbsp;&nbsp;chiral')
            statisticInitial.append(cycle_data['zCHIRAL'][0])
            statisticFinal.append(cycle_data['zCHIRAL'][-1])
        table1.addData(title='', data=statisticNames)
        table1.addData(title='Initial', data=statisticInitial)
        table1.addData(title='Final', data=statisticFinal)
        return table1

    def addTable1_ToDelete(self, xmlnode=None, parent=None, downloadable=False, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self

        # TO DO width and height
        table1 = parent.addTable(xmlnode=xmlnode, style="width:240px;float:left;", downloadable=downloadable,outputXml=self.outputXml,internalId=internalId)
        
        #ResolutionLowNode =xmlnode.findall('.//Overall_stats/resolution_low')
        #ResolutionHighNode =xmlnode.findall('.//Overall_stats/resolution_high')
        #ReflectionsAll =xmlnode.findall('.//Overall_stats/n_reflections_all')
        #ResolutionsFree =xmlnode.findall('.//Overall_stats/n_reflections_free')
        #ReflectionsWorkNode =xmlnode.findall('.//Overall_stats/resolution_high')

        R2WorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/R2work')
        R2FreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/R2free')
        CCIWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCIworkavg')
        CCIFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCIfreeavg')

        RWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/Rwork')
        RFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/Rfree')
        CCFWorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCFworkavg')
        CCFFreeNodes = xmlnode.findall('.//cycle[last()]/data/summary/CCFfreeavg')

        RMSBondsNodes = xmlnode.findall('.//cycle[last()]/geom/summary/rmsd/Bond_distances_non_H')
        RMSAnglesNodes = xmlnode.findall('.//cycle[last()]/geom/summary/rmsd/Bond_angles_non_H')
        # RMS Chiral_centre ?
        
        #MeanBChainNameNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/chain_name')
        #MeanBAllCountNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/number')
        #MeanBAllAverageNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/average')
        
        statisticNames = []
        statisticValues = []
        
        #statisticNames.append('Resolution')
        #if len(ResolutionLowNode)>0: statisticValues.append("{0:.2f}".format(float(ResolutionLowNode[0].text))+'-'+"{0:.2f}".format(float(ResolutionHighNode[0].text)))
        #else: statisticValues.append("-")
        
        #statisticNames.append('No. reflections all/free')
        #if len(ReflectionsAll)>0:
        #    nRefText = ReflectionsAll[0].text+'/'
        #    if len(ResolutionsFree)>0: nRefText += ResolutionsFree[0].text
        #    else: nRefText += "-"
        #    statisticValues.append(nRefText)
        #else: statisticValues.append("-")

        if len(R2WorkNodes) > 0:  # Refinement against intensities
            statisticNames.append('R2work/R2free')
            if len(R2WorkNodes) > 0:
                text = "{:.4f}".format(float(R2WorkNodes[0].text)) + ' / '
                if len(R2FreeNodes) > 0: text += "{:.4f}".format(float(R2FreeNodes[0].text))
                else: text += "-"
                statisticValues.append(text)
            else: statisticValues.append("-")

            statisticNames.append('CCIwork_avg/CCIfree_avg')
            if len(CCIWorkNodes) > 0:
                text = "{:.4f}".format(float(CCIWorkNodes[0].text)) + ' / '
                if len(CCIFreeNodes) > 0: text += "{:.4f}".format(float(CCIFreeNodes[0].text))
                else: text += "-"
                statisticValues.append(text)
            else: statisticValues.append("-")

        else:  # Refinement against amplitudes
            statisticNames.append('Rwork/Rfree')
            if len(RWorkNodes) > 0:
                text = "{:.4f}".format(float(RWorkNodes[0].text)) + ' / '
                if len(RFreeNodes) > 0: text += "{:.4f}".format(float(RFreeNodes[0].text))
                else: text += "-"
                statisticValues.append(text)
            else: statisticValues.append("-")

            statisticNames.append('CCFwork_avg/CCFfree_avg')
            if len(CCFWorkNodes) > 0:
                text = "{:.4f}".format(float(CCFWorkNodes[0].text)) + ' / '
                if len(CCFFreeNodes) > 0: text += "{:.4f}".format(float(CCFFreeNodes[0].text))
                else: text += "-"
                statisticValues.append(text)
            else: statisticValues.append("-")

        if len(self.xmlnode.findall("RigidMode")) == 0:
            statisticNames.append('<i>RMS Deviations</i>')
            statisticValues.append(' ')

            statisticNames.append('Bonds')
            if len(RMSBondsNodes)>0:
                statisticValues.append("{:.4f}".format(float(RMSBondsNodes[0].text)))
            else: statisticValues.append("-")

            statisticNames.append('Angles')
            if len(RMSAnglesNodes)>0:
                statisticValues.append("{:.4f}".format(float(RMSAnglesNodes[0].text)))
            else: statisticValues.append("-")

            #statisticNames.append('<i>Chain B-factors</i>')
            #statisticValues.append('<i>mean B (#atoms)</i>')
            #statisticNames += [MeanBChainNameNodes[i].text for i in range(len(MeanBChainNameNodes))]
            #statisticValues +=  [("{0:.1f}".format(float(MeanBAllAverageNodes[i].text))+'('+MeanBAllCountNodes[i].text+')') for i in range(len(MeanBChainNameNodes))]

        table1.addData(title='', data=statisticNames)
        table1.addData(title='Initial', data=statisticValues)
        
        return table1

def test(xmlFile=None,jobId=None,reportFile=None):
    import sys,os
    print(xmlFile)
    try:
        text = open( xmlFile ).read()
        xmlnode = ET.fromstring( text, PARSER() )
    except:
        print('FAILED loading XML file:', kw['xmlFile'])
    if reportFile is None and xmlFile is not None:
        reportFile = os.path.join(os.path.split(xmlFile)[0],'report.html')
    r = refmac_report(xmlFile=xmlFile,jobId=jobId, xmlnode=xmlnode)
    r.as_html_file(reportFile)

if __name__ == "__main__":
    import sys
    servalcat_xtal_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
