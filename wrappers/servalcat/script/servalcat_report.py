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


class servalcat_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'servalcat'
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
        self.addTablePerCycle(cycle_data, parent=summaryFold, initialFinalOnly=True)
        self.addGraphPerCycle(parent=summaryFold)
        if withTables:
            clearingDiv = parent.addDiv(style="clear:both;")
            perCycleFold = parent.addFold(label='Per cycle statistics', brief='Per cycle', initiallyOpen=False)
            self.addTablePerCycle(cycle_data, parent=perCycleFold, initialFinalOnly=False)
            self.addGraphsVsResolution()
            self.addOutlierAnalysis()

    def addGraphPerCycle(self, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        # Note that when I add the progressgraph, I have to ensure that the select is rooted in my own xmlnode
        progressGraph = parent.addFlotGraph(
            title="Refinement results",
            xmlnode=self.xmlnode,
            style="height:250px;width:400px;float:left;")
        progressGraph.addData(title="Cycle", select=".//cycle/Ncyc") # ycol=1
        progressGraph.addData(title="-LL", select=".//cycle/data/summary/minusLL") # ycol=2
        spa_refinement = False
        intensity_based_refinement = False
        if len(xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')) > 0:
            progressGraph.addData(title="⟨FSCmodel⟩", select=".//cycle/data/summary/FSCaverage", expr="x if float(x)>=0.0 else ''")  # ycol=3
            spa_refinement = True
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/Rwork')) > 0:
            progressGraph.addData(title="Rwork", select=".//cycle/data/summary/Rwork", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="⟨CCFwork⟩", select=".//cycle/data/summary/CCFworkavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
            if len(xmlnode.findall('.//cycle[last()]/data/summary/Rfree')) > 0:
                progressGraph.addData(title="Rfree", select=".//cycle/data/summary/Rfree", expr="x if float(x)>=0.0 else '-'")  # ycol=5
                progressGraph.addData(title="⟨CCFfree⟩", select=".//cycle/data/summary/CCFfreeavg", expr="x if float(x)>=-1.0 else '-'")  # ycol=6
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/R1work')) > 0:
            intensity_based_refinement = True
            progressGraph.addData(title="R1work", select=".//cycle/data/summary/R1work", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="⟨CCIwork⟩", select=".//cycle/data/summary/CCIworkavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
            if len(xmlnode.findall('.//cycle[last()]/data/summary/R1free')) > 0:
                progressGraph.addData(title="R1free", select=".//cycle/data/summary/R1free", expr="x if float(x)>=0.0 else '-'")  # ycol=5
                progressGraph.addData(title="⟨CCIfree⟩", select=".//cycle/data/summary/CCIfreeavg", expr="x if float(x)>=-1.0 else '-'") # ycol=6
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/R')) > 0:
            progressGraph.addData(title="R", select=".//cycle/data/summary/R", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="⟨CCF⟩", select=".//cycle/data/summary/CCFavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/R1')) > 0:
            intensity_based_refinement = True
            progressGraph.addData(title="R1", select=".//cycle/data/summary/R1", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title="⟨CCI⟩", select=".//cycle/data/summary/CCIavg", expr="x if float(x)>=-1.0 else ''")  # ycol=4

        if spa_refinement:
            plotCC = progressGraph.addPlotObject()
            plotCC.append('title', '⟨FSCmodel⟩')
            plotCC.append('plottype', 'xy')
            plotCC.append('xlabel', 'Cycle')
            plotCC.append('ylabel', '⟨FSCmodel⟩')
            plotCC.append('yrange', min=0.0, max=1.0)
            plotCC.append('xintegral', 'true')
            plotCC.append('legendposition', x=0, y=0)
            plotLine = plotCC.append('plotline', xcol=1, ycol=3)
            plotLine.append('colour', 'orange')
            plotLine.append('symbolsize', '0')
        elif intensity_based_refinement:
            addCorrelationProgress(progressGraph)
            addRValuesProgress(progressGraph)
        else:
            addRValuesProgress(progressGraph)
            addCorrelationProgress(progressGraph)

        plotLL = progressGraph.addPlotObject()
        plotLL.append('title', '-LL')
        plotLL.append('plottype', 'xy')
        plotLL.append('xlabel', 'Cycle')
        plotLL.append('ylabel', '-LL')
        plotLL.append('xintegral', 'true')
        plotLL.append('legendposition', x=0, y=0)
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
            progressGraph2.addData(title="Cycle", select=".//cycle/Ncyc")
            progressGraph2.addData(title="RMSD_bond", select=".//cycle/geom/summary/rmsd/Bond_distances_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSD_angle", select=".//cycle/geom/summary/rmsd/Bond_angles_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSZ_bond", select=".//cycle/geom/summary/rmsZ/Bond_distances_non_H", expr="x if float(x)>=0.0 else '-'")
            progressGraph2.addData(title="RMSZ_angle", select=".//cycle/geom/summary/rmsZ/Bond_angles_non_H", expr="x if float(x)>=0.0 else '-'")
            cycles_list_without_zero = range(1, len(xmlnode.findall('.//cycle')))  # weight for the 0th cycle is not defined
            progressGraph2.addData(title="Cycle", data=cycles_list_without_zero)
            progressGraph2.addData(title="Weight", select=".//cycle/weight", expr="x if float(x)>=0.0 else '-'")
            plotRmsd = progressGraph2.addPlotObject()
            plotRmsd.append('title', 'RMS Deviations')
            plotRmsd.append('plottype', 'xy')
            plotRmsd.append('yrange', rightaxis='false')
            plotRmsd.append('xlabel', 'Cycle')
            plotRmsd.append('ylabel', ' ')
            plotRmsd.append('yrange', min=0.0)
            plotRmsd.append('xintegral', 'true')
            plotRmsd.append('legendposition', x=1, y=1)
            plotLine = plotRmsd.append('plotline', xcol=1, ycol=2, rightaxis='false')
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
            plotRmsz.append('yrange', min=0.0)
            plotRmsz.append('xintegral', 'true')
            plotRmsz.append('legendposition', x=1, y=1)
            plotLine = plotRmsz.append('plotline', xcol=1, ycol=4, rightaxis='false')
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')
            plotLine = plotRmsz.append('plotline', xcol=1, ycol=5, rightaxis='false')
            plotLine.append('colour', 'red')
            plotLine.append('symbolsize', '0')

            plotWeight = progressGraph2.addPlotObject()
            plotWeight.append('title', 'Weight')
            plotWeight.append('plottype', 'xy')
            plotWeight.append('xlabel', '')
            plotWeight.append('ylabel', '')
            plotWeight.append('yrange', min=0.0)
            plotWeight.append('xintegral', 'true')
            plotWeight.append('legendposition', x=0, y=0)
            plotLine = plotWeight.append('plotline', xcol=6, ycol=7, rightaxis='false')
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')

        clearingDiv = parent.addDiv(style="clear:both;")
        if len(xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')) > 0:
            return
        noteDiv = parent.addDiv(style="margin-bottom:0;")
        note = ""
        if len(xmlnode.findall('.//cycle[last()]/data/summary/Rwork')) > 0 or \
                len(xmlnode.findall('.//cycle[last()]/data/summary/R')) > 0:
            note = "⟨CC<sub><i>F</i></sub>⟩ is the correlation coefficient between calculated and observed amplitudes averaged over resolution shells."
            # note = "R = &#931; | <i>F</i><sub>obs</sub> &#8722; <i>F</i><sub>calc</sub> | / &#931; <i>F</i><sub>obs</sub><br />" + \
            #     "CC<sub><i>F</i></sub> = (⟨<i>F</i><sub>obs</sub> <i>F</i><sub>calc</sub>⟩ &#8722; ⟨<i>F</i><sub>obs</sub>⟩⟨<i>F</i><sub>calc</sub>⟩ / " + \
            #     "&#8730; (⟨<i>F</i><sub>obs</sub><sup>2</sup>⟩ &#8722; ⟨<i>F</i><sub>obs</sub>⟩<sup>2</sup>) &#8730; (⟨<i>F</i><sub>calc</sub><sup>2</sup>⟩ &#8722; ⟨<i>F</i><sub>calc</sub>⟩<sup>2</sup>)<br />" + \
            #     "⟨CC<sub><i>F</i></sub>⟩ = &#931; <i>N<sub>i</sub></i> CC<sub><i>F,i</i></sub> / &#931; <i>N<sub>i</sub></i> &#160;&#160; where <i>N<sub>i</sub></i> and CC<sub><i>F,i</i></sub> are the number of reflections and the correlation in a resolution bin."
        elif len(xmlnode.findall('.//cycle[last()]/data/summary/R1work')) > 0 or \
                len(xmlnode.findall('.//cycle[last()]/data/summary/R1')) > 0:
            note = "R1 is the R-value calculated from the square root of intensities where <i>I</i>/&#963;(<i>I</i>) > 2.<br />"
            note += "⟨CC<sub><i>I</i></sub>⟩ is the correlation coefficient between calculated and observed intensities averaged over resolution shells."
            # note = "R1 = &#931; | &#8730;<i>I</i><sub>obs</sub> &#8722; <i>F</i><sub>calc</sub> | / &#931; &#8730;<i>I</i><sub>obs</sub> &#160;&#160; " + \
            #     "where <i>I</i><sub>obs</sub>/&#963;(<i>I</i><sub>obs</sub>) >= 2<br />" + \
            #     "CC<sub><i>I</i></sub> = (⟨<i>I</i><sub>obs</sub> <i>I</i><sub>calc</sub>⟩ &#8722; ⟨<i>I</i><sub>obs</sub>⟩⟨<i>I</i><sub>calc</sub>⟩ / " + \
            #     "&#8730; (⟨<i>I</i><sub>obs</sub><sup>2</sup>⟩ &#8722; ⟨<i>I</i><sub>obs</sub>⟩<sup>2</sup>) &#8730; (⟨<i>I</i><sub>calc</sub><sup>2</sup>⟩ &#8722; ⟨<i>I</i><sub>calc</sub>⟩<sup>2</sup>)<br />" + \
            #     "⟨CC<sub><i>I</i></sub>⟩ = &#931; <i>N<sub>i</sub></i> CC<sub><i>I,i</i></sub> / &#931; <i>N<sub>i</sub></i> &#160;&#160; where <i>N<sub>i</sub></i> and CC<sub><i>I,i</i></sub> are the number of reflections and the correlation in a resolution bin."
        noteDiv.append(note)

    def getCycleData(self, xmlnode=None):
        if xmlnode is None: xmlnode = self.xmlnode
        FSCaverageNodes = xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')
        R1WorkNodes = xmlnode.findall('.//cycle[last()]/data/summary/R1work')
        R1Nodes = xmlnode.findall('.//cycle[last()]/data/summary/R1')
        RNodes = xmlnode.findall('.//cycle[last()]/data/summary/R')
        all_cycles = xmlnode.findall('.//cycle')
        ncyc = len(all_cycles)
        cycle_data = {'mode':['-']*ncyc,
                      'cycle':['-']*ncyc,
                      '-LL':['-']*ncyc,
                      'FSCaverage':['-']*ncyc,
                      'Rwork':['-']*ncyc,
                      'Rfree':['-']*ncyc,
                      'R':['-']*ncyc,
                      'CCFworkavg':['-']*ncyc,
                      'CCFfreeavg':['-']*ncyc,
                      'CCFavg':['-']*ncyc,
                      'R1work':['-']*ncyc,
                      'R1free':['-']*ncyc,
                      'R1':['-']*ncyc,
                      'CCIworkavg':['-']*ncyc,
                      'CCIfreeavg':['-']*ncyc,
                      'CCIavg':['-']*ncyc,
                      'rmsBOND':['-']*ncyc,
                      'rmsANGLE':['-']*ncyc,
                      'rmsCHIRAL':['-']*ncyc,
                      'zBOND':['-']*ncyc,
                      'zANGLE':['-']*ncyc,
                      'zCHIRAL':['-']*ncyc,
                      'weight':['-']*ncyc,}
        idx = 0
        for cycle in all_cycles:
            cycle_data['mode'][idx] = 'Restr'
            try: cycle_data['cycle'][idx] = str(int(cycle.findall('Ncyc')[0].text))
            except: pass
            try: cycle_data['-LL'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/-LL')[0].text))
            except: pass
            if idx == 0:  # weight for the 0th cycle is not defined
                cycle_data['weight'][idx] = '-'
            else:
                try: cycle_data['weight'][idx] = "{:.2f}".format(float(cycle.findall('weight')[0].text))
                except: pass
            if len(FSCaverageNodes) > 0: # SPA refinement
                try: cycle_data['FSCaverage'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/FSCaverage')[0].text))
                except: pass
            elif len(R1Nodes) > 0:       # Refinement against intensities without free flags
                try: cycle_data['R1'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R1')[0].text))
                except: pass
                try: cycle_data['CCIavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCIavg')[0].text))
                except: pass
            elif len(RNodes) > 0:       # Refinement against amplitudes without free flags
                try: cycle_data['R'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R')[0].text))
                except: pass
                try: cycle_data['CCFavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCFavg')[0].text))
                except: pass
            elif len(R1WorkNodes) > 0:  # Refinement against intensities with free flags
                try: cycle_data['R1work'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R1work')[0].text))
                except: pass
                try: cycle_data['R1free'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/R1free')[0].text))
                except: pass
                try: cycle_data['CCIworkavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCIworkavg')[0].text))
                except: pass
                try: cycle_data['CCIfreeavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCIfreeavg')[0].text))
                except: pass
            else:                     # Refinement against amplitudes with free flags
                try: cycle_data['Rwork'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/Rwork')[0].text))
                except: pass
                try: cycle_data['Rfree'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/Rfree')[0].text))
                except: pass
                try: cycle_data['CCFworkavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCFworkavg')[0].text))
                except: pass
                try: cycle_data['CCFfreeavg'][idx] = "{:.4f}".format(float(cycle.findall('data/summary/CCFfreeavg')[0].text))
                except: pass
            if cycle_data['mode'][idx] == 'Restr':
                try: cycle_data['rmsBOND'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsd/Bond_distances_non_H')[0].text))
                except: pass
                try: cycle_data['rmsANGLE'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsd/Bond_angles_non_H')[0].text))
                except: pass
                try: cycle_data['rmsCHIRAL'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsd/Chiral_centres')[0].text))
                except: pass
                try: cycle_data['zBOND'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsZ/Bond_distances_non_H')[0].text))
                except: pass
                try: cycle_data['zANGLE'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsZ/Bond_angles_non_H')[0].text))
                except: pass
                try: cycle_data['zCHIRAL'][idx] = "{:.3f}".format(float(cycle.findall('geom/summary/rmsZ/Chiral_centres')[0].text))
                except: pass
            idx += 1
        return cycle_data

    def addTablePerCycle(self, cycle_data, parent=None, initialFinalOnly=False):
        if parent is None: parent = self
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

        clearingDiv = parent.addDiv(style="clear:both;")
        fullTable = None
        fullTable = parent.addTable()
        if len(TLSIdx) > 0 and len(RestrIdx) > 0:
           mode = ['TLS']*len(TLSIdx) + ['Full Atom']*len(RestrIdx)
        fullTable.addData(title="Cycle", data=cycle_data_sel['cycle'])
        if isnumber(cycle_data_sel['FSCaverage'][-1]):  # SPA refinement
            fullTable.addData(title="⟨FSCmodel⟩", data=cycle_data_sel['FSCaverage'])
        elif isnumber(cycle_data_sel['R1'][-1]):        # Refinement against intensities without free flags
            fullTable.addData(title="R1", data=cycle_data_sel['R1'])
            fullTable.addData(title="⟨CCI⟩", data=cycle_data_sel['CCIavg'])
        elif isnumber(cycle_data_sel['R'][-1]):         # Refinement against amplitudes without free flags
            fullTable.addData(title="R", data=cycle_data_sel['R'])
            fullTable.addData(title="⟨CCF⟩", data=cycle_data_sel['CCFavg'])
        elif isnumber(cycle_data_sel['R1work'][-1]):    # Refinement against intensities with free flags
            fullTable.addData(title="R1work", data=cycle_data_sel['R1work'])
            if isnumber(cycle_data_sel['R1free'][-1]):
                fullTable.addData(title="R1free", data=cycle_data_sel['R1free'])
            fullTable.addData(title="⟨CCIwork⟩", data=cycle_data_sel['CCIworkavg'])
            if isnumber(cycle_data_sel['CCIfreeavg'][-1]):
                fullTable.addData(title="⟨CCIfree⟩", data=cycle_data_sel['CCIfreeavg'])
        else:                                # Refinement against amplitudes with free flags
            fullTable.addData(title="Rwork", data=cycle_data_sel['Rwork'])
            if isnumber(cycle_data_sel['Rfree'][-1]):
                fullTable.addData(title="Rfree", data=cycle_data_sel['Rfree'])
            fullTable.addData(title="⟨CCFwork⟩", data=cycle_data_sel['CCFworkavg'])
            if isnumber(cycle_data_sel['CCFfreeavg'][-1]):
                fullTable.addData(title="⟨CCFfree⟩", data=cycle_data_sel['CCFfreeavg'])
        if isnumber(cycle_data_sel['rmsANGLE'][-1]):
            fullTable.addData(title="RMSD (bond/angle/chiral)", subtitle="Bond", data=cycle_data_sel['rmsBOND'])
            fullTable.addData(subtitle="Angle", data=cycle_data_sel['rmsANGLE'])
            fullTable.addData(subtitle="Chiral", data=cycle_data_sel['rmsCHIRAL'])
            fullTable.addData(title="RMSZ (bond/angle/chiral)", subtitle="Bond", data=cycle_data_sel['zBOND'])
            fullTable.addData(subtitle="Angle", data=cycle_data_sel['zANGLE'])
            fullTable.addData(subtitle="Chiral", data=cycle_data_sel['zCHIRAL'])
        fullTable.addData(title="Weight", data=cycle_data_sel['weight'])


    def addGraphsVsResolution(self, parent=None, xmlnode=None, internalIdPrefix=''):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        reportFold = parent.addFold(label='Statistics vs. resolution', brief='Other')

        gallery = reportFold.addObjectGallery(
            height='310px', contentWidth='420px', tableWidth='360px', style='float:left;width:800px;')
        galleryGraphStyle = "width:410px;height:290px;"

        graphCCtitle = "Correlations"
        graphCC = gallery.addFlotGraph(
            xmlnode=xmlnode,
            title=graphCCtitle,
            internalId=graphCCtitle,
            outputXml=self.outputXml,
            label=graphCCtitle,
            style=galleryGraphStyle)
        n_icols = 0
        graphCC.addData(title="Resolution(&Aring;)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
        if len(xmlnode.findall('.//cycle[last()]/data/binned/fsc_FC_full')) > 0:  # SPA refinement
            graphCC.addData(title="fsc_FC_full", select=".//cycle[last()]/data/binned/./fsc_FC_full")
            n_icols += 1
            if len(xmlnode.findall('.//cycle[last()]/data/binned/cc_FC_full')) > 0:
                graphCC.addData(title="CC_FC_full", select=".//cycle[last()]/data/binned/./cc_FC_full")
                n_icols += 1
            if len(xmlnode.findall('.//cycle[last()]/data/binned/mcos_FC_full')) > 0:
                graphCC.addData(title="mcos_FC_full", select=".//cycle[last()]/data/binned/./mcos_FC_full")
                n_icols += 1
        elif len(xmlnode.findall('.//cycle[last()]/data/binned/CCI')) > 0:
            graphCC.addData(title="CCI", select=".//cycle[last()]/data/binned/./CCI")
            n_icols += 1
        elif len(xmlnode.findall('.//cycle[last()]/data/binned/CCF')) > 0:
            graphCC.addData(title="CCF", select=".//cycle[last()]/data/binned/./CCF")
            n_icols += 1
        elif len(xmlnode.findall('.//cycle[last()]/data/binned/CCFwork')) > 0:
            graphCC.addData(title="CCFwork", select=".//cycle[last()]/data/binned/./CCFwork")
            n_icols += 1
            if len(xmlnode.findall('.//cycle[last()]/data/binned/CCFfree')) > 0:
                graphCC.addData(title="CCFfree", select=".//cycle[last()]/data/binned/./CCFfree")
                n_icols += 1
        else:
            graphCC.addData(title="CCIwork", select=".//cycle[last()]/data/binned/./CCIwork")
            n_icols += 1
            if len(xmlnode.findall('.//cycle[last()]/data/binned/CCIfree')) > 0:
                graphCC.addData(title="CCIfree", select=".//cycle[last()]/data/binned/./CCIfree")
                n_icols += 1
        if len(xmlnode.findall('.//cycle[last()]/data/binned/CC')) > 0:
            graphCC.addData(title="CC*", select=".//cycle[last()]/data/binned/./CC")
            n_icols += 1
        plotCC = graphCC.addPlotObject()
        plotCC.append('title', graphCCtitle)
        plotCC.append('plottype', 'xy')
        plotCC.append('xlabel', 'Resolution (&Aring;)')
        plotCC.append('ylabel', 'Correlation')
        plotCC.append('yrange', min=0.0, max=1.0)
        plotCC.append('xscale', 'oneoversqrt')
        plotCC.append('legendposition', x=0, y=0)
        colours = ['orange', 'blue', 'gray']
        for i in range(n_icols):
            plotLine = plotCC.append('plotline', xcol=1, ycol=2 + i)
            plotLine.append('colour', colours[i])
            plotLine.append('symbolsize', '0')

        # R-values vs. resolution - only for servalcat_xtal_norefmac
        if len(xmlnode.findall('.//cycle[last()]/data/binned/R1')) > 0 or \
                len(xmlnode.findall('.//cycle[last()]/data/binned/R')) > 0 or \
                len(xmlnode.findall('.//cycle[last()]/data/binned/Rwork')) > 0 or \
                len(xmlnode.findall('.//cycle[last()]/data/binned/R1work')) > 0:
            graphRtitle = "R-values"
            graphR = gallery.addFlotGraph(
                xmlnode=xmlnode,
                title=graphRtitle,
                internalId=graphRtitle,
                outputXml=self.outputXml,
                label=graphRtitle,
                style=galleryGraphStyle,
                initiallyDrawn=True)
            plotR = graphR.addPlotObject()
            graphR.addData(title="Resolution(&Aring;)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
            if len(xmlnode.findall('.//cycle[last()]/data/binned/Rcmplx_FC_full')) > 0: # SPA refinement
                graphR.addData(title="Rcmplx_FC_full", select=".//cycle[last()]/data/binned/./Rcmplx_FC_full")
            else:
                plotR.append('line', x1=0, x2=1, y1=0.42, y2=0.42, linecolour="red", linestyle="--")
                plotR.append('line', x1=0, x2=1, y1=0.58, y2=0.58, linecolour="red", linestyle="--")
                if len(xmlnode.findall('.//cycle[last()]/data/binned/R1')) > 0:
                    graphR.addData(title="R1", select=".//cycle[last()]/data/binned/./R1")
                elif len(xmlnode.findall('.//cycle[last()]/data/binned/R')) > 0:
                    graphR.addData(title="R", select=".//cycle[last()]/data/binned/./R")
                elif len(xmlnode.findall('.//cycle[last()]/data/binned/Rwork')) > 0:
                    graphR.addData(title="Rwork", select=".//cycle[last()]/data/binned/./Rwork")
                    if len(xmlnode.findall('.//cycle[last()]/data/binned/Rfree')) > 0:
                        graphR.addData(title="Rfree", select=".//cycle[last()]/data/binned/./Rfree")
                elif len(xmlnode.findall('.//cycle[last()]/data/binned/R1work')) > 0:
                    graphR.addData(title="R1work", select=".//cycle[last()]/data/binned/./R1work")
                    if len(xmlnode.findall('.//cycle[last()]/data/binned/R1free')) > 0:
                        graphR.addData(title="R1free", select=".//cycle[last()]/data/binned/./R1free")
            plotR.append('title', graphRtitle)
            plotR.append('plottype', 'xy')
            plotR.append('xlabel', 'Resolution (&Aring;)')
            plotR.append('ylabel', 'R-value')
            plotR.append('yrange', min=0.0)
            plotR.append('xscale', 'oneoversqrt')
            plotR.append('legendposition', x=1, y=0)  # right bottom corner
            plotLine = plotR.append('plotline', xcol=1, ycol=2)
            plotLine.append('colour', 'orange')
            plotLine.append('symbolsize', '0')
            plotLine = plotR.append('plotline', xcol=1, ycol=3)
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')

        # Number of reflections - only for servalcat_xtal_norefmac
        if len(xmlnode.findall('.//cycle[last()]/data/binned/n_obs')) > 0:
            graphNtitle = "Number of reflections"
            graphN = gallery.addFlotGraph(
                xmlnode=xmlnode,
                title=graphNtitle,
                internalId=graphNtitle,
                outputXml=self.outputXml,
                label=graphNtitle,
                style=galleryGraphStyle)
            graphN.addData(title="Resolution(&Aring;)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
            graphN.addData(title="Nobs", select=".//cycle[last()]/data/binned/./n_obs")
            avail_n_work = False
            avail_n_free = False
            avail_n_R1work = False
            avail_n_R1free = False
            avail_n_R1 = False
            if len(xmlnode.findall('.//cycle[last()]/data/binned/n_R1')) > 0:
                graphN.addData(title="N_R1", select=".//cycle[last()]/data/binned/./n_R1")
                avail_n_R1 = True
            else:
                if len(xmlnode.findall('.//cycle[last()]/data/binned/n_work')) > 0:
                    graphN.addData(title="Nwork", select=".//cycle[last()]/data/binned/./n_work")
                    avail_n_work = True
                    if len(xmlnode.findall('.//cycle[last()]/data/binned/n_free')) > 0:
                        graphN.addData(title="Nfree", select=".//cycle[last()]/data/binned/./n_free")
                        avail_n_free = True
                    if len(xmlnode.findall('.//cycle[last()]/data/binned/n_R1work')) > 0:
                        graphN.addData(title="N_R1work", select=".//cycle[last()]/data/binned/./n_R1work")
                        avail_n_R1work = True
                    if len(xmlnode.findall('.//cycle[last()]/data/binned/n_R1free')) > 0:
                        graphN.addData(title="N_R1free", select=".//cycle[last()]/data/binned/./n_R1free")
                        avail_n_R1free = True

            plotN = graphN.addPlotObject()
            plotN.append('title', graphNtitle)
            plotN.append('plottype', 'xy')
            plotN.append('xlabel', 'Resolution (&Aring;)')
            plotN.append('legendposition', x=0, y=1)
            plotN.append('xscale', 'oneoversqrt')
            plotLine = plotN.append('plotline', xcol=1, ycol=2)
            plotLine.append('colour', 'gray')
            plotLine.append('symbolsize', '0')
            if avail_n_R1:
                plotLine = plotN.append('plotline', xcol=1, ycol=3)
                plotLine.append('colour', 'orange')
                plotLine.append('symbolsize', '0')
            else:
                if avail_n_work:
                    plotLine = plotN.append('plotline', xcol=1, ycol=3)
                    plotLine.append('colour', 'orange')
                    plotLine.append('symbolsize', '0')
                if avail_n_free:
                    plotLine = plotN.append('plotline', xcol=1, ycol=4) # , rightaxis='true')
                    plotLine.append('colour', 'blue')
                    plotLine.append('symbolsize', '0')
                if avail_n_R1work:
                    plotLine = plotN.append('plotline', xcol=1, ycol=5)
                    plotLine.append('colour', 'red')
                    plotLine.append('symbolsize', '0')
                if avail_n_R1free:
                    plotLine = plotN.append('plotline', xcol=1, ycol=6)
                    plotLine.append('colour', 'cyan')
                    plotLine.append('symbolsize', '0')

            if avail_n_work and avail_n_free and avail_n_R1work and avail_n_R1free:
                # plot of only n_obs n_work n_free
                plotN2 = graphN.addPlotObject()
                plotN2.append('title', "Number of reflections (only Nobs and Nwork and Nfree)")
                plotN2.append('plottype', 'xy')
                plotN2.append('xlabel', 'Resolution (&Aring;)')
                plotN2.append('legendposition', x=0, y=1)
                plotN2.append('xscale', 'oneoversqrt')
                plotLine = plotN2.append('plotline', xcol=1, ycol=2)
                plotLine.append('colour', 'gray')
                plotLine.append('symbolsize', '0')
                plotLine = plotN2.append('plotline', xcol=1, ycol=3)
                plotLine.append('colour', 'orange')
                plotLine.append('symbolsize', '0')
                plotLine = plotN2.append('plotline', xcol=1, ycol=4) # , rightaxis='true')
                plotLine.append('colour', 'blue')
                plotLine.append('symbolsize', '0')

            if avail_n_R1work and avail_n_R1free:
                plotN3 = graphN.addPlotObject()
                plotN3.append('title', "Number of reflections (only N_R1work and N_R1free)")
                plotN3.append('plottype', 'xy')
                plotN3.append('xlabel', 'Resolution (&Aring;)')
                plotN3.append('legendposition', x=0, y=1)
                plotN3.append('xscale', 'oneoversqrt')
                plotLine = plotN3.append('plotline', xcol=1, ycol=5)
                plotLine.append('colour', 'red')
                plotLine.append('symbolsize', '0')
                plotLine = plotN3.append('plotline', xcol=1, ycol=6)
                plotLine.append('colour', 'cyan')
                plotLine.append('symbolsize', '0')

            if avail_n_work and avail_n_R1work:
                plotN4 = graphN.addPlotObject()
                plotN4.append('title', "Number of reflections (only Nwork and N_R1work)")
                plotN4.append('plottype', 'xy')
                plotN4.append('xlabel', 'Resolution (&Aring;)')
                plotN4.append('legendposition', x=0, y=1)
                plotN4.append('xscale', 'oneoversqrt')
                plotLine = plotN4.append('plotline', xcol=1, ycol=2)
                plotLine.append('colour', 'gray')
                plotLine.append('symbolsize', '0')
                plotLine = plotN4.append('plotline', xcol=1, ycol=3)
                plotLine.append('colour', 'orange')
                plotLine.append('symbolsize', '0')
                plotLine = plotN4.append('plotline', xcol=1, ycol=5)
                plotLine.append('colour', 'red')
                plotLine.append('symbolsize', '0')

            if avail_n_free and avail_n_R1free:
                plotN5 = graphN.addPlotObject()
                plotN5.append('title', "Number of reflections (only Nfree and N_R1free)")
                plotN5.append('plottype', 'xy')
                plotN5.append('xlabel', 'Resolution (&Aring;)')
                plotN5.append('legendposition', x=0, y=1)
                plotN5.append('xscale', 'oneoversqrt')
                plotLine = plotN5.append('plotline', xcol=1, ycol=4)
                plotLine.append('colour', 'blue')
                plotLine.append('symbolsize', '0')
                plotLine = plotN5.append('plotline', xcol=1, ycol=6)
                plotLine.append('colour', 'cyan')
                plotLine.append('symbolsize', '0')

        # Completeness - only for servalcat_xtal_norefmac
        if len(xmlnode.findall('.//cycle[last()]/data/binned/Cmpl')) > 0:
            graphCmplTitle = "Completeness (%)"
            graphCmpl = gallery.addFlotGraph(
                xmlnode=xmlnode,
                title=graphCmplTitle,
                internalId=graphCmplTitle,
                outputXml=self.outputXml,
                label=graphCmplTitle,
                style=galleryGraphStyle)
            graphCmpl.addData(title="Resolution(&Aring;)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
            graphCmpl.addData(title="Completeness(%)", select=".//cycle[last()]/data/binned/./Cmpl")
            plotCmpl = graphCmpl.addPlotObject()
            plotCmpl.append('title', graphCmplTitle)
            plotCmpl.append('plottype', 'xy')
            plotCmpl.append('xlabel', 'Resolution (&Aring;)')
            plotCmpl.append('legendposition', x=0, y=0)
            plotCmpl.append('xscale', 'oneoversqrt')
            plotCmpl.append('yrange', min=0.0, max=100.0)
            plotLine = plotCmpl.append('plotline', xcol=1, ycol=2)
            plotLine.append('colour', 'orange')
            plotLine.append('symbolsize', '0')

        # MnIo & MnIc or MnFo & MnFc - only for servalcat_xtal_norefmac
        if (
            len(xmlnode.findall('.//cycle[last()]/data/binned/MnIo')) > 0
            or len(xmlnode.findall('.//cycle[last()]/data/binned/MnFo')) > 0
        ):
            if len(xmlnode.findall('.//cycle[last()]/data/binned/MnIo')) > 0 and \
                    len(xmlnode.findall('.//cycle[last()]/data/binned/MnIc')) > 0:
                graphMnOCTitle = "Mean Io and Ic"
                MnO = "MnIo"
                MnC = "MnIc"
            elif len(xmlnode.findall('.//cycle[last()]/data/binned/MnFo')) > 0 and \
                    len(xmlnode.findall('.//cycle[last()]/data/binned/MnFc')) > 0:
                graphMnOCTitle = "Mean Fo and Fc"
                MnO = "MnFo"
                MnC = "MnFc"
            graphMnOC = gallery.addFlotGraph(
                xmlnode=xmlnode,
                title=graphMnOCTitle,
                internalId=graphMnOCTitle,
                outputXml=self.outputXml,
                label=graphMnOCTitle,
                style=galleryGraphStyle)
            graphMnOC.addData(title="Resolution(&Aring;)", select=".//cycle[last()]/data/binned/./d_min_4ssqll")
            graphMnOC.addData(title=MnO, select=f".//cycle[last()]/data/binned/./{MnO}")
            graphMnOC.addData(title=MnC, select=f".//cycle[last()]/data/binned/./{MnC}")
            plotMnIoIc = graphMnOC.addPlotObject()
            plotMnIoIc.append('title', graphMnOCTitle)
            plotMnIoIc.append('plottype', 'xy')
            plotMnIoIc.append('xlabel', 'Resolution (&Aring;)')
            plotMnIoIc.append('xscale', 'oneoversqrt')
            plotMnIoIc.append('yrange', min=0.0)
            plotMnIoIc.append('legendposition', x=1, y=1)
            plotLine = plotMnIoIc.append('plotline', xcol=1, ycol=2)
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')
            plotLine = plotMnIoIc.append('plotline', xcol=1, ycol=3)
            plotLine.append('colour', 'red')
            plotLine.append('symbolsize', '0')

        #  MnD0FC0, MnD1FCbulk - only for servalcat_xtal_norefmac
        MnD_parent = ""
        MnD_parents = ["binned", "ml"]
        for p in MnD_parents:
            if len(xmlnode.findall(f'.//cycle[last()]/data/{p}/MnD0FC0')) > 0 and \
                    len(xmlnode.findall(f'.//cycle[last()]/data/{p}/MnD1FCbulk')) > 0:
                MnD_parent = p
        if MnD_parent:
            graphDtitle = "Mean |D0*FC0| and |D1*FCbulk|"
            graphD = gallery.addFlotGraph(
                xmlnode=xmlnode,
                title=graphDtitle,
                internalId=graphDtitle,
                outputXml=self.outputXml,
                label=graphDtitle,
                style=galleryGraphStyle)
            graphD.addData(title="Resolution(&Aring;)", select=f".//cycle[last()]/data/{MnD_parent}/./d_min_4ssqll")
            graphD.addData(title="Mean|D0*FC0|", select=f".//cycle[last()]/data/{MnD_parent}/./MnD0FC0")
            graphD.addData(title="Mean|D1*FCbulk|", select=f".//cycle[last()]/data/{MnD_parent}/./MnD1FCbulk")
            plotD = graphD.addPlotObject()
            plotD.append('title', graphDtitle)
            plotD.append('plottype', 'xy')
            plotD.append('xlabel', 'Resolution (&Aring;)')
            plotD.append('xscale', 'oneoversqrt')
            plotD.append('yrange', min=0.0)
            plotD.append('legendposition', x=1, y=1)
            plotLine = plotD.append('plotline', xcol=1, ycol=2)
            plotLine.append('colour', 'blue')
            plotLine.append('symbolsize', '0')
            plotD.append('yrange', rightaxis='true')
            plotLine = plotD.append('plotline', xcol=1, ycol=3, rightaxis='true')
            plotLine.append('colour', 'red')
            plotLine.append('symbolsize', '0')

        clearingDiv = parent.addDiv(style="clear:both;")

    def addOutlierAnalysis(self, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        # All possible keys are for now: 'bond', 'angle', 'torsion', 'chir', 'plane', 'staca', 'stacd', 'vdw'
        # ADP will be added
        # 'per_atom' not to be used
        outlierFold = parent.addFold(label="Outliers identified by Servalcat", brief='Outliers')
        outBond = xmlnode.findall('.//cycle[last()]/geom/outliers/bond')
        outAngle = xmlnode.findall('.//cycle[last()]/geom/outliers/angle')
        outTorsion = xmlnode.findall('.//cycle[last()]/geom/outliers/torsion')
        outChir = xmlnode.findall('.//cycle[last()]/geom/outliers/chir')
        outPlane = xmlnode.findall('.//cycle[last()]/geom/outliers/plane')
        outStaca = xmlnode.findall('.//cycle[last()]/geom/outliers/staca')
        outStacd = xmlnode.findall('.//cycle[last()]/geom/outliers/stacd')
        outVdw = xmlnode.findall('.//cycle[last()]/geom/outliers/vdw')

        if len(outVdw) > 0:
            n_outliers = len(outVdw)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Van der Waals repulsion outliers indicating close contacts (clashes) between non-bonding atoms:")
            outData = {'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'z_abs': ["-"]*n_outliers,
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
                try:
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    outData['z_abs'][i] = round(abs(z), 2)
                except:
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
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

            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Distance (&Aring;)", data=outData['value'])
            fullTable.addData(title="Critical<br>distance (&Aring;)", data=outData['ideal'])
            fullTable.addData(title="Difference from<br>critical (&Aring;)", data=outData['difference'])
            fullTable.addData(title="Z", data=outData['z'])
            fullTable.addData(title="Type", data=outData['note'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No clashes between atoms (Van der Waals repulsion outliers) observed.")

        if len(outBond) > 0:
            n_outliers = len(outBond)
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("Bond length outliers:")
            outData = {'atom1': ["-"]*n_outliers,
                       'atom2': ["-"]*n_outliers,
                       'value': ["-"]*n_outliers,
                       'ideal': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = round(abs(z), 2)
                    difference = abs(value - ideal)
                    outData['difference'][i] = "{:.2f}".format(difference)
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = round(sigma, 2)
                    percent = 100 * abs(value - ideal) / ideal
                    outData['percent'][i] = round(percent, 2)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['sigma'][i] = '-'
                    outData['percent'][i] = '-'
            outDataZip = list(zip(outData['z_abs'], outData['percent'], outData['sigma'], outData['type'], outData['atom1'], outData['atom2'],
                                  outData['value'], outData['ideal'], outData['difference'], outData['z']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['percent'], outData['sigma'], outData['type'], outData['atom1'], outData['atom2'], \
                outData['value'], outData['ideal'], outData['difference'], outData['z'] = zip(*outDataZip)
    
            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Deviation<br>(in %)", data=outData['percent'])
            fullTable.addData(title="Bond<br>length (&Aring;)", data=outData['value'])
            fullTable.addData(title="Ideal<br>length (&Aring;)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (&Aring;)", data=outData['difference'])
            fullTable.addData(title="Sigma (&Aring;)", data=outData['sigma'])
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
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = round(abs(z), 2)
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = round(sigma, 2)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['sigma'], outData['atom1'], outData['atom2'],
                                  outData['atom3'], outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['sigma'], outData['atom1'], outData['atom2'], \
                outData['atom3'], outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
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
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = round(abs(z), 2)
                    periodicity = int(outlier.findall('per')[0].text)
                    outData['per'][i] = periodicity
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
                    outData['z_abs'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'

            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['label'], outData['atom1'], outData['atom2'], outData['atom3'],
                                  outData['atom4'], outData['value'], outData['ideal'], outData['per'], outData['ideal_per'], outData['z'], outData['difference'], outData['sigma']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['label'], outData['atom1'], outData['atom2'], outData['atom3'], \
                outData['atom4'], outData['value'], outData['ideal'], outData['per'], outData['ideal_per'], outData['z'], outData['difference'], outData['sigma'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
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
                       'signum': ["-"]*n_outliers,
                       'z': ["-"]*n_outliers,
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = round(abs(z), 2)
                    if bool(outlier.findall('both')[0].text) == "True":
                        both = True
                    else:
                        both = False
                    outData['both'][i] = str(both)
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                    if int(sign(value) * sign(ideal)) != 1 and not both:
                        signum = "No"  # "Ideal value has an opposite sign"
                    else:
                        signum = 'Yes'
                    outData['signum'][i] = signum
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
                    outData['both'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
                    outData['signum'][i] = '-'
            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['atomc'], outData['atom1'], outData['atom2'], outData['atom3'],
                                  outData['value'], outData['ideal'], outData['z'], outData['both'], outData['difference'], outData['sigma'], outData['signum']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['atomc'], outData['atom1'], outData['atom2'], outData['atom3'], outData['value'], outData['ideal'], outData['z'], outData['both'], outData['difference'], outData['sigma'], outData['signum'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Chiral atom", data=outData['atomc'])
            fullTable.addData(title="Atom 1", data=outData['atom1'])
            fullTable.addData(title="Atom 2", data=outData['atom2'])
            fullTable.addData(title="Atom 3", data=outData['atom3'])
            fullTable.addData(title="Chiral<br>volume (&Aring;<sup>3</sup>)", data=outData['value'])
            fullTable.addData(title="Ideal<br>value (&Aring;<sup>3</sup>)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (&Aring;<sup>3</sup>)", data=outData['difference'])
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
                       'z': ["-"]*n_outliers,
                       'z_abs': ["-"]*n_outliers}
            for i, outlier in enumerate(outPlane):
                try: outData['label'][i] = str(outlier.findall('label')[0].text)
                except: outData['label'][i] = '-'
                try: outData['atom'][i] = str(outlier.findall('atom')[0].text)
                except: outData['atom'][i] = '-'
                try: outData['dev'][i] = "{:.2f}".format(float(outlier.findall('dev')[0].text))
                except: outData['dev'][i] = '-'
                try:
                    z = float(outlier.findall('z')[0].text)
                    outData['z'][i] = "{:.2f}".format(z)
                    outData['z_abs'][i] = round(abs(z), 2)
                except:
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
            # Does not need to be sorted
            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Label", data=outData['label'])
            fullTable.addData(title="Atom", data=outData['atom'])
            fullTable.addData(title="Deviation (&Aring;)", data=outData['dev'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No planarity outliers observed.")

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
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = round(abs(z), 2)
                    difference = abs(value - ideal)
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                    sigma = abs((value - ideal) / z)
                    outData['sigma'][i] = "{:.2f}".format(sigma)
                except:
                    outData['value'][i] = '-'
                    outData['ideal'][i] = '-'
                    outData['z'][i] = '-'
                    outData['z_abs'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Plane 1", data=outData['plane1'])
            fullTable.addData(title="Plane 2", data=outData['plane2'])
            fullTable.addData(title="Stacking<br>distance (&Aring;)", data=outData['value'])
            fullTable.addData(title="Ideal<br>distance (&Aring;)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (&Aring;)", data=outData['difference'])
            fullTable.addData(title="Sigma (&Aring;)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No stacking distance outliers observed.")

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
                       'z_abs': ["-"]*n_outliers,
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
                    outData['z_abs'][i] = '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'
                    outData['sigma'][i] = '-'
            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['sigma'], outData['plane1'], outData['plane2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv(style="clear:both;")
            styleDiv = outlierFold.addDiv(style="color:navy; text-align: right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Plane 1", data=outData['plane1'])
            fullTable.addData(title="Plane 2", data=outData['plane2'])
            fullTable.addData(title="Stacking<br>angle (°)", data=outData['value'])
            fullTable.addData(title="Ideal<br>angle (°)", data=outData['ideal'])
            fullTable.addData(title="Difference<br>from ideal (°)", data=outData['difference'])
            fullTable.addData(title="Sigma (°)", data=outData['sigma'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style='font-size:110%')
            div.append("No stacking angle outliers observed.")


    def addScrollableDownloadableTable1(self, cycle_data, xmlnode=None, parent=None, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        #create a "shell" div to contain the scrollable table and the hyperlink
        scrollableDownloadableTableDiv = parent.addDiv(style="height:250px; width:315px;float:left;margin-top:2px;")
        #place a scrollable div into the shell: the table will be inserted into this div
        scrollableTableDiv = scrollableDownloadableTableDiv.addDiv(style="height:225px; width:300px;clear:both;overflow:auto;")
        #Put table1 into this (autoscrolling) div
        table1 = self.addTable1(cycle_data, parent=scrollableTableDiv,internalId=internalId)
        #Add a hyperlink to this table in the "tables_as_csv_files" directory
        download = scrollableDownloadableTableDiv.addDownload(jobInfo=self.jobInfo,dataName=table1.id)

    def addTable1(self, cycle_data, xmlnode=None, parent=None, downloadable=False, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        table1 = parent.addTable(xmlnode=xmlnode, style="width:240px;float:left;", downloadable=downloadable, outputXml=self.outputXml, internalId=internalId)
        
        statisticNames = []
        statisticInitial = []
        statisticFinal = []

        if isnumber(cycle_data["FSCaverage"][-1]):# SPA refinement
            statisticNames.append('⟨FSCmodel⟩')
            statisticInitial.append(cycle_data['FSCaverage'][0])
        elif isnumber(cycle_data["R1"][-1]):      # refinement against intensities without free flags
            statisticNames.append('R1')
            statisticInitial.append(cycle_data['R1'][0])
            statisticFinal.append(cycle_data['R1'][-1])
            statisticNames.append('⟨CCI⟩')
            statisticInitial.append(cycle_data['CCIavg'][0])
            statisticFinal.append(cycle_data['CCIavg'][-1])
        elif isnumber(cycle_data["R"][-1]):       # refinement against amplitudes without free flags
            statisticNames.append('R')
            statisticInitial.append(cycle_data['R'][0])
            statisticFinal.append(cycle_data['R'][-1])
            statisticNames.append('⟨CCF⟩')
            statisticInitial.append(cycle_data['CCFavg'][0])
            statisticFinal.append(cycle_data['CCFavg'][-1])
        elif isnumber(cycle_data["R1work"][-1]):  # refinement against intensities with free flags
            statisticNames.append('R1work')
            statisticInitial.append(cycle_data['R1work'][0])
            statisticFinal.append(cycle_data['R1work'][-1])
            if isnumber(cycle_data["R1free"][-1]):
                statisticNames.append('R1free')
                statisticInitial.append(cycle_data['R1free'][0])
                statisticFinal.append(cycle_data['R1free'][-1])
            statisticNames.append('⟨CCIwork⟩')
            statisticInitial.append(cycle_data['CCIworkavg'][0])
            statisticFinal.append(cycle_data['CCIworkavg'][-1])
            if isnumber(cycle_data["CCIfreeavg"][-1]):
                statisticNames.append('⟨CCIfree⟩')
                statisticInitial.append(cycle_data['CCIfreeavg'][0])
                statisticFinal.append(cycle_data['CCIfreeavg'][-1])
        else:                               # refinement against amplitudes with free flags
            statisticNames.append('Rwork')
            statisticInitial.append(cycle_data['Rwork'][0])
            statisticFinal.append(cycle_data['Rwork'][-1])
            if isnumber(cycle_data["Rfree"][-1]):
                statisticNames.append('Rfree')
                statisticInitial.append(cycle_data['Rfree'][0])
                statisticFinal.append(cycle_data['Rfree'][-1])
            statisticNames.append('⟨CCFwork⟩')
            statisticInitial.append(cycle_data['CCFworkavg'][0])
            statisticFinal.append(cycle_data['CCFworkavg'][-1])
            if isnumber(cycle_data["CCFfreeavg"][-1]):
                statisticNames.append('⟨CCFfree⟩')
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


def addCorrelationProgress(progressGraph):
    plot = progressGraph.addPlotObject()
    plot.append('title', 'Correlations')
    plot.append('plottype', 'xy')
    plot.append('xlabel', 'Cycle')
    plot.append('ylabel', 'Correlation')
    plot.append('yrange', min=0.0, max=1.0)
    plot.append('xintegral', 'true')
    plot.append('legendposition', x=0, y=0)
    line = plot.append('plotline', xcol=1, ycol=4)
    line.append('colour', 'orange')
    line.append('symbolsize', '0')
    line = plot.append('plotline', xcol=1, ycol=6)
    line.append('colour', 'blue')
    line.append('symbolsize', '0')


def addRValuesProgress(progressGraph):
    plot = progressGraph.addPlotObject()
    plot.append('title', 'R-values')
    plot.append('plottype', 'xy')
    plot.append('xlabel', 'Cycle')
    plot.append('ylabel', 'R-value')
    plot.append('yrange', min=0.0)
    plot.append('xintegral', 'true')
    plot.append('legendposition', x=0, y=0)
    line = plot.append('plotline', xcol=1, ycol=3)
    line.append('colour', 'orange')
    line.append('symbolsize', '0')
    line = plot.append('plotline', xcol=1, ycol=5)
    line.append('colour', 'blue')
    line.append('symbolsize', '0')


if __name__ == "__main__":
    import sys
    servalcat_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
