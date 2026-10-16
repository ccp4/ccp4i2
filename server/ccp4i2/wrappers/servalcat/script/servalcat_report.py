import json
import os
import sys

from numpy import sign

from ccp4i2.report import Report


def _r_and_cc(xmlnode):
    """The R and CC servalcat reported in the last cycle, or None (SPA).

    Returns (R label, free R label or None, 'I' or 'F'). Servalcat up to
    0.4.149 reports R1work/R1free/R1 when refining against intensities and
    Rwork/Rfree/R against amplitudes; from 0.4.150 it reports Rwork/Rfree/R
    for both. So the R name no longer says which data were used: the CC
    names do (CCI against intensities, CCF against amplitudes).
    """
    last = './/cycle[last()]/data/summary/'
    def has(key):
        return len(xmlnode.findall(last + key)) > 0
    cc = 'I' if has('CCIworkavg') or has('CCIavg') else 'F'
    for work in ('Rwork', 'R1work'):
        if has(work):
            return work, work.replace('work', 'free'), cc
    for overall in ('R', 'R1'):
        if has(overall):
            return overall, None, cc
    return None


_R_AND_CC_KEYS = ('Rwork', 'Rfree', 'R', 'R1work', 'R1free', 'R1',
                  'CCFworkavg', 'CCFfreeavg', 'CCFavg',
                  'CCIworkavg', 'CCIfreeavg', 'CCIavg')


def _r_and_cc_columns(cycle_data):
    """[(column title, cycle_data key), ...] for the R and CC servalcat
    reported, in display order; empty for SPA. As _r_and_cc, but read from
    getCycleData's table, whose last entry is the final cycle."""
    def final(key):
        return isnumber(cycle_data[key][-1])
    cc = 'I' if final('CCIworkavg') or final('CCIavg') else 'F'
    for work in ('Rwork', 'R1work'):
        if final(work):
            free = work.replace('work', 'free')
            columns = [(work, work)]
            if final(free):
                columns.append((free, free))
            columns.append((f'⟨CC{cc}work⟩', f'CC{cc}workavg'))
            if final(f'CC{cc}freeavg'):
                columns.append((f'⟨CC{cc}free⟩', f'CC{cc}freeavg'))
            return columns
    for overall in ('R', 'R1'):
        if final(overall):
            return [(overall, overall), (f'⟨CC{cc}⟩', f'CC{cc}avg')]
    return []


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
    TASKNAME = 'servalcat'
    RUNNING = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        # 'nooutput' mode would be used by another report class that wanted
        # to use some method(s) from this class for its own report
        self.outputXml = jobStatus is not None and jobStatus.lower().count('running')

        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return

        outputJsonPath = os.path.normpath(
            os.path.join(self.jobInfo["fileroot"], "refined_stats.json")
        )
        self.outputJson = json.load(open(outputJsonPath)) if os.path.isfile(outputJsonPath) else None

        self.addDiv()

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
            clearingDiv = parent.addDiv()
            perCycleFold = parent.addFold(label='Per cycle statistics', brief='Per cycle', initiallyOpen=False)
            self.addTablePerCycle(cycle_data, parent=perCycleFold, initialFinalOnly=False)
            self.addGraphsVsResolution()
            self.addTwinningAnalysis(outputJson=self.outputJson)
            self.addOutlierAnalysis()

    def addGraphPerCycle(self, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        
        # Note that when I add the progressgraph, I have to ensure that the select is rooted in my own xmlnode
        progressGraph = parent.addFlotGraph(
            title="Refinement results",
            xmlnode=self.xmlnode,
            style="height:250px;width:100%;")
        progressGraph.addData(title="Cycle", select=".//cycle/Ncyc") # ycol=1
        progressGraph.addData(title="-LL", select=".//cycle/data/summary/minusLL") # ycol=2
        spa_refinement = False
        intensity_based_refinement = False
        if len(xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')) > 0:
            progressGraph.addData(title="⟨FSCmodel⟩", select=".//cycle/data/summary/FSCaverage", expr="x if float(x)>=0.0 else ''")  # ycol=3
            spa_refinement = True
        elif _r_and_cc(xmlnode) is not None:
            r_label, r_free, cc = _r_and_cc(xmlnode)
            intensity_based_refinement = cc == 'I'
            cc_label = 'work' if r_free else ''
            progressGraph.addData(title=r_label, select=f".//cycle/data/summary/{r_label}", expr="x if float(x)>=0.0 else ''")  # ycol=3
            progressGraph.addData(title=f"⟨CC{cc}{cc_label}⟩", select=f".//cycle/data/summary/CC{cc}{cc_label}avg", expr="x if float(x)>=-1.0 else ''")  # ycol=4
            if r_free and len(xmlnode.findall(f'.//cycle[last()]/data/summary/{r_free}')) > 0:
                progressGraph.addData(title=r_free, select=f".//cycle/data/summary/{r_free}", expr="x if float(x)>=0.0 else '-'")  # ycol=5
                progressGraph.addData(title=f"⟨CC{cc}free⟩", select=f".//cycle/data/summary/CC{cc}freeavg", expr="x if float(x)>=-1.0 else '-'")  # ycol=6

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

        clearingDiv = parent.addDiv()
        clearingDiv.append('&#160;&#160;')

        if len(xmlnode.findall('.//cycle/geom/summary/rmsd/Bond_distances_non_H')) > 0:
            progressGraph2 = parent.addFlotGraph(
                title="Deviation from standard geometry",
                xmlnode=self.xmlnode,
                style="height:250px;width:100%;")
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

        clearingDiv = parent.addDiv()
        if len(xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')) > 0:
            return
        noteDiv = parent.addDiv()
        note = ""
        r_and_cc = _r_and_cc(xmlnode)
        if r_and_cc is not None:
            r_label, _, cc = r_and_cc
            if r_label.startswith('R1'):
                note = "R1 is the R-value calculated from the square root of intensities where <i>I</i>/&#963;(<i>I</i>) > 2.<br />"
            elif cc == 'I':
                # Servalcat 0.4.150 and later: a conventional R, on amplitudes
                # derived from the intensities, in place of R1.
                note = "R is calculated against amplitudes derived from the intensities.<br />"
            if cc == 'I':
                note += "⟨CC<sub><i>I</i></sub>⟩ is the correlation coefficient between calculated and observed intensities averaged over resolution shells."
            else:
                note += "⟨CC<sub><i>F</i></sub>⟩ is the correlation coefficient between calculated and observed amplitudes averaged over resolution shells."
        noteDiv.append(note)

    def getCycleData(self, xmlnode=None):
        if xmlnode is None: xmlnode = self.xmlnode
        FSCaverageNodes = xmlnode.findall('.//cycle[last()]/data/summary/FSCaverage')
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
            else:
                # Whichever R and CC servalcat reported: the names changed in
                # 0.4.150 (see _r_and_cc), so read them all and let the
                # tables choose.
                for key in _R_AND_CC_KEYS:
                    try: cycle_data[key][idx] = "{:.4f}".format(float(cycle.findall('data/summary/' + key)[0].text))
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
        if not cycle_data.get('cycle'):  # no <cycle> in the XML: refinement stopped before cycle 0
            return
        clearingDiv = parent.addDiv()
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

        clearingDiv = parent.addDiv()
        fullTable = None
        fullTable = parent.addTable()
        if len(TLSIdx) > 0 and len(RestrIdx) > 0:
           mode = ['TLS']*len(TLSIdx) + ['Full Atom']*len(RestrIdx)
        fullTable.addData(title="Cycle", data=cycle_data_sel['cycle'])
        if isnumber(cycle_data_sel['FSCaverage'][-1]):  # SPA refinement
            fullTable.addData(title="⟨FSCmodel⟩", data=cycle_data_sel['FSCaverage'])
        else:
            for title, key in _r_and_cc_columns(cycle_data_sel):
                fullTable.addData(title=title, data=cycle_data_sel[key])
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

        clearingDiv = parent.addDiv()


    def addTwinningAnalysis(self, outputJson=None, parent=None, xmlnode=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if outputJson is None: outputJson = self.outputJson
        # if len(xmlnode.findall('.//cycle[last()]/twin_alpha')) == 0:
        if not outputJson or outputJson[-1].get('twin_alpha', {}) == {}:
            # No twinning analysis
            return
        twinFold = parent.addFold(label="Twinning analysis", brief='Twinning', initiallyOpen=True)
        divLeft = twinFold.addDiv(style="font-size:110%;")

        try:
            twin_alpha_final = outputJson[-1]['twin_alpha']
            twin_operators = list(twin_alpha_final.keys())
            twin_fraction_final_values = list(twin_alpha_final.values())
            twin_fraction_final_values = ["{:.2f}".format(v) for v in twin_fraction_final_values]
            twin_operators_labels = twin_operators
            twin_alpha_data = {op: [] for op in twin_operators}
            cycles_list = list(range(1, len(outputJson) + 1))
            for c in outputJson:
                for op in twin_operators:
                    twin_alpha_data[op].append(c.get('twin_alpha', {}).get(op, 0.0))

            twinGraph = divLeft.addFlotGraph(
                title="Twin fraction vs cycle",
                xmlnode=self.xmlnode,
                style="height:250px;width:400px;float:left;")
            twinGraph.addData(title="Cycle", data=cycles_list)
            for i, twin_op in enumerate(twin_operators):
                twinGraph.addData(
                    title=twin_operators_labels[i],
                    data=twin_alpha_data[twin_op],
                )
            plotTwin = twinGraph.addPlotObject()
            plotTwin.append('title', 'Twin fraction vs cycle')
            plotTwin.append('plottype', 'xy')
            plotTwin.append('xlabel', 'Cycle')
            plotTwin.append('ylabel', 'Twin fraction vs cycle')
            plotTwin.append('xintegral', 'true')
            plotTwin.append('legendposition', x=0, y=0)
            plotTwin.append('yrange', min=0.0, max=1.0)
            for i, twin_op in enumerate(twin_operators):
                plotLine = plotTwin.append('plotline', xcol=1, ycol=i + 2)
                plotLine.append('symbolsize', '0')

            divRight = twinFold.addDiv(style="font-size:110%;")
            divRight.append("Final twin fractions after refinement:")
            tableTwin = divRight.addTable()
            tableTwin.addData(title="Twin operator", data=twin_operators_labels)
            tableTwin.addData(title="Twin fraction", data=twin_fraction_final_values)
        except:
            divLeft.append("Error: Twinning report was not found.")

        clearingDiv = parent.addDiv()

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
            div = outlierFold.addDiv(style="font-size:110%;")
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

                # Outlier type
                outType = outlier.findall('type')[0].text
                if (isinstance(outType, int)) or (isinstance(outType, str) and outType.isdigit()):
                    # before Servalcat 0.4.123
                    # convert outlier type given as integer to description
                    try:
                        outType = int(outType)
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
                else:
                    # Servalcat 0.4.123 and later
                    try:
                        outData['type'][i] = outType
                        outData['note'][i] = outType
                    except:
                        outData['type'][i] = '-'
                        outData['note'][i] = '-'

                try:
                    # difference = | value - ideal |
                    difference = abs(float(outlier.findall('value')[0].text) - float(outlier.findall('ideal')[0].text))
                    outData['difference_float'][i] = difference
                    outData['difference'][i] = "{:.2f}".format(difference)
                except:
                    # This outlier's entry only: assigning the whole list
                    # replaced every difference with one '-'
                    outData['difference'][i] = '-'
                    outData['difference_float'][i] = '-'

            outDataZip = list(zip(outData['z_abs'], outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'],
                                  outData['value'], outData['ideal'], outData['z'], outData['difference']))
            outDataZip.sort(reverse=True)
            outData['z_abs'], outData['difference_float'], outData['type'], outData['atom1'], outData['atom2'], \
                outData['value'], outData['ideal'], outData['z'], outData['difference'] = zip(*outDataZip)

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No clashes between atoms (Van der Waals repulsion outliers) observed.")

        if len(outBond) > 0:
            n_outliers = len(outBond)
            div = outlierFold.addDiv(style="font-size:110%;")
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
    
            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No bond length outliers observed.")
            
        if len(outAngle) > 0:
            n_outliers = len(outAngle)
            div = outlierFold.addDiv(style="font-size:110%;")
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

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No bond angle outliers observed.")

        if len(outTorsion) > 0:
            n_outliers = len(outTorsion)
            div = outlierFold.addDiv(style="font-size:110%;")
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

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No torsion angle outliers observed.")

        if len(outChir) > 0:
            n_outliers = len(outChir)
            div = outlierFold.addDiv(style="font-size:110%;")
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

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No chirality outliers observed.")

        if len(outPlane) > 0:
            n_outliers = len(outPlane)
            div = outlierFold.addDiv(style="font-size:110%;")
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
            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
            fullTable = None
            fullTable = styleDiv.addTable()
            fullTable.addData(title="Label", data=outData['label'])
            fullTable.addData(title="Atom", data=outData['atom'])
            fullTable.addData(title="Deviation (&Aring;)", data=outData['dev'])
            fullTable.addData(title="Z", data=outData['z'])
        else:
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No planarity outliers observed.")

        if len(outStacd) > 0:
            n_outliers = len(outStacd)
            div = outlierFold.addDiv(style="font-size:110%;")
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

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No stacking distance outliers observed.")

        if len(outStaca) > 0:
            n_outliers = len(outStaca)
            div = outlierFold.addDiv(style="font-size:110%;")
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

            clearingDiv = outlierFold.addDiv()
            styleDiv = outlierFold.addDiv(style="color:navy;text-align:right;")
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
            div = outlierFold.addDiv(style="font-size:110%;")
            div.append("No stacking angle outliers observed.")


    def addScrollableDownloadableTable1(self, cycle_data, xmlnode=None, parent=None, internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        #create a "shell" div to contain the scrollable table and the hyperlink
        scrollableDownloadableTableDiv = parent.addDiv()
        #place a scrollable div into the shell: the table will be inserted into this div
        scrollableTableDiv = scrollableDownloadableTableDiv.addDiv()
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
        else:
            for title, key in _r_and_cc_columns(cycle_data):
                statisticNames.append(title)
                statisticInitial.append(cycle_data[key][0])
                statisticFinal.append(cycle_data[key][-1])
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
