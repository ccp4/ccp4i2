from __future__ import print_function

"""
Qt/matplotlib widget: made using matplotlib.backends.backend_qt4agg
"""
"""
     qtgui/MGQTmatplotlib.py: CCP4MG Molecular Graphics Program
     Copyright (C) 2011 University of York
     Copyright (C) 2012 Science and Technology Facilities Council

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

import sys
import os
import math
import glob
import functools
import traceback
import copy
import html
from collections.abc import Callable

from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, ScalarFormatter
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
import matplotlib
import numpy
from lxml import etree

from matplotlib.backends.backend_pdf import PdfPages


numpy.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore","'Legend' needs 'contains' method")


scaling_functions = { 'oneoversqrt' : lambda x,pos,format='%.1f': (format)%(math.pow(x,-0.5)) if x > 1e-8  else ' ' }
scaling_inverses = {scaling_functions['oneoversqrt'] : lambda x: 1/(x*x) if x*x > sys.float_info.min else sys.float_info.max }
scaling_numerical = {scaling_functions['oneoversqrt'] : lambda x,pos,format='%.1f': (format)%(math.pow(x,-0.5)) if x > sys.float_info.min else ("%e")%(sys.float_info.max) }

def getMatplotlibFontList():
    db = QtGui.QFontDatabase()

    warnings.filterwarnings("error","findfont: Font family")

    items = []
    for fname in db.families():
        try:
            fn = matplotlib.font_manager.findfont(str(fname))
            items.append(fname)
        except:
            pass
            #print "Font",fname," not found"

    fp = os.path.normpath(os.path.join(os.path.dirname(matplotlib.__file__),"mpl-data","fonts","ttf"))
    for root, dirs, files in os.walk(fp):
        for f in files:
            id = db.addApplicationFont(os.path.join(root,f))
            for family in db.applicationFontFamilies(id):
                items.append(family)
            # This is slow !?
            #db.removeApplicationFont(id)

    items = list(set(items))
    items.sort(key=lambda x: str(x).lower())

    warnings.filterwarnings("always","findfont: Font family")

    return items

def openFileToEtree(fileName=None,printout=True):
  from lxml import etree

  # Use this as etree.parse() seg faults on some Linux
  parser = etree.HTMLParser()
  f = open(fileName)
  s = f.read()
  f.close()
  tree = etree.fromstring(s, parser)
  return tree

class StackedWidget(QtWidgets.QStackedWidget):

    DetectPointClickEvent = QtCore.Signal(dict)
    DetectPointMoveEvent = QtCore.Signal(dict)
    RegionSelected = QtCore.Signal(dict)

    def __init__(self,parent=None):
        QtWidgets.QStackedWidget.__init__(self,parent)

class CCP4Table:
    def __str__(self):
        string = ''
        string = string + self.title +'\n'
        string = string + str(self.graphs_list) +'\n'
        string = string + str(self.headers) +'\n'
        string = string + str(self.data_lines) +'\n'
        return string

    def __init__(self,table):
        tname_start = table.find('$TABLE')+len('$TABLE')
        tname_start = table[tname_start:].find(':')+tname_start+1
        
        tname_end = table[tname_start:].find(':')
        if table[tname_start:].find('\n') < tname_end:
            tname_end= table[tname_start:].find('\n')
        self.title = table[tname_start:tname_end+tname_start].strip()

        rest = table[tname_end+tname_start:]

        self.graphs_list = []

        graphs_start = rest.find('$GRAPHS')+len('$GRAPHS')
        graphs_end = rest[graphs_start:].find('$$')
        graphs = rest[graphs_start:graphs_end+graphs_start].strip()
        colCount = 0
        currentHeader = ''
        
        for g in graphs:
            if g == ':':    
                colCount += 1
            if colCount == 4:
                s,name,n,cols = currentHeader.strip().split(':')
                theRange = None
                if n.startswith('N'):
                    theRange = ((None,None),(0,None))
                if not n.startswith('XD') and not n.startswith('N') and not n.startswith('A') and n.find('x')>0 and n.find('|')>0:
                    try:
                        xRange,yRange = n.split('x')
                        xRange = xRange.split('|')
                        yRange = yRange.split('|')
                        theRange = (float(xRange[0]),float(xRange[1])),(float(yRange[0]),float(yRange[1]))
                    except:
                        pass
                self.graphs_list.append((name.strip(), cols, theRange))
                colCount = 0
                currentHeader = ''
            else:
                currentHeader += g
                
                
        rest2 = rest[graphs_end+graphs_start:]

        headers_start = rest2.find('$$')+2
        headers_end = rest2[headers_start:].find('$$')
        self.headers = rest2[headers_start:headers_end+2].replace('\n',' ').split()

        rest3 = rest2[headers_end+2:]
        rest3 = rest3[rest3.find('$$')+2:]
        data_start = rest3.find('$$')+2
        rest3 = rest3[data_start:].lstrip()
        data_end = rest3.find('$$')
        data_start = 0
        data_lines_str =  rest3[data_start:data_end].strip().split('\n')
        self.data_lines = []
        for line_str in data_lines_str:
            vals = []
            vals_str = line_str.split()
            for val_str in vals_str:
                try:
                    vals.append(float(val_str))
                except:
                    vals.append('NA')
            self.data_lines.append(vals)
        self.x_scaling = False
        self.y_scaling = False
        self.custom_x_label = False
        # This is a bit iffy. We need to check if first column of data is one of these columns need this.
        if len(self.headers)>0 and (self.headers[0] == 'M(4SSQ/LL)' or self.headers[0] == '1/d^2' or self.headers[1] == '1/d^2' or self.headers[0] == '1/resol^2' or self.headers[1] == '1/resol^2'):
            self.custom_x_label = u'Resolution / \xc5'
            self.x_scaling = scaling_functions["oneoversqrt"]

    def toEtree(self):
        t = etree.Element("CCP4Table")
        t.attrib["title"] = self.title
        theData = etree.Element("data")
        theDataText = ''
        for tdl in self.data_lines:
            theDataText += ''.join([("%s "%n) for n in tdl[:]]).strip() + '\n'
        theData.text = theDataText
        t.append(theData)
        theHeaders = etree.Element("headers",separator=' ')
        theHeaders.text = ''.join([("%s "%n) for n in self.headers[:]]).strip() + '\n'
        t.append(theHeaders)
        for g in self.graphs_list:
            p = etree.Element('plot')
            theTitle = etree.Element("title")
            theTitle.text = g[0]
            p.append(theTitle)
            t.append(p)
            for s in g[1].split(',')[1:]:
                plotline = etree.Element('plotline', xcol=g[1].split(',')[0], ycol=s)
                p.append(plotline)
        return t

class LogGraphMainWindow(QtWidgets.QMainWindow):
    def __init__(self,parent=None):
        QtWidgets.QMainWindow.__init__(self,parent)
        self.loggraph = LogGraph()
        self.setCentralWidget(self.loggraph)
        menu = QtWidgets.QMenu("File",self)
        menu.addAction(self.loggraph.fileOpen)
        menu.addAction(self.loggraph.fileSave)
        menu.addAction(self.loggraph.fileSaveStatus)
        menu.addAction(self.loggraph.fileSaveAll)
        menu.addSeparator()
        exitAction = menu.addAction("Quit")
        exitAction.setShortcut(QtGui.QKeySequence.Quit)
        exitAction.triggered.connect(self.close)
        self.menuBar().addMenu(menu)
        editMenu = QtWidgets.QMenu("Edit",self)
        editMenu.addAction(self.loggraph.editLegendPos)
        editMenu.addAction(self.loggraph.editPlotStyle)
        editMenu.addAction(self.loggraph.prefAction)
        self.menuBar().addMenu(editMenu)
        self.setWindowTitle("Pimple")

        toolBar = self.addToolBar('')
        toolBar.addAction(self.loggraph.fileOpen)
        toolBar.addAction(self.loggraph.fileSave)
        toolBar.addAction(self.loggraph.fileSaveStatus)
        toolBar.addSeparator()
        toolBar.addAction(self.loggraph.editLegendPos)
        toolBar.addAction(self.loggraph.editPlotStyle)
        if self.loggraph.icons_missing:
            toolBar.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon)

class LogGraph(QtWidgets.QWidget):
    linestyles = ['-','--','-.',':','.']
    linestyles_alias = ['Solid','Dashed','Dash-dot','Dotted','Blank']
    markers = ['o','s','d','^','>','<','p','h','*']
    styles = markers + [
    r'$\lambda$',
    r'$\bowtie$',
    r'$\circlearrowleft$',
    r'$\clubsuit$',
    r'$\checkmark$']
    markers_alias = ['Circle','Square','Diamond','Up arrow','Right arrow','Left arrow','Pentagon','Hexagon','Star','Lambda','Bow tie','Circle arrow left','Clubs','Tick']

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    colors_alias = ['Blue', 'Green', 'Red', 'Cyan', 'Magenta', 'Yellow', 'Black']
    axis_styles_alias = ['Linear','Log','Symmetric Log']
    axis_styles = ['linear','log','symlog']

    currentDataChanged = QtCore.Signal(int)
    currentTableChanged = QtCore.Signal(int)

    def fctocolname(self,col):
        if abs(col[0]-1)<1e-5 and abs(col[1])<1e-5 and abs(col[2])<1e-5:
            return 'r'
        if abs(col[1]-1)<1e-5 and abs(col[0])<1e-5 and abs(col[2])<1e-5:
            return 'g'
        if abs(col[2]-1)<1e-5 and abs(col[0])<1e-5 and abs(col[1])<1e-5:
            return 'b'
        if abs(col[0]-0.75)<1e-5 and abs(col[2]-0.75)<1e-5 and abs(col[1])<1e-5:
            return 'm'
        if abs(col[1]-0.75)<1e-5 and abs(col[2]-0.75)<1e-5 and abs(col[0])<1e-5:
            return 'c'
        if abs(col[0]-0.75)<1e-5 and abs(col[1]-0.75)<1e-5 and abs(col[2])<1e-5:
            return 'y'
        if abs(col[0])<1e-5 and abs(col[1])<1e-5 and abs(col[2])<1e-5:
            return 'k'
        return 'b'

    @QtCore.Slot('QWidget','QWidget','QWidget')
    def applyPlotStyle(self,tab,plot,otherWidget):
        showLegend = otherWidget.showLegend.isChecked()
        if hasattr(tab.currentWidget(),"showPdf"):
            showPdf = tab.currentWidget().showPdf.isChecked()
            hist = self.histograms[tab.currentIndex()]['patches']
            colour = self.colors[tab.currentWidget().colourStyleCombo.currentIndex()]
            for p in hist:
                p.set_facecolor(colour)
            if len(hist)>1:
                if showPdf:
                    hist[0].show_pdf = True
                else:
                    hist[0].show_pdf = False
            doBestFit = showPdf
            normed = tab.currentWidget().normed.isChecked()
            if normed and not self.histograms[tab.currentIndex()]['normed']:
                label = hist[0].get_label()
                ccp4_uuid = hist[0].ccp4_uuid
                plot.canvas.fig.get_axes()[0].patches = []
                mu = self.histograms[tab.currentIndex()]['mu']
                sigma = self.histograms[tab.currentIndex()]['sigma']
                maxx = self.histograms[tab.currentIndex()]['maxx']
                vals = self.histograms[tab.currentIndex()]['vals']
                #n, bins, patches = plot.canvas.ax.hist(vals, int(maxx), normed=1, facecolor=hist[0].get_facecolor(), label=label)
                n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=1, facecolor=hist[0].get_facecolor(), label=label)
                self.histograms[tab.currentIndex()] = {'normed':True,'vals':vals,'maxx':maxx,'patches':patches,'bins':bins,'mu':mu, 'sigma':sigma}
                hist = self.histograms[tab.currentIndex()]['patches']
                hist[0].ccp4_uuid = ccp4_uuid
                plot.canvas.ax[0].set_ylim(0,1)
                plot.canvas.draw()
            if not normed and self.histograms[tab.currentIndex()]['normed']:
                label = hist[0].get_label()
                ccp4_uuid = hist[0].ccp4_uuid
                plot.canvas.fig.get_axes()[0].patches = []
                mu = self.histograms[tab.currentIndex()]['mu']
                sigma = self.histograms[tab.currentIndex()]['sigma']
                maxx = self.histograms[tab.currentIndex()]['maxx']
                vals = self.histograms[tab.currentIndex()]['vals']
                #n, bins, patches = plot.canvas.ax.hist(vals, int(maxx), normed=0, facecolor=hist[0].get_facecolor(), label=label)
                n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=1, facecolor=hist[0].get_facecolor(), label=label)
                self.histograms[tab.currentIndex()] = {'normed':False,'vals':vals,'maxx':maxx,'patches':patches,'bins':bins,'mu':mu, 'sigma':sigma}
                hist = self.histograms[tab.currentIndex()]['patches']
                hist[0].ccp4_uuid = ccp4_uuid
                plot.canvas.ax[0].set_ylim(0,1.1*numpy.max(n))
            lines = plot.canvas.fig.get_axes()[0].get_lines()
            lines = [l for l in lines if not hasattr(l,"ignore")]
            #print lines
            iline = 0
            l = self.histograms[tab.currentIndex()]
            if True:
                        if doBestFit:

                            bins = self.histograms[tab.currentIndex()]['bins']
                            mu = self.histograms[tab.currentIndex()]['mu']
                            sigma = self.histograms[tab.currentIndex()]['sigma']

                            bincenters = 0.5*(bins[1:]+bins[:-1])
                            y = matplotlib.mlab.normpdf( bincenters, mu, sigma)
                            haveBestFit = False
                            for l2 in lines:
                                if l2.get_label() == hist[0].ccp4_uuid+'_matplotlib_line_of_best_fit':
                                    haveBestFit = True
                                    l2.set_xdata(bincenters)
                                    l2.set_ydata(y)
                                    break
                            if not haveBestFit:
                                plot.canvas.plot(bincenters, y, 'r--', linewidth=1,label=hist[0].ccp4_uuid+'_matplotlib_line_of_best_fit')
                        else:
                            haveBestFit = False
                            for l2 in lines:
                                if l2.get_label() == hist[0].ccp4_uuid+'_matplotlib_line_of_best_fit':
                                    haveBestFit = True
                                    print("Need to get rid of", l2)
                                    l2.set_xdata([])
                                    l2.set_ydata([])
                                    break
        elif hasattr(tab.currentWidget(),"lineStyleCombo"):
            linestyle = self.linestyles[tab.currentWidget().lineStyleCombo.currentIndex()]
            linesize = tab.currentWidget().lineSizeSpin.value()
            if str(tab.currentWidget().markerStyleCombo.currentText()) == 'None':
                marker = ''
            else:
                marker = self.styles[tab.currentWidget().markerStyleCombo.currentIndex()]
            colour = self.colors[tab.currentWidget().colourStyleCombo.currentIndex()]
            markersize = tab.currentWidget().sizeSpin.value()
            markeredgewidth = tab.currentWidget().markeredgewidthSpin.value()
            label = str(tab.currentWidget().titleEdit.text())
            plot = self.graph.currentWidget().currentWidget()
            doBestFit = tab.currentWidget().bestFitCombo.currentIndex()
            thisInLegend = tab.currentWidget().showLegend.isChecked()
            isVisible = tab.currentWidget().showVisible.isChecked()
            excludeFromLegend = []
            def apply_to_lines(lines):
                iline = 0
                for l in lines:
                    if not l.get_label().endswith('_matplotlib_line_of_best_fit'):
                        if tab.currentIndex()-len(self.histograms) == iline:
                            if doBestFit:
                                xvals = l.get_xdata()
                                v = l.get_ydata()
                                if doBestFit == 4:
                                    from scipy import optimize
                                    e = lambda v,x,y: sum(((v[1]*numpy.exp(v[0]*x))-y)**2)

                                    v0=[0.3,2.]
                                    vmin = optimize.fmin(e, v0, args=(xvals,v), maxiter=10000, maxfun=10000)
                                    func = lambda x,v0,v1: v1*numpy.exp(v0*x)
                                    vec_func = numpy.vectorize(func)
                                    #print vmin
                                    yfit = vec_func(xvals,vmin[0],vmin[1])
                                else:
                                    m = numpy.polyfit(xvals,v,doBestFit)
                                    #print m
                                    yfit = numpy.polyval(m,xvals)
                                haveBestFit = False
                                for l2 in lines:
                                    if hasattr(l,"ccp4_uuid") and l2.get_label() == l.ccp4_uuid+'_matplotlib_line_of_best_fit':
                                        haveBestFit = True
                                        l2.set_xdata(xvals)
                                        l2.set_ydata(yfit)
                                        break
                                if not haveBestFit and hasattr(l,"ccp4_uuid"):
                                    plot.canvas.plot(xvals,yfit,'k',label=l.ccp4_uuid+'_matplotlib_line_of_best_fit')
                            else:
                                haveBestFit = False
                                for l2 in lines:
                                    if hasattr(l,"ccp4_uuid") and l2.get_label() == l.ccp4_uuid+'_matplotlib_line_of_best_fit':
                                        haveBestFit = True
                                        #print "Need to get rid of", l2
                                        l2.set_xdata([])
                                        l2.set_ydata([])
                                        break
                            l.set_label(label)
                            l.set_marker(marker)
                            l.set_markersize(markersize)
                            l.set_linewidth(linesize)
                            l.set_linestyle(linestyle)
                            l.set_color(colour)
                            if not thisInLegend:
                                l.hide_from_legend = True
                            else:
                                l.hide_from_legend = False
                            if isVisible:
                                l.set_visible(True)
                            else:
                                l.set_visible(False)
                            l.set_markeredgewidth(markeredgewidth)
                            break
                    iline = iline + 1

            for axes in plot.canvas.fig.get_axes():
                lines = axes.get_lines()
                lines = [l for l in lines if not hasattr(l,"ignore")]
                apply_to_lines(lines)

        plot.canvas.ax[0].set_xscale(self.axis_styles[otherWidget.x_axis_style.currentIndex()])
        plot.canvas.ax[0].set_yscale(self.axis_styles[otherWidget.y_axis_style.currentIndex()])

        if sys.version_info >= (3,0):
            plot.canvas.ax[0].set_xlabel(str(otherWidget.x_axis_label.text()))
            plot.canvas.ax[0].set_ylabel(str(otherWidget.y_axis_label.text()))
            plot.canvas.ax[0].set_title(str(otherWidget.plot_title.text()))
            plot.canvas.xlabel = str(otherWidget.x_axis_label.text())
            plot.canvas.ylabel = str(otherWidget.y_axis_label.text())
            plot.canvas.title = str(otherWidget.plot_title.text())
        else:
            plot.canvas.ax[0].set_xlabel(unicode(otherWidget.x_axis_label.text()))
            plot.canvas.ax[0].set_ylabel(unicode(otherWidget.y_axis_label.text()))
            plot.canvas.ax[0].set_title(unicode(otherWidget.plot_title.text()))
            plot.canvas.xlabel = unicode(otherWidget.x_axis_label.text())
            plot.canvas.ylabel = unicode(otherWidget.y_axis_label.text())
            plot.canvas.title = unicode(otherWidget.plot_title.text())

        xminold = plot.canvas.ax[0].get_xlim()[0]
        xmaxold = plot.canvas.ax[-1].get_xlim()[1]
        yminold = plot.canvas.ax[-1].get_ylim()[0]
        ymaxold = plot.canvas.ax[0].get_ylim()[1]
        xminstr, xmaxstr = str(otherWidget.x_axis_range.text()).split(',')
        xmin = float(xminstr)
        xmax = float(xmaxstr)

        if plot.canvas.x_scaling and plot.canvas.x_scaling in scaling_inverses:
            xmin = scaling_inverses[plot.canvas.x_scaling](xmin)
            xmax = scaling_inverses[plot.canvas.x_scaling](xmax)

        if abs(xmin-xminold)>1e-7 or abs(xmax-xmaxold)>1e-7:
            plot.canvas.ax[0].set_xlim(xmin,plot.canvas.ax[0].get_xlim()[1])
            plot.canvas.ax[-1].set_xlim(plot.canvas.ax[-1].get_xlim()[0],xmax)
            plot.canvas.custom_xlim = xmin,xmax

        yminstr, ymaxstr = str(otherWidget.y_axis_range.text()).split(',')
        ymin = float(yminstr)
        ymax = float(ymaxstr)

        if plot.canvas.y_scaling and plot.canvas.y_scaling in scaling_inverses:
            ymin = scaling_inverses[plot.canvas.x_scaling](ymin)
            ymax = scaling_inverses[plot.canvas.x_scaling](ymax)

        if abs(ymin-yminold)>1e-7 or abs(ymax-ymaxold)>1e-7:
            plot.canvas.ax[0].set_ylim(plot.canvas.ax[0].get_ylim()[0],ymax)
            plot.canvas.ax[-1].set_ylim(ymin,plot.canvas.ax[-1].get_ylim()[1])
            plot.canvas.custom_ylim = ymin,ymax

        try:
            xbreak = str(otherWidget.x_axis_break.text()).split(')')
            if len(xbreak)>0:
                new_xbreak = []
                for br in xbreak:
                    theBr = br.strip().lstrip('('). split(',')
                    if len(theBr)==2:
                        lo = float(theBr[0].strip())
                        hi = float(theBr[1].strip())
                        if plot.canvas.x_scaling and plot.canvas.x_scaling in scaling_inverses:
                            lo = scaling_inverses[plot.canvas.x_scaling](lo)
                            hi = scaling_inverses[plot.canvas.x_scaling](hi)
                        new_xbreak.append((lo,hi))
                new_xbreak = tuple(new_xbreak)
                if len(new_xbreak) == 0:
                    new_xbreak = None
            else:
                new_xbreak = None
        except:
            new_xbreak = None
            print("Badly formatted x-axis breaks declaration.")
            print("Should be of form: (xlow,xhigh), (xlow,xhigh), (xlow,xhigh), ...")
            QtWidgets.QMessageBox.critical(self,"Bad x-axis break","Bad x-axis break declaration.<br>Should be of form: (xlow,xhigh), (xlow,xhigh), (xlow,xhigh), ...")

        try:
            ybreak = str(otherWidget.y_axis_break.text()).split(')')
            if len(ybreak)>0:
                new_ybreak = []
                for br in ybreak:
                    theBr = br.strip().lstrip('('). split(',')
                    if len(theBr)==2:
                        lo = float(theBr[0].strip())
                        hi = float(theBr[1].strip())
                        if plot.canvas.y_scaling and plot.canvas.y_scaling in scaling_inverses:
                            lo = scaling_inverses[plot.canvas.y_scaling](lo)
                            hi = scaling_inverses[plot.canvas.y_scaling](hi)
                        new_ybreak.append((lo,hi))
                new_ybreak = tuple(new_ybreak)
                if len(new_ybreak) == 0:
                    new_ybreak = None
            else:
                new_ybreak = None
        except:
            new_ybreak = None
            print("Badly formatted y-axis breaks declaration.")
            print("Should be of form: (ylow,yhigh), (ylow,yhigh), (ylow,yhigh), ...")
            QtWidgets.QMessageBox.critical(self,"Bad y-axis break","Bad y-axis break declaration.<br>Should be of form: (ylow,yhigh), (ylow,yhigh), (ylow,yhigh), ...")

        plot.canvas.xbreak = new_xbreak
        plot.canvas.ybreak = new_ybreak
        plot.canvas.set_breaks(plot.canvas.xbreak,plot.canvas.ybreak)

        if plot.canvas.xbreak and abs(xmin-xminold)<1e-7 or abs(xmax-xmaxold)<1e-7:
            plot.canvas.ax[0].set_xlim(xmin,xmax)
            plot.canvas.ax[-1].set_xlim(xmin,xmax)

        if plot.canvas.ybreak and abs(ymin-yminold)<1e-7 or abs(ymax-ymaxold)<1e-7:
            plot.canvas.ax[0].set_ylim(ymin,ymax)
            plot.canvas.ax[-1].set_ylim(ymin,ymax)

        plot.canvas.fix_aspect_ratio = otherWidget.fixedAspectRatio.isChecked()

        if showLegend:
            plot.canvas.legend()
        else:
            if plot.canvas.ax[0].get_legend():
                plot.canvas.ax[0].legend_.set_visible(False)

        plot.canvas.plot()
        plot.canvas.draw()
        plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))
        plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))

    @QtCore.Slot()
    def editPreferences(self):
        if not hasattr(LogGraph,"matplotlibFontList"):
            LogGraph.matplotlibFontList = getMatplotlibFontList()
        editWidget = QtWidgets.QDialog(self)
        layout = QtWidgets.QGridLayout()
        label = QtWidgets.QLabel(self.tr("Title font"))
        titleFontSel = QtWidgets.QComboBox()
        titleFontSel.addItems(LogGraph.matplotlibFontList)
        titleFontSel.setEditable(True)
        layout.addWidget(label,layout.rowCount(),0)
        titleSize = QtWidgets.QSpinBox()
        titleSize.setRange(8,100)
        layout.addWidget(titleSize,layout.rowCount()-1,1)
        layout.addWidget(titleFontSel,layout.rowCount()-1,2)
        titleItalic = QtWidgets.QCheckBox(self.tr("Italic"))
        layout.addWidget(titleItalic,layout.rowCount()-1,3)
        titleBold = QtWidgets.QCheckBox(self.tr("Bold"))
        layout.addWidget(titleBold,layout.rowCount()-1,4)
        label = QtWidgets.QLabel(self.tr("Axes ticker font"))
        axesFontSel = QtWidgets.QComboBox()
        axesFontSel.addItems(LogGraph.matplotlibFontList)
        axesFontSel.setEditable(True)
        layout.addWidget(label,layout.rowCount(),0)
        axesSize = QtWidgets.QSpinBox()
        axesSize.setRange(8,100)
        layout.addWidget(axesSize,layout.rowCount()-1,1)
        layout.addWidget(axesFontSel,layout.rowCount()-1,2)
        axesItalic = QtWidgets.QCheckBox(self.tr("Italic"))
        layout.addWidget(axesItalic,layout.rowCount()-1,3)
        axesBold = QtWidgets.QCheckBox(self.tr("Bold"))
        layout.addWidget(axesBold,layout.rowCount()-1,4)
        label = QtWidgets.QLabel(self.tr("Axes label font"))
        axesLabelFontSel = QtWidgets.QComboBox()
        axesLabelFontSel.addItems(LogGraph.matplotlibFontList)
        axesLabelFontSel.setEditable(True)
        layout.addWidget(label,layout.rowCount(),0)
        axesLabelSize = QtWidgets.QSpinBox()
        axesLabelSize.setRange(8,100)
        layout.addWidget(axesLabelSize,layout.rowCount()-1,1)
        layout.addWidget(axesLabelFontSel,layout.rowCount()-1,2)
        axesLabelItalic = QtWidgets.QCheckBox(self.tr("Italic"))
        layout.addWidget(axesLabelItalic,layout.rowCount()-1,3)
        axesLabelBold = QtWidgets.QCheckBox(self.tr("Bold"))
        layout.addWidget(axesLabelBold,layout.rowCount()-1,4)
        label = QtWidgets.QLabel(self.tr("Legend font"))
        legendFontSel = QtWidgets.QComboBox()
        legendFontSel.addItems(LogGraph.matplotlibFontList)
        legendFontSel.setEditable(True)
        layout.addWidget(label,layout.rowCount(),0)
        legendSize = QtWidgets.QSpinBox()
        legendSize.setRange(8,100)
        layout.addWidget(legendSize,layout.rowCount()-1,1)
        layout.addWidget(legendFontSel,layout.rowCount()-1,2)
        legendItalic = QtWidgets.QCheckBox(self.tr("Italic"))
        layout.addWidget(legendItalic,layout.rowCount()-1,3)
        legendBold = QtWidgets.QCheckBox(self.tr("Bold"))
        layout.addWidget(legendBold,layout.rowCount()-1,4)
        dialogButtons = QtWidgets.QDialogButtonBox()
        applyButton = dialogButtons.addButton(QtWidgets.QDialogButtonBox.Apply)
        closeButton = dialogButtons.addButton(QtWidgets.QDialogButtonBox.Close)
        closeButton.clicked.connect(editWidget.close)
        layout.addWidget(dialogButtons,layout.rowCount(),0,1,layout.columnCount())

        legendFontSelParams = {'family':legendFontSel,'italic':legendItalic,'weight':legendBold,'size':legendSize}
        titleFontSelParams = {'family':titleFontSel,'italic':titleItalic,'weight':titleBold,'size':titleSize}
        axesFontSelParams = {'family':axesFontSel,'italic':axesItalic,'weight':axesBold,'size':axesSize}
        axesLabelFontSelParams = {'family':axesLabelFontSel,'italic':axesLabelItalic,'weight':axesLabelBold,'size':axesLabelSize}

        for fn,selector in ((self.legendFontSel,legendFontSelParams),(self.titleFontSel,titleFontSelParams),(self.axesFontSel,axesFontSelParams),(self.axesLabelFontSel,axesLabelFontSelParams)):
                    isItalic = False
                    if fn and "slant" in fn and fn['slant'] == "italic":
                        selector['italic'].setChecked(True)
                    else:
                        selector['italic'].setChecked(False)
                    isBold = False
                    if fn and "weight" in fn and fn['weight'] == "bold":
                        selector['weight'].setChecked(True)
                    else:
                        selector['weight'].setChecked(False)
                    if fn and "family" in fn:
                        if sys.version_info >= (3,0):
                            selector['family'].setCurrentIndex(selector['family'].findText(fn['family']))
                        else:
                            selector['family'].setCurrentIndex(selector['family'].findText(unicode(fn['family'],"utf-8")))
                    else:
                        if sys.version_info >= (3,0):
                            selector['family'].setCurrentIndex(selector['family'].findText('Bitstream Vera Sans'))
                        else:
                            selector['family'].setCurrentIndex(selector['family'].findText(unicode('Bitstream Vera Sans',"utf-8")))
                    if fn and "family" in fn:
                        selector['size'].setValue(fn['size'])
                    else:
                        selector['size'].setValue(10)

        applyButton.clicked.connect(functools.partial(self.applyPreferencesFromWidget,titleFontSelParams,axesFontSelParams,axesLabelFontSelParams,legendFontSelParams))
        editWidget.setLayout(layout)
        
        editWidget.exec_()

    def widgetPrefsToPrefs(self,params_in):
        params = {}
        if 'family' in params_in:
            params['family'] = str(params_in['family'].currentText())
        if 'size' in params_in:
            params['size'] = params_in['size'].value()
        if 'weight' in params_in and params_in['weight'].isChecked():
            params['weight'] = 'bold'
        else:
            params['weight'] = 'normal'
        if 'italic' in params_in and params_in['italic'].isChecked():
            params['slant'] = 'italic'
        else:
            params['slant'] = 'normal'
        return params

    @QtCore.Slot(dict,dict,dict,dict)
    def applyPreferencesFromWidget(self,titleFontSel_in,axesFontSel_in,axesLabelFontSel_in,legendFontSel_in):
        titleFontSel = self.widgetPrefsToPrefs(titleFontSel_in)
        axesFontSel = self.widgetPrefsToPrefs(axesFontSel_in)
        axesLabelFontSel = self.widgetPrefsToPrefs(axesLabelFontSel_in)
        legendFontSel = self.widgetPrefsToPrefs(legendFontSel_in)
        self.applyPreferences(titleFontSel,axesFontSel,axesLabelFontSel,legendFontSel)

        
    def fontsToEtree(self):
        tree = etree.Element('Fonts')
        tree.append(etree.Element('titleFont'))
        tree.append(etree.Element('axesTickerFont'))
        tree.append(etree.Element('axesLabelFont'))
        tree.append(etree.Element('legendFont'))

        for p2, sel in zip(['titleFont','axesTickerFont','axesLabelFont','legendFont'],[self.titleFontSel,self.axesFontSel,self.axesLabelFontSel,self.legendFontSel]):
            for p in tree.xpath(p2):
                for k in ["family","size","weight","slant"]:
                    if sel and k in sel:
                        if type(sel[k]) == int:
                            p.set(k,str(sel[k]))
                        else:
                            p.set(k,sel[k])
                    else:
                        if k in p.attrib:
                            p.attrib.pop(k)
        return tree

    def applyPreferences(self,titleFontSel=None,axesFontSel=None,axesLabelFontSel=None,legendFontSel=None):
        self.titleFontSel = titleFontSel
        self.axesFontSel = axesFontSel
        self.axesLabelFontSel = axesLabelFontSel
        self.legendFontSel = legendFontSel

        currenti = self.graph.currentIndex()
        if not self.graph.currentWidget():
            return
        currentj = self.graph.currentWidget().currentIndex()
        for i in range(self.graph.count()):
            self.graph.setCurrentIndex(i)
            sw = self.graph.widget(i)
            for j in range(sw.count()):
                sw.setCurrentIndex(j)
                aplot = sw.widget(j)
                if not hasattr(aplot,"canvas"):
                    continue
                aplot.canvas.applyPreferences(titleFontSel,axesFontSel,axesLabelFontSel,legendFontSel)

        self.graph.setCurrentIndex(currenti)
        self.graph.currentWidget().setCurrentIndex(currentj)

    @QtCore.Slot()
    def editPlotPlotStyle(self):
        plot = self.graph.currentWidget().currentWidget()
        lines = plot.canvas.fig.get_axes()[0].get_lines()
        lines = [l for l in lines if not hasattr(l,"ignore")]
        editWidget = QtWidgets.QDialog(self)
        layout = QtWidgets.QVBoxLayout()
        tab = QtWidgets.QTabWidget()
        for hist in self.histograms:
            if len(hist)>1:
                if not hasattr(hist['patches'][0],"ccp4_uuid"):
                    if sys.platform != 'win32':
                        import uuid
                        uuid._uuid_generate_time = None
                        uuid._uuid_generate_random = None
                        if sys.version_info >= (3,0):
                            uuid_str = uuid.uuid4().hex
                        else:
                            uuid_str = uuid.uuid4().get_hex()
                    else:
                        import msilib
                        uuid_str = msilib.gen_uuid().strip('{').strip('}').replace('-','')
                    hist['patches'][0].ccp4_uuid = uuid_str
                widget = QtWidgets.QWidget()
                widgetLayout = QtWidgets.QGridLayout()
                widget.setLayout(widgetLayout)
                label = QtWidgets.QLabel(self.tr("Colour"))
                widget.colourStyleCombo = QtWidgets.QComboBox()
                widget.colourStyleCombo.addItems(self.colors_alias)
                widget.colourStyleCombo.setCurrentIndex(self.colors.index(self.fctocolname(hist['patches'][0].get_facecolor())))
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.colourStyleCombo,widgetLayout.rowCount()-1,1)
                widget.showPdf = QtWidgets.QCheckBox()
                pdfLabel = QtWidgets.QLabel(self.tr("Show prob. distn. function"))
                widgetLayout.addWidget(pdfLabel,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.showPdf,widgetLayout.rowCount()-1,1)
                if hasattr(hist['patches'][0],"show_pdf") and hist['patches'][0].show_pdf:
                    widget.showPdf.setChecked(True)
                else:
                    widget.showPdf.setChecked(False)
                widget.normed = QtWidgets.QCheckBox()
                normLabel = QtWidgets.QLabel(self.tr("Normalized"))
                widgetLayout.addWidget(normLabel,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.normed,widgetLayout.rowCount()-1,1)
                widget.normed.setChecked(hist['normed'])
                """
                # Do not know how to do this for histograms.
                widget.showLegend = QtWidgets.QCheckBox()
                legendLabel = QtWidgets.QLabel(self.tr("Display in legend"))
                widgetLayout.addWidget(legendLabel,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.showLegend,widgetLayout.rowCount()-1,1)
                print hasattr(hist[0],"hide_from_legend")
                if hasattr(hist[0],"hide_from_legend") and hist[0].hide_from_legend:
                    widget.showLegend.setChecked(False)
                else:
                    widget.showLegend.setChecked(True)
                """
                widgetLayout.setRowStretch(widgetLayout.rowCount(),10)
                tab.addTab(widget,hist['patches'][0].get_label())
        for l in lines:
            if l.get_label() is not None and not l.get_label().endswith('_matplotlib_line_of_best_fit'):
                if not hasattr(l,"ccp4_uuid"):
                    if sys.platform != 'win32':
                        import uuid
                        if sys.version_info >= (3,0):
                            uuid_str = uuid.uuid4().hex
                        else:
                            uuid_str = uuid.uuid4().get_hex()
                    else:
                        import msilib
                        uuid_str = msilib.gen_uuid().strip('{').strip('}').replace('-','')
                    l.ccp4_uuid = uuid_str
                widget = QtWidgets.QWidget()
                widgetLayout = QtWidgets.QGridLayout()
                widget.setLayout(widgetLayout)
                label = QtWidgets.QLabel(self.tr("Line style"))
                widget.lineStyleCombo = QtWidgets.QComboBox()
                widget.lineStyleCombo.addItems(self.linestyles_alias)
                widget.lineStyleCombo.setCurrentIndex(self.linestyles.index(l.get_linestyle()))
                widgetLayout.addWidget(label,0,0)
                widgetLayout.addWidget(widget.lineStyleCombo,0,1)
                label = QtWidgets.QLabel(self.tr("Marker"))
                widget.markerStyleCombo = QtWidgets.QComboBox()
                widget.markerStyleCombo.addItems(self.markers_alias)
                widget.markerStyleCombo.addItem('None')
                marker = l.get_marker()
                if marker and marker != 'None':
                    widget.markerStyleCombo.setCurrentIndex(self.styles.index(l.get_marker()))
                else:
                    widget.markerStyleCombo.setCurrentIndex(widget.markerStyleCombo.count()-1)
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.markerStyleCombo,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Colour"))
                widget.colourStyleCombo = QtWidgets.QComboBox()
                widget.colourStyleCombo.addItems(self.colors_alias)
                if not l.get_color() in self.colors:
                    try:
                        c1 = matplotlib.colors.ColorConverter().to_rgb(l.get_color())
                        for cn in self.colors:
                            c2 = matplotlib.colors.ColorConverter().to_rgb(cn)
                            if c2 == c1:
                                idx = self.colors.index(cn)
                                break
                        else:
                            self.colors.append(l.get_color())
                            self.colors_alias.append(l.get_color().capitalize())
                            widget.colourStyleCombo.addItem(l.get_color())
                            idx = self.colors.index(l.get_color())
                    except:
                        exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                        sys.stderr.write(str(exc_type)+'\n')
                        sys.stderr.write(str(exc_value)+'\n')
                        idx = 0
                else:
                    idx = self.colors.index(l.get_color())
                for c in self.colors:
                    mplc = matplotlib.colors.ColorConverter().to_rgb(c)
                    p = QtGui.QPixmap(20,12)
                    col = QtGui.QColor()
                    col.setRgbF(mplc[0],mplc[1],mplc[2])
                    p.fill(col)
                    widget.colourStyleCombo.setItemIcon(self.colors.index(c),QtGui.QIcon(p))
                widget.colourStyleCombo.setCurrentIndex(idx)
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.colourStyleCombo,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Line width"))
                widget.lineSizeSpin = QtWidgets.QDoubleSpinBox()
                widget.lineSizeSpin.setValue(float(l.get_linewidth()))
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.lineSizeSpin,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Marker size"))
                widget.sizeSpin = QtWidgets.QSpinBox()
                widget.sizeSpin.setValue(l.get_markersize())
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.sizeSpin,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Marker edge width"))
                widget.markeredgewidthSpin = QtWidgets.QDoubleSpinBox()
                widget.markeredgewidthSpin.setValue(l.get_markeredgewidth())
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.markeredgewidthSpin,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Title"))
                widget.titleEdit = QtWidgets.QLineEdit()
                widget.titleEdit.setText(l.get_label())
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.titleEdit,widgetLayout.rowCount()-1,1)
                label = QtWidgets.QLabel(self.tr("Line of best fit"))
                widget.bestFitCombo = QtWidgets.QComboBox()
                fitItems = ['None','Linear','Quadratic','Cubic']
                try:
                    from scipy import optimize
                    fitItems.append('Exponential')
                except:
                    pass
                widget.bestFitCombo.addItems(fitItems)
                widgetLayout.addWidget(label,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.bestFitCombo,widgetLayout.rowCount()-1,1)
                widget.showVisible = QtWidgets.QCheckBox()
                visibleLabel = QtWidgets.QLabel(self.tr("Visible on graph"))
                widgetLayout.addWidget(visibleLabel,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.showVisible,widgetLayout.rowCount()-1,1)
                widget.showLegend = QtWidgets.QCheckBox()
                legendLabel = QtWidgets.QLabel(self.tr("Display in legend"))
                widgetLayout.addWidget(legendLabel,widgetLayout.rowCount(),0)
                widgetLayout.addWidget(widget.showLegend,widgetLayout.rowCount()-1,1)
                if hasattr(l,"hide_from_legend") and l.hide_from_legend:
                    widget.showLegend.setChecked(False)
                else:
                    widget.showLegend.setChecked(True)
                if hasattr(l,"get_visible") and l.get_visible():
                    widget.showVisible.setChecked(True)
                else:
                    widget.showVisible.setChecked(False)
                tab.addTab(widget,l.get_label())
        layout.addWidget(tab)

        otherWidget = QtWidgets.QWidget()
        otherLayout = QtWidgets.QGridLayout()
        otherWidget.setLayout(otherLayout)
        showLegendLabel = QtWidgets.QLabel(self.tr("Show legend"))
        otherWidget.showLegend = QtWidgets.QCheckBox()
        otherWidget.x_axis_style = QtWidgets.QComboBox()
        otherWidget.y_axis_style = QtWidgets.QComboBox()
        otherWidget.x_axis_style.addItems(self.axis_styles_alias)
        otherWidget.y_axis_style.addItems(self.axis_styles_alias)
        xaxisLabel = QtWidgets.QLabel(self.tr("x-axis type"))
        yaxisLabel = QtWidgets.QLabel(self.tr("y-axis type"))
        otherLayout.addWidget(showLegendLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.showLegend,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(xaxisLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.x_axis_style,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(yaxisLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.y_axis_style,otherLayout.rowCount()-1,1)
        otherWidget.x_axis_style.setCurrentIndex(self.axis_styles.index(plot.canvas.ax[0].get_xscale()))
        otherWidget.y_axis_style.setCurrentIndex(self.axis_styles.index(plot.canvas.ax[0].get_yscale()))

        xaxisLabelLabel = QtWidgets.QLabel(self.tr("x-axis label"))
        yaxisLabelLabel = QtWidgets.QLabel(self.tr("y-axis label"))
        otherWidget.x_axis_label = QtWidgets.QLineEdit()
        otherWidget.x_axis_label.setText(plot.canvas.xlabel)
        otherWidget.y_axis_label = QtWidgets.QLineEdit()
        otherWidget.y_axis_label.setText(plot.canvas.ylabel)
        plotTitleLabel = QtWidgets.QLabel(self.tr("Plot title"))
        otherWidget.plot_title = QtWidgets.QLineEdit()
        otherWidget.plot_title.setText(plot.canvas.title)
        
        xRangeLabel = QtWidgets.QLabel(self.tr("x-axis range"))
        otherWidget.x_axis_range = QtWidgets.QLineEdit()
        xmin = plot.canvas.ax[0].get_xlim()[0]
        xmax = plot.canvas.ax[-1].get_xlim()[1]
        if plot.canvas.x_scaling and plot.canvas.x_scaling in scaling_numerical:
            otherWidget.x_axis_range.setText("%s, %s" % (scaling_numerical[plot.canvas.x_scaling](xmin,0,"%.5f"),scaling_numerical[plot.canvas.x_scaling](xmax,0,"%.5f")))
        else:
            otherWidget.x_axis_range.setText("%.3f, %.3f" % (xmin,xmax))
        yRangeLabel = QtWidgets.QLabel(self.tr("y-axis range"))
        otherWidget.y_axis_range = QtWidgets.QLineEdit()
        ymax = plot.canvas.ax[0].get_ylim()[1]
        ymin = plot.canvas.ax[-1].get_ylim()[0]
        if plot.canvas.y_scaling and plot.canvas.y_scaling in scaling_numerical:
            otherWidget.y_axis_range.setText("%s, %s" % (scaling_numerical[plot.canvas.y_scaling](ymin,0,"%.5f"),scaling_numerical[plot.canvas.y_scaling](ymax,0,"%.5f")))
        else:
            otherWidget.y_axis_range.setText("%.3f, %.3f" % (ymin,ymax))

        xBreakLabel = QtWidgets.QLabel(self.tr("x-axis breaks"))
        otherWidget.x_axis_break = QtWidgets.QLineEdit()
        xBreakText = ''
        if plot.canvas.xbreak:
            for b in plot.canvas.xbreak:
                if plot.canvas.x_scaling and plot.canvas.x_scaling in scaling_numerical:
                    xBreakText += ("(%s, %s) " % (scaling_numerical[plot.canvas.x_scaling](b[0],0,"%.5f"),scaling_numerical[plot.canvas.x_scaling](b[1],0,"%.5f")))
                else:
                    xBreakText += ("(%.3f, %.3f) " % (b[0],b[1]))
        otherWidget.x_axis_break.setText(xBreakText)
        yBreakLabel = QtWidgets.QLabel(self.tr("y-axis breaks"))
        otherWidget.y_axis_break = QtWidgets.QLineEdit()
        yBreakText = ''
        if plot.canvas.ybreak:
            for b in plot.canvas.ybreak:
                if plot.canvas.y_scaling and plot.canvas.y_scaling in scaling_numerical:
                    yBreakText += ("(%s, %s) " % (scaling_numerical[plot.canvas.y_scaling](b[0],0,"%.5f"),scaling_numerical[plot.canvas.y_scaling](b[1],0,"%.5f")))
                else:
                    yBreakText += ("(%.3f, %.3f) " % (b[0],b[1]))
        otherWidget.y_axis_break.setText(yBreakText)

        otherLayout.addWidget(xaxisLabelLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.x_axis_label,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(yaxisLabelLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.y_axis_label,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(plotTitleLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.plot_title,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(xRangeLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.x_axis_range,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(yRangeLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.y_axis_range,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(xBreakLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.x_axis_break,otherLayout.rowCount()-1,1)
        otherLayout.addWidget(yBreakLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.y_axis_break,otherLayout.rowCount()-1,1)

        fixedAspectRatioLabel = QtWidgets.QLabel(self.tr("1:1 aspect"))
        otherWidget.fixedAspectRatio = QtWidgets.QCheckBox()
        if plot.canvas.fix_aspect_ratio:
            otherWidget.fixedAspectRatio.setChecked(True)
        else:
            otherWidget.fixedAspectRatio.setChecked(False)
        otherLayout.addWidget(fixedAspectRatioLabel,otherLayout.rowCount(),0)
        otherLayout.addWidget(otherWidget.fixedAspectRatio,otherLayout.rowCount()-1,1)

        if plot.canvas.ax[0].get_legend():
            otherWidget.showLegend.setChecked(True)
        else:
            otherWidget.showLegend.setChecked(False)
        layout.addWidget(otherWidget)

        dialogButtons = QtWidgets.QDialogButtonBox()
        applyButton = dialogButtons.addButton(QtWidgets.QDialogButtonBox.Apply)
        closeButton = dialogButtons.addButton(QtWidgets.QDialogButtonBox.Close)
        closeButton.clicked.connect(editWidget.close)
        layout.addWidget(dialogButtons)
        editWidget.setLayout(layout)
        applyButton.clicked.connect(functools.partial(self.applyPlotStyle,tab,plot,otherWidget))
        editWidget.exec_()

    @QtCore.Slot()
    def editLegendPosition(self):
        plot = self.graph.currentWidget().currentWidget()
        plot.canvas.movingLegend = True
        plot.canvas.setFocus(QtCore.Qt.OtherFocusReason)
        
    def loadFile(self,filename='',hist=False,xsd=None,select=None):
        oldCount = self.table_combo.count()
        if filename.endswith('.log') or filename.endswith('.txt'):
            f = open(filename)
            b = f.read()
            f.close()
            splits = b[b.find("$TABLE"):].split("$TABLE")
            newsplits = []
            gs = []
            #self.table_combo.blockSignals(True)
            self.data_combo.blockSignals(True)
            for ns in splits[1:]:
                if select is None or ns in select:
                    ns = "$TABLE"+ns
                    newsplits.append(ns)
                    table = CCP4Table(ns)
                    gs.append(self.addCCP4Table(table))
                    #print self.table_combo.count(), oldCount
            if self.table_combo.count()>0:
                self.table_combo.setCurrentIndex(oldCount)
                self.data_combo.setCurrentIndex(0)
                self.graph.setCurrentIndex(oldCount)
                self.graph.currentWidget().setCurrentIndex(0)
                self.setCurrentData(0)
            self.table_combo.blockSignals(False)
            self.data_combo.blockSignals(False)
            return gs
        elif filename.endswith('.xml'):
            self.data_combo.blockSignals(True)
            graphs = self.addXMLCCP4TableFile(filename,xsd=xsd)
            if self.table_combo.count()>0:
                self.table_combo.setCurrentIndex(oldCount)
                self.data_combo.setCurrentIndex(0)
                self.graph.setCurrentIndex(oldCount)
                self.graph.currentWidget().setCurrentIndex(0)
                self.setCurrentData(0)
            self.table_combo.blockSignals(False)
            self.data_combo.blockSignals(False)
            return graphs
        elif filename.endswith('.html'):
            graphs = self.addCCP4ReportFile(filename,select=select)
            if self.table_combo.count()>0:
                self.table_combo.setCurrentIndex(oldCount)
                self.data_combo.setCurrentIndex(0)
                self.graph.setCurrentIndex(oldCount)
                self.graph.currentWidget().setCurrentIndex(0)
                self.setCurrentData(0)
            return graphs
        elif filename.endswith('.surf'):
            try:
                ret = self.addSurfaceFile(filename)
                if self.table_combo.count()>0:
                    self.table_combo.setCurrentIndex(oldCount)
                    self.data_combo.setCurrentIndex(0)
                    self.graph.setCurrentIndex(oldCount)
                    self.graph.currentWidget().setCurrentIndex(0)
                return [ret]
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type)+'\n')
                sys.stderr.write(str(exc_value)+'\n')
                traceback.print_tb(exc_tb)
                QtWidgets.QMessageBox.critical(self,"File open failed","Failed to open file "+filename+ \
                               '\n'+str(exc_type)+'\n'+str(exc_value)+'\n')
        elif filename.endswith('.agr') or filename.endswith('.xmgr'):
            try:
                ret = self.addGraceFile(filename)
                if self.table_combo.count()>0:
                    self.table_combo.setCurrentIndex(oldCount)
                    self.data_combo.setCurrentIndex(0)
                    self.graph.setCurrentIndex(oldCount)
                    self.graph.currentWidget().setCurrentIndex(0)
                return [ret]
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type)+'\n')
                sys.stderr.write(str(exc_value)+'\n')
                traceback.print_tb(exc_tb)
                QtWidgets.QMessageBox.critical(self,"File open failed","Failed to open file "+filename+\
                              '\n'+str(exc_type)+'\n'+str(exc_value)+'\n' )
        else:
            try:
                ret = self.addPlainFile(filename,hist=hist)
                if self.table_combo.count()>0:
                    self.table_combo.setCurrentIndex(oldCount)
                    self.data_combo.setCurrentIndex(0)
                    self.graph.setCurrentIndex(oldCount)
                    self.graph.currentWidget().setCurrentIndex(0)
                return [ret]
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type)+'\n')
                sys.stderr.write(str(exc_value)+'\n')
                traceback.print_tb(exc_tb)
                QtWidgets.QMessageBox.critical(self,"File open failed","Failed to open file "+filename+\
                              '\n'+str(exc_type)+'\n'+str(exc_value)+'\n'  )

        self.statusFormat = None

    @QtCore.Slot()
    def loadFileDialog(self):
        filter_list = []
        filter_list.append(self.tr("CCP4 log files (*.log *.html)"))
        filter_list.append(self.tr("CCP4 XML log files (*.xml)"))
        filter_list.append(self.tr("XMGR/Grace files (*.xmgr *.agr)"))
        filter_list.append(self.tr("3D Surface (*.surf)"))
        filter_list.append(self.tr("Any file (*.*)"))
        loadFileDialog = QtWidgets.QFileDialog()
        loadFileDialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        loadFileDialog.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        loadFileDialog.setAcceptMode(QtWidgets.QFileDialog.AcceptOpen)
        if loadFileDialog.exec_():
            for f in loadFileDialog.selectedFiles():
                    self.loadFile(str(f))

    def statusAsXML(self):
        status_xml = ""
        header ="""<?xml version="1.0" encoding="UTF-8" ?>\n"""
        status_xml += header

        NSMAP = {'xsi':"http://www.w3.org/2001/XMLSchema-instance"}
        NS = NSMAP['xsi']
        location_attribute = '{%s}noNamespaceSchemaLocation' % NS
        tree = etree.Element("CCP4ApplicationOutput",nsmap = NSMAP,attrib={location_attribute: 'http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd'})
        tree.append(self.fontsToEtree())
        for surf in self.surface_etrees:
            tree.append(surf)
        for i in range(self.table_combo.count()):
            cw = self.graph.widget(i)
            if cw in self.table_etrees:
                t = str(self.table_combo.itemText(i))
                tableTree = etree.Element('CCP4Table',title=t)
                tree.append(tableTree)
                for data in self.table_etrees[cw].xpath('data'):
                    tableTree.append(copy.deepcopy(data))
                for eltype in ["headers"]:
                    for el in self.table_etrees[cw].xpath(eltype):
                        tableTree.append(copy.deepcopy(el))
                for j in range(cw.count()):
                    cw2 = cw.widget(j)
                    if hasattr(cw2,"canvas"):
                        #print "Have a canvas, already drawn",j
                        el = self.table_etrees[cw].xpath('plot')[j]
                        if hasattr(cw2.canvas,"x_scaling") and cw2.canvas.x_scaling is scaling_functions['oneoversqrt']:
                            el.append(etree.Element('xscale'))
                            el.xpath("xscale")[0].text = 'oneoversqrt'
                        if hasattr(cw2.canvas,"y_scaling") and cw2.canvas.y_scaling is scaling_functions['oneoversqrt']:
                            el.append(etree.Element('yscale'))
                            el.xpath("yscale")[0].text = 'oneoversqrt'
                        if cw2.canvas.custom_xlim:
                            if not len(el.xpath("xrange")) > 0:
                                el.append(etree.Element('xrange'))
                            el.xpath("xrange")[0].set("min",str(cw2.canvas.custom_xlim[0]))
                            el.xpath("xrange")[0].set("max",str(cw2.canvas.custom_xlim[1]))
                        if cw2.canvas.custom_ylim:
                            for yrangeel in el.xpath("yrange"):
                                if "rightaxis" not in yrangeel.attrib:
                                    theEl = yrangeel
                                    break
                                if yrangeel.attrib['rightaxis'] != "true" and yrangeel.attrib['rightaxis'] != "1":
                                    theEl = yrangeel
                                    break
                            else:
                                el.append(etree.Element('yrange'))
                                theEl = el.xpath("yrange")[-1]
                            if len(el.xpath("yrange"))==0:
                                el.append(etree.Element('yrange'))
                                theEl = el.xpath("yrange")[-1]
                            theEl.set("min",str(cw2.canvas.custom_ylim[0]))
                            theEl.set("max",str(cw2.canvas.custom_ylim[1]))
                        if cw2.canvas.title:
                            if not len(el.xpath("title")) > 0:
                                el.append(etree.Element('title'))
                            el.xpath("title")[0].text = cw2.canvas.title
                        if cw2.canvas.xlabel:
                            if not len(el.xpath("xlabel")) > 0:
                                el.append(etree.Element('xlabel'))
                            el.xpath("xlabel")[0].text = cw2.canvas.xlabel
                        if cw2.canvas.ylabel:
                            if not len(el.xpath("ylabel")) > 0:
                                el.append(etree.Element('ylabel'))
                            el.xpath("ylabel")[0].text = cw2.canvas.ylabel

                        if cw2.canvas.xbreak:
                            if not len(el.xpath("xbreaks")) > 0:
                                el.append(etree.Element('xbreaks'))
                            for node in el.xpath("xbreaks") :
                                for n in node:
                                    node.remove(n)
                            for br in cw2.canvas.xbreak:
                                xmlbr = etree.Element('break')
                                xmlbr.set("min",str(br[0]))
                                xmlbr.set("max",str(br[1]))
                                for node in el.xpath("xbreaks") :
                                    node.append(xmlbr)

                        if cw2.canvas.ybreak:
                            if not len(el.xpath("ybreaks")) > 0:
                                el.append(etree.Element('ybreaks'))
                            for node in el.xpath("ybreaks") :
                                for n in node:
                                    node.remove(n)
                            for br in cw2.canvas.ybreak:
                                xmlbr = etree.Element('break')
                                xmlbr.set("min",str(br[0]))
                                xmlbr.set("max",str(br[1]))
                                for node in el.xpath("ybreaks") :
                                    node.append(xmlbr)

                        if hasattr(cw2.canvas,"custom_legend_position"):
                            if not len(el.xpath("legendposition")) > 0:
                                el.append(etree.Element('legendposition'))
                            for node in el.xpath("legendposition") :
                                node.set("x",str(cw2.canvas.custom_legend_position[0]))
                                node.set("y",str(cw2.canvas.custom_legend_position[1]))

                        if cw2.canvas.ax[0].get_legend() and not cw2.canvas.ax[0].legend_.get_visible():
                            if not len(el.xpath("showlegend")) > 0:
                                el.append(etree.Element('showlegend'))
                            for node in el.xpath("showlegend") :
                                node.text = "false"
                        if hasattr(cw2.canvas,"fix_aspect_ratio") and cw2.canvas.fix_aspect_ratio:
                            print("should be saving with fixed aspect ratio")
                            if not len(el.xpath("fixaspectratio")) > 0:
                                el.append(etree.Element('fixaspectratio'))
                            for node in el.xpath("fixaspectratio") :
                                node.text = "true"

                        #print 'plotline length',len(el.xpath("plotline"))
                        for plotline in el.xpath("plotline"):
                            # Just do the simple lines. Histograms will be more complex.
                            if len(self.elementToPlotMap[plotline])==1:
                                if len(plotline.xpath("label")) == 0:
                                    plotline.append(etree.Element('label'))
                                for node in plotline.xpath("label") :
                                    node.text = self.elementToPlotMap[plotline][0].get_label()
                                if len(plotline.xpath("symbol")) == 0:
                                    plotline.append(etree.Element('symbol'))
                                for node in plotline.xpath("symbol") :
                                    node.text = self.elementToPlotMap[plotline][0].get_marker()
                                if len(plotline.xpath("symbolsize")) == 0:
                                    plotline.append(etree.Element('symbolsize'))
                                for node in plotline.xpath("symbolsize") :
                                    node.text = str(int(self.elementToPlotMap[plotline][0].get_markersize()))
                                if len(plotline.xpath("markeredgewidth")) == 0:
                                    plotline.append(etree.Element('markeredgewidth'))
                                for node in plotline.xpath("markeredgewidth") :
                                    node.text = str(int(self.elementToPlotMap[plotline][0].get_markeredgewidth()))
                                if len(plotline.xpath("colour")) == 0:
                                    plotline.append(etree.Element('colour'))
                                for node in plotline.xpath("colour") :
                                    node.text = self.elementToPlotMap[plotline][0].get_color()
                                if len(plotline.xpath("linestyle")) == 0:
                                    plotline.append(etree.Element('linestyle'))
                                for node in plotline.xpath("linestyle") :
                                    node.text = self.elementToPlotMap[plotline][0].get_linestyle()
                                if len(plotline.xpath("linesize")) == 0:
                                    plotline.append(etree.Element('linesize'))
                                for node in plotline.xpath("linesize") :
                                    node.text = str(int(self.elementToPlotMap[plotline][0].get_linewidth()))
                                if not self.elementToPlotMap[plotline][0].get_visible():
                                    if len(plotline.xpath("visible")) == 0:
                                        plotline.append(etree.Element('visible'))
                                    for node in plotline.xpath("visible") :
                                        node.text = "false"
                                if hasattr(self.elementToPlotMap[plotline][0],"hide_from_legend") and self.elementToPlotMap[plotline][0].hide_from_legend:
                                    if len(plotline.xpath("showinlegend")) == 0:
                                        plotline.append(etree.Element('showinlegend'))
                                    for node in plotline.xpath("showinlegend") :
                                        node.text = "false"

                        tableTree.append(copy.deepcopy(el))

                    else:
                        el = self.table_etrees[cw].xpath('plot')[j]
                        if hasattr(cw2,"_do_the_plot_"):
                            if hasattr(cw2._do_the_plot_,"args"):
                                if len(cw2._do_the_plot_.args)>0:
                                    if hasattr(cw2._do_the_plot_.args[0],"x_scaling"):
                                        if cw2._do_the_plot_.args[0].x_scaling is scaling_functions['oneoversqrt']:
                                            if len(el.xpath("xscale"))==0:
                                                el.append(etree.Element('xscale'))
                                                el.xpath("xscale")[0].text = 'oneoversqrt'
                                    if hasattr(cw2._do_the_plot_.args[0],"y_scaling"):
                                        if cw2._do_the_plot_.args[0].y_scaling is scaling_functions['oneoversqrt']:
                                            if len(el.xpath("yscale"))==0:
                                                el.append(etree.Element('yscale'))
                                                el.xpath("yscale")[0].text = 'oneoversqrt'
                        tableTree.append(copy.deepcopy(el))



        if sys.version_info > (3,0):
            status_xml += etree.tostring(tree,encoding='utf-8', pretty_print=True).decode("utf-8")
        else:
            status_xml += etree.tostring(tree,encoding='utf-8', pretty_print=True)
        return status_xml

    @QtCore.Slot(bool,bool)
    def saveFileDialog(self,saveAll=False,saveStatus=False):
        """
        if saveStatus:
            print "Saving status"
            status = self.statusAsXML()
            #print status; sys.stdout.flush()
            return
        """
        filter_list = []
        if not saveAll:
            filter_list.append(self.tr("Postscript (*.ps)"))
        #filter_list.append(self.tr("Encapsulated Postscript (*.eps)"))
        filter_list.append(self.tr("Portable Document Format (*.pdf)"))
        if not saveAll:
            filter_list.append(self.tr("PNG files (*.png)"))

        if saveStatus:
            filter_list = [self.tr("XML files (*.xml)")]

        formats_py = []
        if not saveAll and not saveStatus:
            formats = QtGui.QImageWriter.supportedImageFormats() 
            for format in formats:
                if str(format) != 'png':
                    filter_list.append(self.tr(str(format).upper()+" files (*."+str(format)+")"))
                    formats_py.append(str(format))
        formats_py = tuple(formats_py)

        saveFileDialog = QtWidgets.QFileDialog()
        saveFileDialog.setOption(QtWidgets.QFileDialog.DontUseNativeDialog)
        saveFileDialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        saveFileDialog.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        if saveAll:
              saveFileDialog.setWindowTitle('Save all figures')
        elif saveStatus:
              saveFileDialog.setWindowTitle('Save status')
        else:
              saveFileDialog.setWindowTitle('Save figure')
        saveFileDialog.setNameFilters(filter_list)
        if saveFileDialog.exec_():
            filename = str(saveFileDialog.selectedFiles()[0])
            if not filename.endswith('.xml') and not filename.endswith('.XML') and not filename.endswith('.pdf') and not filename.endswith('.PDF') and not filename.endswith('.png') and not filename.endswith('.PNG') and not filename.endswith('.ps') and not filename.endswith('.PS') and not filename.endswith(formats_py):
                selectedFilter = str(saveFileDialog.selectedNameFilter())
                splitFilters = selectedFilter[selectedFilter.find('(')+1:selectedFilter.find(')')].strip().split(' ')
                if len(splitFilters) == 1 and splitFilters[0].startswith('*.'):
                    filename = filename + splitFilters[0][1:]
                else:
                    filename = filename + 'png'
            dpi = 100

            if saveStatus:
                print("Saving status")
                status = self.statusAsXML()
                f = open(filename,"w")
                f.write(status)
                f.close()
                return
            if saveAll:
                pp = PdfPages(filename)
                currenti = self.graph.currentIndex()
                currentj = self.graph.currentWidget().currentIndex()
                for i in range(self.graph.count()):
                    self.graph.setCurrentIndex(i)
                    cw = self.graph.widget(i)
                    for j in range(cw.count()):
                        cw.setCurrentIndex(j)
                        plot = cw.widget(j)
                        if hasattr(plot,"canvas"):
                            plot.canvas.draw()
                            plot.canvas.fig.savefig(pp, format='pdf')
                        elif hasattr(plot,"_do_the_plot_"):
                            newPlot = plot._do_the_plot_()
                            newPlot.show()
                            newPlot.canvas.draw()
                            newPlot.canvas.fig.savefig(pp, format='pdf')
                pp.close()
                self.graph.setCurrentIndex(currenti)
                self.graph.currentWidget().setCurrentIndex(currentj)
                return

            plot = self.graph.currentWidget().currentWidget()
            if not filename.endswith('.pdf') and not filename.endswith('.PDF') and not filename.endswith('.png') and not filename.endswith('.PNG') and not filename.endswith('.ps') and not filename.endswith('.PS'):
                if filename.endswith(formats_py):
                    import tempfile
                    f = tempfile.NamedTemporaryFile(suffix=".png",prefix="ccp4mg"+str(os.getpid()))
                    fn = f.name
                    f.close()
                    plot.canvas.fig.savefig(fn,dpi=dpi)
                    im = QtGui.QImage(fn)
                    #print "Save ", filename
                    im.save(filename)
                    
            else:
                print("Save ", filename)
                if filename.endswith('.pdf') or filename.endswith('.PDF'):
                    pp = PdfPages(filename)
                    plot.canvas.fig.savefig(pp, format='pdf')
                    pp.close()
                else:
                    plot.canvas.fig.savefig(filename,dpi=dpi)

    def updateStatusFormat(self):
      import math
      def logBase(num):
        try:
          l = math.log10(num)
        except:
          # Probably trying for log(0)
          l = 0.0
        if l<0:
          return int(l)-1
        else:
          return int(l)
      def format(minX,maxX):
        if maxX>3 or minX<-6:
          #Use scientific notation
          xF = '%.3E'
        elif minX < 0:
          if maxX<0:
            xF = '%'+str(abs(minX))+'.'+str(abs(minX))+'f'
          else:
            xF = '%'+str(abs(minX)+abs(maxX))+'.'+str(abs(minX))+'f'
        else:
          xf = '%'+str(abs(maxX))+'.0f'
        return xF

      try:
        plot = self.graph.currentWidget().currentWidget()
        #The range of axes
        xr = math.fabs(plot.canvas.ax[-1].get_xlim()[1]-plot.canvas.ax[0].get_xlim()[0])
        yr = math.fabs(plot.canvas.ax[-1].get_ylim()[1]-plot.canvas.ax[0].get_ylim()[0])
        #The max abs value of axes
        xx = max(math.fabs(plot.canvas.ax[-1].get_xlim()[1]),math.fabs(plot.canvas.ax[0].get_xlim()[0]))
        yy = max(math.fabs(plot.canvas.ax[-1].get_ylim()[1]),math.fabs(plot.canvas.ax[0].get_ylim()[0]))
        #print 'updateStatusFormat',xx,yy,xr,yr
        maxX = logBase(xx)
        minX = logBase(xr) - 3
        xF = format(minX,maxX)
        maxY = logBase(yy)
        minY = logBase(yr) - 3
        yF = format(minY,maxY)
        #print 'updateStatusFormat',minX,maxX,xF,minY,maxY,yF
          
        self.statusFormat = [ xF, yF ]
      except:
        self.statusFormat = [ '%12.6f', '%12.6f' ]



    @QtCore.Slot(float,float,dict)
    def updateStatusWithPickXandY(self,x_scaling,y_scaling,pos):
        plot = self.graph.currentWidget().currentWidget()
        if plot.canvas.movingLegend == True:
            self.statusBar.showMessage(self.tr("Moving legend, press escape to cancel"))
            return
            
        if pos and type(pos)==dict and 'x' in pos and 'y' in pos and 'label' in pos:
            name = pos['label']
            if x_scaling:
                x = "x = " + x_scaling(pos['x'],None,format=self.statusFormat[0])
            else:
                x = "x = "+self.statusFormat[0] % (pos['x'])
            if y_scaling:
                y = "y = " + y_scaling(pos['y'],None,format=self.statusFormat[1])
            else:
                y = "y = "+self.statusFormat[1] % (pos['y'])
            text = name + " " + x + ", " + y
            self.statusBar.showMessage(text)

    @QtCore.Slot(float,float,dict)
    def updateStatusWithXandY(self,x_scaling,y_scaling,pos):
        if self.statusFormat is None: self.updateStatusFormat()
        plot = self.graph.currentWidget().currentWidget()
        if plot.canvas.movingLegend == True:
            self.statusBar.showMessage(self.tr("Moving legend, press escape to cancel"))
            return
            
        if pos and type(pos)==tuple and len(pos)==2:
            if x_scaling:
                x = "x = " + x_scaling(pos[0],None,format=self.statusFormat[0])
            else:
                x = "x = "+self.statusFormat[0] % (pos[0])
            if y_scaling:
                y = "y = " + y_scaling(pos[1],None,format=self.statusFormat[1])
            else:
                y = "y = "+self.statusFormat[1] % (pos[1])
            text = x + ", " + y
            self.statusBar.showMessage(text)

    @QtCore.Slot(int)
    def setCurrentData(self,idx):
        print('setCurrentData', idx)
        #traceback.print_stack()
        if idx < 0: return
        #print self.table_combo.currentIndex(), idx, self.graph.currentIndex(), self.graph.count()
        if not self.graph.currentWidget():
            #print "returning"
            return
        if self.graph.currentWidget().widget(idx) and hasattr(self.graph.currentWidget().widget(idx),"_do_the_plot_"):
            #print self.graph.currentWidget().widget(idx)
            newWidget = self.graph.currentWidget().widget(idx)._do_the_plot_()
            #print newWidget
            oldWidget = self.graph.currentWidget().widget(idx)
            self.graph.currentWidget().removeWidget(oldWidget)
            oldWidget.close()
            self.graph.currentWidget().insertWidget(idx,newWidget)
        self.graph.currentWidget().setCurrentIndex(idx)
        #print self.graph.currentWidget().currentWidget()
        if sys.platform == "darwin":
            # This is a hack to get around a(n almost certainly Qt) bug on OS X.
            # The first widget on each stack seems not to get mouse focus unless
            # the winfow has been resized. This does just that.
            self.resize(self.width()+1,self.height())
            self.repaint()
            self.resize(self.width()-1,self.height())
            self.repaint()
        self.graph.currentWidget().currentWidget().canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(self.graph.currentWidget().currentWidget().canvas.width(),self.graph.currentWidget().currentWidget().canvas.height()),QtCore.QSize(self.graph.currentWidget().currentWidget().canvas.width(),self.graph.currentWidget().currentWidget().canvas.height())))
        self.graph.currentWidget().currentWidget().canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(self.graph.currentWidget().currentWidget().canvas.width(),self.graph.currentWidget().currentWidget().canvas.height()),QtCore.QSize(self.graph.currentWidget().currentWidget().canvas.width(),self.graph.currentWidget().currentWidget().canvas.height())))
        self.statusFormat=None
        self.currentDataChanged.emit(idx)

    @QtCore.Slot(int)
    def setCurrentTable(self,idx):
        print('setCurrentTable',self.graph_lists)
        self.data_combo.clear()
        for g in self.graph_lists[idx]:
            #print "Adding", g[0]
            self.data_combo.addItem(g[0])
            self.comboItemsChanged()
        self.data_combo.setCurrentIndex(0)
        self.graph.setCurrentIndex(idx)
        self.setCurrentData(0)
        self.currentTableChanged.emit(idx)

    def addSurfaceFromEtree(self,tree=None,graph_uuid=None):
        title = tree.attrib.get('title','')
        columns = int(tree.attrib.get('columns',''))
        rows = int(tree.attrib.get('rows',''))
        Z = []
        Zi = []
        theText = tree.text.split()
        while theText:
            while len(Zi)<columns:
                Zi.append(float(theText.pop(0)))
            Z.append(Zi)
            Zi = []
        return self.addSurfaceData(Z,title,graph_uuid)

    def addSurfaceFile(self,filename,graph_uuid=None):
        f = open(filename)
        d = f.readlines()
        f.close()
        Z = []
        surfaceTree = etree.Element('CCP4Surface',title=filename, rows=str(len(d)),columns=str(len(d[0].split())))
        surfaceTree.text = ''
        for l in d:
            surfaceTree.text += l
        self.surface_etrees.append(surfaceTree)
        for l in d:
            Zi = []
            for v in l.strip().split():
                Zi.append(float(v))
            Z.append(Zi)
        return self.addSurfaceData(Z,filename,graph_uuid)

    def addSurfaceData(self,Z,dataname,graph_uuid=None):
        mind = 1e12
        maxd = -1e12

        for zi in Z:
           for v in zi:
               if float(v)> maxd:
                    maxd = float(v)
               if float(v)< mind:
                    mind = float(v)
        levels = []
        diff = float(maxd-mind)

        ncont = 50
        step = diff/ncont

        for i in range(ncont):
            levels.append(mind+i*step)
        graph = StackedWidget()
        if graph_uuid != None:
             graph.setObjectName(graph_uuid)
        self.graph.addWidget(graph)
        graph_lists = []
        graph_lists.append((dataname,''))
        plot = QtMatplotlibWidget(title=dataname)
        plot.setDefaultSymbolSize(self.defaultSymbolSize)
        plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,None,None))
        graph.addWidget(plot)
        from pylab import cm
        plot.canvas.ax[0].contour(Z,80,cmap=cm.bone,linewidths=1,origin='lower')
        CS = plot.canvas.ax[0].contourf(Z,80,linestyles='dotted',origin='lower')
        plot.canvas.fig.colorbar(CS)
        plot.canvas._type = 'surface'
        plot.canvas.ax[0].set_title(dataname)
        self.graph_lists.append(graph_lists)
        self.table_combo.addItem(str(dataname))
        self.comboItemsChanged()
        self.table_combo.setCurrentIndex(self.table_combo.count()-1)
        return graph

    def addGraceFile(self,filename,graph_uuid=None):

        f = open(filename)
        lines = f.readlines()
        f.close()

        ig = 0
        iset = 1
        vals = []
        xvals = []
        graph_lists = []
        graph = StackedWidget()
        if graph_uuid != None:
            graph.setObjectName(graph_uuid)
        self.graph.addWidget(graph)
        t = etree.Element("CCP4Table")
        theDataEl = etree.Element("data")
        self.table_etrees[graph] = t
        inGraph = False
        graph_lists.append((filename+' '+str(iset),''))
        p = etree.Element('plot')
        theTitle = etree.Element("title")
        theTitle.text = filename+' '+str(iset)
        p.append(theTitle)
        t.append(p)
        plot = QtMatplotlibWidget(title=filename+' '+str(iset))
        plot.setDefaultSymbolSize(self.defaultSymbolSize)
        plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,None,None))
        plot.canvas.DetectPointMoveEvent.connect(functools.partial(self.updateStatusWithPickXandY,None,None))
        plot.canvas.DetectPointClickEvent.connect(graph.DetectPointClickEvent.emit)
        plot.canvas.DetectPointMoveEvent.connect(graph.DetectPointMoveEvent.emit)
        plot.canvas.RegionSelected.connect(graph.RegionSelected.emit)
        graph.addWidget(plot)
        currentLineStyle = None
        currentColour = None
        currentMarker = None
        currentLegendString = None
        currentMarkerScaling = 1
        haveLegend = False
        xRange = [None,None]
        yRange = [None,None]

        for l in lines:
            if l.startswith('@') and len(l.split()) == 4 and l.split()[2] == 'type':
                inGraph = True
            elif l.startswith('&'):
                inGraph = False
                label = filename+' '+str(ig)
                if currentLegendString:
                    label = currentLegendString

                markersize = int(self.defaultSymbolSize*currentMarkerScaling)
                if len(vals)>0:
                    if len(xvals)==len(vals):
                        if type(vals[0]) == list:
                            for j in range(len(vals[0])):
                                v = []
                                for i in range(len(xvals)):
                                    v.append(vals[i][j])
                                if currentColour:
                                    color = currentColour
                                else:
                                    color = self.colors[ig % len(self.colors)]
                                if currentLineStyle:
                                    style = currentLineStyle
                                    if currentMarker:
                                        thePlot = plot.canvas.plot(xvals,v,label=label,picker=markersize, marker=currentMarker, linestyle=style, color=color, markersize=markersize)
                                    else:
                                        thePlot = plot.canvas.plot(xvals,v,label=label,picker=markersize, linestyle=style, color=color, markersize=markersize)
                                else:
                                    style = self.styles[ig % len(self.styles)]
                                    thePlot = plot.canvas.plot(xvals,v,label=label,picker=markersize, marker=style, color=color, markersize=markersize)
                                for line in thePlot:
                                    if not currentLegendString: line.hide_from_legend = True
                                ig = ig + 1
                        else:
                            if currentColour:
                                color = currentColour
                            else:
                                color = self.colors[ig % len(self.colors)]
                            if currentLineStyle:
                                style = currentLineStyle
                                if currentMarker:
                                    thePlot = plot.canvas.plot(xvals,vals,label=label,picker=markersize, marker=currentMarker, linestyle=style, color=color, markersize=markersize)
                                else:
                                    thePlot = plot.canvas.plot(xvals,vals,label=label,picker=markersize, linestyle=style, color=color, markersize=markersize)
                            else:
                                style = self.styles[ig % len(self.styles)]
                                thePlot = plot.canvas.plot(xvals,vals,label=label,picker=markersize, marker=style, color=color, markersize=markersize)
                            for line in thePlot:
                                if not currentLegendString: line.hide_from_legend = True
                            ig = ig + 1
                    else:
                        if currentColour:
                            color = currentColour
                        else:
                            color = self.colors[ig % len(self.colors)]
                        if currentLineStyle:
                            style = currentLineStyle
                            if currentMarker:
                                thePlot = plot.canvas.plot(vals,label=label,picker=markersize, marker=currentMarker, linestyle=style, color=color, markersize=markersize)
                            else:
                                thePlot = plot.canvas.plot(vals,label=label,picker=markersize, linestyle=style, color=color, markersize=markersize)
                        else:
                            style = self.styles[ig % len(self.styles)]
                            thePlot = plot.canvas.plot(vals,label=label,picker=markersize, marker=style, color=color, markersize=markersize)
                        for line in thePlot:
                            if not currentLegendString: line.hide_from_legend = True
                        ig = ig + 1

                    #print "adding plotline",label,ig,p
                    plotline = etree.Element('plotline',xcol=str(2*(ig-1)),ycol=str(2*(ig-1)+1))
                    p.append(plotline)
                    self.elementToPlotMap[plotline] = thePlot
                    
                vals = []
                xvals = []
                currentLineStyle = None
                currentColour = None
                currentMarker = None
                currentLegendString = None
                currentMarkerScaling = 1
            elif inGraph and l.startswith('@'):
                if len(l.split()) == 4 and l.split()[2] == 'linestyle':
                    lineStyle = l.split()[3]
                    if int(lineStyle) == 0:
                        currentLineStyle = '.'
                        currentMarker = 'o'
                        print(".")
                    if int(lineStyle) == 1:
                        currentLineStyle = '-'
                elif len(l.split()) == 5 and l.split()[2] == 'symbol' and l.split()[3] == 'size':
                    currentMarkerScaling = 1+float(l.split()[4])
                elif len(l.split()) > 4 and l.split()[1] == 'legend' and l.split()[2] == 'string':
                    currentLegendString = l[l.find('"'):].strip().strip('\'"')
                    haveLegend = True
                elif len(l.split()) == 4 and l.split()[2] == 'color':
                    currentColour = self.colors[int(l.split()[3])-1]
                else:
                    print("IGNORE",l)
            elif inGraph:
                thisVal = l.strip().split()
                if len(thisVal) == 1:
                    vals.append(float(thisVal[0]))
                elif len(thisVal) == 2:
                    xvals.append(float(thisVal[0]))
                    vals.append(float(thisVal[1]))
                elif len(thisVal) > 2:
                    xvals.append(float(thisVal[0]))
                    tvals = []
                    for v in thisVal[1:]:
                        tvals.append(float(v))
                    vals.append(tvals)
            elif l.startswith('@'):
                if len(l.split()) > 2 and l.split()[1] == 'title':
                    title = l[l.find('title')+len('title'):].strip().strip('\'"')
                    plot.canvas.ax[0].set_title(title)
                    plot.canvas.title = title
                elif len(l.split()) > 3 and l.split()[1] == 'world':
                    if l.split()[2] == 'xmin':
                        try:
                            xRange[0] = float(l.split()[3])
                        except:
                            xRange[0] = None
                    if l.split()[2] == 'xmax':
                        try:
                            xRange[1] = float(l.split()[3])
                        except:
                            xRange[1] = None
                    if l.split()[2] == 'ymin':
                        try:
                            yRange[0] = float(l.split()[3])
                        except:
                            yRange[0] = None
                    if l.split()[2] == 'ymax':
                        try:
                            yRange[1] = float(l.split()[3])
                        except:
                            yRange[1] = None
                else:
                    print("IGNORE",l)

        if xRange[0] and xRange[1] and yRange[0] and yRange[1]:
            plot.canvas.ax[0].set_xlim(xRange[0],xRange[1])
            plot.canvas.custom_xlim = xRange[0],xRange[1]
            plot.canvas.ax[0].set_ylim(yRange[0],yRange[1])
            plot.canvas.custom_ylim = yRange[0],yRange[1]
            if abs((xRange[1]-xRange[0])-(yRange[1]-yRange[0]))<1e-5:
                plot.canvas.fix_aspect_ratio = True
                plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))
                plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))

        #print p.xpath('plotline')
        theData = []
        maxXData, maxYData = 0,0
        for plotline in p.xpath('plotline'):
            maxXData, maxYData = len(self.elementToPlotMap[plotline][0].get_xdata()) , len(self.elementToPlotMap[plotline][0].get_ydata())
        theDataStr = ''
        for i,plotline in zip(list(range(len(p.xpath('plotline')))),p.xpath('plotline')):
            theData.append(self.elementToPlotMap[plotline][0].get_xdata())
            theData.append(self.elementToPlotMap[plotline][0].get_ydata())
        #print theData
        theData2 = map(None, *theData)
        theDataText = ''
        for row in theData2:
            theRow = ''
            for c in row:
                if c !=None:
                    theRow += "%.6f " %c
                else:
                    theRow += "%6s " % "-"
            theDataText += theRow + '\n'
        theDataEl.text = theDataText
        t.append(theDataEl)

        if haveLegend: plot.canvas.legend()
        self.graph_lists.append(graph_lists)
        self.table_combo.addItem(str(filename))
        self.comboItemsChanged()
        self.table_combo.setCurrentIndex(self.table_combo.count()-1)
        return graph

    def addPlainFile(self,filename,graph_uuid=None,hist=False):

        f = open(filename)
        lines = f.readlines()
        f.close()
        return self.addPlainData(lines,filename,graph_uuid,hist=hist)

    def addPlainData(self,lines,dataname,graph_uuid=None,hist=False):
        ig = 0
        iset = 1
        vals = []
        xvals = []
        graph_lists = []
        graph = StackedWidget()
        if graph_uuid != None:
            graph.setObjectName(graph_uuid)
        self.graph.addWidget(graph)
        for l in lines:
            thisVal = l.strip().split()
            if len(thisVal) == 1:
                vals.append(float(thisVal[0]))
            elif len(thisVal) == 2:
                xvals.append(float(thisVal[0]))
                vals.append(float(thisVal[1]))
            elif len(thisVal) > 2:
                xvals.append(float(thisVal[0]))
                tvals = []
                for v in thisVal[1:]:
                    tvals.append(float(v))
                vals.append(tvals)
                    
            elif len(thisVal) == 0:
                if len(vals)>0:
                    graph_lists.append((dataname,''))

                    plot = QtMatplotlibWidget(title=dataname)
                    plot.setDefaultSymbolSize(self.defaultSymbolSize)
                    plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,None,None))
                    plot.canvas.DetectPointMoveEvent.connect(functools.partial(self.updateStatusWithPickXandY,None,None))
                    plot.canvas.DetectPointClickEvent.connect(graph.DetectPointClickEvent.emit)
                    plot.canvas.DetectPointMoveEvent.connect(graph.DetectPointMoveEvent.emit)
                    plot.canvas.RegionSelected.connect(graph.RegionSelected.emit)
                    graph.addWidget(plot)

                    if len(xvals)==len(vals):
                        if type(vals[0]) == list:
                            for j in range(len(vals[0])):
                                v = []
                                for i in range(len(xvals)):
                                    v.append(vals[i][j])
                                color = self.colors[ig % len(self.colors)]
                                style = self.styles[ig % len(self.styles)]
                                plot.canvas.plot(xvals,v,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
                                ig = ig + 1
                        else:
                            color = self.colors[ig % len(self.colors)]
                            style = self.styles[ig % len(self.styles)]
                            plot.canvas.plot(xvals,vals,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
                            ig = ig + 1
                    else:
                        if hist:
                            color = self.colors[ig % len(self.colors)]
                            mu, sigma = 0,1
                            if len(vals)>0:
                                mu = numpy.mean(vals)
                                sigma = numpy.std(vals)
                                maxx = numpy.max(vals)
                            #n, bins, patches = plot.canvas.ax.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname+' '+str(ig))
                            if(ig)>0:
                                n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname+' '+str(ig+1))
                            else:
                                n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname)
                            if len(patches)>0:
                                # Maybe need a Histogram class ?
                                self.histograms.append({'normed':False,'vals':vals,'maxx':maxx,'patches':patches,'bins':bins,'mu':mu, 'sigma':sigma})
                        else:
                            color = self.colors[ig % len(self.colors)]
                            style = self.styles[ig % len(self.styles)]
                            plot.canvas.plot(vals,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
                        ig = ig + 1
                    vals = []
                    xvals = []
                    iset = iset + 1
                    ig = 0

        if len(vals)>0:
            if iset>0:
                graph_lists.append((dataname+' '+str(iset+1),''))
            else:
                graph_lists.append((dataname,''))

            plot = QtMatplotlibWidget(title=dataname+' '+str(iset))
            plot.setDefaultSymbolSize(self.defaultSymbolSize)
            plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,None,None))
            plot.canvas.DetectPointMoveEvent.connect(functools.partial(self.updateStatusWithPickXandY,None,None))
            plot.canvas.DetectPointClickEvent.connect(graph.DetectPointClickEvent.emit)
            plot.canvas.DetectPointMoveEvent.connect(graph.DetectPointMoveEvent.emit)
            plot.canvas.RegionSelected.connect(graph.RegionSelected.emit)
            graph.addWidget(plot)

            if len(xvals)==len(vals):
                if type(vals[0]) == list:
                    for j in range(len(vals[0])):
                        v = []
                        for i in range(len(xvals)):
                            v.append(vals[i][j])
                        color = self.colors[ig % len(self.colors)]
                        style = self.styles[ig % len(self.styles)]
                        plot.canvas.plot(xvals,v,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
                        ig = ig + 1
                else:
                    color = self.colors[ig % len(self.colors)]
                    style = self.styles[ig % len(self.styles)]
                    plot.canvas.plot(xvals,vals,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
            else:
                if hist:
                    color = self.colors[ig % len(self.colors)]
                    mu, sigma = 0,1
                    if len(vals)>0:
                        mu = numpy.mean(vals)
                        sigma = numpy.std(vals)
                        maxx = numpy.max(vals)
                    if(ig)>0:
                        n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname+' '+str(ig+1))
                    else:
                        n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname)
                    if len(patches)>0:
                        # Maybe need a Histogram class ?
                        self.histograms.append({'normed':False,'vals':vals,'maxx':maxx,'patches':patches,'bins':bins,'mu':mu, 'sigma':sigma})
                else:
                    color = self.colors[ig % len(self.colors)]
                    style = self.styles[ig % len(self.styles)]
                    plot.canvas.plot(vals,label=dataname+' '+str(ig),picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)

        self.graph_lists.append(graph_lists)
        self.table_combo.addItem(str(dataname))
        self.comboItemsChanged()
        self.table_combo.setCurrentIndex(self.table_combo.count()-1)
        return graph

    def addCCP4ReportFile(self,fname,select=None):
        #print('CCP4Table.addCCP4ReportFile',fname)
        from lxml import etree
        try:
            doc = openFileToEtree(fname)
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')
            traceback.print_tb(exc_tb)
            QtWidgets.QMessageBox.warning(self,"Report HTML parse error","File "+fname+" does not seem to be valid HTML")
            return []
        graphList = []
        #print 'CCP4Table.addCCP4ReportFile doc',etree.tostring(doc,pretty_print=True)
        for tag in ('ccp4_data','{http://www.ccp4.ac.uk/ccp4ns}ccp4_data'):
            for tableEle in doc.iter(tag = tag):
                #print 'CCP4Table.addCCP4ReportFile found ccp4_data',tableEle.tag
                if select is None or tableEle.get('id',None) in select:
                    try:
                        graph = self.addTableFromEtree(tableEle)
                        graphList.append(graph)
                    except:
                        print('ERROR loading table',tableEle.tag)

        # This section deals with data coming from external XML files.
        for tableEle in doc.iter(tag = "div"):
                        dataData = tableEle.attrib.get('data-data','')
                        # We allow for the bizarre situation where the string ends with .xml but is actually a div id.
                        if len(dataData.strip())>0 and dataData.endswith(".xml") and len(doc.findall(".//*[@id='"+dataData+"']"))==0:
                            fnameXML = os.path.join(os.path.dirname(fname),dataData)
                            try:
                                docXML = openFileToEtree(fnameXML)
                                tables = docXML.xpath('/CCP4ApplicationOutput/CCP4Table')
                                for t in tables:
                                    #print CCP4Table.addCCP4ReportFile XML table',etree.tostring(doc,pretty_print=True)
                                    graph = self.addTableFromEtree(t)
                                    graphList.append(graph)
                            except:
                                print('POSSIBLE ERROR loading table file',fnameXML,".")
                                print("This may not be a table file, in which case there is no problem.")
        return graphList
    
    def addXMLCCP4TableFile(self,fname,graph_uuid=None,xsd=None):
        try:
            doc = openFileToEtree(fname)
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')
            traceback.print_tb(exc_tb)
            QtWidgets.QMessageBox.warning(self,"XML parse error","File "+fname+" does not seem to be valid XML")
            return []

        if xsd:
            xmlschema_doc = None
            try:
                    xmlschema_doc = openFileToEtree(xsd)
            except:
                    exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                    sys.stderr.write(str(exc_type)+'\n')
                    sys.stderr.write(str(exc_value)+'\n')
                    traceback.print_tb(exc_tb)
                    QtWidgets.QMessageBox.warning(self,"XML validation error","File "+xsd+ " is not a valid schema")

            if xmlschema_doc is not None:
                try:
                    xmlschema = etree.XMLSchema(xmlschema_doc)
                    xmlschema.assertValid(doc)
                    print("Validated successfully against",xsd)
                except:
                    exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                    sys.stderr.write(str(exc_type)+'\n')
                    sys.stderr.write(str(exc_value)+'\n')
                    traceback.print_tb(exc_tb)
                    QtWidgets.QMessageBox.warning(self,"XML validation error","File "+fname+ " does not conform to schema\n\n"+str(exc_type)+'\n\n'+str(exc_value))

        graphs = []
        surfaces = doc.xpath('/CCP4ApplicationOutput/CCP4Surface')
        for surf in surfaces:
            graph = self.addSurfaceFromEtree(surf)
            graphs.append(graph)
            self.surface_etrees.append(copy.deepcopy(surf))
        tables = doc.xpath('/CCP4ApplicationOutput/CCP4Table')
        for t in tables:
            try:
              graph = self.addTableFromEtree(t)
              graphs.append(graph)
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                sys.stderr.write(str(exc_type)+'\n')
                sys.stderr.write(str(exc_value)+'\n')
                print('ERROR loading table',t.tag)
        fonts = doc.xpath('/CCP4ApplicationOutput/Fonts')
        def elementToDict(element):
            theFont = {}
            if len(element)>0:
                titleFontElement = element[-1]
                family = titleFontElement.attrib.get('family','')
                if family:
                    theFont['family'] = family
                size = titleFontElement.attrib.get('size',0)
                if size:
                    theFont['size'] = int(size)
                weight = titleFontElement.attrib.get('weight')
                if weight:
                    theFont['weight'] = weight
                slant = titleFontElement.attrib.get('slant')
                if slant:
                    theFont['slant'] = slant
            return theFont

        titleFont = None
        axesFont = None
        axesLabelFont = None
        legendFont = None
        for font in fonts:
            titleFonts = font.xpath('titleFont')
            titleFont = elementToDict(titleFonts)
            axesFonts = font.xpath('axesTickerFont')
            axesFont = elementToDict(axesFonts)
            axesLabelFonts = font.xpath('axesLabelFont')
            axesLabelFont = elementToDict(axesLabelFonts)
            legendFonts = font.xpath('legendFont')
            legendFont = elementToDict(legendFonts)

        if titleFont or axesFont or axesLabelFont or legendFont:
            self.applyPreferences(titleFontSel=titleFont,axesFontSel=axesFont,axesLabelFontSel=axesLabelFont,legendFontSel=legendFont)

        return graphs

    def analyzeDistributionOfPoints(self,xdata,ydata,halfx,halfy):
        nTopLeft = 0
        nTopRight = 0
        nBottomLeft = 0
        nBottomRight = 0

        if len(xdata) != len(ydata):
            return

        for i in range(len(xdata)):
            if xdata[i] > halfx:
                if ydata[i]> halfy:
                    nTopRight += 1
                else:
                    nBottomRight += 1
            else:
                if ydata[i]> halfy:
                    nTopLeft += 1
                else:
                    nBottomLeft += 1

        return nTopLeft, nTopRight, nBottomLeft, nBottomRight

    def getDataArray(self,t,theId):
        for data in t.xpath('data'):
            if data.attrib.get('id',None) == theId:
                a = data.text.strip().split()
                a2 = []
                for i in a:
                    j = i.replace("NA","Nan")
                    if j == "*" or j=="-":
                         j = "NaN"
                    a2.append(j)

                array1 = numpy.array(a2,dtype=float)
                array2 = numpy.array_split(array1,len(array1)/len(t.xpath('headers')[0].text.strip().split()))
                array = numpy.array(array2)
                return array

    def addTableFromEtreeAndSetVisible(self,table):
        oldCount = self.table_combo.count()
        self.data_combo.blockSignals(True)
        graphs = []
        graph = self.addTableFromEtree(table)
        graphs.append(graph)
        self.table_combo.setCurrentIndex(oldCount)
        self.data_combo.setCurrentIndex(0)
        self.graph.setCurrentIndex(oldCount)
        self.graph.currentWidget().setCurrentIndex(0)
        self.setCurrentData(0)
        self.table_combo.blockSignals(False)
        self.data_combo.blockSignals(False)
        return graphs

    def addTableFromDataIslandDivs(self,tree=None):
            print("addTableFromEtree"); sys.stdout.flush()
            print("addTableFromEtree",tree); sys.stdout.flush()
            print(etree.tostring(t,pretty_print=True))

    def addTableFromEtree(self,tree=None):
        t = tree
        ttitle = t.attrib.get('title','').strip()
        data = t.xpath('data')[0]
        a = data.text.strip().split()
        a2 = []
        for i in a:
            j = i.replace("NA","Nan")
            if j == "*" or j=="-":
                 j = "NaN"
            a2.append(j)
                
        array1 = numpy.array(a2,dtype=float)
        #jspimple trims by data row, not headers and so works, we'll try the same ...., but headers need proper separators also!
        ncols = len(data.text.strip().split("\n")[0].split())
        #array2 = numpy.array_split(array1,len(array1)/len(t.xpath('headers')[0].text.strip().split()))
        array2 = numpy.array_split(array1,len(array1)/ncols)
        theFirstArray = numpy.array(array2)
        graph_lists = []
        ydata = []
        headers = []
        for i,p in zip(list(range(len(t.xpath('plot')))),t.xpath('plot')):
            if len(p.xpath("title")) > 0:
                ptitle = p.xpath("title")[0].text.strip()
            else:
                ptitle = 'Plot %d' % (i+1)
            graph_lists.append((ptitle,""))
        self.graph_lists.append(graph_lists)
        self.table_combo.addItem(ttitle)
        self.comboItemsChanged()
        graph = StackedWidget()
        self.graph.addWidget(graph)
        self.table_etrees[graph] = t
        for p in t.xpath('plot'):
            def add_plot(p,t):
                if len(p.xpath("title")) > 0:
                    ptitle = p.xpath("title")[0].text.strip()
                else:
                    ptitle = ''
                plot = QtMatplotlibWidget(title=ptitle)
                plot.setDefaultSymbolSize(self.defaultSymbolSize)
                x_scaling = None
                y_scaling = None
                for obj in p.xpath('polygon|circle|line'):
                    if obj.tag == 'polygon':
                        poly = obj
                        array1 = numpy.array(poly.text.strip().split(),dtype=float)
                        array2 = numpy.hsplit(array1,len(array1)/2)
                        polygon = matplotlib.patches.Polygon(array2)
                        polygon.set_linewidth(float(poly.attrib.get('linesize','1')))
                        polygon.set_facecolor(poly.attrib.get('fillcolour','none'))
                        polygon.set_edgecolor(poly.attrib.get('linecolour','black'))
                        polygon.set_linestyle(poly.attrib.get('linestyle','solid'))
                        plot.canvas.ax[0].add_patch(polygon)
                    elif obj.tag == 'circle':
                        circ = obj
                        circle = matplotlib.patches.Circle((float(circ.attrib['xpos']),float(circ.attrib['ypos'])),float(circ.attrib['radius']))
                        circle.set_linewidth(float(circ.attrib.get('linesize','1')))
                        circle.set_facecolor(circ.attrib.get('fillcolour','none'))
                        circle.set_edgecolor(circ.attrib.get('linecolour','black'))
                        circle.set_linestyle(circ.attrib.get('linestyle','solid'))
                        plot.canvas.ax[0].add_patch(circle)
                    elif obj.tag == 'line':
                        l = obj
                        x1 = float(l.attrib['x1'])
                        x2 = float(l.attrib['x2'])
                        y1 = float(l.attrib['y1'])
                        y2 = float(l.attrib['y2'])
                        line = plot.canvas.ax[0].add_line(matplotlib.lines.Line2D([x1,x2],[y1,y2]))
                        line.set_linewidth(float(l.attrib.get('linesize','1')))
                        line.set_linestyle(l.attrib.get('linestyle','solid'))
                        line.set_color(l.attrib.get('linecolour','black'))
                considerFixedAspect = False
                if len(p.xpath("xrange")) > 0:
                    try:
                        try:
                            xmin = float(p.xpath("xrange")[0].attrib['min'])
                        except:
                            xmin = None
                        try:
                            xmax = float(p.xpath("xrange")[0].attrib['max'])
                        except:
                            xmax = None
                        plot.canvas.ax[0].set_xlim(xmin,xmax)
                        plot.canvas.custom_xlim = xmin,xmax
                    except:
                        pass
                if len(p.xpath("yrange")) > 0:
                    try:
                        iax = 0
                        for yrange in p.xpath("yrange"):
                            try:
                                ymin = float(yrange.attrib['min'])
                            except:
                                ymin = None
                            try:
                                ymax = float(yrange.attrib['max'])
                            except:
                                ymax = None
                            try:
                                if yrange.attrib['rightaxis'] == "true" or yrange.attrib['rightaxis'] == '1':
                                    rightaxis = True
                                else:
                                    rightaxis = False
                            except:
                                rightaxis = False
                            if not rightaxis:
                                plot.canvas.ax[iax].set_ylim(ymin,ymax)
                                plot.canvas.custom_ylim = ymin,ymax
                            else:
                                plot.canvas.rax.append(plot.canvas.ax[0].twinx())
                                plot.canvas.rax[-1].set_ylim(ymin,ymax)
                                plot.canvas.custom_rylim = ymin,ymax
                                try:
                                    rhcol = yrange.attrib['colour']
                                    plot.canvas.rax[-1].spines['right'].set_color(rhcol)
                                    plot.canvas.ax[-1].spines['right'].set_color(rhcol)
                                    for tl in plot.canvas.rax[-1].get_yticklabels():
                                            tl.set_color(rhcol)
                                except:
                                    pass
                    except:
                        pass

                    if plot.canvas.custom_ylim and len (plot.canvas.custom_ylim)==2 and (plot.canvas.custom_ylim[0] is not None) and (plot.canvas.custom_ylim[1] is not None) and plot.canvas.custom_xlim and len (plot.canvas.custom_xlim)==2 and (plot.canvas.custom_xlim[0] is not None) and (plot.canvas.custom_xlim[1] is not None):
                        considerFixedAspect = True
                    try:
                        if considerFixedAspect and plot.canvas.ax[iax].get_xlim()[0] is not None and plot.canvas.ax[iax].get_xlim()[1] is not None and plot.canvas.ax[iax].get_ylim()[0] is not None and plot.canvas.ax[iax].get_ylim()[1] is not None:
                            if abs((plot.canvas.ax[iax].get_xlim()[1]-plot.canvas.ax[iax].get_xlim()[0])-(plot.canvas.ax[iax].get_ylim()[1]-plot.canvas.ax[iax].get_ylim()[0]))<1e-5:
                                plot.canvas.fix_aspect_ratio = True
                                plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))
                                plot.canvas.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(plot.canvas.width(),plot.canvas.height()),QtCore.QSize(plot.canvas.width(),plot.canvas.height())))
                    except:
                        pass

                if len(p.xpath("legendposition")) > 0:
                    try:
                        xlegpos = float(p.xpath("legendposition")[0].attrib['x'])
                        ylegpos = float(p.xpath("legendposition")[0].attrib['y'])
                        plot.canvas.initial_legend_position = [xlegpos,ylegpos]
                    except:
                        pass
                if len(p.xpath('xscale')) > 0 :
                    x_scaling = scaling_functions[p.xpath("xscale")[0].text.strip()]
                    plot.canvas.x_scaling = x_scaling
                if len(p.xpath('yscale')) > 0 :
                    y_scaling = scaling_functions[p.xpath("yscale")[0].text.strip()]
                    plot.canvas.y_scaling = y_scaling
                if y_scaling:
                    plot.canvas.ax[0].yaxis.set_major_formatter(FuncFormatter(y_scaling))
                if x_scaling:
                    plot.canvas.ax[0].xaxis.set_major_formatter(FuncFormatter(x_scaling))
                if len(p.xpath('xlabel')) > 0 :
                    theText = plot.canvas.ax[0].set_xlabel(p.xpath("xlabel")[0].text.strip())
                    plot.canvas.xlabel = p.xpath("xlabel")[0].text.strip()
                if len(p.xpath('ylabel')) > 0 :
                    theText = plot.canvas.ax[0].set_ylabel(p.xpath("ylabel")[0].text.strip())
                    plot.canvas.ylabel = p.xpath("ylabel")[0].text.strip()
                if len(p.xpath('rylabel')) > 0 :
                    theText = plot.canvas.rax[0].set_ylabel(p.xpath("rylabel")[0].text.strip())
                    plot.canvas.rylabel = p.xpath("rylabel")[0].text.strip()
                    try:
                        theText.set_color(rhcol)
                    except:
                        pass
                if len(p.xpath('xbreaks')) > 0 :
                    xbreaks = p.xpath('xbreaks')[0]
                    breaks = []
                    for ibr in range(len(xbreaks.xpath('break'))):
                        theBreak = xbreaks.xpath('break')[ibr]
                        breaks.append((float(theBreak.attrib['min']), float(theBreak.attrib['max'])),)
                    breaks = tuple(breaks)
                    plot.canvas.set_breaks(breaks,None)
                elif len(p.xpath('ybreaks')) > 0 :
                    ybreaks = p.xpath('ybreaks')[0]
                    breaks = []
                    for ibr in range(len(ybreaks.xpath('break'))):
                        theBreak = ybreaks.xpath('break')[ibr]
                        breaks.append((float(theBreak.attrib['min']), float(theBreak.attrib['max'])),)
                    breaks = tuple(breaks)
                    plot.canvas.set_breaks(None,breaks)
                plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,x_scaling,y_scaling))
                plot.canvas.DetectPointMoveEvent.connect(functools.partial(self.updateStatusWithPickXandY,x_scaling,y_scaling))
                plot.canvas.DetectPointClickEvent.connect(graph.DetectPointClickEvent.emit)
                plot.canvas.DetectPointMoveEvent.connect(graph.DetectPointMoveEvent.emit)
                plot.canvas.RegionSelected.connect(graph.RegionSelected.emit)
                ig = 0
                theseColours = None; theColour = None
                theseSymbols = None; theSymbol = None
                theseSymbolSizes = None; theSymbolSize = None
                theseLineSizes = None; theLineSize = None
                theseLineStyles = None; theLineStyle = None
                for i in range(len(p.xpath('histogram'))):
                    histo = p.xpath('histogram')[i]
                    vals = []
                    col = int(histo.attrib['col'])
                    array = theFirstArray
                    vals = array[:,int(col)-1]
                    head = t.xpath('headers')[0].text.strip().split()[int(col)-1]
                    if len(histo.xpath('colour'))>0:
                        color = histo.xpath('colour')[0].text.strip()
                    else:
                        color = self.colors[ig % len(self.colors)]
                    if len(vals)>0:
                        mu = numpy.mean(vals)
                        sigma = numpy.std(vals)
                        maxx = numpy.max(vals)
                        #n, bins, patches = plot.canvas.ax.hist(vals, int(maxx), normed=0, facecolor=color, label=dataname+' '+str(ig))
                        #print mu, sigma, maxx, len(vals)
                        if maxx < 0.0:
                            maxx = -maxx
                        if abs(maxx) < 1.0:
                            maxx = 1.0
                        if len(histo.xpath('nbins'))>0:
                            maxx = histo.xpath('nbins')[0].text.strip()
                        n, bins, patches = plot.canvas.hist(vals, int(maxx), normed=0, facecolor=color, label=head)
                        #sys.exit()
                        if len(patches)>0:
                            # Maybe need a Histogram class ?
                            self.histograms.append({'normed':False,'vals':vals,'maxx':maxx,'patches':patches,'bins':bins,'mu':mu, 'sigma':sigma})
                for i in range(len(p.xpath('plotline'))):
                    plotline = p.xpath('plotline')[i]
                    ycol = int(plotline.attrib['ycol'])
                    array = theFirstArray
                    if not plotline.attrib.get('dataid',None) is None:
                        array = self.getDataArray(t,plotline.attrib['dataid'])
                        if array is None:
                            array = theFirstArray
                    ydat2 = array[:,int(ycol)-1]
                    head = t.xpath('headers')[0].text.strip().split()[int(ycol)-1]
                    vals = []
                    if len(plotline.xpath('colour'))>0:
                        color = plotline.xpath('colour')[0].text.strip()
                    else:
                        color = self.colors[ig % len(self.colors)]
                    if len(plotline.xpath('symbol'))>0:
                        style = plotline.xpath('symbol')[0].text.strip()
                    else:
                        style = self.styles[ig % len(self.styles)]
                    if len(plotline.xpath('symbolsize'))>0:
                        markersize = float(plotline.xpath('symbolsize')[0].text.strip())
                    else:
                        markersize = self.defaultSymbolSize
                    if len(plotline.xpath('linestyle'))>0:
                        linestyle = plotline.xpath('linestyle')[0].text.strip()
                    else:
                        linestyle = '-'
                    if len(plotline.xpath('linesize'))>0:
                        linewidth = int(plotline.xpath('linesize')[0].text.strip())
                    else:
                        linewidth = 1
                    if not color in self.colors:
                        try:
                            c1 = matplotlib.colors.ColorConverter().to_rgb(color)
                            for cn in self.colors:
                                c2 = matplotlib.colors.ColorConverter().to_rgb(cn)
                                if c2 == c1:
                                    break
                            else:
                                self.colors.append(color)
                                self.colors_alias.append(color.capitalize())
                        except:
                            color = 'blue'
                            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                            sys.stderr.write(str(exc_type)+'\n')
                            sys.stderr.write(str(exc_value)+'\n')
                            idx = 0
                    xcol = int(plotline.attrib['xcol'])-1
                    xdata = array[:,xcol]

                    if len(p.xpath('xintegral'))>0:
                        if p.xpath('xintegral')[0].text == "true" or p.xpath('xintegral')[0].text == '1':
                            plot.canvas.xintegral=True
                    if len(p.xpath('yintegral'))>0:
                        if p.xpath('yintegral')[0].text == "true" or p.xpath('yintegral')[0].text == '1':
                            plot.canvas.yintegral=True

                    visible=True
                    if len(plotline.xpath('visible'))>0:
                        if plotline.xpath('visible')[0].text == "false" or plotline.xpath('visible')[0].text == '0':
                            visible=False
                    rightaxis = False
                    try:
                        if plotline.attrib['rightaxis'] == "true" or yrange.attrib['rightaxis'] == '1':
                            rightaxis = True
                            if len(plot.canvas.rax)==0:
                                plot.canvas.rax.append(plot.canvas.ax[0].twinx())
                                plot.canvas.rax[-1].set_ylim(None,None)
                                plot.canvas.custom_rylim = None,None
                    except:
                        pass
                    if linestyle == '.':
                        thePlot = plot.canvas.plot(xdata,ydat2,label=head,picker=self.defaultSymbolSize, marker=style, color=color, markersize=markersize,linestyle='None',linewidth=linewidth,visible=visible,rightaxis=rightaxis)
                    else:
                        thePlot = plot.canvas.plot(xdata,ydat2,label=head,picker=self.defaultSymbolSize, marker=style, color=color, markersize=markersize,linestyle=linestyle,linewidth=linewidth,visible=visible,rightaxis=rightaxis)
                    self.elementToPlotMap[plotline] = thePlot
                    if len(plotline.xpath('showinlegend'))>0 and len(thePlot)>0:
                        if plotline.xpath('showinlegend')[0].text == "false" or plotline.xpath('showinlegend')[0].text == '0':
                            thePlot[0].hide_from_legend = True
                    if len(plotline.xpath('label'))>0 and len(thePlot)>0:
                        thePlot[0].set_label(plotline.xpath('label')[0].text)

                    if len(plotline.xpath('markeredgewidth'))>0 and len(thePlot)>0:
                        thePlot[0].set_markeredgewidth(float(plotline.xpath('markeredgewidth')[0].text))
                    
                    ig = ig+1

                if len(p.xpath("legendposition")) == 0:
                    halfx = (plot.canvas.ax[0].get_xlim()[0]+plot.canvas.ax[0].get_xlim()[1])/2.
                    halfy = (plot.canvas.ax[0].get_ylim()[0]+plot.canvas.ax[0].get_ylim()[1])/2.
                    nTopLeft = nTopRight = nBottomLeft = nBottomRight = 0
                    for i in range(len(p.xpath('plotline'))):
                        plotline = p.xpath('plotline')[i]
                        ycol = int(plotline.attrib['ycol'])
                        ydat2 = array[:,int(ycol)-1]
                        xcol = int(plotline.attrib['xcol'])-1
                        xdata = array[:,xcol]
                        nTopLeft_i, nTopRight_i, nBottomLeft_i, nBottomRight_i = self.analyzeDistributionOfPoints(xdata, ydat2, halfx, halfy)
                        nTopLeft += nTopLeft_i
                        nBottomLeft += nBottomLeft_i
                        nTopRight += nTopRight_i
                        nBottomRight += nBottomRight_i
                    legend_pos_idx = (nTopLeft, nTopRight, nBottomLeft, nBottomRight).index(min(nTopLeft, nTopRight, nBottomLeft, nBottomRight))
                    if legend_pos_idx == 0:
                        plot.canvas.initial_legend_position = [0.15,0.75]
                    elif legend_pos_idx == 1:
                        plot.canvas.initial_legend_position = [0.75,0.75]
                    elif legend_pos_idx == 2:
                        plot.canvas.initial_legend_position = [0.15,0.15]
                    elif legend_pos_idx == 3:
                        plot.canvas.initial_legend_position = [0.75,0.15]
                    
                plot.canvas.legend()


                try:
                    plot.canvas.applyPreferences(self.titleFontSel,self.axesFontSel,self.axesLabelFontSel,self.legendFontSel,redraw=False)
                except:
                    sys.stderr.write("Failed to apply font preferences\n")
                    exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                    sys.stderr.write(str(exc_type)+'\n')
                    sys.stderr.write(str(exc_value)+'\n')
                if len(p.xpath('fixaspectratio'))>0:
                    fixaspectratio = p.xpath('fixaspectratio')[0].text.strip()
                    if fixaspectratio =='1' or fixaspectratio =='true':
                        plot.canvas.fix_aspect_ratio = True
                if len(p.xpath('showlegend'))>0:
                    showlegend = p.xpath('showlegend')[0].text.strip()
                    if showlegend =='0' or showlegend =='false':
                        plot.canvas.ax[0].legend_.set_visible(False)
                return plot
            
            #plot = add_plot(p,t)
            #graph.addWidget(plot)
            nullWidget = QtWidgets.QWidget()
            nullWidget._do_the_plot_ = functools.partial(add_plot,p,t)
            graph.addWidget(nullWidget)

            #self.table_combo.setCurrentIndex(self.table_combo.count()-1)
        self.statusFormat = None
        return graph

    def addCCP4Table(self,table,graph_uuid=None):

        self.graph_lists.append(table.graphs_list)
        self.table_combo.addItem(table.title)
        self.comboItemsChanged()
        graph = StackedWidget()
        if graph_uuid != None:
            graph.setObjectName(graph_uuid)
        self.graph.addWidget(graph)
        t = etree.Element("CCP4Table")
        theData = etree.Element("data")
        theDataText = ''
        for tdl in table.data_lines:
            theDataText += ''.join([("%s "%n) for n in tdl[:]]).strip() + '\n'
        theData.text = theDataText
        t.append(theData)
        self.table_etrees[graph] = t
        theHeaders = etree.Element("headers",separator=' ')
        theHeaders.text = ''.join([("%s "%n) for n in table.headers[:]]).strip() + '\n'
        t.append(theHeaders)
        for g in table.graphs_list:
            p = etree.Element('plot')
            theTitle = etree.Element("title")
            theTitle.text = g[0]
            p.append(theTitle)
            t.append(p)
            for s in g[1].split(',')[1:]:
                plotline = etree.Element('plotline', xcol=g[1].split(',')[0], ycol=s)
                p.append(plotline)
            def add_plot(table,graph,g,t,p):
                plot = QtMatplotlibWidget(title=g[0])
                plot.setDefaultSymbolSize(self.defaultSymbolSize)
                plot.canvas.MouseMoveEvent.connect(functools.partial(self.updateStatusWithXandY,table.x_scaling,table.y_scaling))
                plot.canvas.DetectPointMoveEvent.connect(functools.partial(self.updateStatusWithPickXandY,table.x_scaling,table.y_scaling))
                plot.canvas.DetectPointClickEvent.connect(graph.DetectPointClickEvent.emit)
                plot.canvas.DetectPointMoveEvent.connect(graph.DetectPointMoveEvent.emit)
                plot.canvas.RegionSelected.connect(graph.RegionSelected.emit)
                for s in g[1].split(',')[0:1]:
                    xvals = []
                    col = int(s)-1
                    if table.custom_x_label:
                        plot.canvas.ax[0].set_xlabel(table.custom_x_label)
                        plot.canvas.xlabel = table.custom_x_label
                    else:
                        plot.canvas.ax[0].set_xlabel(table.headers[col])
                        plot.canvas.xlabel = table.headers[col]
                    for l in table.data_lines:
                        try:
                            v = l[col]
                            xvals.append(v)
                        except:
                            print("Error",l,col)
                ig = 0
                theRange = None
                if len(g)>2:
                    theRange = g[2]

                theXVals = []
                theYVals = []
                for s in g[1].split(',')[1:]:
                    plotline = p.xpath('plotline')[ig]
                    vals = []
                    col = int(s)-1
                    myxvals = []
                    for l,il in zip(table.data_lines,list(range(len(table.data_lines)))):
                        try:
                            f = float(l[col])
                            vals.append(l[col])
                            myxvals.append(xvals[il])
                        except:
                            pass

                    color = self.colors[ig % len(self.colors)]
                    style = self.styles[ig % len(self.styles)]
                    try:
                        thePlot = plot.canvas.plot(myxvals,vals,label=table.headers[col],picker=self.defaultSymbolSize, marker=style, color=color, markersize=self.defaultSymbolSize)
                        theXVals.append(myxvals)
                        theYVals.append(vals)
                        self.elementToPlotMap[plotline] = thePlot
                        try:
                            if theRange and len(theRange)==2 and theRange[0] and theRange[1] and len(theRange[0])==2 and len(theRange[1])==2:
                                xmin,xmax = theRange[0] 
                                ymin,ymax = theRange[1] 
                                plot.canvas.ax[0].set_xlim(xmin,xmax)
                                plot.canvas.ax[0].set_ylim(ymin,ymax)
                                plot.canvas.custom_xlim = xmin,xmax
                                plot.canvas.custom_ylim = ymin,ymax
                        except:
                            pass
                    except:
                        print("Badly formatted line",l,vals, len(myxvals), len(vals))
                    ig = ig+1

                if table.y_scaling:
                    plot.canvas.ax[0].yaxis.set_major_formatter(FuncFormatter(table.y_scaling))
                    plot.canvas.y_scaling = table.y_scaling
                if table.x_scaling:
                    plot.canvas.ax[0].xaxis.set_major_formatter(FuncFormatter(table.x_scaling))
                    plot.canvas.x_scaling = table.x_scaling

                halfx = (plot.canvas.ax[0].get_xlim()[0]+plot.canvas.ax[0].get_xlim()[1])/2.
                halfy = (plot.canvas.ax[0].get_ylim()[0]+plot.canvas.ax[0].get_ylim()[1])/2.
                nTopLeft = nTopRight = nBottomLeft = nBottomRight = 0
                for i in range(len(theXVals)):
                    nTopLeft_i, nTopRight_i, nBottomLeft_i, nBottomRight_i = self.analyzeDistributionOfPoints(theXVals[i], theYVals[i], halfx, halfy)
                    nTopLeft += nTopLeft_i
                    nBottomLeft += nBottomLeft_i
                    nTopRight += nTopRight_i
                    nBottomRight += nBottomRight_i
                legend_pos_idx = (nTopLeft, nTopRight, nBottomLeft, nBottomRight).index(min(nTopLeft, nTopRight, nBottomLeft, nBottomRight))
                if legend_pos_idx == 0:
                    plot.canvas.initial_legend_position = [0.15,0.75]
                elif legend_pos_idx == 1:
                    plot.canvas.initial_legend_position = [0.75,0.75]
                elif legend_pos_idx == 2:
                    plot.canvas.initial_legend_position = [0.15,0.15]
                elif legend_pos_idx == 3:
                    plot.canvas.initial_legend_position = [0.75,0.15]

                plot.canvas.legend()

                try:
                    plot.canvas.applyPreferences(self.titleFontSel,self.axesFontSel,self.axesLabelFontSel,self.legendFontSel,redraw=False)
                except:
                    sys.stderr.write("Failed to apply font preferences\n")
                    exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                    sys.stderr.write(str(exc_type)+'\n')
                    sys.stderr.write(str(exc_value)+'\n')
                return plot

            #graph.addWidget(add_plot(table,graph,g))
            
            nullWidget = QtWidgets.QWidget()
            nullWidget._do_the_plot_ = functools.partial(add_plot,table,graph,g,t,p)
            graph.addWidget(nullWidget)

        #self.table_combo.setCurrentIndex(self.table_combo.count()-1)
        return graph

    def sizeHint(self):
        return QtCore.QSize(800,600)

    def comboItemsChanged(self):
        if self.table_combo.count()>1 or self.data_combo.count()>1:
            self.table_combo.show()
            self.data_combo.show()
        else:
            self.table_combo.hide()
            self.data_combo.hide()

    def __init__(self,parent=None):
        QtWidgets.QWidget.__init__(self,parent)
        #print 'Loggraph parent',self.parent()
        self.titleFontSel,self.axesFontSel,self.axesLabelFontSel,self.legendFontSel = None,None,None,None
        layout = QtWidgets.QHBoxLayout()
        mainLayout = QtWidgets.QVBoxLayout()
        self.setLayout(layout)
        layout.addLayout(mainLayout)
        self.graph = StackedWidget()
        nullGraph = QtMatplotlibWidget()
        #self.graph.addWidget(nullGraph)
        self.table_combo = QtWidgets.QComboBox()
        self.data_combo = QtWidgets.QComboBox()
        self.table_combo.hide()
        self.data_combo.hide()
        comboLayout = QtWidgets.QHBoxLayout()
        comboLayout.addWidget(self.table_combo)
        comboLayout.addWidget(self.data_combo)
        mainLayout.addWidget(self.graph)
        mainLayout.addLayout(comboLayout)
        self.graph_lists = []
        self.histograms = []
        self.table_etrees = {}
        self.surface_etrees = []
        self.elementToPlotMap = {}
        line = QtWidgets.QHBoxLayout()
        line.setContentsMargins(0,0,0,0)
        line.setSpacing(0)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        comboLayout.setContentsMargins(0,0,0,0)
        comboLayout.setSpacing(0)
        mainLayout.setContentsMargins(0,0,0,0)
        mainLayout.setSpacing(0)
        self.statusBar=QtWidgets.QStatusBar()
        self.statusBar.setSizeGripEnabled(False)
        self.statusFormat = None
        # Add spacing to statusBar line to get nicer layout if add a widget in ccp4i2
        line.addWidget(self.statusBar)
        line.addSpacing(10)
        mainLayout.addLayout(line)
        self.table_combo.currentIndexChanged[int].connect(self.setCurrentTable)
        self.data_combo.currentIndexChanged[int].connect(self.setCurrentData)
        self.setWindowTitle("Pimple")
        self.fileOpen = QtWidgets.QAction("Open..",self)
        self.fileSave = QtWidgets.QAction("Save figure..",self)
        self.fileSaveStatus = QtWidgets.QAction("Save status..",self)
        self.fileSaveAll = QtWidgets.QAction("Save all figures..",self)
        self.prefAction = QtWidgets.QAction("Preferences..",self)
        self.editLegendPos = QtWidgets.QAction("Edit legend position",self)
        self.editPlotStyle = QtWidgets.QAction("Edit plot style",self)
        haveEditIcon = True
        try:
            icon = QtGui.QIcon()
            if getattr(sys, 'frozen', None):
                basedir = sys._MEIPASS
                fbase = os.path.normpath(os.path.join(basedir,"qticons","actions","edit"))
            else:
                fbase = os.path.normpath(os.path.join(os.path.dirname(__file__),"..","qticons","actions","edit"))
            pngs = glob.glob(os.path.join(fbase,"*.png"))
            if len(pngs)==0:
                haveEditIcon = False
            for png in pngs:
                if sys.version_info >= (3,0):
                    icon.addFile(png)
                else:
                    icon.addFile(unicode(png,'utf-8'))
            self.editPlotStyle.setIcon(icon)
        except:
            haveEditIcon = False
            print("Failed to set edit icon")
        self.icons_missing = False
        self.fileOpen.setShortcut(QtGui.QKeySequence.Open)
        self.fileSaveStatus.setShortcut(QtGui.QKeySequence.Save)
        @QtCore.Slot(bool)
        def saveTriggered(dum):
            self.saveFileDialog()
        @QtCore.Slot(bool)
        def saveAllTriggered(dum):
            self.saveFileDialog(True,False)
        @QtCore.Slot(bool)
        def saveStatusTriggered(dum):
            self.saveFileDialog(False,True)
        self.fileSave.triggered.connect(saveTriggered)
        self.fileSaveAll.triggered.connect(saveAllTriggered)
        self.fileSaveStatus.triggered.connect(saveStatusTriggered)
        self.fileOpen.triggered.connect(self.loadFileDialog)
        self.editLegendPos.triggered.connect(self.editLegendPosition)
        self.editPlotStyle.triggered.connect(self.editPlotPlotStyle)
        self.prefAction.triggered.connect(self.editPreferences)
        if getattr(sys, 'frozen', None):
            basedir = sys._MEIPASS
            icondir = os.path.abspath(os.path.join(basedir,'qticons','loggraph'))
        else:
            icondir = os.path.abspath(os.path.join(os.path.dirname(__file__),'..','qticons','loggraph'))
        if os.path.exists(icondir):
            self.fileOpen.setIcon(QtGui.QIcon(os.path.join(icondir,"fileopen.svg")))
            self.fileSave.setIcon(QtGui.QIcon(os.path.join(icondir,"save_picture.svg")))
            self.prefAction.setIcon(QtGui.QIcon(os.path.join(icondir,"preferences.svg")))
            self.editLegendPos.setIcon(QtGui.QIcon(os.path.join(icondir,"moving.svg")))
            self.fileSaveStatus.setIcon(QtGui.QIcon(os.path.join(icondir,"save_status.svg")))
        self.actionBar = QtWidgets.QWidget()
        actionLayout = QtWidgets.QVBoxLayout()
        openButton = QtWidgets.QToolButton()
        openButton.setDefaultAction(self.fileOpen)
        if not os.path.exists(os.path.join(icondir,"fileopen.svg")):
            openButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            self.icons_missing = True
        actionLayout.addWidget(openButton)
        saveStatusButton = QtWidgets.QToolButton()
        saveStatusButton.setDefaultAction(self.fileSaveStatus)
        if not os.path.exists(os.path.join(icondir,"save_status.svg")):
            saveStatusButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            self.icons_missing = True
        saveButton = QtWidgets.QToolButton()
        saveButton.setDefaultAction(self.fileSave)
        if not os.path.exists(os.path.join(icondir,"save_picture.svg")):
            saveButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            self.icons_missing = True
        actionLayout.addWidget(saveButton)
        editLegendButton = QtWidgets.QToolButton()
        editLegendButton.setDefaultAction(self.editLegendPos)
        if not os.path.exists(os.path.join(icondir,"moving.svg")):
            editLegendButton.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            self.icons_missing = True
        actionLayout.addWidget(editLegendButton)
        editPlotStyle = QtWidgets.QToolButton()
        editPlotStyle.setDefaultAction(self.editPlotStyle)
        if not haveEditIcon:
            editPlotStyle.setToolButtonStyle(QtCore.Qt.ToolButtonTextOnly)
            self.icons_missing = True
        actionLayout.addWidget(editPlotStyle)
        self.actionBar.setLayout(actionLayout)
        actionLayout.addStretch(2)
        layout.addWidget(self.actionBar)
        actionLayout.setContentsMargins(3,3,3,3)
        self.actionBar.hide()
        self.defaultSymbolSize = 5

    def showActionBar(self):
        self.actionBar.show()

    def setDefaultSymbolSize(self,size=5):
        self.defaultSymbolSize = size
        
class  QtMatplotlibWidget(QtWidgets.QWidget):

    def __init__(self,parent=None,figsize = (800, 600), dpi=100,facecolor = '#FFFFFF',title='',x_axis_label='',y_axis_label=''):
        QtWidgets.QWidget.__init__(self,parent)
        fig = Figure(figsize=figsize,dpi=dpi,facecolor=facecolor)
        self.canvas = QtMatplotlibCanvas(fig,title=title,x_axis_label=x_axis_label,y_axis_label=y_axis_label)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)
        self.setCursor(QtCore.Qt.CrossCursor)
        layout.setContentsMargins(0,0,0,0)
        self.defaultSymbolSize = 5

    def addPoint(self,point):
        self.canvas.addPoint(point)

    def setDefaultSymbolSize(self,size=5):
        self.defaultSymbolSize = size
        self.canvas.setDefaultSymbolSize(size)

class QtMatplotlibCanvas(FigureCanvas):

    MouseMoveEvent = QtCore.Signal(tuple)
    DetectPointMoveEvent = QtCore.Signal(dict)
    DetectPointClickEvent = QtCore.Signal(dict)
    RegionSelected = QtCore.Signal(dict)

    def setDefaultSymbolSize(self,size=5):
        self.defaultSymbolSize = size

    def applyPreferences(self,titleFontSel=None,axesFontSel=None,axesLabelFontSel=None,legendFontSel=None,redraw=True):
        try:
            params = legendFontSel
            if params:
                if 'family' in params:
                    family = params['family']
                else:
                    family = "Bitstream Vera Sans" 
                if 'size' in params:
                    size = params['size']
                else:
                    size = 8
                newFont = FontProperties(family=family,size=size)
                self.legend_font = newFont
            else:
                newFont = FontProperties(family="Bitstream Vera Sans",size=8)
                self.legend_font = newFont
        except:
            print("error 1")
            pass
        #print 'applyPreferences titleFontSel',titleFontSel
        try:
            params = titleFontSel
            if params:
                if 'family' in params:
                    family = params['family']
                else:
                    family = "Bitstream Vera Sans" 
                if 'size' in params:
                    size = params['size']
                else:
                    size = 8
                newFont = FontProperties(family=family,size=size)
                self.title_font = newFont
            else:
                newFont = FontProperties(family="Bitstream Vera Sans",size=8)
                self.title_font = newFont
        except:
            print("error 2")
            pass
        try:
            params = axesFontSel
            if params:
                if 'family' in params:
                    family = params['family']
                else:
                    family = "Bitstream Vera Sans" 
                if 'size' in params:
                    size = params['size']
                else:
                    size = 8
                newFont = FontProperties(family=family,size=size)
                self.axes_font = newFont
            else:
                newFont = FontProperties(family="Bitstream Vera Sans",size=8)
                self.axes_font = newFont
        except:
            print("error 3")
            pass
        try:
            params = axesLabelFontSel
            if params:
                if 'family' in params:
                    family = params['family']
                else:
                    family = "Bitstream Vera Sans" 
                if 'size' in params:
                    size = params['size']
                else:
                    size = 8
                newFont = FontProperties(family=family,size=size)
                self.label_font = newFont
            else:
                newFont = FontProperties(family="Bitstream Vera Sans",size=8)
                self.label_font = newFont
        except:
            print("error 4")
            pass

        for fn,params in ((self.legend_font,legendFontSel),(self.title_font,titleFontSel),(self.axes_font,axesFontSel),(self.label_font,axesLabelFontSel)):
            if params:
                if 'weight' in params and params['weight']=='bold':
                    fn.set_weight('bold')
                else:
                    fn.set_weight('normal')
                if 'slant' in params and params['slant']=='italic':
                    fn.set_slant('italic')
                else:
                    fn.set_slant('normal')

        self.legend()

        if self.xbreak:
            nbreaks = len(self.xbreak)+1
            theTitle = self.ax[int(nbreaks/2.)].set_title(self.title)
            theXLabel = self.ax[int(nbreaks/2.)].set_xlabel(self.xlabel)
            theYLabel = self.ax[0].set_ylabel(self.ylabel)
            if len(self.rax):
                theRYLabel = self.rax[-1].set_ylabel(self.rylabel)
                theRYLabel.set_fontproperties(self.label_font)
        elif self.ybreak:
            nbreaks = len(self.ybreak)+1
            theTitle = self.ax[0].set_title(self.title)
            theXLabel = self.ax[-1].set_xlabel(self.xlabel)
            theYLabel = self.ax[int(nbreaks/2.)].set_ylabel(self.ylabel)
        else:
            theTitle = self.ax[0].set_title(self.title)
            theXLabel = self.ax[0].set_xlabel(self.xlabel)
            theYLabel = self.ax[0].set_ylabel(self.ylabel)

        theTitle.set_fontproperties(self.title_font)
        theXLabel.set_fontproperties(self.label_font)
        theYLabel.set_fontproperties(self.label_font)

        self.set_tick_labels_font()

        if not redraw:
            return

        self.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(self.width(),self.height()),QtCore.QSize(self.width(),self.height())))
        self.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(self.width(),self.height()),QtCore.QSize(self.width(),self.height())))
        self.draw()

    def legend(self,loc=None,bbox_to_anchor=None):
        xratio = 1.0
        if hasattr(self,"old_xbreak"):
            if (self.old_xbreak and not self.xbreak) or (not self.old_xbreak and self.xbreak) or (self.old_xbreak and self.xbreak and len(self.old_xbreak) != len(self.xbreak)):
                if self.xbreak:
                    nax = len(self.xbreak) + 1
                else:
                    nax = 1
                if self.old_xbreak:
                    naxold = len(self.old_xbreak) + 1
                else:
                    naxold = 1
                xratio = float(nax)/float(naxold)

        haveLegend = False
        handles,labels = self.ax[0].get_legend_handles_labels()
        if len(self.rax)>0:
            rhandles,rlabels = self.rax[0].get_legend_handles_labels()
            handles.extend(rhandles)
            labels.extend(rlabels)
        leg_handles = []
        leg_labels = []
        for handle, label in zip(handles,labels):
            if not label.endswith('_matplotlib_line_of_best_fit'):
                if (hasattr(handle,"hide_from_legend") and handle.hide_from_legend == True) or not handle.get_visible():
                    pass
                else:
                    leg_handles.append(handle)    
                    leg_labels.append(label)    
                    haveLegend = True
        if not haveLegend:
            if self.ax[0].get_legend():
                self.ax[0].legend_.set_visible(False)
        else:
            if loc and bbox_to_anchor.all():
                if self.xbreak:
                    self.custom_legend_position =  bbox_to_anchor[0]/(xratio*1.3*(len(self.xbreak)+1)),bbox_to_anchor[1]
                self.ax[0].legend(leg_handles,leg_labels,fancybox=True, shadow=True, loc = loc, bbox_to_anchor = bbox_to_anchor, prop=self.legend_font)
            elif self.ax[0].get_legend():
                if(hasattr,self.ax[0].get_legend().get_bbox_to_anchor().get_points(),"item") and isinstance(self.ax[0].get_legend().get_bbox_to_anchor().get_points().item, Callable):
                    points = self.ax[0].get_legend().get_bbox_to_anchor().get_points().item(0), self.ax[0].get_legend().get_bbox_to_anchor().get_points().item(1)
                else:
                    points = self.ax[0].get_legend().get_bbox_to_anchor().get_points()[1]
                loc = self.ax[0].get_legend().parent.transAxes.inverted().transform_point((xratio*points[0],points[1]))
                fm = QtGui.QFontMetrics(QtGui.QFont(self.legend_font.get_family()[0],self.legend_font.get_size()))
                maxp = 0
                for lab in leg_labels:
                    maxp = max(maxp,fm.width(lab))
                leg_h = fm.height()*len(leg_labels)
                if sys.platform == "darwin":
                    #No idea why I have to do this hack.
                    leg_h *= 2
                    maxp *= 2
                if loc[0] + float(maxp)/self.parent().width() > 0.8:
                    loc[0] -= float(maxp)/self.parent().width()
                    #print "Adjusting x pos"
                if loc[1] + float(leg_h)/self.parent().height()>0.8:
                    if sys.platform == "darwin":
                        loc[1] -= float(leg_h)/self.parent().height()
                    else:
                        loc[1] -= 1.5*float(leg_h)/self.parent().height()
                    #print "Adjusting y pos"
                self.custom_legend_position =  loc[0],loc[1]
                if self.xbreak:
                    self.custom_legend_position =  loc[0]/(xratio*1.3*(len(self.xbreak)+1)),loc[1]

                self.ax[0].legend(leg_handles,leg_labels,fancybox=True, shadow=True, prop=self.legend_font, loc ='lower left', bbox_to_anchor =loc)
            else:
                self.ax[0].legend(leg_handles,leg_labels,fancybox=True, shadow=True, prop=self.legend_font, loc ='lower left', bbox_to_anchor =self.initial_legend_position)
                if self.xbreak:
                    loc = (xratio*1.0*(len(self.xbreak)+1)*self.initial_legend_position[0],self.initial_legend_position[1])
                    self.ax[0].legend(leg_handles,leg_labels,fancybox=True, shadow=True, prop=self.legend_font, loc = 'lower left', bbox_to_anchor =loc)
        self.old_xbreak = self.xbreak
        self.old_ybreak = self.ybreak

    def keyPressEvent(self,e):
        if e.key() == QtCore.Qt.Key_Escape:
            self.movingLegend = False
            self.repaint()

    def mousePressEvent(self,e):
        self.mouseDownPos = QtCore.QPoint(e.pos())
        FigureCanvas.mousePressEvent(self,e)

    def mouseReleaseEvent(self, e):
        self.mouseUpPos = QtCore.QPoint(e.pos())
        FigureCanvas.mouseReleaseEvent(self,e)
        self.repaint()
        if not hasattr(self,"mouseDownPos"):
            return
        if self.mouseDownPos.x() != self.mouseUpPos.x():
            # FIXME, this needs to cope with broken axes! Hard.
            xlim = self.ax[0].get_xlim()
            figWidth = self.width()*(float(self.right)-self.left)
            xp = (float(self.mouseDownPos.x())-float(self.left)*self.width())/figWidth
            x1 = xp*(xlim[1]-xlim[0])+xlim[0]
            xp = (float(self.mouseUpPos.x())-float(self.left)*self.width())/figWidth
            x2 = xp*(xlim[1]-xlim[0])+xlim[0]
            self.RegionSelected.emit({'plotwidget':self.parent(), 'x1':x1,'x2':x2})

    def mouseMoveEvent(self,e):
        if e.buttons() & QtCore.Qt.LeftButton:
            self.mouseUpPos = QtCore.QPoint(e.pos())
        FigureCanvas.mouseMoveEvent(self,e)
        self.repaint()

    def resizeEvent(self,e):
        if self._type == 'surface':
            return FigureCanvas.resizeEvent(self,e)
        if self.width() == 0:
            self.left = 0.3
            self.right = 1
        else:
            if self.ylabel or self.xlabel:
                if sys.version_info >= (3,0) or type(self.label_font.get_family()[0]) == unicode:
                        labelFontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(self.label_font.get_family()[0],self.label_font.get_size())).height())/self.defaultTextHeight,0.6)
                else:
                        labelFontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(unicode(self.label_font.get_family()[0],"utf-8"),self.label_font.get_size())).height())/self.defaultTextHeight,0.6)
            if self.ylabel:
                self.left = 0.3 / self.width()*200.*labelFontscale
            else:
                self.left = 0.2 / self.width()*200.
            self.right = 1 - 0.1 / self.width()*200.
            if hasattr(self,"rylabel") and self.rylabel:
                self.right -= 0.1 / self.width()*200.*labelFontscale
        if self.height() == 0:
            self.top = 0.95
            self.bottom = 0.1
        else:
            if self.title:
                if sys.version_info >= (3,0) or type(self.label_font.get_family()[0]) == unicode:
                        fontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(self.title_font.get_family()[0],self.title_font.get_size())).height())/self.defaultTextHeight,0.6)
                else:
                        fontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(unicode(self.title_font.get_family()[0],"utf-8"),self.title_font.get_size())).height())/self.defaultTextHeight,0.6)
                self.top = 1 - 0.15 / self.height()*200.*fontscale
            else:
                self.top = 1 - 0.05 / self.height()*200.
            if self.xlabel:
                self.bottom = 0.2 / self.height()*200.*labelFontscale
            else:
                self.bottom = 0.1 / self.height()*200.

        if sys.version_info >= (3,0) or type(self.label_font.get_family()[0]) == unicode:
            axesFontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(self.axes_font.get_family()[0],self.axes_font.get_size())).height())/self.defaultTextHeight,0.6)-1.
        else:
            axesFontscale = math.pow(float(QtGui.QFontMetrics(QtGui.QFont(unicode(self.axes_font.get_family()[0],"utf-8"),self.axes_font.get_size())).height())/self.defaultTextHeight,0.6)-1.
        self.left += 0.3 / self.width()*200.*axesFontscale
        self.bottom += 0.3 / self.width()*200.*axesFontscale
        if len(self.rax)>0:
            self.right -= 50*(self.axes_font.get_size()/self.defaultTextHeight)/self.width()
        if self.fix_aspect_ratio and self.width()>0 and self.height()>0:
            if self.width()>self.height():
                diff = (1-self.right+self.left)/2*self.height()/self.width()
                ratio = (float(self.width())/self.height()-1)/2/(float(self.width())/self.height())
                self.left = ratio
                self.right = 1-ratio
                self.left += (1-self.top+self.bottom)/2*self.height()/self.width()
                self.right -= (1-self.top+self.bottom)/2*self.height()/self.width()
                if self.ylabel:
                    self.left += diff
                    self.bottom += diff
            else:
                ratio = (float(self.height())/self.width()-1)/2/(float(self.height())/self.width())
                self.bottom = ratio
                self.top = 1-ratio
                self.bottom += (1-self.right+self.left)/2*self.width()/self.height()
                self.top -= (1-self.right+self.left)/2*self.width()/self.height()

        if self.left >= self.right:
            self.left = self.right - 0.001
        if self.bottom >= self.top:
            self.bottom = self.top - 0.001
        self.fig.subplots_adjust(self.left, self.bottom, self.right, self.top)
        if self.xbreak:
            nbreaks = len(self.xbreak)+1
            self.gs[0].update(left=self.left, right=(self.left+self.right)/float(nbreaks)-0.5*self.breaksize, top=self.top, bottom=self.bottom, wspace=0.0)
            for i in range(1,len(self.gs)):
                self.gs[i].update(left=i*(self.left+self.right)/float(nbreaks)+0.5*self.breaksize, right=min(self.right,(i+1)*(self.left+self.right)/float(nbreaks)-0.5*self.breaksize), top=self.top, bottom=self.bottom, wspace=0.0)
        elif self.ybreak:
            nbreaks = len(self.ybreak)+1
            bottom = (nbreaks-1)*(self.top+self.bottom)/float(nbreaks)+0.5*self.breaksize
            if bottom >= self.top:
                bottom = self.top - 0.001
            self.gs[0].update(left=self.left, right=self.right, top=self.top, bottom=bottom, wspace=0.0)
            for i in range(1,len(self.gs)):
                bottom = max(self.bottom,(i-1)*(self.top+self.bottom)/float(nbreaks)+0.5*self.breaksize)
                top = i*(self.top+self.bottom)/float(nbreaks)-0.5*self.breaksize
                if bottom >=top:
                    bottom = top - 0.001
                self.gs[nbreaks-i].update(left=self.left, right=self.right, top=top, bottom=bottom, wspace=0.0)
        else:
            self.gs[0].update(left=self.left, right=self.right, top=self.top, bottom=self.bottom, wspace=0.0)
            
        FigureCanvas.resizeEvent(self,e)

    def paintEvent(self,e):
        FigureCanvas.paintEvent(self,e)
        painter = QtGui.QPainter(self)
        painter.setPen(QtCore.Qt.blue)
        pos = self.mapFromGlobal(self.cursor().pos())
        if hasattr(self,"mouseDownPos") and hasattr(self,"mouseUpPos"):
            if self.mouseDownPos.x() != self.mouseUpPos.x():
                if self.mouseDownPos.x() < self.mouseUpPos.x():
                    painter.fillRect(self.mouseDownPos.x(),(1-self.top)*self.height(),self.mouseUpPos.x()-self.mouseDownPos.x(),(1-self.bottom)*self.height()-(1-self.top)*self.height(),QtGui.QBrush(QtGui.QColor(128,128,255,60),QtCore.Qt.SolidPattern))
                else:
                    painter.fillRect(self.mouseUpPos.x(),(1-self.top)*self.height(),self.mouseDownPos.x()-self.mouseUpPos.x(),(1-self.bottom)*self.height()-(1-self.top)*self.height(),QtGui.QBrush(QtGui.QColor(128,128,255,60),QtCore.Qt.SolidPattern))
    
                painter.drawLine(self.mouseDownPos.x(),(1-self.top)*self.height(),self.mouseDownPos.x(),(1-self.bottom)*self.height())
                painter.drawLine(self.mouseUpPos.x(),(1-self.top)*self.height(),self.mouseUpPos.x(),(1-self.bottom)*self.height())


        if self.movingLegend == True:
            if self.ax[0].get_legend():
                bbox = self.ax[0].get_legend().get_window_extent()
                bbox2 = bbox.transformed(self.ax[0].get_legend().axes.transAxes.inverted())
                painter.drawLine(pos.x(),pos.y(),pos.x()+bbox.width,pos.y())
                painter.drawLine(pos.x(),pos.y()-bbox.height,pos.x()+bbox.width,pos.y()-bbox.height)
                painter.drawLine(pos.x(),pos.y(),pos.x(),pos.y()-bbox.height)
                painter.drawLine(pos.x()+bbox.width,pos.y(),pos.x()+bbox.width,pos.y()-bbox.height)
            else:
                self.movingLegend = False
        else:
            if pos.x() > self.left*self.width() and pos.x() < self.right*self.width() and pos.y() > (1-self.top)*self.height() and pos.y() < (1-self.bottom)*self.height():
                painter.drawLine(pos.x(),int((1-self.top)*self.height()),pos.x(),int((1-self.bottom)*self.height()))
                painter.drawLine(int(self.left*self.width()),pos.y(),int(self.right*self.width()),pos.y())
        if self._painted == False:
            self._painted = True
            self.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(self.width(),self.height()),QtCore.QSize(self.width(),self.height())))
            self.legend()
        painter.end()

    def autoscale_based_on(self,ax):
        lines = ax.get_lines()
        ax.dataLim = mtransforms.Bbox([(1,1),(0,0)])
        for line in lines:
            if not hasattr(line,"ignore"):
                xy = numpy.vstack(line.get_data()).T
                ax.dataLim.update_from_data_xy(xy, ignore=False)
        ax.autoscale_view()

    def plot(self,*args,**kwargs):
        axis = self.ax
        if "rightaxis" in kwargs:
            rightaxis = kwargs.pop("rightaxis")
            if rightaxis:
                axis = self.rax

        self.remove_break_lines_x()

        axis[0].set_autoscale_on(True)

        if kwargs or args:
            thePlot = axis[0].plot(*args,**kwargs)
        else:
            thePlot = None

        if self.xbreak or self.ybreak:
            axis[0].set_autoscale_on(False)
            self.show_break_lines()

        self.autoscale_based_on(axis[0])
        axis[0].set_autoscale_on(False)

        if self.custom_xlim: axis[0].set_xlim(self.custom_xlim)
        if self.custom_rylim and axis is self.rax:
            axis[0].set_ylim(self.custom_rylim)
        if self.custom_ylim and axis is not self.rax:
            axis[0].set_ylim(self.custom_ylim)

        for l in axis[0].get_lines():
            if not hasattr(l,"ignore"):
                if not hasattr(l,"ccp4_uuid"):
                    if sys.platform != 'win32':
                        import uuid
                        if sys.version_info >= (3,0):
                            uuid_str = uuid.uuid4().hex
                        else:
                            uuid_str = uuid.uuid4().get_hex()
                    else:
                        import msilib
                        uuid_str = msilib.gen_uuid().strip('{').strip('}').replace('-','')
                    l.ccp4_uuid = uuid_str
                for theAxis in axis[1:]:
                    for l2t in theAxis.get_lines():
                        if hasattr(l2t,"ccp4_uuid") and l2t.ccp4_uuid == l.ccp4_uuid:
                            haveThis = True
                            break
                    else:
                        l2 = matplotlib.lines.Line2D(l.get_xdata(),l.get_ydata(),linewidth=l.get_linewidth(),linestyle=l.get_linestyle(),
                        color=l.get_color(),markersize=l.get_markersize(),marker=l.get_marker(),label=l.get_label(),picker=l.get_picker()
                        )
                        if not l.get_visible():
                            l2.set_visible(0)
                        theAxis.add_line(l2)
                        l2.ccp4_uuid = l.ccp4_uuid

        theXLim = axis[0].get_xlim()
        theYLim = axis[0].get_ylim()
        theXRange = theXLim[1]-theXLim[0]
        theYRange = theYLim[1]-theYLim[0]

        for ax in axis:
            ax.set_xlim(theXLim[0],theXLim[1])
            ax.set_ylim(theYLim[0],theYLim[1])
        if self.xbreak:
            axis[0].set_xlim(theXLim[0],self.xbreak[0][0])
            nbreaks = len(self.xbreak)
            for i in range(nbreaks-1):
                axis[i+1].set_xlim(self.xbreak[i][1],self.xbreak[i+1][0])
                axis[i+1].set_ylim(theYLim[0],theYLim[1])
                fmt = ScalarFormatter()
                fmt.set_useOffset(0)
                axis[i+1].xaxis.set_major_formatter(fmt)
            axis[-1].set_xlim(self.xbreak[-1][1],theXLim[1])
            axis[-1].set_ylim(theYLim[0],theYLim[1])

        elif self.ybreak:
            axis[-1].set_ylim(theYLim[0],self.ybreak[0][0])
            nbreaks = len(self.ybreak)
            for i in range(nbreaks-1):
                axis[nbreaks-i-1].set_ylim(self.ybreak[i][1],self.ybreak[i+1][0])
                axis[nbreaks-i-1].set_xlim(theXLim[0],theXLim[1])
            axis[0].set_ylim(self.ybreak[-1][1],theXLim[1])
            axis[0].set_xlim(theXLim[0],theXLim[1])

        xticks = 5
        if self.xbreak:
            xticks -= len(self.xbreak)-1
        yticks = 5
        if self.ybreak:
            yticks -= len(self.ybreak)-1
        
        if self.xintegral: xticks += 1
        if self.yintegral: yticks += 1

        xticks = max(xticks,1)
        yticks = max(yticks,1)

        for ax in axis:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(xticks,integer=self.xintegral))
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(yticks,integer=self.yintegral))

        for ax in axis:
            if self.x_scaling: ax.xaxis.set_major_formatter(FuncFormatter(self.x_scaling))
            if self.y_scaling: ax.yaxis.set_major_formatter(FuncFormatter(self.y_scaling))

        if self.xbreak and self.custom_xlim:
            xmin,xmax = axis[0].get_xlim()
            axis[0].set_xlim(self.custom_xlim[0],xmax)
            xmin,xmax = axis[-1].get_xlim()
            axis[-1].set_xlim(xmin,self.custom_xlim[1])
        elif self.ybreak and self.custom_ylim:
            ymin,ymax = axis[0].get_ylim()
            axis[0].set_ylim(ymin,self.custom_ylim[1])
            ymin,ymax = axis[-1].get_ylim()
            axis[-1].set_ylim(self.custom_ylim[0],ymax)

        """
        See q.log for why this might be useful.
        if not self.custom_xlim:
            self.ax[-1].set_xlim(self.ax[-1].get_xlim()[0],self.ax[0].get_lines()[0].get_xdata()[-1])
        """

        return thePlot

    def remove_break_lines_x(self):
        for diag in self.diags:
            for l in diag:
                l.remove()
        self.diags = []
            
    def show_break_lines(self):
        d = 0.015
        self.remove_break_lines_x()
        if self.xbreak:
            kwargs = dict(transform=self.ax[0].transAxes, color='k', clip_on=False)
            self.diags.append(self.ax[0].plot((1-d,1+d),(-d,d),**kwargs))
            self.diags[-1][0].ignore = True
            self.diags.append(self.ax[0].plot((1-d,1+d),(1-d,1+d),**kwargs))
            self.diags[-1][0].ignore = True
            kwargs = dict(transform=self.ax[-1].transAxes, color='k', clip_on=False)
            self.diags.append(self.ax[-1].plot((-d,+d),(-d,d),**kwargs))
            self.diags[-1][0].ignore = True
            self.diags.append(self.ax[-1].plot((-d,+d),(1-d,1+d),**kwargs))
            self.diags[-1][0].ignore = True
            for ax in self.ax[1:-1]:
                kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
                self.diags.append(ax.plot((1-d,1+d),(-d,d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((1-d,1+d),(1-d,1+d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((-d,+d),(-d,d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((-d,+d),(1-d,1+d),**kwargs))
                self.diags[-1][0].ignore = True
        elif self.ybreak:
            kwargs = dict(transform=self.ax[0].transAxes, color='k', clip_on=False)
            self.diags.append(self.ax[0].plot((1-d,1+d),(-d,d),**kwargs))
            self.diags[-1][0].ignore = True
            self.diags.append(self.ax[0].plot((-d,d),(-d,d),**kwargs))
            self.diags[-1][0].ignore = True
            kwargs = dict(transform=self.ax[-1].transAxes, color='k', clip_on=False)
            self.diags.append(self.ax[-1].plot((1-d,1+d),(1-d,1+d),**kwargs))
            self.diags[-1][0].ignore = True
            self.diags.append(self.ax[-1].plot((-d,d),(1-d,1+d),**kwargs))
            self.diags[-1][0].ignore = True
            for ax in self.ax[1:-1]:
                kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
                self.diags.append(ax.plot((1-d,1+d),(1-d,1+d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((1-d,1+d),(-d,d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((-d,d),(1-d,1+d),**kwargs))
                self.diags[-1][0].ignore = True
                self.diags.append(ax.plot((-d,d),(-d,d),**kwargs))
                self.diags[-1][0].ignore = True
            
    def hist(self,*args,**kwargs):
        return self.ax[0].hist(*args,**kwargs)

    def set_breaks(self,xbreak,ybreak):

        self.xbreak = xbreak
        self.ybreak = ybreak

        self.breaksize = 0.02

        for ax in self.ax[1:]:
            self.fig.delaxes(ax)
        self.ax = [self.ax[0]]
        self.gs = [self.gs[0]]

        self.gs[0].update(left=self.left, right=self.right, top=self.top, bottom=self.bottom, wspace=0.0)
        self.ax[0].spines['right'].set_visible(True)
        self.ax[0].spines['bottom'].set_visible(True)
        self.ax[0].tick_params(labelright='off')
        self.ax[0].tick_params(labelbottom='on')

        rcoltick = None
        if len(self.rax)>0:
                if len(self.rax[-1].get_yticklabels())>0:
                    rcoltick = self.rax[-1].get_yticklabels()[0].get_color()
                ryrange = self.rax[-1].get_ylim()
                self.rax[-1].set_ylabel('')
                for i in range(len(self.rax)):
                    self.rax[i].spines['left'].set_visible(False)
                    self.rax[i].spines['right'].set_visible(False)
                    self.rax[i].tick_params(labelleft='off',left='off')
                    self.rax[i].tick_params(labelright='off',right='off')

        for i in range(len(self.ax)):
            self.ax[i].spines['left'].set_visible(False)
            self.ax[i].spines['right'].set_visible(False)
            self.ax[i].tick_params(labelleft='off',left='off')
            self.ax[i].tick_params(labelright='off',right='off')

        if self.xbreak:
            nbreaks = len(self.xbreak)+1
            self.gs[0].update(left=self.left, right=(self.left+self.right)/float(nbreaks)-0.5*self.breaksize, top=self.top, bottom=self.bottom, wspace=0.0)
            for i in range(1,nbreaks):
                self.gs.append(gridspec.GridSpec(1,1))
                self.gs[i].update(left=i*(self.left+self.right)/float(nbreaks)+0.5*self.breaksize, right=min(self.right,(i+1)*(self.left+self.right)/float(nbreaks)-0.5*self.breaksize), top=self.top, bottom=self.bottom, wspace=0.0)
                self.ax.append(self.fig.add_subplot(self.gs[i][:]))
                self.ax[i].spines['left'].set_visible(False)
                self.ax[i].spines['right'].set_visible(False)
                self.ax[i].tick_params(labelleft='off',left='off')
                self.ax[i].tick_params(labelright='off',right='off')
                if len(self.rax)>0:
                    #self.rax.append(self.fig.add_subplot(self.gs[i][:]))
                    self.rax.append(self.ax[-1].twinx())
                    self.rax[i].spines['left'].set_visible(False)
                    self.rax[i].spines['right'].set_visible(False)
                    self.rax[i].tick_params(labelleft='off',left='off')
                    self.rax[i].tick_params(labelright='off',right='off')

            self.ax[0].set_title('')
            self.ax[0].set_xlabel('')
            theText = self.ax[int(nbreaks/2.)].set_title(self.title)
            theText.set_fontproperties(self.title_font)
            if nbreaks%2==0:
                theText.set_position((0,1))
            theText = self.ax[int(nbreaks/2.)].set_xlabel(self.xlabel)
            theText.set_fontproperties(self.label_font)
            if nbreaks%2==0:
                theText.set_position((0,1))

            self.ax[0].yaxis.tick_left()
            self.ax[-1].yaxis.tick_right()
            self.ax[0].spines['left'].set_visible(True)
            self.ax[0].spines['right'].set_visible(False)
            self.ax[-1].spines['right'].set_visible(True)
            self.ax[-1].tick_params(labelright='off')
            # It must be a matplotlib bug that this is necessary.
            for i in range(1,len(self.ax)-1):
                self.ax[i].tick_params(labelright='off',right='off')
                self.ax[i].tick_params(labelleft='off',left='off')
            if len(self.rax)>0:
                self.rax[-1].yaxis.tick_right()
                self.rax[0].spines['left'].set_visible(False)
                self.rax[0].spines['right'].set_visible(False)
                self.rax[0].spines['bottom'].set_visible(False)
                self.rax[-1].spines['right'].set_visible(True)
                self.rax[-1].tick_params(labelright='on')
                self.rax[0].tick_params(labelright='off')
                self.rax[-1].tick_params(labelleft='off')
                self.rax[-1].tick_params(right='on')
                self.rax[0].tick_params(right='off')
                self.rax[-1].tick_params(left='off')
                self.rax[0].tick_params(labelbottom='off')
                self.rax[-1].tick_params(labelbottom='off')
                if rcoltick is not None:
                    self.rax[-1].spines['right'].set_color(rcoltick)
                    self.ax[-1].spines['right'].set_color(rcoltick)
                    for tl in self.rax[-1].get_yticklabels():
                            tl.set_color(rcoltick)
                self.rax[0].set_ylabel('')
                if hasattr(self,"rylabel") and len(self.rylabel)>0:
                    if rcoltick is not None:
                        theText = self.rax[-1].set_ylabel(self.rylabel,color=rcoltick)
                    else:
                        theText = self.rax[-1].set_ylabel(self.rylabel)
                    #theText.set_fontproperties(self.label_font)
                    #print self.label_font
                self.rax[-1].set_ylim(ryrange)

        elif self.ybreak:
            nbreaks = len(self.ybreak)+1
            title = self.ax[0].get_title()
            theText = self.ax[0].set_title(title)
            self.gs[0].update(left=self.left, right=self.right, top=self.top, bottom=(nbreaks-1)*(self.top+self.bottom)/float(nbreaks)+0.5*self.breaksize, wspace=0.0)
            for i in range(1,nbreaks):
                self.gs.append(gridspec.GridSpec(1,1))
                self.ax.append(self.fig.add_subplot(self.gs[i][:]))
            for i in range(1,nbreaks):
                self.gs[nbreaks-i].update(left=self.left, right=self.right, top=i*(self.top+self.bottom)/float(nbreaks)-0.5*self.breaksize, bottom=max(self.bottom,(i-1)*(self.top+self.bottom)/float(nbreaks)+0.5*self.breaksize), wspace=0.0)
                self.ax[nbreaks-i].spines['top'].set_visible(False)
                self.ax[nbreaks-i].spines['bottom'].set_visible(False)
                self.ax[nbreaks-i].tick_params(labeltop='off',top='off')
                self.ax[nbreaks-i].tick_params(labelbottom='off',bottom='off')

            self.ax[0].set_ylabel('')
            theText = self.ax[int(nbreaks/2.)].set_ylabel(self.ylabel)
            theText.set_fontproperties(self.label_font)
            if nbreaks%2==0:
                theText.set_position((0,1))

            theText = self.ax[-1].set_xlabel(self.xlabel)
            self.ax[0].set_xlabel('')
            theText.set_fontproperties(self.label_font)

            self.ax[0].xaxis.tick_top()
            self.ax[-1].xaxis.tick_bottom()
            self.ax[0].spines['top'].set_visible(True)
            self.ax[0].spines['bottom'].set_visible(False)
            self.ax[-1].spines['bottom'].set_visible(True)
            self.ax[0].tick_params(labeltop='off')
        else:
            theText = self.ax[0].set_xlabel(self.xlabel)
            theText.set_fontproperties(self.label_font)
            self.ax[0].yaxis.tick_left()
            self.ax[-1].yaxis.tick_right()
            self.ax[0].spines['left'].set_visible(True)
            self.ax[-1].spines['right'].set_visible(True)
            self.ax[0].tick_params(labelleft='on',left='on')
            self.ax[-1].tick_params(labelright='off',right='on')
            if len(self.rax)>0:
                self.rax[-1].yaxis.tick_right()
                self.rax[0].spines['left'].set_visible(False)
                self.rax[0].spines['right'].set_visible(False)
                self.rax[0].spines['bottom'].set_visible(False)
                self.rax[-1].spines['right'].set_visible(True)
                self.rax[-1].tick_params(labelright='on')
                self.rax[0].tick_params(labelright='off')
                self.rax[-1].tick_params(labelleft='off')
                self.rax[-1].tick_params(right='on')
                self.rax[0].tick_params(right='off')
                self.rax[-1].tick_params(left='off')
                self.rax[0].tick_params(labelbottom='off')
                self.rax[-1].tick_params(labelbottom='off')
                if rcoltick is not None:
                    self.rax[-1].spines['right'].set_color(rcoltick)
                    self.ax[-1].spines['right'].set_color(rcoltick)
                    for tl in self.rax[-1].get_yticklabels():
                            tl.set_color(rcoltick)
                self.rax[0].set_ylabel('')
                if hasattr(self,"rylabel") and len(self.rylabel)>0:
                    if rcoltick is not None:
                        theText = self.rax[-1].set_ylabel(self.rylabel,color=rcoltick)
                    else:
                        theText = self.rax[-1].set_ylabel(self.rylabel)
                    #theText.set_fontproperties(self.label_font)
                    #print self.label_font
                self.rax[-1].set_ylim(ryrange)

        # Try to make sure first plot is the top one.
        self.ax[0].set_zorder(1000)
        if len(self.rax)>0:
            self.rax[0].set_zorder(1001)
            """
            for ax in self.rax:
                ax.set_zorder(1001)
            """

        self.fig.subplots_adjust(self.left, self.bottom, self.right, self.top)
        self.set_tick_labels_font()

    def set_tick_labels_font(self):
        labels_x = []
        labels_y = []
        for ax in self.ax:
            labels_x.extend(ax.get_xticklabels())
            labels_y.extend(ax.get_yticklabels())
        for ax in self.rax:
            labels_x.extend(ax.get_xticklabels())
            labels_y.extend(ax.get_yticklabels())

        for xlabel in labels_x:
            xlabel.set_fontproperties(self.axes_font)
        for ylabel in labels_y:
            ylabel.set_fontproperties(self.axes_font)


    def __init__(self,fig=Figure(figsize = (800, 600), dpi=300, facecolor = '#FFFFFF'),title='',x_axis_label='',y_axis_label='',parent=None):
        self.fig = fig
        FigureCanvas.__init__(self,fig)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)

        self._type = 'xy'
        self.diags = []
        self.xbreak = None
        self.ybreak = None
        self.x_scaling = None
        self.y_scaling = None
        self.custom_xlim = None
        self.custom_ylim = None
        self.custom_rylim = None
        self.fix_aspect_ratio = False

        self.initial_legend_position = [0.75,0.75]

        self.gs = [gridspec.GridSpec(1,1)]
        self.ax = [self.fig.add_subplot(self.gs[0][:])]
        self.rax = []
        self.defaultTextHeight = QtGui.QFontMetrics(QtGui.QFont("Bitstream Vera Sans",10)).height()

        self.top = 0.95
        self.left = 0.1
        self.bottom=0.1
        self.right=0.97

        self.xintegral = False
        self.yintegral = False

        #font = FontProperties(family="sans-serif",style="normal",weight="normal",size="10")
        font = FontProperties(family="Bitstream Vera Sans",size=10)

        self.title_font = font.copy()
        #self.title_font.set_weight('bold')
        self.label_font = font.copy()
        #self.label_font.set_weight('bold')
        self.axes_font = font.copy()
        #self.axes_font.set_weight('bold')
        self.legend_font = font.copy()
        #self.legend_font.set_weight('bold')

        self.title = title
        self.xlabel = x_axis_label
        self.ylabel = y_axis_label

        theText = self.ax[0].set_title(self.title,fontproperties=self.title_font)
        self.ax[0].set_ylabel(y_axis_label,fontproperties=self.label_font)
        self.ax[0].set_xlabel(x_axis_label,fontproperties=self.label_font)
        #The matplotlib installed via pip with python3 on os x does not have hold.
        if hasattr(self.ax[0],"hold"):
            self.ax[0].hold(True)


        self.setMouseTracking(True)
        self.mpl_connect('button_press_event', self.onClick)
        self.mpl_connect('motion_notify_event', self.onMove)
        self.mpl_connect('pick_event',self.onpick)
        self.movingLegend = False

        self.fig.subplots_adjust(self.left, self.bottom, self.right, self.top)
        self.set_tick_labels_font()
        self.defaultSymbolSize = 5
        self._painted = False

    def onpick(self,e):
        try:
            if hasattr(e,"mouseevent"):
                mouseevent = e.mouseevent
            elif hasattr(e,"lastevent"):
                mouseevent = e.lastevent
            dumX = float(mouseevent.xdata)
            dumY = float(mouseevent.ydata)
            if len(self.rax)>0:
                yfrac = (self.rax[0].get_ylim()[0]-mouseevent.ydata)/(self.rax[0].get_ylim()[0]-self.rax[0].get_ylim()[1])
                ydatal = self.ax[0].get_ylim()[0] - yfrac*(self.ax[0].get_ylim()[0]-self.ax[0].get_ylim()[1])
                self.MouseMoveEvent.emit((mouseevent.xdata,ydatal))
            else:
                self.MouseMoveEvent.emit((mouseevent.xdata,mouseevent.ydata))
            if len(self.rax)>0:
                hl = self.fig.hitlist(mouseevent)
                xlim = self.ax[0].get_xlim()
                ylim = self.ax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol
                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-mouseevent.xdata)).argmin()
                        #print mouseevent.xdata-xdata[x_idx], xtol
                        if abs(mouseevent.xdata-xdata[x_idx])<xtol and abs(ydatal-ydata[x_idx])<ytol:
                            self.DetectPointClickEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
                xlim = self.rax[0].get_xlim()
                ylim = self.rax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol
                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-mouseevent.xdata)).argmin()
                        #print mouseevent.xdata-xdata[x_idx], xtol
                        if abs(mouseevent.xdata-xdata[x_idx])<xtol and abs(mouseevent.ydata-ydata[x_idx])<ytol:
                            self.DetectPointClickEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
            else:
                hl = self.fig.hitlist(mouseevent)
                xlim = self.ax[0].get_xlim()
                ylim = self.ax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol
                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-mouseevent.xdata)).argmin()
                        #print mouseevent.xdata-xdata[x_idx], xtol
                        if abs(mouseevent.xdata-xdata[x_idx])<xtol and abs(mouseevent.ydata-ydata[x_idx])<ytol:
                            self.DetectPointClickEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
        except:
            pass
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            sys.stderr.write(str(exc_type)+'\n')
            sys.stderr.write(str(exc_value)+'\n')

    def onMove(self,e):
        try:
            dumX = float(e.xdata)
            dumY = float(e.ydata)
            if len(self.rax)>0:
                yfrac = (self.rax[0].get_ylim()[0]-e.ydata)/(self.rax[0].get_ylim()[0]-self.rax[0].get_ylim()[1])
                ydatal = self.ax[0].get_ylim()[0] - yfrac*(self.ax[0].get_ylim()[0]-self.ax[0].get_ylim()[1])
                self.MouseMoveEvent.emit((e.xdata,ydatal))
            else:
                self.MouseMoveEvent.emit((e.xdata,e.ydata))

            if len(self.rax)>0:
                hl = self.fig.hitlist(e)
                xlim = self.ax[0].get_xlim()
                ylim = self.ax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol

                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-e.xdata)).argmin()
                        if abs(e.xdata-xdata[x_idx])<xtol and abs(ydatal-ydata[x_idx])<ytol:
                            self.DetectPointMoveEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
                xlim = self.rax[0].get_xlim()
                ylim = self.rax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol

                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-e.xdata)).argmin()
                        if abs(e.xdata-xdata[x_idx])<xtol and abs(e.ydata-ydata[x_idx])<ytol:
                            self.DetectPointMoveEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
            else:
                hl = self.fig.hitlist(e)
                xlim = self.ax[0].get_xlim()
                ylim = self.ax[0].get_ylim()
                xtol = self.defaultSymbolSize*abs(xlim[1]-xlim[0])/self.width()
                ytol = self.defaultSymbolSize*abs(ylim[1]-ylim[0])/self.height()
                #print xtol, ytol

                for hit in hl:
                    if type(hit)==matplotlib.lines.Line2D:
                        xdata = hit.get_xdata()
                        ydata = hit.get_ydata()
                        x_idx = (numpy.abs(xdata-e.xdata)).argmin()
                        if abs(e.xdata-xdata[x_idx])<xtol and abs(e.ydata-ydata[x_idx])<ytol:
                            self.DetectPointMoveEvent.emit({'plotwidget':self.parent(),'label':hit.get_label(), 'x':xdata[x_idx],'y':ydata[x_idx]})
        except:
            pass

    def onClick(self,event):
        if self.movingLegend == True:
            self.movingLegend = False
            if 1:#try:
                loc = self.ax[0].get_legend().parent.transAxes.inverted().transform_point((event.x, event.y))
                self.legend(loc = 'lower left', bbox_to_anchor = loc)
            else:#except:
                print("Failed to move legend")
            self.draw()
            self.repaint()
        try:
            dumX = float(event.xdata)
            dumY = float(event.ydata)
            #print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(event.button, event.x, event.y, event.xdata, event.ydata)
            self.onpick(event)
        except:
            pass

    def addPoint(self,point):
        x,y = point
        self.plot(x, y, 'o')
        self.draw()

class MyCanvas(QtMatplotlibCanvas):
    def __init__(self,fig=Figure(figsize = (800, 600), dpi=300, facecolor = '#FFFFFF'),parent=None):
        QtMatplotlibCanvas.__init__(self,fig,parent)

    def mousePressEvent(self,e):
        xlim = self.ax[0].get_xlim()
        figWidth = self.width()*(float(self.right)-self.left)
        xp = (float(e.pos().x())-float(self.left)*self.width())/figWidth
        x = xp*(xlim[1]-xlim[0])+xlim[0]

        ylim = self.ax[0].get_ylim()
        figHeight = self.height()*(float(self.top)-self.bottom)
        yp = (self.height()-float(e.pos().y())-float(self.bottom)*self.height())/figHeight
        y = yp*(ylim[1]-ylim[0])+ylim[0]

        self.plot(x, y, 'o')
        self.draw()
        self.fig.savefig("figure.png")
        self.fig.savefig("figure.ps")
        self.fig.savefig("figure.pdf")

def CCP4LogToEtree(b):
    splits = b[b.find("$TABLE"):].split("$TABLE")
    newsplits = []
    gs = []

    bigtree = etree.Element("CCP4ApplicationOutput")
    for ns in splits[1:]:
        ns = "$TABLE"+ns
        newsplits.append(ns)
        table = CCP4Table(ns)
        tree = table.toEtree()
        bigtree.append(tree)
    return bigtree

def CCP4LogToXML(b):
    splits = b[b.find("$TABLE"):].split("$TABLE")
    newsplits = []
    gs = []

    status_xml = ""
    header ="""<?xml version="1.0" encoding="UTF-8" ?>\n"""
    status_xml += header
    
    NSMAP = {'xsi':"http://www.w3.org/2001/XMLSchema-instance"}
    NS = NSMAP['xsi']
    location_attribute = '{%s}noNamespaceSchemaLocation' % NS
    bigtree = etree.Element("CCP4ApplicationOutput",nsmap = NSMAP,attrib={location_attribute: 'http://www.ysbl.york.ac.uk/~mcnicholas/schema/CCP4ApplicationOutput.xsd'})
    for ns in splits[1:]:
        ns = "$TABLE"+ns
        newsplits.append(ns)
        table = CCP4Table(ns)
        tree = table.toEtree()
        bigtree.append(tree)

    status_xml += etree.tostring(bigtree,encoding='utf-8', pretty_print=True)
    return status_xml

def CCP4LogFileNameToXML(f):
    if f.endswith('.log') or f.endswith('.txt'):
        fo = open(f)
        b = fo.read()
        fo.close()
        return CCP4LogToXML(b)
    else:
        return CCP4LogToXML("")

def CCP4LogFileNameToEtree(f):
    if f.endswith('.log') or f.endswith('.txt'):
        fo = open(f)
        b = fo.read()
        fo.close()
        return CCP4LogToEtree(b)
    else:
        return CCP4LogToEtree("")

if __name__ == "__main__":
    if "-quit" in sys.argv:
        sys.exit()

    app = QtWidgets.QApplication(sys.argv)
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
        app.addLibraryPath(os.path.join(basedir,"qt4_plugins"))
    else:
        #app.addLibraryPath(os.path.join(os.environ['CCP4MG'],"QtPlugins"))
        app.addLibraryPath(os.path.join(os.path.dirname(__file__),'..',"QtPlugins"))

    """
    t = 0
    timer = QtCore.QTimer()
    win = QtMatplotlibWidget()

    def addRandomPoint():
        import random
        global t, timer, win
        t = t + 1
        v = random.randint(1,20)
        print t*timer.interval()/1000., v
        win.addPoint((t,v))
    """

    @QtCore.Slot(dict)
    def regionTest(vals):
        print("Region VALS",vals)

    @QtCore.Slot(dict)
    def clickTest(vals):
        print("Click VALS",vals)

    xsd = None
    hist = False
    nomainwin = False
    testla = False
    select = []

    import getopt

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hnx:s:g:G:", ["hist", "nomainwin","testla", "xsd=", "select=","graph=","graphTotal="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        sys.exit(2)

    theGraph = None
    theGraphTotal = None

    for o, a in opts:
        if o in ("-h", "--hist"):
            hist = True
        elif o in ("-n", "--testla"):
            testla = True
        elif o in ("-n", "--nomainwin"):
            nomainwin = True
        elif o in ("-x", "--xsd"):
            xsd = a
        elif o in ("-s", "--select"):
            items = a.strip('[').strip(']').split(',')
            for item in items: select.append(item.strip("'").strip('"'))
        elif o in ("-g", "--graph"):
            items = a.strip('[').strip(']').split(',')
            theGraph = items
        elif o in ("-G", "--graphTotal"):
            theGraphTotal = int(a)
        else:
            assert False, "unhandled option: "+o

    if len(select)==0: select = None
    if theGraph is not None and len(theGraph)!=2: theGraph = None

    if nomainwin:

        win = LogGraph()
        win.layout().setContentsMargins(0,0,0,0)
        win.showActionBar()
        win.show()
        win.raise_()
        if len(sys.argv)>1:
            for f in args:
                graphs = win.loadFile(f,hist=hist,xsd=xsd)
                if len(graphs)>0:
                    graph = graphs[0]
                    if graph:
                        graph.DetectPointClickEvent.connect(clickTest)
                        graph.RegionSelected.connect(regionTest)
    elif testla:
        owin = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()
        owin.setLayout(layout)
        win = LogGraph()
        layout.addWidget(win)
        owin.show()
        owin.raise_()
        if len(sys.argv)>1:
            for f in args:
                graphs = win.loadFile(f,hist=hist,xsd=xsd)
                if len(graphs)>0:
                    graph = graphs[0]
                    if graph:
                        graph.DetectPointClickEvent.connect(clickTest)
                        graph.RegionSelected.connect(regionTest)
        selectorWidget = QtWidgets.QWidget()
        selectorLayout = QtWidgets.QVBoxLayout()
        selectorWidget.setLayout(selectorLayout)
        layout.addWidget(selectorWidget)
        @QtCore.Slot()
        def drawLineSelector():
            import sip
            child = selectorWidget.layout().takeAt(0)
            while child:
                if child:
                    if hasattr(child,"widget") and child.widget() is not None:
                        sip.delete(child.widget())
                    sip.delete(child)
                    child = selectorWidget.layout().takeAt(0)

            if win.graph and win.graph.currentWidget() and win.graph.currentWidget().currentWidget():
                currentGraph = win.graph.currentWidget().currentWidget()
                if hasattr(currentGraph,"canvas") and currentGraph.canvas:
                    if hasattr(currentGraph.canvas,"fig") and currentGraph.canvas.fig:
                        if len(currentGraph.canvas.fig.get_axes())>0:
                            lines = currentGraph.canvas.fig.get_axes()[0].get_lines()
                            lines = [l for l in lines if not hasattr(l,"ignore")]
                            for l in lines:
                                cb = QtWidgets.QCheckBox(l.get_label())
                                if l.get_visible():
                                    cb.setChecked(True)
                                selectorWidget.layout().addWidget(cb)
                                cb.stateChanged.connect(l.set_visible)
                                cb.stateChanged.connect(currentGraph.canvas.draw)
        win.currentTableChanged.connect(drawLineSelector)
        win.currentDataChanged.connect(drawLineSelector)
        drawLineSelector()
        
    else:
        win = LogGraphMainWindow()
        win.layout().setContentsMargins(0,0,0,0)
        win.show()
        win.raise_()
        if len(sys.argv)>1:
            for f in args:
                graphs = win.loggraph.loadFile(f,hist=hist,xsd=xsd,select=select)
                if graphs is not None and len(graphs)>0:
                    graph = graphs[0]
                    if graph:
                        graph.DetectPointClickEvent.connect(clickTest)
                        graph.RegionSelected.connect(regionTest)
                if theGraph is not None:
                    try:
                        win.loggraph.table_combo.setCurrentIndex(int(theGraph[0]))
                        win.loggraph.data_combo.setCurrentIndex(int(theGraph[1]))
                    except:
                        print("Failed to set graph",theGraph)
                if theGraphTotal is not None:
                                        dataTot = 0
                                        for i in range(win.loggraph.table_combo.count()):
                                                win.loggraph.table_combo.setCurrentIndex(i)
                                                dataTot += win.loggraph.data_combo.count()
                                                if theGraphTotal < dataTot:
                                                        diff = dataTot - theGraphTotal
                                                        win.loggraph.data_combo.setCurrentIndex(win.loggraph.data_combo.count()-diff)
                                                        break

    sys.exit(app.exec_())
