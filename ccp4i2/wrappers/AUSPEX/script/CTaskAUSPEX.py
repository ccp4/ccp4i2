from ....qtgui import CCP4TaskWidget


class CTaskAUSPEX(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'AUSPEX'
    TASKVERSION = 0.0
    TASKMODULE ='data_reduction'
    TASKTITLE = 'Graphical diagnostics by AUSPEX plots'
    SHORTTASKTITLE = "AUSPEX"
    DESCRIPTION = 'Use AUSPEX, generate graphical diagnostics for data set'
    WHATNEXT = []

    def __init__(self,parent):
        CCP4TaskWidget.CTaskWidget.__init__(self,parent)

    def drawContents(self):
        folder = self.openFolder(folderFunction='inputData', title='Input Data')
        self.createLine(['subtitle', 'Select input data', 'Observed intensties (or amplitudes) are required'])
        self.openSubFrame(frame=[True])
        self.createLine(['widget', 'F_SIGF'])
        
        self.closeSubFrame()
        self.createLine( ['subtitle' , 'Plot generation parameters' ] )
       
        self.openSubFrame(frame=[True])
        
        self.createLine(['label', 'Range along y axis:', 'widget', 'YLIM', 'advice', '', 'tip', 'Specify the maximum resolution to show in the plots.'])
        
        self.createLine(['label', 'High resolution cut-off:', 'widget', 'DLIM', 'advice', 'No cut-off includes all data.','tip', 'Specify the maximum resolution to show in the plots.'])
        
        self.createLine(['label', 'Put all plots in one figure:', 'widget', 'SINGFIG', 'advice', '', 'tip', 'Should the images be generated in separate png files? (Default: No.)'])
        
        self.createLine(['label', 'Flag suspected ice rings red:', 'widget', 'FLAGICE', 'advice', 'Automatic ice ring detection.','tip', 'If set, suspected ice rings will be flagged by red bars.'])
        
        self.closeSubFrame()
