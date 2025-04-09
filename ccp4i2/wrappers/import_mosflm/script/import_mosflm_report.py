import xml.etree.ElementTree as etree

from ....report.CCP4ReportParser import Graph, Report


class import_mosflm_report(Report):

    TASKNAME = 'import_mosflm'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        if jobStatus.lower() == "nooutput":
            return
        if len(self.xmlnode.findall('.//integration_result'))==0:
            self.append('No integration data found from iMosflm')
            return
        self.append("<p><h2><b>Mosflm Integration Results</b></h2></p>")
        clearingDiv = self.addDiv(style="clear:both;")
        detectorFold = self.addFold(label="Detector and beam", initiallyOpen=True)
        self.beamGraph(parent=detectorFold)
        self.centralSpotProfile(parent=detectorFold)
        clearingDiv = self.addDiv(style="clear:both;")
        crystalFold = self.addFold(label="Crystal", initiallyOpen=True)
        self.crystalGraph(parent=crystalFold)
        self.blockProfiles(parent=crystalFold)
        clearingDiv = self.addDiv(style="clear:both;")
        dataFold = self.addFold(label="Data", initiallyOpen=True)
        self.dataGraph(parent=dataFold)
        self.histograms(parent=dataFold)
    
    def beamGraph(self, parent=None):
        if parent is None:
            parent = self
        tableData = [("beam_x", "Beam x"), ("beam_y", "Beam y"), ("distance", "Distance"),
                     ("yscale", "Y-scale"), ("tilt", "Tilt"), ("twist", "Twist"),
                     ("tangential_offset", "Tangential offset"), ("radial_offset", "Radial offset"),
                     ("global_absolute_rms_residual", "RMS residual"),
                     ("central_absolute_rms_residual", "RMS res. (central)"),
                     ("global_weighted_rms_residual", "RMS res. (weighted)")]
        dataTuples = []
        imageCount = 0
        integrationNode = self.xmlnode.findall('.//integration_result')[-1]
        for key, label in tableData:
            correspondingData = integrationNode.get(key).split()
            imageCount = len(correspondingData)
            dataTuple = (label,correspondingData)
            dataTuples += (dataTuple,)
        dataTuples.insert(0, ('Image_number', [i for i in range(imageCount)]))
        graph = parent.addFlotGraph(title="Beam Graphs", select=".//integration_result",
                                    style="width:340px;height:290px; float:none; border:0px; margin-left:21px")
        for title, data in dataTuples:
            graph.addData(title=title, data=data)
        plot = graph.addPlotObject()
        plot.append('title', 'Beam positional variation across run')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'pos(mm)')
        l1 = plot.append('plotline', xcol=1, ycol=2)
        l1.append('label','x')
        l1.append('colour','red')
        l2 = plot.append('plotline', xcol=1, ycol=3)
        l2.append('label','y')
        l2.append('colour','orange')
        plot = graph.addPlotObject()
        plot.append('title', 'Detector Tilt and Twist')
        plot.append('plottype', 'xy')
        plot.append('xintegral', 'true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'deg.')
        l2 = plot.append('plotline', xcol=1, ycol=6)
        l2.append('label', 'Tilt')
        l2.append('colour', 'blue')
        l2 = plot.append('plotline', xcol=1, ycol=7)
        l2.append('label', 'Twist')
        l2.append('colour', 'green')
        plot = graph.addPlotObject()
        plot.append('title', 'RMS Beam spot residuals')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'mm')
        l3 = plot.append('plotline', xcol=1, ycol=10)
        l3.append('label', 'res.')
        l3.append('colour','blue')
        l3 = plot.append('plotline', xcol=1, ycol=11)
        l3.append('label', 'res.(central)')
        l3.append('colour','green')

    def centralSpotProfile(self, parent=None):
        if parent is None:
            parent = self
        gallery = parent.addObjectGallery(height='300px', contentWidth='300px', style='float:left;width:420px;')
        profileNodes = self.xmlnode.findall('.//profile')
        firstOne = True
        for iProfile, profileNode in enumerate(profileNodes):
            dataDiv = gallery.addDrawnDiv(data=etree.tostring(profileNode), renderer="mosflm.profileRenderer", 
                                          require='mosflm',initiallyDrawn=firstOne, title="Image "+str(iProfile), label="Image "+str(iProfile))
            firstOne = False

    def crystalGraph(self, parent=None):
        if parent is None:
            parent = self
        tableData = [("phi_x","\u03D5(x)"), ("phi_y","\u03D5(y)"), ("phi_z","\u03D5(z)"), ("cell_a","a"), 
                     ("cell_b","b"),("cell_c","c"),("cell_alpha","\u03B1"),("cell_beta","\u03B2"),
                     ("cell_gamma","\u0263"),("mosaicity","Mosaicity")]
        dataTuples = []
        imageCount = 0
        integrationNode = self.xmlnode.findall('.//integration_result')[-1]
        for key, label in tableData:
            try: 
                correspondingData = integrationNode.get(key).split()
            except:
                print(key)
                break
            imageCount = len(correspondingData)
            dataTuple = (label,correspondingData)
            dataTuples += (dataTuple,)
        dataTuples.insert(0, ('Image_number',[i for i in range(imageCount)]))
        graph = parent.addFlotGraph(title="Crystal Graphs", select=".//integration_result",
                                    style="width:340px;height:290px; float:none; border:0px; margin-left:21px")
        for title, data in dataTuples:
            graph.addData(title=title, data=data)
        plot = graph.addPlotObject()
        plot.append('title', 'Miss-setting angles and Mosaicity vs Image')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'degrees')
        l1 = plot.append('plotline', xcol=1, ycol=2)
        l1.append('label','Miss-set phi(x)')
        l1.append('colour','blue')
        l2 = plot.append('plotline', xcol=1, ycol=3)
        l2.append('label','Miss-set phi(y)')
        l2.append('colour','green')
        l3 = plot.append('plotline', xcol=1, ycol=4)
        l3.append('label','Miss-set phi(z)')
        l3.append('colour','orange')
        l4 = plot.append('plotline', xcol=1, ycol=11)
        l4.append('label','mosaicity')
        l4.append('colour','red')

    def blockProfiles(self, parent=None):
        if parent is None:
            parent = self
        gallery = parent.addObjectGallery(height='300px', contentWidth='300px', style='float:left;width:420px;')
        profileGridNodes = self.xmlnode.findall('.//profile_grid')
        firstOne = True
        for iProfileGridNode, profileGridNode in enumerate(profileGridNodes):
            dataDiv = gallery.addDrawnDiv(data=etree.tostring(profileGridNode), renderer='mosflm.profileGridRenderer', 
                                          require='mosflm', initiallyDrawn=firstOne, title=str(iProfileGridNode),
                                          label="Block "+str(iProfileGridNode))
            firstOne = False

    def dataGraph(self, parent=None):
        if parent is None:
            parent = self
        tableData = [("mean_profile_fitted_fulls","\u003C\u00A0I\u002F\u03c3(I)\u003E(prf,full)"),
                     ("mean_profile_fitted_partials","\u003C\u00A0I\u002F\u03c3(I)\u003E(prf,partials)"),
                     ("mean_summation_integration_fulls","\u003C\u00A0I\u002F\u03c3(I)\u003E(sum,full)"),
                     ("mean_summation_integration_partials","\u003C\u00A0I\u002F\u03c3(I)\u003E(sum,partials)"),
                     ("total_spot_count_fulls","Reflections(full)"),
                     ("total_spot_count_partials","Reflections(partials)"),
                     ("outer_profile_fitted_fulls","\u003C\u00A0I\u002F\u03c3(I)\u003E(prf,full)\u00A0HR"),
                     ("outer_profile_fitted_partials","\u003C\u00A0I\u002F\u03c3(I)\u003E(prf,partials)\u00A0HR"),
                     ("outer_summation_integration_fulls","\u003C\u00A0I\u002F\u03c3(I)\u003E(sum,full)\u00A0HR"),
                     ("outer_summation_integration_partials","\u003C\u00A0I\u002F\u03c3(I)\u003E(sum,partials)\u00A0HR"),
                     ("outer_spot_count_fulls","Reflections(full)\u00A0HR"),
                     ("outer_spot_count_partials","Reflections(partials)\u00A0HR"),
                     ]
        '''
            ("overloads_fulls",u"Overloads(full)\\u00A0HR"),
            ("overloads_partials",u"Overloads(partials)\\u00A0HR"),
            ("badspots_fulls",u"Bad\\u00A0spots(full)\\u00A0HR"),
            ("badspots_partials",u"Bad\\u00A0spots(partials)\\u00A0HR"),
            ("soverlaps_fulls",u"Spatial\\u00A0overlaps(full)\\u00A0HR"),
            ("soverlaps_partials",u"Spatial\\u00A0overlaps(partials)\\u00A0HR"),
        '''
        dataTuples = []
        imageCount = 0
        integrationNode = self.xmlnode.findall('.//integration_result')[-1]
        for key, label in tableData:
            try:
                correspondingData = integrationNode.get(key).split()
            except:
                print(key)
            imageCount = len(correspondingData)
            dataTuple = (label,correspondingData)
            dataTuples += (dataTuple,)
        dataTuples.insert(0, ('Image_number',[i for i in range(imageCount)]))
        graph = parent.addFlotGraph(title="Reflection Info Graphs", select=".//integration_result",
                                    style="width:340px;height:290px; float:none; border:0px; margin-left:21px")
        for title, data in dataTuples:
            graph.addData(title=title, data=data)
        plot = graph.addPlotObject()
        plot.append('title', 'I/sigma(I) Fitted')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=2)
        l1.append('label','full')
        l1.append('colour','blue')
        l2 = plot.append('plotline', xcol=1, ycol=3)
        l2.append('label','partial')
        l2.append('colour','green')
        plot = graph.addPlotObject()
        plot.append('title', 'I/sigma(I) Integrated')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=4)
        l1.append('label','full')
        l1.append('colour','green')
        l2 = plot.append('plotline', xcol=1, ycol=5)
        l2.append('label','partial')
        l2.append('colour','orange')
        plot = graph.addPlotObject()
        plot.append('title', 'Number of Reflections')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=6)
        l1.append('label','full')
        l1.append('colour','red')
        l2 = plot.append('plotline', xcol=1, ycol=7)
        l2.append('label','partial')
        l2.append('colour','orange')
        # Outer Reflections
        plot = graph.addPlotObject()
        plot.append('title', 'I/sigma(I) Fitted (hi-res bin)')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=8)
        l1.append('label','full')
        l1.append('colour','blue')
        l2 = plot.append('plotline', xcol=1, ycol=9)
        l2.append('label','partial')
        l2.append('colour','green')
        plot = graph.addPlotObject()
        plot.append('title', 'I/sigma(I) Integrated (hi-res bin)')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=10)
        l1.append('label','full')
        l1.append('colour','green')
        l2 = plot.append('plotline', xcol=1, ycol=11)
        l2.append('label','partial')
        l2.append('colour','orange')
        plot = graph.addPlotObject()
        plot.append('title', 'Number of Reflections (hi-res bin)')
        plot.append('plottype', 'xy')
        plot.append('xintegral','true')
        plot.append('xlabel', 'image')
        plot.append('ylabel', 'Value')
        l1 = plot.append('plotline', xcol=1, ycol=12)
        l1.append('label','full')
        l1.append('colour','red')
        l2 = plot.append('plotline', xcol=1, ycol=13)
        l2.append('label','partial')
        l2.append('colour','orange')

    def histograms(self, parent=None):
        if parent is None:
            parent = self
        gallery = parent.addObjectGallery(height='300px', contentWidth='300px', style='float:left;width:420px;')
        integrationNode = self.xmlnode.findall('.//integration_result')[-1]
        histogramsText = integrationNode.get('histograms')
        histogramDataRows = histogramsText.split('}')
        collectedData = {}
        for histogramDataRow in histogramDataRows:
            try:
                imageName, histogramData = histogramDataRow.split(',')
                histogramType, dataText = histogramData.split('{')
                imageName = imageName.strip()
                histogramType=histogramType.strip()
                dataText=dataText.strip()
                #print '['+imageName+']['+histogramType+']['+dataText+']'
                if not imageName in collectedData:
                    collectedData[imageName]={}
                collectedData[imageName][histogramType] = dataText
                #print imageName, ' has ', data
            except Exception as e:
                pass
    
        imageNames = sorted([imageName for imageName in collectedData])
        firstOne = True
        for iImageName, imageName in enumerate(imageNames):
            imageHistograms = collectedData[imageName]
            graph = Graph(title=str(iImageName))
            resolutionNumbers = imageHistograms['bin_resolution_limits'].split()[1:-1]
            graph.addData(title='resolution', data=resolutionNumbers)
            coordinate = 2
            colors=['red','green','blue','cyan','magenta','yellow']
            for histogramLabels, plotTitle in [
                                               (['summation_integration_partials','summation_integration_fulls'],'Summation_integration'),
                                               (['profile_fitted_partials','profile_fitted_fulls'],'profile_fitted'),
                                               (['overloads_partials','overloads_fulls'],'overloads'),
                                               (['spot_count_partials','spot_count_fulls'],'spot_count'),
                                               (['soverlaps_partials','soverlaps_fulls'],'soverlaps'),
                                               ]:
                plotObject = graph.addPlotObject()
                plotObject.append('title',plotTitle)
                plotObject.append('plottype','xy')
                plotObject.append('xscale','oneoversqrt')
                plotObject.append('xintegral','false')
                for iHistogramLabel in range(len(histogramLabels)):
                    histogramLabel = histogramLabels[iHistogramLabel]
                    if histogramLabel in imageHistograms:
                        graph.addData(title = histogramLabel, data=imageHistograms[histogramLabel].split()[:-1])
                        plotLine = plotObject.append('plotline',xcol=1,ycol=coordinate,colour=colors[iHistogramLabel%len(colors)])
                        coordinate += 1
            graph.makeTableText()
            drawnDiv = gallery.addDrawnDiv(data=etree.tostring(graph.data_as_etree()), renderer="CCP4i2Widgets.CCP4i2Graph",
                                           require='CCP4i2Widgets', initiallyDrawn=firstOne, width='285px',  height='285px', 
                                           title=str(iImageName),  label="Image "+str(iImageName))
            firstOne = False
