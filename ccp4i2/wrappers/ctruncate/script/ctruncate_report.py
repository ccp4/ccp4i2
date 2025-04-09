import math
import os
import sys

from ....report.CCP4ReportParser import GenericElement, Report


class ctruncate_report(Report):
    TASKNAME='ctruncate'
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=None,**kw)
        
        # 'nooutput' mode would be used by another report class that wanted
        # to use some method(s) from this class for its own report
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
       
        try:
          fail = self.Errors(self)
        except:
          fail = False

        wDiv = self.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")

        if 'outputfiles' in jobInfo and len(jobInfo['outputfiles'])>0:
            wDiv.addText(text='Intensity data has been converted to structure factors')
        else:
            wDiv.addText(text='No structure factor data has been created')

        try:
            self.fileroot = self.jobInfo['fileroot']
        except:
            self.fileroot = None
        self.displayFile(self.fileroot, wDiv, ['job_3/log.txt', './log.txt'], 'Show log file')

        if fail == False:
            self.addWarnings(self)
            self.addReflectionData(self)
            self.addQuality(self)
            self.addAnomalous(self)
            self.addNCS(self)
            self.addAnisotropy(self)
            self.addIceRings(self)
            self.addTwin(self)
            self.addCtruncateGraphs(self)

    def displayFile(self, fileroot, parent, filenames, text):
        ''' display message with link to open file
        fileroot   root of file names, None if not set
        parent     where to put it
        filenames  list of possible names
        text       to display
        '''
        p = GenericElement('p')
        filefound = False
        fnames = filenames

        if fileroot != None:
            # prepend fileroot on list
            fnames = []
            for filename in filenames:
                fnames.append(fileroot+filename)

        for filename in fnames:
            if os.path.isfile(filename):
                filefound = True
                break

        if filefound:
            p.append(GenericElement('a',text,href=filename))
            parent.append(p)

    def keyText(self, parent, size='+0', phasertncs=None, phasertwin=None):
        # The most important warnings for the Aimless pipeline
        self.addWarningTwin(parent,size, phasertwin=phasertwin)
        self.addWarningTNCS(parent,size, phasertncs=phasertncs)
        self.addWarningAnisotropy(parent,size)
        self.addWarningQuality(parent,size)
        self.addWarningWilson(parent,size)
        self.addWarningIce(parent,size)

    def addWarningWilson(self,parent,size='+2',check=False):
        if len(self.xmlnode.findall("DataStatistics/Comment[@id='CompletenessQuality']"))>0:
            badrings = False;
            x = self.xmlnode.findall("OutlierRingsAnalysis/present")[0].text
            if check:
                if x.strip() != "no": 
                    return False
                return True
            if x.strip() != "no":
                if len(self.xmlnode.findall("OutlierRingsAnalysis/percentage"))>0:
                    y = float(self.xmlnode.findall("OutlierRingsAnalysis/percentage")[0].text)
                    if y > 3.0:
                        parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: Severe deviation from Wilson plot, "+str(y)+"% of bins deviate, please look at it</font></b></p>")
                    else:
                        parent.append("<p><b><font size='"+size+"'  color='orange'>&#9679; Warning: Some deviation from Wilson plot, "+str(y)+"% of bins deviate, please look at it</font></b></p>")
                else:
                    parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: Deviation in WIlson plot, please look at it</font></b></p>")
                return True
        elif check:
            return False

            
    def addWarningQuality(self,parent,size='+2',check=False):
        if len(self.xmlnode.findall("DataStatistics/Comment[@id='CompletenessQuality']"))>0:
            x = self.xmlnode.findall("DataStatistics/Comment[@id='CompletenessQuality']")[0]
            if check:
                # Just check for problem and return True if there is
                if ("poor" in x.text) or "some" in x.text:
                    return True
                return False

            if "poor" in x.text:
                parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning:"+x.text+"</font></b></p>")
            elif "some" in x.text:
                parent.append("<p><b><font size='"+size+"' color='orange'>&#9679; Warning:"+x.text+"</font></b></p>")
            else:
                parent.append("<p><b><font size='"+size+"' color='green'>&#10004; "+x.text+"</font></b></p>")
        elif check:
            return False

    def addWarningTNCS(self,parent,size='+2',check=False, phasertncs=None):
        if check:
            # Just check for problem and return True if there is
            if len(self.xmlnode.findall('translationalNCS/detected'))>0:
                if self.xmlnode.findall('translationalNCS/detected')[0].text != "No":
                    return True
            return False

        if len(self.xmlnode.findall('translationalNCS/detected'))>0:
            if self.xmlnode.findall('translationalNCS/detected')[0].text != "No":
                s = 'Warning: Possible translational non-crystallographic symmetry at<br/>'
                peaks = self.xmlnode.findall('translationalNCS/PeakList')
                #print "peaks", peaks, type(peaks)
                coordinates = []
                Qscores = []
                peakRatios = []
                for peak in peaks:
                    #print "peak",peak, type(peak)
                    coordinates.append(peak.findall('Peak/coordinate')[0].text)
                    Qscores.append(peak.findall('Peak/Q-score')[0].text)
                    peakRatios.append(peak.findall('Peak/height')[0].text)
                first = True
                for cdQ in zip(coordinates, Qscores, peakRatios):
                    if not first:
                        s += ', at '
                        first = False
                    s += 'coordinate ('+cdQ[0]+'), chance probability '+cdQ[1]+\
                         ', height '+cdQ[2]+' of origin'

                s = "<font size='"+size+"' color='red'>&#10007; "+s+"</font>"
                #parent.append("<p><b><font size='"+size+"' color='red'>&#10007; "+s+"</font></b>")

                # Did Phaser agree?
                #print("Ph tNCS", phasertncs)
                if phasertncs is not None:
                    s1 = ''
                    if phasertncs == 'True':
                        s1 += "<br/>Phaser also finds possible tNCS"
                        colour = 'red'
                    elif phasertncs == 'False':
                        s1 += "<br/>Note that Phaser does NOT find possible tNCS"
                        colour = 'blue'
                    if s1 != '':
                        s += "<font size='"+size+"' color='"+colour+"'>"+s1+"</font>"
            else:
                s = \
                  "<font size='"+size+"' color='green'> &#10004; No evidence of possible translational non-crystallographic symmetry</font>"
                if phasertncs is not None:
                    s1 = ''
                    if phasertncs == 'True':
                        s1 = " <br/> Note however that Phaser DOES find possible tNCS"
                        colour = 'red'
                    elif phasertncs == 'False':
                        s1 = " <br/> Phaser also finds no evidence of tNCS"
                        colour = 'green'
                    if s1 != '':
                        s += "<font size='"+size+"' color='"+colour+"'>"+s1+"</font>"

            parent.append("<p><b>"+s+"</b></p>")


    def addWarningTwin(self,parent,size='+2',check=False, phasertwin=None):
        if check:
            # Just check for problem and return True if there is
            if len(self.xmlnode.findall('Twinning/L-test/Twinned'))>0:
                if self.xmlnode.findall('Twinning/L-test/Twinned')[0].text != "No":
                    return True
            return False

        if len(self.xmlnode.findall('Twinning/L-test/Twinned'))>0:
            alphaBritton = 0.0
            alphaHtest = 0.0
            if len(self.xmlnode.findall('Twinning/TwinOps'))>0:
                oplist = self.xmlnode.findall("Twinning/TwinOps/Operator")
                if len(oplist)>0:
                    # Twin fraction estimate from ML-Britton
                    alist = []
                    y = self.xmlnode.findall("Twinning/ML-Britton/TwinFraction")
                    for x in y:
                        alist.append(float(x.findall('.')[0].text))
                    alphaBritton = max(alist)
                    # Twin fraction estimate from H-test
                    alist = []
                    y = self.xmlnode.findall("Twinning/H-test/TwinFraction")
                    for x in y:
                        alist.append(float(x.findall('.')[0].text))
                    alphaHtest = max(alist)
                        
            if self.xmlnode.findall('Twinning/L-test/Twinned')[0].text != "No":
                # Apperently twinned
                if alphaBritton > 0.1 or alphaHtest > 0.1:
                    s = 'Warning: Possible twinning, twin fraction estimates from Britton plot '+\
                        str(alphaBritton)+', from H-test '+str(alphaHtest)
                    parent.append("<p><b><font size='"+size+"' color='red'>&#10007; "+s+"</font></b></p>")
                else:
                    alphaLtest = self.xmlnode.findall("Twinning/L-test/TwinFraction")[0].text
                    s = 'Warning: Possibly twinned and merged in a higher symmetry than the true spacegroup'
                    s += ', twin fraction estimate from L-test ' + alphaLtest
                # Did Phaser agree? ctruncate thinks it is twinned
                if phasertwin is not None:
                    s1 = ''
                    if phasertwin > 0:
                        # Phaser agrees
                        s1 += "<br/>Phaser also finds possible twinning"
                        colour = 'red'
                    else:
                        s1 += "<br/>Note that Phaser does NOT find possible twinning"
                        colour = 'blue'
                    if s1 != '':
                        s += "<font size='"+size+"' color='"+colour+"'>"+s1+"</font>"
                        
                parent.append(
                    "<p><b><font size='"+size+"' color='red'>&#10007;"+s+"</font></b></p>")
            else:
                # No twinning
                s = "No evidence of twinning"
                # Did Phaser agree? ctruncate thinks it is NOT twinned
                if phasertwin is not None:
                    s1 = ''
                    if phasertwin > 0:
                        # Phaser disagrees
                        s1 += "<br/>But Phaser DOES find possible twinning"
                        colour = 'red'
                    else:
                        s1 += "<br/>Phaser also finds no evidence of twinning"
                        colour = 'blue'
                    if s1 != '':
                        s += "<font size='"+size+"' color='"+colour+"'>"+s1+"</font>"
                parent.append("<p><b><font size='"+size+"' color='green'>&#10004; "+s+"</font></b></p>")

    def addWarningAnisotropy(self,parent,size='+2',check=False):
        if len(self.xmlnode.findall('AnisotropyAnalysis'))>0:
            if len(self.xmlnode.findall('AnisotropyAnalysis/Eigenvalue'))>0:
                paths = self.xmlnode.findall('AnisotropyAnalysis/Eigenvalue')
                y = []
                for x in paths:
                    y.append(x.findall('.')[0].text)
                ratio = min(list(map(float,y)))/max(list(map(float,y)))
                if check:
                    # Just check for problem and return True if there is
                    if ratio < 0.9:
                        return True
                    return False

                if ratio > 0.9:
                    parent.append("<p><b><font size='"+size+"' color='green'>&#10004; Little or no anisotropy detected</font></b></p>")
                elif ratio > 0.5:
                    parent.append("<p><b><font size='"+size+"' color='orange'>&#9679; Warning: Some anisotropy detected. This may have an effect on statistics.</font></b></p>")
                else:
                    parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: data are highly anisotropic.</font></b></p>")

            elif (len(self.xmlnode.findall('AnisotropyAnalysis/Comment'))>0):
                message = self.xmlnode.findall('AnisotropyAnalysis/Comment')[0].text
                #print "message", message
                parent.append("<font size='"+size+"'>"+message+"</font>")

        elif check:
            return False

    def addWarningIce(self,parent,size='+2',check=False):
        if len(self.xmlnode.findall('IceRingsAnalysis'))>0:
            ringslist = self.xmlnode.findall('IceRingsAnalysis/Ring/Reject')
            #better in xml
            icerings = False
            for x in ringslist:
                y = x.findall('.')[0].text
                if y.strip() != "no": icerings = True

            if check:
                # Just check for problem and return True if there is
                if icerings is True:
                    return True
                return False

            if icerings is True:
                parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: Possible ice rings found.</font></b></p>")
            else:
                parent.append("<p><b><font size='"+size+"' color='green'>&#10004; No ice rings found.</font></b></p>")
        elif check:
            return False

    def addWarningAnomalous(self,parent,size='+2',check=False):
        if len(self.xmlnode.findall("AnomStatistics/Comment[@id='AnomResult']"))>0:
            x=self.xmlnode.findall("AnomStatistics/Comment[@id='AnomResult']")[0]
            if "Warning" in x.text:
                parent.append("<p><b><font size='"+size+"' color='red'>&#10007; "+x.text+"</font></b></p>")
            elif "Some" in x.text:
                parent.append("<p><b><font size='"+size+"' color='orange'>&#9679; "+x.text+"</font></b></p>")
            else:
                parent.append("<p><b><font size='"+size+"' color='green'>&#10004; "+x.text+"</font></b></p>")
        elif len(self.xmlnode.findall("AnomStatistics/ResolutionRange[@id='Measurability Limit']"))>0:
            f = self.xmlnode.findall("AnomStatistics/ResolutionRange[@id='Measurability Limit']")[0].text
            xmax = float(f.findall('max')[0].text)
            xmin = float(f.findall('min')[0].text)
            f = self.xmlnode.findall("AnomStatistics/ResolutionRange[@id='Signal to Noise Ratio']")[0].text
            ymax = float(f.findall('max')[0].text)
            ymin = float(f.findall('min')[0].text) 
            if math.isnan(xmin) or math.isnan(ymin):
                parent.append("<p><b><font color='red'>&#10007; Warning: no anomalous signal.</font></b></p>")
            elif xmin == xmax and ymin == ymax:
                parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: no anomalous signal.</font></b></p>")
            elif xmin > 10 and xmax < 3  and  ymin > 10 and ymax < 3:
                parent.append("<p><b><font size='"+size+"' color='green'>&#10004; Significant anomalous signal to high resoluiton.</font></b></p>")
            elif  xmin > 10 and xmax < 5  and  ymin > 10 and ymax < 5:
                parent.append("<p><b><font color='green'>&#10004; Significant anomalous signal to medium resoluiton.</font></b></p>")
            elif ( xmin > 10 and xmax < 5 ) or ( ymin > 10 and ymax < 5 ):
                parent.append("<p><b><font size='"+size+"' color='orange'>&#9679; Inconsistent results, there may be some anomalous signal.</font></b></p>")
            elif xmin < 5 or ymin < 5:
                parent.append("<p><b><font size='"+size+"' color='red'>&#10007; Warning: No anomalous signal. Anomalous signal is probably artefact of noise at high resolution.</font></b></p>")
            else:
                parent.append("<p><b><font size='"+size+"' color='orange'>&#9679; Some anomalous signal at intermediate resolution.  Low resolution signal missing.</font></b></p>")

    def addWarnings(self, parent):
        #summary and warnings
        wDiv = parent.addDiv(style="border-width: 1px; border-color: black; clear:both; margin:0px; padding:0px;")
        wDiv.addText(text="Summary",style="font-size: 180%;text-align:left")
        #eg tNCS
        self.addWarningQuality(wDiv)
        self.addWarningWilson(wDiv)
        self.addWarningTNCS(wDiv)
        self.addWarningTwin(wDiv)
        self.addWarningAnisotropy(wDiv)
        self.addWarningIce(wDiv)
        self.addWarningAnomalous(wDiv)
            
    def CtruncateLtest(self, parent):
        # A GraphGroup is a group of graphs displayed in the same graph viewer widget
        graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:300px;margin: 0 auto;border:1px solid black;")
        # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        # Loop over all Graph tables in the program output and add to the GraphGroup
        # The plotting instructions are provided as xml text
        graphlist = self.xmlnode.findall("Twinning/CCP4Table[@id='L-test']")
        #print "CtruncateLtest", self.xmlnode
        #print "CtruncateLtest", graphlist
        for thisgraph in graphlist:
            graph = graphgroup.addFlotGraph(xmlnode=thisgraph,
                                            title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)

    def CtruncateWilson(self, parent):
        # A GraphGroup is a group of graphs displayed in the same graph viewer widget
        graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:270px;margin: 0 auto;border:1px solid black;")
        # Add a Graph to the GraphGroup - add a table of data and plot instructions to the graph
        # Loop over all Graph tables in the program output and add to the GraphGroup
        # The plotting instructions are provided as xml text
        graphlist = self.xmlnode.findall("CCP4Table[@id='WilsonB']")
        #print("CtruncateWilson", len(graphlist), self.xmlnode)
        for thisgraph in graphlist:
            graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)

    def addCtruncateReport(self, parent, ndatasets=1):
        # The detailed report for one dataset
        datasetid = self.xmlnode.findall('ReflectionFile/CrystalDatasetId')[0].text
        # print("**Dataset", datasetid)
        reportDivHeader = parent.addDiv(style="font-size:120%;")
        reportDivHeader.append('<span style="font-style: bold">Detailed reports on intensity statistics etc</span><br/>')
        #reportDivHeader.append('<a href="#twinreport">* Twinning</a><br/>')
        #reportDivHeader.append('<a href="#dataqualityreport">* Data quality</a><br/>')
        #reportDivHeader.append('<a href="#ncsreport">* Non-crystallographic symmetry</a><br/>')
        #reportDivHeader.append('<a href="#anomalousreport">* Anomalous signal</a><br/>')
        #reportDivHeader.append('<a href="#anisotropyreport">* Anisotropy</a><br/>')
        #reportDivHeader.append('<a href="#iceringsreport">* Ice rings</a><br/>')

        if ndatasets > 1:
            label = "Statistics for dataset: "+datasetid
            startopen = False
            dsetstats = parent.addFold(label=label, brief='datasetid',
                                      initiallyOpen=startopen)
        else:
            dsetstats = parent

        reportDiv = dsetstats.addDiv()
        self.addTwin(reportDiv)
        self.addQuality(reportDiv)
        self.addNCS(reportDiv)
        self.addAnomalous(reportDiv)
        self.addAnisotropy(reportDiv)
        self.addIceRings(reportDiv)

    def addReflectionData(self,parent=None):
        #parent.append('<a name="reflectiondata"></a>')
        startopen = True
        fold = parent.addFold(label="Details of reflection data",
                              brief='InputData',initiallyOpen=startopen)
        fold.addText(text="Summary of input reflection data",style="font-size:130%;text-align:center")
        fold.append("<br/>")
        f = self.xmlnode.findall('ReflectionFile')
        fold.append("<b>File:</b> "+f[0].get("name"))
        fold.append("<b>Crystal/Dataset:</b> "+f[0].findall('CrystalDatasetId')[0].text)
        fold.append("<b>Cell:</b> "+f[0].findall('cell/a')[0].text+", "+f[0].findall('cell/b')[0].text+", "+f[0].findall('cell/c')[0].text+", "+f[0].findall('cell/alpha')[0].text+", "+f[0].findall('cell/beta')[0].text+", "+f[0].findall('cell/gamma')[0].text)
        fold.append("<b>Spacegroup:</b> "+f[0].findall('SpacegroupName')[0].text)
        fold.append("<p><b>Reflection data summary:</b></p>")
        table = fold.addTable( select = "ReflectionData", transpose=False, class_="center" )
        taglist = \
                    [ [ "Data type", "ReflectionDataType" ],
                      [ "Min resolution", "ResolutionLow"],
                      [ "Max resolution", "ResolutionHigh" ],
                      [ "Nreflections", "NumberReflections" ],
                      [ "NObservations", "NumberObservations" ],
                      [ "NCentric", "NumberCentric" ] ]
        self.addToTable(taglist, "ReflectionData", table)
        
        fold.append("<p><b>Maximum indices in reflection data:</b></p>")
        path=self.xmlnode.findall('ReflectionData/ReflectionIndexMax')
        table1 = fold.addTable( select = "ReflectionData", transpose=False, class_="center" )
        tablist = []
        for x in path:
            tablist.append(x.findall("[@id]")[0].attrib["id"])
        table1.addData(title='Direction',data=tablist)
        taglist = \
            [["Index (file)", "FileReflectionIndexMax/Index"],
             ["Max Resolution (file)", "FileReflectionIndexMax/Resolution"],
             ["Index (P1)", "ReflectionIndexMax/Index"],
             ["Max Resolution (P1)", "ReflectionIndexMax/Resolution"]]
        self.addToTable(taglist, "", table1)

    def addNCS(self, parent=None):
        if len(self.xmlnode.findall('translationalNCS'))>0:
            #parent.append('<a name="ncsreport"></a>')
            startopen = self.addWarningTNCS(parent, check=True)
            path=self.xmlnode.findall('translationalNCS')
            fold = parent.addFold(label="tNCS",initiallyOpen=startopen)
            fold.addText(text='tNCS',style="font-size:130%;")
            x = self.xmlnode.findall("translationalNCS/Comment[@id='tncs']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            x = self.xmlnode.findall("translationalNCS/Comment[@id='tncsExplain']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            x = self.xmlnode.findall("translationalNCS/Reference[@id='tncsExplain']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            if path[0].findall('detected')[0].text != "No": fold.append("<p><b>Warning:</b> flat prior used for truncate procedure</p>")
            
    def addAnisotropy(self, parent=None):
        if len(self.xmlnode.findall('AnisotropyAnalysis'))>0 and \
               len(self.xmlnode.findall("AnisotropyAnalysis/Comment[@id='AnisotropyResult']"))>0:
            #parent.append('<a name="anisotropyreport"></a>')
            startopen = self.addWarningAnisotropy(parent, check=True)
            path=self.xmlnode.findall('AnisotropyAnalysis')
            fold = parent.addFold(label="Anisotropy analysis",brief="Anisotropy",
                                  initiallyOpen=startopen)
            fold.addText(text='Anisotropy Analysis',style="font-size:130%;")
            x = self.xmlnode.findall("AnisotropyAnalysis/Comment[@id='AnisotropyResult']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            #fold.append("Anisotropy will be corrected for in the truncation procedure.")
            if len(self.xmlnode.findall('AnisotropyAnalysis/BAniso'))>0:
                fold.append("<p><b>ML Anisotropy matrices:</b></p>")
                aDiv = fold.addDiv(style="width:50%;float:left;border-width: 0px; border-color: black;text-align:left; margin:0px; padding:0px;")
                bDiv = fold.addDiv(style="width:50%;float:left;border-width: 0px; border-color: black;text-align:left; margin:0px; padding:0px;")
                result = path[0].findall("BAniso/Orthogonal")[0].text
                list = result.split()
                list.insert(6,'<br/>')
                list.insert(3,'<br/>')
                aDiv.append("<p><b>Anisotropy B matrix:</b></p>"+"  ".join(str(x) for x in list))
                result = path[0].findall("UAniso[@id='obs']/Orthogonal")[0]
                list = result.text.split()
                list.insert(6,'<br/>')
                list.insert(3,'<br/>')
                bDiv.append("<p><b>Anisotropy U matrix:</b></p>"+"  ".join(str(x) for x in list))
                fold.append("<p><b>Eigenvalues and Eigenvectors of anisotropy matrix</b></p>")
                table = fold.addTable( select = "AnisotropyAnalysis", transpose=False, class_="center" )
                eigenvaluedata = []
                eigenvectordata = []
                eigenvalues = self.xmlnode.findall('AnisotropyAnalysis/Eigenvalue')
                eigenvectors = self.xmlnode.findall('AnisotropyAnalysis/EigenVector')
                for i in range(3):
                    eigenvaluedata.append(eigenvalues[i].text)
                    eigenvectordata.append(eigenvectors[i].text)
                table.addData( title="Eigenvalue", data=eigenvaluedata)
                table.addData( title="Eigenvector", data=eigenvectordata)
                
                fold.append("<p><b>Directional plots:</b></p>")
                x = self.xmlnode.findall("AnisotropyAnalysis/Comment[@id='YorgoModisPlot']")[0]
                graphlist = self.xmlnode.findall("CCP4Table[@id='AnisotropyAnalysis']")
                if len(graphlist) > 0:
                    # not stable, so should use graph-aniso as grouID
                    graph = fold.addFlotGraph( xmlnode=graphlist[0], title=graphlist[0].get("title") )
                    graph = graph.addPimpleData(xmlnode=graphlist[0])
                fold.append("<p><font color='black'>"+x.text+"</font></p>")
                result = path[0].findall("UAniso[@id='correction']/Orthogonal")[0]
                list = result.text.split()
                list.insert(6,'<br/>')
                list.insert(3,'<br/>')
                fold.append("<p><b>Anisotropic correction:</b></p>"+"  ".join(str(x) for x in list))
                fold.append("<p>This will be used in the remaining statistics calculations.</p>")
                #fold.append("<p><font color='red'>effect of correction on data.</font></p>")
#######

    def addQuality(self, parent=None):
        #completeness table name...
        #parent.append('<a name="dataqualityreport"></a>')
        startopen = self.addWarningQuality(parent, check=True)
        fold = parent.addFold(label="Data quality",brief="Quality",initiallyOpen=startopen)
        fold.addText(text='Data quality',style="font-size:130%;")
        if len(self.xmlnode.findall('DataStatistics/ResolutionRange'))>0:
            path=self.xmlnode.findall('DataStatistics/ResolutionRange')
            fold.append("<p><b>Resolution range statistics:</b></p>")
            table = fold.addTable( select = ".//DataStatistics/ResolutionRange", transpose=False, style="margin: 0 auto;border:1px solid orange;", class_="center" )
            tablist = []
            minlist = []
            maxlist = []
            percentlist = []
            for x in path:
                try:
                    tablist.append(x.findall("[@title]")[0].attrib["title"])
                except:
                    tablist.append("Untitled")
                minlist.append(x.findall('min')[0].text)
                maxlist.append(x.findall('max')[0].text)
                if len(self.xmlnode.findall('DataStatistics/ResolutionRange/percentage'))>0:
                    percentlist.append(x.findall('percentage'))
            table.addData(title='Statistic',data=tablist)
            table.addData(title="Min Resolution", data=minlist)
            table.addData(title="Max Resolution", data=maxlist)
            if len(self.xmlnode.findall('DataStatistics/ResolutionRange/percentage'))>0:
                table.addData(title="Percentage refln", data=percentlist)

            x1 = self.xmlnode.findall("DataStatistics/Comment[@id='CompletenessReso']")[0]
            x2 = self.xmlnode.findall("DataStatistics/Comment[@id='Completeness3']")[0]
            x3 = self.xmlnode.findall("DataStatistics/Comment[@id='ActiveRange']")[0]
            fold.append("<p><font color='black'>"+x1.text+" "+x2.text+"<ul>"+x3.text+"</ul></font></p>")

        if len(self.xmlnode.findall("CCP4Table[@id='completeness']"))>0:
            fold.append("<p><b>Completeness and R-standard plots</b></p>")
            graphlist = self.xmlnode.findall("CCP4Table[@id='completeness']")
            #not stable, so should use graph-aniso as groupID
            graph = fold.addFlotGraph( xmlnode=graphlist[0], title=graphlist[0].get("title") )
            graph = graph.addPimpleData(xmlnode=graphlist[0])
            x = self.xmlnode.findtext("DataStatistics/Comment[@id='Completeness']")
            if x is not None:
                x = x.replace("<", "&lt;");
                x = x.replace(">", "&gt;");
                fold.append("<p><font color='black'>"+x+"</font></p>")
        if len(self.xmlnode.findall("DataStatistics/LowResoCompleteness"))>0:
            fold.append("<p><b>Low Resolution Completeness</b></p>")
            table = fold.addTable( select = "DataStatistics/LowResoCompleteness", transpose=False, style="margin: 0 auto;border:1px solid orange;", class_="center" )
            resolist1 = self.xmlnode.findall("DataStatistics/LowResoCompleteness/ResolutionRange/min")
            resolist2 = self.xmlnode.findall("DataStatistics/LowResoCompleteness/ResolutionRange/max")
            complist1 = self.xmlnode.findall("DataStatistics/LowResoCompleteness/Completeness")
            complist2 = self.xmlnode.findall("DataStatistics/LowResoCompleteness/nObs")
            complist3 = self.xmlnode.findall("DataStatistics/LowResoCompleteness/nTot")
            resotab = []
            comptab = []
            for x,y in zip(resolist1,resolist2):
                resotab.append(x.findall(".")[0].text+" - "+y.findall(".")[0].text )
            for x,y,z in zip(complist1,complist2,complist3):
                comptab.append(x.findall(".")[0].text+"  ["+y.findall(".")[0].text+" : "+z.findall(".")[0].text+"]" )
            table.addData(title='Resolution',data=resotab)
            table.addData(title='Completeness',data=comptab)
            x = self.xmlnode.findall("DataStatistics/Comment[@id='CompletenessLow']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")

        if len(self.xmlnode.findall('DataStatistics/WilsonB'))>0:
            path=self.xmlnode.findall('DataStatistics/WilsonB')
            fold.append("<p><b>Wilson scaling:</b></p>")
            fold.append("<p><font color='black'>Estimate of overall wilson B-factor: "+self.xmlnode.findall('DataStatistics/WilsonB')[0].text+" A^(-2)</font></p>")
            fold.append("<p><font color='black'>Estimate wilson scale factor: "+self.xmlnode.findall('DataStatistics/WilsonScale')[0].text+"</font></p>")
            x = self.xmlnode.findall("DataStatistics/Comment[@id='WilsonB']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            if len(self.xmlnode.findall("CCP4Table[@id='WilsonB']"))>0:
                graphlist = self.xmlnode.findall("CCP4Table[@id='WilsonB']")
                #not stable, so should use graph-aniso as grouID
                graph = fold.addFlotGraph( xmlnode=graphlist[0], title=graphlist[0].get("title") )
                graph = graph.addPimpleData(xmlnode=graphlist[0])
            x = self.xmlnode.findall("DataStatistics/Comment[@id='WilsonPlot']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            x = self.xmlnode.findall("DataStatistics/Reference[@id='Best']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            #fold.append("<p><font color='red'>ML calc<br/>outlier reflections</font></p>")

        if len(self.xmlnode.findall('OutlierRingsAnalysis'))>0:
            fold.append("<p><b>Check for bad rings:</b></p>")
            badrings = False;
            x = self.xmlnode.findall("OutlierRingsAnalysis/present")[0].text
            if x.strip() != "no": badrings = True
            
            if badrings is True:
                if len(self.xmlnode.findall("OutlierRingsAnalysis/percentage"))>0:
                    y = float(self.xmlnode.findall("OutlierRingsAnalysis/percentage")[0].text)
                    fold.append("<p>Possible bad resolution rings found, consisting of "+str(y)+"% of bins. Data in these may be of lower quality.</p>")
                else:
                    fold.append("<p>Possible bad resolution rings found. Data in these may be of lower quality.</p>")
                table = fold.addTable( select = "OutlierRingsAnalysis", transpose=False, style="margin: 0 auto;border:1px solid orange;", class_="center" )
                taglist = \
                        [["Resolution", "Resolution"],
                         ["Imean", "Imean"],
                         ["Z-score", "Z-score"],
                         ["Completeness", "Completeness"],
                         ["Ave. Completeness", "ExpectCompleteness"]]
                self.addToTable(taglist,"OutlierRingsAnalysis/Ring", table)
            else:
                fold.append("<p>No issues found.</p>")
            
            x = self.xmlnode.findall("OutlierRingsAnalysis/Comment[@id='OutlierRingsAnalysis']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")

    def addAnomalous(self, parent=None):
        if len(self.xmlnode.findall('AnomStatistics'))>0:
            #parent.append('<a name="anomalousreport"></a>')
            startopen = False
            fold = parent.addFold(label="Anomalous",initiallyOpen=startopen)
            fold.addText(text='Anomalous signal',style="font-size:130%;")
            #fold.append("<p><font color='red'>Use resolution range to comment on anom signal</font></p>")
            path=self.xmlnode.findall('AnomStatistics/ResolutionRange')
            fold.append("<p><b>Anomalous signal statistics:</b></p>")
            table = fold.addTable( select = "AnomStatistics/ResolutionRange", transpose=False, style="margin: 0 auto;border:1px solid orange;", class_="center" )
            tablist = []
            minlist = []
            maxlist = []
            for x in path:
                tablist.append(x.findall("[@id]")[0].attrib["id"])
                minlist.append(x.findall('min')[0].text)
                maxlist.append(x.findall('max')[0].text)
            table.addData(title='Statistic',data=tablist)
            table.addData(title='Min Resolution',data=minlist)
            table.addData(title='Max Resolution',data=maxlist)
            if len(self.xmlnode.findall("AnomStatistics/Comment[@id='Anom']"))>0:
                x=self.xmlnode.findall("AnomStatistics/Comment[@id='Anom']")[0]
                fold.append("<p><font color='black'>"+x.text+"</font></p>")
                fold.append("<p><font color='green'>The resolution limits are derived from I/sigI of 1.2 for the signal to noise, and measurability above 5%.</font></p>")
            if len(self.xmlnode.findall("AnomStatistics/Comment[@id='AnomResult']"))>0:
                x=self.xmlnode.findall("AnomStatistics/Comment[@id='AnomResult']")[0]
                fold.append("<p><font color='black'>"+x.text+"</font></p>")
                x=self.xmlnode.findall("AnomStatistics/Comment[@id='AnomDescription']")[0]
                fold.append("<p><font color='black'>"+x.text+"</font></p>")
            if len(self.xmlnode.findall("CCP4Table[@id='anomalous intensity plot']"))>0:
                graphlist = self.xmlnode.findall("CCP4Table[@id='anomalous intensity plot']")
                #not stable, so should use graph-aniso as groupID
                graph = fold.addFlotGraph( xmlnode=graphlist[0], title=graphlist[0].get("title") )
                graph = graph.addPimpleData(xmlnode=graphlist[0])

    def addIceRings(self, parent=None):
        if len(self.xmlnode.findall('IceRingsAnalysis'))>0:
            #parent.append('<a name="iceringsreport"></a>')
            startopen = self.addWarningIce(parent, check=True)
            fold = parent.addFold(label="Ice Rings",initiallyOpen=startopen)
            fold.addText(text='Ice Rings',style="font-size:130%;")
            icerings = False;
            ringslist = self.xmlnode.findall('IceRingsAnalysis/Ring/Reject')
            #better in xml
            for x in ringslist:
                y = x.findall('.')[0].text
                if y.strip() != "no": icerings = True
            
            if icerings is True:
                fold.append("<p>Possible ice rings found. Data in these may be of lower quality.</p>")
                fold.append("Ice rings will be ignored in the truncate procedure.")
            else:
                fold.append("<p>No ice rings found.</p>")
            
            table = fold.addTable( select = "IceRingsAnalysis", transpose=False, style="margin: 0 auto;border:1px solid orange;", class_="center" )
                    
            #tablelist = self.xmlnode.findall("Ring")
            taglist = \
                    [["Resolution", "Resolution"],
                     ["Ice Ring", "Reject"],
                     ["Imean", "Imean"],
                     ["Z-score", "Z-score"],
                     ["Completeness", "Completeness"],
                     ["Ave. Completeness", "ExpectCompleteness"]]

            self.addToTable(taglist,"OutlierRingsAnalysis/Ring", table)
            x = self.xmlnode.findall("IceRingsAnalysis/Comment[@id='IceRingsAnalysis']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")

    def addTwin(self, parent=None):
        havetwinblock = False
        if len(self.xmlnode.findall('Twinning'))>0:
            havetwinblock = True
        else:
            return
        
        # Twinning block for this dataset
        twinblock = self.xmlnode.findall('Twinning')[0]

        #parent.append('<a name="twinreport"></a>')
        # Should this be initially open?
        startopen = self.addWarningTwin(parent, check=True)
        fold = parent.addFold(label="Twinning statistics",
                              brief="Twinning",initiallyOpen=startopen)
        fold.addText(text='Twinning statistics',style="font-size:130%;")
        x=twinblock.findall("Comment[@id='TwinningSummary']")[0]
        fold.append("<p><font color='black'>"+x.text+"</font></p>")
        x1 = x.text
        y1 = x1.strip()
        if len(y1) > 0 and y1[0] != "N": fold.append("<p><b>Warning:</b> flat prior used for truncate procedure</p>")
        if len(twinblock.findall('L-test'))>0:
            fold.append("<p><b>Global twinning estimates:</b></p>")
            x=twinblock.findall("Comment[@id='TwinningGlobal']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
            table=fold.addTable(class_="center")
            testtab = []
            resulttab = []
            scoretab = []
            alphatab = []
            if len(twinblock.findall('L-test'))>0:
                testtab.append('L-test')
                resulttab.append(twinblock.findall('L-test/Twinned')[0].text)
                scoretab.append(twinblock.findall('L-test/Result')[0].text)
                alphatab.append(twinblock.findall('L-test/TwinFraction'))
            table.addData(title='Test',data=testtab)
            table.addData(title='Twinnned?',data=resulttab)
            table.addData(title='Score',data=scoretab)
            table.addData(title='alpha',data=alphatab)
                    
            x = twinblock.findall(".//Comment[@id='L-test']")[0]
            fold.append("<p>"+x.text+"</p>")
                    
            if len(self.xmlnode.findall('TwinOps'))>0:
                fold.append("<p><b>Operator based twinning estimates:</b></p>")
                if len(twinblock.findall("Rtwin"))>0:
                    x = twinblock.findtext("Comment[@id='TwinningOperatorBased']")
                    x = x.replace("<", "&lt;");
                    x = x.replace(">", "&gt;");
                    fold.append("<p><font color='black'>"+x+"</font></p>")
                    oplist = twinblock.findall("TwinOps/Operator")
                    if len(oplist)>0:
                        rtwinlist = twinblock.findall("Rtwin")
                        rhlist = twinblock.findall("H-test")
                        rmllist = twinblock.findall("ML-Britton")
                        table=fold.addTable(class_="center")
                        optab = []
                        rtwintab = []
                        rhtab = []
                        rmltab = []
                        for x in rtwinlist:
                            # optab.append(x.get("operator") )
                            optab.append(x.findall("[@operator]")[0].attrib["operator"])
                            rtwintab.append(x.findall(".")[0].text )
                        for x in rhlist:
                            rhtab.append(x.findall('TwinFraction')[0].text)
                        for x in rmllist:
                            rmltab.append(x.findall('TwinFraction')[0].text)
                        table.addData(title='Operator',data=optab)
                        table.addData(title='Rtwin</br>score',data=rtwintab)
                        table.addData(title='H-test</br>Twin Fraction',data=rhtab)
                        table.addData(title='ML-Britton</br>Twin Fraction',data=rmltab)
                else:
                    fold.append("<p>"+twinblock.findall("Comment[@id='TwinOps']")[0].text+"</p>")

        fold.append("<p><b>Twinning plots based on global statistics</b></p>")
        fold.append("<p><b><i>L-test</i></b></p>")
        ltestlist = twinblock.findall("CCP4Table[@id='L-test']")
        for thisgraph in ltestlist:
            graph = fold.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)
        x = twinblock.findall(".//Comment[@id='LTestDesc']")[0]
        fold.append("<p><font color='black'>"+x.text+"</font></p>")
        x = twinblock.findall(".//Reference[@id='LTest']")[0]
        fold.append("<p><font color='black'>"+x.text+"</font></p>")
        fold.append("<p><b><i>N(Z) plot</i></b></p>")
        #problem is have amplitudes too
        cumullist = twinblock.findall("CCP4Table[@id='Cumulative intensity distribution']")
        for thisgraph in cumullist:
            graph = fold.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
            graph = graph.addPimpleData(xmlnode=thisgraph)
            x = twinblock.findall(".//Comment[@id='NZ']")[0]
            fold.append("<p><font color='black'>"+x.text+"</font></p>")
        #fold.append("<p><font color='red'>any stats</font></p>")
        if len(twinblock.findall('Moments'))>0:
            fold.append("<p><b><i>Wilson ratios and moments plot</i></b></p>")
            aDiv = fold.addDiv(style="width:50%;float:left;border-width: 0px; border-color: black;text-align:left; margin:0px; padding:0px;")
            bDiv = fold.addDiv(style="width:50%;float:left;border-width: 0px; border-color: black;text-align:left; margin:0px; padding:0px;")
            #problem is have amplitudes too
            mlist = twinblock.findall("Moments[@id='acentric']/Moment")
            if len(mlist)>0:
                aDiv.append("<p>Acentric reflections</p>")
                table=aDiv.addTable(class_="center")
                optab = []
                vtab = []
                utab = []
                ptab = []
                for x in mlist:
                    y = x.findall("[@id]")[0].attrib["id"]
                    y = y.replace("<", "&lt;");
                    y = y.replace(">", "&gt;");
                    optab.append(y)
                    vtab.append(x.findall("value")[0].text)
                    ptab.append(x.findall("twinned")[0].text)
                    utab.append(x.findall("untwinned")[0].text)
                table.addData(title='Operator',data=optab)
                table.addData(title='Value',data=vtab)
                table.addData(title='Untwinned',data=utab)
                table.addData(title='Perfect</br>twin',data=ptab)
                momlist = twinblock.findall("CCP4Table[@id='acentricMoments']")
                for thisgraph in momlist:
                    graph = aDiv.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
                    graph = graph.addPimpleData(xmlnode=thisgraph)
                mlist = twinblock.findall("Moments[@id='centric']/Moment")
                if len(mlist)>0:
                    bDiv.append("<p>Centric reflections</p>")
                    table=bDiv.addTable(class_="center")
                    optab = []
                    vtab = []
                    utab = []
                    ptab = []
                    for x in mlist:
                        y = x.findall("[@id]")[0].attrib["id"]
                        y = y.replace("<", "&lt;");
                        y = y.replace(">", "&gt;");
                        optab.append(y)
                        vtab.append(x.findall("value")[0].text)
                        ptab.append(x.findall("twinned")[0].text)
                        utab.append(x.findall("untwinned")[0].text)
                    table.addData(title='Operator',data=optab)
                    table.addData(title='Value',data=vtab)
                    table.addData(title='Untwinned',data=utab)
                    table.addData(title='Perfect</br>twin',data=ptab)
                    momlist = twinblock.findall("CCP4Table[@id='centricMoments']")
                    for thisgraph in momlist:
                        graph = bDiv.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
                        graph = graph.addPimpleData(xmlnode=thisgraph)
                        
                x = twinblock.findall(".//Comment[@id='MomentsReso']")[0]
                fold.append("<p><font color='black'>"+x.text+"</font></p>")
                # fold.append("<p><font color='red'>moments of E</font></p>")

        if len(twinblock.findall('TwinOps'))>0:
            if len(twinblock.findall("Rtwin"))>0:
                oplist = twinblock.findall("TwinOps")
                #any twinops
                if len(oplist)>0:
                    fold.append("<p><b>Calculation of possible twinning operators:</b></p>")
                    table = fold.addTable( select = "TwinOps", transpose=False, class_="center" )
                    for title,select in [ [ "Operator", "Operator/Description" ],
                             [ "Score/Obliquenss", "Operator/Score" ],
                             [ "Type(m/pm)", "Operator/Type" ] ]:
                        table.addData( title=title , select = select )
                    x = twinblock.findall("TwinOps/Comment[@id='TwinOps']")[0]
                    fold.append("<p>"+x.text+"</p>")
                    x = twinblock.findall("Comment[@id='IfTwinOpsPresent']")[0]
                    fold.append("<p><font color='black'>"+x.text+"</font></p>")
                    fold.append("<p><b>Twinning operator plots:</b></p>")
                    fold.append("<p><i>H plot</i></p>")
                    x = twinblock.findall("Comment[@id='HTest']")[0]
                    fold.append("<p><font color='black'>"+x.text+"</font></p>")
                    hlist = twinblock.findall("CCP4Table[@id='H-test']")
                    for thisgraph in hlist:
                        graph = fold.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
                        graph = graph.addPimpleData(xmlnode=thisgraph)
                    fold.append("<p><i>ML Britton</i></p>")
                    x = twinblock.findall("Comment[@id='MLBritton']")[0]
                    fold.append("<p><font color='black'>"+x.text+"</font></p>")
                    mllist = twinblock.findall("CCP4Table[@id='ML-Britton']")
                    for thisgraph in mllist:
                        graph = fold.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
                        graph = graph.addPimpleData(xmlnode=thisgraph)


    def addCtruncateGraphs(self, parent):
        reportDiv = parent.addDiv()
        fold = parent.addFold(label="Graphs")
        #summaryNode = self.xmlnode.findall('CCP4Summary')[0]
        #summaryDiv = reportDiv.addDiv()
        #summaryDiv.text = '<pre>'+summaryNode.text + '</pre>'
        graphTableList = self.xmlnode.findall("CCP4Table")
        graphgroup = fold.addFlotGraphGroup(style="width:500px;  height:400px;")
        for graphTableNode in graphTableList:
            #print "CTruncateGraphs",graphTableNode
            graph = graphgroup.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title") )
            graph = graph.addPimpleData(xmlnode=graphTableNode)


    def acentricMoments(self, parent, drawgraphs=True):
        if len(self.xmlnode.findall('Twinning/Moments'))>0:
            parent.addText(text="Acentric intensity moments",
                           style="font-weight:bold; font-size:130%;")
            parent.append('<br/>')
            if drawgraphs:  # option to suppress graphs, if already done, just draw table
                momlist = self.xmlnode.findall("Twinning/CCP4Table[@id='acentricMoments']")
                # A GraphGroup is a group of graphs displayed in the same graph viewer widget
                graphgroup = parent.addFlotGraphGroup(style="width:300px;  height:300px;margin: 0 auto;border:1px solid black;")
                for thisgraph in momlist:
                    graph = graphgroup.addFlotGraph( xmlnode=thisgraph, title=thisgraph.get("title") )
                    graph = graph.addPimpleData(xmlnode=thisgraph)

            mlist = self.xmlnode.findall("Twinning/Moments[@id='acentric']/Moment")
            dsetidsnodes = self.xmlnode.findall("ReflectionFile/CrystalDatasetId")
            #print("CTdts",dsetidsnodes)
            ndts = len(dsetidsnodes)  # number of datasets
            dnames = []
            for dsn in dsetidsnodes:
                dnames.append(dsn.text.split('/')[2])
                dnames.append("")  # dummies for 3rd & 4th moments
                dnames.append("")

            if len(mlist)>0:
                s = "<p>Values for these data, and for ideal data (untwinned or twinned)"
                if (self.addWarningTNCS(parent, check=True)):
                    s += "<br/><span style='color:orange; font-style:italic;'>"
                    s += "Since translation NCS is present, these ideal data values are inappropriate</span>"
                s += "</p>"
                parent.append(s)
                    
                table=parent.addTable(class_="center")
                dntab = []
                optab = []
                vtab = []
                utab = []
                ptab = []
                for x,dn in zip(mlist,dnames):
                    y = x.findall("[@id]")[0].attrib["id"]
                    y = y.replace("<", "&lt;");
                    y = y.replace(">", "&gt;");
                    optab.append(y)
                    vtab.append(x.findall("value")[0].text)
                    ptab.append(x.findall("twinned")[0].text)
                    utab.append(x.findall("untwinned")[0].text)
                    dntab.append(dn)

                if ndts > 1: table.addData(title='Dataset',data=dntab)
                table.addData(title='Operator',data=optab)
                table.addData(title='Value',data=vtab)
                table.addData(title='Untwinned',data=utab)
                table.addData(title='Perfect</br>twin',data=ptab)

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
    def addToTable(self, taglist, nodetag, table):
        # add data to table from taglist,
        # as the AddData(select) option does not work properly
        # when there are multiple entries
        rdata = {}
        for title,tag in taglist:
            rdata[tag] = []
        if(nodetag):
            paths = self.xmlnode.findall(nodetag)
            for path in paths:
                for title,tag in taglist:
                    d = path.findall(tag)[0].text
                    rdata[tag].append(d)
        for title,tag in taglist:
            table.addData( title=title, data=rdata[tag])

    def getDatasetid(self):
        return self.xmlnode.findall('ReflectionFile/CrystalDatasetId')[0].text


############################################################################
if __name__ == "__main__":
  report = ctruncate_report(xmlFile = sys.argv[1] )
  tree= report.as_etree()
  report.as_html_file(fileName='./test-ctruncate.html')
  if len(report.errorReport())>0: print('ERRORS:',r)
