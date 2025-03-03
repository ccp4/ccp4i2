from __future__ import print_function

from report.CCP4ReportParser import *
import sys
from lxml import etree

class refmac_martin_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'refmac_martin'
    RUNNING = True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        # 'nooutput' mode would be used by another report class that wanted
        # to use some method(s) from this class for its own report
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
                
        if jobStatus == 'Running':
            #print ' '
            #print 'In running report for refmac_martin generator'
            #print ' '
            
            if xmlnode.haspath(".//RefmacInProgress"):
                progressGraph = self.addFlotGraph(title="Running refmac", select=".//RefmacOptimiseWeight/RefmacInProgress")
                progressGraph.addData(title="Cycle",    select="Cycle/number")
                progressGraph.addData(title="R_Factor", select="Cycle/r_factor")
                progressGraph.addData(title="R_Free",  select="Cycle/r_free")
                progressGraph.addPlot(plot='''<plot>
                    <title>R Factors</title>
                    <plottype>xy</plottype>
                    <plotline xcol="1" ycol="2">
                    <colour>blue</colour>
                    </plotline>
                    <plotline xcol="1" ycol="3">
                    <colour>red</colour>
                    </plotline>
                    </plot>''')
    
        else:            
            #Raid the refmac report of finished jobs
            summaryFold = self.addSummary()

    def addSummary(self, xmlnode=None, parent=None):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        
        summaryFold = parent.addFold(label='Summary of refinement', initiallyOpen=True)
        self.addScrollableDownloadableTable1(parent=summaryFold)
        self.addProgressGraph(parent=summaryFold)
        self.addTables(parent=summaryFold)

    def addProgressGraph(self, parent=None,xmlnode=None):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        #I *do not know* why This is needed
        
        #Note that when I add the progressgraph, I have to ensure that the select is rooted in my own xmlnode
        progressGraph = parent.addFlotGraph( title="Refinement results", xmlnode=self.xmlnode, select = ".//Overall_stats/stats_vs_cycle/new_cycle",style="height:250px; width:400px;float:left;")
        progressGraph.addData(title="Cycle",   select=".//cycle")
        progressGraph.addData(title="R-free",   select=".//r_free", expr="x if float(x)>=0.0 else ''")
        progressGraph.addData(title="R-factor", select=".//r_factor", expr="x if float(x)>=0.0 else ''")
        #For  lines that dont have a value for each point, the trick is to replace missing values with '-'.
        #Out of refmac, they are flagged with a value of -999.
        progressGraph.addData(title="rmsBONDx100", select=".//rmsBOND", expr="str(100.*float(x)) if float(x)>=0.0 else '-'")
        progressGraph.addData(title="rmsANGLE", select=".//rmsANGLE", expr="x if float(x)>=0.0 else '-'")
        
        plot = progressGraph.addPlotObject()
        plot.append('title','Refmac R-factors')
        plot.append('plottype','xy')
        plot.append('yrange', rightaxis='false')
        for coordinate, colour in [(2,'blue'),(3,'green')]:
            plotLine = plot.append('plotline',xcol=1,ycol=coordinate,rightaxis='false',colour=colour)
        plot.append('yrange', rightaxis='true')
        for coordinate, colour in [(4,'red'),(5,'purple')]:
            plotLine = plot.append('plotline',xcol=1,ycol=coordinate,rightaxis='true',colour=colour)

    def addTables(self, parent=None, xmlnode=None):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        #I *do not know* why This is needed
        clearingDiv = parent.addDiv(style="clear:both;")
        
        perCycleFold = parent.addFold(label='Per cycle statistics', initiallyOpen=False)
        tlsRefinementDiv = perCycleFold.addDiv(style="float:left; width:250px;")
        tlsRefinementDiv.addText(text="TLS refinement")
        table = tlsRefinementDiv.addTable(select=".//Overall_stats/stats_vs_cycle/new_cycle[rmsBOND=-999.0000][1] | //Overall_stats/stats_vs_cycle/new_cycle[rmsBOND=-999.0000][last()] ",style="width:250px;")
        table.addData(title="Cycle", select="cycle")
        table.addData(title="R-factor", select="r_factor")
        table.addData(title="R-free",   select="r_free",   expr="x if float(x)>=0.0 else '-'")
        
        refineAllDiv = perCycleFold.addDiv(style="float:left;")
        refineAllDiv.addText(text="XYZ, (Occ) and B refinement")
        refineAllTable = refineAllDiv.addTable(select="//Overall_stats/stats_vs_cycle/new_cycle[rmsBOND>-999.0000][1] | //Overall_stats/stats_vs_cycle/new_cycle[rmsBOND>-999.0000][last()] ")
        refineAllTable.addData(title="Cycle", select="cycle")
        refineAllTable.addData(title="R-factor", select="r_factor")
        refineAllTable.addData(title="R-free",   select="r_free",   expr="x if float(x)>=0.0 else '-'")
        refineAllTable.addData(title="RMS Deviations", subtitle="Bond", select="rmsBOND", expr="x if float(x)>=0.0 else '-'")
        refineAllTable.addData(subtitle="Angle deviations", select="rmsANGLE", expr="x if float(x)>=0.0 else '-'")
        
        clearingDiv = perCycleFold.addDiv(style="clear:both;")
        detailFold = perCycleFold.addFold(label='All cycles',initiallyOpen=False,style="clear:both;")
        fullTable = detailFold.addTable(select="//Overall_stats/stats_vs_cycle/new_cycle")
        for title,select in  [[ "Cycle" ,"cycle" ], [ "R-factor"  , "r_factor" ],[ "R-free" , "r_free" ],["RMS Bond","rmsBOND"],["RMS Angle","rmsANGLE"]]:
            fullTable.addData(title=title,select=select, expr="x if float(x)>=0.0 else '-'")

    def addSymmetryAnalysis(self, parent=None):
        if parent is None:
            parent=self
        equivalenceNodes = self.xmlnode.xpath('//NCS/equivalence')
        if len(equivalenceNodes) == 0: return
        symmetryFold = parent.addFold(label='Non-crystallographic symmetry identified by Refmac')
        table = symmetryFold.addTable(select='//NCS/equivalence')
            
        table.addData(title='',data=['Selection 1','Selection 2','Residues aligned','Score','RMS','Average(RMS loc)'])
        for equivalenceNode in equivalenceNodes:
            table.addData(title='',data=[equivalenceNode.get(property) for property in ['selection1','selection2','nEquivalent','score','rms','ave_rmsLoc']])
        return

    def addTwinningAnalysis(self, parent=None):
        if parent is None:
            parent=self
        twinningNode = self.xmlnode.xpath0('//Twinning')
        if twinningNode is None: return
        twinningFold = parent.addFold(label='Refmac analysis of twinning')
        criteriaText = ""
        rMergeNode = self.xmlnode.xpath0('//Twinning/RmergeFilter')
        if rMergeNode is not None: criteriaText += ('Twins included if Rmerge < %s.'%rMergeNode.text)
        fractionNode = self.xmlnode.xpath0('//Twinning/FractionFilter')
        if fractionNode is not None: criteriaText += ('Twins included if Fraction > %s.'%fractionNode.text)
        if len(criteriaText)>0:
            textDiv = twinningFold.addDiv()
            newText = textDiv.addText(text=criteriaText)
        table = twinningFold.addTable(select='//Twinning/TwinningSymmetryOperator',style="width:250px;float:left;")
        for title,select in  [[ "Symmetry operator" ,"SymmetryOperator" ], [ "Rmerge"  , "Rmerge" ]]:
            table.addData(title=title,select=select)
        table = twinningFold.addTable(select='//Twinning/TwinOperator',style="width:250px;float:left;")
        for title,select in  [[ "Twinning operator" ,"TwinOperator" ], [ "Fraction"  , "Fraction" ]]:
            table.addData(title=title,select=select)
        clearingDiv = parent.addDiv(style="clear:both;")
        return

    def addSmartieGraphs(self, parent=None):
        if parent is None:
            parent=self
        reportFold = parent.addFold(label='Other plots from log file')
        graphTableList = self.xmlnode.xpath('// CCP4ApplicationOutput/CCP4Table')
        graphgroup = reportFold.addFlotGraphGroup(style="width:500px;  height:400px;")
        for graphTableNode in graphTableList:
            graph = graphgroup.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title") )
            graph = graph.addPimpleData(xmlnode=graphTableNode)

    def addOutlierAnalysis(self, parent=None):
        if parent is None:
            parent=self
        
        outlierFold = parent.addFold(label="Outliers identified by Refmac")
        #Check to see if outliers are captured in the program XML
        outliersByCriteriaNodes = self.xmlnode.xpath('//OutliersByCriteria')
        if len(outliersByCriteriaNodes) == 0 or len(outliersByCriteriaNodes[0])==0:
            outlierFold.append('<span style="font-size:110%">No outliers observed </span>')
        else:
            outlierFold.append('<span style="font-size:110%">Residues failing one or more outlier test are flagged below with corresponding Z-score </span>')
            #identify a unique list of amino acids flagged by Refmac, and their associated set of deviations, flagged by the most deviant score observed
            naughtyBits = {}
            from sets import Set
            criteriaThatFail = Set()
            for outlierDictNode in outliersByCriteriaNodes[0]:
                outliers = outlierDictNode.xpath('Outlier')
                criteriaThatFail.add(outlierDictNode.tag)
                for outlier in outliers:
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
            criteriaFold = outlierFold.addFold(label='Criteria used to spot outliers')
            criteriaTable = criteriaFold.addTable(title='Ooh',label='Aah',select='//OutliersByCriteria')
            criteriaTable.addData(title='Interaction type',data=[outliersByCriteriaNodes[0][i].tag for i in range(len(outliersByCriteriaNodes[0]))])
            criteriaTable.addData(title='Criteria',select='//Criteria')

    def appendMolDisp(self, molDataNode=None, selectionText='all', carbonColour='yellow', othersByElement=True,style='CYLINDER',bondOrder=False):
        if molDataNode is None: return
        molDispNode = etree.SubElement(molDataNode,'MolDisp')
        if bondOrder:
            drawingStyleNode = etree.fromstring('''<show_multiple_bonds>1</show_multiple_bonds>
                                                <deloc_ring>1</deloc_ring>
                                                </drawing_style>''')
            molDispNode.append(drawingStyleNode)
        
        selectNode = etree.SubElement(molDispNode,'select')
        selectNode.text = selectionText
        '''
            selectionParametersNode = etree.SubElement(molDispNode,'selection_parameters')
            selectNode = etree.SubElement(selectionParametersNode,'select')
            selectNode.text = 'cid'
            
            cidNode = etree.SubElement(selectionParametersNode,'cid')
            cidNode.text = selectionText
            '''
        colourParametersNode = etree.SubElement(molDispNode,'colour_parameters')
        
        colourModeNode = etree.SubElement(colourParametersNode,'colour_mode')
        colourModeNode.text = 'one_colour'
        if othersByElement:
            nonCNode = etree.SubElement(colourParametersNode,'non_C_atomtype')
            nonCNode.text='1'
        oneColourNode = etree.SubElement(colourParametersNode,'one_colour')
        oneColourNode.text=carbonColour
        #
        styleNode = etree.SubElement(molDispNode,'style')
        styleNode.text=style

    def addMap(self, mapDataNode = None, fPhiObjectName='FPHIOUT', fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=1.0, isDifmap=False):
        if mapDataNode is None: return
        filenameNode = etree.SubElement(mapDataNode,'filename',database='filenames/'+fPhiObjectName)
        columnsNode = None
        if isDifmap:
            columnsNode = etree.SubElement(mapDataNode,'difColumns')
            isDifferenceMapNode = etree.SubElement(mapDataNode,'isDifferenceMap')
            isDifferenceMapNode.text='True'
        else:
            columnsNode = etree.SubElement(mapDataNode,'columns')
            isDifferenceMapNode = etree.SubElement(mapDataNode,'isDifferenceMap')
            isDifferenceMapNode.text='False'
        fColNode = etree.SubElement(columnsNode,'F')
        fColNode.text = fCol
        phiColNode = etree.SubElement(columnsNode,'PHI')
        phiColNode.text = phiCol
        
        modelNode = etree.SubElement(mapDataNode,'model')
        modelNode.text = model
        gridSizeNode = etree.SubElement(mapDataNode,'gridSize')
        gridSizeNode.text=str(gridSize)
        contourUnitsNode = etree.SubElement(mapDataNode,'contourUnits')
        contourUnitsNode.text = contourUnits
        mapDispNode = etree.SubElement(mapDataNode,'MapDisp')
        contourLevelNode = etree.SubElement(mapDispNode,'contourLevel')
        contourLevelNode.text = str(contourLevel)
        differenceNode = etree.SubElement(mapDispNode,'difference')
        colourNode = etree.SubElement(mapDispNode,'colour')
        if isDifmap:
            differenceNode.text='1'
            colourNode.text='green'
            colourNode2 = etree.SubElement(mapDispNode,'second_colour')
            colourNode2.text='red'
        else:
            differenceNode.text='0'
            colourNode.text='blue'

    def addRefinementPictures(self, xmlnode=None, jobInfo=None, parent=None):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        
        from lxml import etree
        #I *do not know* why This is needed
        clearingDiv = parent.addDiv(style="clear:both;")
        
        pictureFold = parent.addFold(label='Picture')
        pictureGallery = pictureFold.addObjectGallery(style='float:left;',height='450px', tableWidth='260px', contentWidth='450px')
        clearingDiv = parent.addDiv(style="clear:both;")
        jobDirectory = jobInfo['fileroot']
        from core import CCP4Utils
        ccp4i2_root = CCP4Utils.getCCP4I2Dir()
        import os
        baseScenePath = os.path.join(ccp4i2_root,'pipelines','prosmart_refmac','script','prosmart_refmac_1.scene.xml')
        monomerNodes = xmlnode.xpath('//RefmacWeight[1]/REFMAC/ModelComposition/Monomer')
        
        #Subsets for snapshotting are each observed monomer and the whole molecule
        interestingBits = [monomerNode.get('id') for monomerNode in monomerNodes] + ['all']
        
        iMonomer = 1
        for interestingBit in interestingBits:
            baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath)
            sceneNode = baseSceneXML.xpath('//scene')[0]
            
            #Define data and associated display objects
            dataNode = etree.SubElement(sceneNode,'data')
            
            #Load model and draw representations
            molDataNode = etree.SubElement(dataNode,'MolData',id='id1')
            fileNode = etree.SubElement(molDataNode,'filename',database='filenames/XYZOUT')
            
            #Define view (possibly some advantage in doing this after loading the moelcule from which it is inferred)
            viewNode = etree.SubElement(sceneNode,'View')
            autoScaleNode = etree.SubElement(viewNode,'scale_auto')
            autoScaleNode.text = '1'
            autoScaleNode = etree.SubElement(viewNode,'slab_enabled')
            autoScaleNode.text = '1'
            centreMolDataNode = etree.SubElement(viewNode,'centre_MolData')
            centreMolDataNode.text = 'id1'
            centreSelectionNode = etree.SubElement(viewNode,'centre_selection')
            centreSelectionNode.text = interestingBit
            
            orientationAutoNode = etree.SubElement(viewNode,'orientation_auto')
            orientationAutoMolDataNode = etree.SubElement(orientationAutoNode,'molData')
            orientationAutoMolDataNode.text = 'id1'
            orientationAutoSelectionNode = etree.SubElement(orientationAutoNode,'selection')
            orientationAutoSelectionNode.text = interestingBit
            
            #Draw stuff around the interesting bit
            selectionText='neighb cid="'+interestingBit+'" maxd=10.0 group=all excl=central'
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='green', othersByElement=True,style='CYLINDERS')
            
            #Draw the interesting bit
            selectionText=interestingBit
            self.appendMolDisp(molDataNode=molDataNode, selectionText=selectionText, carbonColour='yellow', othersByElement=True,style='BALLSTICK')
            
            #Load and draw 2Fo-Fc
            mapDataNode = etree.SubElement(dataNode,'MapData',id='id3')
            self.addMap(mapDataNode, fPhiObjectName='FPHIOUT', fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=1.0, isDifmap=False)
            
            #Load and draw Fo-Fc
            mapDataNode = etree.SubElement(dataNode,'MapData',id='id4')
            self.addMap(mapDataNode, fPhiObjectName='DIFFPHIOUT', fCol='F', phiCol='PHI', gridSize=0.5, contourUnits='sigma', model='id1', contourLevel=3.0, isDifmap=True)
            
            # Dump out the XML
            et = etree.ElementTree(baseSceneXML)
            sceneFilePath = os.path.join(jobDirectory,'monomer'+str(iMonomer)+'.scene.xml')
            et.write(sceneFilePath,pretty_print=True)
            # And add the picture
            pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of selection "'+interestingBit+'"')

            iMonomer+=1

    def addScrollableDownloadableTable1(self, xmlnode=None, parent=None):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        #create a "shell" div to contain the scrollable table and the hyperlink
        scrollableDownloadableTableDiv = parent.addDiv(style="height:250px; width:260px;float:left;outline:black solid thin;margin-top:2px;")
        #place a scrollable div into the shell: the table will be inserted into this div
        scrollableTableDiv = scrollableDownloadableTableDiv.addDiv(style="height:225px; width:260px;clear:both;overflow:auto;")
        #Put table1 into this (autoscrolling) div
        table1 = self.addTable1(xmlnode=xmlnode, parent=scrollableTableDiv)
        #Add a hyperlink to this table in the "tables_as_csv_files" directory
        #print '\n\nJobInfo',self.jobInfo
        download = scrollableDownloadableTableDiv.addDownload(jobInfo=self.jobInfo,dataName=table1.id)

    def addTable1(self, xmlnode=None, parent=None, downloadable=False):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        table1 = parent.addTable(xmlnode=xmlnode, style="width:240px;float:left;", downloadable=downloadable)
        
        ResolutionLowNode =xmlnode.xpath('.//Overall_stats/resolution_low')
        ResolutionHighNode =xmlnode.xpath('.//Overall_stats/resolution_high')
        ReflectionsAll =xmlnode.xpath('.//Overall_stats/n_reflections_all')
        ResolutionsFree =xmlnode.xpath('.//Overall_stats/n_reflections_free')
        ReflectionsWorkNode =xmlnode.xpath('.//Overall_stats/resolution_high')
        RFactorNodes = xmlnode.xpath('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_factor')
        RFreeNodes = xmlnode.xpath('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_free')
        RMSBondsNodes = xmlnode.xpath('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/rmsBOND')
        RMSAnglesNodes = xmlnode.xpath('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/rmsANGLE')
        
        MeanBChainNameNodes = xmlnode.xpath('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/chain_name')
        MeanBAllCountNodes = xmlnode.xpath('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/number')
        MeanBAllAverageNodes = xmlnode.xpath('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/average')
        
        
        table1.addData(title='Statistic',data=['Resolution','No. reflections all/free','R-factor/R-free','RMS Deviations','Bonds','Angles','Chain mean B (No. atoms)']+[MeanBChainNameNodes[i].text for i in range(len(MeanBChainNameNodes))])
        table1.addData(title='Value',data=[ResolutionLowNode[0].text+'-'+ResolutionHighNode[0].text,
                                           ReflectionsAll[0].text+'/'+ResolutionsFree[0].text,
                                           RFactorNodes[0].text+'/'+RFreeNodes[0].text,
                                           ' ',
                                           RMSBondsNodes[0].text,
                                           RMSAnglesNodes[0].text,
                                           ' ',
                                           ] + [(MeanBAllAverageNodes[i].text+'('+MeanBAllCountNodes[i].text+')') for i in range(len(MeanBChainNameNodes))])
        return table1

def test(xmlFile=None,jobId=None,reportFile=None):
    import sys,os
    print(xmlFile)
    try:
        text = open( xmlFile ).read()
        xmlnode = etree.fromstring( text, PARSER() )
    except:
        print('FAILED loading XML file:', kw['xmlFile'])
    if reportFile is None and xmlFile is not None:
        reportFile = os.path.join(os.path.split(xmlFile)[0],'report.html')
    r = refmac_martin_report(xmlFile=xmlFile,jobId=jobId, xmlnode=xmlnode)
    r.as_html_file(reportFile)

if __name__ == "__main__":
    import sys
    refmac_martin_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
