from __future__ import print_function

from report.CCP4ReportParser import *
import sys
import xml.etree.ElementTree as etree

class refmac_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'refmac'
    TASKTITLE = 'REFMAC5 - Macromolecular refinement'
    RUNNING = True
    SEPARATEDATA=True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        # 'nooutput' mode would be used by another report class that wanted
        # to use some method(s) from this class for its own report
        self.outputXml = jobStatus is not None and jobStatus.lower().count('running')
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        
        self.addDiv(style='clear:both;')

        if jobStatus.lower().count('running'):
            self.addRunningProgressGraph(self)
        else:
            #Raid the refmac report of finished jobs
            summaryFold = self.addSummary()

    def addSummary(self, xmlnode=None, parent=None, withTables=True):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        
        summaryFold = parent.addFold(label='Summary of refinement', brief='Summary', initiallyOpen=True)
        self.addScrollableDownloadableTable1(parent=summaryFold)
        self.addProgressGraph(parent=summaryFold)
        if withTables: self.addTables(parent=summaryFold)
    
    def addRunningProgressGraph(self, parent):
        if len(self.xmlnode.findall("Cycle"))>0:
            progressGraph = parent.addFlotGraph(title="Running refmac",select="Cycle",style="height:250px; width:400px;float:left;",outputXml=self.outputXml,internalId="SummaryGraph")
            progressGraph.addData(title="Cycle",    select="number")
            progressGraph.addData(title="R_Factor", select="r_factor")
            progressGraph.addData(title="R_Free",  select="r_free")
            plot = progressGraph.addPlotObject()
            plot.append('title','Running refmac R-factors')
            plot.append('plottype','xy')
            plot.append('yrange', rightaxis='false')
            plot.append( 'xlabel', 'Cycle' )
            plot.append( 'xintegral', 'true' )
            plot.append( 'ylabel', 'R-factor' )
            plot.append( 'rylabel', 'Geometry' )
            for coordinate, colour in [(2,'blue'),(3,'green')]:
                plotLine = plot.append('plotline',xcol=1,ycol=coordinate,rightaxis='false',colour=colour)
            
            rmsBonds = self.xmlnode.findall('Cycle/rmsBonds')
            if len(rmsBonds)> 0:
                plot.append('yrange', rightaxis='true')
                cycleNodes = self.xmlnode.findall('Cycle')
                data = []
                for cycleNode in cycleNodes:
                    try: data.append(cycleNode.findall('rmsBonds')[0].text)
                    except: data.append(None)
                progressGraph.addData(title="rmsBonds",  data=data)
                plotLine = plot.append('plotline',xcol=1,ycol=4,rightaxis='true',colour='red')


    def addProgressGraph(self, parent=None,xmlnode=None):
        if parent is None: parent=self
        if xmlnode is None: xmlnode = self.xmlnode
        #I *do not know* why This is needed
        
        #Note that when I add the progressgraph, I have to ensure that the select is rooted in my own xmlnode
        progressGraph = parent.addFlotGraph( title="Refinement results", xmlnode=self.xmlnode, select = ".//Overall_stats/stats_vs_cycle/new_cycle",style="height:250px; width:400px;float:left;")
        progressGraph.addData(title="Cycle",   select=".//cycle")
        progressGraph.addData(title="R-free",   select=".//r_free", expr="x if float(x)>=0.0 else '-'")
        progressGraph.addData(title="R-factor", select=".//r_factor", expr="x if float(x)>=0.0 else ''")
        #For  lines that dont have a value for each point, the trick is to replace missing values with '-'.
        #Out of refmac, they are flagged with a value of -999.
        progressGraph.addData(title="rmsBONDx100", select=".//rmsBOND", expr="str(100.*float(x)) if float(x)>=0.0 else '-'")
        progressGraph.addData(title="rmsANGLE", select=".//rmsANGLE", expr="x if float(x)>=0.0 else '-'")
        
        plot = progressGraph.addPlotObject()
        plot.append('title','Refmac R-factors')
        plot.append('plottype','xy')
        plot.append('yrange', rightaxis='false')
        plot.append( 'xlabel', 'Cycle' )
        plot.append( 'ylabel', 'R-factor' )
        plot.append( 'rylabel', 'Geometry' )
        plot.append('xintegral','true')
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
        perCycleFold = parent.addFold(label='Per cycle statistics', brief='Per cycle', initiallyOpen=False)
        
        all_cycles = xmlnode.findall('.//Overall_stats/stats_vs_cycle/new_cycle')
        ncyc = len(all_cycles)
        cycle_data = {'mode':['-']*ncyc,'cycle':['-']*ncyc,'r_factor':['-']*ncyc,'r_free':['-']*ncyc,'fscAver':['-']*ncyc,'fscAverFree':['-']*ncyc,'rmsBOND':['-']*ncyc,'rmsANGLE':['-']*ncyc,'rmsCHIRAL':['-']*ncyc}
        
        idx = 0
        for cycle in all_cycles:
           try:
              if len(xmlnode.findall("RigidMode"))>0:
                 cycle_data['mode'][idx] = 'Rigid'
              elif float(cycle.findall('rmsBOND')[0].text) >= 0.0:
                 cycle_data['mode'][idx] = 'Restr'
              elif str(float(cycle.findall('rmsBOND')[0].text)) == '-999.0' and len(xmlnode.findall("TLSMode"))>0:
                 cycle_data['mode'][idx] = 'TLS'
              else: raise
           except:
              print("*** ERROR - UNABLE TO INTERPRET REFMAC5 LOG FILE TO CONSTRUCT REFINEMENT TABLES - PLEASE CONTACT CCP4 WITH THIS ERROR (ccp4@ccp4.ac.uk) ***")
              return None

           try: cycle_data['cycle'][idx] = str(int(cycle.findall('cycle')[0].text))
           except: pass
           try: cycle_data['r_factor'][idx] = str(float(cycle.findall('r_factor')[0].text))
           except: pass
           try: cycle_data['r_free'][idx] = str(float(cycle.findall('r_free')[0].text))
           except: pass
           try: cycle_data['fscAver'][idx] = str(float(cycle.findall('fscAver')[0].text))
           except: pass
           try: cycle_data['fscAverFree'][idx] = str(float(cycle.findall('fscAverFree')[0].text))
           except: pass
           if cycle_data['mode'][idx] == 'Restr':
              try: cycle_data['rmsBOND'][idx] = str(float(cycle.findall('rmsBOND')[0].text))
              except: pass
              try: cycle_data['rmsANGLE'][idx] = str(float(cycle.findall('rmsANGLE')[0].text))
              except: pass
              try: cycle_data['rmsCHIRAL'][idx] = str(float(cycle.findall('rmsCHIRAL')[0].text))
              except: pass
           idx = idx+1
        
        #print(cycle_data)
        
        RigidIdx = [idx for idx, val in enumerate([mode == 'Rigid' for mode in cycle_data['mode']]) if val]
        TLSIdx = [idx for idx, val in enumerate([mode == 'TLS' for mode in cycle_data['mode']]) if val]
        RestrIdx = [idx for idx, val in enumerate([mode == 'Restr' for mode in cycle_data['mode']]) if val]
        
        # Only display first and last in the summary tables
        if len(RigidIdx)>2: RigidIdx = [RigidIdx[0],RigidIdx[-1]]
        if len(TLSIdx)>2: TLSIdx = [TLSIdx[0],TLSIdx[-1]]
        if len(RestrIdx)>2: RestrIdx = [RestrIdx[0],RestrIdx[-1]]
        start_end = ['Start','End']
        
        rigid_data = {}
        tls_data = {}
        restr_data = {}
        for key,val in cycle_data.items():
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

        clearingDiv = perCycleFold.addDiv(style="clear:both;")
        detailFold = perCycleFold.addFold(label='All cycles',initiallyOpen=False,brief='All cycles')
        fullTable = detailFold.addTable()
        if len(TLSIdx)>0 and len(RestrIdx)>0:
           mode = ['TLS']*len(TLSIdx) + ['Full Atom']*len(RestrIdx)
        fullTable.addData(title="Cycle", data=cycle_data['cycle'])
        fullTable.addData(title="R-factor", data=cycle_data['r_factor'])
        fullTable.addData(title="R-free", data=cycle_data['r_free'])
        #fullTable.addData(title="FSCAv", data=cycle_data['fscAver'])
        #fullTable.addData(title="FSCAv-free", data=cycle_data['fscAverFree'])
        if len(RestrIdx)>0:
           fullTable.addData(title="RMSD (bond/angle/chiral)", subtitle="Bond", data=cycle_data['rmsBOND'])
           fullTable.addData(subtitle="Angle", data=cycle_data['rmsANGLE'])
           fullTable.addData(subtitle="Chiral", data=cycle_data['rmsCHIRAL'])
        
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

    def addSmartieGraphs(self, parent=None,internalIdPrefix=''):
        if parent is None:
            parent=self
        reportFold = parent.addFold(label='Other plots from log file',brief='Other')
        
        """
        #MN Here explicitly correct the absence of annotation of xscale as oneoversqrt in plots versus resln
        for fixmeNode in self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table[contains(@title,'resln')]/plot"):
            fixmeNode.append(etree.fromstring('<xscale>oneoversqrt</xscale>'))
        """
        #SJM - A python xml library version of above.    
        for tableNode in self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table"):
            if "title" in tableNode.attrib and "resln" in tableNode.attrib["title"]:
                plots = tableNode.findall("./plot")
                for plot in plots:
                    plot.append(etree.fromstring('<xscale>oneoversqrt</xscale>'))
        
        graphTableList = self.xmlnode.findall('SmartieGraphs/CCP4ApplicationOutput/CCP4Table')
        gallery = reportFold.addObjectGallery(height='300px',contentWidth='420px',tableWidth='360px',style='float:left;width:800px;')
        
        isFirstGraph = True
        for iGraph, graphTableNode in enumerate(graphTableList):
            graph = gallery.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title"), internalId=internalIdPrefix+'SmartiePlot'+str(iGraph), outputXml=self.outputXml, label=graphTableNode.get("title"),style="width:410px;height:290px;", initiallyDrawn=isFirstGraph )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
            isFirstGraph = False

    def addOutlierAnalysis(self, parent=None):
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
        molDispNode = etree.SubElement(molDataNode,'MolDisp')
        if bondOrder:
            drawingStyleNode = etree.fromstring('<drawing_style><show_multiple_bonds>1</show_multiple_bonds><deloc_ring>1</deloc_ring></drawing_style>')
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
            colourNode.text='ice blue'

    def addRefinementPictures(self, xmlnode=None, jobInfo=None, parent=None, objectNameMap={}, initiallyOpen=False):
        if parent is None: parent = self
        if xmlnode is None: xmlnode = self.xmlnode
        if jobInfo is None: jobInfo = self.jobInfo
        
        #Provide default filenames to identify database object
        coordObjectName = objectNameMap.get('XYZ','XYZOUT')
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

        # Ligands are now pre-filtered by CPdbDataComposition:
        # - Only NonPolymer entities (not polymer chains)
        # - Not water molecules
        # - Not single metal ions
        # - At least 5 atoms (significant ligands)
        # Format: "chain:resname:seqnum" e.g. "A:ATP:501"
        interestingBits = [monomerNode.get('id') for monomerNode in monomerNodes] + ['all']
        
        for iMonomer, interestingBit in enumerate(interestingBits):
            baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath,useLXML=False) #This is lxml, not xml ...
            sceneNode = baseSceneXML#.findall('.//scene')[0]
            
            #Define data and associated display objects
            dataNode = etree.SubElement(sceneNode,'data')
            
            #Load model and draw representations
            molDataNode = etree.SubElement(dataNode,'MolData',id='id1')
            fileNode = etree.SubElement(molDataNode,'filename',database='filenames/'+coordObjectName)
            
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
                                customResNode = etree.fromstring('''<customResCIFFiles>
                                    <cifmonomer>
                                    <name>'''+tlc+'''</name>
                                    <filename>'''+dictPath+'''</filename>
                                    </cifmonomer>
                                    </customResCIFFiles>''')
                                molDataNode.append(customResNode)
            except BaseException as err:
                print('Failed analysing DICT'+"Base error: {0}".format(err))
            
            #Define view (possibly some advantage in doing this after loading the moelcule from which it is inferred)
            viewNode = etree.SubElement(sceneNode,'View')
            autoScaleNode = self.subElementWithText(viewNode,'scale_auto',1)
            autoScaleNode = self.subElementWithText(viewNode,'slab_enabled',1)
            centreMolDataNode = self.subElementWithText(viewNode,'centre_MolData','id1')
            centreSelectionNode = self.subElementWithText(viewNode,'centre_selection',interestingBit)
            
            orientationAutoNode = etree.SubElement(viewNode,'orientation_auto')
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

            # Dump out the XML with deterministic filename based on ligand ID
            et = etree.ElementTree(baseSceneXML)
            # Sanitize ligand ID for filename (replace : with _ to avoid filesystem issues)
            safe_ligand_id = interestingBit.replace(':', '_')
            sceneFilePath = os.path.join(jobDirectory, f'ligand_{safe_ligand_id}.scene.xml')
            et.write(sceneFilePath)
            # And add the picture
            pic = pictureGallery.addPicture(sceneFile=sceneFilePath, label=f'Ligand {interestingBit}')

        #And finally the full monty picture :-)
        baseSceneXML = CCP4Utils.openFileToEtree(baseScenePath,useLXML=False)
        sceneNode = baseSceneXML#.findall('scene')[0]
        
        #Define data and associated display objects
        dataNode = etree.SubElement(sceneNode,'data')
        
        #Load model and draw representations
        molDataNode = etree.SubElement(dataNode,'MolData',id='id1')
        fileNode = etree.SubElement(molDataNode,'filename',database='filenames/'+coordObjectName)
        
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
                            customResNode = etree.fromstring('''<customResCIFFiles>
                                <cifmonomer>
                                <name>'''+tlc+'''</name>
                                <filename>'''+dictPath+'''</filename>
                                </cifmonomer>
                                </customResCIFFiles>''')
                            molDataNode.append(customResNode)
        except BaseException as err:
            print('Failed analysing DICT'+"Base error: {0}".format(err))
        
        #Define view (possibly some advantage in doing this after loading the moelcule from which it is inferred)
        viewNode = etree.SubElement(sceneNode,'View')
        autoScaleNode = self.subElementWithText(viewNode,'scale_auto',1)
        autoScaleNode = self.subElementWithText(viewNode,'slab_enabled',1)
        centreMolDataNode = self.subElementWithText(viewNode,'centre_MolData','id1')
        centreSelectionNode = self.subElementWithText(viewNode,'centre_selection',interestingBit)
        
        orientationAutoNode = etree.SubElement(viewNode,'orientation_auto')
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

        # Dump out the XML
        et = etree.ElementTree(baseSceneXML)
        sceneFilePath = os.path.join(jobDirectory,'fullMonty.scene.xml')
        et.write(sceneFilePath)
        # And add the picture
        pic = pictureGallery.addPicture(sceneFile=sceneFilePath,label='Picture of Full Monty')


    def subElementWithText(self, parent, tag, text):
        result = etree.SubElement(parent, tag)
        result.text = str(text)
        return result

    def addScrollableDownloadableTable1(self, xmlnode=None, parent=None,internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        #create a "shell" div to contain the scrollable table and the hyperlink
        scrollableDownloadableTableDiv = parent.addDiv(style="height:250px; width:315px;float:left;margin-top:2px;")
        #place a scrollable div into the shell: the table will be inserted into this div
        scrollableTableDiv = scrollableDownloadableTableDiv.addDiv(style="height:225px; width:300px;clear:both;overflow:auto;")
        #Put table1 into this (autoscrolling) div
        table1 = self.addTable1(xmlnode=xmlnode, parent=scrollableTableDiv,internalId=internalId)
        #scrollableDownloadableTableDiv.addDiv(style="height:10px;width:15px; float:right;")
        #Add a hyperlink to this table in the "tables_as_csv_files" directory
        #print '\n\nJobInfo',self.jobInfo
        download = scrollableDownloadableTableDiv.addDownload(jobInfo=self.jobInfo,dataName=table1.id)

    def addTable1(self, xmlnode=None, parent=None, downloadable=False,internalId='Table1'):
        if xmlnode is None: xmlnode = self.xmlnode
        if parent is None: parent = self
        
        table1 = parent.addTable(xmlnode=xmlnode, style="width:240px;float:left;", downloadable=downloadable,outputXml=self.outputXml,internalId=internalId)
        
        ResolutionLowNode =xmlnode.findall('.//Overall_stats/resolution_low')
        ResolutionHighNode =xmlnode.findall('.//Overall_stats/resolution_high')
        ReflectionsAll =xmlnode.findall('.//Overall_stats/n_reflections_all')
        ResolutionsFree =xmlnode.findall('.//Overall_stats/n_reflections_free')
        ReflectionsWorkNode =xmlnode.findall('.//Overall_stats/resolution_high')
        RFactorNodes = xmlnode.findall('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_factor')
        RFreeNodes = xmlnode.findall('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/r_free')
        RMSBondsNodes = xmlnode.findall('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/rmsBOND')
        RMSAnglesNodes = xmlnode.findall('.//Overall_stats/stats_vs_cycle/new_cycle[last()]/rmsANGLE')
        
        MeanBChainNameNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/chain_name')
        MeanBAllCountNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/number')
        MeanBAllAverageNodes = xmlnode.findall('.//Overall_stats/bvalue_stats/chain_by_chain/new_chain/all/average')
        
        statisticNames = []
        statisticValues = []
        
        statisticNames.append('Resolution')
        if len(ResolutionLowNode)>0: statisticValues.append("{0:.2f}".format(float(ResolutionLowNode[0].text))+'-'+"{0:.2f}".format(float(ResolutionHighNode[0].text)))
        else: statisticValues.append("-")
        
        statisticNames.append('No. reflections all/free')
        if len(ReflectionsAll)>0:
            nRefText = ReflectionsAll[0].text+'/'
            if len(ResolutionsFree)>0: nRefText += ResolutionsFree[0].text
            else: nRefText += "-"
            statisticValues.append(nRefText)
        else: statisticValues.append("-")
        
        statisticNames.append('R-factor/R-free')
        if len(RFactorNodes)>0:
            text = RFactorNodes[0].text+'/'
            if len(RFreeNodes)>0: text += RFreeNodes[0].text
            else: text += "-"
            statisticValues.append(text)
        else: statisticValues.append("-")

        if len(self.xmlnode.findall("RigidMode"))==0:
           statisticNames.append('<i>RMS Deviations</i>')
           statisticValues.append(' ')

           statisticNames.append('Bonds')
           if len(RMSBondsNodes)>0:
               statisticValues.append(RMSBondsNodes[0].text)
           else: statisticValues.append("-")
           
           statisticNames.append('Angles')
           if len(RMSAnglesNodes)>0:
               statisticValues.append(RMSAnglesNodes[0].text)
           else: statisticValues.append("-")
           
           statisticNames.append('<i>Chain B-factors</i>')
           statisticValues.append('<i>mean B (#atoms)</i>')
           
           statisticNames += [MeanBChainNameNodes[i].text for i in range(len(MeanBChainNameNodes))]
           statisticValues +=  [("{0:.1f}".format(float(MeanBAllAverageNodes[i].text))+'('+MeanBAllCountNodes[i].text+')') for i in range(len(MeanBChainNameNodes))]
        
        table1.addData(title='Statistic',data=statisticNames)
        table1.addData(title='Value',data=statisticValues)
        
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
    r = refmac_report(xmlFile=xmlFile,jobId=jobId, xmlnode=xmlnode)
    r.as_html_file(reportFile)

if __name__ == "__main__":
    import sys
    refmac_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
