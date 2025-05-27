import xml.etree.ElementTree as etree

from ....report.CCP4ReportParser import Report
from ...import_mosflm.script import import_mosflm_report


class mosflm_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'mosflm'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        if jobStatus.lower()=="nooutput": return
        
        integrationResultNode=self.makeIntegrationResult()
        #print etree.tostring(integrationResultNode)
        importMosflmReport = import_mosflm_report.import_mosflm_report(xmlnode=integrationResultNode,jobInfo=jobInfo, jobStatus="nooutput")

        detectorFold = self.addFold(label="Detector and beam", initiallyOpen=True)
        detectorDiv = detectorFold.addDiv(style="border:2px solid black;display:inline-block;")
        graphLineChooser = detectorDiv.addGraphLineChooser(tableWidth='220px',contentWidth='260px',height='270px')
        importMosflmReport.beamGraph(parent=graphLineChooser)
        importMosflmReport.centralSpotProfile(parent=detectorDiv)

        clearingDiv = self.addDiv(style="clear:both;")
        
        crystalFold = self.addFold(label="Crystal", initiallyOpen=True)
        crystalDiv = crystalFold.addDiv(style="border:2px solid black;display:inline-block;")
        otherGraphLineChooser = crystalDiv.addGraphLineChooser(tableWidth='220px',contentWidth='260px',height='270px')
        importMosflmReport.crystalGraph(parent=otherGraphLineChooser)
        importMosflmReport.blockProfiles(parent=crystalDiv)
        
        clearingDiv = self.addDiv(style="clear:both;")

        dataFold = self.addFold(label="Data", initiallyOpen=True)
        dataDiv = dataFold.addDiv(style="border:2px solid black;display:inline-block;")
        thirdGraphLineChooser = dataDiv.addGraphLineChooser(tableWidth='220px',contentWidth='260px',height='270px')
        importMosflmReport.dataGraph(parent=thirdGraphLineChooser)
        importMosflmReport.histograms(parent=dataDiv)

    def makeIntegrationResult(self):
        collectedData = {}
        integrationNodes = self.xmlnode.findall('.//integration_positional_refinement')
        integration_positional_refinement_tags = []
        for integrationNode in integrationNodes:
            imageNameNodes = integrationNode.findall('image_name')
            imageName = imageNameNodes[0].text
            if not imageName in collectedData: collectedData[imageName] = {}
            for child in integrationNode:
                if child.tag != 'image_name':
                    if child.text is not None:
                        collectedData[imageName][child.tag] = child.text.strip()
                        if child.tag not in integration_positional_refinement_tags: integration_positional_refinement_tags.append (child.tag)
    
        integrationNodes = self.xmlnode.findall('.//integration_postrefinement')
        integration_post_refinement_tags = []
        for integrationNode in integrationNodes:
            imageNameNodes = integrationNode.findall('image_name')
            imageName = imageNameNodes[0].text
            if not imageName in collectedData: collectedData[imageName] = {}
            for child in integrationNode:
                if child.tag != 'image_name':
                    if child.text is not None:
                        collectedData[imageName][child.tag] = child.text.strip()
                        if child.tag not in integration_post_refinement_tags: integration_post_refinement_tags.append (child.tag)
    
        spotProfiles = {}
        spotProfileNodes = self.xmlnode.findall('.//spot_profile')
        spot_profile_tags = []
        for spotProfileNode in spotProfileNodes:
            imageNameNodes = spotProfileNode.findall('image_name')
            imageName = imageNameNodes[0].text
            if not imageName in spotProfiles: spotProfiles[imageName] = {}
            for child in spotProfileNode:
                if child.text is not None:
                    if child.tag == 'raw_data' or child.tag == 'mask':
                        spotProfiles[imageName][child.tag] = " ".join(child.text.split())
                    else:
                        spotProfiles[imageName][child.tag] = child.text.strip()
                if child.tag not in spot_profile_tags: spot_profile_tags.append (child.tag)
        
        regionalProfiles = []
        regionalProfileNodes = self.xmlnode.findall('.//regional_spot_profile_response')
        for regionalProfileNode in regionalProfileNodes:
            regionalProfile = ET.Element('profile_grid')
            regionalProfiles += [regionalProfile]
            regionalProfile.set('block', str(len(regionalProfiles)))
            for child in regionalProfileNode:
                if child.tag == 'number_of_profiles_x':
                    regionalProfile.set('num_x',child.text)
                elif child.tag == 'number_of_profiles_y':
                    regionalProfile.set('num_y', child.text)
                elif child.tag == 'profile':
                    subProfile = ET.SubElement(regionalProfile,'regional_profile')
                    for grandChild in child:
                        if grandChild.tag in ['box','width','height']:
                            subProfile.set(grandChild.tag, grandChild.text)
                        elif grandChild.tag == 'original':
                            for greatGrandChild in grandChild:
                                if greatGrandChild.tag == 'data':
                                    subProfile.set('raw_data', " ".join(greatGrandChild.text.split()))
                                    subProfile.set('max_data_value',str(max([int(value) for value in greatGrandChild.text.strip().split()])))
                                elif greatGrandChild.tag == 'mask': subProfile.set('mask', " ".join(greatGrandChild.text.split()))
                        elif grandChild.tag == 'averaged':
                            for greatGrandChild in grandChild:
                                if greatGrandChild.tag == 'data':
                                    subProfile.set('alt_raw_data', " ".join(greatGrandChild.text.split()))
                                elif greatGrandChild.tag == 'mask': subProfile.set('alt_mask', " ".join(greatGrandChild.text.split()))
                    if 'alt_raw_data' in subProfile.attrib:
                        if 'raw_data' in subProfile.attrib:
                            subProfile.set('type','dual')
                        else:
                            subProfile.set('type','averaged_only')
                            subProfile.set('raw_data','')
                            subProfile.set('mask','')
                    else:
                        subProfile.set('type','original_only')
                        subProfile.set('alt_raw_data','')
                        subProfile.set('alt_mask','')

        # Need to add image by image data analyses and histograms.  First off pass through and digest into the growing memory data structure
        integrationNodes = self.xmlnode.findall('.//integration_response')
        for iIntegrationNode in range(len(integrationNodes)):
            integrationNode = integrationNodes[iIntegrationNode]
            imageName = integrationNode.findall('image_name')[0].text
            if not imageName in collectedData:
                print(imageName)
                collectedData[imageName] = {}
            for child in integrationNode:
                if child.tag == 'image_name':
                    #no action
                    pass
                elif child.tag == 'status':
                    collectedData[imageName][child.tag] = {}
                    for grandChild in child:
                        collectedData[imageName][child.tag][grandChild.tag] = grandChild.text
                elif child.tag.startswith('bin'):
                    if not 'bins' in collectedData[imageName]: collectedData[imageName]['bins'] = []
                    collectedData[imageName]['bins'].append({})
                    for grandChild in child:
                        collectedData[imageName]['bins'][-1][grandChild.tag] = grandChild.text


        #And then turn this into additional XML
        #Get sorted list of image names
        imageNames = sorted([imageName for imageName in collectedData])
        integrationResultNode = ET.Element('integration_result')
        integrationResultNode.set('image_files_being_processed',' '.join([imageName for imageName in imageNames]))
        for tag in integration_positional_refinement_tags:
            integrationResultNode.set(tag,' '.join([collectedData[imageName][tag] for imageName in imageNames]))
        for tag in integration_post_refinement_tags:
            dataArray = []
            for imageName in imageNames:
                if tag in collectedData[imageName]: dataArray += [collectedData[imageName][tag]]
                else:
                    dataArray += ['{}']
            integrationResultNode.set(tag,' '.join(dataArray))
        nameMappingTotal = {"profile_fitted_fulls":"mean_profile_fitted_fulls",
            "profile_fitted_partials":"mean_profile_fitted_partials",
            "summation_integration_fulls":"mean_summation_integration_fulls",
            "summation_integration_partials":"mean_summation_integration_partials",
            "spot_count_fulls":"total_spot_count_fulls",
            "spot_count_partials":"total_spot_count_partials"}
        for tag in nameMappingTotal:
            integrationResultNode.set(nameMappingTotal[tag],' '.join([collectedData[imageName]['bins'][-1][tag] for imageName in imageNames]))

        nameMappingOuter = {"profile_fitted_fulls":"outer_profile_fitted_fulls",
            "profile_fitted_partials":"outer_profile_fitted_partials",
            "summation_integration_fulls":"outer_summation_integration_fulls",
            "summation_integration_partials":"outer_summation_integration_partials",
            "spot_count_fulls":"outer_spot_count_fulls",
            "spot_count_partials":"outer_spot_count_partials"}
        for tag in nameMappingTotal:
            integrationResultNode.set(nameMappingOuter[tag],' '.join([collectedData[imageName]['bins'][-2][tag] for imageName in imageNames]))

        histogramText = ''
        for imageName in collectedData:
            nBins = len(collectedData[imageName]['bins'])
            resoHistogramText = imageName.strip()+',bin_resolution_limits {\u221e'+" ".join([collectedData[imageName]['bins'][i]['upper_limit'].strip() for i in range(nBins)]).strip()+'} '
            histogramText += resoHistogramText
            #print resoHistogramText
            for tag in nameMappingOuter:
                if tag in collectedData[imageName]['bins'][0]:
                    thisHistogramText = imageName.strip()+','+tag.strip()+' '+'{'+" ".join([collectedData[imageName]['bins'][i][tag].strip() for i in range(nBins)]).strip()+'} '
                    #print thisHistogramText
                    histogramText += thisHistogramText
        integrationResultNode.set('histograms', histogramText)
        #Get sorted list of spot profiles
        imageNames = sorted([imageName for imageName in spotProfiles])
        profileTagLookup = {'width':'width','height':'height', 'profile':'raw_data', 'mask':'mask','image_name':'label'}
        for imageName in imageNames:
            newProfile = ET.SubElement(integrationResultNode,'profile')
            for tag in spot_profile_tags:
                newProfile.set(profileTagLookup[tag], spotProfiles[imageName][tag].strip())
            newProfile.set('max_data_value',str(max([int(value) for value in newProfile.get('raw_data').strip().split()])))
        
        #Append Block Profiles
        for regionalProfile in regionalProfiles:
            integrationResultNode.append(regionalProfile)


                                                                       
        return integrationResultNode
