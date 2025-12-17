
import os
import shutil
import xml.etree.ElementTree as etree

import numpy

from ccp4i2.core import CCP4Utils
from ccp4i2.report import Report


class privateer_report(Report):
  TASKNAME = 'privateer'
  CSS_VERSION = '0.2.0'

  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)

    projectid = self.jobInfo.get("projectid", None)
    jobNumber = self.jobInfo.get("jobnumber", None)

    imageFile = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=pyranose_pseudorotational.svg?jobNumber="
                + jobNumber
            )
    
    imageFileSrc = os.path.join(CCP4Utils.getCCP4I2Dir(), 'wrappers/privateer/script/pyranose_pseudorotational.svg' )
    imageFileJob = os.path.join(jobInfo['fileroot'],"pyranose_pseudorotational.svg")
    shutil.copyfile(imageFileSrc,imageFileJob)

    background_pyranoses = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=mercator_pyranoses.png?jobNumber="
                + jobNumber
            )
    background_pyranosesSrc = os.path.join(CCP4Utils.getCCP4I2Dir(), 'wrappers/privateer/script/mercator_pyranoses.png' )
    background_pyranosesJob = os.path.join(jobInfo['fileroot'], 'mercator_pyranoses.png' )
    shutil.copyfile(background_pyranosesSrc,background_pyranosesJob)

    results = self.addResults()
    results.append('The Cremer-Pople analysis (Cremer and Pople, 1975, JACS 97:1354-58) is used to determine sugar ring conformation. Below is a 2D plot of the conformational'+\
            ' parameters (Q, Phi, Theta for pyranoses; Q and Theta for furanoses) along with a depiction of the conformational sphere for pyranoses:')


    htmlCode = '<img src="' + imageFile + '" alt="Diagram" style="width:250px; height:250px; margin-top:20px; margin-right: 30px; float:left;" />'
    results.append(htmlCode)

    # What follows is an attempt to support a "traffic lights" representation system.
    # We will need to re-use Privateer MKIII XML until this is replaced with the newer
    #   python-wrapped XML generation in MKIV.
    #
    # This change to graphical representation was suggested by Dave Brown and Paul Emsley
    #   at a CCP4 WG2 meeting at Birkbeck in London (January 2019). Users nodded in agreement.

    new_tree = etree.Element('privateer_report')
    pyranose_tree = etree.SubElement ( new_tree, 'pyranoses' )
    furanose_tree = etree.SubElement ( new_tree, 'furanoses')

    low_energy_pyranoses   = etree.SubElement ( pyranose_tree, 'low_energy'   )
    high_energy_pyranoses  = etree.SubElement ( pyranose_tree, 'high_energy'  )
    other_issues_pyranoses = etree.SubElement ( pyranose_tree, 'other_issues' )

    for sugar in self.xmlnode.findall ( ".//ValidationData/Pyranose" ):
        if "conformation" in sugar.find("SugarDiagnostic").text.lower() :
            high_energy_pyranoses.append ( sugar )
        elif "Ok" in sugar.find("SugarDiagnostic").text :
            low_energy_pyranoses.append ( sugar )
        else : other_issues_pyranoses.append ( sugar )

    # Append new_tree to self.xmlnode, checking if it has an append method
    if hasattr(self.xmlnode, 'append') and callable(getattr(self.xmlnode, 'append')):
        self.xmlnode.append(new_tree)
    else:
        # Assume it's an ElementTree and get the root
        self.xmlnode.getroot().append(new_tree)

    graph = results.addFlotGraph(title="graphical representations", select=".//privateer_report/pyranoses", style="height:300px; width:450px; border:0px; float:left;" )
    graph.addData ( title="phi_pyranose_he"  , select="high_energy/Pyranose/SugarPhi"      )
    graph.addData ( title="theta_pyranose_he", select="high_energy/Pyranose/SugarTheta"    )
    graph.addData ( title="phi_pyranose_oi"  , select="other_issues/Pyranose/SugarPhi"     )
    graph.addData ( title="theta_pyranose_oi", select="other_issues/Pyranose/SugarTheta"   )
    graph.addData ( title="phi_pyranose_le"  , select="low_energy/Pyranose/SugarPhi"       )
    graph.addData ( title="theta_pyranose_le", select="low_energy/Pyranose/SugarTheta"     )
    graph.addData ( title="rscc_pyranose_he" , select="high_energy/Pyranose/SugarRSCC"     )
    graph.addData ( title="bfact_pyranose_he", select="high_energy/Pyranose/SugarBFactor"  )
    graph.addData ( title="rscc_pyranose_oi" , select="other_issues/Pyranose/SugarRSCC"    )
    graph.addData ( title="bfact_pyranose_oi", select="other_issues/Pyranose/SugarBFactor" )
    graph.addData ( title="rscc_pyranose_le" , select="low_energy/Pyranose/SugarRSCC"      )
    graph.addData ( title="bfact_pyranose_le", select="low_energy/Pyranose/SugarBFactor"   )

    # We're gonna need a separate plot for furanoses
    #graph.addData ( title="PhiFuranoses", select="Furanose/SugarPhi" )
    #graph.addData ( title="QFuranoses", select="Furanose/SugarQ" )
    #graph.addData ( title="Furanose", select="Furanose/SugarRSCC" )
    #graph.addData ( title="BfacFuranoses", select="Furanose/SugarBFactor" )


    ########## Pyranose conformation ###########
    p = graph.addPlotObject()
    p.append( 'description', 'This graph is a 2D projection of the Cremer-Pople sphere, with chairs at the poles, boats '+\
              'and skew-boats at the equator, and envelopes and half-chairs elsewhere.')
    p.append( 'title', 'Conformational landscape for pyranoses' )
    p.append( 'plottype', 'xy' )
    p.append( 'showlegend', 'false' )
    p.append( 'background' ).text = background_pyranoses
    p.append( 'xintegral', 'false' )
    p.append( 'yintegral', 'false' )
    p.append( 'xlabel', 'Phi' )
    p.append( 'ylabel', 'Theta' )
    p.append( 'xrange', min = '360.0', max = '0.0' )
    p.append( 'yrange', rightaxis='false', max = '0.0', min = '180.0' )

    l = p.append( 'customYLabels' )
    l.text = '0,45,90,135,180'

    l = p.append( 'customYTicks' )
    l.text = '0,45,90,135,180'

    l = p.append( 'customXLabels' )
    l.text = '0,30,60,90,120,150,180,210,240,270,300,330,360'

    l = p.append( 'customXTicks' )
    l.text = '0,30,60,90,120,150,180,210,240,270,300,330,360'

    l = p.append( 'plotline', xcol = 5, ycol = 6 )
    l.append( 'colour', 'olivedrab' )
    l.append( 'linestyle', '.' )
    l.append( 'symbol', '.' )
    l.append( 'symbolsize', '3' )
    l.append( 'fillcolour', 'darkgreen' )
    l = p.append( 'plotline', xcol = 3, ycol = 4 )
    l.append( 'colour', 'orange' )
    l.append( 'linestyle', '.' )
    l.append( 'symbolsize', '5' )
    l.append( 'symbol', 'o' )
    l.append( 'fillcolour', 'gold' )
    l = p.append( 'plotline', xcol = 1, ycol = 2 )
    l.append( 'colour', 'red' )
    l.append( 'linestyle', '.' )
    l.append( 'symbol', 'x' )
    l.append( 'symbolsize', '9' )

    pyranoses = self.xmlnode.findall(".//Pyranose")
    furanoses = self.xmlnode.findall(".//Furanose")

    ########## Furanose conformation VS Q ###########
    # if furanoses != "" :
    #     p = graph.addPlotObject()
    #     p.append( 'title', 'Furanose conformation VS puckering amplitude' )
    #     p.append( 'description', 'A plot of conformation (Theta) VS puckering amplitude (Q). Q measures how much atoms separate from the mean ring plane.')
    #     p.append( 'plottype', 'xy' )
    #     p.append( 'showlegend', 'false' )
    #     p.append( 'xintegral', 'true' )
    #     p.append( 'yintegral', 'true' )
    #     p.append( 'xrange', min = '0.0', max = '360.0' )
    #     p.append( 'yrange', rightaxis='false', min = '0.0', max = '1.0' )
    #     l = p.append( 'plotline', xcol = 3, ycol = 4 )
    #     l.append( 'colour', 'blue' )
    #     l.append( 'linestyle', '.' )

    ########## B-factor VS RSCC ###########

    list_of_bfac_elems = self.xmlnode.findall(".//SugarBFactor")
    list_of_bfac = [ ]

    for elem in list_of_bfac_elems :
        list_of_bfac.append( elem.text )

    bfac_array = numpy.array (list_of_bfac, dtype=float)
    max_bfac = numpy.amax(bfac_array)

    rounded_max = 100.0

    while rounded_max < max_bfac :
        rounded_max += 100.0

    p = graph.addPlotObject()
    p.append( 'title', '<B-factor> VS Real Space CC' )
    p.append( 'plottype', 'xy' )
    p.append( 'showlegend', 'false' )

    p.append( 'xintegral', 'false' )
    p.append( 'yintegral', 'false' )
    p.append( 'yrange', rightaxis='false', min = '0.0', max = '1.0' )
    p.append( 'xrange', min='0.0', max=str(rounded_max) )

    l = p.append( 'customYLabels' )
    l.text = '0.0,0.5,0.7,1.0'

    l = p.append( 'customYTicks' )
    l.text = '0.0,0.5,0.7,1.0'

    p.append( 'xlabel', '<Isotropic B-factor>' )
    p.append( 'ylabel', 'Real Space CC' )

    l = p.append( 'plotline', xcol = 12, ycol = 11 )
    l.append( 'colour', 'olivedrab' )
    l.append( 'linestyle', '.' )
    l.append( 'symbol', '.' )
    l.append( 'symbolsize', '3' )
    l.append( 'fillcolour', 'darkgreen' )
    l = p.append( 'plotline', xcol = 10, ycol = 9 )
    l.append( 'colour', 'orange' )
    l.append( 'linestyle', '.' )
    l.append( 'symbolsize', '5' )
    l.append( 'symbol', 'o' )
    l.append( 'fillcolour', 'gold' )
    l = p.append( 'plotline', xcol = 8, ycol = 7 )
    l.append( 'colour', 'red' )
    l.append( 'linestyle', '.' )
    l.append( 'symbol', 'x' )
    l.append( 'symbolsize', '9' )
    l = p.append ( 'polygon', linecolour="#cccccc", fillcolour="#222222", alpha = 0.2 )
    l.text = ( "0.0 0.0 0.0 0.7 %f 0.7 %f 0.0" % (rounded_max, rounded_max) )

    l = p.append( 'text', xpos=5, ypos=0.72, font='bold 1.25em Helvetica', colour="#AAAAAA" )
    l.text = "Low CC"

    results.append('<br/>')
    clearingDiv = results.addDiv(style="clear:both;")


    ########## Glycan structures ########
    trees = ""
    if len(self.xmlnode.findall ( ".//ValidationData/Glycan" ))>0:
        self.xmlnode.findall ( ".//ValidationData/Glycan" )[0].text

    if trees != "" :
        treesFold = results.addFold ( label="N- and O-glycan structure 2D descriptions", initiallyOpen=False )

        index = 0
        glycanindexforpermutation = 0
        permutationindex = 0 

        chain = 'empty'

        treesFold.append ( 'Below are graphical plots of the detected glycan trees. Placing your mouse pointer' +\
                          ' over any of the sugars will display a tooltip containing its residue name and number from the PDB file.' )
        directory = jobInfo['fileroot']
        for glycan in self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanSVG" ) :

            if chain != self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanChain" )[index].text :
                chain = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanChain" )[index].text
                treesFold.append ( '<p style="font-size:130%; padding:2px; margin-top:20px; margin-bottom:0px; font-weight:bold; margin-left:15px; clear:both"> Chain ' + chain + '</p>' )
                chaindiv = treesFold.addDiv(style="border-width: 1px; padding-top: 10px; padding-bottom:10px; border-color:black; border-style:solid; border-radius:15px;")
#                treesFold.append ( '<p style="background-color:rgb(145,203,219); padding:2px;' +\
#                        'margin-top:20px; margin-bottom:12px; clear:both"> Chain ' + chain + '</p>' )
            

            svg_filename = os.path.join ( directory, glycan.text )

            svg_file = open(svg_filename, 'r')
            svg_string = svg_file.read()
            svg_file.close()

            svg_string_partitioned = svg_string.partition("width=\"")
            svg_width = ''
            for symbol in svg_string_partitioned[2]:
                if symbol != "\"":
                    svg_width = svg_width + symbol
                else:
                    break

            #chaindiv.append ('<div style="float:right; padding-top: 5px; padding-bottom: 5px; margin-right:30px;>' + svg_string + '</div>' )
            svg_div = chaindiv.addDiv ( style="float:right; padding-right:10px; " )
            

            svg_div.append ( svg_string )

            wurcs = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanWURCS" )[index].text
            svg_div.append ( '<p style="font-size:130%; max-width:400px; font-weight:bold"> ' + wurcs + ' </p> ')

            try:
                if len(self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanGTCID" )[index].text) != 0:
                    gtcid = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanGTCID" )[index].text
                    glyconnectid = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanGlyConnectID" )[index].text
                    if gtcid != "Unable to find GlyTouCan ID":
                        svg_div.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold">GlyTouCan ID:<a href="https://glytoucan.org/Structures/Glycans/' + gtcid + '">' + gtcid + '</a></p> ')
                    else:
                        svg_div.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> GlyTouCan ID: ' + 'Not Found' + ' </p> ')
                    
                    if glyconnectid != "Unable to find GlyConnect ID":
                        svg_div.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold">GlyConnect ID:<a href="https://glyconnect.expasy.org/browser/structures/' + glyconnectid + '">' + glyconnectid + '</a></p> ')
                    else:
                        svg_div.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold">GlyConnect ID: ' + 'Not Found' + ' </p> ')
                        
                        if self.xmlnode.findall ( ".//ValidationData/Glycan" )[index].find('GlycanPermutations') != None:
                            permutationsFold = svg_div.addFold ( label="Closest permutations detected on GlyConnect database", initiallyOpen=False )
                            permutationsdiv = permutationsFold.addDiv(style="border-width: 1px; padding-top: 10px; padding-bottom:10px; border-color:black; border-style:solid; border-radius:15px;")
                            
                            for permutation in self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations" )[glycanindexforpermutation]:
                                wurcspermutation = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/PermutationWURCS" )[permutationindex].text
                                gtcidpermutation = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/PermutationGTCID" )[permutationindex].text
                                glyconnectidpermutation = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/PermutationGlyConnectID" )[permutationindex].text
                                permutationsvg = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/PermutationSVG" )[permutationindex].text
                                permutationscore = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/PermutationScore" )[permutationindex].text
                                anomerpermutations = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/anomerPermutations" )[permutationindex].text
                                residuepermutations = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/residuePermutations" )[permutationindex].text
                                residuedeletions = self.xmlnode.findall ( ".//ValidationData/Glycan/GlycanPermutations/GlycanPermutation/residueDeletions" )[permutationindex].text

                                svg_filename_permutation = os.path.join ( directory, permutationsvg )

                                svg_file_permutation = open(svg_filename_permutation, 'r')
                                svg_string_permutation = svg_file_permutation.read()
                                svg_file_permutation.close()

                                # svg_string_partitioned_permutation = svg_string_permutation.partition("width=\"")
                                # svg_width_permutation = ''
                                # for symbol in svg_string_partitioned_permutation[2]:
                                #     if symbol != "\"":
                                #         svg_width_permutation = svg_width_permutation + symbol
                                #     else:
                                #         break

                                permutationsdiv.append ( svg_string_permutation )
                                
                                permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> ' + wurcspermutation + ' </p> ')
                                if float(permutationscore) <= 1.00:
                                    permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> Permutation Score(out of 100): <span style="color: #00ff00">' + permutationscore + '</span></p> ')
                                elif float(permutationscore) > 1.00 and float(permutationscore) <= 10.00:
                                    permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> Permutation Score(out of 100): <span style="color: #ffa500">' + permutationscore + '</span></p> ')
                                elif float(permutationscore) > 10.00:
                                    permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> Permutation Score(out of 100): <span style="color: #ff3300">' + permutationscore + '</span></p> ')
                                permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold"> Anomer Permutations: ' + anomerpermutations + '<br>Residue Permutations: ' + residuepermutations + '<br>Residue Deletions: ' + residuedeletions + '</br></br></p> ')
                                permutationsdiv.append ( '<p style="font-size:130%; max-width:' + svg_width + 'px; font-weight:bold">GlyTouCan ID:<a href="https://glytoucan.org/Structures/Glycans/' + gtcidpermutation + '">' + gtcidpermutation + '</a>' 
                                + '<br>GlyConnect ID:' + '<a href="https://glyconnect.expasy.org/browser/structures/' + glyconnectidpermutation + '">' + glyconnectidpermutation + '</a></br></p>')
                                
                                svg_div_permutation = permutationsdiv.addDiv ( style="float:right; padding-right:10px; " )
                                permutationindex = permutationindex + 1
                    
                            glycanindexforpermutation = glycanindexforpermutation + 1
            except IndexError:
                pass
            
            index = index + 1
                    


#        treeshtmlCode = '<img src="' + treesName + '" alt="Diagram" />'#width="415" height="183" />'
#        treesFold.append ( treeshtmlCode )
            chaindiv.append ( '<div style="clear:both" />' )

        treesFold.append ( '<div style="clear:both; margin-bottom:30px;" />' )

    ######### Detailed validation data ######
    tableFold = results.addFold(label='Detailed monosaccharide validation data', initiallyOpen=True)

    if pyranoses != "" :

        tableFold.append( 'Validation results for pyranose sugars: <br/>' )
        table1 = tableFold.addTable(select=".//Pyranose")

        for title,select in [ [ "Chain"  , "SugarChain" ],
                              [ "Name" ,"SugarName" ],
                              [ "Q<sup>1</sup>" , "SugarQ"],
                              [ "Phi" , "SugarPhi" ],
                              [ "Theta" , "SugarTheta" ],
                              [ "Anomer" , "SugarAnomer" ],
                              [ "D/L<sup>2</sup>", "SugarHand" ],
                              [ "Conformation" ,"SugarConformation" ],
                              [ "RSCC<sup>3</sup>" , "SugarRSCC" ],
                              [ "&lt;Bfactor&gt;", "SugarBFactor" ],
                              [ "Diagnostic" , "SugarDiagnostic" ] ] :

            table1.addData(title=title,select=select)

    if furanoses != "" :

        tableFold.append( 'Validation results for furanose sugars: <br/>' )
        table2 = tableFold.addTable(select=".//ValidationData/Furanose")

        for title,select in [ [ "Chain"  , "SugarChain" ],
                              [ "Name" ,"SugarName" ],
                              [ "Q" , "SugarQ"],
                              [ "Phi" , "SugarPhi" ],
                              [ "Anomer" , "SugarAnomer" ],
                              [ "D/L", "SugarHand" ],
                              [ "Conformation" ,"SugarConformation" ],
                              [ "RSCC" , "SugarRSCC" ],
                              [ "&lt;Bfactor&gt;", "SugarBFactor" ],
                              [ "Diagnostic" , "SugarDiagnostic" ] ] :

            table2.addData(title=title,select=select)

    tableFold.append('<sup>1</sup>Q is the total puckering amplitude, measured in Angstroems.<br/>'+\
                '<sup>2</sup>Whenever N is displayed in the D/L column, it means that Privateer has been unable to determine the handedness based solely on the structure.<br/>'+\
        '<sup>3</sup>RSCC, short for Real Space Correlation Coefficient, measures the agreement between model and positive omit density. A RSCC below 0.8 is tipically considered poor.')


    ############ Summary ##############
    summaryFold = self.addFold ( label='Summary for publications', initiallyOpen=True )

    ngeom = nconf = nanom = ndl = nq = 0

    for diagnostics in self.xmlnode.findall(".//ValidationData/Pyranose/SugarDiagnostic") :
        ngeom += diagnostics.text.count ( "geometry" )
        nconf += diagnostics.text.count ( "mistaken" )
        nanom += diagnostics.text.count ( "anomer" )
        ndl   += diagnostics.text.count ( "configuration" )
        nq    += diagnostics.text.count ( "Q=" )

    for diagnostics in self.xmlnode.findall(".//ValidationData/Furanose/SugarDiagnostic") :
        ngeom += diagnostics.text.count ( "geometry" )
        nconf += diagnostics.text.count ( "mistaken" )
        nanom += diagnostics.text.count ( "anomer" )
        ndl   += diagnostics.text.count ( "configuration" )
        nq    += diagnostics.text.count ( "Q=" )

    tableSummary = summaryFold.addTable ( transpose=True )
    #tableSummary.addData ( title = "Check ring geometry", data = [ ngeom ] )
    tableSummary.addData ( title = "Stereochemical problems", data = [ nanom + ndl ] )
    tableSummary.addData ( title = "Unphysical puckering amplitude", data = [ nq ] )
    tableSummary.addData ( title = "In unlikely ring conformation", data = [ nconf ] )

#    svgFile = '<svg width="100" height="100"><circle cx="50" cy="50" r="40" stroke="green" stroke-width="4" fill="yellow" /></svg>'

#    results.append(svgFile)

    totalIssues = ngeom + nconf + nanom + ndl + nq

    if totalIssues > 0 :
        summaryFold.append ( "There are " + str(totalIssues) + " issues with your model, and they may need correcting. <br/>" )
        summaryFold.append ( "Please select 'Manual model rebuilding' at the bottom of this window, so Coot can guide you through all the problems raised.<br/>" )
    else :
        summaryFold.append ( "No issues were found while analysing the data." )

    self.addTaskReferences()


#    def unescape(s):
#        s = s.replace("&lt;", "<")
#        s = s.replace("&gt;", ">")
#        s = s.replace("&amp;", "&")
#        return s

if __name__ == "__main__":
  import sys
  privateer_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
