import xml.etree.ElementTree as etree

from ccp4i2.core.CCP4ClipperUtils import is_aminoacid
from ccp4i2.report import Report


class edstats_report(Report):
  TASKNAME = 'edstats'
  CSS_VERSION = '0.1.0'

  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)


    results = self.addResults()

    results.append ( "The ZD scores are accuracy metrics, <i>i.e.</i> at least in theory they can be improved by adjusting the model (by eliminating the " +\
                     "obvious difference density). The ZO scores are precision metrics and will be strongly correlated with the B<sub>iso</sub>s " +\
                     "(since that is also a precision metric), <i>i.e.</i> assuming you've fixed any issues with accuracy of that residue there's " +\
                     "nothing you can do about the precision, short of re-collecting the data. <br/>The advisable rejection limits and indeed default values " +\
                     "for this task are &#60; -3&#963; and &#62; 3&#963; for the residue " +\
                     "ZD metrics respectively, and &#60; 1&#963; for the residue ZO metrics." )

    list_of_chains = [ ] # We are going to separate results by chains

    new_tree = etree.Element('edstats_report')
    protein_fragments = etree.SubElement ( new_tree, 'protein' )
    ligands = etree.SubElement ( new_tree, 'ligands')
    waters = etree.SubElement ( new_tree, 'waters' )

    for residue in self.xmlnode.findall ( ".//Residue" ):
        if is_aminoacid ( residue.find("Name").text ) :
            chain_id = residue.find ( "Chain" ).text
            residue.remove ( residue.find ( "Chain" ) )
            if chain_id not in list_of_chains :
                list_of_chains.append ( chain_id )
                chain = etree.SubElement ( protein_fragments, 'chain', id=chain_id )
                chain.append(residue)
            else :
                chain.append(residue)
        elif residue.find("Name").text == "HOH" :
            waters.append ( residue )
        else :
            ligands.append ( residue )

    # Append new_tree to self.xmlnode, checking if it has an append method
    if hasattr(self.xmlnode, 'append') and callable(getattr(self.xmlnode, 'append')):
        self.xmlnode.append(new_tree)
    else:
        # Assume it's an ElementTree and get the root
        self.xmlnode.getroot().append(new_tree)

    for chain in list_of_chains :

        n_res = len(self.xmlnode.findall ( ".//edstats_report/protein/chain[@id='%s']/Residue" % chain ))
        plot_width = '10'

        if n_res < 150 :
            plot_width = '1'
        elif n_res < 500 :
            plot_width = '4'

        residues = self.xmlnode.findall ( ".//edstats_report/protein/chain[@id='%s']/Residue" % chain )
        bad_mc_positives = []
        bad_mc_negatives = []
        bad_sc_positives = []
        bad_sc_negatives = []
        bad_mc_precision = []
        bad_sc_precision = []
        for e in residues:
           ZDpm = e.findall("ZDpm")[0]
           ZDmm = e.findall("ZDmm")[0]
           ZDps = e.findall("ZDps")[0]
           ZDms = e.findall("ZDms")[0]
           ZOm = e.findall("ZOm")[0]
           ZOs = e.findall("ZOs")[0]
           if ZDpm.text != "n/a" and float(ZDpm.text)>3.0:
               bad_mc_positives.append(e)
           if ZDmm.text != "n/a" and float(ZDmm.text)<-3.0:
               bad_mc_negatives.append(e)
           if ZDps.text != "n/a" and float(ZDps.text)>3.0:
               bad_sc_positives.append(e)
           if ZDms.text != "n/a" and float(ZDms.text)<-3.0:
               bad_sc_negatives.append(e)
           if ZOm.text != "n/a" and float(ZOm.text)>1.0:
               bad_mc_precision.append(e)
           if ZOs.text != "n/a" and float(ZOs.text)>1.0:
               bad_sc_precision.append(e)

        results.append ( '<br/><h2>Chain %s</h2>' % chain )

        graph = results.addFlotGraph ( title="Chain %s" % chain,
                                       select=".//edstats_report/protein/chain[@id='%s']/Residue" % chain,
                                       style="height:350px; width:600px; border:0px; float:right;" )

        graph.addData ( title="Residue", select="Number" )
        graph.addData ( title="ZO&nbsp;(right&nbsp;axis)",      select="ZOm"    )
        graph.addData ( title="ZD-",     select="ZDmm"   )
        graph.addData ( title="ZD+",     select="ZDpm"   )
        graph.addData ( title="ZO",      select="ZOs"    )
        graph.addData ( title="ZD-",     select="ZDms"   )
        graph.addData ( title="ZD+",     select="ZDps"   )

        p = graph.addPlotObject()
        p.append( 'title', 'Main chain (%s)' % chain )
        p.append( 'plottype', 'xy' )
        p.append( 'showlegend', 'true' )
        p.append( 'xintegral', 'true' )
        p.append( 'yintegral', 'false' )
        p.append( 'xlabel', 'Residue' )
        p.append( 'ylabel', 'Z-scores' )
        p.append( 'yrange',  rightaxis='true', min = '0.0', max = '20.0' )
        p.append( 'yrange',  rightaxis='false', min = '-8.0', max = '6.0' )
        #plot_tree = p.as_etree()
        #etree.SubElement( plot_tree, 'polygon', linecolour="#cccccc", fillcolour="#aaaaaa", alpha="0.2" ).text = "0 5 10 15 20 25 30 35"
        #rect = p.append( 'polygon', linecolour="#cccccc", fillcolour="#aaaaaa", alpha="0.2" )

        l = p.append( 'plotline', xcol = 1, ycol = 2, rightaxis=True  ) # ZOm
        l.append( 'colour', 'black' )
        l.append( 'linestyle', '.' )
        l.append( 'symbol', 'x' )
        l.append( 'symbolsize', '2' )

        l = p.append( 'barchart', col = 1, tcol = 3 ) # ZDmm
        l.append( 'colour', 'red' )
        l.append ( 'width' ).text = plot_width

        l = p.append( 'barchart', col = 1, tcol = 4 ) # ZDpm
        l.append('colour','green')
        l.append ( 'width' ).text = plot_width

        p = graph.addPlotObject()
        p.append( 'title', 'Side chains (%s)' % chain )
        p.append( 'plottype', 'xy' )
        p.append( 'showlegend', 'true' )
        p.append( 'xintegral', 'true' )
        p.append( 'yintegral', 'false' )
        p.append( 'xlabel', 'Residue' )
        p.append( 'ylabel', 'Z-scores' )
        p.append( 'xrange',  min = '0.0' )
        p.append( 'yrange',  rightaxis='true', min = '0.0', max = '20.0' )
        p.append( 'yrange',  rightaxis='false', min = '-8.0', max = '6.0' )

        l = p.append( 'plotline', xcol = 1, ycol = 5, rightaxis=True  ) # ZOm
        l.append( 'colour', 'black' )
        l.append( 'linestyle', '.' )
        l.append( 'symbol', 'x' )
        l.append( 'symbolsize', '2' )

        l = p.append( 'barchart', col = 1, tcol = 6 ) # ZDmm
        l.append( 'colour', 'red' )
        l.append ( 'width' ).text = plot_width

        l = p.append( 'barchart', col = 1, tcol = 7 ) # ZDpm
        l.append('colour','green')
        l.append ( 'width' ).text = plot_width

        tab_div = results.addDiv (style="float:left;")
        tab_div.append ( "<p >Main chain: </p>" )
        bad_mc_negatives+bad_mc_positives
        table = tab_div.addTable(selectNodes = bad_mc_negatives+bad_mc_positives,
                                        id='bad_mc_positives' )

        for title,subtitle,select in [ [ "Residue",    "Name" ,  "Name"   ],
                                       [ "Residue",    "Number", "Number" ],
                                       [ "Main chain", "ZO",     "ZOm"    ],
                                       [ "Main chain", "ZD-",    "ZDmm"   ],
                                       [ "Main chain", "ZD+",    "ZDpm"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)

        tab_div.append ( "<p>Side chains: </p>" )

        table = tab_div.addTable(selectNodes = bad_sc_negatives+bad_sc_positives,
                                        id='bad_sc_positives' )

        for title,subtitle,select in [ [ "Residue",    "Name" ,  "Name"   ],
                                       [ "Residue",    "Number", "Number" ],
                                       [ "Side chain", "ZO",     "ZOs"    ],
                                       [ "Side chain", "ZD-",    "ZDms"   ],
                                       [ "Side chain", "ZD+",    "ZDps"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)


        results.addDiv ( style="clear:both;" )

    if len(waters) > 0 :

        waters_folder = results.addFold(label='Correlation analysis: water molecules', initiallyOpen=False)

        n_res = len(self.xmlnode.findall (".//edstats_report/waters/Residue" ))

        plot_width = '10'

        if n_res < 150 :
            plot_width = '1'
        elif n_res < 500 :
            plot_width = '4'

        graph = waters_folder.addFlotGraph ( title="Waters", select=".//edstats_report/waters/Residue",
                                             style="height:320px; width:550px; border:0px; float:right;" )
        graph.addData ( title="Residue", select="Number" )
        graph.addData ( title="ZO&nbsp;(right&nbsp;axis)",     select="ZOa"    )
        graph.addData ( title="ZD-",    select="ZDma"   )
        graph.addData ( title="ZD+",    select="ZDpa"   )

        p = graph.addPlotObject()
        p.append( 'title', 'Individual waters (all chains, average metrics)' )
        p.append( 'plottype', 'xy' )
        p.append( 'showlegend', 'true' )
        p.append( 'xintegral', 'true' )
        p.append( 'yintegral', 'false' )
        p.append( 'xlabel', 'Water molecule' )
        p.append( 'ylabel', 'Z-scores' )
        p.append( 'xrange',  min = '0.0' )
        p.append( 'yrange',  rightaxis='true', min = '0.0', max = '20.0' )
        p.append( 'yrange',  rightaxis='false', min = '-8.0', max = '6.0' )

        l = p.append( 'plotline', xcol = 1, ycol = 2, rightaxis=True  ) # ZOm
        l.append( 'colour', 'black' )
        l.append( 'linestyle', '.' )
        l.append( 'symbol', 'x' )
        l.append( 'symbolsize', '2' )

        l = p.append( 'barchart', col = 1, tcol = 3 ) # ZDmm
        l.append( 'colour', 'red' )
        l.append ( 'width' ).text = plot_width

        l = p.append( 'barchart', col = 1, tcol = 4 ) # ZDpm
        l.append('colour','green')
        l.append ( 'width' ).text = plot_width

        tab_div = waters_folder.addDiv (style="float:left; width:28%; padding-right:50px;")
        tab_div.append ( "<p >The entries in the following table may not be water molecules but part of a bigger, undetermined structure (a ligand, perhaps?): </p>" )
        table = tab_div.addTable(select=".//edstats_report/waters/Residue[ZDpa>3.0]",
                                        id='bad_water_positives' )

        for title,subtitle,select in [ [ "Residue",    "Chain" ,  "Chain"  ],
                                       [ "Residue",    "Number",  "Number" ],
                                       [ "Main chain", "ZO",      "ZOm"    ],
                                       [ "Main chain", "ZD-",     "ZDma"   ],
                                       [ "Main chain", "ZD+",     "ZDpa"   ],
                                       [ "Main chain", "RSCC",    "CCPa"   ],
                                       [ "Main chain", "Biso",    "BAa"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)

        tab_div.append ( "<p>The following entries sit in weak density or next to noise peaks. In the former case, removing them might improve R-free</p>" )

        table = tab_div.addTable(select=".//edstats_report/waters/Residue[ZDma<-3.0]",
                                        id='bad_water_negatives' )

        for title,subtitle,select in [ [ "Residue",    "Chain" ,  "Chain"  ],
                                       [ "Residue",    "Number",  "Number" ],
                                       [ "Main chain", "ZO",      "ZOm"    ],
                                       [ "Main chain", "ZD-",     "ZDma"   ],
                                       [ "Main chain", "ZD+",     "ZDpa"   ],
                                       [ "Main chain", "RSCC",    "CCPa"   ],
                                       [ "Main chain", "Biso",    "BAa"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)

        results.addDiv ( style="clear:both;" )


    if len(ligands) > 0 :

        ligands_folder = results.addFold(label='Correlation analysis: ligands and other solutes', initiallyOpen=False)

        n_res = len(self.xmlnode.findall (".//edstats_report/ligands/Residue" ))

        plot_width = '1'

        graph = ligands_folder.addFlotGraph ( title="Ligands", select=".//edstats_report/ligands/Residue",
                                              style="height:300px; width:500px; border:0px; float:right;" )

        graph.addData ( title="ZO",     select="ZOa"     )
        graph.addData ( title="ZD-",    select="ZDma"    )
        graph.addData ( title="ZD+",    select="ZDpa"    )
        graph.addData ( title="Biso",   select="BAa"     )
        graph.addData ( title="RSCC",   select="CCPa"    )

        p = graph.addPlotObject()
        p.append( 'title', 'Ligands and solutes - accuracy vs precision (all chains, average metrics)' )
        p.append( 'plottype', 'xy' )
        p.append( 'showlegend', 'true' )
        p.append( 'xintegral', 'false' )
        p.append( 'yintegral', 'false' )
        p.append( 'xlabel',
                  '(density is less precise)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' +\
                  'ZO&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(density is more precise)' )
        p.append( 'ylabel', 'ZD scores' )
        p.append( 'xrange',  min = '0.0' )
        p.append( 'yrange',  rightaxis='true', min = '0.0', max = '250.0' )
        p.append( 'yrange',  rightaxis='false', min = '-30.0')#, max = '20.0' )

        l = p.append( 'barchart', col = 1, tcol = 2 ) # ZDmm
        l.append( 'colour', 'red' )
        l.append ( 'width' ).text = plot_width

        l = p.append( 'barchart', col = 1, tcol = 3 ) # ZDpm
        l.append('colour','green')
        l.append ( 'width' ).text = plot_width

        l = p.append( 'plotline', xcol = 1, ycol = 4, rightaxis=True ) # BAa
        l.append( 'colour', 'blue')
        l.append( 'linestyle', '.' )
        l.append( 'symbolsize', '3' )
        l.append( 'symbol', 'o' )


        p = graph.addPlotObject()
        p.append( 'title', 'Ligands and solutes - RSCC vs Biso (all chains, average metrics)' )
        p.append( 'plottype', 'xy' )
        p.append( 'showlegend', 'true' )
        p.append( 'xintegral', 'false' )
        p.append( 'yintegral', 'false' )
        p.append( 'xlabel', 'Isotropic B factor' )
        p.append( 'ylabel', 'Real Space Correlation Coefficient' )
        p.append( 'xrange',  min = '0.0' )
        #p.append( 'yrange',  rightaxis='true', min = '0.0', max = '300.0' )
        p.append( 'yrange',  rightaxis='false', min = '0.0', max = '1.0' )

        l = p.append( 'plotline', xcol = 4, ycol = 5 ) # BAa
        l.append( 'colour', 'blue')
        l.append( 'linestyle', '.' )
        l.append( 'symbolsize', '3' )
        l.append( 'symbol', 'o' )

        tab_div = ligands_folder.addDiv (style="float:left; width:30%; padding-right:50px;")
        tab_div.append ( "<p >The ligand models in the following table may be incomplete or fit density poorly: </p>" )
        table = tab_div.addTable(select=".//edstats_report/ligands/Residue[ZDpa>3.0]",
                                        id='bad_ligand_positives' )

        for title,subtitle,select in [ [ "Residue",    "Chain" ,  "Chain"  ],
                                       [ "Residue",    "Name",    "Name"   ],
                                       [ "Residue",    "Number",  "Number" ],
                                       [ "Main chain", "ZO",      "ZOm"    ],
                                       [ "Main chain", "ZD-",     "ZDma"   ],
                                       [ "Main chain", "ZD+",     "ZDpa"   ],
                                       [ "Main chain", "RSCC",    "CCPa"   ],
                                       [ "Main chain", "Biso",    "BAa"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)

        tab_div.append ( "<p>The following ligands sit in weak density or next to noise peaks. In the former case, removing them might improve R-free</p>" )

        table = tab_div.addTable(select=".//edstats_report/ligands/Residue[ZDma<-3.0]",
                                        id='bad_ligand_negatives' )

        for title,subtitle,select in [ [ "Residue",    "Chain" ,  "Chain"  ],
                                       [ "Residue",    "Name",    "Name"   ],
                                       [ "Residue",    "Number",  "Number" ],
                                       [ "Main chain", "ZO",      "ZOm"    ],
                                       [ "Main chain", "ZD-",     "ZDma"   ],
                                       [ "Main chain", "ZD+",     "ZDpa"   ],
                                       [ "Main chain", "RSCC",    "CCPa"   ],
                                       [ "Main chain", "Biso",    "BAa"   ]
                                     ] :


            table.addData(title=subtitle, subtitle=subtitle, select=select)




        results.addDiv ( style="clear:both;" )



    for chain in list_of_chains :

        tableFolder = self.addFold(label='Chain %s: complete per-residue metrics' % chain, initiallyOpen=False)
        table1 = tableFolder.addTable(select=".//edstats_report/protein/chain[@id='%s']/Residue" % chain, downloadable=True, id='edstats_table')

        for title,subtitle,select in [ [ "Residue",    "Name" ,  "Name"   ],
                                       [ "Residue",    "Number", "Number" ],
                                       [ "Main chain", "RSCC",   "CCPa"   ],
                                       [ "Main chain", "Biso",   "BAa"    ],
                                       [ "Main chain", "ZCCm",   "ZCCm"   ],
                                       [ "Main chain", "ZOm",    "ZOm"    ],
                                       [ "Main chain", "ZD-m",   "ZDmm"   ],
                                       [ "Main chain", "ZD+m",   "ZDpm"   ],
                                       [ "Side chain", "ZCCs",   "ZCCs"   ],
                                       [ "Side chain", "ZOs",    "ZOs"    ],
                                       [ "Side chain", "ZD-s",   "ZDms"   ],
                                       [ "Side chain", "ZD+s",   "ZDps"   ],
                                       [ "Average",    "ZCCa",   "ZCCa"   ],
                                       [ "Average",    "ZOa",    "ZOa"    ],
                                       [ "Average",    "ZD-a",   "ZDma"   ],
                                       [ "Average",    "ZD+a",   "ZDpa"   ] ] :


            table1.addData(title=subtitle, subtitle=subtitle, select=select)

    self.addTaskReferences()

if __name__ == "__main__":
  import sys
  edstats_report(xmlFile=sys.argv[1],jobId=sys.argv[2])
