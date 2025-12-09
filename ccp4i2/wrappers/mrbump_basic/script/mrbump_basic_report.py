from report.CCP4ReportParser import *

class mrbump_basic_report(Report):
  # Specify which gui task and/or pluginscript this applies to
  TASKNAME = 'mrbump_basic'
  RUNNING = True
  CSS_VERSION = '0.1.0'

  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)
  
    from ccp4i2.core import CCP4Utils
    import os
    
    results = self.addResults()
    results.append( 'MrBUMP is a pipeline to trial many search models in molecular replacement. \
                     It will find and prepare possible search models based on a sequence alignment \
                     between your target sequence and that of known structures in the PDB (using \
                     Phmmer by default). It then goes on to run molecular replacement (MR) on the \
                     best of these. Post MR it refine each resulting solution and do some initial \
                     model building to assess the likely success or otherwise of the molecular \
                     replacement search.' )

    #tableFoldsearch = results.addFold(label='sequence-based search results', initiallyOpen=True)
    #tableFoldsearch.append('These are the results of the sequence based search (using Phmmer). Below is the alignment of the best matches to the target sequences.<br/>')

    jobDirectory = jobInfo['fileroot']
    
    #if not os.path.isfile(os.path.join(jobDirectory, "search_mrbump_1", "logs", "alignment_report.log")):
    #    results.append("Warning: Cannot find alignment report file!")
    #else: 
    #    alog=open(os.path.join(jobDirectory, "search_mrbump_1", "logs", "alignment_report.log"), "r")
    #    lines="".join(alog.readlines())
    #    alog.close()
    #    tableFoldsearch.append("<pre>%s</pre>" % lines)

    #stuff=tableFoldsearch.addTable(select=".//MRBUMP/target_details/Sequence_align_list")
    #for title, select in [ ["Chain ID", "alignment/chainID"],
    #                       ["Alignment", "alignment/sequence"] ] :
    #    stuff.addData(title=title, select=select)
#
#    tableFoldx = results.addFold(label='detailed results for MR model search and preparation', initiallyOpen=True)
#    tableFoldx.append('Target Sequence:')
#    tableFoldx.append('Scores:')
#    tableFoldx.append('List of search results:')
#    tableFoldx.append('List of prepared models:')
#    tablex =  tableFoldx.addTable(select=".//MRBUMP/model_name/mrprogram")
#    for title, select in [ ["Model Name", "SearchModel_name"] ] :
#        tablex.addData(title=title, select=select)

    tableFoldsearch = results.addFold(label='search model preparation', initiallyOpen=True)
    tableFoldsearch.append('These are the search models that have been found and prepared for use in Molecular Replacement.<br/>')

    if not os.path.isfile(os.path.join(jobDirectory, "search_mrbump_1", "results", "models.txt")):
        results.append("Models will be listed shortly..")
    else: 
        alog=open(os.path.join(jobDirectory, "search_mrbump_1", "results", "models.txt"), "r")
        lines="".join(alog.readlines())
        alog.close()
        tableFoldsearch.append("<pre>%s</pre>" % lines)

    tableFoldmr = results.addFold(label='detailed results for molecular replacement and refinement', initiallyOpen=True)

    tableFoldmr.append('Model names have the following format: PDB ID_chain/domain ID_Search Source_MR preparation method_sequence ID_residue range in target \
                      Each model is used in Phaser to do molecular replacement. \
                      The resulting MR solution is then refined with Refmac (Final R, Final R)<br/>')

    if not os.path.isfile(os.path.join(jobDirectory, "search_mrbump_1", "results", "results.txt")):
        results.append("Molecular replacement results will appear here soon...")
    else: 
        alog=open(os.path.join(jobDirectory, "search_mrbump_1", "results", "results.txt"), "r")
        lines="".join(alog.readlines())
        #lines=alog.readlines()
        alog.close()

        tableFoldmr.append("<pre>%s</pre>" % lines)

#    table1 =  tableFoldmr.addTable(select=".//MRBUMP/model_name/mrprogram")
#
#    for title, select in [ ["Model Name", "SearchModel_name"],
#                           ["Phaser LLG", "PHASER_LLG"],
#                           ["Phaser TFZ", "PHASER_TFZ"],
#                           ["Final R<sub>fact</sub>", "final_Rfact"],
#                           ["Final R<sub>free</sub>", "final_Rfree"],
#                           ["SHELXE CC", "SHELXE_CC"] ] :
#
#        table1.addData(title=title, select=select)

    if os.path.isfile(os.path.join(jobDirectory, "search_mrbump_1", "logs", "programs.json")):
        tableFoldreferences = results.addFold(label='references', initiallyOpen=True)
        tableFoldreferences.append('The following programs were used in this run:') 
        
        import json
        with open(os.path.join(jobDirectory, "search_mrbump_1", "logs", "programs.json")) as json_file:
            programsUsed = json.load(json_file)

        from mrbump.initialisation import MRBUMP_master
        references=MRBUMP_master.References()
        for program in programsUsed:
            rprog=references.getReference(program)
            tableFoldreferences.append("<b>" + rprog.name + "</b> : " + rprog.paper)
            #tableFoldreferences.append(rprog.paper)

    self.addTaskReferences()

if __name__ == "__main__":
  import sys
  mrbump_basic_report(xmlFile=sys.argv[1],jobId=sys.argv[2])


