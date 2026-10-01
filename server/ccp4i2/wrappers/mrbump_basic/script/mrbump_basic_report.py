import os
import re

from ccp4i2.report import Report


def final_solutions(results_text):
    """The rows of MrBUMP's "Final MR solution from Phaser" table, as dicts.

    MrBUMP writes its results as fixed-width text too wide for the report's
    box, which cut off the numbers that decide whether it worked (TFZ, the
    space group Phaser chose, R and R-free). Columns are grouped by "|":
    "# name copy rid | eLLG seqid cover | res expt | RFZ TFZ LLG SG | R Rfree".
    """
    if 'Final MR solution' not in (results_text or ''):
        return []
    block = results_text.split('Final MR solution', 1)[1]
    rows = []
    for line in block.splitlines():
        groups = [g.split() for g in line.split('|')]
        if len(groups) != 5 or not groups[0] or not groups[0][0].isdigit():
            continue
        try:
            rows.append({'model': groups[0][1], 'rfz': groups[3][0], 'tfz': groups[3][1],
                         'llg': groups[3][2], 'sg': groups[3][3],
                         'r': groups[4][0], 'rfree': groups[4][1]})
        except IndexError:
            continue
    return rows


class mrbump_basic_report(Report):
  TASKNAME = 'mrbump_basic'
  RUNNING = True
  # The report reads MrBUMP's results/*.txt; MrBUMP writes no program.xml.
  USEPROGRAMXML = False
  CSS_VERSION = '0.1.0'

  def __init__(self,xmlnode=None,jobInfo={},**kw):
    Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)

    results = self.addResults()
    results.append( 'MrBUMP is a pipeline to trial many search models in molecular replacement. \
                     It will find and prepare possible search models based on a sequence alignment \
                     between your target sequence and that of known structures in the PDB (using \
                     Phmmer by default). It then goes on to run molecular replacement (MR) on the \
                     best of these. Post MR it refine each resulting solution and do some initial \
                     model building to assess the likely success or otherwise of the molecular \
                     replacement search.' )

    jobDirectory = jobInfo['fileroot']

    results_txt = os.path.join(jobDirectory, "search_mrbump_1", "results", "results.txt")
    rows = []
    if os.path.isfile(results_txt):
        with open(results_txt) as f:
            rows = final_solutions(f.read())
    if rows:
        finalFold = results.addFold(label='Final solution', initiallyOpen=True)
        finalFold.append('The best placement after Phaser has also tried the alternative '
                         'space groups, refined with Refmac. The tables further down give '
                         'every search model and every stage.')
        table = finalFold.addTable()
        table.addData(title='Search model', data=[r['model'] for r in rows])
        table.addData(title='Space group', data=[r['sg'] for r in rows])
        table.addData(title='TFZ', data=[r['tfz'] for r in rows])
        table.addData(title='LLG', data=[r['llg'] for r in rows])
        table.addData(title='R', data=[r['r'] for r in rows])
        table.addData(title='R-free', data=[r['rfree'] for r in rows])

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
                      The resulting MR solution is then refined with Refmac (Final R, Final R-free)<br/>')

    if not os.path.isfile(os.path.join(jobDirectory, "search_mrbump_1", "results", "results.txt")):
        results.append("Molecular replacement results will appear here soon...")
    else: 
        alog=open(os.path.join(jobDirectory, "search_mrbump_1", "results", "results.txt"), "r")
        lines="".join(alog.readlines())
        #lines=alog.readlines()
        alog.close()

        tableFoldmr.append("<pre>%s</pre>" % lines)

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
