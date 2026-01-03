import json
import os

from ccp4i2.report import Report


class slicendice_report(Report):
    TASKNAME = 'slicendice'
    USEPROGRAMXML = True
    RUNNING = True
    CSS_VERSION = '0.1.0'

    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,cssVersion=self.CSS_VERSION,**kw)
        
        jobDirectory = jobInfo['fileroot']
        jsdat = None
        with open(os.path.join(jobDirectory, "slicendice_0", "slicendice_results.json")) as json_file:
                jsdat = json.load(json_file)
        sbid = str(self.xmlnode.findall('.//RunInfo/Best/bid')[0].text)
        bestR = float(self.xmlnode.findall('.//RunInfo/Best/R')[0].text)
        bestRFree = float(self.xmlnode.findall('.//RunInfo/Best/RFree')[0].text)
        sbestR = f'{round(bestR, 3):.3g}'
        sbestRFree = f'{round(bestRFree, 3):.3g}'
        results = self.addResults()
        results.append('SliceNDice is an automated pipeline for processing models for use in Molecular Replacement (MR). It can correct B-factor \
                columns in predicted models as well "slicing" search models into rigid regions that may achieve better placement in MR. \
                Following this, it will attempt MR using all models created.')
    
        # New stuff
        results.append('The best MR solution found by Slice\' n\' Dice using %s splits, resulted in an R-Factor of %s and an R-Free \
                        of %s after refinement (all solutions are shown in the table below).'%(sbid,sbestR,sbestRFree))
        # Table
        table = results.addTable(select=".//RunInfo/Sol", style="float:left;")
        table_contents = [["No. of Splits", "SolID"], ["Log-liklihood gain", "llg"], ["tfz score", "tfz"],
                          ["R-Factor", "srf"], ["R-Free", "sre"]]
        for title, solres in table_contents:
            table.addData(title=title, select=solres)

        results.append('The output data from the Refmac refinement, for the best solution, are given below in the Output Data Section')

        tableFoldmr = results.addFold(label='Results for molecular replacement and refinement', initiallyOpen=False)
        tableFoldmr.append('Each model is used in Phaser to do molecular replacement. \
                          The resulting MR solution is then refined with Refmac (Final R, Final R)<br/>')


        # This is the log file
        if not os.path.isfile(os.path.join(jobDirectory, "slicendice_0", "slicendice.log")):
            results.append("Molecular replacement results will appear here soon...")
        else: 
            alog=open(os.path.join(jobDirectory, "slicendice_0", "slicendice.log"), "r")
            lines="".join(alog.readlines())
            #lines=alog.readlines()
            alog.close()
    
            tableFoldmr.append("<pre>%s</pre>" % lines)
    
    
        self.addTaskReferences()
