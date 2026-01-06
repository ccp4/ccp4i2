import os
import re

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript


class molrep_mr(CPluginScript):

    TASKTITLE='Molecular replacement using Molrep'
    TASKNAME = 'molrep_mr'
    TASKCOMMAND = 'molrep'
    TASKMODULE = 'wrappers'
    MAINTAINER = 'andrey.lebedev@stfc.ac.uk'

    ERROR_CODES = {
        101: {'severity': SEVERITY_WARNING, 'description': 'Failed extracting tables from molrep.doc'},
        102: {'severity': SEVERITY_WARNING, 'description': 'Failed writing tables to program.xml'},
        103: {'severity': SEVERITY_WARNING, 'description': 'No tables extracted from molrep.doc'},
        104: {'description': 'No output coordinate file from Molrep'},
        105: {'description': 'No output log file from Molrep'},
        106: {'description': 'Error parsing log file from Molrep'},
        107: {'severity': SEVERITY_WARNING, 'description': 'molrep.doc file not found - molrep may have failed'},
        108: {'severity': SEVERITY_WARNING, 'description': 'Failed to convert input observations to F_SIGF'},
        109: {'severity': SEVERITY_WARNING, 'description': 'Failed to process input coordinate file'},
        110: {'severity': SEVERITY_WARNING, 'description': 'Failed to create scratch directory'},
        111: {'severity': SEVERITY_WARNING, 'description': 'Failed to prepare model for Molrep'},
        112: {'severity': SEVERITY_WARNING, 'description': 'Failed to write sequence file'},
        113: {'severity': SEVERITY_WARNING, 'description': 'Failed to process output coordinate file'},
        114: {'severity': SEVERITY_WARNING, 'description': 'Failed to read molrep.doc file'},
        115: {'severity': SEVERITY_WARNING, 'description': 'Failed to save program XML'},
        116: {'severity': SEVERITY_WARNING, 'description': 'Failed to extract Laue data from log'},
    }

    def processInputFiles(self):
        # Ensure the obs data is in form of F_SIGF
        if self.container.guiParameters.PERFORM != 'den':
            try:
                from ccp4i2.core import CCP4XtalData

                # Using CObsDataFile.convert() did not work for input of anomalous Is (as from aimless)
                self.F_SIGF_hklin, errReport = self.makeHklin([['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
                if self.F_SIGF_hklin is None:
                    self.appendErrorReport(108, 'makeHklin returned None')
                    return CPluginScript.FAILED
            except Exception as e:
                self.appendErrorReport(108, str(e))
                return CPluginScript.FAILED

        if self.container.guiParameters.PERFORM != 'srf':
            # Check if XYZIN is actually set before processing
            if self.container.inputData.XYZIN.isSet():
                try:
                    if self.container.inputData.XYZIN.isSelectionSet():
                        self.selectedXYZIN = os.path.join(self.workDirectory, 'XYZIN_selected_atoms.pdb')
                        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedXYZIN)
                    elif self.container.inputData.XYZIN.getExt() != '.pdb':
                        self.selectedXYZIN = os.path.join(self.workDirectory, 'XYZIN_pdb.pdb')
                        self.container.inputData.XYZIN.convertFormat('pdb', self.selectedXYZIN)
                    else:
                        self.selectedXYZIN = str(self.container.inputData.XYZIN)
                except Exception as e:
                    self.appendErrorReport(109, str(e))
                    return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.controlParameters
        gui = self.container.guiParameters

        self.path_wrk = str(self.getWorkDirectory())
        self.path_scr = os.path.join(self.path_wrk, 'scratch')

        try:
            if not os.path.exists(self.path_scr):
                os.mkdir(self.path_scr)
        except Exception as e:
            self.appendErrorReport(110, str(e))
            return CPluginScript.FAILED

        if str(gui.PERFORM) != 'srf' and inp.XYZIN.isSet():
            try:
                pdbin = os.path.join(self.path_wrk, 'molrep_input.pdb')
                inp.XYZIN.replaceMSE(pdbin)
                self.appendCommandLine(['-m', pdbin])
            except Exception as e:
                self.appendErrorReport(111, str(e))
                return CPluginScript.FAILED

        if inp.ASUIN.isSet() and par.SEQ.__str__() != 'n':
            try:
                seqin = os.path.join(self.workDirectory, 'SEQIN.fasta')
                inp.ASUIN.writeFasta(seqin)
                self.appendCommandLine(['-s', seqin])
            except Exception as e:
                self.appendErrorReport(112, str(e))

        if str(gui.PERFORM) == 'pat':
            if inp.XYZIN_FIX.isSet():
                self.appendCommandLine(['-mx', str(inp.XYZIN_FIX.fullPath)])

        if str(gui.PERFORM) == 'den':
            self.appendCommandLine(['-mx', str(inp.XYZIN_FIX.fullPath)])

        self.appendCommandLine(['-ps', str(self.path_scr) + '/'])
        self.appendCommandLine(['-po', str(self.path_wrk) + '/'])
        self.appendCommandLine(['-i'])

        if str(gui.PERFORM) == 'den':
            self.appendCommandLine(['-f', str(inp.F_PHI_MAP.fullPath)])
            self.appendCommandScript("labin F=F PH=PHI")
            self.appendCommandScript("diff m")
        else:
            self.appendCommandLine(['-f', self.F_SIGF_hklin])
            self.appendCommandScript("labin F=F SIGF=SIGF")

        if str(gui.PERFORM) == 'den':
            if str(par.PRF) == 'n':
                pass
            elif str(par.PRF) == 'y':
                self.appendCommandScript("prf y")
            elif str(par.PRF) == 's':
                self.appendCommandScript("prf s")

        if str(par.NP) != 'Auto':
            self.appendCommandScript("np %s" % str(par.NP))

        if str(par.NMON) != 'Auto':
            self.appendCommandScript("nmon %s" % str(par.NMON))

        if par.SCORE.isSet():
            if str(par.SCORE) == 'y':
                pass
            elif str(par.SCORE) == 'n':
                self.appendCommandScript("score n")
            elif str(par.SCORE) == 'c':
                self.appendCommandScript("score c")

        if par.ANISO.isSet():
            if str(par.ANISO) == 'y':
                pass
            elif str(par.ANISO) == 'n':
                self.appendCommandScript("aniso n")
            elif str(par.ANISO) == 'k':
                self.appendCommandScript("aniso k")

        if par.SEQ.isSet():
            if str(par.SEQ) == 'y':
                pass
            elif str(par.SEQ) == 'd':
                self.appendCommandScript("seq d")
            elif str(par.SEQ) == 'n':
                self.appendCommandScript("seq n")

        if par.SURF.isSet():
            if str(par.SURF) == 'y':
                pass
            elif str(par.SURF) == 'c':
                self.appendCommandScript("surf c")
            elif str(par.SURF) == 'n':
                self.appendCommandScript("surf n")
            elif str(par.SURF) == '2':
                self.appendCommandScript("surf 2")
            elif str(par.SURF) == 'a':
                self.appendCommandScript("surf a")

        if par.NMON_EXP.isSet():
            if str(par.NMON_EXP) != 'Auto':
                self.appendCommandScript("nmon_exp %s" % str(par.NMON_EXP))

        if par.RESMIN.isSet():
            self.appendCommandScript("resmin %s" % str(par.RESMIN))
        if par.RESMAX.isSet():
            self.appendCommandScript("resmax %s" % str(par.RESMAX))

        if par.SG_OPTIONS == 'specify':
            if par.SG.isSet():
                self.appendCommandScript("sg %s" % str(par.SG))
        elif par.SG_OPTIONS == 'laue':
            self.appendCommandScript("sg all")

        return 0

    def processOutputFiles(self):
        import os
        out = self.container.outputData
        gui = self.container.guiParameters
        par = self.container.controlParameters

        # Rename output file or report failure
        if par.SG_OPTIONS != 'laue' and gui.PERFORM != 'srf':
            fileName = os.path.join(self.path_wrk, 'molrep.pdb')
            if os.path.isfile(fileName):
                try:
                    with open(fileName) as istream:
                        content = istream.read()
                    content = re.sub(r'\n#MOLECULE\s+[0-9]+\s*', '\n', content)
                    with open(str(out.XYZOUT.fullPath), 'w') as ostream:
                        ostream.write(content)
                except Exception as e:
                    self.appendErrorReport(113, str(e))
                    return CPluginScript.FAILED
            else:
                self.appendErrorReport(104, str(out.XYZOUT), stack=False)
                return CPluginScript.FAILED

            if gui.PERFORM == 'pat':
                out.XYZOUT.annotation = 'Model from molrep MR'
            elif gui.PERFORM == 'den':
                out.XYZOUT.annotation = 'Model from molrep density search'

        docout = os.path.join(self.path_wrk, 'molrep.doc')
        xmlout = str(self.makeFileName('PROGRAMXML'))

        if os.path.exists(docout):
            try:
                self.saveProgramXml(docout, xmlout)
                os.rename(docout, docout + '.txt')
            except Exception as e:
                self.appendErrorReport(115, str(e))
        else:
            self.appendErrorReport(107, docout, stack=False)

        if gui.PERFORM == 'srf':
            psfile = os.path.join(self.getWorkDirectory(), 'molrep_rf.ps')
            if os.path.exists(psfile):
                out.PSOUT.set(psfile)
                out.PSOUT.annotation.set('Self-rotation function')

        return CPluginScript.SUCCEEDED

    def saveProgramXml(self, docFileName, programXmlFileName):
        from ccp4i2.core import CCP4Utils
        from lxml import etree

        titles = []
        status = 0
        rf = None

        results = etree.Element('MolrepResult')
        tf = etree.Element('MR_TF')
        results.append(tf)

        for key, value in [['err_level', '0'],
                           ['err_message', 'normal termination'],
                           ['n_solution', '1'],
                           ['mr_score', '0.0000']]:
            e = etree.Element(key)
            e.text = value
            tf.append(e)

        try:
            docfileText = CCP4Utils.readFile(docFileName)
            docfileList = docfileText.split('\n')
        except Exception as e:
            self.appendErrorReport(114, str(e))
            return

        for line in docfileList:
            if status == 0:
                lstr = line.strip()
                lst1 = "--- Translation function ---"
                lst2 = "--- phased translation function ---"
                if lstr == lst1 or lstr == lst2:
                    status = 1
            elif status == 1:
                if line.strip().startswith('RF '):
                    titles = line.replace("(", " ").replace(")", "").replace("/", "_").split()
                    rf = etree.Element('RFpeaks')
                    results.append(rf)
                else:
                    words = line.replace("(", " ").replace(")", "").replace("-", " -").split()
                    if len(words) == len(titles):
                        try:
                            for i in (0, 1):
                                int(words[i])
                            for i in range(2, len(words)):
                                float(words[i])
                            peak = etree.Element('RFpeak')
                            for key, value in zip(titles, words):
                                e = etree.Element(key)
                                e.text = str(float(value))
                                peak.append(e)
                            if rf is not None:
                                rf.append(peak)
                        except (ValueError, IndexError):
                            pass

        if self.container.controlParameters.SG_OPTIONS == 'laue':
            data = self.extractLaueDataFromLog()
            if data:
                try:
                    scores = ['%6.4f' % float(e1.findall('RFpeak')[-1].findtext('for'))
                              for e1 in results.findall('RFpeaks')]
                    if len(scores) != len(data):
                        scores = None
                except (ValueError, IndexError, AttributeError):
                    scores = None

                eleNames = ['space_group', 'score', 'contrast']
                eLaue = etree.Element('laue_group_alternatives')
                for d in data:
                    eTest = etree.Element('test')
                    eLaue.append(eTest)
                    for ii in range(3):
                        e = etree.Element(eleNames[ii])
                        e.text = str(d[ii])
                        eTest.append(e)
                    if scores:
                        eTest.find('score').text = scores.pop(0)
                results.append(eLaue)

        CCP4Utils.saveEtreeToFile(results, programXmlFileName)

    def extractLaueDataFromLog(self):
        from ccp4i2.core import CCP4Utils
        try:
            text = CCP4Utils.readFile(self.makeFileName('LOG'))
        except Exception as e:
            self.appendErrorReport(105, f"{self.makeFileName('LOG')}: {str(e)}")
            return []

        m = re.search(r'(.*)Space Group Checking(.*)', text, re.DOTALL)
        if m is None:
            self.appendErrorReport(116, 'Space Group Checking section not found in log')
            return []

        lines = m.groups()[1].split('\n')[4:]
        data = []
        for line in lines:
            if line.find('+-') >= 0:
                break
            try:
                spg = line[14:30].strip()
                score = line[44:49].strip()
                cntr = line[50:56].strip()
                if spg:
                    data.append([spg, score, cntr])
            except IndexError:
                continue
        return data
