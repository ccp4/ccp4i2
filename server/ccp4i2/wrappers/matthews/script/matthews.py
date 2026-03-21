import shutil
from lxml import etree
from ccp4i2.core.CCP4PluginScript import CPluginScript


class matthews(CPluginScript):

    TASKNAME = 'matthews'
    TASKVERSION = 0.1
    WHATNEXT = ['ProvideAsuContents']

    ERROR_CODES = {
        201: {'description': 'No molecular weight could be determined from inputs'},
        202: {'description': 'Failed to load reflection file'},
        203: {'description': 'Matthews coefficient calculation failed'},
    }

    def startProcess(self):
        # Load the reflection file to get cell parameters and space group
        try:
            self.container.inputData.HKLIN.loadFile()
        except Exception as e:
            self.appendErrorReport(202, str(e))
            return CPluginScript.FAILED

        fileContent = self.container.inputData.HKLIN.fileContent

        # Determine molecular weight from the available input mode
        mode = str(self.container.inputData.MODE) if self.container.inputData.MODE.isSet() else 'asu_components'
        molWt = None
        nRes = None
        polymerMode = ""

        if mode == 'asu_components' and self.container.inputData.ASUIN.isSet():
            # Calculate total weight from ASU content sequences
            try:
                self.container.inputData.ASUIN.loadFile()
                seqList = self.container.inputData.ASUIN.fileContent.seqList
                totWeight = 0.0
                for seqObj in seqList:
                    if seqObj.nCopies > 0:
                        if seqObj.polymerType == "PROTEIN":
                            if polymerMode == "D":
                                polymerMode = "C"
                            elif polymerMode == "":
                                polymerMode = "P"
                        if seqObj.polymerType in ["DNA", "RNA"]:
                            if polymerMode == "P":
                                polymerMode = "C"
                            elif polymerMode == "":
                                polymerMode = "D"
                    totWeight += seqObj.molecularWeight(seqObj.polymerType)
                molWt = totWeight
            except Exception:
                pass
        elif mode == 'nres' and self.container.inputData.NRES.isSet():
            nRes = int(self.container.inputData.NRES)
        elif mode == 'molwt' and self.container.inputData.MOLWT.isSet():
            molWt = float(self.container.inputData.MOLWT)

        if molWt is None and nRes is None:
            self.appendErrorReport(201, f'Mode={mode}')
            return CPluginScript.FAILED

        # Run matthews_coef via the MTZ fileContent method
        try:
            if nRes is not None:
                rv = fileContent.matthewsCoeff(nRes=nRes, polymerMode=polymerMode)
            else:
                rv = fileContent.matthewsCoeff(molWt=molWt, polymerMode=polymerMode)
        except Exception as e:
            self.appendErrorReport(203, str(e))
            return CPluginScript.FAILED

        # Build XML output with results table
        xmlroot = etree.Element('MatthewsAnalysis')

        vol = rv.get('cell_volume', None)
        if vol is not None:
            volTag = etree.SubElement(xmlroot, 'cellVolume')
            volTag.text = '{0:.1f}'.format(float(vol))

        results = rv.get('results', [])
        if results:
            compositions = etree.SubElement(xmlroot, 'matthewsCompositions')
            best_prob = -1
            best_nmol = 1
            for result in results:
                comp = etree.SubElement(compositions, 'composition')
                for tag, key, fmt in [
                    ('nMolecules', 'nmol_in_asu', '{0}'),
                    ('solventPercentage', 'percent_solvent', '{0:.2f}'),
                    ('matthewsCoeff', 'matth_coef', '{0:.2f}'),
                    ('matthewsProbability', 'prob_matth', '{0:.2f}'),
                ]:
                    el = etree.SubElement(comp, tag)
                    el.text = fmt.format(result.get(key, 0))

                prob = float(result.get('prob_matth', 0))
                if prob > best_prob:
                    best_prob = prob
                    best_nmol = result.get('nmol_in_asu', 1)

            # Summary element
            summary = etree.SubElement(xmlroot, 'summary')
            nmolTag = etree.SubElement(summary, 'mostLikelyNmol')
            nmolTag.text = str(best_nmol)
            best_result = next(
                (r for r in results if r.get('nmol_in_asu') == best_nmol), results[0]
            )
            solvTag = etree.SubElement(summary, 'solventContent')
            solvTag.text = '{0:.1f}'.format(float(best_result.get('percent_solvent', 0)))

        # Write program XML
        newXml = etree.tostring(xmlroot, pretty_print=True)
        with open(self.makeFileName('PROGRAMXML') + '_tmp', 'w') as f:
            f.write(newXml.decode('utf-8'))
        shutil.move(self.makeFileName('PROGRAMXML') + '_tmp', self.makeFileName('PROGRAMXML'))

        return CPluginScript.SUCCEEDED
