import json
import os
import shutil
from pathlib import Path

from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.lib.utils.logs.mrparse_log import search_program_failure


class mrparse(CPluginScript):
    TASKNAME = 'mrparse'
    TASKCOMMAND = 'mrparse'
    WHATNEXT = ['phaser_simple_phil', 'phaser_pipeline_phil', 'molrep_pipe']
    PERFORMANCECLASS = 'CExpPhasPerformance'

    ERROR_CODES = {
        201: {'description': 'MrParse could not run a search program'},
        202: {'description': 'MrParse found no homologues or models',
              'severity': CCP4ErrorHandling.SEVERITY_WARNING},
        203: {'description': 'A MrParse search program failed; results are incomplete',
              'severity': CCP4ErrorHandling.SEVERITY_WARNING},
    }

    def __init__(self, *args, **kwargs):
        self.seqin = None
        self.hklin = None
        CPluginScript.__init__(self, *args, **kwargs)

    def processInputFiles(self):
        self.seqin = self.container.inputData.SEQIN
        if self.container.inputData.F_SIGF.isSet():
            self.hklin, error = self.makeHklin([['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
            if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        pdb_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "homologs.json")
        af_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "af_models.json")
        esm_json = os.path.join(self.getWorkDirectory(), "mrparse_0", "esm_models.json")

        # Each hit's model, best first. A hit without a model file is skipped:
        # the loop used to carry the previous hit's paths over, registering a
        # duplicate under the wrong name (or NameError on the first hit).
        def register(json_path, label, key):
            if not os.path.exists(json_path):
                return
            with open(json_path, 'r') as f:
                hits = sorted(json.load(f), key=lambda k: k[key], reverse=True)
            for hit in hits:
                if not hit.get('pdb_file'):
                    continue
                xyz_in = os.path.join(self.getWorkDirectory(), "mrparse_0", hit['pdb_file'])
                if not os.path.isfile(xyz_in):
                    continue
                xyz_out = os.path.join(self.getWorkDirectory(), os.path.basename(hit['pdb_file']))
                shutil.copy(xyz_in, xyz_out)
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(xyz_out)
                self.container.outputData.XYZOUT[-1].annotation = "{} hit: {}".format(label, hit['name'])

        register(pdb_json, "PDB", 'ellg' if self.hklin else 'seq_ident')
        register(af_json, "AFDB", 'seq_ident')
        register(esm_json, "ESM", 'seq_ident')
        complexes = self._register_complexes(pdb_json)
        self._write_program_xml(complexes)

        # MrParse logs a failed search program and then writes an empty report,
        # so an unrunnable binary looks exactly like an honest "nothing found".
        # Say which of the two happened rather than reporting a bare success.
        failure = search_program_failure(
            Path(self.getWorkDirectory()) / "mrparse_0" / "mrparse.log"
        )
        found_any = len(self.container.outputData.XYZOUT) > 0
        if failure and not found_any:
            self.appendErrorReport(201, failure)
            return CPluginScript.FAILED
        if failure:
            self.appendErrorReport(203, failure)
        elif not found_any:
            self.appendErrorReport(
                202, 'The searches ran but matched no homologues or predicted models'
            )

        return CPluginScript.SUCCEEDED

    def _targets(self):
        """[(name, sequence)] searched for, as MrParse merged SEQIN: one is
        one kind of molecule, several a complex."""
        from ccp4i2.lib.utils.formats.mrparse_complexes import read_targets
        try:
            return read_targets(str(self.container.inputData.SEQIN.getFullPath()))
        except OSError:
            return []

    def _register_complexes(self, pdb_json):
        """A PDB entry whose hits match more than one target is a template for
        the complex (docs/multi-component-mr.md, route A): its matching chains
        are written as one more model, after the single-chain hits, to be
        placed as one rigid body. CDK4/cyclin D1 was solved with such a file
        and nothing here offered it. Returns [(entry record, XYZOUT index)]."""
        from ccp4i2.lib.utils.formats.mrparse_complexes import (
            describe, find_complexes, write_complex)
        targets = self._targets()
        if len(targets) < 2 or not os.path.exists(pdb_json):
            return []
        with open(pdb_json) as stream:
            hits = json.load(stream)
        cut_dir = os.path.join(self.getWorkDirectory(), "mrparse_0", "homologs")
        registered = []
        for entry in find_complexes(hits, targets):
            path = write_complex(entry, cut_dir, self.getWorkDirectory())
            if path is None:
                continue
            out = self.container.outputData.XYZOUT
            out.append(out.makeItem())
            out[-1].setFullPath(path)
            out[-1].annotation = "PDB complex template: " + describe(entry)
            registered.append((entry, len(out) - 1))
        return registered

    def _write_program_xml(self, complexes):
        """program.xml: how many sequences were searched for, and the complex
        templates found. The judgement routes on both: one target is one kind
        of molecule (phaser_simple_phil); several, with a template, is that
        template as one search model; several without is one ensemble per
        component (phaser_pipeline_phil). Haiku, given a CDK4/cyclin D1 FASTA,
        took the one-model route and fixed the second component."""
        from lxml import etree
        targets = self._targets()
        root = etree.Element("MrParse")
        node = etree.SubElement(root, "Targets")
        node.set("count", str(max(len(targets), 1)))
        for name, _seq in targets:
            etree.SubElement(node, "Target").text = name
        node = etree.SubElement(root, "Complexes")
        node.set("count", str(len(complexes)))
        for entry, index in complexes:
            element = etree.SubElement(node, "Complex")
            element.set("entry", entry["entry"])
            element.set("index", str(index))
            element.set("file", os.path.basename(str(self.container.outputData.XYZOUT[index])))
            for chain in entry["chains"]:
                component = etree.SubElement(element, "Component")
                component.set("chain", str(chain["chain"]))
                component.set("target", str(chain["target"]))
                component.set("identity", f"{chain['identity']:.2f}")
                component.set("hit", str(chain["hit"]))
        with open(self.makeFileName("PROGRAMXML"), "wb") as stream:
            stream.write(etree.tostring(root, pretty_print=True))

    def makeCommandAndScript(self):
        self.appendCommandLine("--seqin")
        self.appendCommandLine(str(self.seqin))
        if self.hklin:
            self.appendCommandLine("--hklin")
            self.appendCommandLine(str(self.hklin))
        if self.container.options.MAXHITS:
            self.appendCommandLine("--max_hits")
            self.appendCommandLine(str(self.container.options.MAXHITS))
        if self.container.options.DATABASE:
            self.appendCommandLine("--database")
            self.appendCommandLine((str(self.container.options.DATABASE)).lower())
        if bool(self.container.options.USEAPI):  # a CBoolean never equals the string 'True'
            self.appendCommandLine("--use_api")
        if self.container.options.PDBLOCAL.isSet():
            self.appendCommandLine("--pdb_local")
            self.appendCommandLine(str(self.container.options.PDBLOCAL))
        if self.container.options.PDBSEQDB.isSet():
            self.appendCommandLine("--pdb_seqdb")
            self.appendCommandLine(str(self.container.options.PDBSEQDB))
        if self.container.options.AFDBSEQDB.isSet():
            self.appendCommandLine("--afdb_seqdb")
            self.appendCommandLine(str(self.container.options.AFDBSEQDB))
        if self.container.options.NPROC:
            self.appendCommandLine("--nproc")
            self.appendCommandLine(str(self.container.options.NPROC))
        if self.container.options.DO_CLASSIFY == 'True':
            self.appendCommandLine("--do_classify")
        self.appendCommandLine("--ccp4cloud")
        return CPluginScript.SUCCEEDED
