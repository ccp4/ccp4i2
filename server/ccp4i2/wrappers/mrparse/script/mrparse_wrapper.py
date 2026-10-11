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
        # What MrParse knew about each model, kept with it: the judgement
        # says to choose by identity and coverage (and eLLG when reflections
        # were given), and those numbers were only in homologs.json, which
        # an agent cannot see. Haiku aligned the sequences itself to guess
        # the identity (Eleanor's case, 2026-10-10).
        models = []

        def register(json_path, label, key):
            if not os.path.exists(json_path):
                return
            with open(json_path, 'r') as f:
                hits = sorted(json.load(f), key=lambda k: k.get(key) or 0, reverse=True)
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
                # The annotation is what project_jobs and the app show: say
                # the identity and eLLG there too, not only in program.xml.
                facts = []
                if isinstance(hit.get('seq_ident'), (int, float)):
                    facts.append("%.0f%% id" % (100 * hit['seq_ident']))
                if isinstance(hit.get('ellg'), (int, float)) and hit['ellg']:
                    facts.append("eLLG %.0f" % hit['ellg'])
                if isinstance(hit.get('query_start'), int) and isinstance(hit.get('query_stop'), int):
                    facts.append("%d-%d" % (hit['query_start'], hit['query_stop']))
                self.container.outputData.XYZOUT[-1].annotation = "{} hit: {}{}".format(
                    label, hit['name'], " (%s)" % ", ".join(facts) if facts else "")
                models.append((label, os.path.basename(hit['pdb_file']), hit))

        register(pdb_json, "PDB", 'ellg' if self.hklin else 'seq_ident')
        register(af_json, "AFDB", 'seq_ident')
        register(esm_json, "ESM", 'seq_ident')
        complexes = self._register_complexes(pdb_json)
        self._write_program_xml(complexes, models)

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

    def _write_program_xml(self, complexes, models=()):
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
        # One Model per XYZOUT model, in XYZOUT order, with what MrParse
        # recorded of it: identity (a fraction), the target range it covers
        # and the coverage of a single target, eLLG (reflections given),
        # MrParse's score, the RMSD it expects, pLDDT for a prediction.
        target_length = len(targets[0][1]) if len(targets) == 1 and targets[0][1] else None
        node = etree.SubElement(root, "Models")
        node.set("count", str(len(models)))
        for index, (source, filename, hit) in enumerate(models):
            element = etree.SubElement(node, "Model")
            element.set("index", str(index))
            element.set("file", filename)
            element.set("source", source)
            element.set("name", str(hit.get("name", "")))
            for key, attr, fmt in (("seq_ident", "identity", "{:.2f}"), ("ellg", "ellg", "{:.1f}"),
                                   ("score", "score", "{:.1f}"), ("rmsd", "rmsd", "{:.2f}"),
                                   ("avg_plddt", "plddt", "{:.1f}"), ("resolution", "resolution", "{:.2f}"),
                                   ("length", "length", "{:d}"), ("ncopies", "copies", "{:d}")):
                value = hit.get(key)
                if isinstance(value, (int, float)):
                    try:
                        element.set(attr, fmt.format(value if fmt != "{:d}" else int(value)))
                    except (ValueError, TypeError):
                        pass
            start, stop = hit.get("query_start"), hit.get("query_stop")
            if isinstance(start, int) and isinstance(stop, int):
                element.set("range", f"{start}-{stop}")
                if target_length:
                    element.set("coverage", f"{min(1.0, (stop - start + 1) / target_length):.2f}")
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
