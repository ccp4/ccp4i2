"""Multi-domain NCS averaging via `dm`.

Improves phases by NCS averaging where DIFFERENT domains follow DIFFERENT NCS
transformations -- the case `parrot` (whole-monomer NCS) cannot handle. The
user supplies a model and per-domain residue ranges; the wrapper derives the
per-domain operators (by superposition) and averaging masks (gemmi), then runs
`dm`. Masks and operators are pipeline-internal -- not modelled as CData yet.

Operator convention (validated against the AHIR driver case): operators map the
reference/masked copy onto each NCS copy (x' = R x + t, identity first).
"""
import os
import re

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript

from . import dm_ncs_lib


class dm_multidomain(CPluginScript):
    TASKNAME = 'dm_multidomain'
    TASKCOMMAND = 'dm'
    PERFORMANCECLASS = 'CExpPhasPerformance'

    def validity(self):
        # CCP4-free: only inspects container parameters.
        error = super(dm_multidomain, self).validity()
        ctrl = self.container.controlParameters
        phase_source = (str(ctrl.PHASE_SOURCE)
                        if ctrl.PHASE_SOURCE.isSet() else 'input')
        if phase_source == 'input' and \
                not self.container.inputData.ABCD.isSet():
            error.append(
                klass=self.TASKNAME, code=201,
                details='Starting phases (ABCD) are required unless phases '
                        'are calculated from the model',
                name=f'{self.TASKNAME}.container.inputData.ABCD',
                severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def processInputFiles(self):
        import gemmi

        inp = self.container.inputData
        ctrl = self.container.controlParameters

        # 1. build the dm input MTZ + LABIN, from either supplied phases (ABCD)
        #    or phases CALCULATED from the model with servalcat sigmaa (bulk
        #    solvent + sigmaA weighting). dm consumes phase + FOM, never map
        #    coefficients, so the sigmaA weighting is delivered via FOMO.
        phase_source = (str(ctrl.PHASE_SOURCE)
                        if ctrl.PHASE_SOURCE.isSet() else 'input')
        try:
            if phase_source == 'model':
                self.hklin, self._labin, self._has_free = self._phases_from_model()
            else:
                cols = [['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],
                        'ABCD']
                if inp.FREERFLAG.isSet():
                    cols.append('FREERFLAG')
                self.hklin, _, error = self.makeHklin0(cols)
                if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                    return CPluginScript.FAILED
                inp.ABCD.setContentFlag()
                if inp.ABCD.contentFlag == \
                        CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
                    phase = "HLA=ABCD_HLA HLB=ABCD_HLB HLC=ABCD_HLC HLD=ABCD_HLD"
                else:
                    phase = "PHIO=ABCD_PHI FOMO=ABCD_FOM"
                self._labin = f"FP=F_SIGF_F SIGFP=F_SIGF_SIGF {phase}"
                self._has_free = bool(inp.FREERFLAG.isSet())
                if self._has_free:
                    self._labin += " FREE=FREERFLAG_FREER"
        except Exception as exc:
            import traceback
            traceback.print_exc()
            self.appendErrorReport(203, f'Phase preparation failed: {exc}')
            return CPluginScript.FAILED

        # 2. derive per-domain operators + averaging masks from the model.
        #    Working state is stored under _-prefixed names so the base object's
        #    "smart" __setattr__ (which turns dict/list into CData) is bypassed.
        try:
            structure = gemmi.read_structure(str(inp.XYZIN.fullPath))
            structure.setup_entities()
            model = structure[0]
            cell = structure.cell
            sg = structure.find_spacegroup()

            domains = []
            for i, d in enumerate(ctrl.DOMAINS, start=1):
                bounds = d.residue_bounds()
                if bounds is None:
                    raise ValueError(f"domain {i}: firstRes/lastRes not set")
                lo, hi = bounds
                prefix = (str(d.chainId) + '_') if d.chainId.isSet() else ''
                domains.append(dict(name=f"{prefix}{lo}_{hi}", lo=lo, hi=hi,
                                    mode=d.averaging_mode()))
            ref, copies = dm_ncs_lib.detect_reference_and_copies(
                model,
                reference=str(ctrl.REFERENCE_CHAIN) if ctrl.REFERENCE_CHAIN.isSet() else None,
                copies=str(ctrl.COPY_CHAINS).split(',') if ctrl.COPY_CHAINS.isSet() else None)
            print(f"dm_multidomain: reference {ref}, copies {copies}")

            radius = float(ctrl.MASK_RADIUS) if ctrl.MASK_RADIUS.isSet() else 2.5
            operators_by_domain = {}
            ncsin = []   # ordered [(domain_name, mask_path)] for non-excluded
            for d in domains:
                if d['mode'] == 'exclude':
                    continue
                ops, rmsds = dm_ncs_lib.domain_operators(
                    model, ref, copies, d['lo'], d['hi'])
                operators_by_domain[d['name']] = ops
                rmsd_str = ", ".join(f"{c}:{r:.2f}" for c, r in rmsds.items())
                print(f"  domain {d['name']} ({d['lo']}-{d['hi']}, {d['mode']}): "
                      f"RMSD A {rmsd_str}")
                mask = os.path.join(self.workDirectory, f"mask_{d['name']}.msk")
                nset = dm_ncs_lib.write_domain_mask(
                    model, ref, d['lo'], d['hi'], cell, sg, mask, radius=radius)
                print(f"    mask {os.path.basename(mask)}: {nset} points")
                ncsin.append((d['name'], mask))

            if not ncsin:
                self.appendErrorReport(201, 'No domains to average (all excluded)')
                return CPluginScript.FAILED

            # solvent content: explicit override or Matthews estimate
            if ctrl.SOLVENT_CONTENT.isSet():
                solc = float(ctrl.SOLVENT_CONTENT)
            else:
                solc = dm_ncs_lib.estimate_solvent_fraction(structure) or 0.5
                print(f"  estimated solvent fraction: {solc}")

            self._domains = domains
            self._operators_by_domain = operators_by_domain
            self._ncsin = ncsin
            self._solc = solc
        except Exception as exc:
            import traceback
            traceback.print_exc()
            self.appendErrorReport(202, f'NCS preparation failed: {exc}')
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def _phases_from_model(self):
        """Calculate sigmaA-weighted, bulk-solvent-corrected starting phases
        from XYZIN with `servalcat sigmaa`. Returns (hklin, labin, has_free).

        dm gets the OBSERVED amplitudes (FP/SIGFP) plus the sigmaA map phase
        (PHWT, == model phase) and the sigmaA FOM -- so the weighting enters
        through FOMO, exactly the channel dm uses.
        """
        import subprocess
        inp = self.container.inputData
        f_mtz, _, error = self.makeHklin0(
            [['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            raise RuntimeError('could not prepare F/SIGF for servalcat')
        prefix = os.path.join(self.workDirectory, 'sigmaa')
        cmd = ['servalcat', 'sigmaa', '--hklin', f_mtz,
               '--labin', 'F_SIGF_F,F_SIGF_SIGF',
               '--model', str(inp.XYZIN.fullPath),
               '--source', 'xray', '-o', prefix]
        print('dm_multidomain: ' + ' '.join(cmd))
        result = subprocess.run(cmd, cwd=self.workDirectory,
                                capture_output=True, text=True, timeout=600)
        sa_mtz = prefix + '.mtz'
        if result.returncode != 0 or not os.path.exists(sa_mtz):
            raise RuntimeError(
                f'servalcat sigmaa failed (rc={result.returncode}): '
                f'{result.stderr[-400:]}')
        # FP/SIGFP = observed; PHWT = sigmaA map phase; FOM = sigmaA weight
        return sa_mtz, 'FP=FP SIGFP=SIGFP PHIO=PHWT FOMO=FOM', False

    def makeCommandAndScript(self):
        ctrl = self.container.controlParameters
        self.hklout = os.path.join(self.workDirectory, "hklout.mtz")

        # command line: HKLIN/HKLOUT + one NCSIN<n> per averaged domain
        self.appendCommandLine(['HKLIN', self.hklin, 'HKLOUT', self.hklout])
        for i, (_, mask) in enumerate(self._ncsin, start=1):
            self.appendCommandLine([f'NCSIN{i}', mask])

        labin = self._labin
        labout = "PHIDM=PHIDM FOMDM=FOMDM FCDM=FCDM PHICDM=PHICDM"
        # cross-validate (free-R per cycle) when a free set is available
        ncross = 2 if self._has_free else 1

        ncycle = int(ctrl.NCYCLES) if ctrl.NCYCLES.isSet() else 10
        for line in dm_ncs_lib.build_keyword_script(
                self._domains, self._operators_by_domain, self._solc, ncycle,
                mode_solv=bool(ctrl.MODE_SOLVENT),
                mode_hist=bool(ctrl.MODE_HISTOGRAM),
                labin=labin, labout=labout, ncross=ncross):
            self.appendCommandScript(line)

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        # PHIDM/FOMDM are PHI/FOM phases; FCDM/PHICDM are map coefficients
        self.container.outputData.ABCDOUT.contentFlag.set(
            CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM)
        self.container.outputData.ABCDOUT.annotation = \
            self.jobNumberString() + ' Phases from multi-domain NCS averaging'
        self.container.outputData.FPHIOUT.contentFlag.set(1)
        self.container.outputData.FPHIOUT.subType.set(1)
        self.container.outputData.FPHIOUT.annotation = \
            self.jobNumberString() + ' Map coefficients from multi-domain dm'

        # build the result XML: loggraph tables + per-cycle + per-domain NCS
        # correlations, as substrate for a graphically rich report.
        rootNode = etree.Element("DmMultidomainResult")
        logtext = self._read_log()
        self._add_smartie_graphs(etree.SubElement(rootNode, 'SmartieGraphs'))
        self._add_per_cycle(rootNode, logtext)
        self._add_ncs_correlations(rootNode, logtext)
        final_fom = self._mean_final_fom(logtext)

        self.xmlout = self.makeFileName('PROGRAMXML')
        with open(self.xmlout, 'w') as xmlFile:
            CCP4Utils.writeXML(xmlFile, etree.tostring(rootNode, pretty_print=True))

        if final_fom is not None:
            self.container.outputData.PERFORMANCE.FOM = final_fom

        return self.splitHklout(
            ['FPHIOUT', 'ABCDOUT'],
            ['FCDM,PHICDM', 'PHIDM,FOMDM'])

    # -- log scraping ---------------------------------------------------------
    def _read_log(self):
        # NB: do NOT strip HTML tags -- dm embeds each $TABLE inside a multi-line
        # <param value="..."> applet tag, and a tag-strip would delete the table
        # content. The plain-text per-cycle/NCS lines parse fine from raw text.
        try:
            with open(self.makeFileName('LOG'), encoding='utf-8', errors='replace') as fh:
                return fh.read()
        except OSError:
            return ''

    def _add_smartie_graphs(self, smartieNode):
        """Native dm $TABLE/$GRAPHS (completeness, mean FOM & dphi vs resolution,
        Free-R vs cycle) -> directly renderable loggraphs."""
        from ccp4i2.smartie import smartie
        from ccp4i2.pimple.logtable import CCP4LogToEtree
        logfile = smartie.parselog(self.makeFileName('LOG'))
        for table in logfile.tables():
            if table.ngraphs() > 0:
                smartieNode.append(CCP4LogToEtree(table.rawtable()))

    def _add_per_cycle(self, rootNode, logtext):
        """Per-cycle metrics for graphing: perturbation gamma, mean combined
        FOM, and the mean NCS correlation of each masked domain. Emitted as
        child elements (Number/Gamma/FOM/Corr_<n>) so the report can plot them.
        """
        text = re.sub(r'<[^>]+>', '', logtext)   # cycle blocks are plain text
        # segment by the "Cycle    N" section headers
        marks = list(re.finditer(r'\n\s*Cycle\s+(\d+)\s*\n', text))
        if not marks:
            return
        node = etree.SubElement(rootNode, 'PerCycle')
        node.set('title', 'Per-cycle statistics')
        for k, m in enumerate(marks):
            start = m.end()
            end = marks[k + 1].start() if k + 1 < len(marks) else len(text)
            seg = text[start:end]
            row = etree.SubElement(node, 'Cycle')
            etree.SubElement(row, 'Number').text = m.group(1)
            gm = re.search(r'Overall value\s+([-\d.]+)', seg)
            if gm:
                etree.SubElement(row, 'Gamma').text = gm.group(1)
            fom = self._cycle_fom(seg)
            if fom is not None:
                etree.SubElement(row, 'FOM').text = f"{fom:.4f}"
            for dom, corr in self._cycle_domain_correlations(seg).items():
                etree.SubElement(row, f'Corr_{dom}').text = f"{corr:.4f}"

    @staticmethod
    def _cycle_fom(seg):
        """NREFLS-weighted mean of the combined FOM (col 8) of the per-cycle
        sigma-a table (cols: RMIN RMAX S^2 NREFLS SIGMAA FOMobs FOMcalc
        FOMcomb DPHI*3)."""
        num, den = 0.0, 0.0
        for line in seg.splitlines():
            f = line.split()
            if len(f) == 11:
                try:
                    nref = float(f[3])
                    num += nref * float(f[7])
                    den += nref
                except ValueError:
                    pass
        return num / den if den else None

    @staticmethod
    def _cycle_domain_correlations(seg):
        """Mean off-diagonal of each domain's per-cycle NCS correlation matrix
        (the last matrix printed for each domain in the cycle)."""
        out = {}
        for m in re.finditer(
                r'CORRELATIONS BETWEEN REGIONS IN DOMAIN\s+(\d+)(.*?)'
                r'(?=CORRELATIONS BETWEEN REGIONS|\Z)', seg, re.S):
            rows = []
            for line in m.group(2).splitlines():
                vals = re.findall(r'[-\d]+\.\d+', line)
                if len(vals) >= 3 and all(
                        abs(float(v)) <= 1.5 for v in vals):
                    rows.append([float(v) for v in vals])
                elif rows:
                    break
            n = len(rows)
            if n >= 2 and all(len(r) == n for r in rows):
                off = [rows[i][j] for i in range(n) for j in range(n) if i != j]
                if off:
                    out[m.group(1)] = sum(off) / len(off)
        return out

    def _add_ncs_correlations(self, rootNode, logtext):
        """Per-domain NCS averaging correlation, initial vs final (+ dm's
        OK/WARNING verdict) -- the headline 'did averaging work' signal."""
        m = re.search(r'NCS correlations between related density regions:(.*?)'
                      r'(?:Refined NCS matrices|\Z)', logtext, re.S)
        if not m:
            return
        node = etree.SubElement(rootNode, 'NCSCorrelations')
        for line in m.group(1).splitlines():
            mm = re.match(r'\s*(\d+)\s+([-\d.]+)\s+([-\d.]+)\s*(\S.*)?$', line)
            if mm:
                d = etree.SubElement(node, 'Domain')
                d.set('number', mm.group(1))
                d.set('initial', mm.group(2))
                d.set('final', mm.group(3))
                d.set('status', (mm.group(4) or 'OK').strip())

    def _mean_final_fom(self, logtext):
        """Mean final FOM (FOMdm column of the 'Phase and weight statistics'
        table) for the performance indicator."""
        block = logtext.split('Phase and weight statistics', 1)
        if len(block) < 2:
            return None
        seg = block[1].split('$TABLE', 1)[0]   # bound to this table
        # the data block is the $$-segment with the most 4-column rows
        best = []
        for part in seg.split('$$'):
            rows = re.findall(r'^\s*[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)\s*$',
                              part, re.M)
            if len(rows) > len(best):
                best = rows
        foms = [float(x) for x in best]
        return round(sum(foms) / len(foms), 4) if foms else None
