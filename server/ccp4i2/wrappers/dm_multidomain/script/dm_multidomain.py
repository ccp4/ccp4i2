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

    # -- validation -----------------------------------------------------------
    #
    # The assembly and the rigid bodies are cross-referenced by role name, and
    # a role that matches nothing is the failure this task is most prone to.
    # Catching it here puts the complaint on the field instead of 200 lines
    # into a run. validity() is polled during editing, so it stays on strings;
    # anything that needs to read the model goes in runTimeValidity().

    _ASSEMBLY_PATH = 'container.controlParameters.ASSEMBLY'
    _DOMAINS_PATH = 'container.controlParameters.DOMAINS'

    def _err(self, error, code, details, path,
             severity=CCP4ErrorHandling.SEVERITY_ERROR):
        error.append(klass=self.TASKNAME, code=code, details=details,
                     name=f'{self.TASKNAME}.{path}', severity=severity)

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

        instances = self._assembly_instances()
        if instances is not None:
            self._check_assembly(error, instances)
        self._check_bodies(error, instances)
        return error

    def _assembly_instances(self):
        """The ASSEMBLY rows as instances, or None when the list is empty
        (empty means 'detect from the model', which validity cannot do)."""
        ctrl = self.container.controlParameters
        if not ctrl.ASSEMBLY.isSet() or len(ctrl.ASSEMBLY) == 0:
            return None
        return dm_ncs_lib.parse_assembly_rows([str(r) for r in ctrl.ASSEMBLY])

    def _check_assembly(self, error, instances):
        if len(instances) < 2:
            self._err(error, 210,
                      'NCS averaging needs at least two copies; the assembly '
                      'has %d. Add a row per copy, or clear the list to '
                      'detect the copies from the model.' % len(instances),
                      self._ASSEMBLY_PATH)
        seen = {}
        for i, inst in enumerate(instances):
            for role, chain in inst.items():
                if chain in seen:
                    self._err(error, 211,
                              f'Chain {chain} appears in more than one copy '
                              f'(rows {seen[chain] + 1} and {i + 1}); each '
                              f'chain belongs to exactly one copy',
                              self._ASSEMBLY_PATH)
                seen[chain] = i

    def _check_bodies(self, error, instances):
        """Parse every rigid body, and cross-check its roles against the
        assembly. `instances` is None when the assembly is auto-detected, in
        which case only the single implicit role can be verified."""
        ctrl = self.container.controlParameters
        if len(ctrl.DOMAINS) == 0:
            self._err(error, 212,
                      'Define at least one rigid body (the whole copy is a '
                      'reasonable starting point)', self._DOMAINS_PATH)
            return

        ref_roles = set(instances[0]) if instances else None
        spans = []          # (role, lo, hi, body index) for the overlap check
        n_averaged = 0
        for i, d in enumerate(ctrl.DOMAINS, start=1):
            spec = d.segments_spec()
            if not spec:
                self._err(error, 213, f'Rigid body {i} has no residue ranges',
                          self._DOMAINS_PATH)
                continue
            try:
                segments = dm_ncs_lib.parse_segments(spec)
            except ValueError as exc:
                self._err(error, 214, f'Rigid body {i}: {exc}',
                          self._DOMAINS_PATH)
                continue
            if d.averaging_mode() != 'exclude':
                n_averaged += 1
            for role, lo, hi in segments:
                if lo > hi:
                    self._err(error, 215,
                              f'Rigid body {i}: residue range {lo}-{hi} runs '
                              f'backwards', self._DOMAINS_PATH)
                if ref_roles is not None and role not in ref_roles:
                    shown = ('the unnamed entity' if role == '_'
                             else f'"{role}"')
                    self._err(error, 216,
                              f'Rigid body {i} refers to {shown}, which is '
                              f'not in the reference copy '
                              f'({", ".join(sorted(ref_roles)) or "empty"}). '
                              f'A body can only span entities the assembly '
                              f'defines.', self._DOMAINS_PATH)
                spans.append((role, lo, hi, i))

        if n_averaged == 0 and len(ctrl.DOMAINS) > 0:
            self._err(error, 217,
                      'Every rigid body is excluded, so there is nothing to '
                      'average', self._DOMAINS_PATH)

        for a in range(len(spans)):
            for b in range(a + 1, len(spans)):
                role_a, lo_a, hi_a, body_a = spans[a]
                role_b, lo_b, hi_b, body_b = spans[b]
                if body_a == body_b or role_a != role_b:
                    continue
                if lo_a <= hi_b and lo_b <= hi_a:
                    # Not fatal -- write_partitioned_masks resolves contested
                    # points by nearest atom -- but two bodies claiming the
                    # same residues is almost always a typo, and the operator
                    # for each is fitted to residues the other also owns.
                    self._err(error, 218,
                              f'Rigid bodies {body_a} and {body_b} both claim '
                              f'residues {max(lo_a, lo_b)}-{min(hi_a, hi_b)}; '
                              f'their masks will be split between them',
                              self._DOMAINS_PATH,
                              severity=CCP4ErrorHandling.SEVERITY_WARNING)

    def runTimeValidity(self):
        """Expensive checks that need the model: do the named chains and
        residue ranges exist, and does each body actually have copies to
        average against."""
        error = super(dm_multidomain, self).runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        preview = self.ncs_preview()
        if not preview.get('ok'):
            reason = preview.get('error')
            if reason:
                self._err(error, 219, reason, self._ASSEMBLY_PATH)
            return error
        for message in preview.get('messages', []):
            self._err(error, message.get('code', 220), message['text'],
                      message.get('path', self._ASSEMBLY_PATH),
                      severity=(CCP4ErrorHandling.SEVERITY_WARNING
                                if message.get('severity') == 'warning'
                                else CCP4ErrorHandling.SEVERITY_ERROR))
        return error

    # -- interactive preview --------------------------------------------------

    def ncs_preview(self):
        """Everything the task interface needs to show the user their own
        model, reachable from the client through `plugin_method`.

        Returns the chains and entities found in XYZIN, the assembly that is
        either set or would be detected, and -- the point of the exercise --
        for every rigid body the number of matched CA atoms and the
        superposition RMSD against each copy. Those two numbers are what tell
        a user whether the ranges they typed really do move as one unit,
        before they spend a run finding out.

        Pure gemmi (via dm_ncs_lib), so it runs on the slim server, and it
        never raises: every failure comes back as data for the interface to
        show.
        """
        try:
            return self._ncs_preview()
        except Exception as exc:               # pragma: no cover - guard rail
            return {'ok': False, 'error': f'{type(exc).__name__}: {exc}'}

    def _ncs_preview(self):
        import gemmi

        inp = self.container.inputData
        ctrl = self.container.controlParameters
        if not inp.XYZIN.isSet():
            return {'ok': False, 'error': 'No model set'}
        path = str(inp.XYZIN.fullPath)
        if not os.path.exists(path):
            return {'ok': False, 'error': f'Model file not found: {path}'}

        structure = gemmi.read_structure(path)
        structure.setup_entities()
        model = structure[0]

        entities = dm_ncs_lib.group_chains_by_entity(model)
        chains = []
        for name in dm_ncs_lib.protein_chains(model):
            bounds = dm_ncs_lib.residue_bounds(model, name)
            chains.append({
                'id': name,
                'first': bounds[0] if bounds else None,
                'last': bounds[1] if bounds else None,
                'nResidues': len(dm_ncs_lib.chain_residue_map(model, name)),
                'entity': next((i for i, g in enumerate(entities)
                                if name in g), None),
            })

        detected, detected_roles = dm_ncs_lib.detect_assembly(model)
        suggestion = {
            'assembly': dm_ncs_lib.format_assembly_rows(detected,
                                                        detected_roles),
            'segments': dm_ncs_lib.suggest_segments(model, detected),
        }

        rows = ([str(r) for r in ctrl.ASSEMBLY]
                if ctrl.ASSEMBLY.isSet() else [])
        rows = [r for r in rows if r.strip()]
        if rows:
            instances = dm_ncs_lib.parse_assembly_rows(rows)
            source = 'parameter'
        else:
            instances, source = detected, 'detected'

        messages = []
        if len(instances) < 2:
            messages.append({
                'code': 210, 'severity': 'error', 'path': self._ASSEMBLY_PATH,
                'text': 'NCS averaging needs at least two copies of the '
                        'assembly; %d found in the model.' % len(instances)})

        known = set(dm_ncs_lib.protein_chains(model))
        for i, inst in enumerate(instances, start=1):
            for role, chain in inst.items():
                if chain not in known:
                    messages.append({
                        'code': 221, 'severity': 'error',
                        'path': self._ASSEMBLY_PATH,
                        'text': f'Copy {i} names chain {chain}, which the '
                                f'model does not have.'})

        bodies = self._preview_bodies(model, instances, messages)
        return {
            'ok': True,
            'model': {'chains': chains, 'entities': entities,
                      'nCopiesDetected': len(detected)},
            'suggestion': suggestion,
            'assembly': {
                'source': source,
                'rows': (rows if source == 'parameter'
                         else suggestion['assembly']),
                'instances': [{'label': dm_ncs_lib.instance_label(inst),
                               'roles': dict(inst)} for inst in instances],
            },
            'bodies': bodies,
            'messages': messages,
        }

    def _preview_bodies(self, model, instances, messages):
        """Per rigid body: its parsed segments, and against each copy the
        number of matched CA atoms and the achieved superposition RMSD."""
        ctrl = self.container.controlParameters
        bodies = []
        ref = instances[0] if instances else {}
        for i, d in enumerate(ctrl.DOMAINS, start=1):
            spec = d.segments_spec()
            entry = {'index': i, 'spec': spec, 'mode': d.averaging_mode(),
                     'segments': [], 'copies': [], 'error': None}
            bodies.append(entry)
            if not spec:
                entry['error'] = 'no residue ranges'
                continue
            try:
                segments = dm_ncs_lib.parse_segments(spec)
            except ValueError as exc:
                entry['error'] = str(exc)
                continue
            entry['segments'] = [{'role': role, 'lo': lo, 'hi': hi}
                                 for role, lo, hi in segments]
            if entry['mode'] == 'exclude' or not ref:
                continue

            roles = dm_ncs_lib.body_segment_roles(segments)
            missing = [r for r in roles if not ref.get(r)]
            if missing:
                entry['error'] = ('reference copy has no chain for '
                                  + ', '.join(sorted(missing)))
                continue
            n_ref = len(dm_ncs_lib.body_ca_positions(model, ref, segments))
            entry['nReferenceCA'] = n_ref
            if n_ref == 0:
                entry['error'] = 'no CA atoms in these ranges'
                messages.append({
                    'code': 222, 'severity': 'error',
                    'path': self._DOMAINS_PATH,
                    'text': f'Rigid body {i} ({spec}) covers no CA atoms of '
                            f'the reference copy.'})
                continue

            for copy in instances[1:]:
                label = dm_ncs_lib.instance_label(copy)
                if not all(copy.get(r) for r in roles):
                    entry['copies'].append({'label': label, 'skipped': True})
                    continue
                try:
                    _, rmsd = dm_ncs_lib.operator_ref_to_copy_body(
                        model, ref, copy, segments)
                except ValueError as exc:
                    entry['copies'].append({'label': label, 'error': str(exc)})
                    continue
                ref_pts, _ = dm_ncs_lib.body_ca_correspondence(
                    model, ref, copy, segments)
                entry['copies'].append({'label': label, 'nCA': len(ref_pts),
                                        'rmsd': round(rmsd, 3)})
            matched = [c for c in entry['copies'] if 'rmsd' in c]
            if not matched:
                messages.append({
                    'code': 223, 'severity': 'error',
                    'path': self._DOMAINS_PATH,
                    'text': f'Rigid body {i} ({spec}) has no copy to average '
                            f'against.'})
            elif max(c['rmsd'] for c in matched) > 3.0:
                worst = max(matched, key=lambda c: c['rmsd'])
                messages.append({
                    'code': 224, 'severity': 'warning',
                    'path': self._DOMAINS_PATH,
                    'text': f'Rigid body {i} ({spec}) superposes on copy '
                            f'{worst["label"]} to only {worst["rmsd"]:.1f} A '
                            f'-- these residues may not move as one unit.'})
        return bodies

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

            # NCS instances (role->chain per copy), index 0 == reference. From
            # ASSEMBLY if supplied, else auto-detect a homomer (first protein
            # chain reference, the rest its copies).
            if ctrl.ASSEMBLY.isSet() and len(ctrl.ASSEMBLY) > 0:
                instances = dm_ncs_lib.parse_assembly_rows(
                    [str(r) for r in ctrl.ASSEMBLY])
            else:
                # Same detection the task interface offers as its suggestion,
                # so a job run without touching ASSEMBLY does what the
                # interface said it would (hetero-complexes included, not
                # just homomers).
                instances, _ = dm_ncs_lib.detect_assembly(model)
            if not instances:
                self.appendErrorReport(
                    201, 'No NCS instances (ASSEMBLY empty and auto-detect '
                         'found no protein chains)')
                return CPluginScript.FAILED
            ref_inst = instances[0]
            print(f"dm_multidomain: {len(instances)} instances, reference "
                  f"{dm_ncs_lib.instance_label(ref_inst)}")

            # rigid bodies: each a list of (role, lo, hi) segments + a mode.
            domains = []
            for i, d in enumerate(ctrl.DOMAINS, start=1):
                spec = d.segments_spec()
                if not spec:
                    raise ValueError(f"domain {i}: no segments set")
                segments = dm_ncs_lib.parse_segments(spec)
                missing = [r for r in dm_ncs_lib.body_segment_roles(segments)
                           if not ref_inst.get(r)]
                if missing:
                    raise ValueError(
                        f"domain {i} ({spec}): reference instance "
                        f"{dm_ncs_lib.instance_label(ref_inst)} lacks role(s) "
                        f"{missing} -- add them to the reference ASSEMBLY row")
                domains.append(dict(name=dm_ncs_lib.body_name(segments),
                                    segments=segments, mode=d.averaging_mode()))

            radius = float(ctrl.MASK_RADIUS) if ctrl.MASK_RADIUS.isSet() else 2.5
            # 1. operators per averaged body (skip excluded; skip bodies no copy
            #    shares -- nothing to average against).
            operators_by_domain = {}
            rmsds_by_domain = {}
            averaged = []   # bodies that survive, in command-line order
            for d in domains:
                if d['mode'] == 'exclude':
                    continue
                ops, rmsds = dm_ncs_lib.body_operators(
                    model, instances, d['segments'], cell=cell)
                rmsds_by_domain[d['name']] = rmsds
                if len(ops) < 2:
                    print(f"  body {d['name']} ({d['mode']}): no copies share "
                          f"its roles -- skipped")
                    continue
                operators_by_domain[d['name']] = ops
                rmsd_str = ", ".join(f"{c}:{r:.2f}" for c, r in rmsds.items())
                print(f"  body {d['name']} ({d['mode']}): RMSD {rmsd_str}")
                averaged.append(d)

            # 2. masks for the surviving bodies, written together so a nearest-
            #    atom competition makes them DISJOINT (dm needs non-overlapping
            #    NCS masks; the radius dilation otherwise bleeds across abutting
            #    boundaries).
            ncsin = []   # ordered [(domain_name, mask_path)], in step with above
            if averaged:
                # .map (not .msk) so the gleaned MASKOUT is an unambiguous CCP4
                # map for the viewer; dm reads it as NCSIN regardless.
                mask_paths = [os.path.join(self.workDirectory,
                                           f"mask_{d['name']}.map")
                              for d in averaged]
                nsets, n_contested = dm_ncs_lib.write_partitioned_masks(
                    model, ref_inst, [d['segments'] for d in averaged],
                    cell, sg, mask_paths, radius=radius)
                if n_contested:
                    print(f"  resolved {n_contested} overlapping mask points "
                          f"(nearest-atom partition) -> disjoint masks")
                for d, path, nset in zip(averaged, mask_paths, nsets):
                    print(f"    mask {os.path.basename(path)}: {nset} points")
                    ncsin.append((d['name'], path))

            # keep the script's domain list in step with the NCSIN<n> order,
            # but remember every body the user asked for (excluded ones
            # included) so the report can draw the whole partition.
            self._all_domains = domains
            self._instances = instances
            self._rmsds_by_domain = rmsds_by_domain
            domains = averaged
            if not ncsin:
                self.appendErrorReport(
                    201, 'No domains to average (all excluded, or no copies '
                         'share any body\'s roles)')
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

        # capture the per-body averaging masks (disjoint mode-0 CCP4 maps) as
        # outputs, so they can be gleaned and overlaid (e.g. in Moorhen) to
        # inspect/verify the domain partition. self._ncsin is [(name, path)] in
        # NCSIN order; the files were written by write_partitioned_masks.
        maskout = self.container.outputData.MASKOUT
        for name, path in getattr(self, '_ncsin', []):
            if os.path.exists(path):
                maskout.append(maskout.makeItem())
                maskout[-1].setFullPath(path)
                # Mark as a mask (not a density map): gleaned to File.sub_type
                # so the Moorhen viewers / scene format render it as a mask.
                maskout[-1].subType.set(maskout[-1].SUBTYPE_MASK)
                maskout[-1].annotation = \
                    self.jobNumberString() + f' NCS averaging mask: {name}'

        # build the result XML: loggraph tables + per-cycle + per-domain NCS
        # correlations, as substrate for a graphically rich report.
        rootNode = etree.Element("DmMultidomainResult")
        logtext = self._read_log()
        self._add_smartie_graphs(etree.SubElement(rootNode, 'SmartieGraphs'))
        self._add_per_cycle(rootNode, logtext)
        self._add_ncs_correlations(rootNode, logtext)
        self._add_body_map(rootNode)
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

    def _add_body_map(self, rootNode):
        """The rigid-body partition itself: which residues of which entity each
        body claimed, and how well it superposed.

        The report otherwise says how averaging went without ever saying what
        was averaged, which is the one thing a reader coming back to a job
        needs in order to judge the rest. It is the same picture the task
        interface draws while the job is being set up.
        """
        import gemmi

        domains = getattr(self, '_all_domains', None)
        instances = getattr(self, '_instances', None)
        if not domains or not instances:
            return
        try:
            structure = gemmi.read_structure(
                str(self.container.inputData.XYZIN.fullPath))
            structure.setup_entities()
            model = structure[0]
        except Exception:
            return

        node = etree.SubElement(rootNode, 'BodyMap')
        reference = instances[0]
        for role, chain in reference.items():
            bounds = dm_ncs_lib.residue_bounds(model, chain)
            if bounds is None:
                continue
            track = etree.SubElement(node, 'Track')
            track.set('role', role)
            track.set('chain', chain)
            track.set('first', str(bounds[0]))
            track.set('last', str(bounds[1]))

        rmsds = getattr(self, '_rmsds_by_domain', {})
        for i, d in enumerate(domains, start=1):
            body = etree.SubElement(node, 'Body')
            body.set('number', str(i))
            body.set('mode', d['mode'])
            body.set('name', d['name'])
            for role, lo, hi in d['segments']:
                segment = etree.SubElement(body, 'Segment')
                segment.set('role', role)
                segment.set('lo', str(lo))
                segment.set('hi', str(hi))
            for label, rmsd in sorted(rmsds.get(d['name'], {}).items()):
                fit = etree.SubElement(body, 'Fit')
                fit.set('copy', label)
                fit.set('rmsd', f'{rmsd:.2f}')

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
