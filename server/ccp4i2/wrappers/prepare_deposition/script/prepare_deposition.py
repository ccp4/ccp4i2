"""Prepare the files a wwPDB OneDep deposition takes: the model mmCIF with
the data statistics, one entity per crystallised molecule and the descriptive
metadata written into it, and the reflections mmCIF; then the wwPDB
validation report.

Everything is done with gemmi (ccp4i2.lib.deposition); no external package
rewrites the file. The model is sanitised first (the things BUSTER, Refmac
and servalcat each leave in their mmCIF), and every change is reported.
"""
import os
import sys

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class prepare_deposition(CPluginScript):
    TASKNAME = 'prepare_deposition'
    ERROR_CODES = {
        201: {'description': 'The model is not fit to deposit (see the problems listed)'},
        202: {'description': 'The reflections mmCIF could not be written'},
        203: {'description': 'The scaling statistics could not be merged'},
        204: {'description': 'The metadata could not be written'},
        205: {'description': 'No sequence selected in the AU contents: the deposited model '
                             'needs the full sequence of every chain'},
        206: {'description': 'The written categories do not conform to the PDBx dictionary'},
        207: {'description': 'A title, at least one author, and a contact name and email are '
                             'needed before the files can be prepared'},
    }

    def __init__(self, *args, **kws):
        super().__init__(*args, **kws)
        self.xmlroot = etree.Element('prepare_deposition')

    # -- interactive helpers (plugin_method) -------------------------------

    def trace_xyzin_lineage(self, file_uuid):
        """What the task interface calls when a model is chosen: the files and
        settings of the refinement and scaling jobs, plus the refinement job
        directory and the program versions found."""
        from ccp4i2.lib.deposition.lineage import trace_xyzin_lineage
        from ccp4i2.db.models import Job

        result = trace_xyzin_lineage(file_uuid)
        job_id = result.get('refinement_job_id')
        if job_id:
            try:
                result['params']['REFINEMENT_JOB_DIR'] = str(Job.objects.get(id=job_id).directory)
            except Job.DoesNotExist:
                pass
        return result

    def suggest_metadata(self):
        """What can be pre-filled from the inputs: one molecule per AU-contents
        sequence (name, length, copies) and the wavelength. Never raises."""
        try:
            return self._suggest_metadata()
        except Exception as exc:  # pragma: no cover - guard rail
            return {'ok': False, 'error': '%s: %s' % (type(exc).__name__, exc)}

    def _suggest_metadata(self):
        inp = self.container.inputData
        out = {'ok': True, 'entities': [], 'wavelength': None}
        if inp.ASUCONTENT.isSet():
            asu = inp.ASUCONTENT
            asu.loadFile()
            for seq in asu.fileContent.seqList:
                if not asu.isSelected(seq):
                    continue
                sequence = ''.join(str(seq.sequence).split())
                out['entities'].append({
                    'SEQUENCE_NAME': str(seq.name),
                    'DESCRIPTION': str(getattr(seq, 'description', '') or seq.name),
                    'length': len(sequence),
                    'copies': int(getattr(seq, 'nCopies', 1) or 1),
                })
        if inp.F_SIGF.isSet():
            import gemmi
            try:
                mtz = gemmi.read_mtz_file(str(inp.F_SIGF.fullPath))
                for ds in mtz.datasets:
                    if ds.wavelength and 0.1 < ds.wavelength < 5.0:
                        out['wavelength'] = round(ds.wavelength, 5)
                        break
            except Exception as err:
                out['wavelength_error'] = str(err)
        return out

    # -- validation ----------------------------------------------------------

    def validity(self):
        error = super().validity()
        ctrl = self.container.controlParameters
        missing = []
        if not str(ctrl.TITLE).strip():
            missing.append('a title')
        if not any(_text(_field(a, 'NAME')) for a in ctrl.AUTHORS):
            missing.append('at least one author')
        if not str(ctrl.CONTACT_LAST_NAME).strip() or not str(ctrl.CONTACT_EMAIL).strip():
            missing.append('the contact author name and email')
        if not _text(ctrl.KEYWORDS) and not self._all_designed():
            missing.append('a classification (one term; only a wholly designed entry defaults '
                           'to DE NOVO PROTEIN)')
        if missing:
            error.append(klass=self.TASKNAME, code=207,
                         details='Still needed: ' + ', '.join(missing),
                         name='%s.container.controlParameters.TITLE' % self.TASKNAME,
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def _all_designed(self):
        """True when every segment of every molecule is designed (tags and
        linkers aside), so DE NOVO PROTEIN is the classification."""
        from ccp4i2.lib.deposition.metadata import default_keywords
        try:
            return default_keywords(self._metadata_spec().get('entities')) is not None
        except Exception:
            return False

    # -- the work --------------------------------------------------------

    def startProcess(self):
        import gemmi

        from ccp4i2.lib.deposition import mmcif_edit, sanitise, statistics
        from ccp4i2.lib.deposition.entities import rebuild_entities_in
        from ccp4i2.lib.deposition.metadata import write_metadata

        inp, ctrl, out = self.container.inputData, self.container.controlParameters, self.container.outputData

        sequences = self._sequences()
        if not sequences:
            self.appendErrorReport(205, 'no sequences selected in the AU contents')
            return CPluginScript.FAILED

        # 1. Read and sanitise what the refinement program wrote.
        from ccp4i2.lib.deposition.entities import read_model_document
        doc = read_model_document(str(inp.XYZIN.fullPath))
        job_dir = str(ctrl.REFINEMENT_JOB_DIR) if ctrl.REFINEMENT_JOB_DIR.isSet() else None
        if not job_dir and inp.REFMACINPUTPARAMSXML.isSet():
            job_dir = os.path.dirname(str(inp.REFMACINPUTPARAMSXML.fullPath))
        program, fixes, problems, _st = sanitise.sanitise(
            doc, refinement_job_dir=job_dir, keep_hydrogens=bool(ctrl.KEEP_HYDROGENS))
        self._report_fixes(program, fixes)
        cell_problem = sanitise.cell_matches_reflections(doc[0], str(inp.F_SIGF.fullPath))
        if cell_problem:
            problems.append(cell_problem)

        # 2. One entity per crystallised molecule, label_seq_id on the construct.
        entity_warnings = []
        st = rebuild_entities_in(doc, sequences, warnings=entity_warnings)
        for w in entity_warnings:
            print('WARNING: %s' % w)
        block = doc[0]

        # 3. Data reduction statistics.
        stats_path = self._statistics_path()
        if stats_path:
            try:
                written = statistics.merge_statistics(
                    block, stats_path, all_shells=bool(ctrl.KEEP_ALL_SHELLS))
                print('Merged scaling statistics: %s' % ', '.join(written))
                self.xmlroot.append(statistics.cifstats_to_xml(stats_path))
            except Exception as err:
                self.appendErrorReport(203, str(err))
                return CPluginScript.FAILED
        elif inp.AIMLESSXML.isSet() and bool(ctrl.USEAIMLESSXML):
            self._aimless_xml_to_report()
            print('WARNING: the scaling job has no DataStatistics.cif; its Aimless XML statistics '
                  'are shown in the report but not merged into the model (enter them in OneDep)')

        # 4. The descriptive metadata.
        spec = self._metadata_spec()
        try:
            warnings = write_metadata(block, st, spec)
        except ValueError as err:
            self.appendErrorReport(204, str(err))
            return CPluginScript.FAILED
        self._report_metadata(spec, entity_warnings + warnings)
        if mmcif_edit.value_of(block, '_diffrn_radiation_wavelength.wavelength') is None:
            problems.append('no wavelength: the model has none and neither the Experiment page '
                            'nor the reflection data gave one')

        # 5. The gates: no split category, no problem left, dictionary clean.
        split = mmcif_edit.split_categories(block)
        if split:
            problems.append('categories written in two places: ' + ', '.join(split))
        if bool(ctrl.DICTIONARY_CHECK):
            messages = self._dictionary_check(block)
            if messages:
                problems.append('PDBx dictionary: ' + '; '.join(messages[:10]))
        if problems:
            self._report_problems(problems)
            self.flushXml()
            self.appendErrorReport(201, '; '.join(problems))
            return CPluginScript.FAILED

        model_out = str(out.MMCIFOUT.fullPath)
        doc.write_file(model_out, gemmi.cif.Style.Pdbx)
        print('Wrote %s' % model_out)

        # 6. The reflections.
        rv = self._write_reflections()
        if rv is not None:
            return rv

        # 7. Sequence check for the report, then validation.
        try:
            from ccp4i2.wrappers.modelASUCheck.script.modelASUCheck import sequenceAlignment
            self.xmlroot.append(sequenceAlignment(str(inp.XYZIN.fullPath), inp.ASUCONTENT))
        except Exception as err:
            print('Sequence alignment check failed (non-fatal): %s' % err)
        if bool(ctrl.SENDTOVALIDATIONSERVER):
            etree.SubElement(self.xmlroot, 'Validation').set('sent', 'true')
            self.flushXml()
            self._validate(model_out, str(out.CIFREFLECTIONS.fullPath))
        else:
            self.flushXml()
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        out = self.container.outputData
        # The Coot script reads the validation XML: only when there is one.
        if os.path.exists(str(out.VALIDATIONXML.fullPath)):
            with open(str(out.COOTSCRIPTOUT.fullPath), 'w') as script:
                script.write(COOT_SCRIPT.format(xml=str(out.VALIDATIONXML.fullPath)))
            out.COOTSCRIPTOUT.annotation = 'Coot script to review the validation issues'
        return CPluginScript.SUCCEEDED

    # -- pieces ----------------------------------------------------------

    def _sequences(self):
        asu = self.container.inputData.ASUCONTENT
        asu.loadFile()
        return [(str(seq.name), str(seq.sequence), str(seq.polymerType))
                for seq in asu.fileContent.seqList
                if asu.isSelected(seq) and str(seq.sequence).strip()]

    def _statistics_path(self):
        inp, ctrl = self.container.inputData, self.container.controlParameters
        if not bool(ctrl.USEAIMLESSXML):
            return None
        if inp.CIFSTATSOUT.isSet() and os.path.exists(str(inp.CIFSTATSOUT.fullPath)):
            return str(inp.CIFSTATSOUT.fullPath)
        return None

    def _aimless_xml_to_report(self):
        path = str(self.container.inputData.AIMLESSXML.fullPath)
        try:
            root = etree.parse(path).getroot()
            nodes = root.xpath('.//AIMLESS_PIPE') or root.xpath('.//AIMLESS')
            if nodes:
                self.xmlroot.append(nodes[-1])
        except Exception as err:
            print('Aimless XML not readable: %s' % err)

    def _metadata_spec(self):
        """The plain dict ccp4i2.lib.deposition.metadata.write_metadata takes,
        from the task's parameters."""
        ctrl = self.container.controlParameters
        inp = self.container.inputData

        def s(item):
            return _text(item)

        def n(item):
            text = _text(item)
            try:
                return float(text) if text is not None else None
            except (TypeError, ValueError):
                return None

        def f(obj, name):
            return _field(obj, name)

        authors = [{'name': s(f(a, 'NAME')), 'orcid': s(f(a, 'ORCID'))} for a in ctrl.AUTHORS
                   if s(f(a, 'NAME'))]
        cit_authors = [{'name': s(f(a, 'NAME')), 'orcid': s(f(a, 'ORCID'))}
                       for a in ctrl.CITATION_AUTHORS if s(f(a, 'NAME'))]
        sequences = self._sequences()
        entities = []
        for i, ent in enumerate(ctrl.ENTITIES):
            entity_id = str(i + 1)
            name = s(f(ent, 'SEQUENCE_NAME'))
            for j, (seq_name, _seq, _ptype) in enumerate(sequences):
                if name and seq_name == name:
                    entity_id = str(j + 1)
            segments = []
            for seg in (f(ent, 'SEGMENTS') or []):
                beg, end = n(f(seg, 'BEG')), n(f(seg, 'END'))
                segments.append({
                    'beg': int(beg) if beg else 1,
                    'end': int(end) if end else None,
                    'kind': s(f(seg, 'KIND')) or 'natural',
                    'organism': s(f(seg, 'ORGANISM')), 'taxid': _positive(f(seg, 'TAXID')),
                    'strain': s(f(seg, 'STRAIN')), 'gene': s(f(seg, 'GENE')),
                    'uniprot': s(f(seg, 'UNIPROT')), 'pdb_id': s(f(seg, 'PDB_ID')),
                    'engineered': s(f(seg, 'ENGINEERED')), 'host': s(f(seg, 'HOST')),
                    'host_taxid': _positive(f(seg, 'HOST_TAXID')), 'details': s(f(seg, 'DETAILS')),
                })
            entities.append({
                'entity_id': entity_id, 'description': s(f(ent, 'DESCRIPTION')),
                'common_name': s(f(ent, 'COMMON_NAME')), 'ec': s(f(ent, 'EC')),
                'mutation': s(f(ent, 'MUTATION')), 'details': s(f(ent, 'DETAILS')),
                'segments': segments,
            })
        software = [{'name': s(f(sw, 'NAME')), 'version': s(f(sw, 'VERSION')),
                     'classification': s(f(sw, 'CLASSIFICATION')) or 'other'}
                    for sw in ctrl.SOFTWARE if s(f(sw, 'NAME'))]
        software.extend(self._lineage_software())
        software.append({'name': 'CCP4i2', 'version': _ccp4i2_version(),
                         'classification': 'data processing'})
        wavelength = n(ctrl.WAVELENGTH)
        if wavelength is None and inp.F_SIGF.isSet():
            try:
                import gemmi
                for ds in gemmi.read_mtz_file(str(inp.F_SIGF.fullPath)).datasets:
                    if ds.wavelength and 0.1 < ds.wavelength < 5.0:
                        wavelength = round(ds.wavelength, 5)
                        break   # the first plausible one; HKL_base carries a dummy
            except Exception:
                wavelength = None
        method = s(ctrl.PHASING_METHOD)
        starting_models = [{'type': s(f(m, 'TYPE')) or 'experimental model',
                            'source': s(f(m, 'SOURCE')) or 'PDB',
                            'accession': s(f(m, 'ACCESSION')), 'details': s(f(m, 'DETAILS')),
                            'entity_ids': s(f(m, 'ENTITY_IDS'))}
                           for m in ctrl.STARTING_MODELS
                           if s(f(m, 'ACCESSION')) or s(f(m, 'DETAILS'))]
        return {
            'title': s(ctrl.TITLE), 'keywords': s(ctrl.KEYWORDS), 'keyword_text': s(ctrl.KEYWORD_TEXT),
            'authors': authors,
            'citation': {'published': bool(ctrl.CITATION_PUBLISHED), 'title': s(ctrl.CITATION_TITLE),
                         'journal': s(ctrl.JOURNAL), 'volume': s(ctrl.VOLUME),
                         'first_page': s(ctrl.FIRST_PAGE), 'last_page': s(ctrl.LAST_PAGE),
                         'year': _positive(ctrl.YEAR), 'doi': s(ctrl.DOI), 'pubmed': s(ctrl.PUBMED),
                         'authors': cit_authors},
            'contact': {'salutation': s(ctrl.CONTACT_SALUTATION), 'first_name': s(ctrl.CONTACT_FIRST_NAME),
                        'last_name': s(ctrl.CONTACT_LAST_NAME), 'role': s(ctrl.CONTACT_ROLE),
                        'organization_type': s(ctrl.CONTACT_ORG_TYPE),
                        'address_1': s(ctrl.CONTACT_ADDRESS_1), 'address_2': s(ctrl.CONTACT_ADDRESS_2),
                        'city': s(ctrl.CONTACT_CITY), 'state': s(ctrl.CONTACT_STATE),
                        'postal_code': s(ctrl.CONTACT_POSTCODE), 'country': s(ctrl.CONTACT_COUNTRY),
                        'email': s(ctrl.CONTACT_EMAIL), 'phone': s(ctrl.CONTACT_PHONE),
                        'orcid': s(ctrl.CONTACT_ORCID)},
            'related': [{'db_name': s(f(r, 'DB_NAME')) or 'PDB', 'db_id': s(f(r, 'DB_ID')),
                         'content_type': s(f(r, 'CONTENT_TYPE')) or 'unspecified',
                         'details': s(f(r, 'DETAILS'))}
                        for r in ctrl.RELATED if s(f(r, 'DB_ID'))],
            'entities': entities,
            'design_details': s(ctrl.DESIGN_DETAILS),
            'experiment': {'source': s(ctrl.SOURCE) or 'SYNCHROTRON', 'site': s(ctrl.SYNCHROTRON_SITE),
                           'beamline': s(ctrl.BEAMLINE), 'detector': s(ctrl.DETECTOR),
                           'detector_type': s(ctrl.DETECTOR_TYPE),
                           'collection_date': s(ctrl.COLLECTION_DATE), 'temperature': n(ctrl.TEMPERATURE),
                           'wavelength': wavelength, 'monochromator': s(ctrl.MONOCHROMATOR),
                           'protocol': s(ctrl.PROTOCOL)},
            'crystal': {'method': s(ctrl.CRYSTAL_METHOD), 'temp': n(ctrl.CRYSTAL_TEMP),
                        'ph': n(ctrl.CRYSTAL_PH), 'details': s(ctrl.CRYSTAL_DETAILS)},
            'software': software,
            'phasing': {'method': method, 'starting_model': s(ctrl.STARTING_MODEL),
                        'starting_models': starting_models,
                        'details': s(ctrl.PHASING_DETAILS)} if method else {'starting_models': starting_models},
            'ligand_of_interest': s(ctrl.LIGAND_OF_INTEREST),
            'assembly': {'mode': 'author' if s(ctrl.ASSEMBLY_MODE) == 'AUTHOR' else 'asu',
                         'oligomeric_count': int(ctrl.OLIGOMERIC_COUNT) if ctrl.OLIGOMERIC_COUNT.isSet() else None,
                         'operators': s(ctrl.ASSEMBLY_OPERATORS) or '',
                         'evidence': s(ctrl.ASSEMBLY_EVIDENCE) or 'none',
                         'details': s(ctrl.ASSEMBLY_DETAILS)},
        }

    def _lineage_software(self):
        """Program versions from the scaling job's XML, when the interface
        filled AIMLESSXML or the statistics file sits in that job."""
        inp = self.container.inputData
        rows = []
        for item in (inp.AIMLESSXML,):
            if not item.isSet() or not os.path.exists(str(item.fullPath)):
                continue
            try:
                root = etree.parse(str(item.fullPath)).getroot()
            except etree.XMLSyntaxError:
                continue
            for program, label in (('AIMLESS', 'Aimless'), ('POINTLESS', 'Pointless')):
                node = root.find('.//%s' % program)
                if node is not None and node.get('version'):
                    rows.append({'name': label, 'version': node.get('version'),
                                 'classification': 'data scaling'})
        return rows

    def _dictionary_check(self, block):
        from ccp4i2.lib.deposition import mmcif_edit
        from ccp4i2.lib.deposition.metadata import WRITTEN_CATEGORIES

        dictionary = _dictionary_path()
        if dictionary is None:
            print('WARNING: no PDBx dictionary found in $CCP4; the dictionary check is skipped')
            return []
        return mmcif_edit.check_categories(block, dictionary, WRITTEN_CATEGORIES)

    def _write_reflections(self):
        from ccp4i2.lib.deposition.reflections import write_reflections_cif

        inp, ctrl, out = self.container.inputData, self.container.controlParameters, self.container.outputData
        if bool(ctrl.USEANOMALOUS):
            refmac_flag = 1 if bool(ctrl.USE_TWIN) else 2
        else:
            refmac_flag = 3 if bool(ctrl.USE_TWIN) else 4
        merged = str(inp.F_SIGF.fullPath)
        if int(inp.F_SIGF.contentFlag) != refmac_flag or inp.FREERFLAG.isSet():
            mini = [['F_SIGF', int(inp.F_SIGF.contentFlag)], ['F_SIGF', refmac_flag]]
            mini += [[name, 1] for name in ('FREERFLAG', 'FPHIOUT', 'DIFFPHIOUT')
                     if getattr(inp, name).isSet()]   # an unset file makes makeHklin0 raise
            self.hklin, self.columns, error = self.makeHklin0(
                miniMtzsIn=mini, hklin='mergedForMakingCif')
            if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
                self.appendErrorReport(202, 'could not merge the reflection files')
                return CPluginScript.FAILED
            merged = os.path.join(self.getWorkDirectory(), 'mergedForMakingCif.mtz')
        unmerged = None
        if bool(ctrl.INCLUDEUNMERGED) and inp.SCALEDUNMERGED.isSet():
            unmerged = str(inp.SCALEDUNMERGED.fullPath)
        try:
            write_reflections_cif(merged, str(out.CIFREFLECTIONS.fullPath), unmerged, log=print)
        except Exception as err:
            self.appendErrorReport(202, str(err))
            return CPluginScript.FAILED
        return None

    def _validate(self, model_path, sf_path):
        from ccp4i2.lib.deposition.onedep_validation import run_validation

        out = self.container.outputData
        work = self.getWorkDirectory()
        try:
            status = run_validation(model_path, sf_path, str(out.PDFOUT.fullPath),
                                    str(out.VALIDATIONXML.fullPath),
                                    os.path.join(work, 'validation.svg'),
                                    os.path.join(work, 'validation.log'), log=print)
        except Exception as err:
            print('wwPDB validation failed (non-fatal): %s' % err)
            return
        if status == 'completed':
            out.PDFOUT.annotation.set('wwPDB validation report')
            try:
                with open(str(out.VALIDATIONXML.fullPath), 'r') as handle:
                    self.xmlroot.append(etree.fromstring(handle.read().encode('utf-8')))
            except Exception as err:
                print('validation XML not readable: %s' % err)
            self.flushXml()

    # -- report XML -------------------------------------------------------

    def _report_fixes(self, program, fixes):
        node = etree.SubElement(self.xmlroot, 'Sanitise')
        etree.SubElement(node, 'Program').text = program
        for fix in fixes:
            f = etree.SubElement(node, 'Fix')
            etree.SubElement(f, 'Found').text = fix.found
            etree.SubElement(f, 'Did').text = fix.did
            print('[%s] found %s; %s' % (fix.program, fix.found, fix.did))
        self.flushXml()

    def _report_problems(self, problems):
        node = etree.SubElement(self.xmlroot, 'Problems')
        for p in problems:
            etree.SubElement(node, 'Problem').text = p
            print('PROBLEM: %s' % p)

    def _report_metadata(self, spec, warnings):
        node = etree.SubElement(self.xmlroot, 'Metadata')

        def item(label, text):
            if text:
                etree.SubElement(node, 'Item', label=label).text = str(text)

        item('Title', spec.get('title'))
        item('Classification', spec.get('keywords'))
        item('Authors', ', '.join(a['name'] for a in spec.get('authors', [])))
        item('Contact', '%s %s <%s>' % (spec['contact'].get('first_name') or '',
                                        spec['contact'].get('last_name') or '',
                                        spec['contact'].get('email') or ''))
        for ent in spec.get('entities', []):
            item('Molecule %s' % ent['entity_id'],
                 '%s; segments: %s' % (ent.get('description'), '; '.join(
                     '%s-%s %s' % (sg.get('beg'), sg.get('end'), sg.get('kind'))
                     for sg in ent.get('segments', []))))
        exp = spec.get('experiment', {})
        item('Data collection', ' '.join(str(v) for v in (
            exp.get('site'), exp.get('beamline'), exp.get('detector_type'),
            exp.get('collection_date'), '%s K' % exp.get('temperature') if exp.get('temperature') else None,
            '%s A' % exp.get('wavelength') if exp.get('wavelength') else None) if v))
        item('Software', ', '.join('%s %s' % (sw['name'], sw.get('version') or '')
                                   for sw in spec.get('software', [])))
        item('Assembly', '%s%s; evidence: %s' % (
            spec['assembly'].get('mode'),
            ' x%s' % spec['assembly'].get('oligomeric_count') if spec['assembly'].get('oligomeric_count') else '',
            spec['assembly'].get('evidence')))
        item('Starting models', '; '.join(
            '%s / %s / %s' % (m.get('type'), m.get('source'), m.get('accession') or m.get('details'))
            for m in spec.get('phasing', {}).get('starting_models', [])))
        item('Ligand of interest', spec.get('ligand_of_interest') or
             'none: tick "No ligand of interest" on OneDep\'s Ligands page')
        for w in warnings:
            etree.SubElement(node, 'Warning').text = w
            print('WARNING: %s' % w)

    def flushXml(self):
        with open(self.makeFileName('PROGRAMXML'), 'w') as handle:
            CCP4Utils.writeXML(handle, etree.tostring(self.xmlroot, pretty_print=True))
        sys.stdout.flush()


def _field(obj, name):
    """A field of a list item, whether it is the typed CData the def.xml
    makes or a plain dict (what a scripted caller may append)."""
    if isinstance(obj, dict):
        return obj.get(name)
    return getattr(obj, name, None)


def _text(item):
    """The value of a CData scalar (or a plain value) as stripped text, or None when unset."""
    if item is None:
        return None
    if hasattr(item, 'isSet'):
        if not item.isSet():
            return None
        text = str(item).strip()
    else:
        text = str(item).strip()
    return text if text not in ('', 'None') else None


def _positive(item):
    """An integer field as text, or None when unset -- an unset CInt comes
    back from params.xml as 0, which is not a taxonomy id or a year."""
    text = _text(item)
    try:
        return str(int(float(text))) if text is not None and int(float(text)) > 0 else None
    except (TypeError, ValueError):
        return None


def _dictionary_path():
    ccp4 = os.environ.get('CCP4')
    if not ccp4:
        return None
    for candidate in ('share/libcifpp/mmcif_pdbx.dic',):
        path = os.path.join(ccp4, candidate)
        if os.path.exists(path):
            return path
    return None


def _ccp4i2_version():
    try:
        from ccp4i2 import __version__
        return __version__
    except Exception:
        return '?'


COOT_SCRIPT = '''
try:
    xml_string = open('{xml}').read()
    valid_inf = parse_wwpdb_validation_xml(xml_string)
    if valid_inf:
        entry_validation_info = valid_inf[0]
        subgroups = valid_inf[1]
        ss = sort_subgroups(subgroups)
        validation_to_gui(entry_validation_info, ss, 0)
    else:
        print('problem with valid_inf')
except IOError as e_mess:
    print('load_validate_xml IOERROR: {{}}'.format(e_mess))
except Exception as err:
    print("Exception {{}}".format(err))
'''
