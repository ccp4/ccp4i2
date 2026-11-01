"""Every controlled value the deposition task offers, checks and writes comes
from one source: pdb_extract's cif_enumerate.list (the best offline proxy for
OneDep's DepUI menus) and, for the few items it lacks, the PDBx dictionary.
These tests pin that, and the live-test findings from GammaSAD_093825 job 18.
"""
import os

import gemmi
import pytest

from ccp4i2.lib.deposition import enumerations, metadata, segments, vocab
from ccp4i2.lib.deposition._enumerations_data import ENUMERATIONS as VENDORED
from ccp4i2.tests.unit.plugins.test_prepare_deposition_metadata import (
    BASE_SPEC,
    _fetch,
    _prepared,
    _rows,
)

GAMMA = ("MHHHHHHLVPRGSHMIPSITAYSKNGLKIEFTFERSNTNPSVTVITIQASNSTELDMTDFVFQAAVPKTFQLQLLSPSSSV"
         "VPAFNTGTITQVIKVLNPQKQQLRMRIKLTYNHKGSAMQDLAEVNNFPPQSWQ")
AHIR3252_START = ("MTGTHHHHHHSDLGKKLLEAARAGQDDEVRILMANGADVNAKDEYGLTPLYLATAHGHLEIVEVLLKNGADVNAVDAIGF"
                  "TPLHLAAFIGHLEIAEVLLKHGADVNAQDKFGKTAFDISIGNGNEDLAEILQKLGGSGGSGGLQIYYEHDADLGF")
CCP4_LIST = os.path.join(os.environ.get("CCP4", ""), "share", "pdb-extract-data", "cif_enumerate.list")


@pytest.fixture
def plugin(tmp_path):
    from ccp4i2.core.tasks import get_plugin_class

    return get_plugin_class("prepare_deposition")(workDirectory=str(tmp_path))


# --- one source ------------------------------------------------------------

@pytest.mark.skipif(not os.path.exists(CCP4_LIST), reason="no pdb_extract list in $CCP4")
def test_the_vendored_copy_agrees_with_the_ccp4_list():
    from ccp4i2.lib.deposition.make_enumerations import parse_enumerate_list

    with open(CCP4_LIST) as handle:
        listed = parse_enumerate_list(handle.read())
    for item, values in VENDORED.items():
        if item in listed:
            assert set(values) <= set(listed[item]), item


def test_every_menu_of_the_form_is_a_subset_of_the_source(plugin):
    from ccp4i2.wrappers.prepare_deposition.script import prepare_deposition as task
    from ccp4i2.wrappers.prepare_deposition.script import prepare_deposition_types as types

    ctrl = plugin.container.controlParameters
    for field, item in task.FIELD_ITEMS.items():
        menu = getattr(ctrl, field).get_qualifier("enumerators")
        assert menu, field
        assert set(menu) <= set(enumerations.allowed(item)), field
        assert getattr(ctrl, field).get_qualifier("onlyEnumerators") is True, field
    for cls, field, item in ((types.CSoftwareEntry, "CLASSIFICATION", "_software.classification"),
                             (types.CRelatedEntry, "DB_NAME", "_pdbx_database_related.db_name"),
                             (types.CRelatedEntry, "CONTENT_TYPE", "_pdbx_database_related.content_type"),
                             (types.CStartingModel, "TYPE", "_pdbx_initial_refinement_model.type"),
                             (types.CStartingModel, "SOURCE", "_pdbx_initial_refinement_model.source_name")):
        menu = getattr(cls(), field).get_qualifier("enumerators")
        assert menu and set(menu) <= set(enumerations.allowed(item)), (cls.__name__, field)
        default = getattr(cls(), field).get_qualifier("default")
        assert default in menu, (cls.__name__, field, default)


def test_the_menus_survive_loading_saved_parameters(plugin, tmp_path):
    xml = tmp_path / "input_params.xml"
    plugin.container.controlParameters.CONTACT_COUNTRY.set("United Kingdom")
    plugin.container.saveDataToXml(str(xml))
    from ccp4i2.core.tasks import get_plugin_class

    again = get_plugin_class("prepare_deposition")(workDirectory=str(tmp_path))
    again.container.loadDataFromXml(str(xml))
    country = again.container.controlParameters.CONTACT_COUNTRY
    assert str(country) == "United Kingdom"
    assert len(country.get_qualifier("enumerators")) == 249


def test_values_onedep_rejects_are_not_offered():
    assert not {"SIR", "MIRAS", "AB INITIO", "OTHER", "MR"} & set(
        enumerations.allowed("_refine.pdbx_method_to_determine_struct"))
    assert "SAD" not in enumerations.allowed("_diffrn_radiation.pdbx_diffrn_protocol")
    assert "LIQUID ANODE" not in enumerations.allowed("_diffrn_source.source")
    assert len(enumerations.allowed("_diffrn_source.source")) == 7
    assert "other" not in enumerations.allowed("_software.classification")
    assert "electron microscopy" not in vocab.assembly_evidence()
    assert not {"DNA/RNA", "RECOMBINATION", "SURFACTANT PROTEIN", "TRANSLOCASE"} & set(vocab.classifications())
    assert len(enumerations.allowed("_pdbx_database_related.content_type")) == 16
    assert "VAPOR DIFFUSION, SITTING DROP" in enumerations.allowed("_exptl_crystal_grow.method")
    assert "COUNTER-DIFFUSION" not in enumerations.allowed("_exptl_crystal_grow.method")


# --- mapping to the PDB's spelling --------------------------------------------

def test_obvious_spellings_are_mapped():
    canon = enumerations.canonical
    assert canon("_diffrn_source.pdbx_synchrotron_site", "Diamond Light Source") == "Diamond"
    assert canon("_diffrn_source.pdbx_synchrotron_site", "diamond") == "Diamond"
    assert canon("_diffrn_detector.type", "Pilatus 6M-F") == "DECTRIS PILATUS 6M-F"
    assert canon("_diffrn_source.pdbx_synchrotron_beamline", "i04-1") == "I04-1"
    assert canon("_pdbx_contact_author.country", "united kingdom") == "United Kingdom"
    assert canon("_diffrn_source.pdbx_synchrotron_site", "Atlantis") is None
    assert enumerations.source_type("Diamond", "I02") == "DIAMOND BEAMLINE I02"
    assert enumerations.source_type("Diamond", "I99") is None
    assert vocab.phasing_method("MR") == "MOLECULAR REPLACEMENT"
    assert vocab.phasing_method("sad") == "SAD"
    assert vocab.phasing_method("MIRAS") is None


def test_program_names_floor_live_and_unverified():
    assert vocab.software_name("refmacat") == ("REFMAC", vocab.LISTED)
    assert vocab.software_name("parrot") == ("PARROT", vocab.LISTED)
    assert vocab.software_name("Pointless") == ("pointless", vocab.LIVE)   # DepUI menu 2026-11-01
    assert vocab.software_name("Ctruncate") == ("CTRUNCATE", vocab.UNKNOWN)
    assert vocab.software_name("Servalcat") == ("Servalcat", vocab.UNKNOWN)
    assert vocab.software_name("CCP4i2") == (None, None)
    for alias, spelling in vocab._LISTED_ALIASES.items():
        if vocab.software_name(alias)[1] == vocab.LISTED:
            assert enumerations.is_allowed("_software.name", spelling), alias


# --- validity: one error per field, on that field --------------------------------

def _messages(error):
    return [(str(r.get("name", "")), str(r.get("details", "")), r.get("severity"))
            for r in getattr(error, "_reports", getattr(error, "reports", []))]


def _filled(plugin):
    ctrl = plugin.container.controlParameters
    ctrl.TITLE.set("Gamma")
    ctrl.CONTACT_FIRST_NAME.set("Martin")
    ctrl.CONTACT_LAST_NAME.set("Noble")
    ctrl.CONTACT_EMAIL.set("m@example.org")
    author = ctrl.AUTHORS.addItem()
    author.NAME.set("Noble, M.E.M.")
    ctrl.KEYWORDS.set("OXIDOREDUCTASE")
    return ctrl


def _ours(plugin):
    from ccp4i2.core import CCP4ErrorHandling

    out = []
    for r in plugin.validity().getErrors() if hasattr(plugin.validity(), "getErrors") else []:
        out.append(r)
    return out


def _report(plugin):
    """``{field: (details, is_error)}`` for the task's own validity messages."""
    from ccp4i2.core import CCP4ErrorHandling

    error = plugin.validity()
    out = {}
    for rep in error._reports if hasattr(error, "_reports") else error.reports:
        name = rep.get("name") if isinstance(rep, dict) else getattr(rep, "name", "")
        details = rep.get("details") if isinstance(rep, dict) else getattr(rep, "details", "")
        severity = rep.get("severity") if isinstance(rep, dict) else getattr(rep, "severity", 0)
        if name and "controlParameters." in str(name):
            out.setdefault(str(name).split("controlParameters.", 1)[1], []).append(
                (str(details), severity >= CCP4ErrorHandling.SEVERITY_ERROR))
    return out


def test_a_missing_classification_is_reported_on_the_classification(plugin):
    ctrl = _filled(plugin)
    ctrl.KEYWORDS.unSet()
    report = _report(plugin)
    assert "KEYWORDS" in report and report["KEYWORDS"][0][1]
    assert "TITLE" not in report


def test_each_missing_item_on_its_own_field(plugin):
    report = _report(plugin)
    for field in ("TITLE", "AUTHORS", "CONTACT_LAST_NAME", "CONTACT_EMAIL", "KEYWORDS"):
        assert field in report, field


def test_author_without_initials_is_named(plugin):
    ctrl = _filled(plugin)
    ctrl.AUTHORS[0].NAME.set("Noble")
    cit = ctrl.CITATION_AUTHORS.addItem()
    cit.NAME.set("Noble, M.E.M.")
    report = _report(plugin)
    assert "AUTHORS[0].NAME" in report and "Noble, M.E.M." in report["AUTHORS[0].NAME"][0][0]
    assert "CITATION_AUTHORS[0].NAME" not in report
    assert vocab.author_problem("Noble, M.E.M.") is None
    assert vocab.author_problem("Morris, G.M.") is None
    assert vocab.author_problem("Noble, M. E. M.") is not None
    assert vocab.author_problem("Martin Noble") is not None


def test_a_free_text_site_onedep_does_not_know_is_an_error_and_a_long_name_is_mapped(plugin):
    ctrl = _filled(plugin)
    ctrl.SYNCHROTRON_SITE.set("Diamond Light Source")   # what the user typed
    assert "SYNCHROTRON_SITE" not in _report(plugin)       # mapped to "Diamond"
    ctrl.SYNCHROTRON_SITE.set("Atlantis Light Source")
    report = _report(plugin)
    assert "SYNCHROTRON_SITE" in report and "Diamond" in report["SYNCHROTRON_SITE"][0][0]
    ctrl.SYNCHROTRON_SITE.set("Diamond")
    ctrl.BEAMLINE.set("I24")
    assert "BEAMLINE" not in _report(plugin)


def test_home_source_wavelength_at_a_synchrotron_is_a_warning(plugin):
    ctrl = _filled(plugin)
    ctrl.SOURCE.set("SYNCHROTRON")
    ctrl.WAVELENGTH.set(1.54179)
    report = _report(plugin)
    assert "WAVELENGTH" in report and not report["WAVELENGTH"][0][1]
    assert vocab.wavelength_consistency("ROTATING ANODE", 0.9795) is not None
    assert vocab.wavelength_consistency("ROTATING ANODE", 1.5418) is None
    assert vocab.wavelength_consistency("SYNCHROTRON", 0.9795) is None


def test_duplicate_records_are_flagged(plugin):
    ctrl = _filled(plugin)
    for _ in range(2):
        row = ctrl.SOFTWARE.addItem()
        row.NAME.set("Pointless")
        row.CLASSIFICATION.set("data scaling")
        ent = ctrl.ENTITIES.addItem()
        ent.SEQUENCE_NAME.set("gamma")
        ent.DESCRIPTION.set("gamma adaptin ear domain")
    report = _report(plugin)
    assert "SOFTWARE[1].NAME" in report and not report["SOFTWARE[1].NAME"][0][1]
    assert "ENTITIES[1].SEQUENCE_NAME" in report and report["ENTITIES[1].SEQUENCE_NAME"][0][1]


def test_an_unverified_program_is_a_warning_not_an_error(plugin):
    ctrl = _filled(plugin)
    row = ctrl.SOFTWARE.addItem()
    row.NAME.set("Servalcat")
    row.CLASSIFICATION.set("refinement")
    report = _report(plugin)
    assert "SOFTWARE[0].NAME" in report and not report["SOFTWARE[0].NAME"][0][1]
    assert "Communication tab" in report["SOFTWARE[0].NAME"][0][0]


# --- fill: idempotent, with proposed segments ------------------------------------

def test_segments_are_proposed():
    gamma = segments.propose_segments(GAMMA, "gamma", "gamma adaptin ear domain")["segments"]
    assert len(GAMMA) == 134
    assert gamma == [{"BEG": 1, "END": 15, "KIND": "tag"}, {"BEG": 16, "END": 134, "KIND": "natural"}]
    ahir = segments.propose_segments(AHIR3252_START, "Ahir-Darpin-3252")["segments"]
    assert ahir[0] == {"BEG": 1, "END": 11, "KIND": "tag"}
    assert ahir[1] == {"BEG": 12, "END": len(AHIR3252_START), "KIND": "natural"}
    plain = "ACDEFGHIKLMNPQRSTVWY" * 3
    assert segments.propose_segments(plain)["segments"] == [{"BEG": 1, "END": 60, "KIND": "natural"}]
    tagged = segments.propose_segments(plain + "LEHHHHHH", "designed binder")["segments"]
    assert tagged == [{"BEG": 1, "END": 60, "KIND": "designed"}, {"BEG": 61, "END": 68, "KIND": "tag"}]
    tev = segments.propose_segments("MGSSHHHHHHSSGENLYFQG" + plain)["segments"]
    assert tev[0] == {"BEG": 1, "END": 20, "KIND": "tag"}
    assert segments.guess_kind("DARPin H10-2-G3") == "scaffold"
    assert segments.guess_kind("Ahir-Darpin-3252") == "fusion"
    assert segments.guess_kind("designed binder") == "designed"


def test_fill_is_idempotent(plugin, monkeypatch):
    def suggestion():
        return {"ok": True, "wavelength": 0.9795, "entities": [{
            "SEQUENCE_NAME": "gamma", "DESCRIPTION": "gamma adaptin ear domain",
            "SEGMENTS": segments.propose_segments(GAMMA)["segments"], "notes": ["n"],
            "sequence": GAMMA, "length": 134, "copies": 1}]}

    monkeypatch.setattr(plugin, "_suggest_metadata", suggestion)
    lineage = [{"name": "Aimless", "version": "0.8.3", "classification": "data scaling"},
               {"name": "Pointless", "version": "1.15.2", "classification": "data scaling"},
               {"name": "Ctruncate", "version": "1.17.29", "classification": "data scaling"}]
    first = plugin.fill_form(software=lineage)
    second = plugin.fill_form(software=lineage + lineage)
    ctrl = plugin.container.controlParameters
    assert first["entities"] == ["gamma"] and second["entities"] == []
    assert second["software"] == [] and second["wavelength"] is None
    assert len(ctrl.ENTITIES) == 1
    assert [str(r.NAME) for r in ctrl.SOFTWARE] == ["Aimless", "pointless", "CTRUNCATE"]
    segs = [(int(s.BEG), int(s.END), str(s.KIND)) for s in ctrl.ENTITIES[0].SEGMENTS]
    assert segs == [(1, 15, "tag"), (16, 134, "natural")]
    assert (plugin.getWorkDirectory() and
            os.path.exists(os.path.join(plugin.getWorkDirectory(), "input_params.xml")))


def test_fill_refuses_a_job_that_has_run(plugin):
    open(os.path.join(plugin.getWorkDirectory(), "params.xml"), "w").close()
    assert plugin.fill_form()["ok"] is False


def test_reference_search_offline_and_mocked():
    class Down:
        def post(self, *a, **k):
            raise OSError("no network")

    assert segments.find_references(GAMMA, session=Down())["ok"] is False

    class Response:
        def __init__(self, code, payload):
            self.status_code, self._payload = code, payload

        def json(self):
            return self._payload

    class Up:
        def post(self, url, json=None, timeout=None):
            return Response(200, {"result_set": [{"identifier": "1GYU_1", "services": [{"nodes": [
                {"match_context": [{"sequence_identity": 0.98, "query_length": 119,
                                    "query_aligned_seq": "A" * 119}]}]}]}]})

        def get(self, url, timeout=None):
            return Response(200, {"rcsb_polymer_entity_container_identifiers": {"uniprot_ids": ["O43747"]},
                                  "rcsb_polymer_entity": {"pdbx_description": "AP-1 complex subunit gamma-1"}})

    result = segments.find_references(GAMMA[15:], session=Up())
    assert result["ok"] and result["hits"][0]["uniprot"] == "O43747"
    assert result["hits"][0]["identity"] == 0.98 and result["hits"][0]["coverage"] == 1.0


# --- the writer and the gate ---------------------------------------------------

def test_written_values_are_the_pdb_spellings_and_pass_the_gate(tmp_path):
    doc, block, st = _prepared(tmp_path)
    spec = dict(BASE_SPEC, experiment={"source": "SYNCHROTRON", "site": "Diamond Light Source",
                                       "beamline": "i04-1", "detector": "PIXEL",
                                       "detector_type": "Pilatus 6M-F", "protocol": "SAD",
                                       "temperature": 100},
                crystal={"method": "vapor diffusion, sitting drop"},
                contact=dict(BASE_SPEC["contact"], country="united kingdom"),
                phasing={"method": "SAD"},
                software=[{"name": "phaser", "classification": "other"}])
    warnings = metadata.write_metadata(block, st, spec, fetch=_fetch)
    value = lambda tag: gemmi.cif.as_string(block.find_value(tag))
    assert value("_diffrn_source.pdbx_synchrotron_site") == "Diamond"
    assert value("_diffrn_source.type") == "DIAMOND BEAMLINE I04-1"
    assert value("_diffrn_detector.type") == "DECTRIS PILATUS 6M-F"
    assert value("_diffrn_radiation.pdbx_diffrn_protocol") == "SINGLE WAVELENGTH"
    assert value("_exptl_crystal_grow.method") == "VAPOR DIFFUSION, SITTING DROP"
    assert value("_pdbx_contact_author.country") == "United Kingdom"
    assert value("_refine.pdbx_method_to_determine_struct") == "SAD"
    assert _rows(block, "_software", ["name", "classification"])[0] == ("PHASER", "data processing")
    assert any("'other'" in w for w in warnings)
    assert enumerations.violations(block, metadata.WRITTEN_CATEGORIES) == []


def test_the_gate_names_a_value_onedep_rejects(tmp_path):
    doc, block, st = _prepared(tmp_path)
    metadata.write_metadata(block, st, dict(BASE_SPEC, crystal={"method": "COUNTER-DIFFUSION"}),
                            fetch=_fetch)
    assert ("_exptl_crystal_grow.method", "COUNTER-DIFFUSION") in \
        enumerations.violations(block, metadata.WRITTEN_CATEGORIES)


def test_refinement_flags_are_found_in_a_pipeline_shaped_job(tmp_path):
    """prosmart_refmac keeps refmac's settings below the top level (GammaSAD
    job 436): the reader searches the saved parameters, never the plugin."""
    from ccp4i2.lib.deposition.lineage import refinement_flags

    job = tmp_path / "job_436"
    sub = job / "job_436_1"
    sub.mkdir(parents=True)
    (job / "input_params.xml").write_text(
        "<ccp4i2><ccp4i2_body><container id='prosmartRefmac'><controlParameters>"
        "<USE_TWIN>False</USE_TWIN></controlParameters></container></ccp4i2_body></ccp4i2>")
    (sub / "params.xml").write_text(
        "<ccp4i2><ccp4i2_body><controlParameters><USEANOMALOUS>True</USEANOMALOUS>"
        "</controlParameters></ccp4i2_body></ccp4i2>")
    assert refinement_flags(job) == {"USE_TWIN": False, "USEANOMALOUS": True}
    (job / "input_params.xml").write_text("not xml")
    assert refinement_flags(job) == {"USEANOMALOUS": True}
    assert refinement_flags(tmp_path / "missing") == {}


# --- the Molecules tab's sequence picture ----------------------------------------

def test_sequence_analysis_gamma_gapped_and_overlapping():
    from ccp4i2.lib.deposition.segments import analyse_segments

    gamma = analyse_segments(GAMMA, [
        {"BEG": "1", "END": "15", "KIND": "tag"},
        {"BEG": "16", "END": "134", "KIND": "natural", "UNIPROT": "O43747"}])
    assert gamma["length"] == 134 and gamma["uncovered"] == [] and gamma["overlaps"] == []
    assert gamma["problems"] == []
    assert [(s["KIND"], s["BEG"], s["END"], s["length"], s["reference"]) for s in gamma["segments"]] == \
        [("tag", 1, 15, 15, "none"), ("natural", 16, 134, 119, "O43747")]
    assert gamma["segments"][0]["first"] == "MHHHH" and gamma["segments"][1]["last"] == "PQSWQ"

    gapped = analyse_segments(GAMMA, [{"BEG": 1, "END": 15, "KIND": "tag"},
                                      {"BEG": 21, "END": 140, "KIND": "natural"}])
    assert gapped["uncovered"] == [(16, 20)]
    joined = " | ".join(gapped["problems"])
    assert "past the sequence end" in joined and "needs a UniProt accession" in joined
    assert "covered by no segment: 16-20" in joined

    overlapping = analyse_segments(GAMMA, [{"BEG": 1, "END": 20, "KIND": "tag"},
                                           {"BEG": 16, "END": 134, "KIND": "designed"}])
    assert overlapping["overlaps"] == [(16, 20, [0, 1])]
    assert overlapping["uncovered"] == []
    assert any("16-20 are in segments 1 and 2" in p for p in overlapping["problems"])


def test_sequence_view_reads_the_saved_records_and_never_writes(plugin, monkeypatch, tmp_path):
    ctrl = _filled(plugin)
    ent = ctrl.ENTITIES.addItem()
    ent.SEQUENCE_NAME.set("gamma")
    for beg, end, kind in ((1, 15, "tag"), (16, 134, "natural")):
        seg = ent.SEGMENTS.addItem()
        seg.BEG.set(beg)
        seg.END.set(end)
        seg.KIND.set(kind)

    class Seq:
        def __init__(self, name, sequence):
            self.name, self.sequence, self.polymerType = name, sequence, "PROTEIN"

    class Asu:
        def isSet(self):
            return True

        def loadFile(self):
            pass

        class fileContent:
            seqList = [Seq("gamma", GAMMA)]

    monkeypatch.setattr(plugin.container.inputData, "ASUCONTENT", Asu(), raising=False)
    view = plugin.sequence_view("gamma")
    assert view["ok"] and view["sequence"] == GAMMA
    assert [(s["KIND"], s["BEG"], s["END"]) for s in view["segments"]] == [("tag", 1, 15), ("natural", 16, 134)]
    assert view["uncovered"] == [] and view["modelled"] is None
    assert view["problems"] == ["segment 2 (16-134 natural): a natural segment needs a UniProt accession"]
    missing = plugin.sequence_view("beta")
    assert missing["ok"] is False and "no sequence found for beta" in missing["error"]
    assert not os.path.exists(os.path.join(str(tmp_path), "input_params.xml"))


def test_modelled_positions_mark_the_unmodelled_tag(tmp_path):
    from ccp4i2.lib.deposition.segments import modelled_positions
    from ccp4i2.tests.unit.plugins.test_prepare_deposition_metadata import _toy_model

    construct = GAMMA[:40]
    path = _toy_model(tmp_path / "m.cif", construct, start=16, skip=(30,))
    covered = modelled_positions(construct, path)
    assert covered[0] == 16 and 30 not in covered and covered[-1] == 40
    assert modelled_positions(construct, tmp_path / "missing.cif") is None
