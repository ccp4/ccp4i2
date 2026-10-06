"""UniProt from a name: how typed text is read, how candidates rank, and the
sequence cut to a construct. The network is faked: no test talks to UniProt.

An agent asked to solve "CDK4/cyclin D" had names and no sequences; nothing
in CCP4i2 turned one into the other.
"""
import io
import json
import urllib.parse

import pytest

from ccp4i2.lib.utils.sequences import uniprot


@pytest.mark.parametrize("typed, kind, term, taxid", [
    ("Human CDK2", "gene", "CDK2", 9606),
    ("human-CDK2", "gene", "CDK2", 9606),
    ("CDK2 from Human", "gene", "CDK2", 9606),
    ("CDK2 (human)", "gene", "CDK2", 9606),
    ("CDK2, Homo sapiens", "gene", "CDK2", 9606),
    ("CDK2 human", "gene", "CDK2", 9606),
    ("human cyclin dependent kinase 2", "protein_name", "cyclin dependent kinase 2", 9606),
    ("CDK4 mouse", "gene", "CDK4", 10090),
    ("E. coli DnaK", "gene", "DnaK", 83333),
    ("P24941", "accession", "P24941", None),
    ("CDK2_HUMAN", "entry_name", "CDK2_HUMAN", None),
    ("hCDK2", "gene", "hCDK2", None),   # a one-letter prefix is too ambiguous to guess
    ("cyclin D1", "protein_name", "cyclin D1", None),
])
def test_typed_text_is_read(typed, kind, term, taxid):
    reading = uniprot.read_query(typed)
    assert (reading["kind"], reading["term"]) == (kind, term)
    assert (reading["organism"] or {}).get("taxid") == taxid


def test_an_organism_given_wins_over_one_in_the_text():
    assert uniprot.read_query("CDK2 human", organism="mouse")["organism"]["taxid"] == 10090


def _entry(acc, entry_name, name, gene, taxid, reviewed=True, sequence="MENFQKVEKIGEGTYGVVYKA"):
    return {"primaryAccession": acc, "uniProtkbId": entry_name,
            "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
            "organism": {"scientificName": "Homo sapiens" if taxid == 9606 else "Mus musculus",
                         "taxonId": taxid},
            "proteinDescription": {"recommendedName": {"fullName": {"value": name}}},
            "genes": [{"geneName": {"value": gene}}],
            "sequence": {"value": sequence, "length": len(sequence)}}


class _Opener:
    """Answers UniProt URLs from a table; records what was asked."""

    def __init__(self, search=(), entries=None):
        self.search, self.entries, self.asked = list(search), entries or {}, []

    def __call__(self, request, timeout=None):
        url = request.full_url
        self.asked.append(urllib.parse.unquote_plus(url))
        if "/search?" in url:
            body = {"results": self.search}
        else:
            key = url.rsplit("/", 1)[1].replace(".json", "")
            body = self.entries[key]
        return io.BytesIO(json.dumps(body).encode())


def test_candidates_rank_reviewed_exact_names_first_and_none_is_chosen():
    opener = _Opener(search=[
        _entry("Q9BW66", "CINP_HUMAN", "Cyclin-dependent kinase 2-interacting protein", "CINP", 9606),
        _entry("B0X", "B0X_HUMAN", "Cyclin-dependent kinase 2", "CDK2", 9606, reviewed=False),
        _entry("P24941", "CDK2_HUMAN", "Cyclin-dependent kinase 2", "CDK2", 9606),
    ])
    out = uniprot.search("human cyclin dependent kinase 2", opener=opener)
    assert [c["accession"] for c in out["candidates"]] == ["P24941", "B0X", "Q9BW66"]
    assert out["read_as"]["organism"]["taxid"] == 9606
    assert any("organism_id:9606" in url for url in opener.asked)


def test_a_family_name_finds_the_family_members_first():
    opener = _Opener(search=[
        _entry("Q9Y222", "DMTF1_HUMAN", "Cyclin-D-binding Myb-like transcription factor 1", "DMTF1", 9606),
        _entry("P24385", "CCND1_HUMAN", "G1/S-specific cyclin-D1", "CCND1", 9606),
        _entry("P30279", "CCND2_HUMAN", "G1/S-specific cyclin-D2", "CCND2", 9606),
    ])
    out = uniprot.search("cyclin D human", opener=opener)
    assert [c["gene"] for c in out["candidates"][:2]] == ["CCND1", "CCND2"]
    assert any("cyclin AND D*" in url for url in opener.asked)  # the prefix query


def test_an_accession_is_looked_up_directly():
    opener = _Opener(entries={"P24941": _entry("P24941", "CDK2_HUMAN", "Cyclin-dependent kinase 2",
                                               "CDK2", 9606)})
    out = uniprot.search("P24941", opener=opener)
    assert out["candidates"][0]["entry_name"] == "CDK2_HUMAN" and len(opener.asked) == 1


def test_fetch_cuts_the_construct_and_records_where_it_came_from():
    seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(300))
    opener = _Opener(entries={"P14635": _entry("P14635", "CCNB1_HUMAN", "G2/mitotic-specific cyclin-B1",
                                               "CCNB1", 9606, sequence=seq)})
    out = uniprot.fetch("P14635", "175-200", opener=opener)
    assert out["sequence"] == seq[174:200] and out["range"] == "175-200"
    assert out["fasta"].startswith(">sp|P14635|CCNB1_HUMAN G2/mitotic-specific cyclin-B1 "
                                   "OS=Homo sapiens residues 175-200\n")
    with pytest.raises(uniprot.UniProtError, match="has 300 residues"):
        uniprot.fetch("P14635", "175-400", opener=opener)


@pytest.mark.parametrize("bad", ["175", "432-175", "0-10", "a-b"])
def test_a_bad_range_is_refused(bad):
    with pytest.raises(uniprot.UniProtError):
        uniprot.parse_range(bad)


def test_a_name_is_not_fetched_and_being_offline_is_said():
    with pytest.raises(uniprot.UniProtError, match="search by name first"):
        uniprot.fetch("cyclin B1")

    def offline(request, timeout=None):
        raise OSError("no network")
    with pytest.raises(uniprot.UniProtError, match="could not be reached"):
        uniprot.search("CDK2", opener=offline)
