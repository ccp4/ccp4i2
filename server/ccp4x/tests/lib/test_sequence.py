import pytest
from ...lib.sequence import code1, SequenceType, L_PEPTIDE, D_PEPTIDE, DNA, RNA


@pytest.mark.parametrize(
    ("resname", "expected"),
    [
        ("", None),
        ("A", "A"),
        ("A1LU6", None),
        ("ALA", "A"),
        ("ASX", "B"),
        ("DA", "A"),
        ("DAL", "A"),
        ("DN", "X"),
        ("DU", "U"),
        ("GLX", "Z"),
        ("HOH", None),
        ("I", "I"),
        ("N", "X"),
        ("PYL", "O"),
        ("SEC", "U"),
        ("UNK", "X"),
        ("XXX", None),
    ],
)
def test_code1(resname: str, expected: str | None):
    assert code1(resname) == expected


@pytest.mark.parametrize(
    ("sequence_type", "sequence", "expected"),
    [
        (L_PEPTIDE, "AA", ["ALA", "ALA"]),
        (L_PEPTIDE, "AA*", ["ALA", "ALA"]),
        (L_PEPTIDE, " A A ", ["ALA", "ALA"]),
        (L_PEPTIDE, "-A-A-", ["ALA", "ALA"]),
        (L_PEPTIDE, "\nAA\n", ["ALA", "ALA"]),
        (L_PEPTIDE, "A(A)A", ["ALA", "A", "ALA"]),
        (L_PEPTIDE, "A(ALA)A", ["ALA", "ALA", "ALA"]),
        (L_PEPTIDE, "BOUXZ", ["ASX", "PYL", "SEC", "UNK", "GLX"]),
        (D_PEPTIDE, "ACDEF", ["DAL", "DCY", "DAS", "DGL", "DPN"]),
        (DNA, "ACGITUX", ["DA", "DC", "DG", "DI", "DT", "DU", "DN"]),
        (RNA, "ACGIUX", ["A", "C", "G", "I", "U", "N"]),
    ],
)
def test_parse(sequence_type: SequenceType, sequence: str, expected: list[str]):
    assert sequence_type.parse(sequence) == expected


@pytest.mark.parametrize(
    ("sequence_type", "sequence"),
    [
        (L_PEPTIDE, "J"),
        (L_PEPTIDE, "A/A"),
        (L_PEPTIDE, "A,A"),
        (L_PEPTIDE, "A()A"),
        (D_PEPTIDE, "X"),
        (DNA, "J"),
        (RNA, "T"),
    ],
)
def test_parse_error(sequence_type: SequenceType, sequence: str):
    with pytest.raises(ValueError):
        sequence_type.parse(sequence)
