from .utils import demoData, i2run


def test_paste_pir_alignment():
    """Test ProvideAlignment with pasted PIR-format alignment."""
    pir_text = (
        ">P1;target\n"
        "target sequence\n"
        "MKYLLPTAAAGLLLLAAQPAMA*\n"
        ">P1;template\n"
        "template sequence\n"
        "MKYLLPTAAAGLLLLAAQPAMA*\n"
    )
    args = ["ProvideAlignment"]
    args += ["--PASTEORREAD", "PASTE"]
    args += ["--SEQUENCETEXT", pir_text]
    with i2run(args) as job:
        alignout = job / "ALIGNMENTFILE.aln"
        assert alignout.exists(), f"No alignment file: {list(job.iterdir())}"


def test_read_pir_file():
    """Test ProvideAlignment by reading gamma PIR file."""
    args = ["ProvideAlignment"]
    args += ["--PASTEORREAD", "ALIGNIN"]
    args += ["--ALIGNIN", demoData("gamma", "gamma.pir")]
    with i2run(args) as job:
        alignout = job / "ALIGNMENTFILE.aln"
        assert alignout.exists(), f"No alignment file: {list(job.iterdir())}"
