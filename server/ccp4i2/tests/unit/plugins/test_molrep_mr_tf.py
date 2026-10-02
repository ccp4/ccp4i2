"""molrep_mr and molrep_den put MOLREP's own result in program.xml.

They wrote ``n_solution 1`` and ``mr_score 0.0000`` whatever MOLREP found,
so a reader (an agent judging the job, or a person) could not tell a good
search from a failed one; MOLREP's real score and z-score were only in its
molrep.xml. Now that is where they are read from.
"""
import pytest

pytest.importorskip("lxml")

from ccp4i2.wrappers.molrep_mr.script.molrep_mr import mr_tf_element  # noqa: E402

# As MOLREP 11.9.02 writes it (Gamma, docs scenario): mixed content, values padded.
MOLREP_XML = """<?xml version="1.0" encoding="ASCII" standalone="yes"?>
<!-- MR_TF output  -->
<MR_TF>
Error <err_level>         0</err_level>
Message <err_message>normal termination</err_message>
Job <job>MR_TF</job>
vers <vers>11.9.02</vers>
nmon_solution <n_solution>         1</n_solution>
mr_resmin <mr_resmin>   26.6370</mr_resmin>
mr_resmax <mr_resmax>    1.9006</mr_resmax>
mr_score <mr_score>    0.7854</mr_score>
mr_score_previous <mr_score_previous>    0.0000</mr_score_previous>
mr_zscore <mr_zscore>   18.9449</mr_zscore>
mr_zscore_previous <mr_zscore_previous>    0.0000</mr_zscore_previous>
<solution sol_file="/x/molrep.pdb"/>
<found dimer_file="none"/>
Time_elapsed <tel>     0h  0m  5s</tel>
Time_remained <trm>null</trm>
</MR_TF>
"""


def test_molreps_result_is_copied(tmp_path):
    path = tmp_path / "molrep.xml"
    path.write_text(MOLREP_XML)
    tf = mr_tf_element(path)
    assert tf.tag == "MR_TF"
    assert {child.tag: child.text for child in tf} == {
        "err_level": "0", "err_message": "normal termination", "n_solution": "1",
        "mr_resmin": "26.6370", "mr_resmax": "1.9006",
        "mr_score": "0.7854", "mr_zscore": "18.9449",
    }


@pytest.mark.parametrize("content", [None, "", "<MR_TF><mr_score>", "<MR_TF/>"])
def test_nothing_is_invented_when_molrep_wrote_nothing(tmp_path, content):
    path = tmp_path / "molrep.xml"
    if content is not None:
        path.write_text(content)
    tf = mr_tf_element(path)
    assert tf.tag == "MR_TF" and len(tf) == 0
