"""Parsing a fileUse reference: the rule, with no database in sight.

parse_file_use is pure on purpose (the models import is function-local), so the
syntax can be tested wherever CI runs.
"""
import pytest

from ccp4i2.lib.utils.files.file_use import (
    FileUseError,
    parse_file_use,
)


class TestTheForms:
    def test_a_bare_positive_index_is_a_job_number(self):
        ref = parse_file_use("3.XYZOUT")
        assert ref.job_number == "3"
        assert ref.index is None
        assert ref.param_name == "XYZOUT"
        assert ref.task_name is None

    def test_brackets_mean_the_same_thing(self):
        assert parse_file_use("[3].XYZOUT").job_number == "3"

    def test_a_sub_job_can_be_named_outright(self):
        """Job.number is a CharField holding "1" or "1.1", and inside brackets
        a dotted number is unambiguous."""
        assert parse_file_use("[3.1].XYZOUT").job_number == "3.1"

    def test_a_bracketed_negative_is_relative(self):
        ref = parse_file_use("[-1].XYZOUT")
        assert ref.index == -1
        assert ref.job_number is None

    def test_a_task_name_scopes_the_index(self):
        ref = parse_file_use("prosmart_refmac[-1].XYZOUT")
        assert ref.task_name == "prosmart_refmac"
        assert ref.index == -1

    def test_a_parameter_index_can_be_given(self):
        ref = parse_file_use("prosmart_refmac[0].XYZOUT[1]")
        assert (ref.task_name, ref.index, ref.param_index) == (
            "prosmart_refmac",
            0,
            1,
        )

    def test_the_default_parameter_index_is_the_last(self):
        assert parse_file_use("[3].XYZOUT").param_index == -1


class TestWhatItRefuses:
    def test_a_bare_negative_says_to_use_brackets(self):
        """The reason brackets are canonical: argparse rejects a bare negative
        outright ("expected at least one argument"), because a token starting
        with '-' is an option string and its numeric escape hatch only exempts
        '-1' and '-1.5'. The old parser accepted this form happily, which meant
        the code looked like it supported something unreachable."""
        with pytest.raises(FileUseError) as caught:
            parse_file_use("-1.XYZOUT")
        message = str(caught.value)
        assert "[-1].XYZOUT" in message, message
        assert "argparse" in message, message

    def test_nonsense_says_what_was_expected(self):
        with pytest.raises(FileUseError) as caught:
            parse_file_use("rubbish")
        assert "[jobNumber].PARAM" in str(caught.value)

    def test_an_empty_reference_is_refused(self):
        with pytest.raises(FileUseError):
            parse_file_use("")

    @pytest.mark.parametrize(
        "text", ["[].X", "[1].", "[1.].X", "[a].X", "1..X", "[1][2].X"]
    )
    def test_malformed_indices_raise_the_documented_error(self, text):
        """Never a bare ValueError from int(): the old parser only caught
        AttributeError, so '[3].XYZOUT' escaped as an unhandled
        'invalid literal for int()'."""
        with pytest.raises(FileUseError):
            parse_file_use(text)
