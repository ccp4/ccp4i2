"""Parsing a fileUse reference: the rule, with no database in sight.

parse_file_use is pure on purpose (the models import is function-local), so the
syntax can be tested wherever CI runs.
"""
import pytest

from ccp4i2.lib.utils.files.file_use import (
    FileUseError,
    _and_these_exist,
    _did_you_mean,
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


class TestTheSuggestions:
    """A reference that does not resolve must say why and offer the near miss.

    Printing i2run's usage instead is no help: servalcat_pipe has 219
    arguments, so the answer would be buried. The relevant list is always short
    -- the parameters of one job, or the registered task names.
    """

    def test_a_near_miss_is_offered(self):
        assert _did_you_mean("freer_flag", ["freerflag", "refmac"]) == (
            " Did you mean 'freerflag'?"
        )

    def test_wrong_case_alone_is_offered(self):
        """difflib is strict about case and about short strings, and case is
        the typo people actually make."""
        assert _did_you_mean("freerout", ["FREEROUT", "F_SIGF"]) == (
            " Did you mean 'FREEROUT'?"
        )

    def test_an_exact_hit_suggests_nothing(self):
        assert _did_you_mean("FREEROUT", ["FREEROUT"]) == ""

    def test_nothing_remotely_close_suggests_nothing(self):
        assert _did_you_mean("QQQQQQ", ["FREEROUT", "F_SIGF"]) == ""

    def test_no_candidates_suggests_nothing(self):
        assert _did_you_mean("anything", []) == ""

    def test_the_inventory_is_sorted_and_deduplicated(self):
        assert _and_these_exist("It has", ["B", "A", "B", None]) == (
            " It has: A, B."
        )

    def test_a_long_inventory_is_capped(self):
        text = _and_these_exist("It has", [f"P{i}" for i in range(20)], limit=3)
        assert text.startswith(" It has: P0, P1, P10")
        assert "and 17 more" in text

    def test_an_empty_inventory_says_nothing(self):
        assert _and_these_exist("It has", []) == ""


def test_the_messages_are_ascii_only():
    """These strings reach print() in the i2run management command's error
    path. On Windows that console is cp1252, and a UnicodeEncodeError raised
    from inside a try/except is the documented way to lose a job silently
    (CLAUDE.md). An em dash slipped in here and was caught this way."""
    import inspect

    from ccp4i2.lib.utils.files import file_use

    source = inspect.getsource(file_use)
    offenders = sorted({char for char in source if ord(char) > 127})
    assert not offenders, f"non-ASCII in file_use.py: {offenders}"
