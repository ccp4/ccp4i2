"""The pandda_campaign report renders from whatever program.xml holds.

Its stats arrive as XML text, so a value used in arithmetic or ``:.2f``
formatting has to be turned into a number first. It was not: the resolution
warning did ``best_input + 0.5`` on a string, and the whole Analysis section
became a REPORT_GENERATION_FAILED the moment a run had both resolutions --
which, now that a dispatched run refreshes its report mid-run, is every
campaign with a resolution spread.
"""
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_report import (
    pandda_campaign_report as Report,
)

DRAGGED = """
<pandda_campaign>
  <state>dispatched</state>
  <analysis><stats>
    <n_analysed>19</n_analysed>
    <median_resolution>3.10</median_resolution>
    <best_input_resolution>1.98</best_input_resolution>
    <n_input_better_than_processing>140</n_input_better_than_processing>
    <n_events>0</n_events>
  </stats></analysis>
</pandda_campaign>
"""


def _render(xml: str):
    return Report(xmlnode=ET.fromstring(xml), jobInfo={}, jobStatus="Running remotely")


def test_a_run_dragged_to_a_worse_resolution_renders_its_warning():
    """The case that raised: both resolutions present and datasets dragged."""
    assert _render(DRAGGED) is not None


@pytest.mark.parametrize("stats", [
    "<median_resolution>3.10</median_resolution><best_input_resolution>1.98</best_input_resolution>",
    "<median_resolution>2.00</median_resolution><best_input_resolution>1.98</best_input_resolution>",
    "<median_resolution></median_resolution><best_input_resolution>1.98</best_input_resolution>",
    "<median_resolution>None</median_resolution><best_input_resolution>1.98</best_input_resolution>",
    "<best_input_resolution>1.98</best_input_resolution>",
    "",
])
def test_partial_or_odd_stats_still_render(stats):
    """Mid-run the tables are partial, so a stat can be absent, empty or the
    string 'None'. None of those may stop the page rendering."""
    assert _render(f"<pandda_campaign><state>dispatched</state><analysis><stats>"
                   f"<n_input_better_than_processing>140</n_input_better_than_processing>"
                   f"{stats}</stats></analysis></pandda_campaign>") is not None


def test_a_report_with_no_analysis_at_all_renders():
    assert _render("<pandda_campaign><state>dispatched</state></pandda_campaign>") is not None


def test_the_number_helper_is_total():
    assert Report._n("2.5") == 2.5
    assert Report._n("") is None
    assert Report._n(None) is None
    assert Report._n("None") is None
    assert Report._n("n/a") is None


MID_RUN = """
<pandda_campaign>
  <state>dispatched</state>
  <analysis><stats>
    <events_table>False</events_table>
    <n_analysed>35</n_analysed>
    <n_events>26</n_events>
    <n_datasets_with_events>13</n_datasets_with_events>
  </stats></analysis>
</pandda_campaign>
"""

FINISHED_EMPTY = """
<pandda_campaign>
  <state>finished</state>
  <analysis><stats>
    <events_table>True</events_table>
    <n_analysed>35</n_analysed>
    <n_events>0</n_events>
  </stats></analysis>
</pandda_campaign>
"""


def _texts(report):
    return " ".join(e.text or "" for e in report.as_data_etree().iter()
                    if (e.text or "").strip())


def test_a_run_still_going_does_not_claim_no_events_were_found():
    """The events table is written at the end, so mid-run `events` is empty
    however many the datasets have found. Saying PanDDA found none then sat
    above a per-dataset column plainly showing several."""
    text = _texts(_render(MID_RUN))
    assert "found no events" not in text
    assert "still going" in text
    assert "26 event(s) so far" in text


def test_a_finished_run_with_nothing_still_says_so():
    text = _texts(_render(FINISHED_EMPTY))
    assert "found no events" in text
