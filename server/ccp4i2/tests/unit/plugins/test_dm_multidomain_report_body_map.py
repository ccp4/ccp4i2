"""The dm_multidomain report must say WHAT was averaged, not just how it went.

Per-domain correlations numbered 1, 2, 3 leave a reader to remember which
domain was which. The body map draws the partition on the reference copy's
residues -- the same picture the task interface draws while the job is being
set up -- so a job can be judged months later without reopening its
parameters.

No binary and no Report machinery: the renderer reads self.xmlnode and appends
HTML, so a stub parent is enough to test it.
"""
import xml.etree.ElementTree as ET

from ccp4i2.wrappers.dm_multidomain.script.dm_multidomain_report import (
    dm_multidomain_report,
)


class _Parent:
    def __init__(self):
        self.html = []

    def append(self, text):
        self.html.append(text)


def _render(xml_text):
    """Render just the body map. Report.__init__ wants a job on disk, and the
    renderer needs nothing but xmlnode and a parent to append to."""
    report = object.__new__(dm_multidomain_report)
    report.xmlnode = ET.fromstring(xml_text)
    parent = _Parent()
    report._add_body_map(parent)
    return "".join(parent.html)


_TWO_BODIES = """
<DmMultidomainResult>
  <BodyMap>
    <Track role="CDK" chain="A" first="1" last="298"/>
    <Track role="cyclin" chain="B" first="175" last="432"/>
    <Body number="1" mode="average" name="CDK_1_298">
      <Segment role="CDK" lo="1" hi="298"/>
      <Fit copy="C+D" rmsd="0.84"/>
    </Body>
    <Body number="2" mode="refine" name="cyclin_175_432">
      <Segment role="cyclin" lo="175" hi="432"/>
      <Fit copy="C+D" rmsd="0.58"/>
    </Body>
  </BodyMap>
</DmMultidomainResult>
"""


def test_draws_a_track_per_entity_and_a_block_per_body():
    html = _render(_TWO_BODIES)
    assert "<svg" in html and "</svg>" in html
    assert "CDK (A)" in html and "cyclin (B)" in html
    assert "body 1: 1-298" in html
    assert "body 2: 175-432" in html


def test_lists_the_fit_of_each_body_against_each_copy():
    html = _render(_TWO_BODIES)
    assert "C+D 0.84" in html and "C+D 0.58" in html
    assert "refine" in html


def test_overlapping_bodies_are_drawn_as_a_clash():
    """Two bodies claiming the same residues is the mistake the picture
    exists to make visible."""
    html = _render("""
    <DmMultidomainResult>
      <BodyMap>
        <Track role="_" chain="A" first="1" last="485"/>
        <Body number="1" mode="average" name="a"><Segment role="_" lo="1" hi="300"/></Body>
        <Body number="2" mode="average" name="b"><Segment role="_" lo="200" hi="485"/></Body>
      </BodyMap>
    </DmMultidomainResult>
    """)
    assert "url(#dmclash)" in html
    assert "chain A" in html        # the implicit role is shown as its chain


def test_an_excluded_body_is_drawn_faint_not_omitted():
    html = _render("""
    <DmMultidomainResult>
      <BodyMap>
        <Track role="_" chain="A" first="1" last="485"/>
        <Body number="1" mode="exclude" name="a"><Segment role="_" lo="1" hi="139"/></Body>
      </BodyMap>
    </DmMultidomainResult>
    """)
    assert "fill-opacity='0.35'" in html
    assert "exclude" in html


def test_a_job_without_a_body_map_renders_nothing_rather_than_failing():
    """Older jobs have no BodyMap in their program.xml."""
    assert _render("<DmMultidomainResult/>") == ""
