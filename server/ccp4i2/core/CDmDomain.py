"""CDmDomain - a rigid body (dm "domain") for multi-domain NCS averaging.

A rigid body is a set of SEGMENTS that move together, plus an averaging
``mode``. Each segment is a ``role:firstRes-lastRes`` range; a token with no
``role:`` prefix uses the single implicit role, so the homomer case stays
terse (``segments = "340-485"``). Cross-chain bodies -- the CDK C-helix that
travels with the cyclin N-lobe -- are expressed by mixing roles in one
``segments`` string (``"cyclin:10-95,CDK:45-60"``). The role->chain mapping for
each NCS copy lives in the task's ASSEMBLY parameter, not here.

The segment string is parsed by the wrapper (dm_ncs_lib.parse_segments); core
deliberately holds only the data, so it stays free of any wrapper/binary
dependency. Following the codebase convention: the *Stub* class defines the data
model (content fields, via the @cdata_class ``attributes=``), and the derived
class carries methods.

Resolvable by the def.xml class-name lookup via ``ccp4i2.core.CDmDomain``.
"""
from ccp4i2.core.base_object.class_metadata import (
    cdata_class, attribute, AttributeType)
from ccp4i2.core.cdata_stubs.CCP4Data import CData


@cdata_class(
    attributes={
        "segments": attribute(AttributeType.STRING),
        "mode": attribute(AttributeType.STRING),
    },
    contents_order=['segments', 'mode'],
    content_qualifiers={
        "segments": {
            'default': '',
            'toolTip': ("Comma-separated role:first-last residue ranges that "
                        "move as one rigid body, e.g. '340-485' or "
                        "'cyclin:10-95,CDK:45-60'"),
        },
        "mode": {
            # enumerators/menuText must be LISTS at runtime (the def.xml parser
            # splits "<enumerators>a,b,c</enumerators>" into a list; the
            # frontend autocomplete calls .map on them). Passing comma-strings
            # here reaches the UI as a string -> "enumerators.map is not a
            # function".
            'onlyEnumerators': True,
            'enumerators': ['average', 'refine', 'exclude'],
            'menuText': ['average', 'refine', 'exclude'],
            'default': 'average',
            'toolTip': 'NCS averaging treatment for this rigid body',
        },
    },
)
class CDmDomainStub(CData):
    """Data model: a rigid body's segment string plus an NCS averaging mode."""
    pass


class CDmDomain(CDmDomainStub):
    """A rigid body for multi-domain NCS averaging: one or more residue-range
    segments (possibly spanning several chains via roles) plus how it should be
    averaged (average / refine / exclude).
    """

    def averaging_mode(self):
        m = str(self.mode) if self.mode.isSet() else 'average'
        return m if m in ('average', 'refine', 'exclude') else 'average'

    def segments_spec(self):
        """The raw 'role:lo-hi,...' segment string (empty if unset). Parsing
        into (role, lo, hi) tuples is the wrapper's job (dm_ncs_lib)."""
        return str(self.segments) if self.segments.isSet() else ''
