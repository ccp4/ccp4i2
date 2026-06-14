"""CDmDomain - a multi-domain NCS averaging domain.

Composed from the existing CResidueRange data model (chainId / firstRes /
lastRes) plus an averaging ``mode``. Following the codebase convention:
the *Stub* class defines the data model (content fields, via the @cdata_class
``attributes=`` merged across the MRO), and the derived class carries methods.

Resolvable by the def.xml class-name lookup via ``ccp4i2.core.CDmDomain``.
"""
from ccp4i2.core.base_object.class_metadata import (
    cdata_class, attribute, AttributeType)
from ccp4i2.core.cdata_stubs.CCP4ModelData import CResidueRangeStub


@cdata_class(
    attributes={
        # chainId / firstRes / lastRes are inherited from CResidueRangeStub's
        # metadata (merged across the MRO); we only add the averaging mode.
        "mode": attribute(AttributeType.STRING),
    },
    contents_order=['chainId', 'firstRes', 'lastRes', 'mode'],
    content_qualifiers={
        "chainId": {'default': ''},
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
            'toolTip': 'NCS averaging treatment for this domain',
        },
    },
)
class CDmDomainStub(CResidueRangeStub):
    """Data model: a residue range plus an NCS averaging mode."""
    pass


class CDmDomain(CDmDomainStub):
    """A domain for multi-domain NCS averaging: a residue range over the
    reference copy plus how it should be averaged (average / refine / exclude).
    """

    def averaging_mode(self):
        m = str(self.mode) if self.mode.isSet() else 'average'
        return m if m in ('average', 'refine', 'exclude') else 'average'

    def residue_bounds(self):
        """(lo, hi) integers for the residue range, or None if not set."""
        if not (self.firstRes.isSet() and self.lastRes.isSet()):
            return None
        return int(str(self.firstRes)), int(str(self.lastRes))
