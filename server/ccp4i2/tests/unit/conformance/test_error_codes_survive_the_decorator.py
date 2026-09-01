"""A class's own ERROR_CODES must survive the decorator that merges them.

`cdata_class` merges error codes along the MRO so a subclass can declare only
what it adds. It does that by assigning `cls.ERROR_CODES`, and a decorator runs
*after* the class body --- so a class that also writes `ERROR_CODES = {...}` by
hand had it discarded.

The discard is silent. `ERROR_CODES` is still populated, so nothing looks
wrong until a code is looked up, and `validity()` reads codes by key:

    details=self.ERROR_CODES[201]['description']    ->  KeyError: 201

CRangeSelection is the case that found it, and it is a nasty one: the merged
declaration carries "201" as a *string* while the hand-written body carries
201 as an *int*, so even a len() of the dict looks healthy.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


@pytest.fixture(scope="module", autouse=True)
def _django():
    import os
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    django.setup()


def test_a_hand_written_error_code_is_not_discarded():
    from ccp4i2.core.base_object.cdata import CData
    from ccp4i2.core.base_object.class_metadata import cdata_class

    @cdata_class(error_codes={"201": {"description": "from the decorator"}})
    class Declared(CData):
        ERROR_CODES = {201: {"description": "written by hand"}}

    assert 201 in Declared.ERROR_CODES, \
        f"the class's own code was discarded: {sorted(map(str, Declared.ERROR_CODES))}"
    assert Declared.ERROR_CODES[201]["description"] == "written by hand"


def test_the_class_that_found_it_still_validates():
    """CRangeSelection.validity() reads ERROR_CODES[201] by integer key."""
    from ccp4i2.core.CCP4Data import CRangeSelection

    bad = CRangeSelection()
    bad.set("not-a-range")
    assert list(bad.validity().entries()), "a malformed range was accepted"


def test_inherited_codes_are_still_merged():
    """The body must win without switching inheritance back off."""
    from ccp4i2.core.base_object.cdata import CData
    from ccp4i2.core.base_object.class_metadata import cdata_class

    @cdata_class(error_codes={1: {"description": "base"}})
    class Base(CData): pass

    @cdata_class(error_codes={2: {"description": "child"}})
    class Child(Base):
        ERROR_CODES = {3: {"description": "by hand"}}

    assert {1, 2, 3} <= set(Child.ERROR_CODES)
