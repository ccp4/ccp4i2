"""A stub's annotations are its field declaration.

Every generated stub declares its fields twice:

    @cdata_class(attributes={"structure": attribute(CUSTOM, "CPdbDataFileStub"), ...})
    class CPdbEnsembleItemStub(CData):
        structure: Optional[CPdbDataFileStub] = None

Only the decorator was load-bearing; the annotation was decoration, and the
`= None` is what the runtime walk overwrites. But the annotation is the
declaration a dataclass uses, and it carries strictly more --- the class
itself, rather than its name as a string to be looked up later.

Measured across the stubs the two agree exactly: 621 fields; for the 278
naming a custom class the annotation names the same one, and the other 343
map one-to-one (STRING->CString, INT->CInt, FLOAT->CFloat, BOOLEAN->CBoolean).

The runtime now reads the annotations, falling back to the decorator for
hand-written classes such as CDataFile which carry one but no annotations.
That makes the annotation the single declaration, and the duplicate in the
decorator removable --- which is what a dataclass conversion needs, and why
this comes before it.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


def _stub_classes():
    import ccp4i2.core.CCP4ModelData          # noqa: F401  (registers subclasses)
    import ccp4i2.core.CCP4XtalData           # noqa: F401
    import ccp4i2.core.CCP4PerformanceData    # noqa: F401
    import ccp4i2.core.CCP4Data               # noqa: F401
    from ccp4i2.core.base_object.cdata import CData

    def walk(cls):
        for sub in cls.__subclasses__():
            yield sub
            yield from walk(sub)

    return [k for k in set(walk(CData))
            if "_metadata" in vars(k) and vars(k)["_metadata"].attributes]


def test_annotations_reproduce_the_decorator():
    """The two declarations must agree, or reading one changes behaviour."""
    from ccp4i2.core.base_object.class_metadata import attributes_from_annotations

    checked, disagreeing = 0, []
    for klass in _stub_classes():
        declared = vars(klass)["_metadata"].attributes
        from_annotations = attributes_from_annotations(klass)
        if not from_annotations:
            continue          # hand-written, no annotations: the fallback case
        checked += 1
        a = {n: (d.attr_type, d.custom_class) for n, d in declared.items()}
        b = {n: (d.attr_type, d.custom_class) for n, d in from_annotations.items()}
        if a != b:
            disagreeing.append(
                f"{klass.__name__}: {[n for n in a if b.get(n) != a[n]][:3]}")

    assert checked > 100, f"only {checked} stubs checked --- proves little"
    assert not disagreeing, (
        f"{len(disagreeing)} classes where the annotation and the decorator "
        f"declare different fields: {disagreeing[:5]}"
    )


def test_a_class_with_no_annotations_still_works():
    """CDataFile is hand-written: decorator, no annotations. The fallback."""
    from ccp4i2.core.base_object.cdata_file import CDataFile

    obj = CDataFile()
    names = {c._name for c in obj.children()}
    assert {"baseName", "relPath", "project", "annotation"} <= names, \
        f"the decorator fallback did not build the fields: {sorted(names)}"


def test_the_declaration_actually_builds_the_object():
    """The end-to-end statement: annotations produce the right children."""
    from ccp4i2.core.CCP4ModelData import CPdbEnsembleItem

    item = CPdbEnsembleItem()
    names = [c._name for c in item.children()]
    assert names[:3] == ["structure", "identity_to_target", "rms_to_target"], names
    assert type(item.identity_to_target).__name__ == "CFloat"
