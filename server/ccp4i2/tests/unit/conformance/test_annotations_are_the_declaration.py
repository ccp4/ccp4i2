"""A stub's annotations are its field declaration.

Every generated stub declares its fields twice:

    @cdata_class(attributes={"structure": attribute(CUSTOM, "CPdbDataFileStub"), ...})
    class CPdbEnsembleItemStub(CData):
        structure: Optional[CPdbDataFileStub] = None

Only the decorator was load-bearing; the annotation was decoration, and the
`= None` is what the runtime walk overwrites. But the annotation is the
declaration a dataclass uses, and it carries strictly more --- the class
itself, rather than its name as a string to be looked up later.

Measured across the stubs the two agree exactly. At the time this was written
that was 621 fields over 122 classes; removing the inherited restatements
later reduced it to 319 fields over 82, because a class that adds nothing now
declares nothing. The agreement is the invariant, not the count.

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


def _built_class(definition):
    """The class a definition builds, rather than how it spells it.

    The two encodings are equivalent and both appear: a generated stub writes
    `attribute(AttributeType.STRING)` while `CDataFile`, hand-written, writes
    `attribute(CUSTOM, custom_class="CString")`. Comparing the spelling would
    report a difference where there is none.
    """
    from ccp4i2.core.base_object.class_metadata import AttributeType
    if definition.custom_class:
        return definition.custom_class
    return {AttributeType.STRING: "CString", AttributeType.INT: "CInt",
            AttributeType.FLOAT: "CFloat", AttributeType.BOOLEAN: "CBoolean",
            }.get(definition.attr_type, definition.attr_type)


def test_annotations_reproduce_the_decorator():
    """The two declarations must agree, or reading one changes behaviour."""
    from ccp4i2.core.base_object.class_metadata import attributes_from_annotations

    checked, disagreeing = 0, []
    for klass in _stub_classes():
        declared = vars(klass)["_metadata"].attributes
        from_annotations = attributes_from_annotations(klass)
        if not from_annotations:
            continue          # no annotations at all: the fallback case
        checked += len(from_annotations)
        a = {n: _built_class(d) for n, d in declared.items()}
        b = {n: _built_class(d) for n, d in from_annotations.items()}
        if a != b:
            disagreeing.append(
                f"{klass.__name__}: {[n for n in a if b.get(n) != a[n]][:3]}")

    # Counted in fields, not classes. The class count is not a measure of
    # coverage and fell from 122 to 82 when the restatements were removed ---
    # a class that declares nothing and inherits everything is the intended
    # end state, not a weaker test. The fields are still compared wherever
    # they are declared.
    assert checked > 250, f"only {checked} fields checked --- proves little"
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
