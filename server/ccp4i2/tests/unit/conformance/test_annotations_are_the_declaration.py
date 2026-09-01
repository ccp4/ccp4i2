"""A class's annotations are its field declaration --- the only one.

Generated stubs used to declare their fields twice, once as an annotation and
once in the decorator:

    @cdata_class(attributes={"structure": attribute(CUSTOM, "CPdbDataFile")})
    class CPdbEnsembleItem(CData):
        structure: Optional[CPdbDataFile] = None

Only the decorator was load-bearing; the annotation was decoration. That was
inverted first --- the runtime reads the annotations --- and then the duplicate
was removed, because two declarations of one fact can disagree and this pair
had no way to notice if they did.

There was a test here asserting the two agreed, over 319 fields. It has been
deleted rather than adapted: it existed to license removing the decorator, and
once `attributes=` is gone there is nothing to compare a declaration against.
Keeping it would have meant keeping the thing it was written to retire. What
survives are the two statements that do not depend on there being two sources
--- that a declaration builds the object it describes, and that inheritance
supplies what a class does not declare.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


@pytest.fixture(scope="module", autouse=True)
def _django():
    import os
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    django.setup()


def test_a_class_inherits_what_it_does_not_declare():
    """`CDataFile` declares these; its subclasses no longer restate them.

    This used to read "CDataFile is hand-written: decorator, no annotations",
    which stopped being true when it was given annotations, and the test kept
    passing while describing the opposite of the arrangement it checked. The
    assertion was always the useful part.
    """
    from ccp4i2.core.base_object.cdata_file import CDataFile
    from ccp4i2.core.CCP4ModelData import CPdbDataFile

    declared = {c._name for c in CDataFile().children()}
    assert {"baseName", "relPath", "project", "annotation"} <= declared, \
        f"the declaration did not build the fields: {sorted(declared)}"

    # The subclass declares none of these and must still have them.
    inherited = {c._name for c in CPdbDataFile().children()}
    assert {"baseName", "relPath", "project", "annotation"} <= inherited, \
        f"inheritance did not supply the base fields: {sorted(inherited)}"


def test_the_declaration_actually_builds_the_object():
    """The end-to-end statement: annotations produce the right children."""
    from ccp4i2.core.CCP4ModelData import CPdbEnsembleItem

    item = CPdbEnsembleItem()
    names = [c._name for c in item.children()]
    assert names[:3] == ["structure", "identity_to_target", "rms_to_target"], names
    assert type(item.identity_to_target).__name__ == "CFloat"


def test_project_keeps_its_type_through_inheritance():
    """`CProjectId`, not `CUUID` --- the qualifier that restatement concealed.

    47 stub classes restated `project` as a CProjectId while the base declared
    CUUID, so removing the restatements would have widened the type in every
    one of them and dropped `allowUnfound`. The base was corrected first; this
    pins it, because nothing else would notice.
    """
    from ccp4i2.core.CCP4ModelData import CPdbDataFile

    project = CPdbDataFile().project
    assert type(project).__name__ == "CProjectId", type(project).__name__
    assert project.qualifiers().get("allowUnfound") is True
