"""What order children come out in, when nobody has said.

`contents_order` orders an unordered bag. It never contracted to span every
field: the ones it names come first, in the order given, and the rest follow.
So the interesting question is what "the rest" means, and the answer used to
be alphabetical --- not because anyone chose it, but because `CONTENTS_ORDER`
is a *property* whose last resort was `sorted(self.CONTENTS)`, and
`dataOrder()` stage 1 asks `hasattr(self, 'CONTENTS_ORDER')`, which is always
true for a property. Stage 1 therefore always won, and stage 3's MRO walk ---
which already computed declaration order correctly --- was unreachable.

The visible cost was borne by `CDataFile`: declared `project, baseName,
relPath, dbFileId, annotation, subType, contentFlag`, rendered `annotation,
baseName, contentFlag, dbFileId, project, relPath, subType`, in every task and
for every file parameter.

These tests pin the three properties that matter, in the order they matter:
declaration order is the fallback, an explicit `contents_order` still wins, and
fields it does not name still appear.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


@pytest.fixture(scope="module", autouse=True)
def _django():
    import os
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")
    django.setup()


def test_a_data_file_renders_in_the_order_it_is_declared():
    """The real-world case, and the one that moved 164 tasks.

    Named explicitly rather than derived, so that reordering the declaration
    has to be a deliberate act that fails this test, not a silent drift.
    """
    from ccp4i2.core.base_object.cdata_file import CDataFile

    assert list(CDataFile().dataOrder()) == [
        "project", "baseName", "relPath", "dbFileId",
        "annotation", "subType", "contentFlag",
    ]


def test_declaration_order_is_not_alphabetical_order():
    """The property under test, on a class the fallback actually governs.

    This deliberately does not use `CCell`, though a unit cell rendered
    `a, alpha, b, beta, c, gamma` is the vivid illustration. `CCell` sets
    `contents_order`, so it escaped the alphabetical fallback and a test
    written on it passes with the defect present --- which is what the first
    draft of this test did.

    `CSequenceMeta` sets no `contents_order`, so it is governed by the
    fallback and nothing else.
    """
    from ccp4i2.core.CCP4ModelData import CSequenceMeta

    meta = CSequenceMeta()
    assert meta._metadata.contents_order in (None, [], ()), \
        "pick a class the fallback actually governs, or this proves nothing"

    order = list(meta.dataOrder())
    assert order == ["uniprotId", "organism", "expressionSystem"]
    assert order != sorted(order), "alphabetical is not an ordering anyone chose"


def test_an_explicit_contents_order_still_wins():
    """The override has to keep working, or the fallback is a regression.

    This test has now been re-pointed twice, which is itself the lesson: it
    used `CCell`, then `CImportUnmerged`, and both later stopped declaring a
    `contents_order` because their declaration order says the same thing. A
    test of an override must use a case where the override is *structurally*
    required, not merely present.

    `CAsuDataFile` is that case. It declares `selection` and inherits
    `project`, `baseName` and the rest from `CDataFile`. Inherited fields come
    first in MRO order, so no arrangement of its own declarations can put
    `selection` in front --- only `contents_order` can.
    """
    from ccp4i2.core.CCP4ModelData import CAsuDataFile

    obj = CAsuDataFile()
    assert obj._metadata.contents_order, "this class no longer proves anything"

    order = list(obj.dataOrder())
    assert order[0] == "selection", order[:3]
    assert "project" in order, "an inherited field was dropped, not reordered"
    assert order.index("selection") < order.index("project"), \
        "the override did not hoist the class's own field"


def test_fields_absent_from_contents_order_are_still_rendered():
    """`contents_order` is a partial hint, not a whitelist.

    `CAsuDataFile` names only `selection` --- hoisting its own field ahead of
    the inherited ones, which is the one thing declaration order cannot express,
    since MRO order puts a parent's fields first. The inherited fields must
    still follow rather than vanish.
    """
    from ccp4i2.core.CCP4ModelData import CAsuDataFile

    order = list(CAsuDataFile().dataOrder())
    assert order[0] == "selection"
    for inherited in ("project", "baseName", "relPath", "annotation"):
        assert inherited in order, f"{inherited} was dropped, not merely reordered"
