"""
The order a campaign's sites are shown in (lib/campaign_site_order).

The client sorts with Intl.Collator("en", {numeric: true, sensitivity:
"base"}); the server must agree, so that its "first site" is the one the user
sees first. ``EXPECTED`` is the collator's own order for ``NAMES``, produced
by Node -- the same list pins client/renderer/__tests__/site-order.test.ts.
Change one, change both.

Pure Python -- no Django, no CCP4.
"""
from types import SimpleNamespace

from ccp4i2.lib.campaign_site_order import ordered_sites, site_sort_key

NAMES = [
    "Site 10", "Pocket_A_Tyr127", "site 1", "Pocket-B", "Interface_AB_top",
    "x0010", "Pocket", "Hinge_B", "Pocket_A-Met", "Site 2", "Pocket B",
    "pocket_a", "Back_A", "Interface_AB_low", "Pocket_A_Trp", "x0009",
    "Pocket_A_Tyr12", "Pocket_A",
]

EXPECTED = [
    "Back_A", "Hinge_B", "Interface_AB_low", "Interface_AB_top", "Pocket",
    "Pocket B", "pocket_a", "Pocket_A", "Pocket_A_Trp", "Pocket_A_Tyr12",
    "Pocket_A_Tyr127", "Pocket_A-Met", "Pocket-B", "site 1", "Site 2",
    "Site 10", "x0009", "x0010",
]


def _sites(names):
    return [SimpleNamespace(name=name, id=i + 1) for i, name in enumerate(names)]


def test_matches_the_client_collator():
    assert [s.name for s in ordered_sites(_sites(NAMES))] == EXPECTED


def test_a_site_and_its_subsites_sit_together():
    names = ["Pocket_B", "Other", "Pocket_A_Trp", "Pocket_A", "Pocket_B_Met", "Pocket_A_Met"]
    assert [s.name for s in ordered_sites(_sites(names))] == [
        "Other", "Pocket_A", "Pocket_A_Met", "Pocket_A_Trp", "Pocket_B", "Pocket_B_Met",
    ]


def test_numbers_compare_as_numbers():
    assert [s.name for s in ordered_sites(_sites(["Site 10", "Site 9", "Site 100"]))] == [
        "Site 9", "Site 10", "Site 100",
    ]


def test_names_equal_but_for_case_keep_creation_order():
    late = SimpleNamespace(name="pocket", id=7)
    early = SimpleNamespace(name="POCKET", id=3)
    assert ordered_sites([late, early]) == [early, late]
    assert site_sort_key(early) < site_sort_key(late)


def test_empty_and_missing_names_sort_first_without_error():
    sites = [SimpleNamespace(name="A", id=1), SimpleNamespace(name="", id=2),
             SimpleNamespace(name=None, id=3)]
    assert [s.id for s in ordered_sites(sites)] == [2, 3, 1]
