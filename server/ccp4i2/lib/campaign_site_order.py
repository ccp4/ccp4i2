"""The order a campaign's sites are shown in: by name, alphabetically.

Names are how people group sites -- a site and its subsites share a stem
(``Pocket``, ``Pocket_A``, ``Pocket_A_Trp``) -- so sorting by name puts like
with like. The client sorts every list of sites with ``Intl.Collator("en",
{numeric: true, sensitivity: "base"})`` (``client/renderer/lib/site-matrix.ts``,
``compareSites``); this key reproduces that order on the server, so that the
"first site" the server picks (the summary scene's camera) is the first site
the user sees. Both sides are pinned to the same examples in their tests.

The collator's order, for the names sites actually get:

  * case is ignored;
  * a run of digits compares as a number ("Site 2" before "Site 10");
  * at any one position, space < punctuation < digits < letters, and
    punctuation in ICU's order (``PUNCTUATION_ORDER``, read from the
    collator itself) -- so "Site_A" comes before "Site-A";
  * a name that is a prefix of another comes first.

Ties (names equal but for case) fall back to creation order (the id).
"""
import re

#: ASCII punctuation and symbols in the order the client's collator puts
#: them, space first. Read from Node's Intl.Collator, not written by hand.
PUNCTUATION_ORDER = " _-,;:!?.'\"()[]{}@*/\\&#%`^+<=>|~$"

_RUNS = re.compile(r"\d+|\D")


def _name_key(name: str) -> tuple:
    key = []
    for run in _RUNS.findall(name or ""):
        if run.isdigit():
            key.append((2, int(run), ""))
        elif run in PUNCTUATION_ORDER:
            key.append((1, PUNCTUATION_ORDER.index(run), ""))
        elif run.isalpha():
            key.append((3, 0, run.casefold()))
        else:
            # Anything else (non-ASCII symbols): after the letters, by
            # codepoint. Rare in a site name; deterministic is enough.
            key.append((4, ord(run), ""))
    return tuple(key)


def site_sort_key(site) -> tuple:
    """Sort key for a ``CampaignSite`` (anything with ``name`` and ``id``)."""
    return (_name_key(site.name), site.id)


def ordered_sites(sites) -> list:
    """The sites, in display order."""
    return sorted(sites, key=site_sort_key)
