"""Editing an mmCIF block so that every category is written exactly once.

OneDep's upload check answers "Bad model file" -- within seconds, with no
per-file detail -- when a category's items appear in two separate places in
the block. gemmi's ``Block.set_pair`` appends a *new* item at the end of the
block, so adding ``_refine.pdbx_starting_model`` to a model whose ``_refine``
pairs sit near the top splits the category and sinks the upload. Learned the
hard way (eight bisecting uploads, 2026-11-01).

Everything here therefore rewrites a category as a whole: the old items are
erased, the new ones are written, and they are moved to where the old ones
were (or, for a category the block did not have, ahead of ``_atom_site`` so
the header stays readable).
"""
import gemmi

ATOM_SITE = "_atom_site"


def quote(value):
    """A CIF token for ``value``: ``?`` for None/empty, quoted otherwise."""
    if value is None:
        return "?"
    if isinstance(value, float):
        text = repr(value) if value != int(value) or abs(value) > 1e15 else str(value)
    else:
        text = str(value)
    if text == "":
        return "?"
    if text in ("?", "."):
        return text
    return gemmi.cif.quote(text)


def category_of(item):
    """The category name (``_refine``) of a block item, or None for an erased item."""
    if item.pair is not None:
        return item.pair[0].split(".")[0]
    if item.loop is not None and len(item.loop.tags):
        return item.loop.tags[0].split(".")[0]
    return None


def category_runs(block):
    """``{category: number of separate runs of items}`` for the block.

    A category whose items are contiguous has one run; a split category has
    more. ``split_categories`` is what a deposition file must not have.
    """
    runs = {}
    previous = None
    for item in block:
        cat = category_of(item)
        if cat is None:
            continue
        if cat != previous:
            runs[cat] = runs.get(cat, 0) + 1
        previous = cat
    return runs


def split_categories(block):
    """The categories that appear in more than one place in the block."""
    return sorted(cat for cat, n in category_runs(block).items() if n > 1)


def category_rows(block, category):
    """``(tags, rows)`` for a category of a block, whether written as pairs
    or as a loop; tags are bare item names, values raw CIF tokens. Empty
    when the block lacks the category."""
    pair_tags, pair_values = [], []
    for item in block:
        if category_of(item) != category:
            continue
        if item.pair is not None:
            pair_tags.append(item.pair[0].split(".", 1)[1])
            pair_values.append(item.pair[1])
        elif item.loop is not None:
            tags = [t.split(".", 1)[1] for t in item.loop.tags]
            width = item.loop.width()
            return tags, [[item.loop[i, j] for j in range(width)] for i in range(item.loop.length())]
    return pair_tags, ([pair_values] if pair_values else [])


def value_of(block, tag):
    """The first value of ``tag`` as a Python string (pairs or loops), or None."""
    column = block.find_values(tag)
    if not column:
        return None
    first = column[0]
    if first in ("?", "."):
        return None
    return gemmi.cif.as_string(first)


def _indices(block, category):
    return [i for i, item in enumerate(block) if category_of(item) == category]


def _insertion_point(block, category):
    """Where a category the block does not have should go: ahead of the atom
    table when there is one, else at the end."""
    atoms = _indices(block, ATOM_SITE)
    if atoms:
        return atoms[0]
    return sum(1 for _ in block)


def _erase(block, category):
    for i, item in enumerate(block):
        if category_of(item) == category:
            item.erase()


def _move_tail_to(block, position, count):
    """Move the last ``count`` items of the block to ``position``, keeping their order."""
    total = sum(1 for _ in block)
    for j in range(count):
        current = total - count + j
        target = position + j
        # After each move the moved item sits at position+j and the rest of
        # the tail is still at the end. A target at or past the item's own
        # place means it is already where it belongs (the category was last).
        if target < current:
            block.move_item(current, target)


def replace_loop(block, category, tags, rows):
    """Write ``category`` as one loop with ``tags`` (bare item names) and
    ``rows`` (sequences of Python values), in place of whatever the block had.

    An empty ``rows`` removes the category.
    """
    existing = _indices(block, category)
    position = existing[0] if existing else _insertion_point(block, category)
    _erase(block, category)
    if not rows:
        return
    loop = block.init_loop(category + ".", list(tags))
    for row in rows:
        loop.add_row([quote(v) for v in row])
    _move_tail_to(block, position, 1)


def replace_pairs(block, category, items):
    """Write ``category`` as pairs (``items`` is a sequence of (name, value)),
    in place of whatever the block had."""
    existing = _indices(block, category)
    position = existing[0] if existing else _insertion_point(block, category)
    _erase(block, category)
    names = []
    for name, value in items:
        block.set_pair("%s.%s" % (category, name), quote(value))
        names.append(name)
    _move_tail_to(block, position, len(names))


def rewrite_category(block, category, tags, rows):
    """Write a category back from raw CIF tokens, in place: as pairs for one
    row, as a loop otherwise."""
    if len(rows) == 1:
        replace_pairs(block, category, [(t, None) for t in tags])
        for tag, token in zip(tags, rows[0]):
            block.set_pair("%s.%s" % (category, tag), token)
    else:
        replace_loop(block, category, tags, [["?"] * len(tags)])
        loop = block.find_loop("%s.%s" % (category, tags[0])).get_loop()
        loop.set_all_values([[row[j] for row in rows] for j in range(len(tags))])


def set_pairs_in_place(block, category, items):
    """Set pairs of an existing pair-style category without moving it: an
    existing tag is replaced where it is, a new tag is appended at the end of
    the block and then moved in beside the category's last item.

    ``_refine`` is the case: the refinement program wrote it and the new items
    (starting model, method) must join it, not trail the atom table. A
    category written as a one-row loop is first rewritten as pairs (gemmi's
    set_pair on a looped tag would replace the whole loop with that one pair).
    """
    tags, rows = category_rows(block, category)
    if rows and block.find_loop("%s.%s" % (category, tags[0])):
        if len(rows) > 1:
            raise ValueError("%s has %d rows; cannot set single values on it" % (category, len(rows)))
        rewrite_category(block, category, tags, rows)
    for name, value in items:
        tag = "%s.%s" % (category, name)
        if block.find_pair_item(tag) is not None:
            block.set_pair(tag, quote(value))
            continue
        existing = _indices(block, category)
        if not existing:
            replace_pairs(block, category, [(name, value)])
            continue
        block.set_pair(tag, quote(value))
        _move_tail_to(block, existing[-1] + 1, 1)


def validate_with_dictionary(doc, dictionary_path, full=False):
    """Validate a gemmi Document against a DDL2 dictionary.

    Returns the list of messages the dictionary check produced. With
    ``full`` the check is the whole of ``gemmi validate -d`` (mandatory items,
    parent keys, linked groups), which a refinement program's mmCIF fails in
    hundreds of harmless ways (no ``_pdbx_poly_seq_scheme``, refmac's own
    ``_ccp4_*`` tags, the monomer-library blocks it appends). The default is
    the narrower check a deposition file must pass: unknown tags, item types,
    enumerations and unique keys.
    """
    messages = []
    ddl = gemmi.cif.Ddl(lambda text: messages.append(str(text)),
                        print_unknown_tags=True, use_regex=True, use_context=full,
                        use_linked_groups=full, use_mandatory=full, use_unique_keys=True)
    ddl.read_ddl(gemmi.cif.read(str(dictionary_path)))
    # Reading the dictionary itself logs a few complaints about the DDL
    # (regex quirks, double parents); those are about the dictionary, not
    # the document, so only what the document check adds is reported.
    read_noise = len(messages)
    ddl.validate_cif(doc)
    return [m for m in messages[read_noise:] if not m.startswith("Bad DDL2")]


def check_categories(block, dictionary_path, categories):
    """Dictionary messages about ``categories`` (names such as ``_struct_ref``)
    in ``block``: the guard a writer runs on what it wrote, ignoring whatever
    the refinement program's own tables provoke."""
    doc = gemmi.cif.Document()
    copy = doc.add_new_block(block.name)
    wanted = set(categories)
    for item in block:
        cat = category_of(item)
        if cat in wanted:
            copy.add_item(item)
    messages = validate_with_dictionary(doc, dictionary_path)
    return [m for m in messages if "EMDB_PUBLIC" not in m and "WWPDB_LOCAL" not in m
            and "for pdb internal use" not in m]
