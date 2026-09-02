"""A Sphinx directive that documents a CData class's fields from the class.

    .. cdata-fields:: ccp4i2.core.CCP4XtalData.CCell

renders a table of the fields the class declares and inherits. Nothing is
generated into the tree and nothing needs regenerating: the table is built from
the live class when the docs are built, so it cannot drift from the code.

The descriptions are not written twice either. A field's `toolTip` and
`guiLabel` are declared for the GUI, and this reads the same declarations:

    a = content(CCellLength, toolTip='Cell length a in A', guiLabel='a')

Qualifiers are read *as declared*, walking the MRO, rather than as resolved on
an instance. Resolved, a class's own qualifiers are inherited by every child,
so `CObsDataFile.project` reports the file's guiLabel ("Reflections") and its
toolTip ("Observed structure factors or intensities") --- confidently wrong for
a field that holds a project id.
"""
from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.util.docutils import SphinxDirective

COLUMNS = ("field", "type", "declared by", "label", "description")


def _declared(cls):
    """Per-field (type, qualifiers, declaring class), walking the MRO."""
    from ccp4i2.core.base_object.class_metadata import contents_from_declarations
    out = {}
    for klass in reversed(cls.__mro__):
        for name, decl in contents_from_declarations(klass).items():
            out[name] = (decl.cls if isinstance(decl.cls, str) else decl.cls.__name__,
                         dict(decl.qualifiers), klass.__name__)
    return out


def _import(path):
    module, _, name = path.rpartition(".")
    import importlib
    return getattr(importlib.import_module(module), name)


class CDataFields(SphinxDirective):
    required_arguments = 1
    option_spec = {}

    def run(self):
        cls = _import(self.arguments[0])
        declared = _declared(cls)
        try:
            order = list(cls().CONTENTS)
        except Exception:
            order = list(declared)
        if not order:
            return [nodes.paragraph(text="This class declares no fields.")]

        rows = []
        for name in order:
            kind, quals, owner = declared.get(name, ("", {}, ""))
            rows.append((name, kind, owner, str(quals.get("guiLabel", "")),
                         str(quals.get("toolTip", ""))))

        table = nodes.table()
        group = nodes.tgroup(cols=len(COLUMNS))
        table += group
        # width by content, or every column is squeezed to the same few
        # characters and the table becomes unreadable in the text builder
        for i, heading in enumerate(COLUMNS):
            group += nodes.colspec(
                colwidth=max([len(heading)] + [len(r[i]) for r in rows]))
        head = nodes.thead()
        group += head
        head += self._row(COLUMNS)
        body = nodes.tbody()
        group += body
        for row in rows:
            body += self._row(row)
        return [table]

    @staticmethod
    def _row(cells):
        row = nodes.row()
        for cell in cells:
            entry = nodes.entry()
            entry += nodes.paragraph(text=str(cell))
            row += entry
        return row


def setup(app):
    app.add_directive("cdata-fields", CDataFields)
    return {"version": "1.0", "parallel_read_safe": True}
