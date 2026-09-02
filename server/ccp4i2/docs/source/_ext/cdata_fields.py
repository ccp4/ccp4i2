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


class CDataErrors(SphinxDirective):
    """Render a class's error codes as a table rather than a dict dump.

        .. cdata-errors:: ccp4i2.core.CCP4File.CDataFile

    Codes merge along the MRO, so most classes resolve to scores of them and
    autodoc prints the lot as one unreadable line. Only the codes a class
    *adds* are interesting, so the inherited ones are counted rather than
    listed.
    """
    required_arguments = 1
    option_spec = {}

    SEVERITY = {0: 'OK', 1: 'info', 2: 'warning', 3: 'error', 4: 'error'}

    def run(self):
        cls = _import(self.arguments[0])
        # cdata_class merges codes along the MRO and assigns the result, so a
        # class's own __dict__ holds the whole inherited set. What it *adds*
        # is what no ancestor already had.
        ancestral = {}
        for base in cls.__mro__[1:]:
            ancestral.update(vars(base).get('ERROR_CODES', {}) or {})
        resolved = getattr(cls, 'ERROR_CODES', {}) or {}
        own = {k: v for k, v in resolved.items() if k not in ancestral}
        inherited = {k: v for k, v in resolved.items() if k in ancestral}
        out = []
        if own:
            table = nodes.table(); group = nodes.tgroup(cols=3); table += group
            rows = [(str(k), self.SEVERITY.get((v or {}).get('severity', 3), ''),
                     str((v or {}).get('description', ''))) for k, v in sorted(own.items(), key=lambda kv: str(kv[0]))]
            for i, h in enumerate(('code', 'severity', 'meaning')):
                group += nodes.colspec(colwidth=max([len(h)] + [len(r[i]) for r in rows]))
            head = nodes.thead(); group += head; head += CDataFields._row(('code', 'severity', 'meaning'))
            body = nodes.tbody(); group += body
            for r in rows:
                body += CDataFields._row(r)
            out.append(table)
        if inherited:
            out.append(nodes.paragraph(
                text=f'{len(inherited)} further codes are inherited and shared with '
                     f'every CData; they cover undefined values, type errors and '
                     f'comparison outcomes.'))
        return out or [nodes.paragraph(text='This class declares no error codes.')]


class ErrorCodeIndex(SphinxDirective):
    """Every error code a task can report, indexed by the number itself.

        .. error-code-index::

    For when a code reaches a user without its message. Codes are *per class*:
    only 40 distinct numbers are in use across 79 tasks, so 301 means one thing
    in molrep_mr and another in servalcat_pipe. Sorting by number and listing
    the tasks against it is therefore the way round that answers the question
    someone actually has.

    Read from the source rather than by importing, so a task whose plugin
    cannot be imported without CCP4 still appears.
    """
    required_arguments = 0
    option_spec = {}

    def run(self):
        import ast, pathlib, collections
        root = pathlib.Path(__file__).resolve().parents[3]      # ccp4i2/
        found = collections.defaultdict(list)
        for area in ('wrappers', 'pipelines'):
            for path in sorted((root / area).rglob('*.py')):
                try:
                    tree = ast.parse(path.read_text(errors='ignore'))
                except Exception:
                    continue
                for node in ast.walk(tree):
                    if not isinstance(node, ast.ClassDef):
                        continue
                    for stmt in node.body:
                        if not (isinstance(stmt, ast.Assign)
                                and any(getattr(t, 'id', None) == 'ERROR_CODES' for t in stmt.targets)
                                and isinstance(stmt.value, ast.Dict)):
                            continue
                        for k, v in zip(stmt.value.keys, stmt.value.values):
                            try:
                                code = int(ast.literal_eval(k))
                                entry = ast.literal_eval(v)
                            except Exception:
                                continue
                            desc = (entry or {}).get('description', '') if isinstance(entry, dict) else str(entry)
                            found[code].append((node.name, str(desc)))

        rows = [(str(code), task, desc)
                for code in sorted(found)
                for task, desc in sorted(found[code])]
        if not rows:
            return [nodes.paragraph(text='No task error codes found.')]

        table = nodes.table()
        group = nodes.tgroup(cols=3)
        table += group
        heads = ('code', 'task', 'meaning')
        for i, h in enumerate(heads):
            group += nodes.colspec(colwidth=max([len(h)] + [len(r[i]) for r in rows]))
        head = nodes.thead(); group += head; head += CDataFields._row(heads)
        body = nodes.tbody(); group += body
        for r in rows:
            body += CDataFields._row(r)
        return [table]


def setup(app):
    app.add_directive("cdata-fields", CDataFields)
    app.add_directive("cdata-errors", CDataErrors)
    app.add_directive("error-code-index", ErrorCodeIndex)
    return {"version": "1.0", "parallel_read_safe": True}
