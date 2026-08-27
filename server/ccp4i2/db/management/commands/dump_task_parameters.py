# Copyright (C) 2026 Newcastle University
"""Emit the parameters every registered task declares, for the GUI coverage test.

A task interface is only as good as the parameters it lets a user reach, and
nothing tested that: ``arp_warp_classic`` shipped as two empty accordions
against a 47-parameter def.xml and no tier noticed, because the i2run suite
drives the backend through the CLI and the api/e2e suite through REST --- both
bypass React entirely.

This command writes the *declared* side of that comparison to a JSON fixture.
The frontend test (``renderer/__tests__/task-parameter-coverage.test.ts``)
supplies the *exposed* side by reading the ``.tsx`` interfaces, so the check
itself needs no Python, no CCP4 and no database, and runs in the CCP4-free CI
job.

Enumeration goes through :func:`parse_def_xml_file`, the same loader the server
uses, so composed definitions resolve the way they do in a running job ---
``servalcat_pipe`` inherits ``servalcat.def.xml``, ``prosmart_refmac`` inherits
``refmac``, and a naive read of the outer file sees a handful of parameters
where a job sees a hundred.

The parsed tree cannot be taken at face value, though. Its second pass
re-parses every nested container onto the root, so phaser_rnp_pipeline's four
root containers come back as twenty-seven, and 1197 parameters across 15 tasks
appear twice -- once in their real place under ``keywords``, once as a phantom
scope at the top. :func:`_file_root_containers` reads the composed XML to settle
which containers are genuinely at root, so nothing is counted twice here
whatever the parser does.

One thing this deliberately does not do: treat the declared set as a target. It
is an upper bound. crank2 declares 238 parameters and its Qt GUI exposed 111, so
the question a reviewer wants answered is "does our interface reach what the Qt
one reached", which is what ``--qt-baseline`` records.

That Qt figure has to account for auto-generation or it libels the port. A Qt
``CTaskWidget`` did not have to name a parameter to show it: ``autoGenerate``
(directly, or via ``openFolder(autoSelection=...)``) walks a container's
``dataOrder()`` and lays out everything matching a selection. ``prosmart`` draws
its entire GUI that way and names three parameters explicitly, so counting only
``createLine`` calls would score a 65-parameter task as a 3-parameter one and
make our interface look finished. Fifteen Qt GUIs use it. This module reproduces
the selection logic --- including that the wildcards are ``re.match`` patterns
rather than globs --- against the declared set.

Usage::

    ccp4-python manage.py dump_task_parameters
    ccp4-python manage.py dump_task_parameters --qt-baseline ccp4/main
"""
import ast
import io
import json
import re
import subprocess
import tokenize
import xml.etree.ElementTree as ET
from pathlib import Path

from django.core.management.base import BaseCommand

from ccp4i2.core.tasks import TASKS, locate_def_xml
from ccp4i2.core.task_manager.def_xml_handler import parse_def_xml_file
from ccp4i2.lib.utils.parameters.load_xml import load_nested_xml

#: Root containers kept out of the coverage comparison, but still reported so a
#: check for names that resolve to nothing can see them.
#:
#: ``outputData`` is named by the plugin rather than chosen in the GUI.
#: ``guiParameters`` is a mixed bag across the twenty tasks that declare one:
#: mostly the GUI's own state (``OPEN_*`` subframe toggles, ``EXPERT_LEVEL``)
#: and cached digest results (``MMCIF_*``, ``HKLINISMTZ``) -- but also a handful
#: of genuine user controls, ``ProvideTLS.EDIT_MODE``, ``molrep_*.PERFORM``,
#: ``mrbump_basic.USE_FIXED``. Counting the whole container as parameters would
#: manufacture work, because React expresses an ``OPEN_*`` toggle as a
#: ``visibility`` prop and rightly has no field for it; counting none of it
#: calls those genuine controls dead. So: reported, not counted.
#:
#: Everything else at root does count, which is wider than the obvious two: the
#: phaser family keeps its keywords in a dozen sibling containers
#: (``macmr``, ``ellg``, ``tncs``, ...) that ``phaser_EP_LLG`` and friends draw
#: with ``autoGenerate``. Walking only inputData/controlParameters would hide
#: them from both sides of the comparison at once, which looks like agreement.
NON_PARAMETER_SECTIONS = ("outputData", "guiParameters")

#: Sections listed first when present, so the fixture reads in a familiar order.
PREFERRED_SECTION_ORDER = ("inputData", "controlParameters")

#: Tasks with no Qt GUI of their own that should be judged against another
#: task's, because they are the same machinery under a different name.
#:
#: ``shelx`` is the crank2 pipeline driven through SHELX C/D/E. It arrived
#: after the Qt era so nothing on ``main`` is called that, and without this it
#: would be exempt from the comparison entirely -- which is how it came to
#: reach 56 of the 108 parameters crank2's GUI offers while looking unmeasured
#: rather than incomplete. Qt's crank2 GUI was curated by the pipeline's author,
#: so it is the authority on what this machinery should put in front of a user,
#: whichever task name the user reached it by.
#:
#: Applied only when the two tasks declare exactly the same parameters, so a
#: definition drifting apart later withdraws the borrowed reference rather than
#: silently judging one task by another's contents.
QT_REFERENCE_ALIASES = {"shelx": "crank2"}

#: Guard against a cyclic definition turning the walk into a hang.
MAX_DEPTH = 8

#: .../<root>/server/ccp4i2/db/management/commands/<this file>
REPO_ROOT = Path(__file__).resolve().parents[5]

DEFAULT_OUT = (
    REPO_ROOT
    / "client"
    / "renderer"
    / "__tests__"
    / "fixtures"
    / "task-parameters.json"
)


def _parameter_paths(container, prefix="", depth=0):
    """Dotted paths of every leaf under *container*, in declaration order.

    ``dataOrder()`` is the accessor that answers "what did the def.xml
    declare", and it preserves file order --- unlike ``children()``, which is
    backed by a weakref set hashed on identity and so comes back in a different
    order run to run.
    """
    paths = []
    if depth > MAX_DEPTH:
        return paths
    for name in container.dataOrder() or []:
        child = getattr(container, name, None)
        if child is None:
            continue
        if type(child).__name__ == "CContainer":
            paths += _parameter_paths(child, f"{prefix}{name}.", depth + 1)
        else:
            paths.append(f"{prefix}{name}")
    return paths


def _is_phil_driven(task):
    """Does this task build its parameters from a PHIL tree at runtime?

    A ``PhilPluginScript`` reads its parameters out of the program it drives
    (xia2, phaser, phasertng) when the job is set up, so its def.xml declares
    almost nothing and its interface names fields no static reading of the
    def.xml can know about. Those names are correct, not dead, and a checker
    that does not know this reports fifteen phantom fields in
    ``xia2_ssx_reduce`` alone.

    Asked of the plugin class rather than by grepping for the base class name,
    so a task that inherits it at one remove -- ``xia2_xds`` from
    ``xia2_dials`` -- answers correctly.
    """
    try:
        from ccp4i2.core.PhilPluginScript import PhilPluginScript
        from ccp4i2.core.tasks import get_plugin_class

        plugin = get_plugin_class(task)
        return plugin is not None and issubclass(plugin, PhilPluginScript)
    except Exception:                                # noqa: BLE001 - absent is an answer
        return False


def _file_root_containers(xml_path):
    """Container ids that genuinely sit at the top of the composed def.xml.

    Needed because the parser puts more there than the file does. Its second
    pass collects ``.//container[@id]`` -- every container at any depth -- and
    re-parses onto the root any whose id it has not already seen, so a scope
    nested inside another appears twice: once in its real place and once at the
    top. phaser_rnp_pipeline declares four root containers and comes back with
    twenty-seven. Across the registry that is 222 phantom containers holding
    1197 duplicate parameters in 15 tasks.

    Counting those would inflate every total and make one interface field look
    like two satisfied parameters. Reading the file settles it, and keeps this
    correct whichever way the parser behaves later.
    """
    try:
        root = load_nested_xml(ET.parse(xml_path).getroot())
    except Exception:                                # noqa: BLE001 - fall back below
        return None
    # Most def.xml files identify the body, but cmapcoeff, comit, cpatterson,
    # fft and validate_protein write a bare <ccp4i2_body>. Looking only for the
    # identified form leaves those five trusting whatever the parser produced.
    body = root.find(".//ccp4i2_body[@id]")
    if body is None:
        body = root.find(".//ccp4i2_body")
    if body is None:
        return None

    names, stack = set(), [body]
    while stack:
        element = stack.pop()
        for child in element:
            if child.tag == "container" and child.get("id"):
                names.add(child.get("id"))
            elif child.tag == "ccp4i2_body":
                # Composition splices an included def.xml in as a nested body;
                # its containers are root-level for our purposes.
                stack.append(child)
    return names


def _sections(root, xml_path):
    """The root containers whose contents a user is meant to set, in order."""
    genuine = _file_root_containers(xml_path)
    present = [
        name for name in (root.dataOrder() or [])
        if name not in NON_PARAMETER_SECTIONS
        and type(getattr(root, name, None)).__name__ == "CContainer"
        and (genuine is None or name in genuine)
    ]
    preferred = [n for n in PREFERRED_SECTION_ORDER if n in present]
    return preferred + [n for n in present if n not in preferred]


def _without_comments(source):
    """*source* with Python comments blanked out.

    Commenting a ``createLine`` out is how these GUIs retire a parameter, and
    they are full of it -- crank2 alone keeps ten commented-out widgets,
    ``KEYWORDS_PHAS`` and ``PRESENT_STYLE`` among them. Reading those as
    exposed would hold the port to parameters the GUI's own author decided not
    to show, which is the opposite of taking that GUI as authoritative.

    Tokenised rather than regexed, so a ``#`` inside a string stays put.
    """
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        # A file we cannot tokenise (Python 2 leftovers) keeps its comments;
        # over-reporting a few parameters beats dropping the whole GUI.
        return source

    lines = source.splitlines(keepends=True)
    for token in tokens:
        if token.type != tokenize.COMMENT:
            continue
        row, col = token.start[0] - 1, token.start[1]
        if 0 <= row < len(lines):
            line = lines[row]
            newline = "\n" if line.endswith("\n") else ""
            lines[row] = line[:col] + newline
    return "".join(lines)


def _withheld_by_taskname(source):
    """Parameters a shared GUI shows to every task *except* a named one.

    ``CTaskCrank2`` serves both crank2 and shelx and branches on its own
    ``TASKNAME``: ``if self.TASKNAME != 'shelx':`` guards the FA-estimation and
    detection program choosers, the partial-model inputs, and the SHELXC/D/E
    switch -- all fixed or meaningless when the pipeline *is* SHELX. Borrowing
    crank2's parameter list for shelx without honouring those branches would
    manufacture five parameters its author decided not to show.

    Returns ``{task name: parameters withheld from it}``.
    """
    withheld = {}
    lines = source.splitlines()
    for i, line in enumerate(lines):
        match = re.search(r"if\s+self\.TASKNAME\s*!=\s*['\"](\w+)['\"]\s*:", line)
        if not match:
            continue
        indent = len(line) - len(line.lstrip())
        block = []
        for following in lines[i + 1:]:
            if following.strip() and len(following) - len(following.lstrip()) <= indent:
                break
            block.append(following)
        found = _createline_parameters("\n".join(block))
        if found:
            withheld.setdefault(match.group(1), set()).update(found)
    return withheld


def _createline_token_lists(source):
    """The string tokens of each ``createLine([...])`` definition list.

    Parsed rather than regexed. A naive scan for quoted runs loses alignment on
    the first escaped quote, and these GUIs are full of them --- "Use GYRE
    option when running Phaser\\'s Rotation Function" hid both BORGES_GYRE and
    BORGES_GYRE_T from the Qt reference, along with the GIMBLE and LLG pairs, so
    arcimboldo's interface looked complete where it was not.

    Falls back to the regex only for a file Python cannot parse, where an
    approximate answer beats no answer.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        out = []
        for call in re.finditer(r"createLine\s*\(\s*\[(.*?)\]", source, re.S):
            out.append(re.findall(r"['\"]([^'\"]*)['\"]", call.group(1)))
        return out

    lists = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = node.func
        called = getattr(name, "attr", None) or getattr(name, "id", None)
        if called != "createLine":
            continue
        # The definition is the first positional argument, or `definition=`.
        argument = next(iter(node.args), None)
        for keyword in node.keywords:
            if keyword.arg == "definition":
                argument = keyword.value
        if not isinstance(argument, (ast.List, ast.Tuple)):
            continue
        lists.append(
            [e.value for e in argument.elts
             if isinstance(e, ast.Constant) and isinstance(e.value, str)]
        )
    return lists


def _createline_parameters(source):
    """Parameters named explicitly in ``createLine`` calls.

    An entry looks like ``['widget', '-guiMode', 'multiLine', 'PARAM']``:
    options come in ``-key value`` pairs, and the parameter is the first plain
    token after them.
    """
    params = set()
    for tokens in _createline_token_lists(source):
        i = 0
        while i < len(tokens):
            if tokens[i] != "widget":
                i += 1
                continue
            j = i + 1
            while j < len(tokens) and tokens[j].startswith("-"):
                j += 2
            if j < len(tokens) and tokens[j] not in ("widget", "label", "stretch"):
                params.add(tokens[j])
            i = j + 1
    return params


def _balanced(source, start, opener, closer):
    """Substring from *start* (which must be *opener*) to its matching *closer*."""
    depth = 0
    for i in range(start, len(source)):
        if source[i] == opener:
            depth += 1
        elif source[i] == closer:
            depth -= 1
            if depth == 0:
                return source[start:i + 1]
    return None


def _literal_selection(call):
    """The ``selection=``/``autoSelection=`` dict inside one call, if literal.

    Anything not a plain literal is skipped: ``xia2_dials`` builds its selection
    at runtime from a helper, and a selection that cannot be resolved without
    running the GUI would have to be guessed at. Guessing would put parameters
    into the baseline that Qt may never have shown, which is the one error this
    file must not make --- it would manufacture work rather than record it.
    """
    match = re.search(r"\b(?:auto)?[Ss]election\s*=\s*\{", call)
    if not match:
        return None
    text = _balanced(call, match.end() - 1, "{", "}")
    if text is None:
        return None
    try:
        value = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return None
    return value if isinstance(value, dict) else None


def _autogenerate_calls(source):
    """``(container name, selection)`` for each auto-generating call in *source*.

    Two spellings reach the same code. ``self.autoGenerate(container=...)`` names
    its container --- often ``controlParameters``, but ``phaser_EP_LLG`` draws
    ``self.container.keywords`` this way. ``openFolder(autoSelection=...)`` takes
    no container and always applies to ``controlParameters``
    (``CCP4TaskWidget.openFolder`` hard-codes that).
    """
    calls, unresolved = [], 0
    for match in re.finditer(r"\bautoGenerate\s*\(", source):
        call = _balanced(source, match.end() - 1, "(", ")")
        if call is None:
            continue
        container = re.search(r"container\s*=\s*self\.container\.(\w+)", call)
        selection = _literal_selection(call)
        if selection is None:
            unresolved += 1
            continue
        calls.append((container.group(1) if container else "controlParameters", selection))
    for match in re.finditer(r"\bopenFolder\s*\(", source):
        call = _balanced(source, match.end() - 1, "(", ")")
        if call is None or "autoSelection" not in call:
            continue
        selection = _literal_selection(call)
        if selection is None:
            unresolved += 1
        else:
            calls.append(("controlParameters", selection))
    return calls, unresolved


def _autogenerated_parameters(source, declared):
    """Parameters an ``autoGenerate`` call would lay out, given the declared set.

    Mirrors ``CCP4TaskWidget.autoGenerate``: with ``includeParameters``, a name
    is shown if it is listed or if some listed pattern ``re.match``es it (the
    ``*`` in ``'CLUSTER_*'`` is a regex quantifier on the underscore, not a
    glob); without it, everything not in ``excludeParameters`` is shown.

    ``autoGenerate`` recurses into sub-containers, so a dotted path matches on
    any of its segments.

    Returns the parameters found and how many calls could not be read
    statically, so a task whose Qt figure is knowingly incomplete can say so
    rather than passing quietly.
    """
    shown = set()
    calls, unresolved = _autogenerate_calls(source)
    for container, selection in calls:
        paths = declared.get(container, [])
        include = selection.get("includeParameters", [])
        exclude = set(selection.get("excludeParameters", []))
        patterns = [p for p in include if "*" in p]

        for path in paths:
            segments = path.split(".")
            if include:
                hit = any(s in include for s in segments) or any(
                    re.match(p, s) for p in patterns for s in segments
                )
            else:
                hit = not any(s in exclude for s in segments)
            if hit:
                shown.add(path)
    return shown, unresolved


def _qt_exposed(ref, declared):
    """Parameters each Qt ``CTaskWidget`` put in front of a user, by task name.

    The Qt GUIs are the only record of which of a task's declared parameters
    were judged worth showing, so they are the practical target for a port. We
    read them out of *ref* (a git ref, normally ``ccp4/main``) rather than
    vendoring a copy, so refreshing the hint is one flag rather than an edit.

    *declared* maps task name to its section lists, and is needed to resolve
    what an ``autoGenerate`` call actually put on screen.
    """
    listing = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-tree", "-r", "--name-only", ref],
        capture_output=True, text=True, check=True,
    ).stdout.splitlines()
    candidates = [
        f for f in listing
        if f.endswith(".py") and f.split("/")[0] in ("wrappers", "wrappers2", "pipelines")
    ]

    exposed, auto_generated, unreadable = {}, set(), set()
    withheld = {}
    for path in candidates:
        source = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "show", f"{ref}:{path}"],
            capture_output=True, text=True,
        ).stdout
        if "CTaskWidget" not in source:
            continue
        source = _without_comments(source)
        name = re.search(r"^\s*TASKNAME\s*=\s*['\"]([^'\"]+)['\"]", source, re.M)
        if not name:
            continue
        task = name.group(1)

        params = _createline_parameters(source)
        auto, unresolved = _autogenerated_parameters(source, declared.get(task, {}))
        if auto:
            auto_generated.add(task)
        if unresolved:
            unreadable.add(task)
        params |= auto

        if params:
            exposed.setdefault(task, set()).update(params)

        for other, hidden in _withheld_by_taskname(source).items():
            withheld.setdefault(other, set()).update(hidden)

    return (
        {task: sorted(params) for task, params in exposed.items()},
        sorted(auto_generated),
        sorted(unreadable),
        {task: sorted(names) for task, names in withheld.items()},
    )


class Command(BaseCommand):
    help = "Write the declared-parameter fixture used by the GUI coverage test"
    requires_system_checks = []

    def add_arguments(self, parser):
        parser.add_argument(
            "--out", type=Path, default=DEFAULT_OUT,
            help=f"Fixture to write (default: {DEFAULT_OUT})",
        )
        parser.add_argument(
            "--qt-baseline", metavar="GIT_REF",
            help="Refresh the Qt-exposed hint from this ref, e.g. ccp4/main. "
                 "Omitted, whatever the fixture already records is kept.",
        )

    def handle(self, **options):
        out = options["out"]
        previous = {}
        if out.exists():
            previous = json.loads(out.read_text()).get("tasks", {})

        tasks, skipped = {}, {}
        for task in sorted(TASKS):
            def_xml = locate_def_xml(task)
            if not def_xml:
                skipped[task] = "no def.xml"
                continue
            try:
                root = parse_def_xml_file(def_xml)
            except Exception as exc:                    # noqa: BLE001 - report, don't fail
                skipped[task] = f"{type(exc).__name__}: {exc}"
                continue

            declared = {}
            for section in _sections(root, def_xml):
                paths = _parameter_paths(getattr(root, section))
                if paths:
                    declared[section] = paths
            if not declared:
                skipped[task] = "declares no input or control parameters"
                continue
            entry = {"declared": declared}
            # Reported but not counted -- see NON_PARAMETER_SECTIONS. An
            # interface may legitimately render either (mergeMtz and mtzutils
            # both show HKLOUT; ProvideTLS shows EDIT_MODE), so a check for
            # itemNames that name nothing has to know they exist.
            for section in NON_PARAMETER_SECTIONS:
                container = getattr(root, section, None)
                if container is not None:
                    paths = _parameter_paths(container)
                    if paths:
                        entry[section] = paths
            if _is_phil_driven(task):
                entry["philDriven"] = True
            tasks[task] = entry

        if options["qt_baseline"]:
            self.stdout.write(f"Reading Qt GUIs from {options['qt_baseline']}...")
            hint, auto, unreadable, withheld = _qt_exposed(
                options["qt_baseline"],
                {t: v["declared"] for t, v in tasks.items()},
            )
            self.stdout.write(
                f"  {len(hint)} Qt task GUIs, {len(auto)} auto-generating "
                f"({', '.join(auto)}), {len(unreadable)} with a selection that "
                f"cannot be read statically ({', '.join(unreadable)})"
            )
            for borrower, source in QT_REFERENCE_ALIASES.items():
                if borrower in hint or source not in hint:
                    continue
                if tasks.get(borrower, {}).get("declared") != tasks.get(source, {}).get(
                    "declared"
                ):
                    self.stdout.write(
                        f"  {borrower} no longer declares the same parameters as "
                        f"{source}; not borrowing its Qt reference"
                    )
                    continue
                hidden = set(withheld.get(borrower, ()))
                hint[borrower] = sorted(set(hint[source]) - hidden)
                self.stdout.write(
                    f"  {borrower} judged against {source}'s Qt GUI "
                    f"({len(hint[borrower])} parameters; identical def.xml"
                    + (
                        f"; {len(hidden)} withheld from it there: "
                        f"{', '.join(sorted(hidden))}"
                        if hidden
                        else ""
                    )
                    + ")"
                )

            for task, entry in tasks.items():
                if task in hint:
                    entry["qtExposed"] = hint[task]
                if task in QT_REFERENCE_ALIASES and task in hint:
                    entry["qtReferenceFrom"] = QT_REFERENCE_ALIASES[task]
                if task in auto:
                    entry["qtAutoGenerated"] = True
                if task in unreadable:
                    entry["qtExposedIncomplete"] = True
        else:
            for task, entry in tasks.items():
                for key in ("qtExposed", "qtAutoGenerated", "qtExposedIncomplete"):
                    if key in previous.get(task, {}):
                        entry[key] = previous[task][key]

        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(
            json.dumps({"tasks": tasks, "skipped": skipped}, indent=1, sort_keys=True) + "\n"
        )

        total = sum(
            sum(len(v) for v in t["declared"].values()) for t in tasks.values()
        )
        self.stdout.write(f"{len(tasks)} tasks, {total} declared parameters -> {out}")
        if skipped:
            self.stdout.write(f"{len(skipped)} skipped:")
            for task, why in sorted(skipped.items()):
                self.stdout.write(f"  {task}: {why}")
