"""
A nested container belongs where the def.xml put it, and nowhere else.

The parser attached every container it could find under `<ccp4i2_body>` to the
root a second time. What it tracked as "already placed" was a set of *id
strings* of the body's *direct* children, so anything nested deeper --- the 29
keyword sub-containers of every phaser task, servalcat_pipe's prosmart ones ---
failed that test and was re-parsed onto the root.

The result was a parallel ghost tree of 2,304 parameters across 15 tasks that
nothing reads. i2run's bare keyword flags addressed the ghost, so they set a
value no program ever saw: silent no-ops. The ghosts were also serialised into
params.xml.

The trap, which cost a first attempt at this fix: five def.xml files open with
a `<ccp4i2_body>` carrying no id, and for those the root-level loop is the only
thing that builds a container at all. Deleting the loop strips inputData,
outputData and controlParameters from them entirely.

Pure Python -- no CCP4 binaries needed.
"""
import tempfile
from pathlib import Path

import pytest

from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

NESTED = """<?xml version='1.0' encoding='ASCII'?>
<ccp4:ccp4i2 xmlns:ccp4="http://www.ccp4.ac.uk/ccp4ns">
    <ccp4i2_header><function>DEF</function><pluginName>toy</pluginName></ccp4i2_header>
    <ccp4i2_body id="toy">
        <container id="inputData">
            <content id="XYZIN"><className>CPdbDataFile</className><qualifiers/></content>
        </container>
        <container id="keywords">
            <container id="tncs">
                <content id="RESTRICT"><className>CString</className><qualifiers/></content>
            </container>
        </container>
    </ccp4i2_body>
</ccp4:ccp4i2>
"""

ID_LESS_BODY = NESTED.replace('<ccp4i2_body id="toy">', '<ccp4i2_body>')


def _parse(text):
    directory = Path(tempfile.mkdtemp())
    path = directory / "toy.def.xml"
    path.write_text(text, encoding="utf-8")
    return DefXmlParser().parse_def_xml(path)


def test_a_nested_container_is_not_also_at_the_root():
    container = _parse(NESTED)
    top = list(container.dataOrder())
    assert "keywords" in top
    assert "tncs" not in top, "the nested container was copied to the root"


def test_the_nested_container_is_still_where_it_belongs():
    container = _parse(NESTED)
    assert "tncs" in list(container.keywords.dataOrder())
    assert "RESTRICT" in list(container.keywords.tncs.dataOrder())


def test_a_body_with_no_id_still_yields_its_containers():
    """The five-file trap: without this the task loses everything."""
    container = _parse(ID_LESS_BODY)
    top = list(container.dataOrder())
    assert "inputData" in top
    assert "keywords" in top
    assert "tncs" not in top


# --- the real tasks ---------------------------------------------------------

def _plugin(task):
    """The plugin, not just its container.

    A container's children are reached through weak references, so a container
    outliving its plugin comes back empty. Returning the plugin keeps it alive
    for as long as the caller needs it.
    """
    from ccp4i2.core.tasks import get_plugin_class

    plugin_class = get_plugin_class(task)
    if plugin_class is None:
        pytest.skip(f"{task} is not available in this environment")
    return plugin_class(workDirectory=tempfile.mkdtemp(), name="t")


@pytest.mark.parametrize("task", [
    "phaser_MR", "phaser_MR_AUTO", "phaser_pipeline", "phaser_rnp_pipeline",
    "servalcat_pipe",
])
def test_a_real_task_has_no_ghost_sections(task):
    """No id that the def.xml uses *only* below body level may sit at the root.

    Not "no name appears twice": servalcat_pipe's sub-task containers each
    legitimately hold their own inputData, and inputData is a body-level
    section too. What must never appear at the root is an id the def.xml only
    ever nests --- phaser's `tncs`, `solparameters`, `translation`.
    """
    import xml.etree.ElementTree as ET

    from ccp4i2.core.tasks import locate_def_xml
    from ccp4i2.core.task_manager.def_xml_handler import load_nested_xml

    plugin = _plugin(task)
    root = load_nested_xml(ET.parse(locate_def_xml(task)).getroot())
    body = root.find(".//ccp4i2_body[@id]")
    if body is None:
        body = root.find(".//ccp4i2_body")
    nested_body = body.find("./ccp4i2_body")
    if nested_body is not None:
        body = nested_body

    at_body_level = {c.get("id") for c in body.findall("./container[@id]")}
    anywhere = {c.get("id") for c in body.findall(".//container[@id]")}
    only_nested = anywhere - at_body_level

    ghosts = only_nested & set(plugin.container.dataOrder())
    assert not ghosts, \
        f"{task} carries {len(ghosts)} nested sections at its root: {sorted(ghosts)[:6]}"


def test_the_sections_a_task_declares_are_all_there():
    """Removing the ghosts must not remove anything real."""
    import xml.etree.ElementTree as ET

    from ccp4i2.core.tasks import locate_def_xml
    from ccp4i2.core.task_manager.def_xml_handler import load_nested_xml

    plugin = _plugin("phaser_pipeline")
    root = load_nested_xml(ET.parse(locate_def_xml("phaser_pipeline")).getroot())
    body = root.find(".//ccp4i2_body[@id]")
    declared = {c.get("id") for c in body.findall("./container[@id]")}

    assert declared <= set(plugin.container.dataOrder())


@pytest.mark.parametrize("task", [
    "cmapcoeff", "comit", "cpatterson", "fft", "validate_protein",
])
def test_the_id_less_body_tasks_keep_their_containers(task):
    plugin = _plugin(task)
    top = set(plugin.container.dataOrder())
    assert {"inputData", "outputData", "controlParameters"} <= top, \
        f"{task} lost sections: has {sorted(top)}"
