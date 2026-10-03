"""import_merged does not demand crystal and dataset names: they come from
the file. Every agent in the 2026-10-03 trials was stopped by "Required
value not set" on both and had to invent names; the interface fills them
from the file, a job set up any other way did not."""
from ccp4i2.core.tasks import get_plugin_class


def test_names_are_not_required(tmp_path):
    plugin = get_plugin_class("import_merged")(workDirectory=str(tmp_path), name="imp")
    named = [str(r.get("name", "")) for r in plugin.validity()._reports]
    assert not any(n.endswith(".CRYSTALNAME") or n.endswith(".DATASETNAME") for n in named)
