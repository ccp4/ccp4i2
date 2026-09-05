"""A job saved before a .multiple scope became a CList still loads.

Before, a repeated PHIL scope was one container and a repeated definition
one leaf, and params.xml recorded them that way. Such an element now lands
on a CList and must be read as one item -- not dropped, and not misread as
one garbage item per field.
"""
import os
import tempfile

import pytest

parse = pytest.importorskip("libtbx.phil", reason="needs libtbx (CCP4/cctbx)").parse

from ccp4i2.core.PhilPluginScript import PhilPluginScript

PHIL = parse("""
    composition {
      solvent = None
        .type = float
      chain
        .optional = True
        .multiple = True
      {
        nres = None
          .type = int
        num = 1
          .type = int
      }
    }
    model = None
      .type = path
      .multiple = True
""")


class OldJobPlugin(PhilPluginScript):
    TASKNAME = "old_job_plugin"
    TASKCOMMAND = "echo"

    def get_master_phil(self):
        return PHIL

    def get_phil_exclude_scopes(self):
        return []

    def get_command_target(self):
        return "echo"


OLD_PARAMS = """<?xml version='1.0' encoding='utf-8'?>
<ccp4i2><ccp4i2_header><function>PARAMS</function><pluginName>old_job_plugin</pluginName></ccp4i2_header>
<ccp4i2_body><inputData/><outputData/><controlParameters>
  <composition>
    <composition__chain>
      <composition__chain__nres>120</composition__chain__nres>
      <composition__chain__num>2</composition__chain__num>
    </composition__chain>
  </composition>
  <model>/x/a.pdb</model>
</controlParameters></ccp4i2_body></ccp4i2>"""


def test_a_single_scope_and_a_single_leaf_load_as_one_item_each():
    with tempfile.TemporaryDirectory() as tmp:
        plugin = OldJobPlugin(workDirectory=tmp)
        path = os.path.join(tmp, "params.xml")
        with open(path, "w") as f:
            f.write(OLD_PARAMS)
        plugin.loadDataFromXml(path)
        chain = plugin.container.controlParameters.composition.composition__chain
        assert len(chain) == 1
        assert chain[0].composition__chain__nres.value == 120
        assert chain[0].composition__chain__num.value == 2
        model = plugin.container.controlParameters.model
        assert len(model) == 1 and str(model[0]) == "/x/a.pdb"
        assert plugin.extract_phil_lines() == [
            "composition.chain {", "  nres = 120", "  num = 2", "}", "model = /x/a.pdb"]


def test_the_new_format_still_round_trips():
    with tempfile.TemporaryDirectory() as tmp:
        plugin = OldJobPlugin(workDirectory=tmp)
        chain = plugin.container.controlParameters.composition.composition__chain
        for nres in (10, 20):
            chain.append(chain.makeItem())
            chain[-1].composition__chain__nres.value = nres
        plugin.container.controlParameters.model.append("/x/a.pdb")
        path = os.path.join(tmp, "params.xml")
        plugin.saveDataToXml(path)
        again = OldJobPlugin(workDirectory=tmp)
        again.loadDataFromXml(path)
        chain2 = again.container.controlParameters.composition.composition__chain
        assert [c.composition__chain__nres.value for c in chain2] == [10, 20]
        assert [str(m) for m in again.container.controlParameters.model] == ["/x/a.pdb"]
