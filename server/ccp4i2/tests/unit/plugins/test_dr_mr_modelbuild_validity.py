"""dr_mr_modelbuild_pipeline required F_SIGF, FREERFLAG and UNMERGEDFILES on
every job, though the pipeline fills F_SIGF and FREERFLAG itself and each
input type needs only its own data, so no job could pass validation. Each
type now requires what it reads, on the field the interface shows for it."""
from ccp4i2.core.tasks import get_plugin_class


def _errors(tmp_path, data_type, **files):
    plugin = get_plugin_class("dr_mr_modelbuild_pipeline")(
        workDirectory=str(tmp_path), name="dr_mr")
    plugin.container.controlParameters.MERGED_OR_UNMERGED.set(data_type)
    for name, path in files.items():
        path.write_text("")
        getattr(plugin.container.inputData, name).setFullPath(str(path))
    return [(e.get("code"), e.get("name", "")) for e in plugin.validity().getErrors()]


def test_internal_inputs_never_required(tmp_path):
    names = [name for _, name in _errors(tmp_path, "MERGED")]
    assert not any(n.endswith((".F_SIGF", ".FREERFLAG")) for n in names)


def test_each_type_requires_its_own_data(tmp_path):
    assert (203, "dr_mr_modelbuild_pipeline.container.inputData.F_SIGF_IN") in \
        _errors(tmp_path, "MERGED")
    assert (203, "dr_mr_modelbuild_pipeline.container.inputData.HKLIN") in \
        _errors(tmp_path, "MERGED_F")
    assert (203, "dr_mr_modelbuild_pipeline.container.inputData.UNMERGEDFILES") in \
        _errors(tmp_path, "UNMERGED")
    codes = [c for c, _ in _errors(tmp_path, "MERGED", F_SIGF_IN=tmp_path / "f.mtz")]
    assert 203 not in codes
