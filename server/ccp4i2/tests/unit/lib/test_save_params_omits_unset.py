"""A parameter the user never touched must not be written as though they had.

The regression this pins: `save_params_for_job` used to default to
`exclude_unset=False`, and only job creation overrode it. Creating an
aimless_pipe job therefore wrote a clean input_params.xml, and then *uploading
a file* --- which rewrites the same file through the same function, without
the override --- added every untouched integer at its zero value:

    <NPROC>0</NPROC>  <SCALES_NTILEX>0</SCALES_NTILEX>  ...

Those read back as explicitly set. NPROC declares min=1, so run-time
validation then blocked the job on a value the user had never chosen.

Pure Python --- no CCP4 binaries, no database (the function touches only a
handful of attributes on the job, so a stand-in carries them).
"""
import types
import uuid
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.tests.unit.conformance import harness

save_params = pytest.importorskip(
    'ccp4i2.lib.utils.parameters.save_params',
    reason='needs the Django app importable',
)

# aimless_pipe is where this was found; NPROC is the parameter that bit.
TASK = 'aimless_pipe'
BOUNDED = 'NPROC'


def _plugin():
    plugin = harness.build(TASK)
    if plugin is None:
        pytest.skip(f'{TASK} is not available in this environment')
    return plugin


def _stub_job(tmp_path):
    """Just the attributes save_params_for_job reads off a Job."""
    return types.SimpleNamespace(
        directory=tmp_path,
        uuid=uuid.uuid4(),
        number='1',
        project=types.SimpleNamespace(uuid=uuid.uuid4(), name='a_project'),
    )


def _written(tmp_path, plugin, **kwargs):
    save_params.save_params_for_job(plugin, _stub_job(tmp_path), **kwargs)
    written = list(tmp_path.glob('*.xml'))
    assert written, 'save_params_for_job wrote nothing'
    return ET.parse(str(written[0])).getroot()


def _named(root, name):
    return [e for e in root.iter() if e.tag == name]


def test_an_untouched_bounded_parameter_is_not_written(tmp_path):
    root = _written(tmp_path, _plugin())
    found = _named(root, BOUNDED)
    assert not found, (
        f'{BOUNDED} was written as {found[0].text!r} without anyone setting it; '
        'reloading that blocks the job on its own minimum'
    )


def test_nothing_untouched_is_written(tmp_path):
    """Not just the one that had a lower bound to trip over."""
    plugin = _plugin()
    unset = {
        path.rsplit('.', 1)[-1]
        for path, obj in harness.walk(plugin, TASK).leaves
        if callable(getattr(obj, 'isSet', None)) and not obj.isSet()
    }
    root = _written(tmp_path, plugin)
    leaked = sorted({e.tag for e in root.iter()} & unset)
    assert not leaked, f'{len(leaked)} unset parameters were written: {leaked[:8]}'


def test_what_was_actually_chosen_still_survives(tmp_path):
    """The fix must not cost the file the values that belong in it."""
    plugin = _plugin()
    plugin.container.controlParameters.NPROC.set(4)

    root = _written(tmp_path, plugin)
    found = _named(root, BOUNDED)
    assert found and found[0].text.strip() == '4', \
        'an explicitly chosen value did not reach the file'


def test_the_written_file_reloads_without_inventing_choices(tmp_path):
    """The end of the user-visible failure: reload and ask isSet()."""
    root = _written(tmp_path, _plugin())

    read = _plugin()
    read.container.setEtree(root, ignore_missing=True)
    assert not read.container.controlParameters.NPROC.isSet(), \
        f'{BOUNDED} came back set from a file nobody set it in'


def test_a_caller_can_still_ask_for_everything(tmp_path):
    """exclude_unset=False stays available for a deliberate full record."""
    root = _written(tmp_path, _plugin(), exclude_unset=False)
    assert _named(root, BOUNDED), 'exclude_unset=False no longer writes unset values'
