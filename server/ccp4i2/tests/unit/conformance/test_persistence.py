"""
A job's parameters must survive the round trip, and outlive a container change.

Tier 4. A params.xml is not written and read by the same tree: it is written
by whatever the code was then, and read back into a tree built by the code
now, from a def.xml that may since have changed. Defect A already relied on
that --- old phaser params files carry 29 ghost sections that must be dropped
without disturbing the real keywords.* values.

Two directions, and the second is the one that protects users:

  forwards   written today, read into a tree built today, values identical
  backwards  written before a container change, still readable afterwards,
             unknown elements dropped and known ones intact

Value *state* is checked alongside value, because isSet() decides whether a
parameter reaches a program at all and 719 uses depend on it. A round trip
that restores values but flattens EXPLICITLY_SET into DEFAULT would silently
change which keywords a re-run sends.

Pure Python -- no CCP4 binaries needed.
"""
import xml.etree.ElementTree as ET

import pytest

from ccp4i2.tests.unit.conformance import harness

# Tasks chosen for shape rather than importance: a plain wrapper, one with
# lists, one with deep keyword nesting, one that inherits another's def.xml.
SAMPLE = ['freerflag', 'gesamt', 'phaser_MR', 'servalcat_pipe', 'aimless_pipe']


def _build(task):
    plugin = harness.build(task)
    if plugin is None:
        pytest.skip(f'{task} is not available in this environment')
    return plugin


def _values(plugin):
    """path -> (str value, isSet) for every leaf, as the tree sees it now."""
    out = {}
    for path, obj in harness.walk(plugin, 'x').leaves:
        try:
            is_set = obj.isSet() if callable(getattr(obj, 'isSet', None)) else None
        except Exception:
            is_set = 'raised'
        try:
            out[path] = (str(obj), is_set)
        except Exception:
            out[path] = ('<unstringifiable>', is_set)
    return out


@pytest.mark.parametrize('task', SAMPLE)
def test_a_tree_round_trips_through_its_own_xml(task):
    written = _build(task)
    before = _values(written)

    etree = written.container.getEtree(excludeUnset=True)
    read = _build(task)
    read.container.setEtree(etree, ignore_missing=True)

    after = _values(read)
    differing = {p: (before[p], after.get(p)) for p in before
                 if after.get(p) != before[p]}
    assert not differing, (
        f'{len(differing)} parameters changed across a round trip:\n  '
        + '\n  '.join(f'{p}: {b} -> {a}' for p, (b, a) in list(differing.items())[:10]))


@pytest.mark.parametrize('task', SAMPLE)
def test_the_round_trip_reaches_the_same_parameters(task):
    written = _build(task)
    etree = written.container.getEtree(excludeUnset=True)
    read = _build(task)
    read.container.setEtree(etree, ignore_missing=True)
    assert set(_values(read)) == set(_values(written))


# --- value state, not merely value ------------------------------------------

def test_an_explicitly_set_value_is_still_explicitly_set_afterwards():
    """A value equal to the default, chosen deliberately, must stay chosen.

    phaser skips keywords that are at their default: "don't set keywords set
    to their default values". So a user who types the default and a user who
    leaves it alone are meant to differ, and only the state records that.
    """
    plugin = _build('freerflag')
    target = plugin.container.controlParameters
    name = next(n for n in target.dataOrder())
    obj = getattr(target, name)
    obj.set(str(obj))                       # set it to what it already shows
    assert obj.isSet()

    etree = plugin.container.getEtree(excludeUnset=True)
    read = _build('freerflag')
    read.container.setEtree(etree, ignore_missing=True)

    assert getattr(read.container.controlParameters, name).isSet(), \
        'an explicit choice was flattened back into a default by the round trip'


def test_an_untouched_value_does_not_become_explicitly_set():
    """The other direction: loading must not invent choices the user never made."""
    plugin = _build('freerflag')
    before = {p: v[1] for p, v in _values(plugin).items()}

    read = _build('freerflag')
    read.container.setEtree(plugin.container.getEtree(excludeUnset=True),
                            ignore_missing=True)
    after = {p: v[1] for p, v in _values(read).items()}

    became_set = [p for p in before if not before[p] and after.get(p)]
    assert not became_set, \
        f'{len(became_set)} parameters became set by being loaded: {became_set[:8]}'


# --- backwards compatibility: a params.xml older than the container ---------

def test_a_params_file_naming_things_the_container_lost_still_loads():
    """Defect A removed 1,358 ghost paths; jobs written before it still exist."""
    plugin = _build('freerflag')
    etree = plugin.container.getEtree(excludeUnset=True)
    ghost = ET.SubElement(etree, 'aSectionThatNoLongerExists')
    ET.SubElement(ghost, 'SOMEPARAM').text = 'value'

    read = _build('freerflag')
    read.container.setEtree(etree, ignore_missing=True)      # must not raise

    assert not hasattr(read.container, 'aSectionThatNoLongerExists')


def test_the_real_values_survive_alongside_the_dropped_ones():
    plugin = _build('freerflag')
    plugin.container.inputData.F_SIGF.setFullPath('/tmp/some_reflections.mtz')
    etree = plugin.container.getEtree(excludeUnset=True)
    ET.SubElement(etree, 'aSectionThatNoLongerExists').text = 'junk'

    read = _build('freerflag')
    read.container.setEtree(etree, ignore_missing=True)

    assert 'some_reflections.mtz' in str(read.container.inputData.F_SIGF)


def test_a_params_file_missing_a_newly_added_parameter_still_loads():
    """The other half of backwards compatibility: the container gained a field."""
    plugin = _build('freerflag')
    etree = plugin.container.getEtree(excludeUnset=True)
    section = etree.find('inputData')
    if section is not None and len(section):
        section.remove(section[0])           # an old file that predates a field

    read = _build('freerflag')
    read.container.setEtree(etree, ignore_missing=True)      # must not raise
    assert read.container.inputData is not None


# --- a property worth pinning, not a defect ---------------------------------

def test_serialising_everything_marks_everything_explicit():
    """`getEtree()` without excludeUnset writes unset values, and loading them
    back marks them EXPLICITLY_SET --- five parameters in freerflag alone go
    from unset to set, with the same value.

    Production never does this: saveParams and saveDataToXml both default to
    exclude_unset=True, so a params.xml holds only what was chosen. The pin is
    here because the reevaluation could easily change which of the two is the
    default, and the difference is invisible until something asks isSet().
    """
    plugin = _build('freerflag')
    before = {p: v[1] for p, v in _values(plugin).items()}

    read = _build('freerflag')
    read.container.setEtree(plugin.container.getEtree(), ignore_missing=True)
    after = {p: v[1] for p, v in _values(read).items()}

    became_set = [p for p in before if not before[p] and after.get(p)]
    assert became_set, \
        'a full serialisation no longer marks unset values as set -- good, ' \
        'but the reevaluation should say so deliberately'


# --- the two ways a value gets set must agree -------------------------------

def test_the_gui_route_and_the_python_route_reach_the_same_state():
    """A job built in the GUI and the same job built by i2run must agree.

    The GUI goes through CContainer.set_parameter(); a plugin or i2run calls
    .set() directly. If one marked EXPLICITLY_SET and the other DEFAULT for
    the same visible value, the two jobs would send different keywords to the
    program while looking identical in every view we have --- phaser skips
    keywords that are at their default.
    """
    by_python = _build('freerflag')
    section = by_python.container.controlParameters
    name = next(n for n in section.dataOrder())
    getattr(section, name).set('1')

    by_gui = _build('freerflag')
    by_gui.container.set_parameter(f'controlParameters.{name}', '1', skip_first=False)

    python_obj = getattr(by_python.container.controlParameters, name)
    gui_obj = getattr(by_gui.container.controlParameters, name)

    assert str(python_obj) == str(gui_obj)
    assert python_obj.isSet() == gui_obj.isSet() is True
    assert python_obj.getValueState() == gui_obj.getValueState()
