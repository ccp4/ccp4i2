"""
An ensemble list that asks phaser to place nothing.

A job was submitted with every ensemble asking for 0 copies: phaser ran and
found nothing, because there was nothing to look for. The copies count now
defaults to 1, an implausible count on one ensemble is a warning, a list in
which no ensemble asks for copies is a warning while editing, and the tasks
that search refuse to run in that state.
"""

import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4ModelData import CEnsemble, CEnsembleList
from ccp4i2.core.tasks import get_plugin_class
from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR.script import phaser_MR

WARNING = CCP4ErrorHandling.SEVERITY_WARNING
ERROR = CCP4ErrorHandling.SEVERITY_ERROR


def _reports(report, code=None, name_ends=None, klass=None):
    out = []
    for err in report.getErrors():
        if code is not None and err.get('code') != code:
            continue
        if klass is not None and err.get('class') != klass:
            continue
        if name_ends is not None and not err.get('name', '').endswith(name_ends):
            continue
        out.append(err)
    return out


def _task_with_ensembles(name, copies):
    """A task whose ENSEMBLES holds one ensemble per entry in copies.

    An entry may be an int (copies, in use) or a (copies, use) pair.
    """
    task = get_plugin_class(name)()
    ensembles = task.container.inputData.ENSEMBLES
    for entry in copies:
        n, use = entry if isinstance(entry, tuple) else (entry, True)
        item = ensembles.makeItem()
        item.number.set(n)
        item.use.set(use)
        ensembles.append(item)
    return task


# ---------------------------------------------------------------- default

def test_a_new_ensemble_asks_for_one_copy_and_is_in_use():
    ensembles = CEnsembleList(name='ENSEMBLES')
    item = ensembles.makeItem()
    assert int(item.number) == 1
    assert item.isUsed()
    assert item.copiesToPlace() == 1


def test_the_default_survives_the_json_template_the_gui_builds_rows_from():
    """The GUI builds a new row from the encoder's _subItem: its copies must be 1."""
    import json
    from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder
    ensembles = CEnsembleList(name='ENSEMBLES')
    encoded = json.loads(json.dumps(ensembles, cls=CCP4i2JsonEncoder))
    assert encoded['_subItem']['_value']['number']['_value'] == 1


# ------------------------------------------------- one ensemble's copies

@pytest.mark.parametrize("copies", [0, 13, 40])
def test_an_implausible_copy_count_is_a_warning(copies):
    ensemble = CEnsemble(name='ens')
    ensemble.number.set(copies)
    found = _reports(ensemble.validity(), code=101, klass='CEnsemble')
    assert len(found) == 1
    assert found[0]['severity'] == WARNING
    assert found[0]['class'] == 'CEnsemble'


@pytest.mark.parametrize("copies", [1, 2, 12])
def test_a_plausible_copy_count_is_not_remarked_on(copies):
    ensemble = CEnsemble(name='ens')
    ensemble.number.set(copies)
    assert _reports(ensemble.validity(), code=101, klass='CEnsemble') == []


# ------------------------------------------------------------ use flag

def test_an_ensemble_not_in_use_places_nothing_whatever_its_count():
    ensemble = CEnsemble(name='ens')
    ensemble.number.set(3)
    ensemble.use.set(False)
    assert not ensemble.isUsed()
    assert ensemble.copiesToPlace() == 0


def test_zero_copies_on_an_unused_ensemble_is_not_a_warning():
    ensemble = CEnsemble(name='ens')
    ensemble.number.set(0)
    ensemble.use.set(False)
    assert _reports(ensemble.validity(), code=101, klass='CEnsemble') == []


def test_an_unset_use_flag_counts_as_in_use():
    ensemble = CEnsemble(name='ens')
    ensemble.number.set(2)
    ensemble.use.unSet()
    assert ensemble.isUsed()
    assert ensemble.copiesToPlace() == 2


# --------------------------------------------------------- the whole list

def test_a_list_of_zeros_is_a_warning_while_editing():
    task = _task_with_ensembles('phaser_pipeline', [0, 0])
    report = task.validity()
    found = _reports(report, code=103)
    assert len(found) == 1
    assert found[0]['severity'] == WARNING
    assert found[0]['name'] == 'phaser_pipeline.container.inputData.ENSEMBLES'
    # ...and each zero is flagged on its own row too
    assert len(_reports(report, code=101, klass='CEnsemble')) == 2


def test_one_ensemble_with_copies_is_enough():
    task = _task_with_ensembles('phaser_pipeline', [0, 2])
    assert _reports(task.validity(), code=103) == []
    assert task.container.inputData.ENSEMBLES.copiesToPlace() == 2


def test_copies_on_an_ensemble_not_in_use_do_not_count():
    task = _task_with_ensembles('phaser_pipeline', [(0, True), (2, False)])
    report = task.validity()
    assert len(_reports(report, code=103)) == 1
    # only the zero on the ensemble in use is flagged on its own row
    assert len(_reports(report, code=101, klass='CEnsemble')) == 1
    assert task.container.inputData.ENSEMBLES.copiesToPlace() == 0


def test_an_empty_list_is_reported_as_empty_not_as_zeros():
    task = get_plugin_class('phaser_pipeline')()
    report = task.validity()
    assert _reports(report, code=103) == []
    assert _reports(report, code=101, klass='CEnsembleList', name_ends='ENSEMBLES')  # listMinLength


# ------------------------------------------------------------- run time

def _run_time_ensemble_errors(task):
    return _reports(task.runTimeValidity(), code=111, name_ends='.ENSEMBLES')


@pytest.mark.parametrize("name", ["phaser_pipeline", "phaser_MR_AUTO",
                                  "phaser_MR_FRF", "phaser_MR_FTF"])
def test_a_searching_task_refuses_to_run_with_nothing_to_place(name):
    task = _task_with_ensembles(name, [0, 0])
    found = _run_time_ensemble_errors(task)
    assert len(found) == 1
    assert found[0]['severity'] == ERROR
    assert found[0]['name'] == f'{name}.container.inputData.ENSEMBLES'
    # while editing it is only a warning
    assert _reports(task.validity(), code=111) == []


@pytest.mark.parametrize("name", ["phaser_pipeline", "phaser_MR_AUTO"])
def test_a_searching_task_runs_when_one_ensemble_has_copies(name):
    task = _task_with_ensembles(name, [0, 1])
    assert _run_time_ensemble_errors(task) == []


@pytest.mark.parametrize("name", ["phaser_pipeline", "phaser_MR_AUTO"])
def test_a_searching_task_refuses_when_the_only_copies_are_switched_off(name):
    task = _task_with_ensembles(name, [(0, True), (3, False)])
    assert len(_run_time_ensemble_errors(task)) == 1


def test_the_search_skips_ensembles_not_in_use():
    """addSearches hands phaser only the ensembles in use with copies."""
    task = _task_with_ensembles('phaser_MR_AUTO', [(2, True), (3, False), (0, True)])
    for i, ens in enumerate(task.container.inputData.ENSEMBLES):
        ens.label.set(f'E{i}')
    searches = []

    class FakeInput:
        def addSEAR_ENSE_NUM(self, label, n):
            searches.append((label, n))

    task.addSearches(FakeInput())
    assert searches == [('E0', 2)]


@pytest.mark.parametrize("name", ["phaser_MR_PAK", "phaser_MR_RNP"])
def test_modes_that_place_nothing_do_not_care(name):
    task = _task_with_ensembles(name, [0])
    assert _run_time_ensemble_errors(task) == []


@pytest.mark.parametrize("mode", ["MR_PAK", "MR_RNP"])
def test_the_pipeline_follows_its_mode(mode):
    task = _task_with_ensembles('phaser_pipeline', [0])
    task.container.inputData.MODE_TY.set(mode)
    assert _run_time_ensemble_errors(task) == []


def test_the_helper_leaves_an_empty_list_alone():
    """phaser_simple fills ENSEMBLES in process(); an empty list is not a zero list."""
    task = get_plugin_class('phaser_pipeline')()
    error = CCP4ErrorHandling.CErrorReport()
    phaser_MR.ensembles_run_time_validity('phaser_pipeline', task.container.inputData, error)
    assert len(error.getErrors()) == 0
