"""KPIs are labelled in words, never by their field names (issue #596).

The job list once showed ``spaceGroup``, ``highResLimit`` and ``rMeas`` beside
an aimless job's numbers. The keys stay -- they are the database's -- but what
is shown is ``KPI_LABELS``, and every performance-indicator field anywhere must
have an entry there that agrees with its ``guiLabel``.
"""

import importlib
import inspect

import pytest

from ccp4i2.core import CCP4PerformanceData
from ccp4i2.core.base_object.class_metadata import contents_from_declarations
from ccp4i2.lib.kpi_labels import KPI_LABELS, humanise_key, kpi_label, kpi_labels

# Task-local performance types, declared outside core/ (Task.dataTypes).
TASK_TYPE_MODULES = (
    "ccp4i2.wrappers.pandda_campaign.script.pandda_campaign_types",
    "ccp4i2.wrappers.pandda_events.script.pandda_events_types",
)


def _performance_classes():
    modules = [CCP4PerformanceData]
    for name in TASK_TYPE_MODULES:
        modules.append(importlib.import_module(name))
    found = {}
    for module in modules:
        for _, cls in inspect.getmembers(module, inspect.isclass):
            if issubclass(cls, CCP4PerformanceData.CPerformanceIndicator):
                found[cls.__name__] = cls
    return sorted(found.values(), key=lambda c: c.__name__)


def _declared_fields():
    """(class name, field name, guiLabel) for every declared KPI field."""
    rows = []
    for cls in _performance_classes():
        for name, declaration in contents_from_declarations(cls).items():
            rows.append((cls.__name__, name, declaration.qualifiers.get("guiLabel")))
    return rows


FIELDS = _declared_fields()


def test_the_classes_are_found():
    names = {cls for cls, _, _ in FIELDS}
    assert "CDataReductionPerformance" in names
    assert "CRefinementPerformance" in names
    assert "CPanddaRunPerformance" in names
    assert len(names) >= 15


@pytest.mark.parametrize("cls, field, gui_label", FIELDS,
                         ids=[f"{c}.{f}" for c, f, _ in FIELDS])
def test_every_field_has_a_label_that_matches_its_gui_label(cls, field, gui_label):
    assert field in KPI_LABELS, f"{cls}.{field} has no entry in KPI_LABELS"
    assert gui_label == KPI_LABELS[field], (
        f"{cls}.{field}: guiLabel {gui_label!r} != KPI_LABELS {KPI_LABELS[field]!r}"
    )


@pytest.mark.parametrize("key, label", [
    ("spaceGroup", "Space group"),
    ("highResLimit", "High resolution (Å)"),
    ("rMeas", "Rmeas"),
    ("ccHalf", "CC1/2"),
    ("RFactor", "R"),
    ("RFree", "R-free"),
    ("FOM", "FOM"),
    ("CC", "CC"),
    ("LLG", "LLG"),
    ("TFZ", "TFZ"),
])
def test_the_labels_asked_for(key, label):
    assert kpi_label(key) == label


@pytest.mark.parametrize("key, words", [
    ("highResLimit", "High res limit"),
    ("nEvents", "N events"),
    ("someNewMetric", "Some new metric"),
    ("mean_phase_error", "Mean phase error"),
    ("XMLFile", "XML file"),
    ("Hand1Score", "Hand1 score"),
    ("rmsd", "Rmsd"),
    ("TFZ", "TFZ"),
    ("", ""),
])
def test_an_unlabelled_key_is_turned_into_words(key, words):
    assert humanise_key(key) == words


def test_no_label_is_camel_case():
    for key, label in KPI_LABELS.items():
        for word in label.split():
            # A capital after a lower-case letter inside a word is camelCase;
            # symbols (CC1/2, R-free, FSC) do not have one.
            assert not any(a.islower() and b.isupper() for a, b in zip(word, word[1:])), (
                f"{key}: {label!r}"
            )


def test_kpi_labels_covers_each_key():
    assert kpi_labels(["RFree", "madeUpKey"]) == {
        "RFree": "R-free", "madeUpKey": "Made up key",
    }
