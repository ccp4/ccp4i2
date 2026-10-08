"""What a KPI is called on screen.

A job's key performance indicators are stored under the name of the field
they came from (``JobValueKey.name``: ``spaceGroup``, ``highResLimit``,
``RFree``). Those names are the database's and the params files' and must not
change; but they are not words, and a job list that reads ``highResLimit: 1.8``
is showing a user the plumbing (issue #596).

``KPI_LABELS`` is the one mapping from key to label. The key namespace is
global -- every task that reports ``RFree`` stores it under the same
``JobValueKey`` -- so a label belongs to the name, not to the class that
declared it, and each name has exactly one. The ``guiLabel`` of every
performance-indicator field is held equal to its entry here by a unit test,
so the task interface and the job list say the same thing.

Labels are sentence case and keep the field's conventional symbols as
written (Rmeas, R-free, CC1/2, FOM); a resolution or a distance carries its
unit, so the value beside it needs none.

A key with no entry -- a new field nobody has labelled yet, or a key imported
from an older database -- still never shows as camelCase: ``kpi_label``
turns it into words (``highResLimit`` -> ``High res limit``). The client
carries the same fallback (``client/renderer/lib/format-kpi.ts``).
"""

from __future__ import annotations

import re
from typing import Dict, Iterable

KPI_LABELS: Dict[str, str] = {
    # The base indicator's own fields
    "value": "Value",
    "annotation": "Annotation",
    # Data reduction
    "spaceGroup": "Space group",
    "highResLimit": "High resolution (Å)",
    "rMeas": "Rmeas",
    "ccHalf": "CC1/2",
    # Paired refinement
    "cutoff": "Resolution cutoff (Å)",
    # Refinement and model fit
    "RFactor": "R",
    "RFree": "R-free",
    "R": "R (all)",
    "R1Factor": "R1",
    "R1Free": "R1-free",
    "R1": "R1 (all)",
    "RMSBond": "RMS bond (Å)",
    "RMSAngle": "RMS angle (°)",
    "weightUsed": "Weight",
    "FSCaverage": "FSC average",
    "CCFwork_avg": "CCF work",
    "CCFfree_avg": "CCF free",
    "CCF_avg": "CCF",
    "CCIwork_avg": "CCI work",
    "CCIfree_avg": "CCI free",
    "CCI_avg": "CCI",
    # Experimental phasing and molecular replacement
    "FOM": "FOM",
    "CFOM": "Combined FOM",
    "Hand1Score": "Hand 1 score",
    "Hand2Score": "Hand 2 score",
    "CC": "CC",
    "LLG": "LLG",
    "TFZ": "TFZ",
    "phaseError": "Phase error (°)",
    "weightedPhaseError": "Weighted phase error (°)",
    "reflectionCorrelation": "Reflection correlation",
    # Model building and model preparation
    "completeness": "Completeness",
    "nAtoms": "Atoms",
    "nResidues": "Residues",
    # Superposition
    "RMSxyz": "RMS deviation (Å)",
    "QScore": "Q-score",
    # Format conversion tests
    "columnLabelsString": "Columns",
    # PanDDA (types in wrappers/pandda_campaign and wrappers/pandda_events)
    "nDatasets": "Datasets staged",
    "nDatasetsProcessed": "Datasets processed",
    "nDatasetsAnalysed": "Datasets analysed",
    "nEvents": "Events",
    "nSites": "Sites",
    "wallSeconds": "Wall time (s)",
    "nCreated": "Receipts created",
    "nSkipped": "Already had a receipt",
    "nAbsent": "Not in the tree",
    "nFailed": "Failed",
    "nNoProject": "Project not in this database",
    "nEventsExpected": "Events expected",
    "nEventsDelivered": "Event maps delivered",
    "nPosesExpected": "Poses expected",
    "nPosesDelivered": "Poses delivered",
    "bestBuildScore": "Best build score",
    "STATE": "State",
    "STDERR": "stderr",
}

# A lower-case letter or digit followed by a capital; a run of capitals
# followed by a capitalised word (the run is an acronym: "XMLFile").
_WORD_BOUNDARY = re.compile(r"(?<=[a-z0-9])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])")


def humanise_key(key: str) -> str:
    """Turn a field name into words, for a key with no label.

    ``highResLimit`` -> ``High res limit``; ``nEvents`` -> ``N events``;
    ``mean_phase_error`` -> ``Mean phase error``. Words in capitals are taken to be
    symbols and left as they are; other words after the first are lower-cased.
    Must agree with ``humaniseKpiKey`` in ``client/renderer/lib/format-kpi.ts``.
    """
    words = [w for w in _WORD_BOUNDARY.sub(" ", key.replace("_", " ")).split() if w]
    if not words:
        return key
    out = []
    for i, word in enumerate(words):
        if word[1:].lower() == word[1:] and not word.isupper():
            word = word.lower()
        out.append(word)
    out[0] = out[0][0].upper() + out[0][1:]
    return " ".join(out)


def kpi_label(key: str) -> str:
    """The on-screen label for a KPI key: its entry, or the key in words."""
    return KPI_LABELS.get(key) or humanise_key(key)


def kpi_labels(keys: Iterable[str]) -> Dict[str, str]:
    """``{key: label}`` for each of `keys`, for an API response to carry."""
    return {key: kpi_label(key) for key in keys}
