"""Pre-flight cell check for tasks that merge reflection files.

Merging observations with a FreeR set is the common case, and the two files
routinely come from different crystals -- a fragment campaign shares one free
set across a whole series, whose cells drift by a percent or more between
soaks. When the merge refuses that on a cell comparison, the failure surfaces
from inside ``merge_mtz_files`` as ``MtzMergeError: Incompatible unit cells``,
part-way through a job that has already done real work, and says nothing about
what to do next.

This puts the same question at the Run dialog instead, where it can name the
two cells, say how far apart they are, and propose the remedy.

On the remedy, and why it is not "re-cell the FreeR set with freerflag":
matching in a permissive merge is by Miller index, and a free-R flag belongs
to a reflection, not to a cell. So carrying the flags across and stamping the
observations' cell IS the re-celling operation -- there is nothing left for a
separate freerflag run to do. What the user has to decide is only whether the
two files really are the same crystal form, which no automatic check can
answer for them.
"""

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4XtalData import cells_are_compatible

#: Clipper's default, and what merge_mtz_files uses when the caller asks for a
#: strict comparison.
DEFAULT_TOLERANCE = 1.0


def _cell_of(data_file):
    """The unit cell of a CDataFile, or None if it cannot be read.

    Deliberately forgiving: this runs during validation, where a file that
    cannot be read yet is somebody else's error to report, and raising here
    would replace a useful message with a stack trace.
    """
    try:
        if not data_file.isSet():
            return None
        data_file.loadFile()
        content = data_file.getFileContent()
        cell = getattr(content, "cell", None)
        if cell is None or not cell.isSet():
            return None
        return (
            float(cell.a), float(cell.b), float(cell.c),
            float(cell.alpha), float(cell.beta), float(cell.gamma),
        )
    except Exception:
        return None


def describe_difference(first_cell, second_cell, tolerance=DEFAULT_TOLERANCE):
    """Compare two cells. Returns the result dict, or None if either is unknown."""
    if first_cell is None or second_cell is None:
        return None
    return cells_are_compatible(first_cell, second_cell, tolerance=tolerance)


def format_cell(cell):
    return "%.2f %.2f %.2f  %.1f %.1f %.1f" % cell


def check_merge_cells(
    error,
    task_name,
    observations,
    free_r,
    *,
    parameter_name="FREERFLAG_IN",
    tolerance=DEFAULT_TOLERANCE,
    code=220,
):
    """Say, before the job runs, that the free-R set will be reconciled.

    A cell difference here is not a fault: it is what a free-R set shared
    across a fragment campaign looks like, and the task handles it by running
    freerflag in COMPLETE mode to produce a set carrying this dataset's cell
    and reaching its resolution. This is an advisory so that the user knows
    that happened and can sanity-check the one thing no code can decide for
    them -- whether the two files really are the same crystal form.

    Appends to ``error`` and returns True if anything was reported.
    """
    first = _cell_of(observations)
    second = _cell_of(free_r)
    result = describe_difference(first, second, tolerance=tolerance)

    if result is None or result["validity"]:
        return False

    error.append(
        klass=task_name,
        code=code,
        details=(
            "The observations and the free-R set have different unit cells:\n"
            "  observations  %s\n"
            "  free-R set    %s\n"
            "\nThe free-R set will be reconciled with this data before "
            "refinement: freerflag joins the two by reflection index, so the "
            "existing flags stay with the reflections they were assigned to, "
            "stamps this dataset's cell, and extends the set to this data's "
            "resolution. The result is saved as this job's FREERFLAG_OUT, so "
            "later jobs on this dataset can use it directly.\n\n"
            "Worth checking only that the two files really are the same "
            "crystal form -- nothing here can tell a legitimate soak-to-soak "
            "drift from a free-R set picked from the wrong crystal."
            % (format_cell(first), format_cell(second))
        ),
        name=f"{task_name}.container.inputData.{parameter_name}",
        severity=CCP4ErrorHandling.SEVERITY_WARNING,
    )
    return True
