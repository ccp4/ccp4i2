import pathlib
import pandas as pd
import numpy
import gemmi


def analyze_mtz(file_path, with_reflections: bool = False):
    mtz = gemmi.read_mtz_file(file_path)
    result = {}

    result["header_info"] = {
        "title": mtz.title,
        "spacegroup": mtz.spacegroup.hm,
        "sort_order": mtz.sort_order,
        "history": mtz.history,
        "resolution_low": mtz.resolution_low(),
        "resolution_high": mtz.resolution_high(),
        "cell": {
            "a": mtz.cell.a,
            "b": mtz.cell.b,
            "c": mtz.cell.c,
            "alpha": mtz.cell.alpha,
            "beta": mtz.cell.beta,
            "gamma": mtz.cell.gamma,
        },
        "datasets": [
            {
                "id": dataset.id,
                "project_name": dataset.project_name,
                "crystal_name": dataset.crystal_name,
            }
            for dataset in mtz.datasets
        ],
        "columns": [
            {
                "label": col.label,
                "type": col.type,
                "min": col.min_value,
                "max": col.max_value,
                "dataset_id": col.source,
            }
            for col in mtz.columns
        ],
    }

    binner = gemmi.Binner()
    binner.setup(20, gemmi.Binner.Method.EqualCount, mtz)

    reflection_counts = []
    for i in range(binner.size):
        limits = {
            "number": i + 1,
            "dmax": binner.dmax_of_bin(i),
            "dmin": binner.dmin_of_bin(i),
        }
        reflection_counts.append(limits)

    uniques = gemmi.make_miller_array(
        mtz.cell, mtz.spacegroup, mtz.resolution_high(), mtz.resolution_low()
    )
    unique_counts = numpy.bincount(binner.get_bins(uniques))
    for i in range(binner.size):
        reflection_counts[i]["Unique"] = unique_counts[i]

    present_counts = numpy.bincount(binner.get_bins(mtz))
    for i in range(binner.size):
        reflection_counts[i]["Present"] = present_counts[i]

    result["counts"] = reflection_counts

    if with_reflections:
        df = pd.DataFrame(data=mtz.array, columns=mtz.column_labels())
        result["reflections"] = df.to_dict(orient="dict")
    return result


def mtz_as_dict(file_path, with_reflections: bool = False):
    return analyze_mtz(str(pathlib.Path(file_path)), with_reflections)
