import pathlib
import logging
import gemmi
from typing import List
import re
from ccp4i2.core import CCP4XtalData
from ..files.available_name import available_file_name_based_on

logger = logging.getLogger(f"ccp4x:{__name__}")


def gemmi_split_mtz(
    input_file_path: pathlib.Path = None,
    input_column_path: str = None,
    preferred_dest: pathlib.Path = None,
    return_metadata: bool = False,
):
    """
    Splits an MTZ file based on specified columns and writes the result to a new file.

    Args:
        input_file_path (pathlib.Path): Path to the input MTZ file. Must be provided and must exist.
        input_column_path (str): Path to the columns to be extracted, e.g., '/*/*/[F,SIGFP]'. Must be provided.
        preferred_dest (pathlib.Path): Preferred destination path for the split file. Must be provided.
        return_metadata (bool): If True, return a dict with path and metadata instead of just the path.

    Raises:
        Exception: If any of the required arguments are not provided or if the input file does not exist.
        Exception: If the input column path is invalid or if the specified columns cannot be found in the input file.

    Returns:
        pathlib.Path: Path to the newly created MTZ file with the selected columns.
        OR (if return_metadata=True):
        dict: {"path": pathlib.Path, "metadata": {"crystal_name": str, "dataset_name": str, "column_labels": list}}
    """
    logger.warning("Starting gemmi_split_mtz with input_file_path=%s, input_column_path=%s, preferred_dest=%s",
                   input_file_path, input_column_path, preferred_dest)  
    
    if input_file_path is None:
        raise Exception("smartSplitMTZ Exception:", "Must provide an input file")
    if not input_file_path.is_file():
        raise Exception(
            "smartSplitMTZ Exception:", "inputFile must exist" + str(input_file_path)
        )
        
    logger.warning("Input file exists: %s", input_file_path)
    if input_column_path is None:
        raise Exception(
            "smartSplitMTZ Exception:",
            "Must provide an input columnPath e.g. '/*/*/[F,SIGFP]'",
        )
    logger.warning("Input column path provided: %s", input_column_path)
    
    if preferred_dest is None:
        raise Exception(
            "smartSplitMTZ Exception:",
            "Provide first guess for destination of split file",
        )
    logger.warning("Preferred destination provided: %s", preferred_dest)
    
    mtzin = gemmi.read_mtz_file(str(input_file_path))
    mtzin.ensure_asu()

    provided_column_names = mtzin.column_labels()
    if input_column_path.startswith("/"):
        input_column_path = input_column_path[1:]
    if len(input_column_path.split("/")) not in [1, 3]:
        raise Exception("smartSplitMTZ Exception:", "Invalid input columnPath")

    selected_columns = re.sub(r"[\[\] ]", "", input_column_path.split("/")[-1]).split(
        ","
    )
    output_columns, type_signature = get_output_columns(
        selected_columns, provided_column_names, mtzin, input_column_path
    )

    logger.info("Output columns are %s", output_columns)
    final_dest = available_file_name_based_on(preferred_dest)

    mtzout = gemmi.Mtz()
    mtzout.spacegroup = mtzin.spacegroup
    mtzout.cell = mtzin.cell
    # Create a base dataset object with the same spacegroup and cell as the input file and populate with a full set of unique reflections
    hkl_base = mtzout.add_dataset("HKL_base")
    mtzout.add_column("H", "H")
    mtzout.add_column("K", "H")
    mtzout.add_column("L", "H")
    uniques = gemmi.make_miller_array(
        mtzout.cell, mtzout.spacegroup, mtzin.resolution_high(), mtzin.resolution_low()
    )
    mtzout.set_data(uniques)

    dataset = hkl_base
    if len(mtzin.datasets) > 1:
        dataset = output_columns[-1].dataset
        ds = mtzout.add_dataset(dataset.project_name)
        ds.crystal_name = dataset.crystal_name
        ds.dataset_name = dataset.dataset_name
        ds.wavelength = dataset.wavelength

    output_column_labels = []
    logger.warning("Type signature: %s", type_signature)
    labels_dict = {
        "FQ": {"cls": CCP4XtalData.CObsDataFile, "contentType": 4},
        "JQ": {"cls": CCP4XtalData.CObsDataFile, "contentType": 3},
        "GLGL": {"cls": CCP4XtalData.CObsDataFile, "contentType": 2},
        # surely not
        "FQFQ": {"cls": CCP4XtalData.CObsDataFile, "contentType": 2},
        "KMKM": {"cls": CCP4XtalData.CObsDataFile, "contentType": 1},
        # surely not
        "JQJQ": {"cls": CCP4XtalData.CObsDataFile, "contentType": 1},
        "AAAA": {"cls": CCP4XtalData.CPhsDataFile, "contentType": 1},
        "PW": {"cls": CCP4XtalData.CPhsDataFile, "contentType": 2},
        "I": {"cls": CCP4XtalData.CFreeRDataFile, "contentType": 1},
        # map coefficient files
        "FP": {"cls": CCP4XtalData.CMapCoeffsDataFile, "contentType": 1},
    }
    output_column_labels.extend(
        getattr(labels_dict[type_signature]["cls"], "CONTENT_SIGNATURE_LIST")[
            labels_dict[type_signature]["contentType"] - 1
        ]
    )
    logger.warning("Output column labels: %s", str(output_column_labels))
    for i, column in enumerate(output_columns):
        new_column = mtzout.copy_column(-1, column)
        new_column.label = output_column_labels[i]
        new_column.dataset_id = dataset.id

    mtzout.history = [f"MTZ file created from {input_file_path.name} using gemmi."]
    mtzout.write_to_file(str(final_dest))

    if return_metadata:
        # Extract metadata for annotation
        metadata = {
            "crystal_name": dataset.crystal_name if dataset != hkl_base else None,
            "dataset_name": dataset.dataset_name if dataset != hkl_base else None,
            "column_labels": output_column_labels,
            "original_columns": [col.label for col in output_columns],
        }
        return {"path": final_dest, "metadata": metadata}

    return final_dest


def get_output_columns(
    selected_columns: List[str],
    provided_column_names: List[str],
    mtzin: gemmi.Mtz,
    input_column_path: str,
):
    output_columns = []
    type_signature = ""
    for columnLabel in selected_columns:
        if provided_column_names.count(columnLabel) == 1:
            column = mtzin.column_with_label(columnLabel)
            output_columns.append(column)
            type_signature += column.type
        else:
            if len(input_column_path.split("/")) != 3:
                raise Exception(
                    "smartSplitMTZ Exception:",
                    "Input file requires full input columnPath e.g. '/crystal/dataset/[F,SIGFP]'",
                )
            for dataset in mtzin.datasets:
                if (
                    dataset.crystal_name == input_column_path.split("/")[-3]
                    and dataset.dataset_name == input_column_path.split("/")[-2]
                ):
                    column = mtzin.column_with_label(
                        columnLabel, mtzin.dataset(dataset.id)
                    )
                    output_columns.append(column)
                    type_signature += column.type

    if len(output_columns) != len(selected_columns):
        raise Exception(
            "smartSplitMTZ Exception:",
            f"Unable to select columns from input file - options are {provided_column_names}",
        )

    return output_columns, type_signature
