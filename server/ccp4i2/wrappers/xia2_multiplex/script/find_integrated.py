"""Find integrated.{expt,refl} file pairs for DIALS integration output,
recursively under a root directory"""

import glob
import os
import os.path
from pathlib import Path
from dxtbx.model.experiment_list import ExperimentList


def is_xia2_dir(path):
    p = Path(path)
    if os.path.isdir(p / "DataFiles") and os.path.exists(p / "xia2.txt"):
        return True
    else:
        return False


def find_integrated(root_dir):

    if not root_dir:
        return []

    # Find directories to search in
    search = [dirpath for dirpath, _, _ in os.walk(root_dir)]

    # xia2 keeps multiple copies of integrated output. Filter out anything
    # that looks like a xia2 directory
    xia2_dirs = [d for d in search if is_xia2_dir(d)]
    for xia2_dir in xia2_dirs:
        search = [s for s in search if not xia2_dir in s]

    # Now add the xia2 DataFiles directories only back in
    for xia2_dir in xia2_dirs:
        search.append(Path(xia2_dir) / "DataFiles")

    # Iterate over all .expt files
    integrated = []
    identifiers = set()
    for pth in search:
        for expt in glob.iglob(str(Path(pth) / "*.expt")):

            # Only keep files for which a prefix.{expt,refl} pair exists
            prefix = os.path.splitext(expt)[0]
            refl = prefix + ".refl"
            if not os.path.exists(refl):
                continue

            # Load the experiments
            el = ExperimentList.from_file(expt, check_format=False)

            # Reject file if any experiments don't have a profile model
            if not all((e.profile for e in el)):
                continue

            # Reject file if any experiments have a scaling model
            if any((e.scaling_model for e in el)):
                continue

            # Reject file if it contains experimental identifiers we've already seen
            ids = [e.identifier for e in el]
            if any((i in identifiers for i in ids)):
                continue
            for i in ids:
                identifiers.add(i)

            # Otherwise, looks good. Keep it
            integrated.append(prefix)

    return sorted(integrated)
