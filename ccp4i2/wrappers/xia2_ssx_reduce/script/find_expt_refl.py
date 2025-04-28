"""Find integrated*.{expt,refl} or scaled*.{expt,refl} file pairs without duplicates,
recursively under a root directory"""

from pathlib import Path
import glob
import os
import sys

from dxtbx.model.experiment_list import ExperimentList


def is_xia2_dir(path):
    p = Path(path)
    if os.path.isdir(p / "DataFiles") and (os.path.exists(p / "xia2.ssx_reduce.txt") or os.path.exists(p / "xia2.ssx.txt")):
        return True
    else:
        return False


def find_expt_refl(root_dir, preference="integrated"):

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

    search = sorted(search)
    print("Searched directories:")
    print(str(search))
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

            print("--- It has .refl " + str(expt))
            # Load the experiments
            el = ExperimentList.from_file(expt, check_format=False)

            # Reject file if any experiments don't have a profile model
            #if not all((e.profile for e in el)): #                      MM
            #    print("--- Rejected - no profile model " + str(expt))
            #    continue

            if preference == "integrated":
                # Reject file if any experiments have a scaling model
                if any((e.scaling_model for e in el)):
                    continue
            elif preference == "scaled":
                if any((not e.scaling_model for e in el)):
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


if __name__ == "__main__":
    for pth in find_expt_refl(sys.argv[1]):
        print(pth)
