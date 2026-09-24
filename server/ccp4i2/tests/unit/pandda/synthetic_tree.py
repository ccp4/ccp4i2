"""A small PanDDA 2 output tree, shaped as PanDDA writes one, built to order.

Real trees live on an external volume; this builds the same shape in a
temporary directory so the reader and the receipt can be tested anywhere,
including the cases no real tree offers on demand -- a declared event map
that is missing, a build with no pose, a run that never wrote its events
table. Maps are tiny gemmi grids; models are one-atom PDBs.
"""
from pathlib import Path

PDB = ("ATOM      1  CA  GLY A   1      {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
       "END\n")
#: What PanDDA writes for a built pose: residue LIG, whatever the dictionary.
POSE_PDB = ("HETATM    1  C1  LIG 0   1      {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
            "END\n")
DICT_CIF = """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
{code} {code}
data_comp_{code}
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
{code} C1 C
"""


def write_map(path, value=1.0):
    import gemmi
    m = gemmi.Ccp4Map()
    m.grid = gemmi.FloatGrid(4, 4, 4)
    m.grid.fill(value)
    m.grid.set_unit_cell(gemmi.UnitCell(10, 10, 10, 90, 90, 90))
    m.grid.spacegroup = gemmi.find_spacegroup_by_name('P1')
    m.update_ccp4_header()
    m.write_ccp4_map(str(path))


def write_pdb(path, xyz=(1.0, 2.0, 3.0), template=PDB):
    Path(path).write_text(template.format(x=xyz[0], y=xyz[1], z=xyz[2]))


def event_record(idx, bdc=0.8, score=0.3, centroid=(15.0, 40.0, 30.0), build=True,
                 build_score=0.57, rscc=0.30, contour=0.98, build_path=None):
    record = {
        "Score": score, "BDC": bdc, "Size": 5.2,
        "Centroid": list(centroid), "Local Strength": 24.0,
        "Position Array": [list(centroid)], "Point Array": [[1, 2, 3]],
    }
    if build:
        record["Build"] = {
            "Build Path": build_path or "", "Ligand Key": "dict", "Score": build_score,
            "Centroid": list(centroid), "BDC": bdc, "Build Score": build_score,
            "Noise": 150.0, "Signal": 90.0, "Num. Contacts": 5, "Num. Points": 191.0,
            "Optimal Contour": contour, "RSCC": rscc,
        }
    return record


def make_tree(root, datasets, *, events_table=True, staged_apo=True, ligand_code="MZ0"):
    """Build ``root`` as a PanDDA 2 output tree.

    ``datasets`` maps dtag -> list of ``(record, options)`` where ``record``
    is an events.yaml entry (see ``event_record``) and ``options`` may set
    ``event_map=False`` (declared but absent) or ``pose=False`` (a Build
    block with no file). A dtag mapped to ``[]`` is a zero-event dataset
    (``events.yaml`` is ``{}``). ``pandda_model=True`` in a dataset's
    options dict (given as a trailing dict) writes the merged model.
    """
    import yaml

    root = Path(root)
    processed = root / "processed_datasets"
    processed.mkdir(parents=True)
    staging = root.parent / "staged_inputs"
    staging.mkdir(exist_ok=True)
    rows = []
    for dtag, spec in datasets.items():
        extra = {}
        if spec and isinstance(spec[-1], dict) and "pandda_model" in spec[-1]:
            extra = spec[-1]
            spec = spec[:-1]
        ddir = processed / dtag
        (ddir / "autobuild").mkdir(parents=True)
        (ddir / "modelled_structures").mkdir()
        (ddir / "ligand_files").mkdir()
        if ligand_code:
            (ddir / "ligand_files" / "dict.cif").write_text(DICT_CIF.format(code=ligand_code))
        # PanDDA symlinks the apo input into its tree; a receipt must follow it.
        apo_real = staging / f"{dtag}-final.pdb"
        write_pdb(apo_real, (0.0, 0.0, 0.0))
        apo_link = ddir / f"{dtag}-pandda-input.pdb"
        if staged_apo:
            apo_link.symlink_to(apo_real)
        else:
            apo_link.symlink_to(staging / f"{dtag}-gone.pdb")   # dangling
        write_map(ddir / f"{dtag}-z_map.native.ccp4", 3.0)
        write_map(ddir / f"{dtag}-ground-state-average-map.native.ccp4", 1.0)
        if extra.get("pandda_model"):
            write_pdb(ddir / "modelled_structures" / f"{dtag}-pandda-model.pdb", template=POSE_PDB)
        records = {}
        for n, entry in enumerate(spec, start=1):
            record, options = (entry if isinstance(entry, tuple) else (entry, {}))
            records[n] = record
            if options.get("event_map", True):
                token = f"{round(1 - record['BDC'], 2):g}"
                write_map(ddir / f"{dtag}-event_{n}_1-BDC_{token}_map.native.ccp4", 2.0)
            if "Build" in record and options.get("pose", True):
                pose = ddir / "autobuild" / f"7_{n}_dict_0.pdb"
                write_pdb(pose, tuple(record["Centroid"]), template=POSE_PDB)
                record["Build"]["Build Path"] = str(pose)
                write_pdb(ddir / f"{dtag}_event_{n}_best_autobuild.pdb", tuple(record["Centroid"]), template=POSE_PDB)
            rows.append((dtag, n, record))
        (ddir / "events.yaml").write_text(yaml.safe_dump(records) if records else "{}\n")
        (ddir / "processed_dataset.yaml").write_text("Summary:\n  Processing Resolution: 1.8\n")
    if events_table:
        analyses = root / "analyses"
        analyses.mkdir()
        header = (",dtag,event_idx,bdc,cluster_size,site_idx,x,y,z,z_peak,"
                  "hit_in_site_probability,interesting\n")
        lines = [header]
        for i, (dtag, n, record) in enumerate(rows):
            x, y, z = record["Centroid"]
            lines.append(f"{i},{dtag},{n},{record['BDC']},59,{1 + i % 2},{x},{y},{z},"
                         f"0.8,{0.9 - 0.1 * i:.2f},False\n")
        (analyses / "pandda_analyse_events.csv").write_text("".join(lines))
    return root
