"""The model Crank2 refined, as one file: the built model plus the substructure.

Crank2 (and its SHELX route) hands on XYZOUT, the built model, and
XYZOUT_SUBSTR, the anomalous scatterers, separately; its own final REFMAC
refinement modelled both. Refining XYZOUT alone leaves the scatterers out of
the model, which matters as much as they scatter: for a mercury soak of a
10 kDa protein (HypF, two Hg sites) Servalcat went from R1-free 0.29 with the
sites to 0.44 without them, the protein distorting to explain the mercury.

A site is not always a new atom. Where the scatterer is part of the model it
is already there: S-SAD sites sit on cysteine SG and methionine SD, and a
selenomethionine's Se sits where the model, built as methionine, has an SD.
So each site is
  - left out if the model has an atom of the same element within SAME_ATOM
    of it (any symmetry copy);
  - for Se, made part of the model if it lies on a methionine SD: that
    residue becomes MSE, with its SD replaced by SE;
  - left out, and said so, if any other model atom is within CLASH of it
    (no real scatterer sits that close to a protein atom; it is a misplaced
    site or a misbuilt residue, and adding it would only distort refinement);
  - otherwise added, as a HETATM residue of its own in a new chain, with the
    occupancy and B factor Crank2 refined.
"""
import gemmi

SAME_ATOM = 2.0
CLASH = 1.5

# Monomer-library names for single-atom residues whose code is not the element
ION_RESIDUE = {"I": "IOD"}


def _near(ns, model, cell, pos, radius):
    """(distance, cra) for every model atom within radius of pos, any image."""
    found = []
    for mark in ns.find_atoms(pos, "\0", radius=radius):
        cra = mark.to_cra(model)
        d = cell.find_nearest_image(pos, cra.atom.pos, gemmi.Asu.Any).dist()
        if d <= radius:
            found.append((d, cra))
    return sorted(found, key=lambda x: x[0])


def _free_chain_name(model):
    used = {ch.name for ch in model}
    for name in "WXYZHIJKLMNOPQRSTUVABCDEFG" + "abcdefghijklmnopqrstuvwxyz0123456789":
        if name not in used:
            return name
    raise ValueError("no free single-character chain name")


def complete_model(model_path, substr_path, out_path):
    """Write the model plus its substructure to out_path; return what was done.

    The result is a dict: "added" (sites added as atoms), "on_model" (sites
    the model already has: "<element> <residue>" for each), "converted"
    (methionines made MSE), "clashes" (sites left out for a clash).
    """
    st = gemmi.read_structure(str(model_path))
    sub = gemmi.read_structure(str(substr_path))
    if len(st) == 0 or len(st[0]) == 0:
        raise ValueError(f"{model_path}: no model to complete")
    st.setup_entities()
    model = st[0]
    cell = st.cell if st.cell.is_crystal() else sub.cell
    ns = gemmi.NeighborSearch(model, cell, 5).populate()

    report = {"added": [], "on_model": [], "converted": [], "clashes": []}
    sites = [atom for chain in sub[0] if len(sub) for residue in chain for atom in residue] if len(sub) else []
    additions = []
    to_mse = []
    for site in sites:
        element = site.element.name
        near = _near(ns, model, cell, site.pos, SAME_ATOM)
        same = [(d, cra) for d, cra in near if cra.atom.element.name == element]
        if same:
            d, cra = same[0]
            report["on_model"].append(f"{element} {cra.chain.name}/{cra.residue.name}{cra.residue.seqid}/{cra.atom.name} {d:.2f} A")
            continue
        met_sd = [(d, cra) for d, cra in near if element == "Se" and cra.residue.name == "MET"
                  and cra.atom.name == "SD"]
        if met_sd:
            d, cra = met_sd[0]
            to_mse.append((cra.chain.name, str(cra.residue.seqid)))
            continue
        clash = [(d, cra) for d, cra in near if d < CLASH]
        if clash:
            d, cra = clash[0]
            report["clashes"].append(f"{element} {d:.2f} A from {cra.chain.name}/{cra.residue.name}{cra.residue.seqid}/{cra.atom.name}")
            continue
        additions.append(site)

    for chain_name, seqid in dict.fromkeys(to_mse):
        for residue in model[chain_name]:
            if str(residue.seqid) == seqid and residue.name == "MET":
                residue.name = "MSE"
                residue.het_flag = "H"
                for atom in residue:
                    if atom.name == "SD":
                        atom.name = "SE"
                        atom.element = gemmi.Element("Se")
                report["converted"].append(f"{chain_name}/MET{seqid}")

    if additions:
        chain = gemmi.Chain(_free_chain_name(model))
        for number, site in enumerate(additions, start=1):
            residue = gemmi.Residue()
            residue.name = ION_RESIDUE.get(site.element.name.upper(), site.element.name.upper())
            residue.seqid = gemmi.SeqId(number, " ")
            residue.het_flag = "H"
            atom = site.clone()
            atom.name = site.element.name.upper()
            residue.add_atom(atom)
            chain.add_residue(residue)
            report["added"].append(f"{residue.name} occupancy {site.occ:.2f}")
        model.add_chain(chain)

    st.setup_entities()
    if str(out_path).lower().endswith((".cif", ".mmcif")):
        st.make_mmcif_document().write_file(str(out_path))
    else:
        st.write_pdb(str(out_path))
    return report
