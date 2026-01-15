import gemmi


def servalcat_utils_read_cif_safe(cif_in):
    # with gzip.open(cif_in, "rt") if cif_in.endswith(".gz") else open(cif_in) as ifs:
    with open(cif_in) as ifs:
        s = ifs.read()
    if "\0" in s: # Refmac occasionally writes \0 in some fields..
        print(" WARNING: null character detected. Replacing with '.'")
        s = s.replace("\0", ".")
    doc = gemmi.cif.read_string(s)
    return doc


def servalcat_utils_write_mmcif(st, cif_out, cif_ref=None, cif_ref_doc=None):
    """
    Refmac fails if _entry.id is longer than 80 chars including quotations
    """
    st_new = st.clone()
    print("Writing mmCIF file: {}".format(cif_out))
    if cif_ref or cif_ref_doc:
        if cif_ref:
            print("  using mmCIF metadata from: {}".format(cif_ref))
        groups = gemmi.MmcifOutputGroups(False)
        groups.group_pdb = True
        groups.ncs = True
        groups.atoms = True
        groups.cell = True
        groups.scale = True
        groups.assembly = True
        groups.entity = True
        groups.entity_poly = True
        groups.entity_poly_seq = True
        groups.cis = True
        groups.conn = True
        groups.software = True
        groups.auth_all = True
        # FIXME is this all?
        if cif_ref:
            try:
                cif_ref_doc = servalcat_utils_read_cif_safe(cif_ref)
            except Exception as e:
                # Sometimes refmac writes a broken mmcif file..
                print("Error in mmCIF reading: {}".format(e))
                print("  Give up using cif reference.")
                return servalcat_utils_write_mmcif(st, cif_out)
            
        blocks = list(filter(lambda b: b.find_loop("_atom_site.id"), cif_ref_doc))
        if len(blocks) == 0:
            print("No _atom_site found in reference")
            print("  Give up using cif reference.")
            return servalcat_utils_write_mmcif(st, cif_out)
        block = blocks[0]
        # to remove fract_transf_matrix. maybe we should keep some (like _atom_sites.solution_hydrogens)?
        # we do not want this because cell may be updated
        block.find_mmcif_category("_atom_sites.").erase()
        st_new.update_mmcif_block(block, groups)
        if "_entry.id" in st_new.info: st_new.info["_entry.id"] = st_new.info["_entry.id"][:78]
        cif_ref_doc.write_file(cif_out, options=gemmi.cif.Style.Aligned)
    else:
        st_new.name = st_new.name[:78] # this will become _entry.id
        if "_entry.id" in st_new.info: st_new.info["_entry.id"] = st_new.info["_entry.id"][:78]
        groups = gemmi.MmcifOutputGroups(True, auth_all=True)
        doc = gemmi.cif.Document()
        block = doc.add_new_block("new")
        st_new.update_mmcif_block(block, groups)
        doc.write_file(cif_out, options=gemmi.cif.Style.Aligned)


def fix_tls_cif_hetero_aniso(cif_in, cif_out):

    u_to_b = 8 * 3.141592653589793**2
    b_to_u = 1. / u_to_b

    doc = gemmi.cif.read(cif_in)
    st = gemmi.read_structure(cif_in)
    tls_groups = {int(x.id): x for x in st.meta.refinement[0].tls_groups}
    tls_details = doc[0].find_value("_ccp4_refine_tls.details")
    if tls_groups and gemmi.cif.as_string(tls_details) == "U values: residual only":
        for cra in st[0].all():
            tlsgr = tls_groups.get(cra.atom.tls_group_id)
            if cra.atom.tls_group_id > 0 and tlsgr is not None:
                if not cra.atom.aniso.nonzero():
                    u = cra.atom.b_iso * b_to_u
                    cra.atom.aniso = gemmi.SMat33f(u, u, u, 0, 0, 0)
                u_from_tls = gemmi.calculate_u_from_tls(tlsgr, cra.atom.pos)
                cra.atom.aniso += gemmi.SMat33f(*u_from_tls.elements_pdb())
                cra.atom.b_iso = cra.atom.aniso.trace() / 3. * u_to_b
        doc[0].set_pair("_ccp4_refine_tls.details", gemmi.cif.quote("U values: with tls added"))
        servalcat_utils_write_mmcif(st, cif_out, cif_ref_doc=doc)


if __name__ == "__main__":
    import sys
    cif_in, cif_out = sys.argv[1:]
    fix_tls_cif_hetero_aniso(cif_in, cif_out)
    