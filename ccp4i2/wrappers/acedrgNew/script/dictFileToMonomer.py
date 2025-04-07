import collections
import sys

from ccp4mg import mmdb2
import ccp4srs


def dictFileToMonomer(dictFileName):
    f = mmdb2.File()
    rc = f.ReadMMCIFFile(dictFileName)
    if rc == 0:
      # FIXME - *REALLY* need to deal with multiple compids in the file!
      chemID = ""
      for i in range(f.GetNofData()):
        d = f.GetCIFData(i)
        if d.GetDataName()=="comp_list" and d.GetLoopLength("_chem_comp")>0:
          l = d.GetLoop("_chem_comp")
          chemID = l.GetString("id",0)[0]

      if not chemID:
          try:
              chemID = d.GetDataName()
          except:
              pass

      #print chemID, d.GetDataName(), f.GetNofData()
      if chemID:
        have_atoms = False
        have_bonds = False
        have_angles = False
        have_torsions = False
        for i in range(f.GetNofData()):
          d = f.GetCIFData(i)
          if (d.GetDataName()=="comp_"+chemID and d.GetLoopLength("_chem_comp_atom")>0) or (d.GetDataName()==chemID and d.GetLoopLength("_chem_comp_atom")>0):
            have_atoms = True
          if (d.GetDataName()=="comp_"+chemID and d.GetLoopLength("_chem_comp_bond")>0) or (d.GetDataName()==chemID and d.GetLoopLength("_chem_comp_bond")>0):
            have_bonds = True
          if (d.GetDataName()=="comp_"+chemID and d.GetLoopLength("_chem_comp_angle")>0) or (d.GetDataName()==chemID and d.GetLoopLength("_chem_comp_angle")>0):
            have_angles = True
          if (d.GetDataName()=="comp_"+chemID and d.GetLoopLength("_chem_comp_tor")>0) or (d.GetDataName()==chemID and d.GetLoopLength("_chem_comp_tor")>0):
            have_torsions = True

      #FIXME - Urgh! d is in loop above, but not here!
      if have_bonds and have_atoms:
        monomer = ccp4srs.Monomer()
        # FIXME - we now have enough info to construct a Graph from the dict
        l = d.GetLoop("_chem_comp_atom")
        atsym = {}
        atomIdIndex = []
        for i in range(l.GetLoopLength()):
          at = l.GetString("atom_id",i)[0]
          atomIdIndex.append(at)
          sym = l.GetString("type_symbol",i)[0]
          try:
             charge = float(l.GetString("partial_charge",i)[0])
          except:
             try:
                 charge = float(l.GetString("charge",i)[0])
             except:
                 charge = 0.0
          type_symbol = l.GetString("type_symbol",i)[0]
          atsym[at] = sym
          pat = ccp4srs.SRSAtom()
          try:
              x = float(l.GetString("x",i)[0])
              y = float(l.GetString("y",i)[0])
              z = float(l.GetString("z",i)[0])
              pat.set_x_ccp4_mlib(x)
              pat.set_y_ccp4_mlib(y)
              pat.set_z_ccp4_mlib(z)
              pat.set_x_rcsb_ideal(x)
              pat.set_y_rcsb_ideal(y)
              pat.set_z_rcsb_ideal(z)
              pat.set_x_rcsb_cartn(x)
              pat.set_y_rcsb_cartn(y)
              pat.set_z_rcsb_cartn(z)
          except:
              pass #Do not set coordinates
          try:
              type_energy = l.GetString("type_energy",i)[0]
              pat.set_energy_type(type_energy)
          except:
              pass #Do not set energy type
          pat.set_atom_charge(charge)
          pat.set_element(type_symbol)
          pat.set_name(at)
          pat.thisown = 0
          monomer.add(pat)

        l = d.GetLoop("_chem_comp_bond")
        for i in range(l.GetLoopLength()):
          at1 = l.GetString("atom_id_1",i)[0]
          at2 = l.GetString("atom_id_2",i)[0]
          bond = ccp4srs.Bond()
          try:
              value_dist = float(l.GetString("value_dist",i)[0])
              bond.set_bond_length(value_dist)
          except:
             pass #Do not set distance
          try:
              value_dist_esd = float(l.GetString("value_dist_esd",i)[0])
              bond.set_bond_length_esd(value_dist_esd)
          except:
             pass #Do not set distance sigma
          try:
              bondType = l.GetString("type",i)[0]
          except:
              try:
                  bondType = l.GetString("value_order",i)[0]
              except:
                  pass #Possibly raise, this is bad
          aromatic = False
          try:
             aromatic_cif = l.GetString("aromatic",i)[0]
             if aromatic_cif.lower() == "y":
                 aromatic = True
          except:
              try:
                 aromatic_cif = l.GetString("pdbx_aromatic_flag",i)[0]
                 if aromatic_cif.lower() == "y":
                     aromatic = True
              except:
                 pass
          bond.set_atom_1(atomIdIndex.index(at1))
          bond.set_atom_2(atomIdIndex.index(at2))
          if aromatic:
              bond.set_bond_order(ccp4srs.Bond.Aromatic)
          elif bondType.lower()[:4] == "doub":
              bond.set_bond_order(ccp4srs.Bond.Double)
          elif bondType.lower()[:4] == "trip":
              bond.set_bond_order(ccp4srs.Bond.Triple)
          elif bondType.lower()[:4] == "arom":
              bond.set_bond_order(ccp4srs.Bond.Aromatic)
          elif bondType.lower()[:4] == "delo":
              bond.set_bond_order(ccp4srs.Bond.Deloc)
          elif bondType.lower()[:4] == "cova":
              bond.set_bond_order(ccp4srs.Bond.Covalent)
          elif bondType.lower()[:4] == "meta":
              bond.set_bond_order(ccp4srs.Bond.Metal)
          elif bondType.lower()[:4] == "sing":
              bond.set_bond_order(ccp4srs.Bond.Single)
          else:
              bond.bond_order(ccp4srs.Bond.noOrder)
          bond.thisown = 0
          monomer.add(bond)

        l = d.GetLoop("_chem_comp_angle")
        if have_angles:
            for i in range(l.GetLoopLength()):
              at1 = l.GetString("atom_id_1",i)[0]
              at2 = l.GetString("atom_id_2",i)[0]
              at3 = l.GetString("atom_id_3",i)[0]
              value_angle = float(l.GetString("value_angle",i)[0])
              value_angle_esd = float(l.GetString("value_angle_esd",i)[0])
              angle = ccp4srs.Angle()
              angle.set_atom_1(atomIdIndex.index(at1))
              angle.set_atom_2(atomIdIndex.index(at2))
              angle.set_atom_3(atomIdIndex.index(at3))
              angle.set_angle(value_angle)
              angle.set_angle_esd(value_angle_esd)
              angle.thisown = 0
              monomer.add(angle)

        l = d.GetLoop("_chem_comp_tor")
        if have_torsions:
            for i in range(l.GetLoopLength()):
              at1 = l.GetString("atom_id_1",i)[0]
              at2 = l.GetString("atom_id_2",i)[0]
              at3 = l.GetString("atom_id_3",i)[0]
              at4 = l.GetString("atom_id_4",i)[0]
              value_angle = float(l.GetString("value_angle",i)[0])
              value_angle_esd = float(l.GetString("value_angle_esd",i)[0])
              period = int(l.GetString("period",i)[0])
              torsion_id = l.GetString("id",i)[0]
              torsion = ccp4srs.Torsion()
              torsion.set_atom_1(atomIdIndex.index(at1))
              torsion.set_atom_2(atomIdIndex.index(at2))
              torsion.set_atom_3(atomIdIndex.index(at3))
              torsion.set_atom_4(atomIdIndex.index(at3))
              torsion.set_angle(value_angle)
              torsion.set_angle_esd(value_angle_esd)
              torsion.set_period(period)
              torsion.set_torsion_id(torsion_id)
              torsion.thisown = 0
              monomer.add(torsion)

        l = d.GetLoop("_chem_comp_chir")
        if l is not None:
            for i in range(l.GetLoopLength()):
              at0 = l.GetString("atom_id_centre",i)[0]
              at1 = l.GetString("atom_id_1",i)[0]
              at2 = l.GetString("atom_id_2",i)[0]
              at3 = l.GetString("atom_id_3",i)[0]
              chiid = l.GetString("id",i)[0]
              volume_sign = l.GetString("volume_sign",i)[0]
              chi_centre = ccp4srs.ChiCenter()
              chi_centre.set_atom_center(atomIdIndex.index(at0))
              chi_centre.set_atom_1(atomIdIndex.index(at1))
              chi_centre.set_atom_2(atomIdIndex.index(at2))
              chi_centre.set_atom_3(atomIdIndex.index(at3))
              chi_centre.set_chicenter_id(chiid)
              #Hmm, ccp4srs::ChiCenter sign enum is private.
              if volume_sign == "positiv":
                  #chi_centre.set_volume_sign(ccp4srs.ChiCenter.Positive
                  chi_centre.set_volume_sign(1)
              elif volume_sign == "negativ":
                  #chi_centre.set_volume_sign(ccp4srs.ChiCenter.Negative
                  chi_centre.set_volume_sign(2)
              elif volume_sign == "both":
                  #chi_centre.set_volume_sign(ccp4srs.ChiCenter.Both
                  chi_centre.set_volume_sign(3)
              else:
                  #chi_centre.set_volume_sign(ccp4srs.ChiCenter.noSign
                  chi_centre.set_volume_sign(0)
              chi_centre.thisown = 0
              monomer.add(chi_centre)

        planes = collections.OrderedDict()
        l = d.GetLoop("_chem_comp_plane_atom")
        if l is not None:
            for i in range(l.GetLoopLength()):
              at = l.GetString("atom_id",i)[0]
              plane_id = l.GetString("plane_id",i)[0]
              dist_esd = float(l.GetString("dist_esd",i)[0])
              if not plane_id in planes:
                  planes[plane_id] = []
              planes[plane_id].append((atomIdIndex.index(at),dist_esd))
        for k,v in list(planes.items()):
            plane = ccp4srs.Plane()
            ids = []
            esds = []
            for atid,esd in v:
                ids.append(atid)
                esds.append(esd)
            plane.set_plane_id(k)
            plane.set_n_atoms(len(ids))
            plane.set_plane_atoms(ids)
            plane.set_dist_esds(esds)
            plane.thisown = 0
            monomer.add(plane)

        return monomer
      else:
        return None

if __name__ == "__main__":
    monomer = dictFileToMonomer(sys.argv[1])
    print("atoms",monomer.n_atoms())
    print("bonds",monomer.n_bonds())
    print("angles",monomer.n_angles())
    print("torsions",monomer.n_torsions())
    print("chirals",monomer.n_chicenters())
    print("planes",monomer.n_planes())

    rcSRS = mmdb2.intp()
    gSRS = monomer.getGraph(rcSRS)
    gSRS.Print()
