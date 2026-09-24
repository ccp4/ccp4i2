# pandda_baz2b_mini

Three BAZ2B fragment-screen datasets, exactly as a CCP4i2 dimple + acedrg
pair leaves them (`<label>.dimple.pdb`, `<label>.dimple.mtz`, `ligand.cif`),
so that PanDDA staging can be tested by reproducing the contract's input tree
from real sources without an external volume being mounted.

They are the first three rows of the 201-dataset staged tree on
`/Volumes/LocalStore/pandda/BAZ2B` (`xtal-0000`..`xtal-0002`), which is what
`tests/unit/pandda/test_staging.py` diffs against when that volume is present.

| Label | Cell (C 2 2 21) | d_min | FreeR column |
|---|---|---|---|
| BAZ2BA-x425 | 82.5 97.0 58.1 | 1.72 | `FreeR_flag` |
| BAZ2BA-x427 | 82.6 96.6 58.0 | 1.69 | `FreeR_flag` |
| BAZ2BA-x428 | 82.4 96.9 58.0 | 1.78 | `FreeR_flag` |

All three dictionaries carry `_chem_comp_bond.type` (CCP4 monomer-library
spelling). The `value_order` case that `prepare_dict_for_pandda()` exists for
is manufactured by the tests, deliberately, because no shipped fixture has it.

Source: Diamond Light Source XChem BAZ2B screen (Bradley et al.), via the
PanDDA tutorial data. Small enough to live here: ~1.8 MB per dataset.
