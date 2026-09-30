==================================
Space Group Validation with Zanuda
==================================

`Zanuda <http://scripts.iucr.org/cgi-bin/paper?S1399004714014795>`_ was developed to automate the validation of space group in case of the presence of pseudosymmetry and twinning. The program can be used to restore the correct space group in structures which were intentionally solved in low symmetry space groups including P1 (the highest symmetry = the smallest cell volume). The validation is based on the results of a series of refinements in space groups, which are compatible with the observed unit cell parameters.

If the pseudosymmetry is very close to an exact crystallographic symmetry, the structure can be solved and partially refined in the **wrong space group**. Typically, in such false structures all or some of the *pseudosymmetry operations are treated as crystallographic symmetry operations and vice versa*. Such misassignment is not uncommon when the structure is solved by molecular replacement (MR) and ``it becomes apparent, when the R-free ceases to decrease at about 35% or even at a higher value, and no further model rebuilding and refinement can improve it.`` At this point the electron density map remains imperfect (breaks in the main chain electron density, poor solvent peaks) while does not suggest any particular ways of model improvement.

In case of space group validation with Zanuda was made the next assumptions:

    * The pseudosymmetry operations, if any, are close enough to exact symmetry operations and, therefore, refinement converges to the global minimum when wrong symmetry constraints are removed and correct constraints are imposed.

    * The errors in individual macromolecules do not hinder the difference between pseudosymmetry and crystallographic symmetry, i.e. the model is already refined well enough (R-free around or below 40%). ``However it is not assumed that this refinement has been performed in the true space group``.

The search for a minimum free energy and, as a consequence, the regular packing of molecules in a crystal lattice often leads to a symmetric relationship between the molecules. A characteristic of a crystal is that it has unit translations in three dimensions, also called three-dimensional translational symmetry, corresponding to the repetition of the unit cells. Application of the symmetry operators, such as rotation, reflection, or inversion, leaves the entire crystal unchanged if the space group was chosen right. There are 230 different ways to combine the allowed symmetry operations in a crystal, leading to 230 space groups. They can be found in the `International Tables for Crystallography <https://it.iucr.org/Ac/>`_
Only 65 space groups are “biological”. The reason is that in protein crystals, the application of mirror planes and inversion centers (centers of symmetry) would change the asymmetry of the amino acids: An l-amino acid would become a d-amino acid, but these are never found in proteins.

**How do we deduce the Space Group in practice?**

• We start in reciprocal space (point group)

• We go all way back from symmetry in reciprocal space to crystal space group

    Data processing gives values of the unit cell parameters

    Lattice symmetry is derived from the unit cell parameters

    Comparison of related intensities gives crystal point group

    Systematic absences allow to reduce the number of possible space groups

    Space group is only a hypothesis until structure is complete

-----------
Point Group
-----------

In three-dimensional objects, such as crystals, symmetry elements may be present in several combinations.  In fact, in crystals there are 32 possible combinations of symmetry elements.  These operations can be thus combined to form a group of symmetry operations. These groups of operations are called **point groups** because the symmetry elements of these operations all pass through a single point of the object. Examples for point group symmetry operations are **rotation axes**, **inversion axes**, and **mirror planes**. ``In case of point group you have a point to which all symmetry operations go through`` with a net result that the actual molecule does not change.

More about **Symmetry in Crystallography** you can find `here <http://xrayweb.chem.ou.edu/notes/symmetry.html>`__

--------------
Pseudosymmetry
--------------

Crystals have an ordered internal arrangement of atoms i.e. the atoms are arranged in a symmetrical fashion on a three-dimensional network referred to as a lattice.
Symmetry operator defines a rotation and a translation of the crystal such that each atom in the repositioned copy matches a certain atom in the original. **Pseudo-symmetry operations** are defined similarly, except that the coordinates of matching atoms are not required to coincide exactly. Therefore, it is convenient to define a **pseudo-symmetry space group (PSSG)** which contains both all the operations from the crystal space group and all the pseudo-symmetry operations.

.. note::
    Noncrystallographic symmetry (NCS) and pseudo-symmetry are different concepts. An NCS operation is local and is defined by the best overlap of two NCS-related molecules after applying the NCS operation to one of them. In contrast, the pseudo-symmetry operation is global and is defined by the best match between the entire crystal and its transformed copy. Thus, the NCS operation and the pseudo-symmetry operation relating the same two molecules are in general different operations and may coincide only in special cases.

**Crystallographic symmetry**: symmetry is **global** and **exact**

**Generic Non-Crystallographic Symmetry (NCS)**: symmetry is	**local** and **approximate**

**Pseudosymmetry**: symmetry is **global**	and	**approximate**

In structures with one molecule per asymmetric unit (AU) there is no pseudosymmetry and PSSG coincides with the space group of the crystal. In many cases of NCS the global mapping of the crystal on itself cannot be defined even formally and PSSG remains equal to the crystal space group.

**Zanuda imposes the upper limit of 3 Å for the C-α r.m.s.d. between the structure and its copy generated by an additional global operation**. Global operations with larger values of r.m.s.d. are ignored as they are unlikely to be misinterpreted.




**References**

`Lebedev A.  & Isupov M.  (2014). Acta Cryst. D70, 2430–2443. <http://scripts.iucr.org/cgi-bin/paper?S1399004714014795>`_

`Drenth J. (2007) Principles of Protein X-Ray Crystallography. Springer-Verlag New York  <https://www.springer.com/gp/book/9780387333342>`_


**ACKNOWLEDGEMENTS**

This article uses materials provided by Andrey Lebedev.

More about **space group validation with Zanuda** you can find `here <https://www.ccp4.ac.uk/newsletters/newsletter48/articles/Zanuda/zanuda.html>`__
