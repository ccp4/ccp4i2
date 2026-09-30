"""The edits MakeLink can make to a monomer, as data rather than commands.

AceDRG takes its modifications as an imperative word stream ("DELETE ATOM O4A
2 CHANGE BOND C4A O4A SINGLE 2 CHARGE 2 N1 1"), where every section has its
own argument order and a malformed one can hang it. The task instead holds
what the modified monomer should be -- which atoms are gone, what order a
bond now has, what charge an atom now carries -- and the plugin writes the
words (link_instruction.py). Declared as data, an edit list cannot say two
contradictory things about one bond or one atom without validity() seeing it.

Resolvable by the def.xml class-name lookup because MakeLink declares this
module in its Task.dataTypes (core/tasks.py).
"""
from ccp4i2.core.base_object.class_metadata import content
from ccp4i2.core.CCP4Data import CData


class CMakeLinkBondOrder(CData):
    """A bond within one monomer, and the order it should have."""

    class Meta:
        contents_order = ['ATOM_1', 'ATOM_2', 'ORDER']
        qualifiers = {"toolTip": "A bond of the monomer and its new order"}

    ATOM_1 = content("CString", guiLabel='Atom', toolTip='Dictionary name of one atom of the bond')
    ATOM_2 = content("CString", guiLabel='Atom', toolTip='Dictionary name of the other atom of the bond')
    ORDER = content(
        "CString", guiLabel='New order', toolTip='The order the bond should have',
        enumerators=['SINGLE', 'DOUBLE', 'TRIPLE'], onlyEnumerators=True)


class CMakeLinkCharge(CData):
    """An atom within one monomer, and the formal charge it should carry."""

    class Meta:
        contents_order = ['ATOM', 'CHARGE']
        qualifiers = {"toolTip": "An atom of the monomer and its new formal charge"}

    ATOM = content("CString", guiLabel='Atom', toolTip='Dictionary name of the atom')
    CHARGE = content("CInt", guiLabel='Charge', toolTip='The formal charge the atom should carry',
                     min=-3, max=3)
