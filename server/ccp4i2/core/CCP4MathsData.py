"""
CCP4MathsData.py --- CData classes.

These were once two classes each: a generated stub carrying the data
model, and an implementation carrying the methods. The generator is no
longer run and the split cost more than it saved --- the two halves
interleaved in the MRO, so an implementation could drop out of its own
subclass's ancestry and `isinstance` would say no to a file that plainly
was one. They are one class now.
"""

from typing import Optional, Any
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType, content
from ccp4i2.core.base_object.base_classes import CData
from ccp4i2.core.base_object.fundamental_types import CFloat



class CMatrix33(CData):


    """
    Extends CMatrix33 with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMatrix33.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CXyzBox(CData):


    """
    Extends CXyzBox with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        error_codes = {
            "201": {
                "description": "Maximum x,y or z value less than minimum"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    xMin = content(
        "CFloat",
        guiLabel='x min',
        toolTip='Lower x bound of the box, in Angstroms')
    yMin = content(
        "CFloat",
        guiLabel='y min',
        toolTip='Lower y bound of the box, in Angstroms')
    zMin = content(
        "CFloat",
        guiLabel='z min',
        toolTip='Lower z bound of the box, in Angstroms')
    xMax = content(
        "CFloat",
        guiLabel='x max',
        toolTip='Upper x bound of the box, in Angstroms')
    yMax = content(
        "CFloat",
        guiLabel='y max',
        toolTip='Upper y bound of the box, in Angstroms')
    zMax = content(
        "CFloat",
        guiLabel='z max',
        toolTip='Upper z bound of the box, in Angstroms')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXyzBox.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CAngle(CFloat):


    """
    An angle
    
    Extends CAngle with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "max": None,
            "min": None,
            "enumerators": [],
            "menuText": [],
            "onlyEnumerators": False,
        }
    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAngle.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CXyz(CData):


    """
    Extends CXyz with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        error_codes = {
            "201": {
                "description": "Attempting arithmetic with inappropriate data type"
            },
            "202": {
                "description": "Attempting arithmetic in unset data object"
            },
            "203": {
                "description": "Attempting arithmetic with unset data object as argument"
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    x = content("CFloat", guiLabel='x', toolTip='x coordinate, in Angstroms')
    y = content("CFloat", guiLabel='y', toolTip='y coordinate, in Angstroms')
    z = content("CFloat", guiLabel='z', toolTip='z coordinate, in Angstroms')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXyz.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CEulerRotation(CData):


    """
    Extends CEulerRotation with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    alpha = content("CAngle", guiLabel='Alpha', toolTip='First Euler angle, in degrees')
    beta = content("CAngle", guiLabel='Beta', toolTip='Second Euler angle, in degrees')
    gamma = content("CAngle", guiLabel='Gamma', toolTip='Third Euler angle, in degrees')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEulerRotation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CTransformation(CData):


    """
    Extends CTransformation with implementation-specific methods.

    Provides flat access to rotation angles (alpha, beta, gamma) and
    translation components (x, y, z) for backward compatibility with
    wrappers that expect direct attribute access (e.g. gesamt).
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    translation = content(
        "CXyz",
        guiLabel='Translation',
        toolTip='Translation part of the transformation, in Angstroms')
    rotation = content(
        "CEulerRotation",
        guiLabel='Rotation',
        toolTip='Rotation part of the transformation, as Euler angles')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTransformation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    @property
    def alpha(self):
        return self.rotation.alpha

    @property
    def beta(self):
        return self.rotation.beta

    @property
    def gamma(self):
        return self.rotation.gamma

    @property
    def x(self):
        return self.translation.x

    @property
    def y(self):
        return self.translation.y

    @property
    def z(self):
        return self.translation.z

