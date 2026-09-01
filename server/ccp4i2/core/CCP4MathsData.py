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



@cdata_class(
    error_codes={
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
)
class CMatrix33(CData):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CMatrix33.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CMatrix33 with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
        "201": {
            "description": "Maximum x,y or z value less than minimum"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
)
class CXyzBox(CData):

    xMin = content("CFloat")
    yMin = content("CFloat")
    zMin = content("CFloat")
    xMax = content("CFloat")
    yMax = content("CFloat")
    zMax = content("CFloat")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXyzBox.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CXyzBox with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "max": None,
        "min": None,
        "enumerators": [],
        "menuText": [],
        "onlyEnumerators": False,
    },
    qualifiers_order=[
        'min',
        'max',
        'onlyEnumerators',
        'enumerators',
        'menuText'],
)
class CAngle(CFloat):

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAngle.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    An angle
    
    Extends CAngle with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
        "201": {
            "description": "Attempting arithmetic with inappropriate data type"
        },
        "202": {
            "description": "Attempting arithmetic in unset data object"
        },
        "203": {
            "description": "Attempting arithmetic with unset data object as argument"
        }
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
)
class CXyz(CData):

    x = content("CFloat")
    y = content("CFloat")
    z = content("CFloat")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CXyz.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CXyz with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
)
class CEulerRotation(CData):

    alpha = content("CAngle")
    beta = content("CAngle")
    gamma = content("CAngle")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CEulerRotation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CEulerRotation with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


@cdata_class(
    error_codes={
    },
    qualifiers={
        "allowUndefined": True,
        "guiDefinition": {},
        "saveToDb": False,
    },
    qualifiers_order=[
        'allowUndefined',
        'default',
        'toolTip',
        'guiLabel',
        'guiDefinition',
        'helpFile',
        'saveToDb'],
)
class CTransformation(CData):

    translation = content("CXyz")
    rotation = content("CEulerRotation")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTransformation.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CTransformation with implementation-specific methods.

    Provides flat access to rotation angles (alpha, beta, gamma) and
    translation components (x, y, z) for backward compatibility with
    wrappers that expect direct attribute access (e.g. gesamt).
    """

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

