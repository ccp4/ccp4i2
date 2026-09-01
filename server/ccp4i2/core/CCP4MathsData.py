"""
Implementation classes for CCP4MathsData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType
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

    xMin: Optional["CFloat"] = None
    yMin: Optional["CFloat"] = None
    zMin: Optional["CFloat"] = None
    xMax: Optional["CFloat"] = None
    yMax: Optional["CFloat"] = None
    zMax: Optional["CFloat"] = None

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

    x: Optional["CFloat"] = None
    y: Optional["CFloat"] = None
    z: Optional["CFloat"] = None

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
    contents_order=['alpha', 'beta', 'gamma'],
)
class CEulerRotation(CData):

    alpha: Optional["CAngle"] = None
    beta: Optional["CAngle"] = None
    gamma: Optional["CAngle"] = None

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
    contents_order=['translation', 'rotation'],
)
class CTransformation(CData):

    translation: Optional["CXyz"] = None
    rotation: Optional["CEulerRotation"] = None

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

