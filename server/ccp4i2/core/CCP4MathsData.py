# Copyright (C) 2025-2026 University of York
# Copyright (C) 2025-2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""
Implementation classes for CCP4MathsData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any

from ccp4i2.core.cdata_stubs.CCP4MathsData import CAngleStub, CEulerRotationStub, CMatrix33Stub, CTransformationStub, CXyzStub, CXyzBoxStub


class CAngle(CAngleStub):
    """
    An angle
    
    Extends CAngleStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CEulerRotation(CEulerRotationStub):
    """
    Extends CEulerRotationStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMatrix33(CMatrix33Stub):
    """
    Extends CMatrix33Stub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CTransformation(CTransformationStub):
    """
    Extends CTransformationStub with implementation-specific methods.

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


class CXyz(CXyzStub):
    """
    Extends CXyzStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CXyzBox(CXyzBoxStub):
    """
    Extends CXyzBoxStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

