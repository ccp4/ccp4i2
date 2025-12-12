from __future__ import print_function

"""
     pdbset.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

from ccp4i2.core.CCP4PluginScript import CPluginScript
     
class pdbset(CPluginScript):

    TASKMODULE = 'demo'
    TASKTITLE = 'PDBSet'
    TASKNAME = 'pdbset'
    TASKCOMMAND = 'pdbset'
    TASKVERSION= 0.0
    COMLINETEMPLATE = '''1 XYZIN $XYZIN
1 XYZOUT $XYZOUT'''
    COMTEMPLATE = '''1 CELL $CELL.a $CELL.b $CELL.c $CELL.alpha $CELL.beta $CELL.gamma
1 END'''
