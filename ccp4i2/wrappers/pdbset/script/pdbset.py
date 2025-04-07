"""
Copyright (C) 2010 University of York
"""

from ....core.CCP4PluginScript import CPluginScript


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
