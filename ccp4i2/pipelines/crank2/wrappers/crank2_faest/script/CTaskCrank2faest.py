"""
     CTaskCrank2faest.py
     Copyright (C) 2011 University of York, Leiden University

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

from pipelines.crank2.script import CTaskCrank2

class CTaskCrank2faest(CTaskCrank2.CTaskCrank2):
  TASKNAME = 'crank2_faest'
  TASKVERSION = 0.01
  TASKMODULE='test'
  TASKTITLE='FA estimation'
  SHORTTASKTITLE = None
  DESCRIPTION='Crank2 FA estimation and other phasing preparations'

def whatNext(*args,**kwargs):
  return CTaskCrank2.whatNext(*args,**kwargs)
