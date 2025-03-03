"""
     crank2.py: CCP4 GUI Project
     Copyright (C) 2010 University of York, Leiden University

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

from pipelines.crank2.script import crank2_script

class crank2_refatompick(crank2_script.crank2):

  TASKMODULE = 'test'
  TASKNAME = 'crank2_refatompick'
  TASKVERSION = 0.01
  SHORTTASKTITLE = ''

  out_params = ["XYZOUT_SUBSTR", "XYZOUT", "FPHOUT_HL", "FPHOUT_DIFFANOM"]
