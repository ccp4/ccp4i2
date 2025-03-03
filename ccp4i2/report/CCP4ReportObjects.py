"""
     CCP4ReportObjects.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

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

"""
   Liz Potterton Aug 2011 -Kludged report objects
"""
from core import CCP4Data
from core.CCP4ErrorHandling import *
import types


    
class CReportAnnotation(CCP4Data.CData):
  """Annotation text"""

  QUALIFIERS = { 'text' : '',
                 'program' : '',
                 'tag' : '',
                 'style' : 'normal' } 
  QUALIFIERS_ORDER = ['text','program','tag','style']
  QUALIFIERS_DEFINITION = { 'text' : { 'type' : types.StringType,
                                       'description' : 'Fixed text' },
                            'program' : {  'type' : types.StringType,
                                       'description' : 'Program log file' },
                            'tag' : {  'type' : types.StringType,
                                       'description' : 'Log file tag' },
                            'style' : { 'type' : 'enumerator',
                                       'menu' : ['normal','small','large'],
                                       'description' : 'Text style' }
                            }

class CReportMiniTable(CCP4Data.CData):
  '''Mini table'''
  QUALIFIERS = { 'title' : '',
                 'program' : '',
                 'tag' : '' }
   
  QUALIFIERS_ORDER = ['title','program','tag']
  
  QUALIFIERS_DEFINITION = { 'title' : { 'type' : types.StringType,
                                       'description' : 'Fixed text' },
                            'program' : {  'type' : types.StringType,
                                       'description' : 'Program log file' },
                            'tag' : {  'type' : types.StringType,
                                       'description' : 'Log file tag' } }
