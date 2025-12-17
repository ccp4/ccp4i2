import types

from ccp4i2.core import CCP4Data


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
