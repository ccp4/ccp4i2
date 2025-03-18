"""
Copyright (C) 2014 STFC
Liz Potterton Mar 2014 - performance indicators
"""

from PySide2 import QtCore

from . import CCP4Data, CCP4XtalData
from .CCP4ErrorHandling import CErrorReport, Severity


def performanceIndicatorClasses():
  # This is used by CCP4DbApi.getJobPerformance (and perhaps other methods) to get the
  # appropriate classes from  XData table
  return "('CPerformanceIndicator','CRefinementPerformance','CServalcatPerformance','CModelBuildPerformance','CDataReductionPerformance','CDataReductionCCPerformance','CTestObsConversionsPerformance','CExpPhasPerformance','CPhaseErrorPerformance','CAtomCountPerformance','CPairefPerformance')"

#  ***************  new performance classes also need the keytype to be registered with the database **************************
#   See the definition of KEYTYPELIST in dbapi/CCP4DbApi.py

class CPerformanceIndicator(CCP4Data.CData):
  # Performance indicator expected in wrapper/pipeline output
  # This class should be reimplemented if value is not a float
  CONTENTS_ORDER = ['value','annotation']
  CONTENTS = { 'value' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'annotation' :  { 'class' : CCP4Data.CString } }
  ERROR_CODES = { 300 : { 'description' : 'Passed', 'severity': Severity.OK },
                  301 : { 'description' : 'Data value not set' },
                  302 : { 'description' : 'Performance indicator value difference greater than tolereance' },
                  303 : { 'description' : 'Performance indicator value different' },
                  304 : { 'description' : 'Performance indicator value difference greater than tolereance - but improved', 'severity' : Severity.WARNING },
                  305 : { 'description' : 'Performance indicator not used', 'severity' : Severity.OK }
                  }
  
  def saveToDb(self):
    # Assume performance indicator is saved to db irespective of the saveToDb qualifier
    return [],None,{}



  # To support the QTreeItem.data() calls from CCP4ProjectWidget.CTreeItemJob
  def  data(self,mode=QtCore.Qt.DisplayRole):
    if mode == QtCore.Qt.DisplayRole:
      return self.__str__()
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def scale(self,valueList,absoluteScale=True):
    # Not sure what this does yet.  Probably need to scale a list of PIs of the same type
    # to determine colouring or whatever graphical presentation in the gui
    # Scaling could be absolute of per project
    pass
  
  def assertSame(self,other,testItems=[],diagnostic=False):
    # Return CErrorReport list any failed/warning differences intended for automatic testing
    # other should be another CPerformanceIndicator object
    # It is assumed that self is the original source project
    err = CErrorReport()
    
    # better is integer: 0=null, 1= increasing PI is improvement, -1= decreasing PI is improvement
    for item,tolerance,better in testItems:
      name = self.objectName()+'.'+item
      if diagnostic: print('Comparing ',item,self.__dict__['_value'][item].objectPath(),'with tolerance',tolerance)
      if not self.__dict__['_value'][item].isSet():
        if not other.__dict__['_value'][item].isSet():
            err.append(CPerformanceIndicator,305,name=name,details='in either project',stack=False)
        else:
            err.append(CPerformanceIndicator,301,name=name,details='in original project',stack=False)
      elif not other.__dict__['_value'][item].isSet():
        err.append(CPerformanceIndicator,301,name=name,details='in test project',stack=False)
      else:
        if isinstance(self.__dict__['_value'][item],CCP4Data.CFloat):
          diff = float(other.__dict__['_value'][item])-float(self.__dict__['_value'][item])
          if abs(diff)>tolerance:
            if better!=0 and better*diff>0.0:            
              err.append(CPerformanceIndicator,304,name=name,details= str(self.__dict__['_value'][item])+' : '+str(other.__dict__['_value'][item]),stack=False)
            else:
              err.append(CPerformanceIndicator,302,name=name,details= str(self.__dict__['_value'][item])+' : '+str(other.__dict__['_value'][item]),stack=False)
        elif isinstance(self.__dict__['_value'][item],CCP4Data.CInt):
          diff = int(other.__dict__['_value'][item])-int(self.__dict__['_value'][item])
          if abs(diff)>tolerance:
            if better!=0 and better*diff>0.0:            
              err.append(CPerformanceIndicator,304,name=name,details= str(self.__dict__['_value'][item])+' : '+str(other.__dict__['_value'][item]),stack=False)
            else:
              err.append(CPerformanceIndicator,302,name=name,details= str(self.__dict__['_value'][item])+' : '+str(other.__dict__['_value'][item]),stack=False)
        elif isinstance(self.__dict__['_value'][item],CCP4Data.CString):
          if self.__dict__['_value'][item] != other.__dict__['_value'][item]:
            err.append(CPerformanceIndicator,303,name=name,details= str(self.__dict__['_value'][item])+' : '+str(other.__dict__['_value'][item]),stack=False)
        if len(err) == 0: err.append(CPerformanceIndicator,300,name=name)
        
    return err
  
class CRefinementPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['RFactor','RFree','RMSBond','RMSAngle','weightUsed','annotation']
  CONTENTS = { 'RFactor' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                  'RFree' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                  'RMSBond' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                  'RMSAngle' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                  'weightUsed' : { 'class' : CCP4Data.CFloat },
                  'annotation' :  { 'class' : CCP4Data.CString } }


  def __str__(self):
    text = ''
    if self.__dict__['_value']['RFactor'].__dict__['_value'] is not None:
      text = text + 'R=' + format(self.__dict__['_value']['RFactor'],'.2f')+' '
    if self.__dict__['_value']['RFree'].__dict__['_value'] is not None:
        text = text + 'RFree=' + format(self.__dict__['_value']['RFree'],'.2f')+' '
    return text
    
       
  def saveToDb(self):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(CRefinementPerformance.saveToDb)
    #Return a list of files, string representation of xml, dict of key,value pairs
    #The dict keys must map to those defined in CCP4DbApi.
    ret = {}
    if self.RFactor.isSet(): ret['RFactor'] = self.RFactor.__float__()
    if self.RFree.isSet(): ret['RFree'] = self.RFree.__float__()
    print(ret)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['RFactor',0.01,-1]],diagnostic=diagnostic)


class CServalcatPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['RFactor', 'RFree', 'R', 'R1Factor', 'R1Free', 'R1', 'FSCaverage', 'annotation']
  CONTENTS = { 'RFactor' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'RFree' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'R' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'R1Factor' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'R1Free' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'R1' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'CCFwork_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'CCFfree_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'CCF_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'CCIwork_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'CCIfree_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'CCI_avg' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'FSCaverage' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : -1.0 } },
               'annotation' :  { 'class' : CCP4Data.CString } }

  def __str__(self):
    text = ''
    if self.__dict__['_value']['RFactor'].__dict__['_value'] is not None:
      text = text + 'Rwork=' + format(self.__dict__['_value']['RFactor'], '.4f')+' '
    if self.__dict__['_value']['RFree'].__dict__['_value'] is not None:
        text = text + 'Rfree=' + format(self.__dict__['_value']['RFree'], '.4f')+' '
    if self.__dict__['_value']['R'].__dict__['_value'] is not None:
      text = text + 'R=' + format(self.__dict__['_value']['R'], '.4f')+' '
    if self.__dict__['_value']['R1Factor'].__dict__['_value'] is not None:
      text = text + 'R1work=' + format(self.__dict__['_value']['R1Factor'], '.4f')+' '
    if self.__dict__['_value']['R1Free'].__dict__['_value'] is not None:
        text = text + 'R1free=' + format(self.__dict__['_value']['R1Free'], '.4f')+' '
    if self.__dict__['_value']['R1'].__dict__['_value'] is not None:
      text = text + 'R1=' + format(self.__dict__['_value']['R1'], '.4f')+' '
    if self.__dict__['_value']['CCFwork_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCFwork⟩=' + format(self.__dict__['_value']['CCFwork_avg'], '.4f')+' '
    if self.__dict__['_value']['CCFfree_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCFfree⟩=' + format(self.__dict__['_value']['CCFfree_avg'], '.4f')+' '
    if self.__dict__['_value']['CCF_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCF⟩=' + format(self.__dict__['_value']['CCF_avg'], '.4f')+' '
    if self.__dict__['_value']['CCIwork_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCIwork⟩=' + format(self.__dict__['_value']['CCIwork_avg'], '.4f')+' '
    if self.__dict__['_value']['CCIfree_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCIfree⟩=' + format(self.__dict__['_value']['CCIfree_avg'], '.4f')+' '
    if self.__dict__['_value']['CCI_avg'].__dict__['_value'] is not None:
      text = text + '⟨CCI⟩=' + format(self.__dict__['_value']['CCI_avg'], '.4f')+' '
    if self.__dict__['_value']['FSCaverage'].__dict__['_value'] is not None:
      text = text + '⟨FSCmodel⟩=' + format(self.__dict__['_value']['FSCaverage'], '.4f')+' '
    return text

  def saveToDb(self):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(CServalcatPerformance.saveToDb)
    #Return a list of files, string representation of xml, dict of key,value pairs
    #The dict keys must map to those defined in CCP4DbApi.
    ret = {}
    if self.RFactor.isSet(): ret['RFactor'] = self.RFactor.__float__()
    if self.RFree.isSet(): ret['RFree'] = self.RFree.__float__()
    if self.R.isSet(): ret['R'] = self.R.__float__()
    if self.R1Factor.isSet(): ret['R1Factor'] = self.R1Factor.__float__()
    if self.R1Free.isSet(): ret['R1Free'] = self.R1Free.__float__()
    if self.R1.isSet(): ret['R1'] = self.R1.__float__()
    if self.CCFwork_avg.isSet(): ret['CCFwork_avg'] = self.CCFwork_avg.__float__()
    if self.CCFfree_avg.isSet(): ret['CCFfree_avg'] = self.CCFfree_avg.__float__()
    if self.CCF_avg.isSet(): ret['CCF_avg'] = self.CCF_avg.__float__()
    if self.CCIwork_avg.isSet(): ret['CCIwork_avg'] = self.CCIwork_avg.__float__()
    if self.CCIfree_avg.isSet(): ret['CCIfree_avg'] = self.CCIfree_avg.__float__()
    if self.CCI_avg.isSet(): ret['CCI_avg'] = self.CCI_avg.__float__()
    if self.FSCaverage.isSet(): ret['FSCaverage'] = self.FSCaverage.__float__()
    print(ret)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['RFactor',0.01,-1]],diagnostic=diagnostic)


class CPairefPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['cutoff']
  CONTENTS = { 'cutoff' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } } }

  def __str__(self):
    text = ''
    if self.__dict__['_value']['cutoff'].__dict__['_value'] is not None:
      text = text + 'Cutoff=' + format(self.__dict__['_value']['cutoff'],'.2f')+' '
    return text

  def saveToDb(self):
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(CPairefPerformance.saveToDb)
    #Return a list of files, string representation of xml, dict of key,value pairs
    #The dict keys must map to those defined in CCP4DbApi.
    ret = {}
    if self.cutoff.isSet(): ret['cutoff'] = self.cutoff.__float__()
    print(ret)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['cutoff',0.0001,-1]],diagnostic=diagnostic)

class CModelBuildPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['RFactor','completeness','annotation']
  CONTENTS = { 'RFactor' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'completeness' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'annotation' :  { 'class' : CCP4Data.CString } }

  def  __str__(self):
    text = ''
    if self.__dict__['_value']['RFactor'].__dict__['_value'] is not None:
      text = text + 'R=' + format(self.__dict__['_value']['RFactor'],'0.2f')+' '
    if self.__dict__['_value']['completeness'].__dict__['_value'] is not None:
      text = text + ' %=' + format(int(self.__dict__['_value']['completeness']*100),'d')
    return text

       
  def saveToDb(self):
    #Return a list of files, string reporesentation of xml, dict of key,value pairs
    #The dict keys must map to does defined in CCP4DbApi.
    #print 'CModelBuildPerformance.saveToDb',{'RFactor' : self.RFactor.__float__(), 'completeness' : self.completeness.__float__() }
    ret = {}
    if self.RFactor.isSet(): ret['RFactor'] = self.RFactor.__float__()
    if self.completeness.isSet(): ret['completeness'] = self.completeness.__float__()
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['RFactor',0.01,-1],['completeness',0.01,1]],diagnostic=diagnostic)


class CDataReductionPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = [ 'spaceGroup', 'highResLimit', 'rMeas' ]
  CONTENTS = { 'spaceGroup' : { 'class' : CCP4XtalData.CSpaceGroup },
               'highResLimit' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'rMeas' :  { 'class' : CCP4Data.CFloat }
                }

  def __str__(self):
    text = ''
    if self.__dict__['_value']['spaceGroup'].isSet(): text += 'Sgp='+str(self.__dict__['_value']['spaceGroup']) 
    if self.__dict__['_value']['highResLimit'].isSet():  text += ' res=' +  format(self.__dict__['_value']['highResLimit'],'.2f')
    if self.__dict__['_value']['rMeas'].isSet():  text += ' Rmeas=' + format(self.__dict__['_value']['rMeas'],'.3f')
    return text


  def saveToDb(self):
    ret = { }
    if self.__dict__['_value']['spaceGroup'].isSet(): ret['spaceGroup'] = str(self.__dict__['_value']['spaceGroup'])
    if self.__dict__['_value']['highResLimit'].isSet(): ret['highResLimit'] = self.highResLimit.__float__()
    if self.__dict__['_value']['rMeas'].isSet(): ret['rMeas'] = self.rMeas.__float__() 
    #print 'CDataReductionPerformance.saveToDb',ret
    return  [],None,ret
  
  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['rMeas',0.01,-1],['highResLimit',0.01,0],['spaceGroup','',0]],diagnostic=diagnostic)


class CDataReductionCCPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = [ 'spaceGroup', 'highResLimit', 'ccHalf' ]
  CONTENTS = { 'spaceGroup' : { 'class' : CCP4XtalData.CSpaceGroup },
               'highResLimit' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'ccHalf' :  { 'class' : CCP4Data.CFloat }
                }

  def __str__(self):
    text = ''
    if self.__dict__['_value']['spaceGroup'].isSet(): text += 'Sgp='+str(self.__dict__['_value']['spaceGroup']) 
    if self.__dict__['_value']['highResLimit'].isSet():  text += ' res=' +  format(self.__dict__['_value']['highResLimit'],'.2f')
    if self.__dict__['_value']['ccHalf'].isSet():  text += ' CC1/2=' + format(self.__dict__['_value']['ccHalf'],'.3f')
    return text


  def saveToDb(self):
    ret = { }
    if self.__dict__['_value']['spaceGroup'].isSet(): ret['spaceGroup'] = str(self.__dict__['_value']['spaceGroup'])
    if self.__dict__['_value']['highResLimit'].isSet(): ret['highResLimit'] = self.highResLimit.__float__()
    if self.__dict__['_value']['ccHalf'].isSet(): ret['ccHalf'] = self.ccHalf.__float__() 
    #print 'CDataReductionPerformance.saveToDb',ret
    return  [],None,ret
  
  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['ccHalf',0.01,-1],['highResLimit',0.01,0],['spaceGroup','',0]],diagnostic=diagnostic)


class CAtomCountPerformance(CPerformanceIndicator):
  # Trivial performance test for chainsaw/sculptor for use in project-based-testing
  CONTENTS_ORDER = [ 'nAtoms', 'nResidues']
  CONTENTS = { 'nAtoms' :  { 'class' : CCP4Data.CInt },
               'nResidues' :  { 'class' : CCP4Data.CInt } }
 
  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[ ['nAtoms',0,1 ] ] ,diagnostic=diagnostic)
  
  def __str__(self):
    text = ''
    if self.__dict__['_value']['nResidues'].isSet():
      text = text + 'nRes='+str(self.__dict__['_value']['nResidues'])
    elif self.__dict__['_value']['nAtoms'].isSet():
      text = text + 'nAtoms='+str(self.__dict__['_value']['nAtoms'])
    return text

  def saveToDb(self):
    ret = {}
    if self.__dict__['_value']['nAtoms'].isSet(): ret['nAtoms'] = int(self.__dict__['_value']['nAtoms'])
    if self.__dict__['_value']['nResidues'].isSet(): ret['nResidues'] = int(self.__dict__['_value']['nResidues'])
    print('CAtomCountPerformance.saveToDb',ret)
    return  [],None,ret

  def testComparisonData(self):   
    ret = { }
    if self.__dict__['_value']['nAtoms'].isSet(): ret['nAtoms'] = self.__dict__['_value']['nAtoms']
    if self.__dict__['_value']['nResidues'].isSet(): ret['nResidues'] = self.__dict__['_value']['nResidues']
    return  [],None,ret

  def setFromPdbDataFile(self,fileName):
    if  isinstance(fileName,str):
      from . import CCP4ModelData
      fileName = CCP4ModelData.CPdbDataFile(fileName)
    try:
      fileName.loadFile()
      print('CAtomCountPerformance.setFromPdbDataFile',fileName.fileContent.composition.nAtoms)
      self.__dict__['_value']['nAtoms'].set(fileName.fileContent.composition.nAtoms)
      self.__dict__['_value']['nResidues'].set(fileName.fileContent.composition.nResidues)
    except:
      pass
      

class CTestObsConversionsPerformance(CPerformanceIndicator):
  # Trivial performance test for chainsaw/sculptor for use in project-based-testing
  CONTENTS_ORDER = [ 'columnLabelsString']
  CONTENTS = { 'columnLabelsString' :  { 'class' : CCP4Data.CString } }
 
  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[ ['columnLabelsString',0,1 ] ] ,diagnostic=diagnostic)
  
  def __str__(self):
    text = ''
    if self.__dict__['_value']['columnLabelsString'].isSet():
      text = text + 'labels:'+str(self.__dict__['_value']['columnLabelsString']).replace("'","")
    return text

  def saveToDb(self):
    ret = {}
    if self.__dict__['_value']['columnLabelsString'].isSet(): ret['columnLabelsString'] = str(self.__dict__['_value']['columnLabelsString'])
    print('CTestObsConversionsPerformance.saveToDb',ret)
    return  [],None,ret

  def testComparisonData(self):   
    ret = { }
    if self.__dict__['_value']['columnLabelsString'].isSet(): ret['columnLabelsString'] = self.__dict__['_value']['columnLabelsString']
    return  [],None,ret


class CExpPhasPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['FOM','CFOM','Hand1Score','Hand2Score','CC','RFactor','RFree','annotation']
  CONTENTS = { 'FOM' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'CFOM' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'Hand1Score' : { 'class' : CCP4Data.CFloat }, 'Hand2Score' : { 'class' : CCP4Data.CFloat },
                'CC' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'RFactor' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'RFree' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
                'annotation' :  { 'class' : CCP4Data.CString } }

  def __str__(self):
    text, val = '', self.__dict__['_value']
    if val['FOM'].isSet():
      text += 'FOM={:.2f} '.format(val['FOM'])
    if val['CFOM'].isSet():
      text += 'CFOM={:.1f} '.format(val['CFOM'])
    if val['Hand1Score'].isSet() and val['Hand2Score'].isSet():
      if val['CC'].isSet() and min(abs(val['CC']-val['Hand1Score']),abs(val['CC']-val['Hand2Score']))<1e-4:
        text += 'CC={:.0f}:{:.0f} '.format(val['Hand1Score'],val['Hand2Score'])
      else:
        text += 'Score={:.0f}:{:.0f} '.format(val['Hand1Score'],val['Hand2Score'])
    elif val['CC'].isSet():
      text += 'CC={:.1f} '.format(val['CC'])
    if val['RFactor'].isSet():
      text += 'R={:.2f} '.format(val['RFactor'])
    if val['RFree'].isSet():
      text += 'Rfree={:.2f} '.format(val['RFree'])
    return text

  def saveToDb(self):
    #Return a list of files, string representation of xml, dict of key,value pairs
    #The dict keys must map to those defined in CCP4DbApi.
    ret = {}
    if self.FOM.isSet(): ret['FOM'] = self.FOM.__float__()
    if self.CFOM.isSet(): ret['CFOM'] = self.CFOM.__float__()
    if self.Hand1Score.isSet(): ret['Hand1Score'] = self.Hand1Score.__float__()
    if self.Hand2Score.isSet(): ret['Hand2Score'] = self.Hand2Score.__float__()
    if self.CC.isSet(): ret['CC'] = self.CC.__float__()
    if self.RFactor.isSet(): ret['RFactor'] = self.RFactor.__float__()
    if self.RFree.isSet(): ret['RFree'] = self.RFree.__float__()
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['FOM',0.05,-1],['CFOM',10.,-1],['Hand1Score',10.,-1],['Hand2Score',10.,-1],['CC',0.01,-1],['RFactor',0.01,-1],['RFree',0.01,-1],],diagnostic=diagnostic)


class CPhaseErrorPerformance(CPerformanceIndicator):
  CONTENTS_ORDER = ['phaseError','weightedPhaseError','reflectionCorrelation']
  CONTENTS = { 'phaseError' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'weightedPhaseError' : { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } },
               'reflectionCorrelation' :  { 'class' : CCP4Data.CFloat, 'qualifiers' : { 'min' : 0.0 } } }

  def  __str__(self):
    #if self.__dict__['_value']['phaseError'].__dict__['_value'] is not None:
    text = 'dphi=' + format(self.__dict__['_value']['phaseError'],'.1f')
    return text

  def saveToDb(self):
    #Return a list of files, string reporesentation of xml, dict of key,value pairs
    #The dict keys must map to does defined in CCP4DbApi.
    ret = {}
    if self.phaseError.isSet(): ret['phaseError'] = self.phaseError.__float__()
    if self.weightedPhaseError.isSet(): ret['weightedPhaseError'] = self.weightedPhaseError.__float__()
    if self.reflectionCorrelation.isSet(): ret['reflectionCorrelation'] = self.reflectionCorrelation.__float__()
    return [],None,ret

  def assertSame(self,other,diagnostic=False):
    return CPerformanceIndicator.assertSame(self,other,[['phaseError',5.0,-1],['weightedPhaseError',5.0,-1],['reflectionCorrelation',0.05,1]],diagnostic=diagnostic)



class CSuperposePerformance(CPerformanceIndicator):
  # Trivial performance test for chainsaw/sculptor for use in project-based-testing
  CONTENTS_ORDER = [ 'RMSxyz' ,'nResidues' ]
  CONTENTS = { 'RMSxyz' :  { 'class' : CCP4Data.CFloat } ,
               'nResidues' :  { 'class' : CCP4Data.CInt } }
             #  'QScore' :  { 'class' : CCP4Data.CFloat },
 
  def assertSame(self,other,diagnostic=False):
    #return CPerformanceIndicator.assertSame(self,other,[ ['nResidues',0,1 ], [ 'QScore', 1.0, 1], ['RMSxyz', 0.1, -1 ]  ] ,diagnostic=diagnostic)
    return CPerformanceIndicator.assertSame(self,other,[  ['RMSxyz', 0.1, -1 ],  ['nResidues',0,1 ]] ,diagnostic=diagnostic)
  
  def __str__(self):
    text = ''
    if self.__dict__['_value']['RMSxyz'].isSet():
      text = text + 'RMS='+str(self.__dict__['_value']['RMSxyz']) + ' '
    if self.__dict__['_value']['nResidues'].isSet():
      text = text + 'nRes='+str(self.__dict__['_value']['nResidues']) + ' '
    return text

  def saveToDb(self):
    ret = {}
    if self.__dict__['_value']['RMSxyz'].isSet(): ret['RMSxyz'] = float(self.__dict__['_value']['RMSxyz'])
    if self.__dict__['_value']['nResidues'].isSet(): ret['nResidues'] = int(self.__dict__['_value']['nResidues'])
    #if self.__dict__['_value']['QScore'].isSet(): ret['QScore'] = float(self.__dict__['_value']['QScore'])
    return  [],None,ret
