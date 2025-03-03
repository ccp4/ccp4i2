from __future__ import print_function

"""
     tasks/dummy.py: CCP4 GUI Project
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

"""
     Liz Potterton Jan 2010 - Create demo.py prototype
"""

from PySide2 import QtCore

from qtgui.CCP4TaskWidget import CTaskWidget

#-------------------------------------------------------------------
class CTaskDummy(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'dummy'
  TASKVERSION = 0.0
  TASKMODULE='demo'
  TASKTITLE='Task interface demo code'
  DESCRIPTION = 'Demo code for developers - does not run'
  SHORTTASKTITLE = 'DumMee'

  ERROR_CODES = { 333 : { 'description' : 'Space group in coordinate and reflection file do not match' } }

  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)


  def drawContents(self):
                            
    self.openFolder(folderFunction='inputData')

    self.createLine ( [ 'widget','HKLIN' ] )

    # Make the behaviour of PDBIN dependent on the checkbox of PDBIN_COMPULSARY
    # The method handlePDBIN_COMPULSARY() is defined below to set the qualifiers of PDBIN appropriately
    self.createLine ( [ 'widget' , 'PDBIN_COMPULSARY' , 'label' , 'Model data input is compulsary'])
    self.createLine ( [  'widget','PDBIN' ] )
    self.connectDataChanged('PDBIN_COMPULSARY',self.handlePDBIN_COMPULSARY)
    self.handlePDBIN_COMPULSARY()

    self.createLine ( [ 'advice' , 'Enter cell parameters' ] )
    self.createLine ( [ 'widget','CELL' ] )
    
    
    self.openFolder(title='Simple Options')

    self.createLine ( [ 'label' , 'Run for',
                        'widget', 'NCYCLES' ,
                        'label' , 'cycles and use',
                        'widget' , 'QUICKNDIRTY',
                        'label', 'evaluation method' ] )


    # The WILLITWORK checkbox controls display of the WILLITWORKTIMES widget
    # and the line with TESTINPUT.

    # This line is displayed if WILLITWORK is False
    self.createLine( [ 'help', 'WILLITWORK',
                       'tip','Click to see line toggled',
                       'widget', 'WILLITWORK',
                       'label','Will this work',
                       ], toggle=['WILLITWORK','open', [False]] )
    
    # Example of using openSubFrame to group lines that are displayed/hidden together
    # If WILLITWORK is True then the subframe is shown - this has contains two lines
    # The first line is the preiously defined line plus the WILLITWORKTIMES widget
    self.openSubFrame( toggle=['WILLITWORK','open', [True]] )
    self.createLine( [ 'help', 'WILLITWORK',
                       'tip','Click to see line toggled',
                       'widget', 'WILLITWORK',
                       'label','Will this work, how many times',
                       'widget','WILLITWORKTIMES'
                        ] )

    self.createLine ( [  'help', 'TESTING',
                        'label', 'Testing..testing..',
                        'widget', 'TESTINPUT'
                        ] )

    self.closeSubFrame()
    

    # The DOHARDSTUFF widget controls whether the 'Advanced Options' folder
    # is accessible 
    self.createLine( [ 'tip', 'Click to see folder toggled',
                       'widget', 'DOHARDSTUFF',
                       'label', 'Do something really difficult' ] )
    

    self.openFolder(title='Advanced options',toggle=['DOHARDSTUFF','open',[1]])

    # Set the menu text - do not need to give menuText for all enumerators
    #self.setMenuText('CHOOSEONE',{None:'program decides',1:'one',3:'three'})
    self.createLine( [  'label','Choose one',
                       'widget','CHOOSEONE',
                       'label','or choose another',
                       'widget', 'CHOOSEOTHER'
                        ] )
    
    # Group the refinement options in a framed (ie outlined) sub-frame
    self.openSubFrame(frame=True,title='Customise the refinement')
    # Show REFINE_MODE options as radio buttons with one button per line
    self.setMenuText( 'REFINE_MODE' , { 'simple' : 'Keep it simple', 'custom' : 'Customise' } )
    self.createLine( ['widget','-guiMode','multiLineRadio','REFINE_MODE'])

    # Control display of all REF_PARAM parameters together by putting then
    # in the same sub-frame
    # Set menuText and toolTip for each parameter
    #self.openSubFrame(toggle=['REFINE_MODE','open',['custom']])
    self.setMenuText('REF_PARAM_A',['program decides','normal','few','many'])
    self.setToolTip('REF_PARAM_A','Enter a number')
    self.createLine( ['label','Refinement parameter A','widget','REF_PARAM_A'],toggle=['REFINE_MODE','open',['custom']])
    self.setToolTip('REF_PARAM_B','Choose a method')
    self.setMenuText('REF_PARAM_B',['Fast method','Slow method'])
    self.createLine( ['label','Refinement parameter B','widget','REF_PARAM_B'],toggle=['REFINE_MODE','open',['custom']])
    self.createLine( ['label','Refinement parameter C','widget','REF_PARAM_C'],toggle=['REFINE_MODE','open',['custom']])

    self.closeSubFrame()
   
    self.openFolder(title='Lists')
    self.createLine( [ 'advice' , 'Example of a simple list of strings' ] )
    self.createLine( ['widget' ,'-title','Enter a list of strings - whatever you like!!', 'MRFILELIST'])


    self.createLine( [ 'advice' , 'Example of list of CAsuComponents' ] )
    self.createLine( ['widget' ,'-title','Contents of asymmetric unit',  'ASUCOMPONENTLIST'])

    self.openFolder(title='Split lines')
    # Vsibility of items on following line dependent on IFWHAT and IFEVER
    # Note that in this example we consider these parameters are not used by the pipeline/wrapper
    # so that, in the def file, they are placed in the guiParameters sub-container.
    self.createLine( ['widget' ,'IFWHAT','label','Show what','widget' ,'IFEVER','label','Show ever'] )
    # Create a blank line (or with a spacer item)
    line = self.createLine(['stretch'])
    # Create two line objects that are appended to the preceeding line.
    # Visibility of these line fragments is dependent on  IFWHAT or IFEVER
    self.createLine( ['label','What is What','widget','WHAT'],appendLine=line,toggle=['IFWHAT','open'])
    self.createLine( ['label','What is Ever','widget','EVER'],appendLine=line,toggle=['IFEVER','open'])
    # The following line fragment is displayed if either IFWHAT or IFEVER are true
    # The function whatEverVisibility controls the visibility.
    # The second part of the toggleFunction argument is a list of the parameters that the toggleFunction
    # is dependent on. This list is needed to ensure that things get updated whenever any of these parameters change.
    self.createLine( ['label','and how about Whatever','widget','WHATEVER'],appendLine=line,
                     toggleFunction=[self.whatEverVisibility,['IFEVER','IFWHAT']])

    # Demo of button to launch another 

    self.openFolder(title='New job')
    self.createLine(['label','Enter a coordinate file to test launch job button'])
    self.createLine(['widget','XYZIN'])
    self.createLine(['label','Superpose the coordinates on another model by','launchButton','gesamt'])
    # Update the status of the lauch button dependent on selection of XYZIN
    self.connectDataChanged('XYZIN',self.updateLauchButton)
    # .. and initialise its status
    self.updateLauchButton()
    
    self.openFolder(title='Autogenerated')
    self.autoGenerate(container=self.container.controlParameters,selection={ 'includeParameters' : ['OBSCURE_A','OBSCURE_B'] } )
    
  def whatEverVisibility(self):
    # Evaluate whether the WHATEVER widget should be visible
    # Dependent on two parameters
    if self.container.guiParameters.IFWHAT or self.container.guiParameters.IFEVER:
      return True
    else:
      return False

  @QtCore.Slot()
  def updateLauchButton(self):
    # Enable the launch button dependent on whether XYZIN is set
    self.findWidget('gesamt').setEnabled(self.container.inputData.XYZIN.isSet())

  def handleLaunchedJob(self,jobId=None,status=None,taskWidget=None):
    print('CTaskDummy.handleLaunchedJob',jobId,status,taskWidget)
    # If this is called with status=1=Pending the taskWidget (gesamt window) has just been opened
    # and we can set a value in it
    if status == 1 and taskWidget is not None:
      taskWidget.container.inputData.XYZIN_QUERY.set(self.container.inputData.XYZIN)
    # Status is 6 (Finished)
    elif status == 6:
      # Can not assume that the gesamt widget is still there - must instead query the database for output file
      # Use CDbApi.getJobFilesInfo() which returns a list of dicts containing description of files output by the job
      # The best way to set the file object ot a new value is by setDbFileId()
      from core import CCP4Modules
      gesamtFileList = CCP4Modules.PROJECTSMANAGER().db().getJobFilesInfo(jobId=jobId,jobParamName='XYZOUT')
      #print 'CTaskDummy.handleLaunchedJob ',gesamtFileList
      if len(gesamtFileList)>0:
        self.getWidget('XYZIN').model.setDbFileId(gesamtFileList[0]['fileId'])
        
      
  @QtCore.Slot()
  def handlePDBIN_COMPULSARY(self):
    # set the qualifiers of PDBIN dependent on the value of PDBIN_COMPULSARY
    # see ccp4/docs/developer/task_guis.html#setQualifiers for discussion
    mode = bool(self.container.inputData.PDBIN_COMPULSARY)
    self.container.inputData.PDBIN.setQualifiers({ 'allowUndefined' : (not mode) } )
    self.getWidget('PDBIN').validate()


  def taskValidity(self):
    from core import CCP4ErrorHandling
    rv = CCP4ErrorHandling.CErrorReport()
    # Check the space group is same in MTZ and PDB
    if self.container.inputData.PDBIN.exists():
      pdbSpgp = self.container.inputData.PDBIN.fileContent.mmdbManager.GetSpaceGroup()
    else:
      pdbSpgp = None
    if self.container.inputData.HKLIN.exists():
      mtzSpgp = str(self.container.inputData.HKLIN.fileContent.spaceGroup)
    else:
      mtzSpgp = None

    print('CTaskDummy.taskValidity',pdbSpgp,mtzSpgp)
    if pdbSpgp is not None and mtzSpgp is not None and pdbSpgp != mtzSpgp:
        rv.append(self.__class__,333,details='Coordinate file space group:'+str(pdbSpgp)+' Reflections file space group:'+str(mtzSpgp),stack=False)

    return rv
