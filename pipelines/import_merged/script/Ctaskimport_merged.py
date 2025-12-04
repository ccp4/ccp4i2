from core.CCP4ErrorHandling import *
from qtgui import CCP4XtalWidgets
from qtgui.CCP4TaskWidget import CTaskWidget

import gemmi
from baselayer import QtCore, QtWidgets

from  pipelines.import_merged.script.mmcifutils import *
from  pipelines.import_merged.script.dybuttons import *
from  pipelines.import_merged.script.importutils import ReflectionDataTypes

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  IMPORT_MERGED
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
class CTaskimport_merged(CTaskWidget):
#-------------------------------------------------------------------

# Subclass CTaskWidget to give specific task window
  TASKNAME = 'import_merged'
  TASKVERSION = 0.1
  TASKMODULE='data_entry'
  TASKTITLE='Import merged reflection data'
  SHORTTASKTITLE='Import merged'
  DESCRIPTION = 'Import reflection data in any format, report on contents and create CCP4i2 data objects'

  # -------------------------------------------------------------
  def __init__(self,parent):
    CTaskWidget.__init__(self,parent)

  # -------------------------------------------------------------
  def drawContents(self):
    print("Ctaskimport_merged.py drawContents")
    self.container.guiParameters.HKLIN_HAS_COLUMNS.set(False)

    if self.container.inputData.SPACEGROUPCELL.isSet():
      #print("SPACEGROUPCELL set")
      #   SPACEGROUPCELL object exists in old cloned jobs only
      # copy information to separate objects, new style (PRE February 2024)
      self.container.inputData.UNITCELL.set\
          (self.container.inputData.SPACEGROUPCELL.cell)
      self.container.inputData.SPACEGROUP.set\
          (self.container.inputData.SPACEGROUPCELL.spaceGroup)
      self.container.inputData.SPACEGROUPCELL.unSet()

    self.openFolder(folderFunction='inputData',followFrom=False)

    self.createLine ( [ 'subtitle','Select a merged data file' ] )
    self.createLine ( [ 'widget','HKLIN' ] )

    #  Extra information wanted for non-MTZ files
    self.openSubFrame(toggle=[ 'HKLINISMTZ', 'open', [ False ] ] )
    self.createLine ( [ 'subtitle','Enter additional information' ] )
    self.createLine ( [  'widget', 'SPACEGROUP', 'widget', 'UNITCELL' ] )
    self.createLine ( [ 'label','Wavelength','widget', 'WAVELENGTH' ] )
    self.createLine ( [ 'label', 'Crystal name','widget','CRYSTALNAME','label', 'Dataset name','widget','DATASETNAME' ] )

    self.openSubFrame(toggle=[ 'SHOW_MMCIF_BLOCKS', 'open', [ True ] ] )
    self.setMMCIFframe()
    self.closeSubFrame()  #  MMCIF frame
    self.closeSubFrame()  #  non-MTZ frame


    self.createLine (['label',''])
    # Selected columns
    self.columnAdviceMTZ   = 'reselect input file to change'
    self.columnAdviceMMCIF = ' reselect cif reflection block to change content selection'

    advlabel = ''
    clabel = ''

    self.container.guiParameters.HKLIN_HAS_COLUMNS.set(False)
    if self.container.inputData.HKLIN_FORMAT == 'MTZ':
      self.container.guiParameters.HKLIN_HAS_COLUMNS.set(True)
      clabel =  ' MTZ columns: ' + str(self.container.inputData.HKLIN_OBS_COLUMNS)
      advlabel = self.columnAdviceMTZ

    elif self.container.inputData.HKLIN_FORMAT == 'MMCIF':
      self.container.guiParameters.HKLIN_HAS_COLUMNS.set(True)
      content = str(self.container.inputData.MMCIF_SELECTED_CONTENT)
      block = str(self.container.inputData.MMCIF_SELECTED_BLOCK)
      clabel =  ' MMCIF Block: ' + block + ', content type: ' + content
      advlabel = self.columnAdviceMMCIF
      
    colLine = self.createLine( ['label','Selected data: ',
                                'label', clabel,
                                'advice', advlabel],
                               toggle=['HKLIN_HAS_COLUMNS','open',[True]])
    self.selectedColumnLabels = colLine.layout().itemAt(1).widget()
    self.selectedColumnAdvice = colLine.layout().itemAt(2).widget()
    
    #print(">*>* Format", self.container.inputData.HKLIN_FORMAT,
    #    self.container.guiParameters.HKLIN_HAS_COLUMNS)

    if self.container.inputData.HKLIN_FORMAT == 'MMCIF':
      self.setCIFblocklist()
    
    self.createLine (['label',
                      'This file appears to be from the StarAniso server'],
                     toggle=['STARANISO_DATA', 'closed', [False]])
    self.createLine (['label',
                      'By default the FreeR flag will be imported unchanged without completion'],
                     toggle=['STARANISO_DATA', 'closed', [False]])
    
    # Optional resolution cutoff
    line = self.createLine( ['label','Select resolution range (\xc5)',
                             'widget', 'RESOLUTION_RANGE',
                             'advice',' Highest resolution in file ',
                             'label','range to be determined'],
                            toggle=['CAN_CUT_RESOLUTION', 'closed', [False]])
    # get label widget, note itemAt counts from 0
    self.maximum_resolution_label = line.layout().itemAt(3).widget()
    self.container.inputData.RESOLUTION_RANGE.dataChanged.connect(self.resorangeset)

    self.createLine (['label',
       'Any FreeR in the input data will be automatically read and completed'],
                     toggle=['STARANISO_DATA', 'open', [False]])
    self.createLine (['label',
       'If no FreeR is present, a new set will be generated, or an existing FreeR set may be defined below'],
                     toggle=['STARANISO_DATA', 'open', [False]])
    self.createLine (['label',''])
    self.createLine (['label',
       'You can define a pre-existing FreeR set below to override a FreeR set from the input data'])
    self.createLine (['subtitle','If a FreeR set is neither present in the input, nor given explicitly, then one will be generated'])

    self.createLine( ['widget', 'FREERFLAG'] )
    self.container.inputData.FREERFLAG.dataChanged.connect(self.handleSelectFreeRflag)
        
    self.createLine(['widget','CUTRESOLUTION',
                     'label',
                     'Cut resolution of FreeR set if necessary to match the data'],
                    toggleFunction=[self.toggleCutResolution, ['FREERFLAG']])
    self.createLine(['label', 'Fraction of reflections in generated freeR set',
                     'widget','FREER_FRACTION',
                     'advice', 'Default fraction is 0.05',
                     'advice',
                     'Potential twinning operations will be taken into account'],
                    toggleFunction=[self.toggleFraction, ['FREERFLAG', 'HASFREER']])
    self.createLine(['widget', 'SKIP_FREER',
                     'label', 'Leave input FreeR set unchanged'],
                    toggleFunction=[self.toggleSkip, ['FREERFLAG', 'HASFREER', 'STARANISO_DATA']])
    self.createLine(['widget', 'SKIP_FREER',
                     'label', 'Leave input FreeR set unchanged (recommended for StarAniso data)'],
                    toggleFunction=[self.toggleSkip2, ['FREERFLAG', 'HASFREER', 'STARANISO_DATA']])

    # Parameters for Free R
    self.openSubFrame(frame=True,
                      toggle=['COMPLETE', 'open', [True]])
    self.createLine( ['subtitle','Options for FreeR set extension:  -------'])

    self.createLine( ['advice','If you are extending an existing FreeR set, it must match the observed data in unit cell and Laue group'])
    self.createLine( ['advice','The unit cells should match to the lower resolution of the two datasets'])
    self.createLine([
      'tip','DANGEROUS: only sensible if the unit cells are very similar',
      'widget', 'OVERRIDE_CELL_DIFFERENCE',
      'label', 'allow existing freeR set to have different unit cells'])
    self.createLine(['advice',
          '<span style="color: DarkOrange;font-weight: bold;">Be sure you know what you are doing: the cells must be very similar even if outside the test limits</span>'])
    self.closeSubFrame()

    # We almost certainly want to keep whatever was set when the data was saved so dont call updateFromFile
    #self.updateFromFile(force=False)
    # Beware connecting to any dataChanged from HKLIN got a signal at run time when dbFileId set on the HKLIN
    self.container.inputData.HKLIN.dataChanged.connect(self.updateFromFile)
    # old combination of spacegroup and cell, now separate
    #self.getWidget('SPACEGROUPCELL').validate()
    ###self.updateFromFile()

  # -------------------------------------------------------------
  def getMaximumResolution( self ) :
    isSet = False
    if not self.container.inputData.CAN_CUT_RESOLUTION:
      return isSet
    highRes = float(self.container.inputData.MAXIMUM_RESOLUTION.get())
    print("getMaximumResolution", highRes)
    if highRes > 0.0:
      label = "%5.2f\xc5" % highRes
      isSet = True
    else:
      label = "Unknown"
    if self.maximum_resolution_label: self.maximum_resolution_label.setText(label)
    return isSet
  #  -------------------------------------------------------------
  def resorangeset(self):
    # Is the resolution cutoff range set
    self.container.inputData.RESOLUTION_RANGE_SET.set(False)
    if self.container.inputData.RESOLUTION_RANGE:
      r1 = self.container.inputData.RESOLUTION_RANGE.start
      r2 = self.container.inputData.RESOLUTION_RANGE.end
      if r1.isSet() or r2.isSet():
        self.container.inputData.RESOLUTION_RANGE_SET.set(True)

# -------------------------------------------------------------
  def fix(self):
    # disconnect updateFromFile() cos that and unSetAll() were being called by processing of the  HKLIN file (ef setDbFileId()) (maybe not?)

    print("\n*** CTaskimport_merge, fix")
    self.container.inputData.HKLIN.dataChanged.disconnect(self.updateFromFile)
    if self.container.guiParameters.HKLINISMTZ:
      for item in [ 'WAVELENGTH', 'CRYSTALNAME', 'DATASETNAME' ]:
        self.container.inputData.get(item).setQualifiers( { 'allowUndefined' : True } )
        self.getWidget(item).validate()
    return CErrorReport()
    
  # -------------------------------------------------------------
  @QtCore.Slot(bool)
  def updateFromFile(self,force=True):
    # Explicit call to CGenericREflnDataFile.getFileContent() otherwise CData properties code gets it wrong
    self.unSetAll()
    fc = self.container.inputData.HKLIN.getFileContent()
    print("\n** updateFromFile HKLIN ",self.container.inputData.HKLIN.fileContent)

    if not self.container.inputData.HKLIN.isSet():
      print("HKLIN unset")
      return

    #print("IDRR", self.container.inputData.RESOLUTION_RANGE)
    #self.unSetAll()

    self.container.inputData.HKLIN.loadFile()
    # What is the format?
    fformat = self.container.inputData.HKLIN.getFormat()
    print("Input file format", fformat)
    self.container.guiParameters.HKLINISMTZ = (fformat == 'mtz')  # MTZ
    if fformat == 'mtz':
        self.container.inputData.HKLIN_FORMAT.set('MTZ')
        self.container.inputData.CAN_CUT_RESOLUTION.set(True)
    elif fformat == 'mmcif':
        self.container.inputData.HKLIN_FORMAT.set('MMCIF')
        self.container.inputData.CAN_CUT_RESOLUTION.set(True)
    elif fformat == 'sca':
        self.container.inputData.HKLIN_FORMAT.set('SCA')
        self.container.inputData.CAN_CUT_RESOLUTION.set(True)
    else:
        self.container.inputData.HKLIN_FORMAT.set('OTHER')

    fileContent = self.container.inputData.HKLIN.getFileContent()
    self.container.controlParameters.STARANISO_DATA.set(False)

    message = ''
    if self.container.inputData.HKLIN.getMerged() == 'unmerged':
        if fformat == 'shelx':
            message = 'Shelx files may contain merged or unmerged data - you should  use the Data Reduction task to import it'
        else:
            message = 'This apppears to be an unmerged file - please use the Data Reduction task to import it'
        msgbox = MyMessageBox()
        msgbox.setDefault(1)
        msgbox.setInformativeText(message)
        ret = msgbox.displayText('Importing merged file')
        if ret == 1:
            # Abort button selected, bail out
            self.container.inputData.HKLIN.unSet()
            self.container.inputData.UNITCELL.unSet()
            self.container.inputData.SPACEGROUP.unSet()
            return
        # Open button, continue to try
       
    if self.container.inputData.HKLIN_FORMAT == 'MMCIF':
      if not self.openMmcifFile():
        return   # fail
      self.processMmcifFile()
      return
    else:
      # not mmcif
      self.container.guiParameters.SHOW_MMCIF_BLOCKS.set(False)
      # Not mmcif
      self.container.inputData.UNITCELL.set(fileContent.cell.fix(fileContent.cell.get()))
      if force or fileContent.spaceGroup.isSet():
        self.container.inputData.SPACEGROUP.set(fileContent.spaceGroup)
      #print('***CTaskimport_merged.updateFromFile wavelengths ',fileContent.contents('wavelengths'))
      if fileContent.contents('wavelengths') is not None and fileContent.wavelengths.isSet() and len(fileContent.wavelengths)>0:
        # extract relevant wavelength
        for wavelength, dataset in zip(fileContent.wavelengths, fileContent.datasets):
          if dataset != 'HKL_base':
            break         
        #print("wavelength, dataset", wavelength, dataset)
        self.container.inputData.WAVELENGTH.set(self.container.inputData.WAVELENGTH.fix(wavelength.__float__()))
      elif fileContent.contents('wavelength') is not None and fileContent.wavelength.isSet():
        self.container.inputData.WAVELENGTH.set(self.container.inputData.WAVELENGTH.fix(fileContent.wavelength.__float__()))
      elif force:
        self.container.inputData.WAVELENGTH.unSet()

      if self.container.inputData.HKLIN.getFormat() == 'mtz':
        #   Resolution
        lowres = fileContent.resolutionRange.low
        highres = fileContent.resolutionRange.high
        print("MTZ Resolution range", lowres, highres)
        self.container.inputData.MAXIMUM_RESOLUTION.set(highres)
        self.container.controlParameters.SKIP_FREER.set(False)
        self.container.controlParameters.STARANISO_DATA.set(False)
        self.container.inputData.HASFREER.set(False)
        for column in fileContent.listOfColumns:
          if 'FREE' in str(column.columnLabel).upper():
          # if 'FREER' in str(column.columnLabel).upper():
            self.container.inputData.HASFREER.set(True)
          if str(column.columnLabel) == 'SA_flag':
            # Data comes from StarAniso
            #print('>>** found SA_FLAG')
            self.container.controlParameters.STARANISO_DATA.set(True)
            self.container.controlParameters.SKIP_FREER.set(True)
            self.container.controlParameters.COMPLETE.set(False)
            self.container.inputData.FREERFLAG.unSet()
            self.container.inputData.HASFREER.unSet()

        # MTZ
        self.selectObsColumns()
        self.getMaximumResolution()
      else:
        # not MTZ (and not mmcif)
        #  Assume no freer set
        self.container.inputData.HASFREER.set(False)
        self.getMaximumResolution()

        try:
          if fileContent.datasetName.isSet() and len(fileContent.datasetName)>=1:
            self.container.inputData.DATASETNAME.set(str(fileContent.datasetName))
          elif force:
            self.container.inputData.DATASETNAME.unSet()
        except:
          pass
        try:
          if fileContent.crystalName.isSet() and len(fileContent.crystalName)>=1:
            self.container.inputData.CRYSTALNAME.set(str(fileContent.crystalName))
          elif force:
            self.container.inputData.CRYSTALNAME.unSet()
        except:
          pass
     
  # -------------------------------------------------------------
  def selectObsColumns(self):
    print("selectObsColumns")
    hklin = self.container.inputData.HKLIN
    mtzModel = self.container.inputData.get('HKLIN_OBS')
    mtzModel.setFullPath(str(hklin))
    mtzModel.loadFile()
    
    errors = mtzModel.validColumns()
   
    if errors.maxSeverity()>SEVERITY_WARNING:
      #print errors.report()
      if errors.count(cls=mtzModel.__class__,code=203)>0:
        QtWidgets.QMessageBox.warning(self,'Error in selected MTZ file','Selected MTZ file does not contain correct type of data')
        hklin.unSet()
        return
      elif errors.count(cls=mtzModel.__class__,code=204)>0 or errors.count(cls=mtzModel.__class__,code=205)>0:
        applyNow = (errors.count(cls=mtzModel.__class__,code=205)>0)
        self.selColDialog=CCP4XtalWidgets.CSelectColumnsWidget(parent=self,model=mtzModel,applyNow=applyNow,filename=str(hklin))
        self.selColDialog.applySignal.connect(self.handleSelColDialogApply)
        self.selColDialog.cancelSignal.connect(self.handleSelColDialogCancel)
      else:
        pass

  # -------------------------------------------------------------
  @QtCore.Slot()
  def handleSelColDialogApply(self):
    contentFlag,i2Names,dataset,colLabels = self.selColDialog.getSelection()
    #print('handleSelColDialogApply getSelection',contentFlag,i2Names,dataset,colLabels)
    try:
      self.selColDialog.hide()
    except:
      pass
    jobId = self.jobId()
    projectId = self.projectId()
    mtzModel =  self.container.inputData.HKLIN_OBS
    mtzModel.setFullPath(str(self.container.inputData.HKLIN))
    
    # PRE Feb2022
    # This section seems to just be here to pick some errors - not sure it is needed
    error = mtzModel.splitMtz(jobId=jobId,projectId=projectId,contentFlag=None,i2Labels=i2Names,columnLabels=colLabels)
    sourceFileName = mtzModel.__dict__.get('sourceFileName','')
    if sourceFileName is None: sourceFileName = ''
    #print 'handleSelColDialogApply',error.report(),sourceFileName
    if error.maxSeverity()==SEVERITY_WARNING and error[0]['code']==212:
      mess = QtWidgets.QMessageBox.warning(self,self.windowTitle(),'This data is already imported as\n'+error[0]['details'])
    elif error.maxSeverity()>=SEVERITY_WARNING:
      if error[0]['code']==211:
        mess = QtWidgets.QMessageBox.warning(self,self.windowTitle(),'No column data selected')
      else:
        error.warningMessage(windowTitle='Splitting MTZ: '+sourceFileName,jobId=jobId,parent=self)
      self.container.inputData.HKLIN.unSet()
      
    else:
      mtzModel.loadFile()
      columns = ''
      ok = True
      # selColDialog has returned the things we need
      #  contentFlag, dataset, column names
      for colName in colLabels:
          if colName is None or colName in ['']:
            ok = False
          else:
            columns = columns+str(colName)+','
      
      self.selColDialog.deleteLater()
      if not ok:
        self.unSetAll(ifHklin=True)
        print("Not OK")
      else:
        # MTZ format
        self.container.inputData.HKLIN_OBS_CONTENT_FLAG = contentFlag
        self.container.inputData.HKLIN_OBS_COLUMNS = columns[0:-1]
        self.container.guiParameters.HKLIN_HAS_COLUMNS.set(True)
        colLabel = ' MTZ columns: ' + str(self.container.inputData.HKLIN_OBS_COLUMNS)
        self.selectedColumnLabels.setText(colLabel)
        self.selectedColumnAdvice.setText(self.columnAdviceMTZ)
        #print 'handleSelColDialogApply',mtzModel.fileContent.datasets,mtzModel.fileContent.crystalNames
        # Expect the datasets/cryalNames lists to be HKL_base and possibly the crystal/dataset for the extracted data
        datasets =  mtzModel.fileContent.datasets
        idts = -1
        for idx in range(len(datasets)):
          dtnm = str(datasets[idx])
          if dataset in dtnm:
            idts = idx
            break
        if idts < 0:
          print("Dataset", dataset," not found")
          self.unSetAll(ifHklin=True)
        else:
          self.container.inputData.DATASETNAME.set(str(mtzModel.fileContent.datasets[idts]))
        self.container.inputData.CRYSTALNAME.set(str(mtzModel.fileContent.crystalNames[idts]))
        self.container.inputData.UNITCELL.set(mtzModel.fileContent.datasetCells[idts])
        self.getWidget('DATASETNAME').validate()
        self.getWidget('CRYSTALNAME').validate()

  # -------------------------------------------------------------
  @QtCore.Slot()
  def handleSelColDialogCancel(self):
    self.selColDialog.hide()
    self.selColDialog.deleteLater()
    self.unSetAll(ifHklin=True)

  # -------------------------------------------------------------
  @QtCore.Slot()
  def handleSelectFreeRflag(self):

    if self.container.inputData.FREERFLAG.isSet():
      self.container.controlParameters.COMPLETE.set(True)
      self.container.inputData.HASFREER.set(True)
    else:
      self.container.controlParameters.COMPLETE.set(False)

  # -------------------------------------------------------------
  def toggleCutResolution(self):
    #print('>> toggleCutResolution', self.container.inputData.FREERFLAG.isSet())
    if self.container.inputData.FREERFLAG.isSet():
      return True
    return False
  # -------------------------------------------------------------
  def toggleFraction(self):
    #print('>> toggleFraction FRF', self.container.inputData.FREERFLAG.isSet(),
    #      "HF", self.container.inputData.HASFREER)
    if self.container.inputData.FREERFLAG.isSet():
      return False
    if self.container.inputData.HASFREER:
      return False
    return True
    
  # -------------------------------------------------------------
  def toggleSkip2(self):
    #  Leave FreeR unchanged
    # print('>> toggleSkip2 FRF', self.container.inputData.FREERFLAG.isSet(),
    #       "HF", self.container.inputData.HASFREER)
    if self.container.inputData.FREERFLAG.isSet():
      return False
    if not self.container.inputData.HASFREER:
      return False
    if not self.container.controlParameters.STARANISO_DATA:
      return False    
    return True
  
  # -------------------------------------------------------------
  def toggleSkip(self):
    #  Leave FreeR unchanged
    #print('>> toggleSkip FRF', self.container.inputData.FREERFLAG.isSet(),
    #      "HF", self.container.inputData.HASFREER)
    if self.container.inputData.FREERFLAG.isSet():
      return False
    if not self.container.inputData.HASFREER:
      return False    
    if self.container.controlParameters.STARANISO_DATA:
      return False    
    return True
    
  # -------------------------------------------------------------
  def unSetAll(self,ifHklin=False):
    #print("unSetAll", ifHklin)
    for item in ['HKLIN_OBS','HKLIN_OBS_COLUMNS','HKLIN_OBS_CONTENT_FLAG',
                 'CRYSTALNAME','DATASETNAME',
                 'MAXIMUM_RESOLUTION','CAN_CUT_RESOLUTION',
                 'MMCIF_SELECTED_BLOCK', 'MMCIF_SELECTED_DETAILS',
                 'MMCIF_SELECTED_INFO', 'MMCIF_SELECTED_COLUMNS',
                 'MMCIF_SELECTED_CONTENT',
                 'MMCIF_SELECTED_ISINTENSITY']:
        self.container.inputData.get(item).unSet()       
    for item in ['HKLIN_HAS_COLUMNS', 'SHOW_MMCIF_BLOCKS', 'HKLINISMTZ']:
      self.container.guiParameters.get(item).unSet()
    for item in ['WAVELENGTH','CRYSTALNAME','DATASETNAME']:
      self.getWidget(item).validate()
    self.container.inputData.RESOLUTION_RANGE_SET.set(False)
    if ifHklin:
      self.container.inputData.HKLIN.blockSignals(True)
      self.container.inputData.HKLIN.unSet()
      self.container.inputData.HKLIN.blockSignals(False)

  # -------------------------------------------------------------
  def openMmcifFile(self):
      print("openMmcifFile", self.container.inputData.HKLIN)
      if not self.container.inputData.HKLIN.isSet():
        # No hklin
        self.unSetAll(False)
        return
      else:
        try:
          self.mmcif = gemmi.cif.read(str(self.container.inputData.HKLIN))
        except Exception as e:
          exceptionMessage = str(e)
          print("Read fail\n", exceptionMessage)
          if 'duplicate' in exceptionMessage:
            message = exceptionMessage[exceptionMessage.find('duplicate'):]
          else:
            message = exceptionMessage
          msgbox = MyMessageBox()
          msgbox.setInformativeText(message)
          msgbox.singleButton()
          ret = msgbox.displayText('Cannot import merged file, illegal format')
          self.container.inputData.HKLIN.unSet()
          self.container.inputData.UNITCELL.unSet()
          self.container.inputData.SPACEGROUP.unSet()
          return False
      
        self.rblocks = gemmi.as_refln_blocks(self.mmcif)
        if len(self.rblocks) == 0 or not self.rblocks[0]:
          # not a reflection mmcif file
          mess = 'This does not seem to be a reflection mmcif file: may be coordinates?'
          msgbox.setInformativeText(message)
          msgbox.singleButton()
          ret = msgbox.displayText('Importing merged file')
          return False

        self.cifblockinfo = []
        for rb in self.rblocks:
          self.cifblockinfo.append(CifBlockInfo(rb))   # in mmcifutils
          cbi = self.cifblockinfo[-1]
          printBlockInfo(cbi)

        # check if file is from StarAniso
        self.checkForStarAniso()
        return True

  # -------------------------------------------------------------
  def processMmcifFile(self):
    print("** processMmcifFile", self.container.inputData.HKLIN)
    self.container.guiParameters.SHOW_MMCIF_BLOCKS.set(True) # for mmCIF file

    nblock = len(self.rblocks)

    # for each accepted reflection block
    #  only list "accepted" blocks which have I or F data
    ids = []          #  block names
    detailslist = []  #  _diffrn.details
    infolist = []     #  info, hkl list type  
    columnlist = []   #  column labels
    #  information about blocks not accepted for input
    #    eg unmerged data, or phases etc, not Is or Fs
    otherlist = []
    idxblkinfo = []
    for idx, cifinfo in enumerate(self.cifblockinfo):
      #print(">>>")
      #printBlockInfo(cifinfo)
      if cifinfo.merged_diffn_data():
        # block contains Is or Fs
        idxblkinfo.append(idx)
        ids.append(cifinfo.bname)
        details = cifinfo.details   # _diffrn.details if present
        if details is None: details = ''
        detailslist.append(details)
        info = cifinfo.info   # formatted column info
        info += "\n    hkl list type: "+cifinfo.hklcheckformat
        # Highest resolution
        info += "\n    Highest resolution: "+ "{:.2f}A".format(cifinfo.highres)
        # FreeR status, list of any warning messages 
        freerwarning = cifinfo.freerWarning()  # a list of warnings
        if freerwarning is not None:
          for w in freerwarning:
            info += "\n   " + w
        if info is None: info = ''
        infolist.append(info)
        columns = cifinfo.columnnames
        if columns is None: columns = ''
        columnlist.append(columns)
      else:
        otherlist.append(self.otherinfo(cifinfo))

    #print("detailslist", len(detailslist), detailslist)
    #print("infolist", infolist)
    #print("columnlist", columnlist)
    self.container.guiParameters.MMCIF_INDICES.set(idxblkinfo)
    self.container.guiParameters.MMCIF_BLOCKNAMES.set(ids)
    self.container.guiParameters.MMCIF_BLOCK_DETAILS.set(detailslist)
    self.container.guiParameters.MMCIF_BLOCK_INFO.set(infolist)
    self.setColumnNames(columnlist)
    # non-accepted blocks
    self.container.guiParameters.MMCIF_BLOCK_OTHER.set(otherlist)
    
    mbnl = list(self.container.guiParameters.MMCIF_BLOCKNAMES)
    for b in self.container.guiParameters.MMCIF_BLOCKNAMES:
        s = str(b)
    self.extractMmcifInfo()
    self.setCIFblocklist()

  # -------------------------------------------------------------
  def otherinfo(self, cifinfo):
    # Information about an mmcif not suitable for merged input
    s = ">>> mmCIF block " + cifinfo.bname + ": "
    details = cifinfo.details   # _diffrn.details if present
    if not cifinfo.ismerged():
      s += "Unmerged data"
      if details is not None:
        s += ": "+details
    else:
      info = cifinfo.info   # formatted column info
      info += "\n   hkl list type: "+cifinfo.hklcheckformat
      s += info

      if details is not None:
        s += ": "+details

      #  cifinfo.columnnames is a dictionary
      columnlist = self.formatColumnNames(cifinfo.columnnames)
      s += '\n   '+columnlist

    return s  

  # -------------------------------------------------------------
  def setColumnNames(self, columnlist):
      # columnlist is a dictionary of columns names indexed by content
      columnlists = []
      for clist in columnlist:
          collist = None
          for ctype in ReflectionDataTypes.DATA_PRIORITY:
              if ctype in clist:
                  collist = clist[ctype]
                  fclmns = self.formatColumnlist(collist)
                  columnlists.append(fclmns)
                  break
      self.container.guiParameters.MMCIF_BLOCK_COLUMNS.set(columnlists)

  # -------------------------------------------------------------
  def formatColumnNames(self, columnlist):
      # columnlist is a dictionary of columns names indexed by content
      # Return list of column names
      
      # Possible content types, in order of priority
      s = ''
      for key in columnlist:
          collist = columnlist[key]
          fclmns = self.formatColumnlist(collist)
          s += str(key) + ': [' + fclmns + '], '
      return s[:-2]

  # -------------------------------------------------------------
  def formatColumnlist(self, collist):
      # collist is a list of strings, return string
      s = ''
      if len(collist) == 0: return s
      for label in collist:
          s += label + ", "
      return s[:-2]
  # -------------------------------------------------------------
  def extractMmcifInfo(self, blockname=None):
      print("**extractMmcifInfo", blockname)

      '''
      We want to set:
      self.container.inputData.UNITCELL
      self.container.inputData.SPACEGROUP
      self.container.inputData.WAVELENGTH.
      self.container.inputData.DATASETNAME
      self.container.inputData.CRYSTALNAME
      self.container.inputData.MAXIMUM_RESOLUTION
      '''
      
      nblocks = len(self.rblocks)
      naccepted = len(self.container.guiParameters.MMCIF_BLOCK_DETAILS)
      if nblocks > 0 and naccepted == 0:
        self.processMmcifFile()
        naccepted = len(self.container.guiParameters.MMCIF_BLOCK_DETAILS)

      cifinfo = None
      jblock = 0
      if blockname is None or nblocks == 1:
          # take the first one
          rblock = self.rblocks[0]
          cifinfo = self.cifblockinfo[0]
          blockname = cifinfo.bname
      else:
          for i in range(naccepted):
            j = self.container.guiParameters.MMCIF_INDICES[i]
            if self.cifblockinfo[i].bname == blockname:
              cifinfo = self.cifblockinfo[j]
              jblock = i  # selected block
              break

      if cifinfo is None:
          print("***Failed") # shouldn't happen111111111111111111111111111

      #print("*Rcell:",cifinfo.cell)
      #print("RSG:", cifinfo.spacegroup_name)
      #print("*Rwvl:",cifinfo.wavelength)

      rc = cifinfo.cell
      self.container.inputData.UNITCELL.set(rc)

      self.container.inputData.SPACEGROUP.set(cifinfo.spacegroup_name)
      self.container.inputData.WAVELENGTH.set(cifinfo.wavelength)

      self.container.inputData.DATASETNAME.set(blockname)
      self.container.inputData.CRYSTALNAME.set(blockname)
      self.container.inputData.MMCIF_SELECTED_BLOCK.set(blockname)
      self.container.inputData.MMCIF_SELECTED_DETAILS.set(\
          self.container.guiParameters.MMCIF_BLOCK_DETAILS[jblock])
      self.container.inputData.MMCIF_SELECTED_INFO.set(\
          self.container.guiParameters.MMCIF_BLOCK_INFO[jblock])
      self.container.inputData.MMCIF_SELECTED_COLUMNS.set(\
          self.container.guiParameters.MMCIF_BLOCK_COLUMNS[jblock])
      self.container.inputData.MAXIMUM_RESOLUTION.set(\
        cifinfo.highres)

      #print("M_S_I", self.container.inputData.MMCIF_SELECTED_INFO)
      #print("M_S_C", self.container.inputData.MMCIF_SELECTED_COLUMNS)
      self.getMaximumResolution()
      #  +1 if intensity, -1 if amplitude, 0 if unknown
      isintensity = -1
      if 'I' in str(self.container.inputData.MMCIF_SELECTED_INFO):
          # only intensity contents contains the letter 'I'
          isintensity = +1
      else:
          isintensity = -1

      self.container.inputData.MMCIF_SELECTED_ISINTENSITY.set(isintensity)

      # Does block have FreeR data?
      frvalid = cifinfo.validFreeR()
      # Returns: None if no FreeR; True if valid; False if invalid (eg all the same)
      freervalid = frvalid
      if frvalid is None: freervalid = False  # invalid and missing are treated the same
      self.container.inputData.HASFREER.set(freervalid)
      if freervalid:
        if self.container.controlParameters.STARANISO_DATA:
          self.container.controlParameters.SKIP_FREER.set(True)
      if not freervalid:
        self.container.controlParameters.SKIP_FREER.set(False)
        #print("\n*** SKIP_FREER set to False")

      self.getWidget('DATASETNAME').validate()
      self.getWidget('CRYSTALNAME').validate()

      merged = cifinfo.ismerged()
      if not merged:
          mess = QtWidgets.QMessageBox.warning(self,'Importing merged file','This apppears to be an unmerged data block - please use the Data Reduction task to import it')

      # Check reflection label sets
      if not cifinfo.OK:
          # Probably not a reflection block, eg maybe coordinates
          mess = QtWidgets.QMessageBox.warning(self,'Importing merged file',
     'This does not seem to be a reflection mmcif file: may be coordinates?')
          return
      
      message = ""
      if nblocks == 1:
          message = "This file"
      else:
          message = "This data block"
          
      status, mess = cifinfo.columnsOK()
      message += mess

      if status < 0:
          mess = QtWidgets.QMessageBox.warning(self,
                                               self.windowTitle(), message)
      
  # -------------------------------------------------------------
  def setMMCIFframe(self):
      # just make an empty area for later, if needed
      # and if HKLIN is MMCIF, create file data objects
      print("markCIFframe", self.container.inputData.HKLIN)
      if not self.container.inputData.HKLIN:
        return
      if self.widget.subFrame is not None:
          #print("subFrame")
          self.cifpane = self.widget.subFrame.layout()
      else:
          #print("no subFrame")
          self.cifpane = self.widget.currentFolderLayout

      self.cifbuttons = ChoiceButtons()
      self.cifpane.addWidget(self.cifbuttons)
      self.selectedBlock = ''
      self.cifbuttons.clickedSignal.connect(self.cifblockClicked)

      if self.container.inputData.HKLIN_FORMAT == 'MMCIF':
          # make self.mmcif & .rblocks
          self.openMmcifFile()
          
  # -------------------------------------------------------------
  def infoList(self):
    # Get list of column content types from selected Cif block
    contents = str(self.container.inputData.MMCIF_SELECTED_INFO).split(',')
    # Remove leading spaces, if present
    for idx, ct in enumerate(contents):
      contents[idx] = ct.strip()

    return contents
    
  # -------------------------------------------------------------
  def cifColumnSelect(self):
    # column selection box, where appropriate
    #  returns type found or selected
    if self.container.inputData.MMCIF_SELECTED_CONTENT.isSet():
      return
    
    msc = self.container.inputData.MMCIF_SELECTED_CONTENT
    if msc == 'NONE':
      msc = "None"
    if msc == '':
      msc = 'Blank'
    #print(">*> cifColumnSelect", msc)
    #print("<< ", self.container.inputData.MMCIF_SELECTED_CONTENT)
    #print("<<** ", type(self.container.inputData.MMCIF_SELECTED_CONTENT))
    #print("<*** ", len(self.container.inputData.MMCIF_SELECTED_CONTENT))
   
    self.selectedContent = None

    contents = self.infoList()
    #print("contents", contents)
    #print("Ctypes", ReflectionDataTypes.DATA_PRIORITY)
    ntypefound = 0
    ctypefound = None
    choices = []
    for ctype in ReflectionDataTypes.DATA_PRIORITY:
      if ctype in contents:
        # Possible type found, in order of priority
        ntypefound += 1
        ctypefound = ctype
        choices.append(ctype)

    if ntypefound == 0:
      return
    
    # If only one found, just accept it
    if ntypefound == 1:
      self.selectedContent = ctypefound
      self.setMMCIFcontentColumns()
      block = str(self.container.inputData.MMCIF_SELECTED_BLOCK)
      clabel =  ' MMCIF Block: ' + block + ', content type: ' + self.selectedContent
      self.selectedColumnLabels.setText(clabel)
      self.selectedColumnAdvice.setText(' sole available option')
      return
  
    title = "Select data column group to import"
    subtitle = "shown in order of priority"
    tags = ['  <<    default', '', '', '']

    notes = []
    for ctype in choices:
      l = ReflectionDataTypes.TYPE_LABELS[ctype]
      notes.append(l)
    
    selbox = selectBox(title, choices, tags, notes, subtitle=subtitle)
    selbox.execit()

    self.selectedContent = selbox.getSelected()
    #print(">Column selected", self.selectedContent)
    if self.selectedContent is None:
      self.selectedContent = "Unspecified"
    if self.selectedContent == 'Cancelled':
      self.unSetAll(ifHklin=True)
      self.selectedContent = "Unspecified"
      
    self.setMMCIFcontentColumns()
    block = str(self.container.inputData.MMCIF_SELECTED_BLOCK)
    clabel =  ' MMCIF Block: ' + block + ', content type: ' + self.selectedContent
    self.selectedColumnLabels.setText(clabel)
    self.selectedColumnAdvice.setText(self.columnAdviceMMCIF)

    selbox.deleteLater()
  # -------------------------------------------------------------
  def setMMCIFcontentColumns(self):
    # store content and mmcif columns
    self.container.inputData.MMCIF_SELECTED_CONTENT.set(self.selectedContent)
    if self.selectedContent == "Unspecified":
      return

    cifspecs = ReflectionDataTypes.REFLECTION_DATA[self.selectedContent]
    #print("cifspecs", cifspecs)
    ciflabels = ''
    for spec in cifspecs:
      l = spec.split()
      ciflabels += l[0] + ', '
    
    self.container.inputData.MMCIF_SELECTED_COLUMNS.set(ciflabels[:-2])
    
  # -------------------------------------------------------------
  @QtCore.Slot(str)
  def cifblockClicked(self):
      s = self.cifbuttons.selected
      #print("Clicked", s)
      self.selectedBlock = s
      self.extractMmcifInfo(s)
      self.container.inputData.MMCIF_SELECTED_CONTENT.unSet()
      self.cifColumnSelect()

  # -------------------------------------------------------------
  def setCIFblocklist(self):

      infolist = [''] # One entry for each block
      if self.container.guiParameters.MMCIF_BLOCK_INFO:
          infolist = self.strlist(self.container.guiParameters.MMCIF_BLOCK_INFO)
          # Edit infolist
           #   (this would be easier if I hadn't lumped different things together)
          for i, info in enumerate(infolist):
            info = 'Column content type: '+info
            if not 'FreeR' in info:
              idx = info.find('\n')
              info = info[:idx]+"   NB no FreeR flag"+info[idx:]
            infolist[i] = info
          
          #print("\n>>* infolist ", infolist)
            
      #print("Info", self.container.guiParameters.MMCIF_BLOCK_INFO)

      nblocks = len(self.container.guiParameters.MMCIF_BLOCKNAMES)
      if nblocks > 1:
          title = 'Select which mmCIF reflection block to use (first is default)'
      else:
          title = 'mmCIF file contains one appropriate reflection block'

      self.cifbuttons.setChoices(title,
            self.strlist(self.container.guiParameters.MMCIF_BLOCKNAMES),
            self.strlist(self.container.guiParameters.MMCIF_BLOCK_DETAILS),
                                 infolist)

      if self.container.guiParameters.MMCIF_BLOCK_OTHER.isSet():
        if len(self.container.guiParameters.MMCIF_BLOCK_OTHER) > 0:
          self.cifbuttons.addOtherText(\
          "Reflection blocks not containing importable merged reflection data",
            self.strlist(self.container.guiParameters.MMCIF_BLOCK_OTHER))

      self.container.guiParameters.HKLIN_HAS_COLUMNS.set(True)
          
      self.cifColumnSelect()

  # -------------------------------------------------------------
  def checkForStarAniso(self):
    # current check loop for _software.name == STARANISO in 1st block
    self.container.controlParameters.STARANISO_DATA.set(False)

    block = self.rblocks[0].block
    listofitems = list(block.find_loop('_software.name'))
    #print("checkForStarAniso item", listofitems);
    
    if listofitems is not None:
      if 'STARANISO' in listofitems:
        print(">> STARANISO data")
        self.container.controlParameters.STARANISO_DATA.set(True)

        hasfreer = self.container.inputData.HASFREER
        if hasfreer:
          self.container.controlParameters.SKIP_FREER.set(True)
        else:
          self.container.controlParameters.SKIP_FREER.set(False)
    
  # -------------------------------------------------------------
  def strlist(self, cstringlist):
      # string list from CString list
      sl = []
      for i in range(len(cstringlist)):
          sl.append(str(cstringlist[i]))
      return sl
