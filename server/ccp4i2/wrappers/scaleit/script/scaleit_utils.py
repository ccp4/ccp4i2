from ccp4i2.pipelines.aimless_pipe.script.aimless_pipe_utils import CellCheck, CellFormat


class DatalistCheck:
  # check that all items in the list are comaptible,
  # and there are at least two entries
  #   for interface and run time
  def __init__(self, obsdatalist):
      # obsdatalist is CMiniMtzDataFileList
      self.obslist = obsdatalist

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def checkAll(self):
    # non-interactive call
    validlist = self.compareAllCells()
    OK = True
    if False in validlist:
        OK = False
        cellsgformat = self.formatAllCellSGs()   # List
        self.failtext = ['Inconsistent cells']
        self.failtext += cellsgformat
        print("** failtext", self.failtext)

    return OK
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def formatFail(self):
    s = ''
    for line in self.failtext:
      s += line + "\n"
    return s[:-2]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def compareAllCells(self):
    nfiles = len(self.obslist)
    validlist = [True]*nfiles  # first one is always OK
    if len(self.obslist) <= 1:
      return validlist
    
    content0 = self.obslist[0].fileContent

    # Compare each file to the first one
    i = 1
    for obsfile in self.obslist[1:]:
      if obsfile.isSet():
        # Check cells
        valid = self.compareCells(content0, obsfile.fileContent)
        if not valid:
          print("File ", i, " is NOT compatible with 1st file")
          #print("obsfile", dir(obsfile))
          validlist[i] = False

      return validlist

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def compareCells(self, content0, content1):
    cellcheck = CellCheck(content0, content1)
    valid = cellcheck.isValid()
    #print("compareCells", valid)
    return valid

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  def formatAllCellSGs(self):
    s = []
    cellformat = CellFormat()
    i = 1
    for obsfile in self.obslist:
      cell = obsfile.fileContent.cell
      cf = cellformat.shortformatCell(cell)
      
      sgname = ''
      if obsfile.fileContent.spaceGroup.isSet():
        sgname = str(obsfile.fileContent.spaceGroup).replace(' ','')
        s.append('{:2}:{};{}'.format(i, sgname, cf))
      else:
        s.append('{:2}:{}'.format(i, cf))
      i += 1
    return s
