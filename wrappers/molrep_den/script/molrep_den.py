
"""
     molrep_mr.py: CCP4 GUI Project
     Copyright (C) 2011 STFC

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

import os
import re
import shutil
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
#from lxml import etree
from xml.etree import ElementTree as ET

class molrep_den(CPluginScript):

    TASKTITLE='Molecular replacement with electron density - MOLREP'
    TASKNAME = 'molrep_den'
    TASKCOMMAND = 'molrep'
    TASKMODULE = 'test'
    MAINTAINER = 'andrey.lebedev@stfc.ac.uk'

    ERROR_CODES = {101 : { 'severity': SEVERITY_WARNING, 'description' : 'Failed extracting tables from molrep.doc' },
                   102 : { 'severity': SEVERITY_WARNING, 'description' : 'Failed writing tables to program.xml' },
                   103 : { 'severity': SEVERITY_WARNING, 'description' : 'No tables extracted from molrep.doc' },
                   104 : { 'description' : 'No output coordinate file from Molrep' },
                   105 : { 'description' : 'No output log file from Molrep' },
                   106 : { 'description' : 'Error parsing log file from Molrep' }}


    def processInputFiles(self):
      # Ensure the obs data is in form of F_SIGF
      if self.container.guiParameters.PERFORM != 'den':
        from core import CCP4XtalData
        # Using CObsDataFile.convert() did not work for input of anomalous Is (as from aimless)
        self.F_SIGF_hklin,errReport = self.makeHklin([['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]])
        #self.F_SIGF_hklin,errReport = self.container.inputData.F_SIGF.convert(targetContent=CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN)
        #print 'molrep_mr.processInputFiles',self.F_SIGF_hklin,errReport
        if self.F_SIGF_hklin is None:
            return CPluginScript.FAILED
      if self.container.guiParameters.PERFORM != 'srf':
        if self.container.inputData.XYZIN.isSelectionSet():
          self.selectedXYZIN = os.path.join(self.workDirectory,'XYZIN_selected_atoms.pdb')
          rv = self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.selectedXYZIN)
        elif self.container.inputData.XYZIN.getExt() != '.pdb':
          self.selectedXYZIN = os.path.join(self.workDirectory,'XYZIN_pdb.pdb')
          self.container.inputData.XYZIN.convertFormat('pdb',self.selectedXYZIN )
        else:
          self.selectedXYZIN = str( self.container.inputData.XYZIN)
      return CPluginScript.SUCCEEDED


    @staticmethod
    def replace_MSE(model, pdbin):
      with open(model) as istream:
        content = istream.read()
      inp = '(:?HETATM|ATOM  ).{11}MSE .*'
      out = lambda x: x.group(0)\
        .replace('MSE', 'MET')\
        .replace('SE', ' S')\
        .replace('S   MET', 'SD  MET')\
        .replace('HETATM', 'ATOM  ')
      content = re.sub(inp, out, content)
      with open(pdbin, 'w') as ostream:
        ostream.write(content)

    def makeCommandAndScript(self):
      inp = self.container.inputData
      par = self.container.controlParameters
      gui = self.container.guiParameters
      from core import CCP4Utils
      self.path_wrk = str( self.getWorkDirectory() )
      self.path_scr = os.path.join( self.path_wrk, 'scratch' )
      if not os.path.exists(self.path_scr): os.mkdir( self.path_scr )
      #print "self.path_wrk", self.path_wrk
      #print "self.path_scr", self.path_scr
      if str( gui.PERFORM ) != 'srf' :
         pdbin = os.path.join(self.path_wrk, 'molrep_input.pdb')
         self.replace_MSE(self.selectedXYZIN, pdbin)
         self.appendCommandLine( [ '-m', pdbin ] )
      if inp.ASUIN.isSet() and par.SEQ.__str__() != 'n':
        seqin = os.path.join(self.workDirectory,'SEQIN.fasta')
        inp.ASUIN.writeFasta(seqin)
        self.appendCommandLine( [ '-s', seqin ] )

      if str( gui.PERFORM ) == 'pat' :
         if  inp.XYZIN_FIX.isSet() :
            self.appendCommandLine( [ '-mx', str( inp.XYZIN_FIX.fullPath ) ] )

      if str( gui.PERFORM ) == 'den' :
         self.appendCommandLine( [ '-mx', str( inp.XYZIN_FIX.fullPath ) ] )

      self.appendCommandLine( [ '-ps', str( self.path_scr ) + '/' ] )
      self.appendCommandLine( [ '-po', str( self.path_wrk ) + '/' ] )
      self.appendCommandLine( [ '-i' ] )

      if str( gui.PERFORM ) == 'den' :
         self.appendCommandLine( [ '-f', str( inp.F_PHI_MAP.fullPath ) ] )
         self.appendCommandScript( "labin F=F PH=PHI" )
         self.appendCommandScript( "diff m" )
      else :
         #self.appendCommandLine( [ '-f', str( inp.F_SIGF.fullPath ) ] )
         self.appendCommandLine( [ '-f', self.F_SIGF_hklin ] )
         self.appendCommandScript( "labin F=F SIGF=SIGF" )

      if str( gui.PERFORM ) == 'den' :
         if str( par.PRF ) == 'n' :
            pass
#-          self.appendCommandScript( "prf n" )
         elif str( par.PRF ) == 'y' :
            self.appendCommandScript( "prf y" )
         elif str( par.PRF ) == 's' :
            self.appendCommandScript( "prf s" )
         else :
            raise Exception()

      if str( par.NP ) != 'Auto' :
         self.appendCommandScript( "np %s" %( str( par.NP ) ) )

      if str( par.NMON ) != 'Auto' :
         self.appendCommandScript( "nmon %s" %( str( par.NMON ) ) )

      if par.SCORE.isSet() :
         if str( par.SCORE ) == 'y' :
            pass
#-          self.appendCommandScript( "score y" )
         elif str( par.SCORE ) == 'n' :
            self.appendCommandScript( "score n" )
         elif str( par.SCORE ) == 'c' :
            self.appendCommandScript( "score c" )
         else :
            raise Exception()

      if par.ANISO.isSet() :
         if str( par.ANISO ) == 'y' :
            pass
#-          self.appendCommandScript( "aniso y" )
         elif str( par.ANISO ) == 'n' :
            self.appendCommandScript( "aniso n" )
         elif str( par.ANISO ) == 'k' :
            self.appendCommandScript( "aniso k" )
         else :
            raise Exception()

#-    if par.HIGH_PATH_VAR.isSet() :
#-       if str( gui.HIGH_PATH_VAR ) == 's' :
#-          pass
##          self.appendCommandScript( "high_pass_var s" )
#-       elif str( gui.HIGH_PATH_VAR ) == 'i' :
#-          self.appendCommandScript( "high_pass_var i" )
#-       elif str( gui.HIGH_PATH_VAR ) == 'r' :
#-          self.appendCommandScript( "high_pass_var r" )
#-       elif str( gui.HIGH_PATH_VAR ) == 'b' :
#-          self.appendCommandScript( "high_pass_var b" )
#-       else :
#-          raise Exception()

#-    if par.LOW_PATH_VAR.isSet() :
#-       if str( gui.LOW_PATH_VAR ) == 'c' :
#-          pass
##          self.appendCommandScript( "low_pass_var c" )
#-       elif str( gui.LOW_PATH_VAR ) == 'r' :
#-          self.appendCommandScript( "low_pass_var r" )
#-       elif str( gui.LOW_PATH_VAR ) == 'b' :
#-          self.appendCommandScript( "low_pass_var b" )
#-       else :
#-          raise Exception()

      if par.SEQ.isSet() :
         if str( par.SEQ ) == 'y' :
            pass
#-          self.appendCommandScript( "seq y" )
         elif str( par.SEQ ) == 'd' :
            self.appendCommandScript( "seq d" )
         elif str( par.SEQ ) == 'n' :
            self.appendCommandScript( "seq n" )
         else :
            raise Exception()

      if par.SURF.isSet() :
         if str( par.SURF ) == 'y' :
            pass
#-          self.appendCommandScript( "surf y" )
         elif str( par.SURF ) == 'c' :
            self.appendCommandScript( "surf c" )
         elif str( par.SURF ) == 'n' :
            self.appendCommandScript( "surf n" )
         elif str( par.SURF ) == '2' :
            self.appendCommandScript( "surf 2" )
         elif str( par.SURF ) == 'a' :
            self.appendCommandScript( "surf a" )
         else :
            raise Exception()

      if par.NMON_EXP.isSet() :
         if str( par.NMON_EXP ) != 'Auto' :
            self.appendCommandScript( "nmon_exp %s" %( str( par.NMON_EXP ) ) )

      if par.RESMIN.isSet(): self.appendCommandScript( "resmin %s" %( str( par.RESMIN ) ) )
      if par.RESMAX.isSet(): self.appendCommandScript( "resmax %s" %( str( par.RESMAX ) ) )

      if par.SG_OPTIONS == 'specify':
        if par.SG.isSet():  self.appendCommandScript( "sg %s" %( str( par.SG ) ) )
      elif par.SG_OPTIONS == 'laue':
        self.appendCommandScript( "sg all")
#-    out = self.container.outputData
#-    print str( out.XYZOUT.fullPath )
#-    print
#-    print '------------------------------------------'
#-    print self.commandLine
#-    print self.commandScript
#-    print '------------------------------------------'
#-    print
      return 0


    def processOutputFiles( self ) :
      out = self.container.outputData
      gui = self.container.guiParameters
      par = self.container.controlParameters

      
      # Rename output file or report failure
      if not par.SG_OPTIONS == 'laue' and not gui.PERFORM == 'srf':
        fileName = os.path.join( self.path_wrk, 'molrep.pdb' )
        #print 'molrep_mr.processOutputFiles XYZOUT',fileName,  os.path.isfile( fileName )
        if os.path.isfile( fileName ) :
           with open(fileName) as istream:
             content = istream.read()

           content = re.sub('\n#MOLECULE\s+[0-9]+\s*', '\n', content)
           with open(str(out.XYZOUT.fullPath), 'w') as ostream:
             ostream.write(content)

        else:
          self.appendErrorReport(104,str( out.XYZOUT ) , stack=False)
          return CPluginScript.FAILED
    
      if gui.PERFORM == 'pat':
        out.XYZOUT.annotation = 'Model from molrep MR'
      elif gui.PERFORM == 'den':
        out.XYZOUT.annotation = 'Model from molrep density search'

#-    file = os.path.join( self.path_wrk, 'molrep.xml' )
#-    if os.path.isfile( file ) :
#-       os.rename( file, xmlout )

#-    print str( out.XYZOUT.fullPath )

      docout = os.path.join( self.path_wrk, 'molrep.doc' )
      xmlout = str( self.makeFileName( 'PROGRAMXML' ) )
      self.saveProgramXml(docout,xmlout )
      os.rename(docout,docout+'.txt')
      

      if  gui.PERFORM == 'srf':
        psfile = os.path.join(self.getWorkDirectory(),'molrep_rf.ps')
        if os.path.exists(psfile):
          out.PSOUT.set(psfile)
          out.PSOUT.annotation.set('Self-rotation function')
      
      return CPluginScript.SUCCEEDED



    def saveProgramXml ( self, docFileName, programXmlFileName ) :
      from core import CCP4Utils
      titles = []
      status = 0
      results = ET.Element('MolrepResult')
      tf = ET.Element('MR_TF')
      results.append(tf)
      for key,value in [ ['err_level','0'],
                         ['err_message','normal termination'],
                         ['n_solution','1'],
                         ['mr_score','0.0000'] ]:
          
        e = ET.Element(key)
        e.text = value
        tf.append(e)

      '''
      table = [ '<?xml version="1.0" encoding="ASCII" standalone="yes"?>' ]
      table.append( "<MolrepResult>" )
      table.append( "<MR_TF>" )
      table.append( "Error <err_level>0</err_level>" )
      table.append( "Message <err_message>normal termination</err_message>" )
      table.append( "nmon_solution <n_solution>1</n_solution>" )
      table.append( "mr_score <mr_score>0.0000</mr_score>" )
      table.append( "</MR_TF>" )
      table.append( " <RFpeaks>" )
      '''

      docfileText = CCP4Utils.readFile(docFileName)
      docfileList = docfileText.split('\n')
      
      for line in docfileList :
         #print 'line in docfile',status,line
         if status == 0 :
            lstr = line.strip()
            lst1 = "--- Translation function ---"
            lst2 = "--- phased translation function ---"
            if lstr == lst1 or lstr == lst2 :
              status = 1

         elif status == 1 :
            if line.strip().startswith( 'RF ' ) :
              titles = line.replace( "(", " " ).replace( ")", "" ).replace( "/", "_" ).split()
              #print 'titles',titles
              rf = ET.Element('RFpeaks')
              results.append(rf)

            else:
              words = line.replace( "(", " " ).replace( ")", "" ).replace( "-", " -" ).split()
              if len( words ) == len( titles ) :
                try :
                  for i in (0,1): ii = int( words[i] )
                  for i in range(2,len(words)): ii = float( words[i] )
                  peak = ET.Element('RFpeak')
                  for key,value in zip( titles, words ) :
                    #print ' key,value', key,value
                    e = ET.Element(key)
                    e.text = str(float(value))
                    peak.append(e)
                except :
                  pass
                else:
                  rf.append(peak)
                 
      if self.container.controlParameters.SG_OPTIONS == 'laue':
        data = self.extractLaueDataFromLog()
        eleNames = ['space_group','score','contrast']
        eLaue = ET.Element('laue_group_alternatives')
        for d in data:
          eTest = ET.Element('test')
          eLaue.append(eTest)
          for ii in range(3):
            e = ET.Element(eleNames[ii])
            e.text = str(d[ii])
            eTest.append(e)     
        results.append(eLaue)
          
      with open(programXmlFileName,"w") as programXMLFile:
            CCP4Utils.writeXML(programXMLFile,ET.tostring(results))

    def extractLaueDataFromLog(self):
      from core import CCP4Utils
      try:
        text = CCP4Utils.readFile(self.makeFileName('LOG'))
      except:
        self.appendErrorReport(105,(self.makeFileName('LOG')))
        return                       
                               
      m = re.search(r'(.*)Space Group Checking(.*)',text,re.DOTALL)
      if m is None:
        self.appendErrorReport(106)
        return
      lines = m.groups()[1].split('\n')[4:]
      #print 'lines',lines
      data = []
      for line in lines:
        if line.find('+-')>=0: break
        spg = line[14:30].strip()
        score = line[44:49].strip()
        cntr = line[50:56].strip()
        data.append([spg,score,cntr])
        #print 'spg,score',spg,score
      return data

'''       
0123456789012345678901234567890123456789012345678901234567890
 +-------------------------------------------------------+
 |        Nsg     sg                        Score   Cntr |
 +-------------------------------------------------------+
 |    1   19   P 21 21 21                   0.656 24.887 |
 |    2   16   P 2 2 2                      0.275  5.305 |
 |    3   17   P 2 2 21                     0.289 14.153 |
 |    4 1017   P 21 2 2                     0.288  9.549 |
 |    5 2017   P 2 21 2                     0.293  9.750 |
 |    6   18   P 21 21 2                    0.312  5.950 |
 |    7 2018   P 21 2 21                    0.288  4.563 |
 |    8 3018   P 2 21 21                    0.296 13.246 |
 +-------------------------------------------------------+
'''       

