from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript


class privateer(CPluginScript):
    TASKNAME = "privateer"
    TASKCOMMAND = 'privateer'
    WHATNEXT = [ 'coot_rebuild', 'prosmart_refmac' ]

    def processInputFiles(self):
      from ccp4i2.core import CCP4XtalData

      #print 'taskMakeHklin F_SIGF',self.container.inputData.F_SIGF,type(self.container.inputData.F_SIGF),self.container.inputData.F_SIGF.contentFlag
      self.hklin,error = self.makeHklin ( [ ['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] ] )
      if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING: return CPluginScript.FAILED

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        print(self.WHATNEXT)
        import os
        out = self.container.outputData
        self.path_wrk = str( self.getWorkDirectory() )

        fileName = os.path.join ( self.path_wrk, 'FPHIOUT.mtz' )
        library  = os.path.join ( self.path_wrk, 'privateer-lib.cif' )
        keywords = os.path.join ( self.path_wrk, 'keywords_refmac5.txt' )
        print(keywords)

        if os.path.isfile ( fileName ) :
            out.FPHIOUT.set ( fileName )
            out.FPHIOUT.annotation.set ( '2mFo-DFc map coefficients' )
            out.FPHIOUT.subType = 1

        fileName2 = os.path.join ( self.path_wrk, 'OMITFPHIOUT.mtz' )

        if os.path.isfile ( fileName2 ) :
            out.DIFFPHIOUT.set ( fileName2 )
            out.DIFFPHIOUT.annotation.set ( 'omit mFo-DFc map coefficients' )
            out.DIFFPHIOUT.subType = 2

        fileName3 = os.path.join ( self.path_wrk, 'privateer-results.py' )

        if os.path.isfile ( fileName3 ) :
            out.COOTSCRIPTOUT.set ( fileName3 )
            out.COOTSCRIPTOUT.annotation.set ( 'guided tour on the reported issues' )

        if os.path.isfile (library) :
            out.LIBOUT.set ( library )
            out.LIBOUT.annotation.set ( 'Dictionary containing unimodal torsion restraints for Coot and REFMAC5' )

        if os.path.isfile ( keywords ) :
            out.REFMAC_KEYWORD_FILE.set ( keywords )
            out.REFMAC_KEYWORD_FILE.annotation.set ( 'Keywords for activating unimodal torsion restraints in REFMAC5' )

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
      self.appendCommandLine(['-stdin'])

      self.appendCommandScript( "pdbin %s"%(str(self.container.inputData.XYZIN)))

      # INPUT DATA
      self.appendCommandScript( "mtzin "+self.hklin )
      self.appendCommandScript( "mode ccp4i2" )
      self.appendCommandScript( "colin-fo F,SIGF")

      # CONTROL PARAMETERS

      if self.container.controlParameters.SHOWGEOM:
        self.appendCommandScript("showgeom")

      if self.container.controlParameters.RADIUSIN.isSet():
        self.appendCommandScript("radiusin %s"%(str(self.container.controlParameters.RADIUSIN)))

      if self.container.controlParameters.NEW_SUGAR :
        self.appendCommandScript("valstring %s,%s/%s/%s/%s/%s/%s,%s,%s,%s"%(self.container.controlParameters.CODEIN,self.container.controlParameters.RING_OXYGEN,self.container.controlParameters.RING_C1,self.container.controlParameters.RING_C2,self.container.controlParameters.RING_C3,self.container.controlParameters.RING_C4,self.container.controlParameters.RING_C5,self.container.controlParameters.ANOMER,self.container.controlParameters.HAND,self.container.controlParameters.CONFORMATION_PYRANOSE ))
        self.appendCommandScript("codein %s"%(self.container.controlParameters.CODEIN))

      if self.container.controlParameters.VERTICAL == "vertical" :
        self.appendCommandScript("vertical")

      if self.container.controlParameters.ESSENTIALS == "essentials" :
        self.appendCommandScript("essentials")

      if self.container.controlParameters.INVERT == "white" :
        self.appendCommandScript("invert")

      if self.container.controlParameters.OLDSTYLE == "original":
        self.appendCommandScript("oldstyle")

      if self.container.controlParameters.BLOBS :
        self.appendCommandScript("check-unmodelled")
        self.appendCommandScript("blobs_threshold %s"%(str(self.container.controlParameters.BLOBSLEVEL)))
        
      if self.container.controlParameters.GLYTOUCAN:
        self.appendCommandScript("glytoucan")
      
      if self.container.controlParameters.CLOSESTMATCH:
        self.appendCommandScript("closest_match_disable")

      if self.container.controlParameters.ALLPERMUTATIONS:
        self.appendCommandScript("all_permutations")

      if self.container.controlParameters.NUMTHREADS.isSet():
        self.appendCommandScript("cores %s"%(str(self.container.controlParameters.NUMTHREADS)))

      # Need to a fix bug here where I cant put GLyToucan after expression and expression before glytoucan

      if self.container.controlParameters.EXPRESSION != "undefined" :
        self.appendCommandScript("expression %s" % (self.container.controlParameters.EXPRESSION) )


#      OUTPUT DATA
#      if self.container.outputData.XYZOUT.isSet():
#          self.appendCommandScript("pdbout %s"%(str(self.container.outputData.XYZOUT.fullPath)))
#      self.appendCommandScript("xmlout %s"%(self.makeFileName('PROGRAMXML')))

      return CPluginScript.SUCCEEDED
