from __future__ import print_function


from core.CCP4PluginScript import CPluginScript

class refmac_wrk( CPluginScript ) :

    TASKTITLE = 'Refmac Replica'
    TASKNAME = 'refmac_wrk'
    TASKCOMMAND = 'refmac5'
    TASKVERSION = 0.0
    TASKMODULE = 'test'

    def makeCommandAndScript( self ) :

      inp = self.container.inputData
      par = self.container.controlParameters
      out = self.container.outputData

      from core import CCP4Utils
      import os
#     self.path_wrk = str( self.getWorkDirectory() )
#     self.path_scr = os.path.join( self.path_wrk, 'scratch' )
#     os.mkdir( self.path_scr )
#     print "self.path_wrk", self.path_wrk
#     print "self.path_scr", self.path_scr

#     xmlout = str( out.XMLOUT.fullPath )
      xmlout = self.makeFileName( 'PROGRAMXML' )

      self.appendCommandLine( [ 'XYZIN', inp.XYZIN.fullPath ] )
      self.appendCommandLine( [ 'HKLIN', inp.HKLIN.fullPath ] )
      if par.IFFIXTLS or str(par.REFINE_TYPE) == 'TLS':
          self.appendCommandLine( [ 'TLSIN', inp.TLSIN.fullPath ] )

      if inp.MAKE_LIBRARY.isSet() :
          self.appendCommandLine( [ 'LIBIN', inp.MAKE_LIBRARY.fullPath ] )

      self.appendCommandLine( [ 'XYZOUT', out.XYZOUT.fullPath ] )
      self.appendCommandLine( [ 'HKLOUT', out.HKLOUT.fullPath ] )
      if str(par.REFINE_TYPE) == 'TLS':
          self.appendCommandLine( [ 'TLSOUT', out.TLSOUT.fullPath ] )

      if inp.MAKE_LIBRARY.isSet() :
          self.appendCommandLine( [ 'LIBOUT', out.LIBOUT.fullPath ] )

      self.appendCommandLine( [ 'XMLOUT', xmlout ] )

      lines = []

      if par.TITLE.isSet():
          lines.append("TITLE %s" %(str(par.TITLE)))

      if par.EXCLUDE_RESOLUTION:
          lines.append("refi resolution %s %s" %(str(par.EXCLUDE_RESOLUTION_MAX), str(par.EXCLUDE_RESOLUTION_MIN)))

      if str(par.REFINE_TYPE) == 'REVIEW':
          lines.append("make check %s" %str(par.REVIEW_MAKE_CHECK))
          lines.append("make -")
          lines.append("   hydrogen %s -" %str(par.REVIEW_MAKE_HYDROGEN))
          if par.REVIEW_MAKE_HOUT:
              lines.append("   hout YES -" )

          else:
              lines.append("   hout NO -" )

          if par.REVIEW_MAKE_PEPTIDE:
              lines.append("   peptide YES -" )

          else:
              lines.append("   peptide NO -" )

          if par.REVIEW_MAKE_CISPEPTIDE:
              lines.append("   cispeptide YES -" )

          else:
              lines.append("   cispeptide NO -" )

          if par.REVIEW_MAKE_SSBRIDGE:
              lines.append("   ssbridge YES -" )

          else:
              lines.append("   ssbridge NO -" )

          if par.REVIEW_MAKE_SYMMETRY:
              lines.append("   symmetry YES -" )

          else:
              lines.append("   symmetry NO -" )

          lines.append("   sugar %s -" %str(par.REVIEW_MAKE_SUGAR))
          lines.append("   connectivity %s -" %str(par.REVIEW_MAKE_CONNECTIVITY))
          lines.append("   link %s -" %str(par.REVIEW_MAKE_LINK))
          lines.append("   format f exit Y" )

      else:
          lines.append("make check %s" %str(par.MAKE_CHECK))
          lines.append("make -")
          lines.append("   hydrogen %s -" %str(par.MAKE_HYDROGEN))
          if par.MAKE_HOUT:
              lines.append("   hout YES -" )

          else:
              lines.append("   hout NO -" )

          if par.MAKE_PEPTIDE:
              lines.append("   peptide YES -" )

          else:
              lines.append("   peptide NO -" )

          if par.MAKE_CISPEPTIDE:
              lines.append("   cispeptide YES -" )

          else:
              lines.append("   cispeptide NO -" )

          if par.MAKE_SSBRIDGE:
              lines.append("   ssbridge YES -" )

          else:
              lines.append("   ssbridge NO -" )

          if par.MAKE_SYMMETRY:
              lines.append("   symmetry YES -" )

          else:
              lines.append("   symmetry NO -" )

          lines.append("   sugar %s -" %str(par.MAKE_SUGAR))
          lines.append("   connectivity %s -" %str(par.MAKE_CONNECTIVITY))
          lines.append("   link %s" %str(par.MAKE_LINK))

      if not str(par.REFINE_TYPE) == 'REVIEW':
          if par.IFAUTONCS:
              lines.append("ncsr %s" %str(par.AUTONCS_MODE))

          lines.append("refi -")
          if str(par.REFINE_TYPE) == 'TLS':
              lines.append("   type REST -")

          else:
              lines.append("   type %s -" %str(par.REFINE_TYPE))

          if str(par.INPUT_PHASE) in ["PHASE","HL"]:
              lines.append("   phase scbl %s bblu %s -" %(str(par.PHASE_SCBLUR),str(par.PHASE_BBLUR)))

          lines.append("   resi MLKF meth CGMAT -")
          if str(par.REFINE_TYPE) in ["IDEA","RIGID"]:
              line = '   bref OVER'

          else:
              line = '   bref %s' %str(par.B_REFINEMENT_MODE)

          if str(par.REFINE_TYPE) == 'TLS':
              lines.append(line + ' -')
              line = '   tlsc %s' %str(par.TLS_NCYCLES)

          lines.append(line)

          if str(par.INPUT_PHASE) == 'NO' and not str(par.TWINREF_TYPE) == "NO":
              lines.append('twin')

          if str(par.INPUT_PHASE) == 'SAD':
              if not par.REF_SUBOCC:
                  lines.append('refi oref no')

              if par.WAVELENGTH.isSet():
                  lines.append("anom wave %s" %str(par.WAVELENGTH))

              for atom in par.ANOMALOUS_ATOMS:
                  if atom.atomType.isSet() and atom.Fp.isSet() and atom.Fpp.isSet():
                      lines.append("anom form %s %s %s" %(str(atom.atomType), str(atom.Fp), str(atom.Fpp)))

          if not str(par.REFINE_TYPE) == 'RIGID':
              lines.append("ncyc %s" %str(par.NCYCLES))

          else:
              lines.append("rigid ncycle %s" %str(par.RIGID_NCYCLES))
              ind = 0
              for rigid_group in par.RIGID_GROUP_LIST :
                  if rigid_group.segmentList :
                     ind = ind + 1
                     line = 'rigid group %s' %str(ind)
                     for segment in rigid_group.segmentList :
                         lines.append(line + ' -')
                         str_seg = (str(segment.residue_1), str(segment.chain_id), str(segment.residue_2), str(segment.chain_id))
                         line = '   from %s %s to %s %s' %str_seg

                     lines.append(line)

          if par.BLIM:
              lines.append('blim %s %s' %(str(par.BLIM_MIN),str(par.BLIM_MAX)))

          lines.append('scal -')
          if str(par.BULK_SOLVENT_SCALING) == 'BULK':
              lines.append('   type BULK -')
              if par.BULK_SCALING_RESOLUTION_MIN.isSet() and par.BULK_SCALING_RESOLUTION_MAX.isSet():
                 lines.append('   reso %s %s -' %(str(par.BULK_SCALING_RESOLUTION_MAX),str(par.BULK_SCALING_RESOLUTION_MIN)))

          else:
              lines.append('   type SIMP -')
              if par.SIMPLE_SCALING_RESOLUTION_MIN.isSet() and par.SIMPLE_SCALING_RESOLUTION_MAX.isSet():
                 lines.append('   reso %s %s -' %(str(par.SIMPLE_SCALING_RESOLUTION_MAX),str(par.SIMPLE_SCALING_RESOLUTION_MIN)))

          sublines = []
          subline = '   LSSC'
          if not str(par.REFINE_TYPE) == 'RIGID':
              sublines.append(subline + ' -')
              subline = '   ANISO'

          if par.SCALING_IF_FIXB:
              sublines.append(subline + ' -')
              subline = '   FIXB bvalue %s' %str(par.SCALING_FIXB_BBULK)

          if par.SCALING_EXPE_SIGMA:
              sublines.append(subline + ' -')
              subline = '   EXPE'

          if str(par.SCALING_REF_SET) == 'FREE':
              sublines.append(subline + ' -')
              subline = '   FREE'

          sublines.append(subline)
          if len(sublines) > 1:
              for subline in sublines:
                  lines.append(subline)

          if par.EXCLUDE_FREER and not str(par.EXCLUDE_FREER_VALUE) == '0':
              lines.append("free %s" %str(par.EXCLUDE_FREER_VALUE))

          line = 'solvent NO'
          if par.IF_SOLVENT:
              line = 'solvent YES'

          if par.SOLVENT_VDWPROB.isSet():
              lines.append(line + ' -')
              line = '   VDWProb %s' %str(par.SOLVENT_VDWPROB)

          if par.SOLVENT_IONPROB.isSet():
              lines.append(line + ' -')
              line = '   IONProb %s' %str(par.SOLVENT_IONPROB)

          if par.SOLVENT_RSHRINK.isSet():
              lines.append(line + ' -')
              line = '   RSHRink %s' %str(par.SOLVENT_RSHRINK)

          lines.append(line)

          lines.append('weight -')
          if not par.EXPERIMENTAL_WEIGHTING:
              lines.append('   NOEX -')

          if par.AUTO_WEIGHTING:
              lines.append('   AUTO')

          else:
              lines.append('   MATRIX %s' %str(par.MATRIX_WEIGHT))

          lines.append('monitor %s -' %str(par.MONI_LEVEL))
          lines.append('   torsion %s -' %str(par.MONI_TORSION))
          lines.append('   distance %s -' %str(par.MONI_DISTANCE))
          lines.append('   angle %s -' %str(par.MONI_ANGLE))
          lines.append('   plane %s -' %str(par.MONI_PLANE))
          lines.append('   vanderwaals %s -' %str(par.MONI_VANDERWAALS))
          lines.append('   chiral %s -' %str(par.MONI_CHIRAL))
          lines.append('   bfactor %s -' %str(par.MONI_BFACTOR))
          lines.append('   bsphere %s -' %str(par.MONI_BSPHERE))
          lines.append('   rbond %s -' %str(par.MONI_RBOND))
          lines.append('   ncsr %s' %str(par.MONI_NCSR))

          if not str(par.REFINE_TYPE) == 'IDEA':

              line = 'labin'
              if str(par.INPUT_PHASE) == 'SAD':
                  lines.append(line + ' -')
                  lines.append('   F+=%s -' %str(inp.F_ANO.F_plus))
                  lines.append('   SIGF+=%s -' %str(inp.F_ANO.SIGF_plus))
                  lines.append('   F-=%s -' %str(inp.F_ANO.F_minus))
                  line = '   SIGF-=%s' %str(inp.F_ANO.SIGF_minus)

              elif str(par.INPUT_PHASE) == 'NO' and str(par.TWINREF_TYPE) == 'INTENSITIES':
                  lines.append(line + ' -')
                  lines.append('   IP=%s -' % str(inp.IOBS.IP))
                  line = '   SIGIP=%s' % str(inp.IOBS.Sigma)

              else:
                  lines.append(line + ' -')
                  lines.append('   FP=%s -' % str(inp.FOBS.FP))
                  line = '   SIGFP=%s' % str(inp.FOBS.Sigma)

              if str(par.INPUT_PHASE) == 'HL':
                  lines.append(line + ' -')
                  lines.append('   HLA=%s -' % str(inp.ABCD.HLA))
                  lines.append('   HLB=%s -' % str(inp.ABCD.HLB))
                  lines.append('   HLC=%s -' % str(inp.ABCD.HLC))
                  line = '   HLD=%s -' % str(inp.ABCD.HLD)

              elif str(par.INPUT_PHASE) == 'PHASE':
                  lines.append(line + ' -')
                  lines.append('   PHIB=%s -' % str(inp.PHOBS.PHIB))
                  line = '   FOM=%s' % str(inp.PHOBS.FOM)

              if par.EXCLUDE_FREER and inp.FREERFLAG.isSet():
                  lines.append(line + ' -')
                  line = '   FREE=%s' % str(inp.FREERFLAG.FREE)

              lines.append(line)

          if par.IFTMP or par.IFBFAC_SET and str(par.REFINE_TYPE) == 'TLS':
              line = 'temp'
              if par.IFTMP:
                  line = line + ' %s %s %s %s %s' %(str(par.WBSKAL), str(par.SIGB1), str(par.SIGB2), str(par.SIGB3), str(par.SIGB4))

              if par.IFBFAC_SET and str(par.REFINE_TYPE) == 'TLS':
                  line = line + ' set %s' %str(par.BFAC_SET)

              lines.append(line)

          if par.IFADDU:
              lines.append('tlso addu')

          if par.IFSHARP:
              line = 'mapc shar'
              if par.IFBSHARP:
                  lines.append(line + ' -')
                  line = '   %s' %str(par.B_SHARP)

              if par.IFALSHARP:
                  lines.append(line + ' -')
                  line = '   alpha %s' %str(par.AL_SHARP)

              lines.append(line)


          if par.IFJELLY:
              lines.append('ridg dist sigm %s' %(str(par.JELLY_SIGMA)))

          if par.IFDIST:
              lines.append('dist %s %s %s %s %s %s'
                  %(str(par.WDSKAL), str(par.SIGD1), str(par.SIGD2),
                    str(par.SIGD3), str(par.SIGD4), str(par.SIGD5)))

          if par.IFANGL:
              lines.append('angle %s'
                  %(str(par.ANGLE_SCALE)))

          if par.IFPLAN:
              lines.append('plan %s %s %s'
                  %(str(par.WPSKAL), str(par.SIGPP), str(par.SIGPA)))

          if par.IFCHIR:
              lines.append('chir %s %s'
                  %(str(par.WCSKAL), str(par.SIGC)))

          if par.IFNCSR:
              lines.append('ncsr %s %s %s %s %s %s %s'
                  %(str(par.WSSKAL), str(par.SIGSP1), str(par.SIGSP2),
                    str(par.SIGSP3), str(par.SIGSB1), str(par.SIGSB2),
                    str(par.SIGSB3)))

          if par.IFTORS:
              lines.append('torsion %s %s %s %s %s'
                  %(str(par.WTSKAL), str(par.SIGT1), str(par.SIGT2),
                    str(par.SIGT3), str(par.SIGT4)))

          if par.IFVAND:
              sublines = []
              subline = 'vand'
              if par.WVSKAL.isSet():
                  sublines.append(subline + ' -')
                  subline = 'overall %s' %str(par.WVSKAL)

              if par.WAND_SIGMA_VDW.isSet():
                  sublines.append(subline + ' -')
                  subline = ' sigma vdw %s' %str(par.WAND_SIGMA_VDW)

              if par.WAND_SIGMA_HBOND.isSet():
                  sublines.append(subline + ' -')
                  subline = 'sigma hbond %s' %str(par.WAND_SIGMA_HBOND)

              if par.WAND_SIGMA_METAL.isSet():
                  sublines.append(subline + ' -')
                  subline = ' sigma metal %s' %str(par.WAND_SIGMA_METAL)

              if par.WAND_SIGMA_TORS.isSet():
                  sublines.append(subline + ' -')
                  subline = 'sigma tors %s' %str(par.WAND_SIGMA_TORS)

              if par.WAND_INCR_TORS.isSet():
                  sublines.append(subline + ' -')
                  subline = 'incr tors %s' %str(par.WAND_INCR_TORS)

              if par.WAND_INCR_ADHB.isSet():
                  sublines.append(subline + ' -')
                  subline = 'incr adhb %s' %str(par.WAND_INCR_ADHB)

              if par.WAND_INCR_AHHB.isSet():
                  sublines.append(subline + ' -')
                  subline = 'incr ahhb %s' %str(par.WAND_INCR_AHHB)

              sublines.append(subline)
              if len(sublines) > 1:
                  for subline in sublines:
                      lines.append(subline)

          if par.IFISO and str(par.B_REFINEMENT_MODE) == 'ANIS':
              if par.SPHERICITY.isSet():
                  lines.append('sphericity %s' %str(par.SPHERICITY))

              if par.RBOND.isSet():
                  lines.append('rbond %s' %str(par.RBOND))

      if par.IFPROSMART:
          if par.RESTRAINTFILE.isSet() and par.RESTRAINTFILE.exists():
              lines.append('@%s' %str(par.RESTRAINTFILE.fullPath))

              if par.IFEXTREST_SCALE:
                  lines.append('EXTERNAL WEIGHT SCALE %s' %str(par.EXTREST_SCALE))

              if par.IFEXTREST_USEMAIN:
                  lines.append('EXTERNAL USE MAIN')

              if par.IFEXTREST_GMWT:
                  lines.append('EXTERNAL WEIGHT GMWT %s' %str(par.EXTREST_GMWT))

              if par.IFEXTREST_DMAX:
                  lines.append('EXTERNAL DMAX %s' %str(par.EXTREST_DMAX))

      if par.INCLUDEFILE.isSet() and par.INCLUDEFILE.exists():
              lines.append('@%s' %str(par.INCLUDEFILE.fullPath))

      lines.append('end')
      for line in lines:
          self.appendCommandScript(line)

      return 0

# -----------------------------------------------------------------------------------
import unittest

class testRefmac( unittest.TestCase ) :

#  def setUp( self ) :
#  def tearDown( self ) :
#  def test2( self ) :

   def test1( self ) :
      self._pattern( 1 )

   def test2( self ) :
      self._pattern( 2 )

   def test3( self ) :
      self._pattern( 3 )

   def test4( self ) :
      self._pattern( 4 )

   def test5( self ) :
      self._pattern( 5 )

   def test6( self ) :
      self._pattern( 6 )

   def test7( self ) :
      self._pattern( 7 )

   def test8( self ) :
      self._pattern( 8 )

   def test9( self ) :
      self._pattern( 9 )

   def _pattern( self, no ) :

      from core.CCP4Utils import getCCP4I2Dir
      import os

      xmlInput = os.path.join( getCCP4I2Dir(), 'wrappers', 'refmac_wrk', 'test_data', 'test'+str(no)+'.params.xml' )
      self.wrapper = refmac_wrk( name='job' )
      self.wrapper.container.loadDataFromXml( xmlInput )
      self.wrapper.setWaitForFinished( 1000000 )

      pid = self.wrapper.process()
      self.wrapper.setWaitForFinished( -1 )
      if len(self.wrapper.errorReport) > 0 :
         print(self.wrapper.errorReport.report())

def TESTSUITE() :

   suite = unittest.TestLoader().loadTestsFromTestCase( testRefmac )
   return suite

def testModule() :

   suite = TESTSUITE()
   unittest.TextTestRunner( verbosity=2 ).run( suite )

