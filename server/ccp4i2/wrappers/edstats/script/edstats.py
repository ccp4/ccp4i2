from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4ClipperUtils import is_aminoacid
from ccp4i2.core.CCP4PluginScript import CPluginScript


class edstats(CPluginScript):

    TASKMODULE          = 'validation'          # Where this plugin will appear on the gui
    TASKTITLE           = 'Measure agreement between model and density'  # A short title for gui menu
    TASKNAME            = 'edstats'  # Task name - should be same as class name
    TASKCOMMAND         = 'edstats'  # The command to execute, should be reachable
    DESCRIPTION         = 'Calculates real-space metrics for evaluating the agreement between model and density (Edstats)'
    TASKVERSION         = 0.1                   # Version of this plugin
    WHATNEXT            = ['coot_rebuild']
    PURGESEARCHLIST =  [[ 'fft%*/MAPOUT.map', 1 ] ]
    MAINTAINER = 'jon.agirre@york.ac.uk'

    def processInputFiles(self):
      from ccp4i2.core import CCP4XtalData

      self.cfftPlugin1 = self.makeCfftPlugin1 ( )
      error = self.cfftPlugin1.process ( )
      if error == CPluginScript.FAILED:
          self.reportStatus ( error )
      else :
          self.container.inputData.MAPIN1 = self.cfftPlugin1.container.outputData.MAPOUT

      self.cfftPlugin2 = self.makeCfftPlugin2 ( )
      error = self.cfftPlugin2.process ( )
      if error == CPluginScript.FAILED:
          self.reportStatus ( error )
      else :
          self.container.inputData.MAPIN2 = self.cfftPlugin2.container.outputData.MAPOUT


      #print 'taskMakeHklin F_SIGF',self.container.inputData.F_SIGF,type(self.container.inputData.F_SIGF),self.container.inputData.F_SIGF.contentFlag
      #self.hklin,error = self.makeHklin ( [ ['F_SIGF',CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN ] ] )
      #if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING: return CPluginScript.FAILED

      return CPluginScript.SUCCEEDED

    def processOutputFiles(self):

        import os
        out = self.container.outputData
        self.path_wrk = str( self.getWorkDirectory() )

        fileName3 = os.path.join ( self.path_wrk, 'coot_script.py' )

        if os.path.isfile ( fileName3 ) :
            out.COOTSCRIPTOUT.set ( fileName3 )
            out.COOTSCRIPTOUT.annotation.set ( 'guided tour on the reported issues' )

        from lxml import etree
        xmlRoot = etree.Element('Edstats')
        segmentNode = None

        lines = open(self.makeFileName('LOG')).readlines()

        readingTable = False

        mcBadBits     = []
        scBadBits     = []
        ligandBadBits = []
        waterBadBits  = []

        zdm = float ( self.container.controlParameters.SIGMA_RZ_MINUS )
        zdp = float ( self.container.controlParameters.SIGMA_RZ_PLUS  )
        zo  = float ( self.container.controlParameters.SIGMA_RO )

        for line in lines:

            if readingTable :
                if len( line.strip().split ( ) ) > 20:
                    residueNode = etree.SubElement(xmlRoot,'Residue')

                    resnameNode     = etree.SubElement(residueNode,'Name')
                    resnameNode.text = line.strip().split( )[0]

                    resChainNode      = etree.SubElement(residueNode,'Chain')
                    resChainNode.text = line.strip().split( )[1]

                    resNumberNode      = etree.SubElement(residueNode,'Number')
                    resNumberNode.text = line.strip().split( )[2]

                    resBAmNode    = etree.SubElement(residueNode,'BAm')
                    resBAmNode.text = line.strip().split( )[3]

                    resZCCmNode    = etree.SubElement(residueNode,'ZCCm')
                    resZCCmNode.text = line.strip().split( )[10]

                    resZOmNode    = etree.SubElement(residueNode,'ZOm')
                    resZOmNode.text = line.strip().split( )[11]

                    resZDmmNode    = etree.SubElement(residueNode,'ZDmm')
                    resZDmmNode.text = line.strip().split( )[13]

                    resZDpmNode    = etree.SubElement(residueNode,'ZDpm')
                    resZDpmNode.text = line.strip().split( )[14]

                    resBAsNode    = etree.SubElement(residueNode,'BAs')
                    resBAsNode.text = line.strip().split( )[15]

                    resZCCsNode    = etree.SubElement(residueNode,'ZCCs')
                    resZCCsNode.text = line.strip().split( )[22]

                    resZOsNode    = etree.SubElement(residueNode,'ZOs')
                    resZOsNode.text = line.strip().split( )[23]

                    resZDmsNode    = etree.SubElement(residueNode,'ZDms')
                    resZDmsNode.text = line.strip().split( )[25]

                    resZDpsNode    = etree.SubElement(residueNode,'ZDps')
                    resZDpsNode.text = line.strip().split( )[26]

                    resBAaNode    = etree.SubElement(residueNode,'BAa')
                    resBAaNode.text = line.strip().split( )[27]

                    resCCPaNode    = etree.SubElement(residueNode,'CCPa')
                    resCCPaNode.text = line.strip().split( )[33]

                    resZCCaNode    = etree.SubElement(residueNode,'ZCCa')
                    resZCCaNode.text = line.strip().split( )[34]

                    resZOaNode    = etree.SubElement(residueNode,'ZOa')
                    resZOaNode.text = line.strip().split( )[35]

                    resZDmaNode    = etree.SubElement(residueNode,'ZDma')
                    resZDmaNode.text = line.strip().split( )[37]

                    resZDpaNode    = etree.SubElement(residueNode,'ZDpa')
                    resZDpaNode.text = line.strip().split( )[38]

                    diagnostic = "Nothing"
                    potential_fix = "Refine?"

                    if is_aminoacid ( resnameNode.text ) :

                        # Now we create a list with main chain outliers to be fixed within Coot
                        if resZDmmNode.text != "n/a" and resZDpmNode.text != "n/a" and resZOmNode.text != "n/a" :
                            if (
                                 float ( resZDmmNode.text ) < zdm or
                                 float ( resZDpmNode.text ) > zdp or
                                 float ( resZOmNode.text )  < zo
                               ) :

                                if (
                                     float ( resZDpmNode.text ) > zdp
                                   ) :
                                    diagnostic = "This peptide needs to be flipped in order to fit the density"
                                    potential_fix = "pep_flip"

                                badResidue = { "name"          :resnameNode.text,
                                               "id"            :resNumberNode.text,
                                               "chain"         :resChainNode.text,
                                               "minus"         :resZDmmNode.text,
                                               "plus"          :resZDpmNode.text,
                                               "zo"            :resZOmNode.text,
                                               "rscc"          :resZCCmNode.text,
                                               "bfact"         :resBAmNode.text,
                                               "diagnostic"    :diagnostic,
                                               "potential_fix" :potential_fix
                                             }
                                mcBadBits.append ( badResidue )

                        # And the same for side chains
                        if resZDmsNode.text != "n/a" and resZDpsNode.text != "n/a" and resZOsNode.text != "n/a" :
                            if float ( resZDmsNode.text ) < zdm or float ( resZDpsNode.text ) > zdp or float ( resZOsNode.text ) < zo :
                                badResidue = { "name"          :resnameNode.text,
                                               "id"            :resNumberNode.text,
                                               "chain"         :resChainNode.text,
                                               "minus"         :resZDmsNode.text,
                                               "plus"          :resZDpsNode.text,
                                               "zo"            :resZOsNode.text,
                                               "rscc"          :resZCCsNode.text,
                                               "bfact"         :resBAsNode.text,
                                               "diagnostic"    :diagnostic,
                                               "potential_fix" :potential_fix
                                             }
                                scBadBits.append ( badResidue )
                    # now waters
                    elif resnameNode.text == "HOH" :
                        if float ( resZDmaNode.text ) < zdm or float ( resZDpaNode.text ) > zdp or float ( resZOaNode.text ) < zo :

                            bad_water = { "name"           :resnameNode.text,
                                           "id"            :resNumberNode.text,
                                           "chain"         :resChainNode.text,
                                           "minus"         :resZDmmNode.text,
                                           "plus"          :resZDpmNode.text,
                                           "zo"            :resZOaNode.text,
                                           "rscc"          :resZCCaNode.text,
                                           "bfact"         :resBAaNode.text,
                                           "diagnostic"    :diagnostic,
                                           "potential_fix" :potential_fix
                                        }
                            waterBadBits.append ( bad_water )
                    else : # and (quite possibly) ligands
                        if float ( resZDmaNode.text ) < zdm or float ( resZDpaNode.text ) > zdp or float ( resZOaNode.text ) < zo :
                            bad_ligand = { "name"          :resnameNode.text,
                                           "id"            :resNumberNode.text,
                                           "chain"         :resChainNode.text,
                                           "minus"         :resZDmmNode.text,
                                           "plus"          :resZDpmNode.text,
                                           "zo"            :resZOaNode.text,
                                           "rscc"          :resZCCaNode.text,
                                           "bfact"         :resBAaNode.text,
                                           "diagnostic"    :diagnostic,
                                           "potential_fix" :potential_fix
                                         }
                            ligandBadBits.append ( bad_ligand )

            elif line.strip().startswith('RT  CI'):
                readingTable = True

        with open(self.container.outputData.COOTSCRIPTOUT.fullPath.__str__(),"w") as cootscript:
            if len(mcBadBits) > 0:
                interestingMCBitsDef = 'main_chain_outliers = ['
                for res in mcBadBits:
                    interestingMCBitsDef += ('''{"name": "%s",
                                                 "id": "%s",
                                                 "chain":"%s",
                                                 "minus":"%s",
                                                 "plus":"%s",
                                                 "zo":"%s",
                                                 "rscc":"%s",
                                                 "bfact":"%s",
                                                 "diagnostic":"%s",
                                                 "potential_fix":"%s"
                                                }\n,'''%(res['name'],res['id'],res['chain'],
                                                         res['minus'],res['plus'],res['zo'],
                                                         res['rscc'],res['bfact'],res['diagnostic'],
                                                         res['potential_fix']))

                interestingMCBitsDef += ']\n'

                cootscript.write(interestingMCBitsDef)
                cootscript.write('ccp4i2Interface.add_consolidated_menu(title="Main chain", interesting_bits=main_chain_outliers)\n')

            if len(scBadBits) > 0:
                interestingSCBitsDef = 'side_chain_outliers = ['
                for res in scBadBits:
                    interestingSCBitsDef += ('''{"name": "%s",
                                                 "id": "%s",
                                                 "chain":"%s",
                                                 "minus":"%s",
                                                 "plus":"%s",
                                                 "zo":"%s",
                                                 "rscc":"%s",
                                                 "bfact":"%s",
                                                 "diagnostic":"%s",
                                                 "potential_fix":"%s"
                                                }\n,'''%(res['name'],res['id'],res['chain'],
                                                         res['minus'],res['plus'],res['zo'],
                                                         res['rscc'],res['bfact'],res['diagnostic'],
                                                         res['potential_fix']))

                interestingSCBitsDef += ']\n'

                cootscript.write(interestingSCBitsDef)
                cootscript.write('ccp4i2Interface.add_consolidated_menu(title="Side chain", interesting_bits=side_chain_outliers)\n')

            if len(ligandBadBits) > 0:
                interestingLigandBitsDef = 'ligand_outliers = ['
                for res in ligandBadBits:
                    interestingLigandBitsDef += ('''{"name": "%s",
                                                 "id": "%s",
                                                 "chain":"%s",
                                                 "minus":"%s",
                                                 "plus":"%s",
                                                 "zo":"%s",
                                                 "rscc":"%s",
                                                 "bfact":"%s",
                                                 "diagnostic":"%s",
                                                 "potential_fix":"%s"
                                                }\n,'''%(res['name'],res['id'],res['chain'],
                                                         res['minus'],res['plus'],res['zo'],
                                                         res['rscc'],res['bfact'],res['diagnostic'],
                                                         res['potential_fix']))

                interestingLigandBitsDef += ']\n'

                cootscript.write(interestingLigandBitsDef)
                cootscript.write('ccp4i2Interface.add_consolidated_menu(title="Ligands", interesting_bits=ligand_outliers)\n')

            if len(waterBadBits) > 0:
                interestingWaterBitsDef = 'water_outliers = ['
                for res in waterBadBits:
                    interestingWaterBitsDef += ('''{"name": "%s",
                                                 "id": "%s",
                                                 "chain":"%s",
                                                 "minus":"%s",
                                                 "plus":"%s",
                                                 "zo":"%s",
                                                 "rscc":"%s",
                                                 "bfact":"%s",
                                                 "diagnostic":"%s",
                                                 "potential_fix":"%s"
                                                }\n,'''%(res['name'],res['id'],res['chain'],
                                                         res['minus'],res['plus'],res['zo'],
                                                         res['rscc'],res['bfact'],res['diagnostic'],
                                                         res['potential_fix']))

                interestingWaterBitsDef += ']\n'

                cootscript.write(interestingWaterBitsDef)
                cootscript.write('ccp4i2Interface.add_consolidated_menu(title="Waters", interesting_bits=water_outliers)\n')


        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            xmlString = etree.tostring(xmlRoot, pretty_print=True)
            CCP4Utils.writeXML(xmlFile,xmlString)

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
      import os

      from ccp4i2.core import CCP4XtalData

      self.path_wrk = str( self.getWorkDirectory() )
      edstatsOut = os.path.join ( self.path_wrk, 'edstats.out' )


      self.appendCommandLine([ 'XYZIN', self.container.inputData.XYZIN.fullPath ])

      # INPUT DATA
      self.appendCommandLine([ 'MAPIN1', self.container.inputData.MAPIN1.fullPath ])

      if self.container.inputData.MAPIN2.isSet():
        self.appendCommandLine([ 'MAPIN2', self.container.inputData.MAPIN2.fullPath ])

      if self.container.controlParameters.OUTPUT_PDB_FILE.isSet() :
          self.appendCommandLine ( [ ' XYZOUT EDSTATS-per_atom_metrics.pdb ' ] )

      self.appendCommandLine ( [ ' OUT edstats.out ' ] )

      self.appendCommandScript( "reslo=%s"%( str ( self.container.inputData.RES_LOW ) ))
      self.appendCommandScript( "reshi=%s"%( str ( self.container.inputData.RES_HIGH ) ))

      if self.container.controlParameters.MAIN_AVERAGING.isSet() :
        self.appendCommandScript("main=%s"%(str(self.container.controlParameters.MAIN_AVERAGING)))

      if self.container.controlParameters.SIDE_AVERAGING.isSet() :
        self.appendCommandScript("side=%s"%(str(self.container.controlParameters.SIDE_AVERAGING)))

      if self.container.controlParameters.SCALING :
        self.appendCommandScript("resc=%s"%(str(self.container.controlParameters.SCALING_TYPE)))
      else :
        self.appendCommandScript("resc=none")

      return CPluginScript.SUCCEEDED

    def makeCfftPlugin1 ( self ):
        cfftPlugin1 = self.makePluginObject( 'fft' )
        cfftPlugin1.container.inputData.FPHIIN = self.container.inputData.FPHIIN1
        cfftPlugin1.container.inputData.RESOLUTION = self.container.inputData.RES_HIGH * 3.0 / 4.1
        return cfftPlugin1

# I'd like to keep separate plugins because I plan to do divergent stuff in the future. 3/4.1 seemed to work.

    def makeCfftPlugin2 ( self ):
        cfftPlugin2 = self.makePluginObject( 'fft' )
        cfftPlugin2.container.inputData.FPHIIN = self.container.inputData.FPHIIN2
        cfftPlugin2.container.inputData.RESOLUTION = self.container.inputData.RES_HIGH * 3.0 / 4.1
        return cfftPlugin2
