from ccp4i2.core.CCP4PluginScript import CPluginScript


class MakeLink(CPluginScript):
    TASKNAME = 'MakeLink'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection to PDB file' },
                    203 : {'description' : 'Required input data not set' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                        ['log_mtzjoin.txt', 0]
                       ]

    def __init__(self, *args, **kws):
        super(MakeLink, self).__init__(*args, **kws)
        self.container.inputData.RES_NAME_1_CIF.setQualifier('onlyEnumerators', False)
        self.container.inputData.RES_NAME_2_CIF.setQualifier('onlyEnumerators', False)
        self.container.inputData.ATOM_NAME_1_CIF.setQualifier('onlyEnumerators', False)
        self.container.inputData.ATOM_NAME_2_CIF.setQualifier('onlyEnumerators', False)
        self.container.inputData.ATOM_NAME_1_TLC.setQualifier('onlyEnumerators', False)
        self.container.inputData.ATOM_NAME_2_TLC.setQualifier('onlyEnumerators', False)
        self.container.inputData.DELETE_1_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.DELETE_2_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHARGE_1_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHARGE_2_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHANGE_BOND_1_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHANGE_BOND_2_LIST.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHANGE_BOND_1_TYPE.setQualifier('onlyEnumerators', False)
        self.container.inputData.CHANGE_BOND_2_TYPE.setQualifier('onlyEnumerators', False)
        self.container.controlParameters.MODEL_RES_LIST.setQualifier('onlyEnumerators', False)

    def validity(self):
        """Override to adjust allowUndefined based on MON_TYPE selection.

        MakeLink has conditional field requirements:
        - When MON_1_TYPE='TLC', the TLC fields are required and CIF fields are optional
        - When MON_1_TYPE='CIF', the CIF fields are required and TLC fields are optional
        - LIST fields are populated via GUI dropdowns and are optional for command-line use
        - Toggle-controlled fields (DELETE, CHARGE, CHANGE_BOND) are optional

        This mirrors the logic in MakeLink_gui.py which dynamically sets allowUndefined
        based on user selection, but the def.xml has allowUndefined=False for all.
        """
        inp = self.container.inputData
        ctrl = self.container.controlParameters

        # Determine which mode we're in (defaults to TLC from def.xml)
        mon_1_type = str(inp.MON_1_TYPE) if inp.MON_1_TYPE.isSet() else 'TLC'
        mon_2_type = str(inp.MON_2_TYPE) if inp.MON_2_TYPE.isSet() else 'TLC'

        # Set allowUndefined for fields based on mode
        # For monomer 1
        if mon_1_type == 'TLC':
            # TLC mode: CIF fields are optional
            inp.RES_NAME_1_CIF.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_1_CIF.setQualifier('allowUndefined', True)
        else:
            # CIF mode: TLC fields are optional
            inp.RES_NAME_1_TLC.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_1_TLC.setQualifier('allowUndefined', True)

        # For monomer 2
        if mon_2_type == 'TLC':
            inp.RES_NAME_2_CIF.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_2_CIF.setQualifier('allowUndefined', True)
        else:
            inp.RES_NAME_2_TLC.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_2_TLC.setQualifier('allowUndefined', True)

        # LIST fields are populated by GUI dropdowns - make them optional for CLI use
        inp.DELETE_1_LIST.setQualifier('allowUndefined', True)
        inp.DELETE_2_LIST.setQualifier('allowUndefined', True)
        inp.CHARGE_1_LIST.setQualifier('allowUndefined', True)
        inp.CHARGE_2_LIST.setQualifier('allowUndefined', True)
        inp.CHANGE_BOND_1_LIST.setQualifier('allowUndefined', True)
        inp.CHANGE_BOND_2_LIST.setQualifier('allowUndefined', True)
        ctrl.MODEL_RES_LIST.setQualifier('allowUndefined', True)

        # Now call parent validity() which will use our updated allowUndefined settings
        return super(MakeLink, self).validity()

    def createLinkInstruction(self):
       instruct = "LINK:"

       if not self.container.inputData.ATOM_NAME_1.isSet():
          print("Error - required parameter is not set: ATOM_NAME_1")
          return CPluginScript.FAILED

       if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
          if not self.container.inputData.RES_NAME_1_CIF.isSet():
             print("Error - required parameter is not set: RES_NAME_1_CIF")
             return CPluginScript.FAILED
          if not self.container.inputData.DICT_1.isSet():
             print("Error - required parameter is not set: DICT_1")
             return CPluginScript.FAILED
          instruct += " RES-NAME-1 " + self.container.inputData.RES_NAME_1_CIF.__str__()
          instruct += " ATOM-NAME-1 " + self.container.inputData.ATOM_NAME_1.__str__()
          instruct += " FILE-1 " + self.container.inputData.DICT_1.fullPath.__str__()
       else:
          if not self.container.inputData.RES_NAME_1_TLC.isSet():
             print("Error - required parameter is not set: RES_NAME_1_TLC")
             return CPluginScript.FAILED
          instruct += " RES-NAME-1 " + self.container.inputData.RES_NAME_1_TLC.__str__()
          instruct += " ATOM-NAME-1 " + self.container.inputData.ATOM_NAME_1.__str__()

       if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
          if not self.container.inputData.RES_NAME_2_CIF.isSet():
             print("Error - required parameter is not set: RES_NAME_2_CIF")
             return CPluginScript.FAILED
          if not self.container.inputData.DICT_2.isSet():
             print("Error - required parameter is not set: DICT_2")
             return CPluginScript.FAILED
          instruct += " RES-NAME-2 " + self.container.inputData.RES_NAME_2_CIF.__str__()
          instruct += " ATOM-NAME-2 " + self.container.inputData.ATOM_NAME_2.__str__()
          instruct += " FILE-2 " + self.container.inputData.DICT_2.fullPath.__str__()
       else:
          if not self.container.inputData.RES_NAME_2_TLC.isSet():
             print("Error - required parameter is not set: RES_NAME_2_TLC")
             return CPluginScript.FAILED
          instruct += " RES-NAME-2 " + self.container.inputData.RES_NAME_2_TLC.__str__()
          instruct += " ATOM-NAME-2 " + self.container.inputData.ATOM_NAME_2.__str__()

       if self.container.controlParameters.BOND_ORDER:
          instruct += " BOND-TYPE " + self.container.controlParameters.BOND_ORDER.__str__()

       if self.container.inputData.TOGGLE_DELETE_1:
          if not self.container.inputData.DELETE_1.isSet():
             print("Error - required parameter is not set: DELETE_1")
             return CPluginScript.FAILED
          instruct += " DELETE ATOM " + self.container.inputData.DELETE_1.__str__() + " 1"

       if self.container.inputData.TOGGLE_DELETE_2:
          if not self.container.inputData.DELETE_2.isSet():
             print("Error - required parameter is not set: DELETE_2")
             return CPluginScript.FAILED
          instruct += " DELETE ATOM " + self.container.inputData.DELETE_2.__str__() + " 2"

       if self.container.inputData.TOGGLE_CHANGE_1:
          if not self.container.inputData.CHANGE_BOND_1.isSet():
             print("Error - required parameter is not set: CHANGE_BOND_1")
             return CPluginScript.FAILED
          if not self.container.inputData.CHANGE_1_TYPE.isSet():
             print("Error - required parameter is not set: CHANGE_1_TYPE")
             return CPluginScript.FAILED
          atoms = self.container.inputData.CHANGE_BOND_1.__str__().split(" -- ")
          if len(atoms) != 2:
             print("Error interpreting bond: "+self.container.inputData.CHANGE_BOND_1.__str__()+" : "+str(bonds))
             return CPluginScript.FAILED
          instruct += " CHANGE BOND " + atoms[0] + " " + atoms[1] + " " + self.container.inputData.CHANGE_1_TYPE.__str__() + " 1"

       if self.container.inputData.TOGGLE_CHANGE_2:
          if not self.container.inputData.CHANGE_BOND_2.isSet():
             print("Error - required parameter is not set: CHANGE_BOND_2")
             return CPluginScript.FAILED
          if not self.container.inputData.CHANGE_2_TYPE.isSet():
             print("Error - required parameter is not set: CHANGE_2_TYPE")
             return CPluginScript.FAILED
          atoms = self.container.inputData.CHANGE_BOND_2.__str__().split(" -- ")
          if len(atoms) != 2:
             print("Error interpreting bond: "+self.container.inputData.CHANGE_BOND_2.__str__()+" : "+str(bonds))
             return CPluginScript.FAILED
          instruct += " CHANGE BOND " + atoms[0] + " " + atoms[1] + " " + self.container.inputData.CHANGE_2_TYPE.__str__() + " 2"

       #MN Patch to handle changes in formal charge
       if self.container.inputData.TOGGLE_CHARGE_1:
          if not self.container.inputData.CHARGE_1.isSet():
             print("Error - required parameter is not set: CHARGE_1")
             return CPluginScript.FAILED
          if not self.container.inputData.CHARGE_1_VALUE.isSet():
             print("Error - required parameter is not set: CHARGE_1_VALUE")
             return CPluginScript.FAILED
          instruct += " CHANGE CHARGE 1 " + self.container.inputData.CHARGE_1.__str__() + " " + self.container.inputData.CHARGE_1_VALUE.__str__()

       if self.container.inputData.TOGGLE_CHARGE_2:
          if not self.container.inputData.CHARGE_2.isSet():
             print("Error - required parameter is not set: CHARGE_2")
             return CPluginScript.FAILED
          if not self.container.inputData.CHARGE_2_VALUE.isSet():
             print("Error - required parameter is not set: CHARGE_2_VALUE")
             return CPluginScript.FAILED
          instruct += " CHANGE CHARGE 2 " + self.container.inputData.CHARGE_2.__str__() + " " + self.container.inputData.CHARGE_2_VALUE.__str__()
       ##END MN CHARGE PATCH##


       if self.container.controlParameters.EXTRA_ACEDRG_INSTRUCTIONS.isSet():
            for kwLine in str(self.container.controlParameters.EXTRA_ACEDRG_INSTRUCTIONS).split('\n'):
                kw = kwLine.lstrip().rstrip()
                #print 'kw','['+str(kw)+']'
                if len(kw)>0:
                   if str(kw)[0] != '#':
                      instruct += " " + kw

       return instruct
    
    def createLinkInstructionFile(self,instruct):
       instructFile = self.workDirectory / "link_instruction.txt"
       with instructFile.open("w") as file:
          file.write(instruct)
       self.container.outputData.INSTRUCTION_FILE.setFullPath(instructFile)
       self.container.outputData.INSTRUCTION_FILE.annotation.set('AceDRG instruction file')
    
    def get_link_bond_value(self,cif_file_path):
       print('')
       print("Getting link bond value from dictionary...")
       try:
          from gemmi import cif
          link_dict = cif.read_file(cif_file_path)
          block = link_dict.find_block("link_list")
          if block:
             if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
                rname1 = self.container.inputData.RES_NAME_1_CIF.__str__()
             else:
                rname1 = self.container.inputData.RES_NAME_1_TLC.__str__()
             if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
                rname2 = self.container.inputData.RES_NAME_2_CIF.__str__()
             else:
                rname2 = self.container.inputData.RES_NAME_2_TLC.__str__()
             aname1 = self.container.inputData.ATOM_NAME_1.__str__()
             aname2 = self.container.inputData.ATOM_NAME_2.__str__()
             
             chem_link = block.find('_chem_link.',['id','comp_id_1','comp_id_2'])
             link_ids = []
             for link in chem_link:
                if link[1] == rname1 and link[2] == rname2:
                   link_ids.append(link[0])

             if len(link_ids) == 0:
                raise Exception("Cannot find correct link in dictionary: "+cif_file_path)

             if len(link_ids) > 1:
                print("Warning - multiple link descriptions found between the same residues in the dictionary. Unexpected behaviour may be encountered. Continuing anyway...")

             print("Searching for link between atoms "+aname1+" and "+aname2)
             bond_values = []
             for link_id in link_ids:
                print("Found link ID: "+link_id)
                link_block_id = "link_"+link_id
                link_block = link_dict.find_block(link_block_id)
                if link_block:
                   chem_link_bond = link_block.find('_chem_link_bond.',['link_id','atom_id_1','atom_id_2','value_dist'])
                   for link_bond in chem_link_bond:
                      print("Found link description: "+link_bond[0]+" between atoms "+link_bond[1]+" and "+link_bond[2])
                      if link_bond[0] == link_id and link_bond[1] == aname1 and link_bond[2] == aname2:
                         bond_values.append(link_bond[3])
                         print("Found description for link between "+aname1+" and "+aname2+". Ideal bond value: "+link_bond[3])
              
             if len(bond_values) == 0:
                raise Exception("Cannot find correct link in dictionary: "+cif_file_path)

             if len(bond_values) > 1:
                raise Exception("Multiple matching link descriptions found in dictionary: "+cif_file_path)
             
             return bond_values[0]

          else:
             raise Exception("Cannot find link_list block in dictionary: "+cif_file_path)
       except Exception as e:
          print("Error: %s" % e)
          print("Cannot continue.")
          print("")
       return None
    
    def applyLinksToModel(self,link_bond_value):
       if not link_bond_value:
          return
       if not self.container.controlParameters.TOGGLE_LINK:
          return
       if not self.container.inputData.XYZIN.isSet():
          return
       path = self.container.inputData.XYZIN.fullPath.__str__().rstrip()

       threshold = 0.0
       if self.container.controlParameters.LINK_DISTANCE.isSet():
         threshold = self.container.controlParameters.LINK_DISTANCE * float(link_bond_value)

       if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
          rname1 = self.container.inputData.RES_NAME_1_CIF.__str__()
       else:
          rname1 = self.container.inputData.RES_NAME_1_TLC.__str__()
       if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
          rname2 = self.container.inputData.RES_NAME_2_CIF.__str__()
       else:
          rname2 = self.container.inputData.RES_NAME_2_TLC.__str__()
       aname1 = self.container.inputData.ATOM_NAME_1.__str__()
       aname2 = self.container.inputData.ATOM_NAME_2.__str__()
       link_id = self.container.inputData.LINK_ID.__str__()

       del1 = None
       del2 = None
       if self.container.inputData.TOGGLE_DELETE_1:
          if self.container.inputData.DELETE_1.isSet():
             del1 = self.container.inputData.DELETE_1.__str__()
       if self.container.inputData.TOGGLE_DELETE_2:
          if self.container.inputData.DELETE_2.isSet():
             del2 = self.container.inputData.DELETE_2.__str__()

       print('')
       print("Applying links to model...")
       print("Using detection threshold: "+str(threshold)+" Angstroms")
       try:
         import gemmi
         
         def create_link(conn_list,a1,a2,linkid,ASU):
           con = gemmi.Connection()
           ctr = 1
           con.name = 'link'+str(ctr)
           con_names = [conn.name for conn in conn_list]
           while con.name in con_names:
             ctr += 1
             con.name = 'link'+str(ctr)
           con.type = gemmi.ConnectionType.Covale
           con.partner1 = a1
           con.partner2 = a2
           con.link_id = linkid
           con.asu = ASU
           return con

         def apply_links_to_model(st,model,link_desc):
            res1,atom1,del1,res2,atom2,del2,linkid,max_dist = link_desc
            
            atom1_list = []
            atom2_list = []
            for chain in model:
             for residue in chain:
               if residue.name == res1:
                 for atom in residue:
                   if atom.name == atom1:
                     atom1_list.append([atom,gemmi.AtomAddress(chain.name,residue.seqid,residue.name,atom.name)])
               if residue.name == res2:
                 for atom in residue:
                   if atom.name == atom2:
                     atom2_list.append([atom,gemmi.AtomAddress(chain.name,residue.seqid,residue.name,atom.name)])

            link_list = []
            for a1,addr1 in atom1_list:
             for a2,addr2 in atom2_list:
               FOUND_ASU = None
               if st.cell.find_nearest_image(a1.pos, a2.pos, gemmi.Asu.Same).dist() < max_dist:
                 FOUND_ASU = gemmi.Asu.Same
               elif st.cell.find_nearest_image(a1.pos, a2.pos, gemmi.Asu.Different).dist() < max_dist:
                 FOUND_ASU = gemmi.Asu.Different
               if FOUND_ASU:
                 link_list.append([addr1,addr2,FOUND_ASU])

            if len(link_list) == 0:
               print("No matching links found - no links will be added to the model")
               return

            for a1,a2,asu in link_list:
                st.connections.append(create_link(st.connections,a1,a2,linkid,asu))
                cra1 = model.find_cra(a1)
                cra2 = model.find_cra(a2)
                print("Created link: "+" ".join([str(st.connections[-1].link_id),str(st.connections[-1].name),str(cra1),'-',str(cra2),str(asu)]))
                remove_atoms = []
                if del1:
                  for atom in cra1.residue:
                    if atom.name == del1:
                      remove_atoms.append([cra1.residue,atom])
                if del2:
                  for atom in cra2.residue:
                    if atom.name == del2:
                      remove_atoms.append([cra2.residue,atom])
                for res,atom in remove_atoms:
                  print("Removed atom: "+atom.name+" from residue: "+str(res))
                  res.remove_atom(atom.name,atom.altloc)
       
         link_desc = [rname1,aname1,del1,rname2,aname2,del2,link_id,threshold]
         doc_in = None
         modelOut = None
         try: # try to read CIF file
           doc_in = gemmi.cif.read(path)
         except:
           try: # try to read PDB file
             st = gemmi.read_structure(path)
             if st:
               for model in st:
                 apply_links_to_model(st,model,link_desc)
               modelOut = str(self.workDirectory / "ModelWithLinks.pdb")
               st.write_pdb(modelOut,use_linkr=True)
             else:
               raise Exception("Cannot read input model: "+path)
           except:
             raise Exception("Cannot read input model: "+path)

         if doc_in: # Input is CIF file:
            doc_out = gemmi.cif.Document()
            for block in doc_in:
              st = gemmi.make_structure_from_block(block)
              if st:
                for model in st:
                  apply_links_to_model(st,model,link_desc)
                doc_out.add_copied_block(st.make_mmcif_document().sole_block())
            modelOut = str(self.workDirectory / "ModelWithLinks.cif")
            doc_out.write_file(modelOut)
   
         self.container.outputData.XYZOUT.setFullPath(modelOut)
         self.container.outputData.XYZOUT.annotation.set('Model with links applied')
         print("Completed applying links to model: "+modelOut)
                
       except Exception as e:
         print("Error: %s" % e)
         print("Cannot continue.")
         print("")


    #The startProcess method is where you build in the pipeline logic
    def startProcess(self):
        self.AcedrgLinkPlugins = []
        self.completedPlugins = []

        print("Creating link instruction")
        instruct = self.createLinkInstruction()
        print(instruct)
        self.createLinkInstructionFile(instruct)
        print("Written link instruction file to: ",self.container.outputData.INSTRUCTION_FILE.fullPath.__str__())

        self.AcedrgLinkPlugins.append(self.makePluginObject("AcedrgLink"))
        self.AcedrgLinkPlugins[-1].container.inputData.INSTRUCTION_FILE = self.container.outputData.INSTRUCTION_FILE
        
        link_id = ""
        if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
           link_id += self.container.inputData.RES_NAME_1_CIF.__str__()
        else:
           link_id += self.container.inputData.RES_NAME_1_TLC.__str__()
        link_id += "-"
        if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
           link_id += self.container.inputData.RES_NAME_2_CIF.__str__()
        else:
           link_id += self.container.inputData.RES_NAME_2_TLC.__str__()
        self.container.inputData.LINK_ID.set(link_id)
        self.AcedrgLinkPlugins[-1].container.inputData.LINK_ID.set(link_id)
        
        annotation = ""
        if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
           annotation += self.container.inputData.RES_NAME_1_CIF.__str__()
        else:
           annotation += self.container.inputData.RES_NAME_1_TLC.__str__()
        annotation += "("+self.container.inputData.ATOM_NAME_1.__str__()+")"
        annotation += " - "
        if self.container.inputData.MON_2_TYPE.__str__() == 'CIF':
           annotation += self.container.inputData.RES_NAME_2_CIF.__str__()
        else:
           annotation += self.container.inputData.RES_NAME_2_TLC.__str__()
        annotation += "("+self.container.inputData.ATOM_NAME_2.__str__()+")"
        self.container.inputData.ANNOTATION.set(annotation)
        self.AcedrgLinkPlugins[-1].container.inputData.ANNOTATION.set(annotation)
        
        if self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS.isSet():
           self.AcedrgLinkPlugins[-1].container.controlParameters.EXTRA_ACEDRG_KEYWORDS = self.container.controlParameters.EXTRA_ACEDRG_KEYWORDS

#        self.applyLinksToModel(1.5) # this is just for testing - it's quicker to apply links before running AceDRG, though really it should de done after running AceDRG.
#        return CPluginScript.FAILED
        
        AcedrgLinkResult = self.AcedrgLinkPlugins[-1].process()

        return CPluginScript.SUCCEEDED

    #This method will be called as each plugin completes if the pipeline is run asynchronously
    def pluginFinished(self, whichPlugin):
        self.completedPlugins.append(whichPlugin)
        if len(self.AcedrgLinkPlugins) == len(self.completedPlugins):
            postProcessStaus = super(MakeLink, self).postProcess(processId=self._runningProcessId)
            self.reportStatus(postProcessStatus)
            
    def processOutputFiles(self):
        #Create (dummy) PROGRAMXML
        import shutil

        from lxml import etree

        from ccp4i2.core import CCP4Utils
        pipelineXMLStructure = etree.Element("MakeLink")
        
        for iPlugin, AcedrgLinkPlugin in enumerate(self.AcedrgLinkPlugins):
            self.container.outputData.CIF_OUT.setFullPath(self.workDirectory / (AcedrgLinkPlugin.container.inputData.LINK_ID.__str__()+"_link.cif"))
            shutil.copyfile(AcedrgLinkPlugin.container.outputData.CIF_OUT.fullPath.__str__(), self.container.outputData.CIF_OUT.fullPath.__str__())
            
            #Create link records, if an input model is provided
            link_bond_value = self.get_link_bond_value(self.container.outputData.CIF_OUT.fullPath.__str__())
            self.applyLinksToModel(link_bond_value)
            
            self.container.outputData.CIF_OUT.annotation.set("Link dictionary: "+self.container.inputData.ANNOTATION.__str__())
            self.container.outputData.UNL_PDB = AcedrgLinkPlugin.container.outputData.UNL_PDB.fullPath.__str__()
            self.container.outputData.UNL_CIF = AcedrgLinkPlugin.container.outputData.UNL_CIF.fullPath.__str__()
            
            #Catenate output XMLs
            pluginXMLStructure = CCP4Utils.openFileToEtree(AcedrgLinkPlugin.makeFileName("PROGRAMXML"))
            cycleElement = etree.SubElement(pluginXMLStructure,"Cycle")
            cycleElement.text = str(iPlugin)
            pipelineXMLStructure.append(pluginXMLStructure)
        
        with open(self.makeFileName("PROGRAMXML"),"w") as pipelineXMLFile:
            CCP4Utils.writeXML(pipelineXMLFile,etree.tostring(pipelineXMLStructure))
        
        return CPluginScript.SUCCEEDED
