from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.pipelines.MakeLink.script.link_instruction import (
    MonomerEdits, edit_problems, edit_words, extra_instruction_problems, instruction_words,
)


def _has_atoms(structure):
    """True if the structure holds at least one atom.

    gemmi hands back an empty Structure for a file it cannot make sense of
    instead of raising, so a successful read says nothing about whether there
    is a model in it.
    """
    for model in structure:
        for chain in model:
            for residue in chain:
                for _atom in residue:
                    return True
    return False


class MakeLink(CPluginScript):
    TASKNAME = 'MakeLink'

    # Applying the link to a model is the half of this task that can fail
    # quietly: AceDRG writes the dictionary, and everything after that is our
    # own gemmi code. Every way it can fail to produce the model the user
    # asked for now ends the job as failed, rather than Finished with nothing.
    ERROR_CODES = {
        301: {'description': 'The input model could not be read'},
        302: {'description': 'The input model contains no atoms'},
        303: {'description': 'The generated dictionary does not describe the requested link'},
        304: {'description': 'Failed to apply the link to the input model'},
        305: {'description': 'No residue pair in the model matches the requested link, '
                             'so no link record was added'},
        306: {'description': 'An input model was given but "Apply links to model" is not '
                             'selected, so the model will not be modified',
              'severity': CCP4ErrorHandling.SEVERITY_WARNING},
        307: {'description': '"Apply links to model" is selected but no input model was given'},
        308: {'description': 'The changes to a monomer contradict each other'},
        309: {'description': 'The extra AceDRG instructions would be misread by AceDRG'},
    }

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
        error = super(MakeLink, self).validity()

        # The two ways of asking for a model and not getting one. In the task
        # interface XYZIN is nested inside the TOGGLE_LINK checkbox so neither
        # is reachable, but i2run, a cloned job and the REST API can all set
        # the two independently -- and the job used to finish happily either way.
        if ctrl.TOGGLE_LINK and not inp.XYZIN.isSet():
            error.append(
                klass=self.TASKNAME, code=307,
                details='Select an input model, or turn off "Apply links to model"',
                name=f'{self.TASKNAME}.container.inputData.XYZIN',
                severity=CCP4ErrorHandling.SEVERITY_ERROR)
        if inp.XYZIN.isSet() and not ctrl.TOGGLE_LINK:
            error.append(
                klass=self.TASKNAME, code=306,
                details=self.ERROR_CODES[306]['description'],
                name=f'{self.TASKNAME}.container.controlParameters.TOGGLE_LINK',
                severity=CCP4ErrorHandling.SEVERITY_WARNING)

        # The edits are a description, so a contradiction in one is found
        # here rather than by AceDRG halfway through a run.
        for monomer in (1, 2):
            link_atom = getattr(inp, f'ATOM_NAME_{monomer}')
            for problem in edit_problems(
                    self.monomerEdits(monomer), str(link_atom) if link_atom.isSet() else ''):
                error.append(
                    klass=self.TASKNAME, code=308,
                    details=f'Monomer {monomer}: {problem}',
                    name=f'{self.TASKNAME}.container.inputData.DELETE_ATOMS_{monomer}',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR)
        # A malformed CHANGE or ADD section does not fail in AceDRG: it hangs.
        if ctrl.EXTRA_ACEDRG_INSTRUCTIONS.isSet():
            for problem in extra_instruction_problems(str(ctrl.EXTRA_ACEDRG_INSTRUCTIONS)):
                error.append(
                    klass=self.TASKNAME, code=309, details=problem,
                    name=f'{self.TASKNAME}.container.controlParameters.EXTRA_ACEDRG_INSTRUCTIONS',
                    severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def monomerEdits(self, monomer):
        """The declared edits to monomer 1 or 2, as a MonomerEdits.

        DELETE_ATOMS_n, BOND_ORDERS_n and CHARGES_n are the description. The
        single-edit fields that came before them (TOGGLE_DELETE_n + DELETE_n,
        TOGGLE_CHANGE_n + CHANGE_BOND_n + CHANGE_n_TYPE, TOGGLE_CHARGE_n +
        CHARGE_n + CHARGE_n_VALUE) are folded in, so an older job, a clone of
        one, or an i2run command written for them still asks for the same
        thing. Repeats are harmless; edit_words() writes each edit once.
        """
        inp = self.container.inputData
        n = str(monomer)
        edits = MonomerEdits()
        for name in getattr(inp, 'DELETE_ATOMS_' + n):
            edits.deletes.append(str(name).strip())
        for bond in getattr(inp, 'BOND_ORDERS_' + n):
            edits.bond_orders.append(
                (str(bond.ATOM_1).strip(), str(bond.ATOM_2).strip(), str(bond.ORDER).strip()))
        for charge in getattr(inp, 'CHARGES_' + n):
            value = int(charge.CHARGE) if charge.CHARGE.isSet() else 0
            edits.charges.append((str(charge.ATOM).strip(), value))

        if getattr(inp, 'TOGGLE_DELETE_' + n) and getattr(inp, 'DELETE_' + n).isSet():
            edits.deletes.append(str(getattr(inp, 'DELETE_' + n)).strip())
        if getattr(inp, 'TOGGLE_CHANGE_' + n) and getattr(inp, 'CHANGE_BOND_' + n).isSet():
            atoms = str(getattr(inp, 'CHANGE_BOND_' + n)).split(' -- ')
            order = getattr(inp, 'CHANGE_' + n + '_TYPE')
            edits.bond_orders.append((
                atoms[0].strip(), atoms[1].strip() if len(atoms) == 2 else '',
                str(order).strip() if order.isSet() else ''))
        if getattr(inp, 'TOGGLE_CHARGE_' + n) and getattr(inp, 'CHARGE_' + n).isSet():
            value = getattr(inp, 'CHARGE_' + n + '_VALUE')
            edits.charges.append((
                str(getattr(inp, 'CHARGE_' + n)).strip(), int(value) if value.isSet() else 0))
        return edits

    def normaliseResidueCodes(self):
        """Upper-case (and trim) the monomer-library residue codes.

        AceDRG finds LYS.cif for "Lys" but then looks inside it for a comp
        called "Lys", and fails. The task interface and the monomer-info
        endpoint both upper-case the lookup, so the atom dropdown filled
        normally and the job looked ready to run. Every library code is upper
        case, so normalise once here, before the instruction, the link id and
        the model matching all read the name. CIF-mode names come from the
        user's own dictionary and are left exactly as it spells them.
        """
        for field in (self.container.inputData.RES_NAME_1_TLC,
                      self.container.inputData.RES_NAME_2_TLC):
            if field.isSet():
                code = str(field).strip().upper()
                if code != str(field):
                    field.set(code)

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

       # The edits are written from their declared form, never assembled
       # piecewise: see link_instruction.py for why.
       for monomer in (1, 2):
          edits = self.monomerEdits(monomer)
          link_atom = str(getattr(self.container.inputData, f'ATOM_NAME_{monomer}'))
          problems = edit_problems(edits, link_atom)
          if problems:
             for problem in problems:
                print(f"Error - monomer {monomer}: {problem}")
                self.appendErrorReport(308, f"Monomer {monomer}: {problem}")
             return CPluginScript.FAILED
          words = edit_words(edits, monomer)
          if words:
             instruct += " " + " ".join(words)

       extra = self.container.controlParameters.EXTRA_ACEDRG_INSTRUCTIONS
       if extra.isSet():
          problems = extra_instruction_problems(str(extra))
          if problems:
             for problem in problems:
                print("Error - " + problem)
                self.appendErrorReport(309, problem)
             return CPluginScript.FAILED
          words = instruction_words(str(extra))
          if words:
             instruct += " " + " ".join(words)

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
          # Reported, not just printed: without this the caller silently
          # skipped applying the link and the job still finished successfully.
          print("Error: %s" % e)
          self.appendErrorReport(303, str(e))
       return None
    
    def applyLinksToModel(self,link_bond_value):
       """Add the link record to the input model. Returns a CPluginScript status.

       Nothing here is optional once the user has asked for it: every way this
       can fail to produce the model returns FAILED with a reported error, so
       the job does not finish successfully holding no model.
       """
       if not self.container.controlParameters.TOGGLE_LINK:
          # Not asked for. validity() has already warned if a model was given.
          if self.container.inputData.XYZIN.isSet():
             print("")
             print("An input model was given but 'Apply links to model' is not selected;"
                   " the model will not be modified.")
          return CPluginScript.SUCCEEDED
       if not self.container.inputData.XYZIN.isSet():
          self.appendErrorReport(307, self.ERROR_CODES[307]['description'])
          return CPluginScript.FAILED
       if not link_bond_value:
          # get_link_bond_value() has already reported 303.
          return CPluginScript.FAILED
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
               return 0

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
            return len(link_list)
       
         link_desc = [rname1,aname1,del1,rname2,aname2,del2,link_id,threshold]
         doc_in = None
         modelOut = None
         try: # try to read CIF file
           doc_in = gemmi.cif.read(path)
         except Exception:
           doc_in = None

         links_made = 0
         if doc_in is None: # Input is a PDB file
           try:
             st = gemmi.read_structure(path)
           except Exception as e:
             self.appendErrorReport(301, path+" ("+str(e)+")")
             return CPluginScript.FAILED
           # gemmi returns an EMPTY structure for a file it cannot make sense
           # of rather than raising, so "did it parse" is not the question --
           # "did it contain a model" is. Without this the task wrote out a
           # two-line PDB and called it "Model with links applied".
           if not _has_atoms(st):
             self.appendErrorReport(302, path)
             return CPluginScript.FAILED
           for model in st:
             links_made += apply_links_to_model(st,model,link_desc)
           modelOut = str(self.workDirectory / "ModelWithLinks.pdb")
           st.write_pdb(modelOut,use_linkr=True)
         else: # Input is a CIF file
            doc_out = gemmi.cif.Document()
            seen_atoms = False
            for block in doc_in:
              st = gemmi.make_structure_from_block(block)
              if st:
                if not _has_atoms(st):
                  continue
                seen_atoms = True
                for model in st:
                  links_made += apply_links_to_model(st,model,link_desc)
                doc_out.add_copied_block(st.make_mmcif_document().sole_block())
            if not seen_atoms:
              self.appendErrorReport(302, path)
              return CPluginScript.FAILED
            modelOut = str(self.workDirectory / "ModelWithLinks.cif")
            doc_out.write_file(modelOut)

         # Setting TOGGLE_LINK *and* supplying a model asserts that there is a
         # link here to be found. If none is, that expectation was wrong and
         # the user needs to know which of the residue names, atom names or
         # search distance is at fault -- so this is an error, not a warning.
         # Handing back an unmodified copy of their own model tells them
         # nothing, and is neither of the two outcomes the task promises.
         if links_made == 0:
           self.appendErrorReport(
             305,
             "%s(%s) - %s(%s) within %.3f A" % (rname1,aname1,rname2,aname2,threshold))
           return CPluginScript.FAILED

         self.container.outputData.XYZOUT.setFullPath(modelOut)
         self.container.outputData.XYZOUT.annotation.set(
             'Model with %d link%s applied' % (links_made, '' if links_made == 1 else 's'))
         print("Completed applying links to model: "+modelOut)
         return CPluginScript.SUCCEEDED

       except Exception as e:
         print("Error: %s" % e)
         self.appendErrorReport(304, str(e))
         return CPluginScript.FAILED


    #The startProcess method is where you build in the pipeline logic
    def startProcess(self):
        self.AcedrgLinkPlugins = []
        self.completedPlugins = []

        self.normaliseResidueCodes()
        print("Creating link instruction")
        instruct = self.createLinkInstruction()
        if instruct == CPluginScript.FAILED:
            return CPluginScript.FAILED
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
        from pathlib import Path

        from lxml import etree

        from ccp4i2.core import CCP4Utils
        pipelineXMLStructure = etree.Element("MakeLink")
        linkStatus = CPluginScript.SUCCEEDED
        
        for iPlugin, AcedrgLinkPlugin in enumerate(self.AcedrgLinkPlugins):
            self.container.outputData.CIF_OUT.setFullPath(self.workDirectory / (AcedrgLinkPlugin.container.inputData.LINK_ID.__str__()+"_link.cif"))
            shutil.copyfile(AcedrgLinkPlugin.container.outputData.CIF_OUT.fullPath.__str__(), self.container.outputData.CIF_OUT.fullPath.__str__())
            
            #Create link records, if an input model is provided
            link_bond_value = self.get_link_bond_value(self.container.outputData.CIF_OUT.fullPath.__str__())
            if self.applyLinksToModel(link_bond_value) == CPluginScript.FAILED:
                linkStatus = CPluginScript.FAILED
            
            self.container.outputData.CIF_OUT.annotation.set("Link dictionary: "+self.container.inputData.ANNOTATION.__str__())
            self.container.outputData.UNL_PDB = AcedrgLinkPlugin.container.outputData.UNL_PDB.fullPath.__str__()
            self.container.outputData.UNL_CIF = AcedrgLinkPlugin.container.outputData.UNL_CIF.fullPath.__str__()
            
            #Catenate output XMLs (AcedrgLink may not produce PROGRAMXML)
            programXmlPath = AcedrgLinkPlugin.makeFileName("PROGRAMXML")
            if Path(programXmlPath).exists():
                pluginXMLStructure = CCP4Utils.openFileToEtree(programXmlPath)
                cycleElement = etree.SubElement(pluginXMLStructure,"Cycle")
                cycleElement.text = str(iPlugin)
                pipelineXMLStructure.append(pluginXMLStructure)
        
        with open(self.makeFileName("PROGRAMXML"),"w") as pipelineXMLFile:
            CCP4Utils.writeXML(pipelineXMLFile,etree.tostring(pipelineXMLStructure))
        
        # The dictionary is written either way, but a job asked to update the
        # model and unable to do so is a failed job, not a finished one.
        if linkStatus == CPluginScript.FAILED:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
