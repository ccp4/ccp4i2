
import gemmi
import argparse
import sys
import os
import logging
import xml.etree.ElementTree as ET

## the nexr function is copied from:
## adding_stats_to_mmcif/xml_parsing.py

def parse_xml(xml_file):
    """
    checks input file is XML file, parses it.
    :return: root of the xml or None
    """
    root = None
    try:
        if os.path.exists(xml_file):
            tree = ET.parse(xml_file)
            root = tree.getroot()
        else:
            logging.error('xml file: "{}" does not exist'.format(xml_file))
    except Exception as e:
        logging.error(e)

    return root

## the fragment of code from here to the end
## of the aimlessReport class was copied from
## adding_stats_to_mmcif/aimless_xml_parser.py
## with small modifications

# AL: not obvious decision, keeping as it was.
#all = 'Overall'
_all = ''

# AL:  corrections/additions in: pdbx_chi_squared, number_unique_obs,
# meanI_over_sigI_obs, pdbx_chi_squared, percent_possible_all
stats_remap = {
    'reflns': {'pos_list': ['Overall'],
               'cif_to_xml':
                   {'d_resolution_low': 'ResolutionLow',
                    'd_resolution_high': 'ResolutionHigh',
                    'pdbx_Rmerge_I_obs': 'Rmerge' + _all,
                    'pdbx_Rrim_I_all': 'Rmeas' + _all,
                    'pdbx_Rpim_I_all': 'Rpim' + _all,
                    'number_obs': 'NumberReflections',
                    'pdbx_netI_over_sigmaI': 'MeanIoverSD',
                    'pdbx_CC_half': 'CChalf',
                    'pdbx_chi_squared': 'MeanChiSq',
                    'percent_possible_obs': 'Completeness',
                    'pdbx_redundancy': 'Multiplicity'}},
    'reflns_shell': {'pos_list': ['Inner', 'Outer'],
                     'cif_to_xml':
                         {'d_res_low': 'ResolutionLow',
                          'd_res_high': 'ResolutionHigh',
                          'Rmerge_I_obs': 'Rmerge' + _all,
                          'pdbx_Rrim_I_all': 'Rmeas' + _all,
                          'pdbx_Rpim_I_all': 'Rpim' + _all,
                          'number_unique_obs': 'NumberReflections',
                          'meanI_over_sigI_obs': 'MeanIoverSD',
                          'pdbx_CC_half': 'CChalf',
                          'pdbx_chi_squared': 'MeanChiSq',
                          'percent_possible_all': 'Completeness',
                          'pdbx_redundancy': 'Multiplicity'}}
}

def split_on_separator(separator, value_to_split):
    value_list = list()
    logging.debug('separator: "%s"' % separator)
    if separator:
        if separator == ' ':
            value_list = value_to_split.split()
        else:
            value_list = value_to_split.split(separator)
        logging.debug(value_list)
    return value_list

class aimlessReport:
    extra_cif_items = {'pdbx_ordinal': ''}
    table_keys = {'deposition': 'reflns'}

    def __init__(self, xml_file):
        """
        :param xml_file: input aimless XML file
        """
        self.xml_file = xml_file
        self.tree = None
        self.root = None
        self.stats_dict = dict()
        self.cif_item = None
        self.cif_cat = None
        self.number_of_values = 0

    def parse_xml(self):
        """
            checks input file is XML file, parses it.
            prevents parsing twice by checking if self.tree already exists
        :return: True if a parsed aimless XML file, False if not
        """
        if not self.tree:
            self.root = parse_xml(xml_file=self.xml_file)

        if self.root is not None:
            if self.root.tag == 'AIMLESS_PIPE' or self.root.tag == 'AIMLESS':
                logging.debug('is an aimless xml file')
                return True
        return False

    def add_data_to_stats_dict(self, instance, value):
        self.stats_dict.setdefault(self.cif_cat, {}).setdefault(self.cif_item, [''] * self.number_of_values)[
            instance] = value

    def get_data(self):
        """
        :return: statistics dictionary from Dataset XML tag.
        """
        datasetresultnodes = self.root.findall(".//Result/Dataset")
        data_set_counter = 0
        for datasetresultnode in datasetresultnodes:
            data_set_counter += 1
            for self.cif_cat in stats_remap:

                location_list = stats_remap[self.cif_cat]['pos_list']
                self.number_of_values = len(location_list)

                for instance, location in enumerate(location_list):
                    for self.cif_item in stats_remap[self.cif_cat]['cif_to_xml']:
                        logging.debug(self.cif_item)
                        xml_item = stats_remap[self.cif_cat]['cif_to_xml'][self.cif_item]

                        xml_node = datasetresultnode.find(xml_item)
                        xml_item_for_location = xml_node.find(location)

                        logging.debug(xml_item_for_location)
                        xml_value = xml_item_for_location.text.strip()
                        logging.debug(xml_value)

                        self.add_data_to_stats_dict(instance=instance, value=xml_value)

                    for self.cif_item in self.extra_cif_items:
                        if self.extra_cif_items[self.cif_item]:
                            value = self.extra_cif_items[self.cif_item]
                        else:
                            value = instance + 1
                        self.add_data_to_stats_dict(instance=instance, value=str(value))

        return self.stats_dict

    def get_data_from_table(self):
        """
        :return: dictionary of statistics from CCP4 tables.
        """
        ccp4tables = self.root.findall(".//CCP4Table")
        logging.debug(ccp4tables)
        for table in ccp4tables:
            logging.debug(table.attrib)
            if 'id' in table.attrib:
                if table.attrib['id'] in self.table_keys:
                    # need to set cif category based on the table name.
                    self.cif_cat = self.table_keys[table.attrib['id']]
                    headers = table.find('headers')
                    header_list = split_on_separator(separator=headers.attrib['separator'].text,
                                                     value_to_split=headers.text)
                    data = table.find('data').text
                    logging.debug(data)
                    data_lines = data.strip().split('\n')
                    logging.debug(data_lines)
                    number_of_data_values = len(data_lines)
                    logging.debug('number of data items: %s' % number_of_data_values)
                    for instance, d in enumerate(data_lines):
                        d = split_on_separator(separator=headers.attrib['separator'].text, value_to_split=d)
                        for header_pos, self.cif_item in enumerate(header_list):
                            logging.debug('%s - position %s' % (self.cif_item, header_pos))
                            value = d[header_pos]
                            self.add_data_to_stats_dict(instance=instance, value=value)

        return self.stats_dict

    def return_data(self):
        """
        master function which returns data from aimless XML file.
        :return:
        """
        is_aimless_file = self.parse_xml()
        if is_aimless_file:
            self.get_data_from_table()
            if not self.stats_dict:
                self.get_data()

        return self.stats_dict

class DepoDoc(object):
    seq_width = 60

    def __init__(self, cifin, log):
        self.log = log
        self.doc = gemmi.cif.read(cifin)
        self.block = self.doc[0]

    def write_doc(self, cifout):
        self._fix_cross_val()
        self._fix_resolution_limits()
        self._add_exptl_xray()
        wo = gemmi.cif.WriteOptions()
        wo.prefer_pairs = False
        self.doc.write_file(cifout, wo)

    @staticmethod
    def _change_if_beyond(c1, sign, *cc):
        if len(c1) == 1:
            v1 = c1[0]
            for c2 in cc:
                for i, v2 in enumerate(c2):
                    try:
                        if (float(v2) - float(v1))* sign > 0:
                            c2[i] = v1
                    except:
                        pass

    def _fix_resolution_limits(self):
        t = self.block
        a1 = t.find_values('_reflns.d_resolution_low')
        a2 = t.find_values('_reflns_shell.d_res_low')
        r1 = t.find_values('_refine.ls_d_res_low')
        r2 = t.find_values('_refine_ls_shell.d_res_low')
        self._change_if_beyond(a1, +1, a2, r1, r2)
        self._change_if_beyond(r1, +1, r2)
        a1 = t.find_values('_reflns.d_resolution_high')
        a2 = t.find_values('_reflns_shell.d_res_high')
        r1 = t.find_values('_refine.ls_d_res_high')
        r2 = t.find_values('_refine_ls_shell.d_res_high')
        self._change_if_beyond(a1, -1, a2, r1, r2)
        self._change_if_beyond(r1, -1, r2)

    def _fix_cross_val(self):
        t = self.block
        cc = t.get_mmcif_category('_refine')
        pp = dict([(k, v[0]) for k, v in cc.items()])
        p1 = pp['ls_number_reflns_R_free']
        v2 = 'THROUGHOUT'
        if p1 and p1.isdigit():
            v2 = 'FREE R-VALUE'
        pp['pdbx_ls_cross_valid_method'] = v2
        t.set_pairs('_refine.', pp, False)

    def _add_exptl_xray(self):
        block = self.block
        if not block.find_mmcif_category('_exptl.'):
            block.set_mmcif_category('_exptl.',
                dict(entry_id=[block.name], method=['X-RAY DIFFRACTION']))
            self.log.write("### Experimental method record added.\n\n")

    def add_aimless_data(self, xml_file):
        ar = aimlessReport(xml_file=xml_file)
        xml_data = ar.return_data()
        key = 'reflns'
        self.block.set_mmcif_category('_' + key, xml_data[key])
        key = 'reflns_shell'
        self.block.set_mmcif_category('_' + key, xml_data[key])
    '''
    def _update_block(self, st):
        og = gemmi.MmcifOutputGroups(False)
        og.atoms = True
        og.chem_comp = True
        og.conn = True
        og.entity = True
        og.entity_poly = True
        og.entity_poly_seq = True
        og.group_pdb = True
        og.struct_asym = True
        og.auth_all = True
        st.update_mmcif_block(self.block, og)
    '''
    def _update_block(self, st):
        og = gemmi.MmcifOutputGroups(True)
        og.atom_type = False
        og.auth_all = True
        og.exptl = True
        og.refine = False
        og.software = False
        og.symmetry = False
        og.reflns = False
        st.update_mmcif_block(self.block, og)

    def _del_tab(self, lab):
        tab = self.block.find_mmcif_category(lab)
        if tab:
            tab.erase()

    def _del_col(self, lab):
        col = self.block.find_values(lab)
        if col:
            col.erase()

    def _cln_col(self, lab):
        col = self.block.find_values(lab)
        if col:
            for i in range(len(col)):
                col[i] = '.'

    def _del_entity_tables(self):
        self._del_tab('_struct_asym')
        self._del_tab('_entity')
        self._del_tab('_entity_poly')
        self._del_tab('_entity_poly_seq')

    def clear_entities(self, deepclean = False):
        log = self.log
        self._del_entity_tables()
        self._del_col('_atom_site.label_seq_id')
        self._del_col('_atom_site.label_entity_id')
        self._cln_col('_struct_conn.ptnr1_label_seq_id')
        self._cln_col('_struct_conn.ptnr2_label_seq_id')
        st = gemmi.make_structure_from_block(self.block)
        st.assign_het_flags()
        if not deepclean:
            st.add_entity_types(overwrite=True)
            st.setup_entities()
        log.write(
            '### Sequence association was not requested.\n'
            '\n'
            'If the output mmCIF file is used for PDB deposition,\n'
            'all sequence-to-chain associations will need to be\n'
            'made manually during the deposition.\n'
            '\n'
        )
        self._update_block(st)

    def del_h0(self):
        st = gemmi.make_structure_from_block(self.block)
        cou = 0
        for mo in st:
            for ch in mo:
                for res in ch:
                    nat = len(res)
                    for ind, atom in enumerate(reversed(res)):
                        if atom.element.is_hydrogen and atom.occ == 0.0:
                            cou += 1
                            del res[nat - ind - 1]
        if cou > 0:
            self._update_block(st)
        return cou

    def attach_seq(self, seqin, deepclean = False):
        log = self.log
        width = self.seq_width
        success = True
        if deepclean:
            self._del_entity_tables()
        self._del_col('_atom_site.label_seq_id')
        self._del_col('_atom_site.label_entity_id')
        self._cln_col('_struct_conn.ptnr1_label_seq_id')
        self._cln_col('_struct_conn.ptnr2_label_seq_id')
        st = gemmi.make_structure_from_block(self.block)
        st.assign_het_flags()
        st.add_entity_types(overwrite=True)
        st.setup_entities()
        if not deepclean:
            st.clear_sequences()
        with open(seqin) as f:
            seq_lines = f.read().replace('\t', ' ').split('\n')
        for i, line in enumerate(seq_lines):
            if line.startswith('>'):
                seq_lines[i] = '\t' + line + '\t'
        seq_lines = ''.join(seq_lines).split('\t')
        seq_list = [''.join(s.split()) for s in seq_lines[2::2]]
        seq_lines = seq_lines[1::2]
        ind = 0
        print("### Reading target sequences", file=log)
        print(file=log)
        for head, seq in zip(seq_lines, seq_list):
            ind += 1
            print(ind, head, file=log)
            for i in range(0, len(seq), width):
                print(seq[i:i+width], file=log)
            print(file=log)
        st.assign_best_sequences(seq_list)
        st.setup_entities()
        # need example where sequence changes in entity
        # correct alg - simply find what matches
        seq_dict = dict()
        for ind, seq in enumerate(seq_list):
            if not seq in seq_dict:
                seq_dict[seq] = []
            seq_dict[seq].append(ind + 1)
        if len(seq_dict) < len(seq_list):
            log.write("### Duplicate input sequences detected:\n")
            for ind_list in sorted(seq_dict.values()):
                if len(ind_list) > 1:
                    log.write(" %s\n" %", ".join(
                        [str(ind) for ind in ind_list]))
            log.write("In each group, only the first sequence number\n"
                "will be used in subsequent output and warnings.\n")
            log.write("\n")
        ent_dict = dict()
        for ent in st.entities:
            if ent.entity_type is gemmi.EntityType.Polymer:
                sec1lc = gemmi.one_letter_code(ent.full_sequence)
                assert ent.name not in ent_dict
                ind_list = seq_dict.pop(sec1lc, None)
                ent_dict[ent.name] = ind_list[0] if ind_list else 0
        if len(seq_dict) > 0:
            log.write("### Input"
                " sequences not matching any polymer chain:\n")
            log.write(" %s\n" %", ".join([str(ind_list[0])
                for ind_list in sorted(seq_dict.values())]))
            log.write("These sequences will be ignored.\n")
            log.write("\n")
        rejected = []
        table = [("Chain", "Seqence", "Type",
            "Length","Missing", "Secial", "Mismatches", "Status")]
        chains = dict(matches=[], errors=[], missing=[])
        fmt00 = "### Chain %s does not match any input sequence\n"
        fmt01 = "### Chain %s against input sequence %s:\n"
        fmt10 = "Special amino acids: %s\n"
        fmt11 = "Special amino acids: %s (allowed)\n"
        fmt20 = "Missing residues: %s\n"
        fmt21 = "Missing residues: %s (allowed)\n"
        fmt30 = "Mismatches: %s\n"
        fmt31 = "Mismatches: %s (error)\n"
        for chain in st[0]:
            span = chain.get_polymer()
            if span:
                seq2 = span.extract_sequence()
                ent = st.get_entity_of(span)
                assert ent.entity_type is gemmi.EntityType.Polymer
                seq1 = ent.full_sequence
                seq1_ind = ent_dict[ent.name]
                seq1_type = str(ent.polymer_type).split('.')[1]
                if seq1 and seq1_ind:
                    result = gemmi.align_sequence_to_polymer(
                        seq1, span, ent.polymer_type,
                        gemmi.AlignmentScoring())
                    line1 = result.add_gaps(
                       gemmi.one_letter_code(seq1), 1)
                    line2 = result.add_gaps(
                       gemmi.one_letter_code(seq2), 2)
                    match = result.match_string
                    mm = ''.join([line2[i]
                        for i, x in enumerate(match) if x != '|'])
                    n_spec = mm.count('X')
                    n_gaps = mm.count('-')
                    n_mism = len(mm) - n_spec - n_gaps
                    log.write(fmt01 %(chain.name, seq1_ind))
                    log.write((fmt11 if n_spec else fmt10) %n_spec)
                    log.write((fmt21 if n_gaps else fmt20) %n_gaps)
                    log.write((fmt31 if n_mism else fmt30) %n_mism)
                    log.write("\n")
                    for i in range(0, len(match), width):
                        log.write(line1[i:i+width] + "\n")
                        log.write(match[i:i+width] + "\n")
                        log.write(line2[i:i+width] + "\n")
                        log.write("\n")
                    if n_mism:
                        rejected.append(chain.name)
                        chains['errors'].append(chain.name)
                    else:
                        chains['matches'].append(chain.name)
                    table.append((
                        chain.name, str(seq1_ind), seq1_type,
                        str(len(match)), str(n_gaps), str(n_spec),
                        str(n_mism), "error" if n_mism else "valid"))
                else:
                    line2 = gemmi.one_letter_code(seq2)
                    log.write(fmt00 %chain.name)
                    log.write("\n")
                    for i in range(0, len(line2), width):
                        log.write(line2[i:i+width] + "\n")
                    log.write("\n")
                    rejected.append(chain.name)
                    chains['missing'].append(chain.name)
                    table.append((chain.name,
                        '-', seq1_type, '-', '-', '-', '-', "error"))
        def fmt_ch(ch_list, mch=3):
            nch = len(ch_list)
            short_list = ch_list
            if nch > mch:
                short_list = ch_list[:mch] + ['...']
            msg = '0'
            if short_list:
                msg = '%d (%s)' %(nch, ', '.join(short_list))
            return msg
        error_count = len(chains['errors'] + chains['missing'])
        summary = (
            'Chains with valid target sequence associations: ' + 
                str(len(chains['matches'])),
            'Chains with sequence mismatches: ' +
                fmt_ch(chains['errors']),
            'Chains with no associated target sequence: ' +
                fmt_ch(chains['missing'])
        )
        log.write(
            '### Sequence association completed.\n'
            '\n'
            '%s\n'
            '\n'
            'Mismatches or missing associations may complicate\n'
            'the deposition process. Ideally, these should be\n'
            'corrected before deposition.\n'
            '\n'
            'The task may also be run without sequence association.\n'
            'In this case, sequence-to-chain associations will need\n'
            'to be made manually during PDB deposition.\n'
            '\n' %'\n'.join(summary)
        )
        st.assign_label_seq_id(force=True)
        # the next two lines are only for
        # updating _entity_poly_seq.hetero:
        self._update_block(st)
        st = gemmi.make_structure_from_block(self.block)
        self._update_block(st)
        return error_count, summary, table

def main(args=None):
    parser = argparse.ArgumentParser(
        description='''
            Preparing mmcif file for deposition
        ''',
    )
    parser.add_argument(
        'cifin',
        help='''
            input coordinate file in mmcif format
        ''',
        type=str,
        metavar='file-in',
    )
    parser.add_argument(
        'cifout',
        help='''
            output coordinate file in mmcif format
        ''',
        type=str,
        metavar='file-out',
    )
    parser.add_argument(
        '-s', '--seqin',
        help='''
            single- or multi-sequence fasta file
        ''',
        metavar='file-seq',
    )
    parser.add_argument(
        '-a', '--aimless-xml',
        help='''
            xmlout from aimless
        ''',
        metavar='file-seq',
    )
    parser.add_argument(
        '-d', '--deep-clean',
        help='''
            experimenting: removing more or adding less
        ''',
        action='store_true',
    )
    parser.add_argument(
        '-0', '--del-h0',
        help='''
            removing H-atoms with occupancy 0
        ''',
        action='store_true',
    )
    opt = parser.parse_args(args)
    print('-------------------')
    print('- in:')
    print(opt.cifin)
    print(opt.seqin if opt.seqin else 'no sequence')
    print('- out:')
    print(opt.cifout)
    print('-------------------')
    print()
    dd = DepoDoc(opt.cifin, sys.stdout)
    if opt.aimless_xml and os.path.isfile(opt.aimless_xml):
        dd.add_aimless_data(opt.aimless_xml)
    if opt.del_h0:
        print('Deleted H-0occ:', dd.del_h0())
    if opt.seqin:
        error_count, _, table = dd.attach_seq(opt.seqin, opt.deep_clean)
        print('error_count =', error_count)
        for row in table:
            print(' '.join(row))
        print()
    else:
        dd.clear_entities(opt.deep_clean)
    dd.write_doc(opt.cifout)

if __name__ == '__main__':
    main()

