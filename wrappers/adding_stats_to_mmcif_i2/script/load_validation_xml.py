"""
Crude script to load wwPDB validation XML file in coot.
Run using the Coot option "Calculate"> "Run Script"

The script will then prompt for the XML filename to be entered using
the Coot console keyboard, with the default of the first *.xml found
in the directory.

Oliver Smart, Anglia Ruskin University 2019
"""
import glob


def load_wwpdb_validation_xml_file(xml_file_name, imol=0):
    """
    in Coot loads wwPDB validation report from filename xml_file_name

    Note:

    Adapted from Coot method pdb_validate
    https://github.com/pemsley/coot/blob/a3b95b34bd007dfdd49c8483df790854beee31c5/python/pdbe_validation_data.py#L851
    """
    try:
        xml_string = open(xml_file_name).read()
        valid_inf = parse_wwpdb_validation_xml(xml_string)
        if valid_inf:
            entry_validation_info = valid_inf[0]
            subgroups = valid_inf[1]
            ss = sort_subgroups(subgroups)
            validation_to_gui(entry_validation_info, ss, imol)
        else:
            message = 'problem with valid_inf'
            print(message)
            add_status_bar_text(message)
    except IOError as e_mess:
        print('load_validate_xml IOERROR: {}'.format(e_mess))


def prompt_for_xml_file():
    """ keyboard prompt for xml file name. returns filename or None"""
    xml_files = glob.glob('*.xml')
    if xml_files:
        xml_first = xml_files[0]
    else:
        xml_first = None
    xml_file_name = raw_input('Enter validation XML file name <{}>: '.
                              format(xml_first))
    if len(xml_file_name) == 0:
        xml_file_name = xml_first
    if os.path.isfile(xml_file_name):
        return xml_file_name
    else:
        return None

if __name__ == '__main__':
    xml_file_name = prompt_for_xml_file()
    if xml_file_name is not None:
        load_wwpdb_validation_xml_file(xml_file_name)
