import sys
import xml.etree.ElementTree as etree
import json

def program_xml_to_json(xmlText):

    xmlnode = etree.fromstring(xmlText)

    validation_rama_outliers = xmlnode.findall('.//Validation/Ramachandran/Outliers/Residue')
    mol_rama_outliers = xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier')
    rota_outliers = xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier')
    cbeta_outliers = xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier')
    flips = xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier')
    clashes = xmlnode.findall('.//Molprobity/Clashes/Outlier')

    json_o = {}
    json_o["title"] = "Interesting things from CCP4i2 validate protein"
    json_o["sections"] = []

    if len(validation_rama_outliers)>0 or len(mol_rama_outliers)>0:

        rama_section = {}
        rama_section["title"] = "Ramachandran plot outliers"
        rama_section["items"] = []

        for el in validation_rama_outliers:
            chain = el.attrib["chain"]
            resno = el.attrib["seqnum"]
            item = {}
            item["label"] = "Ramachandran Outlier"
            item["position-type"] = "by-residue-spec"
            item["residue-spec"] = [chain,resno,""]
            item["action"] = ["triple-refinement-action","triple-refinement-with-rama-restraints-action"]
            rama_section["items"].append(item)
    
        for el in mol_rama_outliers:
            chain = el.attrib["chain"]
            resno = el.attrib["seqnum"]
            item = {}
            item["label"] = "Ramachandran Outlier"
            item["position-type"] = "by-residue-spec"
            item["residue-spec"] = [chain,resno,""]
            item["action"] = ["triple-refinement-action","triple-refinement-with-rama-restraints-action"]
            rama_section["items"].append(item)

        json_o["sections"].append(rama_section)

    if len(rota_outliers)>0:

        rota_section = {}
        rota_section["title"] = "Rotamer outliers"
        rota_section["items"] = []

        for el in rota_outliers:
            chain = el.attrib["chain"]
            resno = el.attrib["seqnum"]
            item = {}
            item["label"] = "Rotamer outlier"
            item["position-type"] = "by-residue-spec"
            item["residue-spec"] = [chain,resno,""]
            item["action"] = ["auto-fit-rotamer-action"]
            rota_section["items"].append(item)

        json_o["sections"].append(rota_section)

    if len(flips)>0:

        flip_section = {}
        flip_section["title"] = "Side chain flips"
        flip_section["items"] = []

        for el in flips:
            chain = el.attrib["chain"]
            resno = el.attrib["seqnum"]
            item = {}
            item["label"] = "Side chain flip"
            item["position-type"] = "by-residue-spec"
            item["residue-spec"] = [chain,resno,""]
            item["action"] = ["side-chain-flip-action"]
            flip_section["items"].append(item)

        json_o["sections"].append(flip_section)

    if len(clashes)>0:

        clash_section = {}
        clash_section["title"] = "Molprobity clashes"
        clash_section["items"] = []

        for el in clashes:
            first_atom = el.attrib["first_atom"].split()
            second_atom = el.attrib["second_atom"].split()
            chain_1 = first_atom[0]
            resno_1 = first_atom[1]
            atom_1 = first_atom[3]
            chain_2 = second_atom[0]
            resno_2 = second_atom[1]
            atom_2 = second_atom[3]
            item = {}
            item["label"] = "Molprobity clash"
            item["position-type"] = "by-atom-spec-pair"
            item["atom-1-spec"] = [chain_1,resno_1,"",atom_1,""]
            item["atom-2-spec"] = [chain_2,resno_2,"",atom_2,""]
            item["action"] = ["sphere-refinement-action"]
            clash_section["items"].append(item)

        json_o["sections"].append(clash_section)

    json_string = json.dumps(json_o,indent=4)

    return json_string

if __name__ == "__main__":
    xmlFile = sys.argv[1]
    xmlText = open( xmlFile ).read()
    json_string = program_xml_to_json(xmlText)
    print(json_string)
