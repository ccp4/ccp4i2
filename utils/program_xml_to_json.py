import sys
import xml.etree.ElementTree as etree
import json

xmlFile = sys.argv[1]

text = open( xmlFile ).read()
xmlnode = etree.fromstring(text)

validation_rama_outliers = xmlnode.findall('.//Validation/Ramachandran/Outliers/Residue')
mol_rama_outliers = xmlnode.findall('.//Molprobity/Ramachandran_outliers/Outlier')
rota_outliers = xmlnode.findall('.//Molprobity/Rotamer_outliers/Outlier')
cbeta_outliers = xmlnode.findall('.//Molprobity/CBeta_outliers/Outlier')
flips = xmlnode.findall('.//Molprobity/Side_chain_flips/Outlier')
clashes = xmlnode.findall('.//Molprobity/Clashes/Outlier')

json_o = {}
json_o["title"] = "Interesting things from CCP4i2 validate protein"
json_o["items"] = []

for el in validation_rama_outliers:
    chain = el.attrib["chain"]
    resno = el.attrib["seqnum"]
    item = {}
    item["label"] = "Ramachandran Outlier"
    item["position-type"] = "by-residue-spec"
    item["residue-spec"] = [chain,resno,""]
    item["action"] = ["triple-refinement-action","triple-refinement-with-rama-restraints-action"]
    json_o["items"].append(item)
    
for el in mol_rama_outliers:
    chain = el.attrib["chain"]
    resno = el.attrib["seqnum"]
    item = {}
    item["label"] = "Ramachandran Outlier"
    item["position-type"] = "by-residue-spec"
    item["residue-spec"] = [chain,resno,""]
    item["action"] = ["triple-refinement-action","triple-refinement-with-rama-restraints-action"]
    json_o["items"].append(item)

for el in rota_outliers:
    chain = el.attrib["chain"]
    resno = el.attrib["seqnum"]
    item = {}
    item["label"] = "Rotamer outlier"
    item["position-type"] = "by-residue-spec"
    item["residue-spec"] = [chain,resno,""]
    item["action"] = ["auto-fit-rotamer-action"]
    json_o["items"].append(item)

for el in flips:
    chain = el.attrib["chain"]
    resno = el.attrib["seqnum"]
    item = {}
    item["label"] = "Side chain flip"
    item["position-type"] = "by-residue-spec"
    item["residue-spec"] = [chain,resno,""]
    item["action"] = ["side-chain-flip-action"]
    json_o["items"].append(item)

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
    json_o["items"].append(item)

json_string = json.dumps(json_o,indent=4)

print(json_string)
