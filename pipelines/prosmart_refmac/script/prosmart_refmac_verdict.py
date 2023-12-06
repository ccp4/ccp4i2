from __future__ import print_function

import sys,os

from lxml import etree

import gemmi

sys.path.append(os.path.join(os.environ["CCP4"],"share","jscofe"))
from pycofe.verdicts import verdict_refmac

def getJSCOFERefmac5Verdict(programxml=None,pdbfile=None,refmaclog=None):
    verdict_meta_refmac = verdict_refmac.parseRefmacLog(refmaclog)
    tree = etree.parse(programxml)
    st = gemmi.read_structure(pdbfile)
    st.setup_entities()

    xyz = []
    ligandlen = 0
    macrolen = 0
    waterlen = 0

    for mod in st:
        model = {'model':int(mod.name)}
        chains = []
        for chain in mod:
            if len(chain) > 0 and chain[0].entity_type == gemmi.EntityType.Polymer:
                subchain = chain[0].subchain
                for entity in st.entities:
                    if subchain in entity.subchains:
                        if entity.polymer_type == gemmi.PolymerType.PeptideL or  entity.polymer_type == gemmi.PolymerType.PeptideD:
                            seq = "".join([gemmi.find_tabulated_residue(residue.name).one_letter_code for residue in chain])
                            chains.append({'type':'Protein','id':chain.name,'seq':str(seq),'size':len(chain)})
            for res in chain:
                if res.entity_type == gemmi.EntityType.NonPolymer and (len(res) != 1 or not res[0].element.is_metal):
                    ligandlen += len(res)
                elif res.entity_type == gemmi.EntityType.Polymer:
                    macrolen += len(res)
                elif res.entity_type == gemmi.EntityType.Water:
                    waterlen += len(res)
        model["chains"] = chains
        xyz.append(model)

    totallen = ligandlen + macrolen + waterlen

    resolution = float(tree.xpath("//resolution_high")[0].text)

    verdict_meta = {
        "refmac":verdict_meta_refmac,
        'params': {
            'refmac': verdict_meta_refmac["params"]["refmac"]}, 
        'data': {
            'resolution': resolution}, 
            'xyzmeta': {
                'cryst': {'spaceGroup': st.spacegroup_hm, 'a': st.cell.a, 'b': st.cell.b, 'c': st.cell.c, 'alpha': st.cell.alpha, 'beta': st.cell.beta, 'gamma': st.cell.gamma},
                'xyz': xyz,
                'ligands': []
            }, 
            'molprobity': {
                'clashscore': float(tree.xpath("//Molprobity/Summary/Clashscore")[0].text),
                'molp_score': float(tree.xpath("//Molprobity/Summary/Molprobity_score")[0].text),
                'rms_bonds': float(tree.xpath("//Molprobity/Summary/RMS_bonds")[0].text),
                'rms_angles': float(tree.xpath("//Molprobity/Summary/RMS_angles")[0].text),
                'rama_outliers': float(tree.xpath("//Molprobity/Summary/Ramachandran_outliers")[0].text.replace('%','')),
                'rama_favored': float(tree.xpath("//Molprobity/Summary/Ramachandran_favoured")[0].text.replace('%','')),
                'rota_outliers':  float(tree.xpath("//Molprobity/Summary/Rotamer_outliers")[0].text.replace('%','')),
                #'cbeta_deviations': 1.0,
                #'bfac_overall': 20.9,
                #'bfac_macro': 20.1,
                #'bfac_water': 31.9,
                #'bfac_ligand': 17.2,
                'natoms_overall': totallen,
                'natoms_water': waterlen,
                'natoms_ligand': ligandlen,
                'natoms_macro': macrolen,
            }
    }

    return verdict_refmac.calculate(verdict_meta)

if __name__ == "__main__":
    programxml = "/Users/stuart/CCP4I2_PROJECTS/MGScene/CCP4_JOBS/job_54/program.xml"
    pdbfile = "/Users/stuart/CCP4I2_PROJECTS/MGScene/CCP4_JOBS/job_54/54_mgscene_xyzout_prosmart_refmac.pdb"
    refmaclog = "/Users/stuart/CCP4I2_PROJECTS/MGScene/CCP4_JOBS/job_54/job_1/log.txt"

    result = getJSCOFERefmac5Verdict(programxml=programxml,pdbfile=pdbfile,refmaclog=refmaclog)
    verdict_score = result["score"]
    verdict_message  = result["message"]
    bottomline = result["bottomLine"]
    meanRfree = result["meanRfree"]
    medianClash = result["medianClash"]
    ramaOutliers = result["ramaOutliers"]
    suggestedParameters = result["suggestedParameters"]

    print("verdict_score")
    print(verdict_score)
    print()
    print("verdict_message")
    print(verdict_message)
    print()
    print("bottomline")
    print(bottomline)
    print()
    print("meanRfree")
    print(meanRfree)
    print()
    print("medianClash")
    print(medianClash)
    print()
    print("ramaOutliers")
    print(ramaOutliers)
    print()
    print("suggestedParameters")
    print(suggestedParameters)
    print()
