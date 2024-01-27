const CootWorker = require('./node_modules/moorhen/baby-gru/wasm/moorhen.js')

const fs = require('fs').promises
const path = require('path')
const createCootModule = require('./node_modules/moorhen/baby-gru/wasm/moorhen.js');
const commander = require('commander')
const printToConsole = (thing) => {
    console.log(thing)
}

const errorToConsole = (thing) => {
    console.error(thing)
}

const loadCoot = async () => createCootModule({
    print: printToConsole,
    printErr: errorToConsole,
    locateFile: (file) => `./node_modules/moorhen/baby-gru/wasm/${file}`
})

class MoleculesContainer {
    constructor(cootModule) {
        this.cootModule = cootModule
        this.molecules_container = new cootModule.molecules_container_js(false);
    }
    async read_pdb(filePath) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(fileName)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const molNo = this.molecules_container.read_pdb(tempFilePath)
        this.cootModule.FS_unlink(tempFilePath)
        return molNo
    }

    async read_mtz(filePath, FLabel, PhiLabel, dunno, isDiff, erm) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(fileName)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const molNo = this.molecules_container.read_mtz(tempFilePath, FLabel, PhiLabel, dunno, isDiff, erm)
        this.cootModule.FS_unlink(tempFilePath)
        return molNo
    }

    async import_cif_dictionary(filePath, iMol) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(fileName)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const molNo = this.molecules_container.import_cif_dictionary(tempFilePath, iMol)
        this.cootModule.FS_unlink(tempFilePath)
        return molNo
    }

    async writeCIFASCII(iMol, filePath) {
        const fileName = filePath.split('/').at(-1)
        const tempFilePath = `./${fileName}`
        this.molecules_container.writeCIFASCII(iMol, tempFilePath)
        const fileContent = this.cootModule.FS.readFile(tempFilePath, { encoding: 'utf8' });
        const bytesWritten = fs.writeFile(filePath, fileContent)
        this.cootModule.FS_unlink(tempFilePath)
        return bytesWritten
    }

    interestingPlaceDataToJSArray(interestingPlaceData) {
        let returnResult = [];

        const interestingPlaceDataSize = interestingPlaceData.size()
        for (let ir = 0; ir < interestingPlaceDataSize; ir++) {
            const residue = interestingPlaceData.get(ir)
            const residueSpec = residue.residue_spec
            returnResult.push({
                modelNumber: residueSpec.model_number,
                chainId: residueSpec.chain_id,
                insCode: residueSpec.ins_code,
                resNum: residueSpec.res_no,
                featureType: residue.feature_type,
                featureValue: residue.feature_value,
                buttonLabel: residue.button_label,
                badness: residue.badness,
                coordX: residue.x,
                coordY: residue.y,
                coordZ: residue.z
            })
            residue.delete()
            residueSpec.delete()
        }
        interestingPlaceData.delete()
        return returnResult
    }

    validationDataToJSArray(validationData, chainID) {
        let returnResult = []
        const cviv = validationData.cviv
        const chainSize = cviv.size()
        for (let chainIndex = 0; chainIndex < chainSize; chainIndex++) {
            const chain = cviv.get(chainIndex)
            if (chainID !== null && chain.chain_id !== chainID) {
                // pass
            } else {
                const resInfo = chain.rviv;
                const resInfoSize = resInfo.size()
                for (let ir = 0; ir < resInfoSize; ir++) {
                    const residue = resInfo.get(ir)
                    const residueSpec = residue.residue_spec
                    returnResult.push({
                        chainId: residueSpec.chain_id,
                        insCode: residueSpec.ins_code,
                        seqNum: residueSpec.res_no,
                        restype: "UNK",
                        value: residue.function_value
                    })
                    residue.delete()
                    residueSpec.delete()
                }
                resInfo.delete()
            }
            chain.delete()
        }
        cviv.delete()
        validationData.delete()
        return returnResult
    }

    set_imol_refinement_map(iMol) {
        return this.molecules_container.set_imol_refinement_map(iMol)
    }

    fit_to_map_by_random_jiggle_using_cid(imol, cid, n_trials, translation_scale_factor) {
        return this.molecules_container.fit_to_map_by_random_jiggle_using_cid(imol, cid, n_trials, translation_scale_factor)
    }

    get_monomer_from_dictionary(tlc, iMol, idealized) {
        return this.molecules_container.get_monomer_from_dictionary(tlc, iMol, idealized)
    }

    merge_molecules(iMol, listOfMolecules) {
        return this.molecules_container.merge_molecules(iMol, listOfMolecules)
    }
    unmodelled_blobs(imol_model, imol_map) {
        const resultAsVector = this.molecules_container.unmodelled_blobs(imol_model, imol_map)
        return this.interestingPlaceDataToJSArray(resultAsVector)
    }

    fit_ligand_right_here(imol_protein, imol_map, imol_ligand, x, y, z, n_rmsd, use_conformers, n_conformers) {
        return this.molecules_container.fit_ligand_right_here(imol_protein, imol_map, imol_ligand, x, y, z, n_rmsd, use_conformers, n_conformers)
    }
}

const executeScript = async () => {
    const cootModule = await loadCoot()
    const molecules_container = new MoleculesContainer(cootModule)
    const iMol = await molecules_container.read_pdb('./t2.pdb')
    const iMap = await molecules_container.read_mtz('./3-prosmart_refmac.mtz', 'F', 'PHI', '', false, false)
    molecules_container.set_imol_refinement_map(iMap)
    const readDictResult = await molecules_container.import_cif_dictionary('./NCL-00025200.cif', iMol)
    const iDictMol = molecules_container.get_monomer_from_dictionary('DRG', iMol, false)
    const blobsResults = molecules_container.unmodelled_blobs(iMol, iMap)
    const modelledBlobs = []
    blobsResults.forEach(blobsResult=>{
        const fitResult = molecules_container.fit_ligand_right_here(iMol, iMap, iDictMol, blobsResult.coordX,
            blobsResult.coordY, blobsResult.coordZ
            , 3, true, 1)
        modelledBlobs.push(fitResult.get(0))
    })
    const mergeResult = molecules_container.merge_molecules(iMol, `${modelledBlobs.join(':')}`)
    molecules_container.writeCIFASCII(iMol, `result.pdb`)
}

console.log(process.arg)
