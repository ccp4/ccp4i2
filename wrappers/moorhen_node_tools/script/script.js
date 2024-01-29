
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

const loadCoot = async (codeRoot) => createCootModule({
    print: printToConsole,
    printErr: errorToConsole,
    locateFile: (file) => `${codeRoot}/node_modules/moorhen/baby-gru/wasm/${file}`
})

class CootWrapper {
    constructor(cootModule, PWD) {
        this.cootModule = cootModule
        this.molecules_container = new cootModule.molecules_container_js(false);
        this.PWD = PWD
    }
    async read_pdb(filePath) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(filePath)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const molNo = this.molecules_container.read_pdb(tempFilePath)
        this.cootModule.FS_unlink(tempFilePath)
        return Promise.resolve(molNo)
    }

    async read_mtz(filePath, FLabel, PhiLabel, dunno, isDiff, erm) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(filePath)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const molNo = this.molecules_container.read_mtz(tempFilePath, FLabel, PhiLabel, dunno, isDiff, erm)
        this.cootModule.FS_unlink(tempFilePath)
        return Promise.resolve(molNo)
    }

    async import_cif_dictionary(filePath, iMol) {
        const fileName = filePath.split('/').at(-1)
        const fileContent = await fs.readFile(filePath)
        const tempFilePath = `./${fileName}`
        this.cootModule.FS_createDataFile(".", fileName, fileContent, true, true);
        const dictReadResult = this.molecules_container.import_cif_dictionary(tempFilePath, iMol)
        this.cootModule.FS_unlink(tempFilePath)
        return Promise.resolve(dictReadResult)
    }

    async writeCIFASCII(iMol, filePath) {
        const fileName = filePath.split('/').at(-1)
        const tempFilePath = `./${fileName}`
        this.molecules_container.writeCIFASCII(iMol, tempFilePath)
        const fileContent = this.cootModule.FS.readFile(tempFilePath, { encoding: 'utf8' });
        const bytesWritten = fs.writeFile(filePath, fileContent)
        this.cootModule.FS_unlink(tempFilePath)
        return Promise.resolve(bytesWritten)
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

    fit_ligand(MolHandle_1, MapHandle_1, imol_lig, val1, val2, val3) {
        return this.molecules_container.fit_ligand(MolHandle_1, MapHandle_1, imol_lig, val1, val2, val3)
    }
}

const FIT_LIGAND = async (cootWrapper, args) => {
    try {
        const { XYZIN_0, FPHIIN_0, DICTIN_0, TLC, MAXCOPIES } = args

        mc = cootWrapper.molecules_container
        const iMol = await cootWrapper.read_pdb(XYZIN_0)
        const iMap = await cootWrapper.read_mtz(FPHIIN_0, 'F', 'PHI', '', false, false)
        mc.set_imol_refinement_map(iMap)
        const readDictResult = await cootWrapper.import_cif_dictionary(DICTIN_0, iMol)
        const iDictMol = mc.get_monomer_from_dictionary(TLC, iMol, false)
        const solutions = mc.fit_ligand(iMol, iMap, iDictMol, 1.0, true, 30)
        const nSolutions = solutions.size()
        console.log({nSolutions})
        const solutionMols = []
        for (let iSolution = 0; iSolution < (nSolutions < MAXCOPIES ? nSolutions: MAXCOPIES); iSolution++) {
            solutionMols.push(`${solutions.get(iSolution).imol}`)
        }
        console.log({solutionMols})
        const mergeResult = mc.merge_molecules(iMol, `${solutionMols.join(':')}`)
        console.log({mergeResult})
        const outputFilePath = `${cootWrapper.PWD}/result.pdb`
        console.log({outputFilePath})
        const writeResult = await cootWrapper.writeCIFASCII(iMol, outputFilePath)
        console.log({writeResult})
        return Promise.resolve({ XYZOUT: [{ filePath: outputFilePath, annotation: "Model with ligands fit" }] })
    }
    catch (err) {
        return Promise.reject()
    }
}

const methods = {
    FIT_LIGAND
}

const main = async () => {
    containerContentJson = await fs.readFile(process.argv[2])
    const cootModule = await loadCoot(process.argv[3])
    const PWD = process.argv[4]

    const cootWrapper = new CootWrapper(cootModule, PWD)
    containerContent = JSON.parse(containerContentJson)
    args = {}
    Object.keys(containerContent._value.inputData._value).forEach(inputType => {
        const inputList = containerContent._value.inputData._value[inputType]._value
        inputList.forEach((inputFile, iInputFile) => {
            args[`${inputType}_${iInputFile}`] = inputFile._fullPath
        })
    })
    const methodName = containerContent._value.controlParameters._value.STARTPOINT._value
    const methodParameters = containerContent._value.controlParameters._value[methodName]
    Object.keys(methodParameters._value).forEach(parameterName => {
        args[parameterName] = methodParameters._value[parameterName]._value
    })
    result = await methods[methodName](cootWrapper, args)
    const outputFilePath = `${PWD}/output.json`
    const bytesWritten = fs.writeFile(outputFilePath, JSON.stringify(result))
}

main()