import { useEffect, useState } from "react"
import { CCP4i2Api } from "./CCP4i2Api"
import { MoorhenDraggableModalBase } from 'moorhen'
import { Modal } from '@mui/material'

import { CCP4i2RunningJobPanel } from './CCP4i2RunningJobPanel'
export const CCP4i2MakeAceDRGLinkPanel = (props) => {
    const [runningJobId, setRunningJobId] = useState(null)
    const [showRunningJobModal, setShowRunningJobModal] = useState(null)

    useEffect(() => {

        if (props.linkData && props.urlRoot && props.projectId) {
            const ccp4i2Api = new CCP4i2Api(props.urlRoot)
            const asyncFunc = async () => {
                console.log(props.linkData)

                let createJobResult = await ccp4i2Api.postCommand("createJob", {
                    projectId: props.projectId,
                    taskName: "MakeLink"
                })
                const newJobId = createJobResult.jobId

                //Upload the selected coordinate set to the new Job
                const molZeros = props.molecules.filter(molecule => molecule.molNo === props.linkData.atomOneFormData.selectedMolNo)
                let pdbData = await molZeros[0].getAtoms()
                const atomsBlob = new Blob([pdbData])
                const uploadResult = await ccp4i2Api.uploadFileForJobObject(newJobId, 'inputData.XYZIN', atomsBlob, `${molZeros[0].name}.pdb`)
                console.log(uploadResult)

                const inputData = {
                    MON_1_TYPE: 'TLC',
                    MON_2_TYPE: 'TLC',
                    RES_NAME_1_TLC: /\(([^)]+)\)/.exec(props.linkData.atomOneFormData.selectedAtom)[1],
                    RES_NAME_2_TLC: /\(([^)]+)\)/.exec(props.linkData.atomTwoFormData.selectedAtom)[1],
                    ATOM_NAME_1: props.linkData.atomOneFormData.selectedAtom.split('/').at(-1),
                    ATOM_NAME_1_TLC: props.linkData.atomOneFormData.selectedAtom.split('/').at(-1),
                    ATOM_NAME_2: props.linkData.atomTwoFormData.selectedAtom.split('/').at(-1),
                    ATOM_NAME_2_TLC: props.linkData.atomTwoFormData.selectedAtom.split('/').at(-1),
                }
                if (props.linkData.atomOneFormData.deleteAtom) {
                    inputData.TOGGLE_DELETE_1 = 'True'
                    inputData.DELETE_1 = props.linkData.atomOneFormData.deleteSelectedAtom.split('/').at(-1)
                    inputData.DELETE_1_LIST = props.linkData.atomOneFormData.deleteSelectedAtom.split('/').at(-1)
                }
                if (props.linkData.atomTwoFormData.deleteAtom) {
                    inputData.TOGGLE_DELETE_2 = 'True'
                    inputData.DELETE_2 = props.linkData.atomTwoFormData.deleteSelectedAtom.split('/').at(-1)
                    inputData.DELETE_2_LIST = props.linkData.atomTwoFormData.deleteSelectedAtom.split('/').at(-1)
                }
                if (props.linkData.atomOneFormData.changeAtomCharge) {
                    inputData.TOGGLE_CHARGE_1 = 'True'
                    inputData.CHARGE_1 = props.linkData.atomOneFormData.changeSelectedAtomCharge.split('/').at(-1)
                    inputData.CHARGE_1_LIST = props.linkData.atomOneFormData.changeSelectedAtomCharge.split('/').at(-1)
                    inputData.CHARGE_1_VALUE = props.linkData.atomOneFormData.newAtomCharge
                }
                if (props.linkData.atomOneFormData.changeBondOrder) {
                    inputData.TOGGLE_CHANGE_1 = 'True'
                    inputData.CHANGE_BOND_1 = props.linkData.atomOneFormData.changeSelectedBondOrder
                    inputData.CHANGE_BOND_1_LIST = props.linkData.atomOneFormData.changeSelectedBondOrder
                    inputData.BOND_1_TYPE = props.linkData.atomOneFormData.newBondOrder
                    inputData.CHANGE_BOND_1_TYPE = props.linkData.atomOneFormData.newBondOrder
                }
                if (props.linkData.atomTwoFormData.changeAtomCharge) {
                    inputData.TOGGLE_CHARGE_2 = 'True'
                    inputData.CHARGE_2 = props.linkData.atomTwoFormData.changeSelectedAtomCharge.split('/').at(-1)
                    inputData.CHARGE_2_LIST = props.linkData.atomTwoFormData.changeSelectedAtomCharge.split('/').at(-1)
                    inputData.CHARGE_2_VALUE = props.linkData.atomTwoFormData.newAtomCharge
                }
                if (props.linkData.atomTwoFormData.changeBondOrder) {
                    inputData.TOGGLE_CHANGE_2 = 'True'
                    inputData.CHANGE_BOND_2 = props.linkData.atomTwoFormData.changeSelectedBondOrder
                    inputData.CHANGE_BOND_2_LIST = props.linkData.atomTwoFormData.changeSelectedBondOrder
                    inputData.BOND_2_TYPE = props.linkData.atomTwoFormData.newBondOrder
                    inputData.CHANGE_BOND_2_TYPE = props.linkData.atomTwoFormData.newBondOrder
                }
                for (let key in inputData) {
                    const valueXml = `<${key}>${inputData[key]}</${key}>`
                    let commandResult = await ccp4i2Api.setJobParameterValue(newJobId, `inputData.${key}`, inputData[key])
                    console.log(commandResult)
                }

                let commandResult = await ccp4i2Api.setJobParameterValue(newJobId, 'controlParameters.TOGGLE_LINK', 'True')
                console.log(commandResult)

                let runJobResult = await ccp4i2Api.postCommand("runJob", {
                    jobId: newJobId
                })
                console.log(runJobResult)

                setRunningJobId(newJobId)
                setShowRunningJobModal(true)

            }
            asyncFunc()
        }
    }, [props.linkData, props.urlRoot, props.projectId])

    return <>
        <MoorhenDraggableModalBase>
            <CCP4i2RunningJobPanel {...props} runningJobId={runningJobId} />
        </MoorhenDraggableModalBase>
    </>
}