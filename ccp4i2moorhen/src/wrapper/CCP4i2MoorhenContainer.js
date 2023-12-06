import { useRef, useState, useReducer } from 'react';
import { useEffect } from 'react';
import { MoorhenContainer, MoorhenMolecule, MoorhenMap, MoorhenMoleculeSelect, addMolecule, MoorhenDraggableModalBase, addMap, setActiveMap } from 'moorhen';
import { InputGroup, Modal, NavDropdown } from 'react-bootstrap';
import { Avatar, Button, Link, MenuItem, Typography } from '@mui/material';
import $ from 'jquery';
import { CCP4i2RefmacPanel } from './CCP4i2RefmacPanel';
import { CCP4i2DataBrowser } from './CCP4i2DataBrowser';
import { CCP4i2MakeAceDRGLinkPanel } from './CCP4i2MakeAceDRGLinkPanel';
import { jobFields } from './CCP4i2JobWidgets';
import { useDispatch, useSelector } from 'react-redux'

const saveMolNo = async (urlRoot, molNo, cootJobId, iSave, setISave, molecules, setToastContent) => {
    const formData = new FormData()
    formData.append('jobId', cootJobId)
    formData.append('fileRoot', `COOT_FILE_DROP/output_${iSave.toString().padStart(3, '0')}`)
    formData.append('fileExtension', ".pdb")
    const molZeros = molecules.filter(molecule => molecule.molNo === molNo)
    if (molZeros.length === 1) {
        let pdbData = await molZeros[0].getAtoms()
        const atomsBlob = new Blob([pdbData])
        formData.append('file', atomsBlob)
        fetch(`${urlRoot}/database/uploadFileToJob`, {
            method: "POST",
            body: formData
        })
            .then(response => response.json())
            .then(result => {
                setISave(iSave + 1)
                setToastContent("File saved")
            })
    }
    else {
        setToastContent(`No mols with index ${molNo}`)
    }
}

export const CCP4i2MoorhenContainer = (props) => {
    const glRef = useRef(null)
    const [appTitle, setAppTitle] = useState('CCP4i2Moorhen')
    const moleculesRef = useRef(null)
    const molecules = useSelector((state) => state.molecules)
    const mapsRef = useRef(null)
    const maps = useSelector((state) => state.maps)
    const [currentDropdownId, setCurrentDropdownId] = useState(-1)

    const [showModal, setShowModal] = useState(false)
    const [iSave, setISave] = useState(0)
    const [showMoleculeSelectorModal, setShowMoleculeSelectorModal] = useState(false)
    const [toastContent, setToastContent] = useState("")
    const [showRefmacModal, setShowRefmacModal] = useState(false)
    const [showAceDRGLinkModal, setShowAceDRGLinkModal] = useState(false)
    const [linkData, setLinkData] = useState(null)

    const selectorRef = useRef(null)

    const [currentJob, setCurrentJob] = useState(null)
    const { urlRoot } = props;

    const ccp4i2MenuRef = useRef(null)

    const dispatch = useDispatch()

    const commandCentre = useRef(null)

    const cootInitialized = useSelector(state => state.generalStates.cootInitialized)

    const aceDRGInstance = {
        createCovalentLink: (atomOneFormData, atomTwoFormData) => {
            setLinkData({ atomOneFormData, atomTwoFormData })
            setShowAceDRGLinkModal(true)
        }
    }

    const collectedProps = {
        glRef, commandCentre, currentDropdownId, setCurrentDropdownId,
        showModal, setShowModal, showMoleculeSelectorModal, setShowMoleculeSelectorModal, iSave, setISave,
        toastContent, setToastContent, urlRoot, aceDRGInstance, moleculesRef, mapsRef, appTitle
    }

    const controls = useRef(null);
    //FIXME: hardwired

    useEffect(() => {
        if (props.cootJob) {
            const asyncFunc = async () => {
                const jobResult = await fetch(`${props.urlRoot}/ModelValues?${$.param({
                    __type__: "Jobs",
                    __values__: jobFields,
                })}`).then(response => response.json())
                const theJob = jobResult.results[0]
                setCurrentJob(theJob)
            }
            asyncFunc()
        }
    }, [props.cootJob])

    const makeDbFilePromise = async (fileOfType, mimeType, arg) => {
        const fileDict = {
            filetypeid__filetypename: mimeType,
        }
        const dbFileIdNodes = $(fileOfType).find("dbFileId")
        if (dbFileIdNodes.length === 1 && $(dbFileIdNodes[0]).text().length > 0) {
            fileDict['fileid'] = $(dbFileIdNodes[0]).text()
        }
        const subTypeNodes = $(fileOfType).find("subType")
        if (subTypeNodes.length === 1 && $(subTypeNodes[0]).text().length > 0) {
            fileDict['filesubtype'] = parseInt($(subTypeNodes[0]).text())
        }
        const anotationNodes = $(fileOfType).find("annotation")
        if (anotationNodes.length === 1 && $(anotationNodes[0]).text().length > 0) {
            fileDict['annotation'] = $(anotationNodes[0]).text()
        }
        const fileUrl = `${props.urlRoot}/getFileWithPredicate?${$.param({ fileid: fileDict.fileid })}`
        const annotation = fileDict.annotation ? fileDict.annotation : "A file"
        const subType = fileDict.filesubtype ? fileDict.filesubtype : 1
        return makeFilePromise(fileUrl, mimeType, annotation, subType, arg)
    }

    const makeParamsFilePromise = async (fileOfType, mimeType, arg) => {
        const fileDict = {}
        const projectNodes = $(fileOfType).find("project")
        if (projectNodes.length === 1 && $(projectNodes[0]).text().length > 0) {
            fileDict['projectId'] = $(projectNodes[0]).text()
        }
        const relPathNodes = $(fileOfType).find("relPath")
        if (relPathNodes.length === 1 && $(relPathNodes[0]).text().length > 0) {
            fileDict['relPath'] = $(relPathNodes[0]).text()
        }
        const baseNameNodes = $(fileOfType).find("baseName")
        if (baseNameNodes.length === 1 && $(baseNameNodes[0]).text().length > 0) {
            fileDict['baseName'] = $(baseNameNodes[0]).text()
        }
        let subType = 1
        const subTypeNodes = $(fileOfType).find("subType")
        if (subTypeNodes.length === 1 && $(subTypeNodes[0]).text().length > 0) {
            subType = $(subTypeNodes[0]).text()
        }
        let annotation = "A file"
        const annotationNodes = $(fileOfType).find("annotation")
        if (annotationNodes.length === 1 && $(annotationNodes[0]).text().length > 0) {
            annotation = $(annotationNodes[0]).text()
        }
        const fileUrl = `${props.urlRoot}/getProjectFileData?${$.param(fileDict)}`
        return makeFilePromise(fileUrl, mimeType, annotation, subType, arg)
    }

    const makeFilePromise = async (fileUrl, mimeType, annotation, subType, arg) => {
        if (mimeType === 'application/refmac-dictionary') {
            return fetch(fileUrl)
                .then(response => response.text())
                .then(async fileContent => {
                    if (fileContent.trim().length > 0) {
                        //Associate with molNos if provide in arg, or place in general dictionary
                        const molNos = arg?.molNos?.length > 0 ? arg.molNos : [-999999]
                        const readDictPromises = molNos.map(molNo => commandCentre.current.cootCommand({
                            returnType: "status",
                            command: 'shim_read_dictionary',
                            commandArgs: [fileContent, molNo]
                        }, true))
                        return Promise.all(readDictPromises)
                            .then(async readReturns => {
                                molecules
                                    .filter(molecule => molNos.includes(molecule.molNo) || molNos[0] === -999999)
                                    .forEach(async molecule => {
                                        molecule.setAtomsDirty(true)
                                        await molecule.redraw()
                                    })
                                return Promise.resolve(true)
                            })
                    }
                    else {
                        return Promise.resolve(true)
                    }
                })
        }
        else if (mimeType === 'chemical/x-pdb') {
            const newMolecule = new MoorhenMolecule(commandCentre, glRef)
            return newMolecule.loadToCootFromURL(fileUrl, annotation)
                .then(async reply => {
                    dispatch(addMolecule(reply))
                    newMolecule.centreOn()
                    await reply.fetchIfDirtyAndDraw('CBs')
                    return Promise.resolve(newMolecule.molNo)
                })
        }
        else if (mimeType === 'application/CCP4-mtz-map') {
            const newMap = new MoorhenMap(commandCentre, glRef)
            const selectedColumns = {
                F: 'F', PHI: 'PHI', W: "", useWeight: false,
                isDifference: [2, 3, 4].includes(subType)
            }
            return newMap.loadToCootFromMtzURL(fileUrl,
                annotation,
                selectedColumns
            )
                .then(reply => {
                    if (subType == 1) {
                        dispatch(setActiveMap(reply))
                    }
                    dispatch(addMap(reply))
                    return Promise.resolve(reply)
                })
        }
    }

    const handleCootJob = () => {
        const arg = { molNos: [] }
        fetch(`${props.urlRoot}/getJobFile?jobId=${props.cootJob}&fileName=input_params.xml`)
            .then(response => response.text())
            .then(text => { return Promise.resolve($.parseXML(text)) })
            .then(xmlStruct => {
                const inputData = $(xmlStruct).find("inputData")
                const filesOfType = inputData.find("CPdbDataFile")
                let readFilePromises = []
                filesOfType.toArray().forEach(fileOfType => {
                    const dbFileIdNodes = $(fileOfType).find("dbFileId")
                    const projectNodes = $(fileOfType).find("project")
                    if (dbFileIdNodes.length === 1 && $(dbFileIdNodes[0]).text().length > 0) {
                        readFilePromises.push(makeDbFilePromise(fileOfType, "chemical/x-pdb"))
                    }
                    else if (projectNodes.length === 1 && $(projectNodes[0]).text().length > 0) {
                        readFilePromises.push(makeParamsFilePromise(fileOfType, "chemical/x-pdb"))
                    }
                })
                return Promise.all(readFilePromises)
                    .then(readMolNos => {
                        readMolNos.forEach(molNo => { arg.molNos.push(molNo) })
                        const filesOfType = inputData.find("DICT")
                        let readFilePromises = []
                        filesOfType.toArray().forEach(fileOfType => {
                            const dbFileIdNodes = $(fileOfType).find("dbFileId")
                            const projectNodes = $(fileOfType).find("project")
                            if (dbFileIdNodes.length === 1 && $(dbFileIdNodes[0]).text().length > 0) {
                                readFilePromises.push(makeDbFilePromise(fileOfType, "application/refmac-dictionary", arg))
                            }
                            else if (projectNodes.length === 1 && $(projectNodes[0]).text().length > 0) {
                                readFilePromises.push(makeParamsFilePromise(fileOfType, "application/refmac-dictionary", arg))
                            }
                        })
                        return Promise.all(readFilePromises)
                    })
                    .then(() => {
                        const filesOfType = inputData.find("CMapCoeffsDataFile")
                        let readFilePromises = []
                        filesOfType.toArray().forEach(fileOfType => {
                            const dbFileIdNodes = $(fileOfType).find("dbFileId")
                            const projectNodes = $(fileOfType).find("project")
                            if (dbFileIdNodes.length === 1 && $(dbFileIdNodes[0]).text().length > 0) {
                                readFilePromises.push(makeDbFilePromise(fileOfType, "application/CCP4-mtz-map"))
                            }
                            else if (projectNodes.length === 1 && $(projectNodes[0]).text().length > 0) {
                                readFilePromises.push(makeParamsFilePromise(fileOfType, "application/CCP4-mtz-map"))
                            }
                        })
                        return Promise.all(readFilePromises)
                    })
            })
    }

    useEffect(() => {
        if (cootInitialized) {
            if (props.cootJob) {
                handleCootJob()
            }
        }
    }, [cootInitialized])

    const ccp4i2ExtraMenu = {
        name: "CCP4i2",
        ref: ccp4i2MenuRef,
        icon: <Avatar><img src={`${props.urlRoot}/svgicons/ccp4`} /></Avatar>,
        JSXElement: <MoorhenCCP4i2Menu
            key={"CCP4i2"}
            dropdownId={7}
            {...collectedProps}
            handleJob={(jobId) => { /*handleJob(jobId)*/ }}
            cootJobId={props.cootJob}
            makeFilePromise={makeFilePromise}
            setShowRefmacModal={setShowRefmacModal}
        />
    }

    return <>
        <MoorhenContainer
            {...collectedProps}
            extraNavBarMenus={[ccp4i2ExtraMenu]}
            extraFileMenuItems={[
                <hr></hr>,
                <MenuItem key="Find files" onClick={() => {
                    setShowModal(true)
                }}>Add files from CCP4i2
                </MenuItem>,
                <MenuItem key="Save to CCP4i2" onClick={
                    async () => {
                        saveMolNo(props.urlRoot, 0, props.cootJob, iSave, setISave, molecules, setToastContent)
                    }}>Save first molecule to CCP4i2
                </MenuItem>,
                <MenuItem key="Select molecule to save" onClick={
                    () => {
                        setShowMoleculeSelectorModal(true)
                    }}>Select molecule to save to CCP4i2
                </MenuItem>
            ]}
            controls={controls}
            urlPrefix={""}
        />

        <MoorhenDraggableModalBase key="Browse database modal"
            transparentOnMouseOut={false}
            left={`400px`}
            top={`200px`}
            show={showModal}
            setShow={setShowModal}
            windowHeight={1200}
            windowWidth={1800}
            headerTitle={'CCP4i2 Database'}
            additionalHeaderButtons={[
                <Button variant="white">
                </Button>
            ]}
            body={
                <CCP4i2DataBrowser {...collectedProps} updated={showModal} makeFilePromise={makeFilePromise} />
            }
            footer={null} />

        <MoorhenDraggableModalBase key="Refmac modal"
            transparentOnMouseOut={false}
            left={`400px`}
            top={`200px`}
            show={showRefmacModal}
            setShow={setShowRefmacModal}
            windowHeight={1200}
            windowWidth={1200}
            headerTitle={'Run refmac'}
            additionalHeaderButtons={[
                <Button variant="white">
                </Button>
            ]}
            body={
                <CCP4i2RefmacPanel {...collectedProps} setShowRefmacModal={setShowRefmacModal}
                    molecules={molecules}
                    makeFilePromise={makeFilePromise}
                    cootJob={props.cootJob} urlPrefix="" />
            }
            footer={null} />
        <Modal key="Make aceDRG Link modal" show={showAceDRGLinkModal} onHide={() => {
            setShowAceDRGLinkModal(false)
        }}>
            <Modal.Body>
                {currentJob && <CCP4i2MakeAceDRGLinkPanel {...collectedProps}
                    molecules={molecules}
                    linkData={linkData}
                    setShowAceDRGLinkModal={showAceDRGLinkModal}
                    projectId={currentJob.projectid__projectid}
                    makeFilePromise={makeFilePromise} urlPrefix="" />}
            </Modal.Body>
        </Modal>

        <Modal key="Save selector modal" show={showMoleculeSelectorModal} onHide={() => {
            setShowMoleculeSelectorModal(false)
        }}><Modal.Body>
                {cootInitialized &&
                    <MoorhenMoleculeSelect width="" ref={selectorRef} molecules={molecules} allowAny={false} />
                }
                <Button onClick={() => {
                    saveMolNo(props.urlRoot, parseInt(selectorRef.current.value))
                    setShowMoleculeSelectorModal(false)
                }}>OK</Button>
                <Button onClick={() => {
                    setShowMoleculeSelectorModal(false)
                }}>Cancel</Button>
            </Modal.Body>
        </Modal>
    </>
}

CCP4i2MoorhenContainer.defaultProps = { job: { jobid: null }, cootJob: null }

const MoorhenCCP4i2Menu = (props) => {
    const [popoverIsShown, setPopoverIsShown] = useState(false)
    const menuItemProps = { setPopoverIsShown, ...props }

    return <>
        <MenuItem id='configure-refmac-menu-item' onClick={() => {
            props.setShowRefmacModal(true)
            props.setCurrentDropdownId('-1')
        }}>Run refmac
        </MenuItem>

    </>
}


