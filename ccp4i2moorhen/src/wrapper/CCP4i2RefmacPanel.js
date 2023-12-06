import { useState } from 'react';
import { useEffect } from 'react';
import { useRef } from 'react';
import { MoorhenMoleculeSelect } from 'moorhen';
import { Button, Popover, Paper, Toolbar, List } from '@mui/material';
import $ from 'jquery';
import { ControlPoint } from '@mui/icons-material';
import { CCP4i2JobBrowser, CCP4i2JobListItem, jobFields } from './CCP4i2JobWidgets';
import { CCP4i2RunningJobPanel } from './CCP4i2RunningJobPanel';
import { CCP4i2Api } from './CCP4i2Api'

export const CCP4i2RefmacPanel = (props) => {
    const [predicate, setPredicate] = useState({})
    const moleculeSelectRef = useRef(null)
    const [selectedMolecule, setSelectedMolecule] = useState(null)
    const [refmacJobs, setRefmacJobs] = useState([])
    const [clonedJob, setClonedJob] = useState(null)
    const [open, setOpen] = useState(false);
    const [runningJobId, setRunningJobId] = useState(null)
    const [anchorEl, setAnchorEl] = useState(null);

    const handleClick = (event) => {
        setAnchorEl(event.currentTarget);
        setOpen(true)
    };

    const handleClose = () => {
        setAnchorEl(null);
    };

    const popoverOpen = Boolean(anchorEl);
    const id = popoverOpen ? 'simple-popover' : undefined;

    useEffect(() => {
        const asyncFunc = async () => {
            const jobId = props.cootJob ? props.cootJob : props.jobId ? props.jobId : null;
            if (jobId) {
                const jobFetchResult = await fetch(`${props.urlRoot}/ModelValues?${$.param({
                    __type__: "Jobs", jobid: jobId, __values__: jobFields
                })}`).then(response => response.json())
                const projectFetchResult = await fetch(`${props.urlRoot}/ModelValues?${$.param({
                    __type__: "Projects", projectid: jobFetchResult.results[0].projectid__projectid
                })}`).then(response => response.json())
                setPredicate({
                    projectid__projectid: projectFetchResult.results[0].projectid,
                    projectid__projectname: projectFetchResult.results[0].projectname,
                    status__statustext: "Finished",
                    taskname: "prosmart_refmac",
                    parentjobid__isnull: true,
                });
            }
        }
        asyncFunc()
    }, [])

    //useEffect to fetch the paramsxml of the cloned Job ...could provide some UI to edit some of these
    useEffect(() => {
        if (clonedJob) {
            const asyncFunc = async () => {
                const paramsText = await fetch(`${props.urlRoot}/getJobFile?jobId=${clonedJob.jobid}&fileName=input_params.xml`)
                    .then(response => response.text())
                const paramsXml = $.parseXML(paramsText)
            }
            asyncFunc()
        }
    }, [clonedJob])

    useEffect(() => {
        if (Object.keys(predicate).length > 0) {
            const asyncFunc = async () => {
                const newJobsResult = await fetch(`${props.urlRoot}/ModelValues?${$.param({
                    __type__: "Jobs", ...predicate, __values__: jobFields
                })}`).then(response => response.json())
                setRefmacJobs(newJobsResult.results)
                if (newJobsResult.results.length > 0) {
                    setClonedJob(newJobsResult.results.at(0))
                }
            }
            asyncFunc()
        }
    }, [predicate])

    const createAndRun = async () => {
        //Clone the selected job and recover the new jobId
        const ccp4i2Api = new CCP4i2Api(props.urlRoot)

        const cloneFormData = new FormData()
        cloneFormData.set('jobId', clonedJob.jobid)
        const cloneResult = await fetch(`${props.urlRoot}/cloneJob`, {
            method: 'POST',
            body: cloneFormData
        }).then(response => response.json())

        //Upload the selected coordinate set to the new Job
        const molToRefine = parseInt(moleculeSelectRef.current.value)
        const molZeros = props.molecules.filter(molecule => molecule.molNo === molToRefine)
        let pdbData = await molZeros[0].getAtoms()
        const atomsBlob = new Blob([pdbData])
        const uploadResult = await ccp4i2Api.uploadFileForJobObject(cloneResult.jobId, 'inputData.XYZIN', atomsBlob, `${molZeros[0].name}.pdb`)

        //Run the new job
        const runFormData = new FormData()
        runFormData.set('jobId', cloneResult.jobId)
        const runResult = await fetch(`${props.urlRoot}/runJob`, {
            method: 'POST',
            body: runFormData
        }).then(response => {
            setRunningJobId(cloneResult.jobId)
            return response.json()
        })
    }


    return <>
        <Paper>
            {!runningJobId && <Paper>
                <MoorhenMoleculeSelect {...props} ref={moleculeSelectRef} onChange={(evt) => {
                    setSelectedMolecule(evt.currentTarget.value)
                }} />
                <h6>Refinement based on job:</h6>
                {clonedJob ? <CCP4i2JobListItem dense job={clonedJob}
                    secondaryAction={<Button onClick={handleClick} ><ControlPoint /></Button>}
                /> : "Job not yet chosen"}
                <Toolbar sx={{ justifyContent: "space-between" }}>
                    <div />
                    <Button variant="outlined" onClick={createAndRun}>Create and run</Button>
                </Toolbar>
            </Paper>}
            {runningJobId && <CCP4i2RunningJobPanel
                {...props}
                runningJobId={runningJobId}
                mainChartKey='CCP4i2ReportFlotGraph[key="SummaryGraph"]' />}
        </Paper >
        <Popover
            id={id}
            open={open}
            anchorEl={anchorEl}
            onClose={handleClose}
            anchorOrigin={{
                vertical: 'bottom',
                horizontal: 'left',
            }}
        >
            <CCP4i2JobBrowser {...props} predicate={predicate} onClick={job => {
                setClonedJob(job)
                setOpen(false)
                setAnchorEl(null)
            }} />
        </Popover>
    </>
}