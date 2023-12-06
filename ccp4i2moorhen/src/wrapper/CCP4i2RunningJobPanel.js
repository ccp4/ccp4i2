import { useState, useEffect } from 'react'
import { Paper, List } from '@mui/material'
import { CCP4i2JobListItem, jobFields } from './CCP4i2JobWidgets'
import { CCP4i2FileBrowser } from './CCP4i2FileWidgets'
import { CCP4i2ReportChart, reportXmlToCharts } from './CCP4i2ReportChart'
import $ from 'jquery'

export const CCP4i2RunningJobPanel = (props) => {
    const [runningJob, setRunningJob] = useState(null)
    const [finishedStatus, setFinishedStatus] = useState(null)
    const [mainChart, setMainChart] = useState(null)
    const { runningJobId } = props

    useEffect(() => {
        if (runningJobId) {
            const pulseFunc = async () => {
                const jobDetailsResult = await fetch(`${props.urlRoot}/ModelValues?${$.param({
                    __type__: "Jobs",
                    __values__: jobFields,
                    jobid: runningJobId
                })}`).then(response => response.json())
                setRunningJob(jobDetailsResult.results[0])
                if (["Pending", "Queued", "Running"].includes(jobDetailsResult.results[0].status__statustext)) {
                    setTimeout(() => { pulseFunc() }, 5000)
                }
                else {
                    setFinishedStatus(jobDetailsResult.results[0].status__statustext)
                }
            }
            pulseFunc()
        }
    }, [runningJobId])

    useEffect(() => {
        if (finishedStatus === 'Finished' && props.mainChartKey) {
            const fetchFunc = async () => {
                const jobReportText = await fetch(`${props.urlRoot}/getReportXML?${$.param({
                    jobId: runningJobId
                })}`).then(response => response.text())
                const xmlDoc = $.parseXML(jobReportText)
                const jQueryXMLDoc = $(xmlDoc);
                const mainGraph = jQueryXMLDoc.find(props.mainChartKey).toArray()[0]
                const graphAsStruct = reportXmlToCharts(mainGraph)
                setMainChart(graphAsStruct[0])
            }
            fetchFunc()
        }
    }, [finishedStatus, props.mainChartKey])

    return <>
        <Paper>
            {runningJob && <h1><CCP4i2JobListItem {...props} job={runningJob} /></h1>}
            {mainChart && <CCP4i2ReportChart height={200} width={300} chart={mainChart}></CCP4i2ReportChart>}
            {(runningJob && finishedStatus !== null) && <Paper>
                <h5>Job files</h5>
                <CCP4i2FileBrowser {...props}
                    predicate={{
                        jobid__jobid: runningJob.jobid,
                        filetypeid__filetypename__in: ["application/CCP4-mtz-map", "application/refmac-dictionary", "chemical/x-pdb"]
                    }}
                    onClick={(file) => {
                        const fileUrl = `${props.urlRoot}/getFileWithPredicate?fileid=${file.fileid}`
                        const mimeType = file.filetypeid__filetypename
                        const annotation = `${file.jobid__jobnumber} : ${file.annotation}`
                        const subType = file.filesubtype
                        props.makeFilePromise(fileUrl, mimeType, annotation, subType)
                    }} />
            </Paper>}
        </Paper >
    </>
}
