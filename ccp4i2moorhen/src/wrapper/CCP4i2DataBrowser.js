import { useRef, useState, useReducer } from 'react';
import { Button, Tab, Tabs, Typography } from '@mui/material';
import { ArrowBack } from '@mui/icons-material';
import { Slide, Box, Paper } from '@mui/material';
import { CCP4i2JobBrowser } from './CCP4i2JobWidgets';
import { CCP4i2ProjectBrowser } from './CCP4i2ProjectWidgets';
import { CCP4i2FileBrowser } from './CCP4i2FileWidgets';

function CustomTabPanel(props) {
    const { children, value, index, ...other } = props;

    return (
        <div
            role="tabpanel"
            hidden={value !== index}
            id={`simple-tabpanel-${index}`}
            aria-labelledby={`simple-tab-${index}`}
            {...other}
        >
            {value === index && (
                <Box sx={{ p: 3 }}>
                    <Typography>{children}</Typography>
                </Box>
            )}
        </div>
    );
}

function a11yProps(index) {
    return {
        id: `simple-tab-${index}`,
        'aria-controls': `simple-tabpanel-${index}`,
    };
}

const initialLevel = { levelWas: 'Nothing', levelIs: 'Projects', direction: 'left', fromRef: null }
export const CCP4i2DataBrowser = (props) => {
    const levelReducer = (oldLevel, newLevel) => {
        let newLevelDict = { levelWas: oldLevel.levelIs, levelIs: newLevel }
        if (oldLevel.is === 'Projects') {
            newLevelDict.fromRef = projectsContainerRef
            newLevelDict.direction = 'left'
        }
        else if (oldLevel.is === 'Files') {
            newLevelDict.fromRef = filesContainerRef
            newLevelDict.direction = 'right'
        }
        if (oldLevel.is === 'Jobs') {
            newLevelDict.fromRef = jobsContainerRef
            if (newLevel === 'Projects') {
                newLevelDict.direction = 'right'
            }
            else if (newLevel === 'Files') {
                newLevelDict.direction = 'left'
            }
        }
        return newLevelDict
    }
    const [level, setLevel] = useReducer(levelReducer, initialLevel)
    const [project, setProject] = useState({})
    const [job, setJob] = useState({})

    const projectsContainerRef = useRef(null)
    const jobsContainerRef = useRef(null)
    const filesContainerRef = useRef(null)
    const containerRef = useRef(null)
    const [value, setValue] = useState(0);

    const handleChange = (event, newValue) => {
        setValue(newValue);
    };

    return <Box sx={{
        padding: 2,
        borderRadius: 1,
        bgcolor: (theme) =>
            theme.palette.mode === 'light' ? 'grey.100' : 'grey.900',
        overflow: 'hidden'
    }}>
        <Paper elevation={3} sx={{ width: "100%", height: "100%", overflow: "hidden" }}>
            <Tabs value={value} onChange={handleChange} aria-label="basic tabs example">
                <Tab label="Projects" {...a11yProps(0)} />
                <Tab label="Jobs" {...a11yProps(1)} />
                <Tab label="Files" {...a11yProps(2)} />
            </Tabs>
            <CustomTabPanel value={value} index={0}>
                <Paper elevation={4}  >
                    <h5>Projects</h5>
                    <CCP4i2ProjectBrowser {...props}
                        predicate={{
                            parentprojectid__isnull: true,
                        }}
                        onClick={project => {
                            setProject(project)
                            setLevel('Jobs')
                            setValue(1)
                        }} />
                </Paper>
            </CustomTabPanel>
            <CustomTabPanel value={value} index={1}>
                <Paper elevation={5}>
                    <h5><Button onClick={() => { setValue(0) }}><ArrowBack /></Button> Jobs of {project.projectname}</h5>
                    <CCP4i2JobBrowser {...props}
                        predicate={{
                            parentjobid__isnull: true,
                            projectid__projectid: project.projectid
                        }}
                        project={project}
                        onClick={job => {
                            setJob(job)
                            setValue(2)
                        }} />
                </Paper>
            </CustomTabPanel>
            <CustomTabPanel value={value} index={2}>
                <Paper elevation={6}>
                    <h5><Button onClick={() => { setValue(1) }}><ArrowBack /></Button> Files of {job.jobnumber} : {job.annotation}</h5>
                    <CCP4i2FileBrowser {...props}
                        predicate={{
                            jobid__jobid: job.jobid,
                            filetypeid__filetypename__in: ["application/CCP4-mtz-map", "application/refmac-dictionary", "chemical/x-pdb"]
                        }}
                        onClick={(file) => {
                            const fileUrl = `${props.urlRoot}/getFileWithPredicate?fileid=${file.fileid}`
                            const mimeType = file.filetypeid__filetypename
                            const annotation = `${file.jobid__jobnumber} : ${file.annotation}`
                            const subType = file.filesubtype
                            props.makeFilePromise(fileUrl, mimeType, annotation, subType)
                        }} />
                </Paper>
            </CustomTabPanel>

        </Paper >

    </Box >
}