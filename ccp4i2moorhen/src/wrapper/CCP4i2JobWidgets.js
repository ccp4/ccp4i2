import ListItem from '@mui/material/ListItem';
import ListItemAvatar from '@mui/material/ListItemAvatar';
import ListItemText from '@mui/material/ListItemText';
import Avatar from '@mui/material/Avatar';
import { CCP4i2BaseBrowser } from './CCP4i2BaseBrowser';
import { CircularProgress } from '@mui/material';

export const jobFields = ['projectid__projectname', 'projectid__projectid',
    'projectid__projectdirectory', 'taskname', 'jobnumber', 'jobid', 'jobtitle',
    'status__statustext', 'creationtime', 'finishtime', 'status__statusid',
    'processid', 'parentjobid']

export const CCP4i2JobBrowser = (props) => {
    return <CCP4i2BaseBrowser {...props}
        dbType="Jobs"
        sort={(a, b) => {
            const jnA = parseInt(a.jobnumber.split('.')[0])
            const jnB = parseInt(b.jobnumber.split('.')[0])
            return jnA < jnB ? 1 : jnB < jnA ? -1 : 0
        }}
        fields={jobFields}
        avatar={job => <CCP4i2JobAvatar job={job} />}
        primary={(job) => { return `${job.jobnumber}:${job.annotation ? job.annotation : job.tasktitle ? job.tasktitle : job.taskname}` }}
        secondary=""
        listItem={job => <CCP4i2JobListItem job={job} {...props} />}
    />
}

export const CCP4i2JobAvatar = (props) => {
    const { job } = props;
    return job.status__statustext !== "Running" ? <Avatar style={{
        backgroundColor:
            job.status__statustext === "Finished" ? "#afa" :
                job.status__statustext === "Failed" ? "#faa" :
                    job.status__statustext === "Running" ? "#aaf" :
                        job.status__statustext === "Queued" ? "#ffa" :
                            job.status__statustext === "Pending" ? "#fff" :
                                "black"
    }}>
        <img
            style={{ width: "2.0rem", height: "2.0rem", marginRight: "1rem" }}
            title={job.taskname}
            src={`/database/svgicons/${job.taskname}`} alt={`${job.taskname} `} />
    </Avatar> : <CircularProgress />
}

export const CCP4i2JobListItem = (props) => {
    const { job } = props;
    return <ListItem
        secondaryAction={props.secondaryAction}
        onClick={() => { if (props.onClick) { props.onClick(job) } }}
        dense
    >
        <ListItemAvatar>
            <CCP4i2JobAvatar {...props} />
        </ListItemAvatar>
        <ListItemText primary={`${job.jobnumber}:${job.taskname}`} secondary={job.annotation} />
    </ListItem>
}