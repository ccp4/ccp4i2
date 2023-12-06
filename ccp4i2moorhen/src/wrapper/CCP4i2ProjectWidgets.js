import { ListItem, ListItemAvatar, ListItemText } from '@mui/material';
import { CCP4i2BaseBrowser} from './CCP4i2BaseBrowser';

export const CCP4i2ProjectBrowser = (props) => {
    return <CCP4i2BaseBrowser {...props}
        dbType="Projects"
        sort={(a, b) => {
            const jnA = a.lastaccess
            const jnB = b.lastaccess
            return jnA < jnB ? 1 : jnB < jnA ? -1 : 0
        }}
        avatar={project => "hi"}
        primary={(project) => { return `${project.projectname}` }}
        secondary=""
        listItem={project => <CCP4i2ProjectListItem project={project} {...props} />}
    />
}

export const CCP4i2ProjectListItem = (props) => {
    const { project } = props;
    return <ListItem
        secondaryAction={props.secondaryAction}
        onClick={() => { if (props.onClick) { props.onClick(project) } }}
        dense
    >
        <ListItemText primary={`${project.projectname}`} secondary={project.projectdirectory} />
    </ListItem>
}