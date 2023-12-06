import ListItem from '@mui/material/ListItem';
import ListItemAvatar from '@mui/material/ListItemAvatar';
import ListItemText from '@mui/material/ListItemText';
import Avatar from '@mui/material/Avatar';
import { CCP4i2BaseBrowser } from './CCP4i2BaseBrowser';

export const fileTypeMapping = {
    "application/CCP4-mtz-observed": "ObsDataFile",
    "application/CCP4-mtz-freerflag": "FreeRDataFile",
    "application/CCP4-mtz-map": "MapCoeffsDataFile",
    "application/CCP4-mtz-phases": "PhsDataFile",
    "application/refmac-dictionary": "DictDataFile",
    "application/coot-script": "CootHistoryDataFile",
    "application/CCP4-unmerged-experimental": "UnmergedDataFile",
    "chemical/x-pdb": "PdbDataFile",
    "application/CCP4-asu-content": "AsuDataFile"
}

export const fileFields = ["filetypeid__filetypename", "filetypeid__filetypeid", "annotation", "fileid",
    "jobid__jobnumber", "jobid__jobid", "filesubtype"]

export const CCP4i2FileBrowser = (props) => {
    return <CCP4i2BaseBrowser {...props}
        dbType="Files"
        sort={(a, b) => {
            const jnA = parseInt(a.filetypeid__filetypeid)
            const jnB = parseInt(b.filetypeid__filetypeid)
            return jnA < jnB ? 1 : jnB < jnA ? -1 : 0
        }}
        fields={fileFields}
        avatar={file => <CCP4i2FileAvatar file={file} />}
        primary={(file) => { return `${file.jobid__jobnumber}:${file.annotation}` }}
        secondary=""
        listItem={file => <CCP4i2FileListItem file={file} {...props} />}
    />
}

export const CCP4i2FileAvatar = (props) => {
    const { file } = props;
    return <Avatar>
        <img
            style={{ width: "2.0rem", height: "2.0rem", marginRight: "1rem" }}
            title={file.filetypeid__filetypename}
            src={`${props.urlRoot}/qticons/${fileTypeMapping[file.filetypeid__filetypename]}`} alt={`${fileTypeMapping[file.filetypeid__filetypename]} `} />
        {file.jobid__jobnumber}: {file.annotation}
    </Avatar>
}

export const CCP4i2FileListItem = (props) => {
    const { file } = props;
    return <ListItem
        secondaryAction={props.secondaryAction}
        onClick={() => { if (props.onClick) { props.onClick(file) } }}
        dense
    >
        <ListItemAvatar>
            <CCP4i2FileAvatar {...props} />
        </ListItemAvatar>
        <ListItemText primary={`${file.jobid__jobnumber}:${file.annotation}`} />
    </ListItem>
}