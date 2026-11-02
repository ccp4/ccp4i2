import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CredentialAlert } from "../task-elements/credential-alert";
import { useJob } from "../../../utils";

/**
 * Submit to OneDep.
 *
 * Everything descriptive about the entry (title, authors, sources, citation,
 * data collection) travels inside the model mmCIF the deposition-preparation
 * task wrote, so this form asks only who is depositing. The API key is a
 * credential, not a parameter: it must never reach the project database or an
 * exported project, so it is set through the banner, not a field.
 */
const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, mutateValidation } = useJob(props.job.id);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CredentialAlert name="onedep" onChanged={() => mutateValidation?.()} />
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Files to deposit" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="MODEL_CIF" {...props} />
        <CCP4i2TaskElement itemName="SF_CIF" {...props} />
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Depositor" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="EMAIL" {...props} />
        <CCP4i2TaskElement itemName="ORCIDS" {...props} />
        <CCP4i2TaskElement itemName="COUNTRY" {...props} />
        <CCP4i2TaskElement itemName="EXPERIMENT_TYPE" {...props} />
        <Typography variant="body2" color="text.secondary">
          The deposition belongs to the ORCID the API key was generated for.
          Every ORCID listed can open it in the OneDep interface, where the
          remaining pages (grant, release status, terms) are completed.
        </Typography>
      </CCP4i2ContainerElement>

      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Existing deposition" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement itemName="REPLACE_DEPOSITION_ID" {...props} />
        <Typography variant="body2" color="text.secondary">
          Leave empty to create a new deposition. Give an existing D_ id to
          replace its model and reflections with these files; OneDep then
          reprocesses them and the pages it fills from the model are updated.
        </Typography>
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
