import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { InlineField } from "../task-elements/inline-field";
import { useJob } from "../../../utils";
import { useBoolToggle } from "../task-elements/shared-hooks";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { useTaskItem, container } = useJob(props.job.id);
  const { value: compBy } = useTaskItem("COMP_BY");

  const resolOn = useBoolToggle(useTaskItem, "RESOL_ON");
  const llgComplOn = useBoolToggle(useTaskItem, "LLG_COMPL_ON");
  const llgSigcoOn = useBoolToggle(useTaskItem, "LLG_COMPL_SIGCO_ON");
  const llgAtmsepOn = useBoolToggle(useTaskItem, "LLG_COMPL_ATMSEP_ON");
  const llgMaxcycOn = useBoolToggle(useTaskItem, "LLG_COMPL_MAXCYC_ON");

  const expPackOn = useBoolToggle(useTaskItem, "EXP_PACKCT_ON");
  const expTranSrchOn = useBoolToggle(useTaskItem, "EXP_TRAN_SRCHPK_ON");
  const expPurgeOn = useBoolToggle(useTaskItem, "EXP_PURGE_TRANPK_ON");
  const expLlgCrtOn = useBoolToggle(useTaskItem, "EXP_EXPLLG_CRT_ON");
  const expResolOn = useBoolToggle(useTaskItem, "EXP_RESRAN_HRREFINE_ON");

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs>
        {/* ── Tab 1: Input Data and Run Parameters ── */}
        <CCP4i2Tab key="inputData" label="Input Data and Run Parameters">
          {/* Select input data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Select input data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement itemName="F_SIGF" {...props} />
            <CCP4i2TaskElement itemName="FREERFLAG" {...props} />

            <InlineField
              width="auto"
              after={
                <InlineField label="Angstroms to" width="6rem">
                  <CCP4i2TaskElement
                    itemName="RESOL_HI"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                itemName="RESOL_ON"
                {...props}
                qualifiers={{ guiLabel: "Use resolution range :" }}
                onChange={resolOn.onChange}
                sx={{ width: "auto" }}
              />
              {resolOn.value && (
                <InlineField width="6rem">
                  <CCP4i2TaskElement
                    itemName="RESOL_LO"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              )}
            </InlineField>
          </CCP4i2ContainerElement>

          {/* Composition */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Composition" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="COMP_BY"
              {...props}
              qualifiers={{ guiLabel: "Composition of asymmetric unit:" }}
            />
            {compBy === "ASU" && (
              <>
                <CCP4i2TaskElement
                  itemName="ASUFILE"
                  {...props}
                  qualifiers={{ guiLabel: "AU contents" }}
                />
                <Typography
                  variant="caption"
                  color="text.secondary"
                  sx={{ pl: 1 }}
                >
                  If a suitable ASU is not available above, you can press the
                  cross &amp; then button to quickly create one.
                </Typography>
              </>
            )}
            {compBy === "MW" && (
              <>
                <CCP4i2TaskElement
                  itemName="ASU_PROTEIN_MW"
                  {...props}
                  qualifiers={{ guiLabel: "Protein molecular weight" }}
                />
                <CCP4i2TaskElement
                  itemName="ASU_NUCLEICACID_MW"
                  {...props}
                  qualifiers={{ guiLabel: "Nucleic acid molecular weight" }}
                />
              </>
            )}
          </CCP4i2ContainerElement>

          {/* Define Search */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Define Search" }}
            containerHint="FolderLevel"
          >
            <InlineField
              label="Search for"
              hint="atom(s) of type"
              after={
                <InlineField width="4rem">
                  <CCP4i2TaskElement
                    itemName="SINGLE_ATOM_TYPE"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                itemName="SINGLE_ATOM_NUM"
                {...props}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>

            <CCP4i2TaskElement
              itemName="LLG_COMPL_ON"
              {...props}
              qualifiers={{ guiLabel: "Turn on Log-likelihood Completion" }}
              onChange={llgComplOn.onChange}
            />

            {llgComplOn.value && (
              <>
                <InlineField
                  label="Log-likelihood Gain Map (LLG) Completion with"
                  hint="atoms"
                  width="4rem"
                >
                  <CCP4i2TaskElement
                    itemName="LLG_COMPL_ATOMTYP"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>

                <InlineField width="auto" after={
                  <InlineField width="6rem">
                    <CCP4i2TaskElement
                      itemName="LLG_COMPL_SIGCO"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                }>
                  <CCP4i2TaskElement
                    itemName="LLG_COMPL_SIGCO_ON"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "LLG-Map : sigma cutoff for adding new atom sites",
                    }}
                    onChange={llgSigcoOn.onChange}
                    sx={{ width: "auto" }}
                  />
                </InlineField>

                <InlineField width="auto" after={
                  <InlineField width="6rem">
                    <CCP4i2TaskElement
                      itemName="LLG_COMPL_ATMSEP"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                }>
                  <CCP4i2TaskElement
                    itemName="LLG_COMPL_ATMSEP_ON"
                    {...props}
                    qualifiers={{
                      guiLabel: "LLG-Map : atomic separation distance cut-off",
                    }}
                    onChange={llgAtmsepOn.onChange}
                    sx={{ width: "auto" }}
                  />
                </InlineField>

                <InlineField width="auto" after={
                  <InlineField width="6rem">
                    <CCP4i2TaskElement
                      itemName="LLG_COMPL_MAXCYC"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                }>
                  <CCP4i2TaskElement
                    itemName="LLG_COMPL_MAXCYC_ON"
                    {...props}
                    qualifiers={{
                      guiLabel:
                        "LLG-Map : max number of completion cycles to use",
                    }}
                    onChange={llgMaxcycOn.onChange}
                    sx={{ width: "auto" }}
                  />
                </InlineField>
              </>
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ── Tab 2: Expert Settings ── */}
        <CCP4i2Tab key="inputDataExp" label="Expert Settings">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Expert Inputs" }}
            containerHint="FolderLevel"
          >
            {/* Packing criteria */}
            <InlineField
              width="auto"
              after={
                <>
                  <CCP4i2TaskElement
                    itemName="EXP_PACKCT_TYPE"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                  <InlineField label="at" width="6rem">
                    <CCP4i2TaskElement
                      itemName="EXP_PACKCT_AMT"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                </>
              }
            >
              <CCP4i2TaskElement
                itemName="EXP_PACKCT_ON"
                {...props}
                qualifiers={{ guiLabel: "Set packing criteria :" }}
                onChange={expPackOn.onChange}
                sx={{ width: "auto" }}
              />
            </InlineField>

            {/* Translation search peak criteria */}
            <InlineField
              width="auto"
              after={
                <>
                  <CCP4i2TaskElement
                    itemName="EXP_TRAN_SRCHPK_TYPE"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                  <InlineField label="at" width="6rem">
                    <CCP4i2TaskElement
                      itemName="EXP_TRAN_SRCHPK_AMT"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                </>
              }
            >
              <CCP4i2TaskElement
                itemName="EXP_TRAN_SRCHPK_ON"
                {...props}
                qualifiers={{ guiLabel: "Translation search peak criteria :" }}
                onChange={expTranSrchOn.onChange}
                sx={{ width: "auto" }}
              />
            </InlineField>

            {/* Purge translation peaks */}
            <InlineField
              width="auto"
              after={
                <>
                  <InlineField
                    label="% cutoff set to"
                    width="6rem"
                  >
                    <CCP4i2TaskElement
                      itemName="EXP_PURGE_TRANPK_PER"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                  <InlineField label="but with maximum number" width="6rem">
                    <CCP4i2TaskElement
                      itemName="EXP_PURGE_TRANPK_NUM"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                </>
              }
            >
              <CCP4i2TaskElement
                itemName="EXP_PURGE_TRANPK_ON"
                {...props}
                qualifiers={{
                  guiLabel: "Purge translation peaks with",
                }}
                onChange={expPurgeOn.onChange}
                sx={{ width: "auto" }}
              />
            </InlineField>

            {/* Expected LLG criteria */}
            <InlineField
              width="auto"
              after={
                <InlineField
                  label="where the targeted gain value is"
                  width="6rem"
                >
                  <CCP4i2TaskElement
                    itemName="EXP_EXPLLG_CRT"
                    {...props}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
              }
            >
              <CCP4i2TaskElement
                itemName="EXP_EXPLLG_CRT_ON"
                {...props}
                qualifiers={{
                  guiLabel: "Use expected LLG criteria,",
                }}
                onChange={expLlgCrtOn.onChange}
                sx={{ width: "auto" }}
              />
            </InlineField>

            {/* Resolution range */}
            <InlineField
              width="auto"
              after={
                <>
                  <InlineField width="6rem">
                    <CCP4i2TaskElement
                      itemName="EXP_RESRAN_HRREFINE_LO"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                  <InlineField label="Angstroms to" width="6rem">
                    <CCP4i2TaskElement
                      itemName="EXP_RESRAN_HRREFINE_HI"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </InlineField>
                </>
              }
            >
              <CCP4i2TaskElement
                itemName="EXP_RESRAN_HRREFINE_ON"
                {...props}
                qualifiers={{ guiLabel: "Use resolution range :" }}
                onChange={expResolOn.onChange}
                sx={{ width: "auto" }}
              />
            </InlineField>

            {/* Translational NCS */}
            <CCP4i2TaskElement
              itemName="EXP_TRAN_NCS_ON"
              {...props}
              qualifiers={{
                guiLabel: "Use translational NCS if present",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
