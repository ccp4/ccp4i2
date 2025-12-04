import { useCallback, useContext, useEffect, useMemo } from "react";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useRunCheck } from "../../../providers/run-check-provider";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, validation } = useJob(job.id);
  const { processedErrors, setProcessedErrors } = useRunCheck();

  // Consolidated task values
  const taskValues = useMemo(
    () => ({
      chooseMode: useTaskItem("CHOOSE_MODE").value,
      mode: useTaskItem("MODE").value,
      hklinRef: useTaskItem("HKLIN_REF").value,
      aimlessRef: useTaskItem("REFERENCE_FOR_AIMLESS").value,
      referenceDataset: useTaskItem("REFERENCE_DATASET").value,
    }),
    [useTaskItem]
  );

  // Process validation errors
  const processedValidationErrors = useMemo(() => {
    if (!validation) return null;

    const filtered = Object.keys(validation)
      .filter((key) => !key.startsWith("aimless_pipe.controlParameters.CELL."))
      .reduce((acc, key) => ({ ...acc, [key]: validation[key] }), {});

    // Add custom validation for HKLIN_REF
    if (
      taskValues.mode === "MATCH" &&
      taskValues.aimlessRef &&
      taskValues.referenceDataset === "HKL" &&
      !taskValues.hklinRef?.dbFileId
    ) {
      filtered["aimless_pipe.inputData.HKLIN_REF"] = {
        messages: ["HKLIN_REF must be set when being used for match"],
        maxSeverity: 2,
      };
    }

    return filtered;
  }, [validation, taskValues]);

  // Update processed errors efficiently
  useEffect(() => {
    if (
      JSON.stringify(processedValidationErrors) !==
      JSON.stringify(processedErrors)
    ) {
      setProcessedErrors(processedValidationErrors);
    }
  }, [processedValidationErrors, processedErrors, setProcessedErrors]);

  // Visibility helpers
  const visibility = useMemo(
    () => ({
      isChooseMode: () => taskValues.mode === "CHOOSE",
      isMatchMode: () => taskValues.mode === "MATCH",
      isChooseSolution: () =>
        taskValues.mode === "CHOOSE" && taskValues.chooseMode === "SOLUTION_NO",
      isChooseSpacegroup: () =>
        taskValues.mode === "CHOOSE" &&
        ["SPACEGROUP", "REINDEX_SPACE"].includes(taskValues.chooseMode),
      isReindexSpace: () =>
        taskValues.mode === "CHOOSE" &&
        taskValues.chooseMode === "REINDEX_SPACE",
      isChooseLauegroup: () =>
        taskValues.mode === "CHOOSE" && taskValues.chooseMode === "LAUEGROUP",
      hasAimlessRef: () => taskValues.mode === "MATCH" && taskValues.aimlessRef,
      isHklReference: () =>
        taskValues.mode === "MATCH" &&
        taskValues.aimlessRef &&
        taskValues.referenceDataset === "HKL",
      isXyzReference: () =>
        taskValues.mode === "MATCH" &&
        taskValues.aimlessRef &&
        taskValues.referenceDataset === "XYZ",
    }),
    [taskValues]
  );

  // Element configurations
  const elementConfigs = useMemo(
    () => ({
      fileInputs: [
        { key: "UNMERGEDFILES", label: "Unmerged files" },
        { key: "FREERFLAG", label: "Free R set to use/extend" },
        {
          key: "HKLIN_IS_SCALED",
          label: "Analyse data without determining scales",
        },
      ],
      parameters: [
        { key: "AUTOCUTOFF", label: "Apply auto. data cutoff" },
        { key: "RESOLUTION_RANGE", label: "Resolution" },
        { key: "OVERRIDE_CELL_DIFFERENCE", label: "Override cell difference" },
      ],
      choiceOptions: [
        {
          key: "CHOOSE_MODE",
          label: "Symmetry choice mode",
          visible: visibility.isChooseMode,
        },
        {
          key: "CHOOSE_SOLUTION_NO",
          label: "Solution no. to choose",
          visible: visibility.isChooseSolution,
        },
        {
          key: "CHOOSE_SPACEGROUP",
          label: "Spacegroup to choose",
          visible: visibility.isChooseSpacegroup,
        },
        {
          key: "REINDEX_OPERATOR",
          label: "Reindexing operator",
          visible: visibility.isReindexSpace,
        },
        {
          key: "CHOOSE_LAUEGROUP",
          label: "Lauegroup to choose",
          visible: visibility.isChooseLauegroup,
        },
      ],
      referenceOptions: [
        {
          key: "REFERENCE_FOR_AIMLESS",
          label: "Reference",
          visible: visibility.isMatchMode,
        },
        {
          key: "REFERENCE_DATASET",
          label: "Reference type",
          visible: visibility.hasAimlessRef,
        },
        {
          key: "HKLIN_REF",
          label: "Reference reflections",
          visible: visibility.isHklReference,
        },
        {
          key: "XYZIN_REF",
          label: "Reference coordinates",
          visible: visibility.isXyzReference,
        },
      ],
    }),
    [visibility]
  );

  // Render helpers
  const renderElements = useCallback(
    (elements: typeof elementConfigs.fileInputs) =>
      elements.map(({ key, label }) => (
        <CCP4i2TaskElement
          {...props}
          key={key}
          itemName={key}
          qualifiers={{ guiLabel: label }}
        />
      )),
    [props]
  );

  const renderConditionalElements = useCallback(
    (elements: typeof elementConfigs.choiceOptions) =>
      elements.map(({ key, label, visible }) => (
        <CCP4i2TaskElement
          {...props}
          key={key}
          itemName={key}
          qualifiers={{ guiLabel: label }}
          visibility={visible}
        />
      )),
    [props]
  );

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs" key="1">
        <CCP4i2ContainerElement
          {...props}
          key="Files"
          itemName=""
          containerHint="BlockLevel"
          qualifiers={{ initiallyOpen: true, guiLabel: "File inputs" }}
        >
          {renderElements(elementConfigs.fileInputs)}
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          key="Parameters"
          itemName=""
          containerHint="BlockLevel"
          qualifiers={{ guiLabel: "Parameters", initiallyOpen: true }}
        >
          {renderElements(elementConfigs.parameters)}
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...props}
          key="ChoosingSpace"
          itemName=""
          containerHint="BlockLevel"
          qualifiers={{ guiLabel: "Choosing spacegroup", initiallyOpen: true }}
        >
          <CCP4i2TaskElement
            {...props}
            key="MODE"
            itemName="MODE"
            qualifiers={{ guiLabel: "Pipeline mode" }}
          />

          <CCP4i2ContainerElement
            {...props}
            key="ChoiceOptions"
            itemName=""
            containerHint="BlockLevel"
            qualifiers={{ guiLabel: "Choice options" }}
            visibility={visibility.isChooseMode}
          >
            {renderConditionalElements(elementConfigs.choiceOptions)}
          </CCP4i2ContainerElement>

          <CCP4i2ContainerElement
            {...props}
            key="SpecifyReference"
            itemName=""
            containerHint="BlockLevel"
            qualifiers={{ guiLabel: "Specify reference" }}
            visibility={visibility.isMatchMode}
          >
            {renderConditionalElements(elementConfigs.referenceOptions)}
          </CCP4i2ContainerElement>
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
