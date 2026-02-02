import { useCallback, useEffect, useMemo } from "react";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useRunCheck } from "../../../providers/run-check-provider";
import { FieldRow } from "../task-elements/field-row";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem, validation } = useJob(job.id);
  const { processedErrors, setProcessedErrors } = useRunCheck();

  // Get syncTo for REFERENCE_FOR_AIMLESS to auto-sync it with mode
  const { syncTo: syncToAimlessRef } = useTaskItem("REFERENCE_FOR_AIMLESS");

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

  // Auto-set REFERENCE_FOR_AIMLESS based on MODE:
  // - true when MODE is "MATCH" (always provide reference when matching)
  // - false otherwise
  // Uses syncTo to prevent bouncing loops (handles pending update tracking internally)
  useEffect(() => {
    syncToAimlessRef(taskValues.mode === "MATCH");
  }, [taskValues.mode, syncToAimlessRef]);

  // Process validation errors
  // TODO: This client-side validity filtering should be moved to Python's validity() method
  // in pipelines/aimless_pipe/script/aimless_pipe.py to avoid race conditions and GUI complexity.
  // The Python validity() method can adjust qualifiers (e.g., allowUndefined) before validation runs.
  const processedValidationErrors = useMemo(() => {
    if (!validation) return null;

    const filtered = Object.keys(validation)
      .filter((key) => !key.startsWith("aimless_pipe.container.controlParameters.CELL."))
      .reduce((acc, key) => ({ ...acc, [key]: validation[key] }), {});

    // Add custom validation for HKLIN_REF
    // Note: aimlessRef check removed since it's now auto-synced to true in MATCH mode
    if (
      taskValues.mode === "MATCH" &&
      taskValues.referenceDataset === "HKL" &&
      !taskValues.hklinRef?.dbFileId
    ) {
      filtered["aimless_pipe.container.inputData.HKLIN_REF"] = {
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

  // Visibility and disabled helpers
  // Note: REFERENCE_FOR_AIMLESS is auto-synced with MODE, so we no longer need
  // to check taskValues.aimlessRef - it's always true when mode is "MATCH"
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
      isHklReference: () =>
        taskValues.mode === "MATCH" && taskValues.referenceDataset === "HKL",
      isXyzReference: () =>
        taskValues.mode === "MATCH" && taskValues.referenceDataset === "XYZ",
    }),
    [taskValues]
  );

  // Disabled helper - reference type selector is visible but disabled unless mode is MATCH
  const isReferenceTypeDisabled = useMemo(
    () => () => taskValues.mode !== "MATCH",
    [taskValues.mode]
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
      // Note: REFERENCE_FOR_AIMLESS toggle is auto-synced with MODE and hidden from UI
      // Reference file selectors are rendered directly in the JSX (not through helpers)
      // to support the side-by-side MODE + REFERENCE_DATASET layout with disabled state
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
          {/* MODE and REFERENCE_DATASET side-by-side:
              - MODE is always editable
              - REFERENCE_DATASET is always visible but disabled unless mode is MATCH */}
          <FieldRow>
            <CCP4i2TaskElement
              {...props}
              key="MODE"
              itemName="MODE"
              qualifiers={{ guiLabel: "Pipeline mode" }}
            />
            <CCP4i2TaskElement
              {...props}
              key="REFERENCE_DATASET"
              itemName="REFERENCE_DATASET"
              qualifiers={{ guiLabel: "Reference type" }}
              disabled={isReferenceTypeDisabled}
            />
          </FieldRow>

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

          {/* Reference file selectors - only shown when mode is MATCH */}
          <CCP4i2TaskElement
            {...props}
            key="HKLIN_REF"
            itemName="HKLIN_REF"
            qualifiers={{ guiLabel: "Reference reflections" }}
            visibility={visibility.isHklReference}
          />
          <CCP4i2TaskElement
            {...props}
            key="XYZIN_REF"
            itemName="XYZIN_REF"
            qualifiers={{ guiLabel: "Reference coordinates" }}
            visibility={visibility.isXyzReference}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};

export default TaskInterface;
