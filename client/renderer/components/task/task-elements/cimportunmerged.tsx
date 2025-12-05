import { useCallback, useMemo } from "react";
import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useJob, valueOfItem } from "../../../utils";
import { apiGet } from "../../../api-fetch";
import { Grid2, Typography } from "@mui/material";
import { CSimpleDataFileElement } from "./csimpledatafile";

export const CImportUnmergedElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const { itemName, job } = props;
  const {
    useTaskItem,
    useFileDigest,
    getValidationColor,
    setParameterNoMutate,
    mutateContainer,
    mutateValidation,
    mutateParams_xml,
  } = useJob(job.id);

  const { item } = useTaskItem(itemName);
  // Get current values for comparison (no update functions needed - we use sparse dict)
  const { value: cell } = useTaskItem(`${itemName}.cell`);
  const { value: wavelength } = useTaskItem(`${itemName}.wavelength`);
  const { value: crystalName } = useTaskItem(`${itemName}.crystalName`);
  const { value: dataset } = useTaskItem(`${itemName}.dataset`);

  const fileObjectPath = useMemo(
    () => (item?._objectPath ? `${item._objectPath}.file` : null),
    [item]
  );
  // fileDigest used only for declarative presentation (batch display)
  const { data: fileDigest } = useFileDigest(fileObjectPath || "");

  /**
   * Imperative handler: only called when user explicitly changes the file.
   * Fetches fresh digest and updates all changed fields in a single API call
   * using a sparse dict (only includes fields that actually changed).
   */
  const handleChange = useCallback(
    async (updated: any) => {
      if (!item || !setParameterNoMutate || !updated) return;

      const updatedValue = valueOfItem(updated);
      const digestResponse = await apiGet(
        `files/${updatedValue.dbFileId}/digest_by_uuid`
      );
      const digestData = digestResponse?.data;
      if (!digestData) return;

      // Build sparse update dict - only include fields that changed
      const updates: Record<string, any> = {};

      if (
        digestData.cell &&
        JSON.stringify(digestData.cell) !== JSON.stringify(cell)
      ) {
        updates.cell = digestData.cell;
      }
      if (digestData.wavelength && digestData.wavelength !== wavelength) {
        updates.wavelength = digestData.wavelength;
      }
      if (digestData.crystalName && digestData.crystalName !== crystalName) {
        updates.crystalName = digestData.crystalName;
      }
      if (digestData.datasetName && digestData.datasetName !== dataset) {
        updates.dataset = digestData.datasetName;
      }

      // Single API call with sparse dict (server's update() handles it)
      if (Object.keys(updates).length > 0) {
        await setParameterNoMutate({
          object_path: item._objectPath,
          value: updates,
        });
        await Promise.all([
          mutateContainer(),
          mutateValidation(),
          mutateParams_xml(),
        ]);
      }
    },
    [
      item,
      setParameterNoMutate,
      cell,
      wavelength,
      crystalName,
      dataset,
      mutateContainer,
      mutateValidation,
      mutateParams_xml,
    ]
  );

  // Helper function for object paths
  const getObjectPath = (field: string) =>
    item ? `${item._objectPath}.${field}` : null;

  // Grid items configuration
  const gridItems = [
    { path: getObjectPath("crystalName"), label: "Crystal name" },
    { path: getObjectPath("dataset"), label: "Dataset name" },
    { path: getObjectPath("wavelength"), label: "Wavelength" },
  ];

  const inferredVisibility = useMemo(() => {
    if (!props.visibility) return true;
    return typeof props.visibility === "function"
      ? props.visibility()
      : props.visibility;
  }, [props.visibility]);

  const hasValidationError = useMemo(
    () => (item ? getValidationColor(item) === "error.light" : false),
    [getValidationColor, item]
  );

  if (!inferredVisibility || !fileObjectPath) return null;

  return (
    <CSimpleDataFileElement
      {...props}
      hasValidationError={hasValidationError}
      itemName={fileObjectPath}
      onChange={handleChange}
    >
      {getObjectPath("cell") && item._value["cell"] && (
        <CCP4i2TaskElement
          {...props}
          itemName={getObjectPath("cell")!}
          qualifiers={{ guiLabel: "Cell" }}
        />
      )}

      <Grid2 container rowSpacing={0} sx={{ mt: 2 }}>
        {gridItems.map(({ path, label }) => (
          <Grid2 key={label} size={{ xs: 4 }}>
            <CCP4i2TaskElement
              {...props}
              sx={{ my: 0, py: 0, minWidth: "10rem" }}
              itemName={path!}
              qualifiers={{
                ...props.qualifiers,
                guiLabel: label,
              }}
            />
          </Grid2>
        ))}

        <Grid2 size={{ xs: 4 }}>
          <Typography variant="body1">Batches in file</Typography>
        </Grid2>
        <Grid2 size={{ xs: 8 }}>
          <Typography variant="body1">
            {fileDigest?.batchs && JSON.stringify(fileDigest.batchs)}
          </Typography>
        </Grid2>

        <Grid2 size={{ xs: 12 }}>
          <CCP4i2TaskElement
            {...props}
            sx={{ mt: 1, mb: 0, py: 0, minWidth: "30rem" }}
            itemName={`${itemName}.excludeSelection`}
            qualifiers={{
              ...props.qualifiers,
              guiLabel: "Batch range(s) to exclude",
              guiMode: "multiLine",
            }}
          />
        </Grid2>
      </Grid2>
    </CSimpleDataFileElement>
  );
};
