import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import {
  Box,
  Card,
  CardContent,
  Chip,
  FormControlLabel,
  Stack,
  Switch,
  Typography,
} from "@mui/material";
import { useJob } from "../../../utils";
import { ErrorInfo } from "./error-info";
import { useMemo } from "react";
import TableChartIcon from "@mui/icons-material/TableChart";

// Color mapping for column group types
const TYPE_COLORS: Record<string, "primary" | "secondary" | "success" | "warning" | "info"> = {
  Obs: "primary",
  Phs: "secondary",
  MapCoeffs: "success",
  FreeR: "warning",
};

// Icons/labels for column group types
const TYPE_LABELS: Record<string, string> = {
  Obs: "Observations",
  Phs: "Phases",
  MapCoeffs: "Map Coefficients",
  FreeR: "Free R Flag",
};

/**
 * CColumnGroupElement - A compact widget for displaying MTZ column groups
 *
 * Displays:
 * - Column group type as a colored chip
 * - Dataset name
 * - List of columns in the group
 * - Selection checkbox
 */
export const CColumnGroupElement: React.FC<CCP4i2TaskElementProps> = (props) => {
  const { itemName, job } = props;
  const { useTaskItem, getValidationColor, setParameter } = useJob(job.id);

  const { item } = useTaskItem(itemName);

  // Extract values from the item
  const columnGroupType = useMemo(() => {
    const typeValue = item?._value?.columnGroupType;
    if (typeof typeValue === "object" && typeValue?._value !== undefined) {
      return typeValue._value;
    }
    return typeValue || "Unknown";
  }, [item]);

  const dataset = useMemo(() => {
    const datasetValue = item?._value?.dataset;
    if (typeof datasetValue === "object" && datasetValue?._value !== undefined) {
      return datasetValue._value;
    }
    return datasetValue || "";
  }, [item]);

  const columnList = useMemo(() => {
    const listValue = item?._value?.columnList;
    if (Array.isArray(listValue?._value)) {
      return listValue._value.map((col: any) =>
        typeof col === "object" ? col._value : col
      );
    }
    if (Array.isArray(listValue)) {
      return listValue.map((col: any) =>
        typeof col === "object" ? col._value : col
      );
    }
    return [];
  }, [item]);

  const selected = useMemo(() => {
    const selectedValue = item?._value?.selected;
    if (typeof selectedValue === "object" && selectedValue?._value !== undefined) {
      return Boolean(selectedValue._value);
    }
    return Boolean(selectedValue);
  }, [item]);

  const isEditable = job.status === 1;

  const handleSelectionChange = async (event: React.ChangeEvent<HTMLInputElement>) => {
    if (!item || !isEditable) return;

    await setParameter({
      object_path: `${item._objectPath}.selected`,
      value: event.target.checked,
    });
  };

  const chipColor = TYPE_COLORS[columnGroupType] || "info";
  const typeLabel = TYPE_LABELS[columnGroupType] || columnGroupType;
  const validationColor = getValidationColor(item);

  if (!item) {
    return <Typography color="text.secondary">Loading...</Typography>;
  }

  return (
    <Card
      variant="outlined"
      sx={{
        borderColor: validationColor,
        borderWidth: 2,
        bgcolor: selected ? "action.selected" : "background.paper",
        transition: "all 0.2s ease-in-out",
        "&:hover": {
          boxShadow: 2,
        },
      }}
    >
      <CardContent sx={{ py: 1.5, px: 2, "&:last-child": { pb: 1.5 } }}>
        <Stack spacing={1}>
          {/* Header row: Type chip, dataset, selection toggle */}
          <Stack
            direction="row"
            alignItems="center"
            justifyContent="space-between"
            spacing={1}
          >
            <Stack direction="row" alignItems="center" spacing={1} sx={{ flex: 1 }}>
              <Chip
                icon={<TableChartIcon />}
                label={typeLabel}
                color={chipColor}
                size="small"
                sx={{ fontWeight: 600 }}
              />
              {dataset && (
                <Typography
                  variant="body2"
                  color="text.secondary"
                  sx={{
                    overflow: "hidden",
                    textOverflow: "ellipsis",
                    whiteSpace: "nowrap",
                  }}
                >
                  {dataset}
                </Typography>
              )}
            </Stack>

            <Stack direction="row" alignItems="center" spacing={1}>
              <ErrorInfo {...props} />
              <FormControlLabel
                control={
                  <Switch
                    checked={selected}
                    onChange={handleSelectionChange}
                    disabled={!isEditable}
                    size="small"
                    color="primary"
                  />
                }
                label=""
                sx={{ m: 0 }}
              />
            </Stack>
          </Stack>

          {/* Column list */}
          {columnList.length > 0 && (
            <Box
              sx={{
                display: "flex",
                flexWrap: "wrap",
                gap: 0.5,
              }}
            >
              {columnList.map((col: string, index: number) => (
                <Chip
                  key={`${col}-${index}`}
                  label={col}
                  size="small"
                  variant="outlined"
                  sx={{
                    height: 22,
                    fontSize: "0.75rem",
                    fontFamily: "monospace",
                    bgcolor: "background.default",
                  }}
                />
              ))}
            </Box>
          )}
        </Stack>
      </CardContent>
    </Card>
  );
};

export default CColumnGroupElement;
