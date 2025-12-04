import React, { useCallback, useMemo, useState } from "react";
import {
  Card,
  CardContent,
  CardHeader,
  Collapse,
  FormControl,
  FormLabel,
  Grid2,
  IconButton,
  Stack,
  Typography,
} from "@mui/material";
import { Add, Delete, ExpandMore as ExpandMoreIcon } from "@mui/icons-material";

import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { useApi } from "../../../api";
import { useJob, useProject, valueOfItem } from "../../../utils";
import { MyExpandMore } from "../../expand-more";
import { useCCP4i2Window } from "../../../app-context";
import { Project } from "../../../types/models";
import { ErrorTrigger } from "./error-info";

// Types
interface CListElementProps extends CCP4i2TaskElementProps {
  initiallyOpen?: boolean;
}

interface SetParameterArg {
  object_path: string;
  value: any;
}

interface SetParameterResponse {
  success: boolean;
  data?: {
    updated_item: any;
  };
  error?: string;
}

// Constants
const DEFAULT_VALUES = {
  CAltSpaceGroup: "P1",
  CInt: 0,
  CFloat: 0.0,
  CString: "",
} as const;

const GRID_SIZES = {
  ELEMENT: { xs: 11 },
  DELETE_BUTTON: { xs: 1 },
} as const;

// Custom hooks
const useListState = (initiallyOpen = true) => {
  const [expanded, setExpanded] = useState(initiallyOpen);

  const toggleExpanded = useCallback(() => {
    setExpanded((prev) => !prev);
  }, []);

  return {
    expanded,
    toggleExpanded,
  };
};

const useGuiLabel = (item: any, qualifiers: any): string => {
  return useMemo(() => {
    return (
      qualifiers?.guiLabel ||
      item?._objectPath?.split(".").at(-1) ||
      "Unnamed List"
    );
  }, [item, qualifiers]);
};

const useVisibility = (visibility: any): boolean => {
  return useMemo(() => {
    if (!visibility) return true;
    if (typeof visibility === "function") {
      try {
        return visibility();
      } catch (error) {
        console.error("Error evaluating visibility function:", error);
        return true; // Default to visible on error
      }
    }
    return Boolean(visibility);
  }, [visibility]);
};

const useValidationBorderColor = (
  itemName: string,
  item: any,
  getValidationColor: (item: any) => string
): string => {
  return useMemo(() => {
    if (itemName && item) {
      return getValidationColor(item);
    }
    return "divider"; // Default border color
  }, [itemName, item, getValidationColor]);
};

// Utility functions
const createNewItemValue = (
  taskElement: any,
  project: Project | undefined
): any => {
  let newItemValue = valueOfItem(taskElement);

  if (taskElement._baseClass === "CDataFile" && newItemValue && project) {
    return {
      ...newItemValue,
      project: project.uuid.replace(/-/g, ""),
      baseName: "UNDEFINED",
    };
  }

  // Use lookup table for default values
  const defaultValue =
    DEFAULT_VALUES[taskElement._class as keyof typeof DEFAULT_VALUES] ||
    DEFAULT_VALUES[taskElement._baseClass as keyof typeof DEFAULT_VALUES];

  return defaultValue !== undefined ? defaultValue : newItemValue;
};

const updateObjectPath = (element: any, newIndex: number): any => {
  if (!element) return element;

  const updatedElement = { ...element };
  updatedElement._objectPath = element._objectPath?.replace(
    "[?]",
    `[${newIndex}]`
  );

  // Update nested value elements if they exist
  if (typeof updatedElement._value === "object" && updatedElement._value) {
    updatedElement._value = Object.keys(updatedElement._value).reduce(
      (acc, key) => {
        const valueElement = updatedElement._value[key];
        if (valueElement?._objectPath) {
          acc[key] = {
            ...valueElement,
            _objectPath: valueElement._objectPath.replace(
              "[?]",
              `[${newIndex}]`
            ),
          };
        } else {
          acc[key] = valueElement;
        }
        return acc;
      },
      {} as any
    );
  }

  return updatedElement;
};

const validateListOperation = (item: any, job: any): boolean => {
  if (!item) {
    console.error("List operation failed: No item provided");
    return false;
  }

  if (job.status !== 1) {
    console.warn("List operation blocked: Job is not in editable state");
    return false;
  }

  return true;
};

// Main component
export const CListElement: React.FC<CListElementProps> = ({
  itemName,
  job,
  qualifiers,
  initiallyOpen = true,
  onChange,
  visibility,
  ...restProps
}) => {
  const api = useApi();
  const { useTaskItem, setParameter, getValidationColor } = useJob(job.id);
  const { projectId } = useCCP4i2Window();
  const { project } = projectId
    ? useProject(projectId)
    : { project: undefined };
  const { item } = useTaskItem(itemName);

  // Custom hooks
  const { expanded, toggleExpanded } = useListState(initiallyOpen);
  const guiLabel = useGuiLabel(item, qualifiers);
  const isVisible = useVisibility(visibility);
  const validationBorderColor = useValidationBorderColor(
    itemName,
    item,
    getValidationColor
  );

  // Computed values
  const isEditable = useMemo(() => job.status === 1, [job.status]);
  const listItems = useMemo(() => item?._value || [], [item]);
  const hasItems = useMemo(() => listItems.length > 0, [listItems]);

  // Card styling with validation color
  const cardSx = useMemo(
    () => ({
      mx: 2,
      border: 2,
      borderColor: validationBorderColor,
      "&:hover": {
        borderColor:
          validationBorderColor === "divider"
            ? "primary.light"
            : validationBorderColor,
      },
    }),
    [validationBorderColor]
  );

  // Event handlers
  const handleAddItem = useCallback(async () => {
    if (!validateListOperation(item, job)) return;

    try {
      // Deep clone the sub-item template
      const taskElement = JSON.parse(JSON.stringify(item._subItem));
      const newIndex = listItems.length;

      // Update object paths
      const updatedTaskElement = updateObjectPath(taskElement, newIndex);

      // Get current list value
      const currentListValue = Array.isArray(valueOfItem(item))
        ? [...valueOfItem(item)]
        : [];

      // Create new item value
      const newItemValue = createNewItemValue(updatedTaskElement, project);
      currentListValue.push(newItemValue);

      console.log("Adding list item:", {
        objectPath: updatedTaskElement._objectPath,
        newValue: currentListValue,
        newItemValue,
      });

      // Set parameter
      const setParameterArg: SetParameterArg = {
        object_path: item._objectPath,
        value: currentListValue,
      };

      const result = (await setParameter(
        setParameterArg
      )) as SetParameterResponse;

      if (result?.success && result.data?.updated_item && onChange) {
        await onChange(result.data.updated_item);
      } else if (result && !result.success) {
        console.error("Failed to add list item:", result.error);
      }
    } catch (error) {
      console.error("Error adding list item:", error);
    }
  }, [item, project, job, listItems.length, setParameter, onChange]);

  const handleDeleteItem = useCallback(
    async (deletedItem: any) => {
      if (!validateListOperation(item, job)) return;

      try {
        const currentArray = [...item._value];
        const itemIndex = currentArray.findIndex(
          (arrayItem) => arrayItem === deletedItem
        );

        if (itemIndex === -1) {
          console.warn("Item not found in array for deletion");
          return;
        }

        // Remove item from array
        currentArray.splice(itemIndex, 1);

        console.log("Deleting list item:", {
          index: itemIndex,
          deletedItem,
          remainingItems: currentArray.length,
        });

        // Set parameter
        const setParameterArg: SetParameterArg = {
          object_path: item._objectPath,
          value: currentArray,
        };

        const result = (await setParameter(
          setParameterArg
        )) as SetParameterResponse;

        if (result?.success && result.data?.updated_item && onChange) {
          await onChange(result.data.updated_item);
        } else if (result && !result.success) {
          console.error("Failed to delete list item:", result.error);
        }
      } catch (error) {
        console.error("Error deleting list item:", error);
      }
    },
    [item, job, setParameter, onChange]
  );

  const handleExpandToggle = useCallback(
    (event: React.MouseEvent) => {
      event.stopPropagation();
      toggleExpanded();
    },
    [toggleExpanded]
  );

  // Render helpers
  const renderListActions = () => (
    <Stack direction="row" alignItems="center">
      <MyExpandMore
        sx={{ color: "primary.text" }}
        expand={expanded}
        onClick={handleExpandToggle}
        aria-expanded={expanded}
        aria-label={expanded ? "Collapse list" : "Expand list"}
      >
        <ExpandMoreIcon />
      </MyExpandMore>

      <IconButton
        disabled={!isEditable}
        onClick={handleAddItem}
        size="small"
        sx={{ color: "primary.text" }}
        aria-label="Add new item to list"
      >
        <Add />
      </IconButton>

      <ErrorTrigger item={item} job={job} />
    </Stack>
  );

  const renderEmptyState = () => (
    <Typography
      variant="caption"
      color="text.secondary"
      sx={{
        fontStyle: "italic",
        textAlign: "center",
        py: 2,
      }}
    >
      No elements in this list
    </Typography>
  );

  const renderListItem = (content: any, index: number) => (
    <Grid2
      key={content._objectPath || `item-${index}`}
      container
      spacing={1}
      sx={{ mb: 1 }}
    >
      <Grid2 size={GRID_SIZES.ELEMENT}>
        <CCP4i2TaskElement
          {...restProps}
          itemName={content._objectPath}
          job={job}
          qualifiers={qualifiers}
          onChange={onChange}
        />
      </Grid2>

      <Grid2 size={GRID_SIZES.DELETE_BUTTON}>
        <FormControl fullWidth>
          <FormLabel
            component="legend"
            sx={{
              fontSize: "0.75rem",
              mb: 0.5,
            }}
          >
            Delete
          </FormLabel>
          <IconButton
            disabled={!isEditable}
            onClick={() => handleDeleteItem(content)}
            size="small"
            color="error"
            aria-label={`Delete item ${index + 1}`}
            sx={{
              "&:hover": {
                backgroundColor: "error.lighter",
              },
            }}
          >
            <Delete />
          </IconButton>
        </FormControl>
      </Grid2>
    </Grid2>
  );

  const renderListContent = () => (
    <Collapse in={expanded} timeout="auto" unmountOnExit>
      {!hasItems ? (
        renderEmptyState()
      ) : (
        <Stack spacing={1}>
          {listItems.map((content: any, index: number) =>
            renderListItem(content, index)
          )}
        </Stack>
      )}
    </Collapse>
  );

  // Early return if not visible
  if (!isVisible) {
    return null;
  }

  return (
    <Card sx={cardSx}>
      <CardHeader
        sx={{
          backgroundColor: "background.paper", // Very light grey background
          color: "text.primary", // Black text (uses theme's primary text color)
          "& .MuiCardHeader-title": {
            color: "text.primary", // Ensure title is black
          },
          "&:hover": {
            backgroundColor: "background.paper", // Slightly darker on hover
          },
        }}
        title={
          <Typography variant="body2" component="div">
            {guiLabel}
          </Typography>
        }
        action={renderListActions()}
      />

      <CardContent sx={{ pt: 1 }}>{renderListContent()}</CardContent>
    </Card>
  );
};

// Set display name for debugging
CListElement.displayName = "CListElement";

export default CListElement;
