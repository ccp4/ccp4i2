import React, { useMemo } from "react";

import { CCP4i2TaskElementProps } from "./task-element";
import { CSimpleTextFieldElement } from "./csimple-textfield";
import { CSimpleAutocompleteElement } from "./csimple-autocomplete";
import { useJob } from "../../../utils";
import { FieldSize } from "./field-sizes";

// Types
export interface CCP4i2CSimpleElementProps extends CCP4i2TaskElementProps {
  type: "int" | "float" | "text" | "checkbox";
  /** Optional explicit field size. If not provided, size is inferred from item type/qualifiers. */
  size?: FieldSize;
}

interface ProcessedQualifiers {
  onlyEnumerators?: boolean;
  enumerators?: any[];
  [key: string]: any;
}

interface ComponentConfig {
  shouldUseSelect: boolean;
  isVisible: boolean;
  qualifiers: ProcessedQualifiers | null;
}

// Constants
const INPUT_TYPES = {
  INT: "int",
  FLOAT: "float",
  TEXT: "text",
  CHECKBOX: "checkbox",
} as const;

// Custom hooks
const useProcessedQualifiers = (
  item: any,
  qualifiers: any,
  itemName: string
): ProcessedQualifiers | null => {
  return useMemo(() => {
    // If item has qualifiers, merge with props qualifiers
    if (item?._qualifiers) {
      try {
        return qualifiers
          ? { ...item._qualifiers, ...qualifiers }
          : item._qualifiers;
      } catch (error) {
        console.error(`Error processing qualifiers for ${itemName}:`, error);
        return qualifiers || null;
      }
    }

    // Return props qualifiers if no item qualifiers
    return qualifiers || null;
  }, [item, qualifiers, itemName]);
};

const useComponentConfig = (
  processedQualifiers: ProcessedQualifiers | null,
  visibility: any
): ComponentConfig => {
  return useMemo(() => {
    // Determine if we should use select/autocomplete component
    const shouldUseSelect = Boolean(
      processedQualifiers?.onlyEnumerators ||
        (Array.isArray(processedQualifiers?.enumerators) &&
          processedQualifiers.enumerators.length > 0)
    );

    // Determine visibility
    let isVisible = true;
    if (visibility === false) {
      isVisible = false;
    } else if (typeof visibility === "function") {
      try {
        isVisible = visibility();
      } catch (error) {
        console.error("Error evaluating visibility function:", error);
        isVisible = true; // Default to visible on error
      }
    }

    return {
      shouldUseSelect,
      isVisible,
      qualifiers: processedQualifiers,
    };
  }, [processedQualifiers, visibility]);
};

// Validation functions
const validateProps = (props: CCP4i2CSimpleElementProps): boolean => {
  const { itemName, job, type } = props;

  if (!itemName || typeof itemName !== "string") {
    console.error("CSimpleElement: itemName is required and must be a string");
    return false;
  }

  if (!job || typeof job.id !== "number") {
    console.error("CSimpleElement: job with valid id is required");
    return false;
  }

  if (!Object.values(INPUT_TYPES).includes(type)) {
    console.error(
      `CSimpleElement: invalid type "${type}". Must be one of: ${Object.values(
        INPUT_TYPES
      ).join(", ")}`
    );
    return false;
  }

  return true;
};

// Main component
export const CSimpleElement: React.FC<CCP4i2CSimpleElementProps> = (props) => {
  const { itemName, job, qualifiers, visibility, ...restProps } = props;

  // Early validation
  if (!validateProps(props)) {
    return null;
  }

  // Get task item
  const { useTaskItem } = useJob(job.id);
  const { item } = useTaskItem(itemName);

  // Process qualifiers
  const processedQualifiers = useProcessedQualifiers(
    item,
    qualifiers,
    itemName
  );

  // Get component configuration
  const config = useComponentConfig(processedQualifiers, visibility);

  // Early return if not visible
  if (!config.isVisible) {
    return null;
  }

  // Common props for both components
  const commonProps = {
    ...restProps,
    itemName,
    job,
    qualifiers: config.qualifiers,
    visibility,
  };

  // Render appropriate component based on configuration
  if (config.shouldUseSelect) {
    return (
      <CSimpleAutocompleteElement
        {...commonProps}
        key={`${itemName}-autocomplete`}
      />
    );
  }

  return (
    <CSimpleTextFieldElement {...commonProps} key={`${itemName}-textfield`} />
  );
};

// Default export with display name for debugging
CSimpleElement.displayName = "CSimpleElement";

export default CSimpleElement;
