/*
 * Copyright (C) 2026 Newcastle University
 *
 * This file is part of CCP4i2.
 *
 * CCP4i2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License version 3,
 * modified in accordance with the provisions of the license to address
 * the requirements of UK law.
 *
 * See https://www.ccp4.ac.uk/ccp4license.php for details.
 */
import React from "react";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { isTruthy } from "../../task-elements/shared-hooks";

interface InputDataTabProps extends CCP4i2TaskInterfaceProps {
  handleF_SIGFChange: () => Promise<void>;
  refinementMode: string;
  hydrUse: any;
  addWaters: any;
  useAnomalous: any;
  F_SIGFItem: any;
}

export const InputDataTab: React.FC<InputDataTabProps> = (props) => {
  const {
    handleF_SIGFChange,
    refinementMode,
    hydrUse,
    addWaters,
    useAnomalous,
    F_SIGFItem,
    ...taskProps
  } = props;

  return (
    <>
      {/* Main inputs */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Main inputs" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="container.inputData.XYZIN" />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="F_SIGF"
          onChange={handleF_SIGFChange}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="WAVELENGTH"
          qualifiers={{ guiLabel: "Wavelength" }}
        />
        <CCP4i2TaskElement {...taskProps} itemName="FREERFLAG" />
      </CCP4i2ContainerElement>

      {/* Experimental phase information */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Experimental phase information" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="ABCD" />
      </CCP4i2ContainerElement>

      {/* Additional geometry dictionaries */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Additional geometry dictionaries" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement {...taskProps} itemName="DICT_LIST" />
      </CCP4i2ContainerElement>

      {/* Options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="REFINEMENT_MODE"
          qualifiers={{
            guiLabel: "Refinement mode (restrained or rigid body):",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="NCYCLES"
          qualifiers={{
            guiLabel: "Number of restrained refinement cycles:",
          }}
          visibility={() => refinementMode !== "RIGID"}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="NCYCRIGID"
          qualifiers={{
            guiLabel: "Number of rigid body refinement cycles:",
          }}
          visibility={() => refinementMode === "RIGID"}
        />

        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HYDR_USE"
            qualifiers={{
              guiLabel: "Use riding hydrogens during refinement",
            }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="HYDR_ALL"
            visibility={() => isTruthy(hydrUse)}
          />
        </FieldRow>

        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="ADD_WATERS"
            qualifiers={{ guiLabel: "Add waters" }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="REFPRO_RSR_RWORK_LIMIT"
            qualifiers={{ guiLabel: "or lower" }}
            visibility={() => isTruthy(addWaters)}
          />
        </FieldRow>

        <CCP4i2TaskElement
          {...taskProps}
          itemName="USE_TWIN"
          qualifiers={{ guiLabel: "Crystal is twinned" }}
        />
      </CCP4i2ContainerElement>

      {/* Anomalous signal */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Anomalous signal" }}
        containerHint="FolderLevel"
        visibility={() =>
          F_SIGFItem?.contentFlag && [1, 2].includes(F_SIGFItem.contentFlag)
        }
      >
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="USEANOMALOUS"
            qualifiers={{ guiLabel: "Use anomalous" }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="USEANOMALOUSFOR"
            qualifiers={{ guiLabel: "Use for" }}
            visibility={() => isTruthy(useAnomalous)}
          />
        </FieldRow>
      </CCP4i2ContainerElement>
    </>
  );
};
