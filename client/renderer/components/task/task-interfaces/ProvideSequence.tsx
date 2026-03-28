/*
 * Copyright (C) 2025-2026 Newcastle University
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
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { useJob } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback } from "react";
import { apiText } from "../../../api-fetch";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  const { forceUpdate: forceSetSEQUENCETEXT } = useTaskItem("SEQUENCETEXT");

  // Handle SEQIN file change - read raw file content and populate SEQUENCETEXT
  // Matches Qt behavior: open file, read content, set SEQUENCETEXT
  const handleSeqInChange = useCallback(
    async (updatedItem: any) => {
      // dbFileId is a CData object with nested _value, not a plain string
      const dbFileId =
        updatedItem?._value?.dbFileId?._value?.trim() ||
        updatedItem?.dbFileId;
      if (!dbFileId) return;

      try {
        const content = await apiText(`files/${dbFileId}/download_by_uuid`);
        if (content) {
          await forceSetSEQUENCETEXT(content);
        }
      } catch (error) {
        console.error("Error reading sequence file:", error);
      }
    },
    [forceSetSEQUENCETEXT]
  );

  return (
    <CCP4i2Tabs {...props}>
      <CCP4i2Tab label="Main inputs">
        <CCP4i2ContainerElement
          {...props}
          itemName=""
          qualifiers={{ guiLabel: "Key files" }}
          containerHint="FolderLevel"
          initiallyOpen={true}
        >
          <CCP4i2TaskElement
            {...props}
            itemName="SEQUENCETEXT"
            qualifiers={{ guiLabel: "Sequence" }}
            sx={{ minWidth: "100%", minHeight: "10rem" }}
          />

          <CCP4i2TaskElement
            {...props}
            itemName="SEQIN"
            qualifiers={{ guiLabel: "File from which to extract sequence" }}
            onChange={handleSeqInChange}
          />

          <CCP4i2TaskElement
            {...props}
            itemName="XYZIN"
            qualifiers={{ guiLabel: "Coordinate file (for extracting sequence or Matthews calc)" }}
          />
        </CCP4i2ContainerElement>
      </CCP4i2Tab>
    </CCP4i2Tabs>
  );
};
export default TaskInterface;
