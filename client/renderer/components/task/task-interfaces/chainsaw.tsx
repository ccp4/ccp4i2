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
import { LinearProgress, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";
import { useCallback, useMemo } from "react";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { container, useTaskItem, fetchDigest } = useJob(props.job.id);

  const { item: alignItem } = useTaskItem("ALIGNIN");
  const { item: targetItem, forceUpdate: setTargetIndex } =
    useTaskItem("TARGETINDEX");

  // Extract current digest-populated enumerators from ALIGNIN
  // useFileDigest watches the ALIGNIN object path and re-fetches on file change
  const { data: alignDigest } = useJob(props.job.id).useFileDigest(
    alignItem?._objectPath || ""
  );

  // Build enumerators and labels from the alignment digest identifiers
  const { enumerators, menuText } = useMemo(() => {
    const ids: string[] = alignDigest?.identifiers || [];
    if (ids.length === 0) return { enumerators: [], menuText: [] };
    return {
      enumerators: ids.map((_: string, i: number) => i),
      menuText: ids,
    };
  }, [alignDigest?.identifiers]);

  // When ALIGNIN changes, fetch its digest to update identifiers
  const handleAlignInChange = useCallback(async () => {
    if (!alignItem?._objectPath) return;
    const digest = await fetchDigest(alignItem._objectPath);
    // Reset TARGETINDEX to 0 (first sequence) when alignment changes
    if (digest?.identifiers?.length > 0) {
      await setTargetIndex(0);
    }
  }, [alignItem?._objectPath, fetchDigest, setTargetIndex]);

  if (!container) return <LinearProgress />;

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      {/* Input Data */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Input Data" }}
        containerHint="FolderLevel"
      >
        <Typography variant="subtitle2">
          Enter structure to be edited
        </Typography>
        <CCP4i2TaskElement itemName="XYZIN" {...props} />
        <Typography variant="subtitle2">
          Enter sequence alignment of edited structure and MR target
        </Typography>
        <CCP4i2TaskElement
          itemName="ALIGNIN"
          {...props}
          onChange={handleAlignInChange}
        />
        {enumerators.length > 0 && (
          <CCP4i2TaskElement
            itemName="TARGETINDEX"
            {...props}
            qualifiers={{
              guiLabel: "Identifier of the target sequence in this alignment",
              enumerators,
              menuText,
              onlyEnumerators: true,
            }}
          />
        )}
      </CCP4i2ContainerElement>

      {/* Simple Options */}
      <CCP4i2ContainerElement
        {...props}
        itemName=""
        qualifiers={{ guiLabel: "Simple Options" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2">
          How severe is side-chain truncation for non-conserved residues:
        </Typography>
        <CCP4i2TaskElement itemName="MODE" {...props} />
      </CCP4i2ContainerElement>
    </Paper>
  );
};

export default TaskInterface;
