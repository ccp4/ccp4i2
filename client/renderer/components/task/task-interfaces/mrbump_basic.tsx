import React, { useCallback, useEffect, useState } from "react";
import { Box, Paper, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "./task-container";
import { CCP4i2TaskElement } from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useJob } from "../../../utils";

/** Normalize CBoolean values - server may return boolean or string */
const isTruthy = (val: any): boolean =>
  val === true || val === "True" || val === "true";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const { job } = props;
  const { useTaskItem } = useJob(job.id);

  // --- Boolean toggles for conditional visibility ---
  const { value: SEARCH_PDB_RAW } = useTaskItem("SEARCH_PDB");
  const { value: SEARCH_AFDB_RAW } = useTaskItem("SEARCH_AFDB");
  const { value: INCLUDE_RAW } = useTaskItem("INCLUDE");

  const [searchPdb, setSearchPdb] = useState(() => isTruthy(SEARCH_PDB_RAW));
  const [searchAfdb, setSearchAfdb] = useState(() => isTruthy(SEARCH_AFDB_RAW));
  const [includeLocal, setIncludeLocal] = useState(() =>
    isTruthy(INCLUDE_RAW)
  );

  useEffect(() => setSearchPdb(isTruthy(SEARCH_PDB_RAW)), [SEARCH_PDB_RAW]);
  useEffect(() => setSearchAfdb(isTruthy(SEARCH_AFDB_RAW)), [SEARCH_AFDB_RAW]);
  useEffect(() => setIncludeLocal(isTruthy(INCLUDE_RAW)), [INCLUDE_RAW]);

  const handleSearchPdb = useCallback(async (item: any) => {
    setSearchPdb(isTruthy(item._value));
  }, []);

  const handleSearchAfdb = useCallback(async (item: any) => {
    setSearchAfdb(isTruthy(item._value));
  }, []);

  const handleIncludeLocal = useCallback(async (item: any) => {
    setIncludeLocal(isTruthy(item._value));
  }, []);

  return (
    <Paper sx={{ display: "flex", flexDirection: "column", gap: 1, p: 1 }}>
      <CCP4i2Tabs {...props}>
        {/* ===== Tab 1: Input Data ===== */}
        <CCP4i2Tab label="Input Data">
          {/* Target Sequence */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Target Sequence" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="ASUIN"
              {...props}
              qualifiers={{ guiLabel: "AU contents" }}
            />
            <Typography variant="body2" sx={{ fontStyle: "italic", my: 0.5 }}>
              If a suitable ASU is not available above, you can press the cross
              &amp; then button to quickly create one.
            </Typography>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                The number of monomers to search for
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="NMON"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </CCP4i2ContainerElement>

          {/* Experimental Data */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Experimental Data" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="F_SIGF"
              {...props}
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <CCP4i2TaskElement
              itemName="FREERFLAG"
              {...props}
              qualifiers={{ guiLabel: "Free R set" }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 2: Search Models ===== */}
        <CCP4i2Tab label="Search Models">
          {/* Model databases */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model databases" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="SEARCH_PDB"
              {...props}
              qualifiers={{
                guiLabel: "Search PDB for possible MR search models",
              }}
              onChange={handleSearchPdb}
            />
            {searchPdb && (
              <Box sx={{ pl: 3 }}>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                    Non-redundancy level for homologue search:
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="REDUNDANCYLEVEL"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </Box>
            )}

            <CCP4i2TaskElement
              itemName="SEARCH_AFDB"
              {...props}
              qualifiers={{
                guiLabel: "Search EBI-AFDB for possible MR search models",
              }}
              onChange={handleSearchAfdb}
            />
            {searchAfdb && (
              <Box sx={{ pl: 3 }}>
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                    EBI-AFDB pLDDT residue score cut-off:
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="AFDBLEVEL"
                      {...props}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </Box>
            )}

            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
                mt: 1,
              }}
            >
              <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                Maximum no. of search models to create:
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="MRMAX"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </CCP4i2ContainerElement>

          {/* Optional Settings */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Optional Settings" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="HHPREDIN"
              {...props}
              qualifiers={{ guiLabel: "HHPred hhr file" }}
            />
            <CCP4i2TaskElement
              itemName="PDBLOCAL"
              {...props}
              qualifiers={{ guiLabel: "Path to local PDB mirror" }}
            />
          </CCP4i2ContainerElement>

          {/* Local coordinate files */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{
              guiLabel: "Local coordinate files to be used as search models",
            }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="INCLUDE"
              {...props}
              qualifiers={{ guiLabel: "Include local files" }}
              onChange={handleIncludeLocal}
            />
            {includeLocal && (
              <CCP4i2TaskElement itemName="XYZIN_LIST" {...props} />
            )}
          </CCP4i2ContainerElement>
        </CCP4i2Tab>

        {/* ===== Tab 3: Options ===== */}
        <CCP4i2Tab label="Options">
          {/* Molecular Replacement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Molecular Replacement" }}
            containerHint="FolderLevel"
          >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                Number of cores for Phaser (maximum=10)
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="PJOBS"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </CCP4i2ContainerElement>

          {/* Refinement */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Refinement" }}
            containerHint="FolderLevel"
          >
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1" sx={{ fontStyle: "italic" }}>
                Number of refinement cycles in Refmac
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="NCYC"
                  {...props}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </CCP4i2ContainerElement>

          {/* Model Building */}
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Model Building" }}
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              itemName="BUCC"
              {...props}
              qualifiers={{
                guiLabel: "Run Buccaneer after refinement",
              }}
            />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </Paper>
  );
};

export default TaskInterface;
