import React from "react";
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { isTruthy } from "../../task-elements/shared-hooks";

interface AdvancedTabProps extends CCP4i2TaskInterfaceProps {
  isXtal: boolean;
  isSpa: boolean;
  bfacSetUse: any;
  randomizeUse: any;
  scatteringFactors: any;
  resCustom: any;
  blurUse: any;
  runAdpAnalysis: any;
  runCoordAdpDev: any;
}

export const AdvancedTab: React.FC<AdvancedTabProps> = (props) => {
  const {
    isXtal,
    isSpa,
    bfacSetUse,
    randomizeUse,
    scatteringFactors,
    resCustom,
    blurUse,
    runAdpAnalysis,
    runCoordAdpDev,
    ...taskProps
  } = props;

  return (
    <>
      {/* Resolution and experiment (xtal only) */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Resolution and experiment" }}
        containerHint="FolderLevel"
        visibility={() => isXtal}
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="RES_CUSTOM"
          qualifiers={{
            guiLabel: "Use custom resolution limits",
          }}
        />
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RES_MIN"
            qualifiers={{ guiLabel: "High resolution (dmin):" }}
            visibility={() => isTruthy(resCustom)}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RES_MAX"
            qualifiers={{ guiLabel: "Low resolution (dmax):" }}
            visibility={() => isTruthy(resCustom)}
          />
        </FieldRow>

        <CCP4i2TaskElement
          {...taskProps}
          itemName="FREERFLAG_NUMBER"
          qualifiers={{
            guiLabel: "FreeR flag number for test set:",
          }}
        />

        <CCP4i2TaskElement
          {...taskProps}
          itemName="SCATTERING_FACTORS"
          qualifiers={{ guiLabel: "Diffraction experiment type:" }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="SCATTERING_ELECTRON"
          qualifiers={{ guiLabel: "Form factor calculation:" }}
          visibility={() => scatteringFactors === "electron"}
        />

        <CCP4i2TaskElement
          {...taskProps}
          itemName="USE_WORK_IN_EST"
          qualifiers={{
            guiLabel:
              "Use work reflections in maximum likelihood parameter estimates",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Hydrogen and charge options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Hydrogen and charge options" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="H_OUT"
          qualifiers={{
            guiLabel:
              "Write hydrogen atoms in the output model",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="H_REFINE"
          qualifiers={{
            guiLabel: "Refine hydrogen positions",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="KEEP_CHARGES"
          qualifiers={{
            guiLabel:
              "Keep charges, i.e. use scattering factor for charged atoms where relevant",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Structure model modification before refinement */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel:
            "Structure model modification before refinement",
        }}
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
          <CCP4i2TaskElement
            {...taskProps}
            itemName="BFACSETUSE"
            qualifiers={{
              guiLabel: "Reset all ADPs at start",
            }}
            sx={{ width: "auto" }}
          />
          {isTruthy(bfacSetUse) && (
            <>
              <Typography variant="body1">
                to fixed value:
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="BFACSET"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </>
          )}
        </Box>

        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            gap: 1,
            flexWrap: "wrap",
          }}
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="RANDOMIZEUSE"
            qualifiers={{
              guiLabel: "Shake coordinates at start",
            }}
            sx={{ width: "auto" }}
          />
          {isTruthy(randomizeUse) && (
            <>
              <Typography variant="body1">
                with RMSD:
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="RANDOMIZE"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </>
          )}
        </Box>
      </CCP4i2ContainerElement>

      {/* Additional keywords */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Additional keywords" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="SERVALCAT_KEYWORD_FILE"
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="EXTRA_SERVALCAT_OPTIONS"
          qualifiers={{
            guiLabel: "Extra servalcat command line options:",
          }}
        />
      </CCP4i2ContainerElement>

      {/* Validation and Analysis */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Validation and Analysis" }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_IRIS"
          qualifiers={{
            guiLabel: "Generate Iris validation report",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_RAMACHANDRAN"
          qualifiers={{
            guiLabel: "Generate Ramachandran plots",
          }}
        />
        <CCP4i2TaskElement
          {...taskProps}
          itemName="VALIDATE_MOLPROBITY"
          qualifiers={{
            guiLabel: "Run MolProbity to analyse geometry",
          }}
        />

        <CCP4i2TaskElement
          {...taskProps}
          itemName="RUN_ADP_ANALYSIS"
          qualifiers={{ guiLabel: "Run ADP analysis" }}
        />
        {isTruthy(runAdpAnalysis) && (
          <>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body2">
                Atoms with a B-value lower than{" "}
                <em>the first quartile - factor * interquartile_range</em>
                {" "}or higher than{" "}
                <em>the third quartile + factor * interquartile_range</em>
                {" "}to be reported. Factor:
              </Typography>
              <Box sx={{ width: "6rem" }}>
                <CCP4i2TaskElement
                  {...taskProps}
                  itemName="ADP_IQR_FACTOR"
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
            </Box>
          </>
        )}

        <CCP4i2TaskElement
          {...taskProps}
          itemName="monitor.RUN_COORDADPDEV_ANALYSIS"
          qualifiers={{
            guiLabel:
              "Run analysis of changes in coordinates and ADPs",
          }}
        />
        {isTruthy(runCoordAdpDev) && (
          <>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="monitor.MIN_COORDDEV"
              qualifiers={{
                guiLabel:
                  "Minimum shift of atom coordinates to be reported:",
              }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="monitor.MIN_ADPDEV"
              qualifiers={{
                guiLabel:
                  "Minimum shift of B-values to be reported:",
              }}
            />
          </>
        )}
      </CCP4i2ContainerElement>

      {/* SPA-specific options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "SPA-specific options" }}
        containerHint="FolderLevel"
        visibility={() => isSpa}
      >
        <CCP4i2TaskElement
          {...taskProps}
          itemName="PIXEL_SIZE"
          qualifiers={{
            guiLabel: "Pixel size (\u00C5/pixel):",
          }}
        />
        <FieldRow>
          <CCP4i2TaskElement
            {...taskProps}
            itemName="POINTGROUP"
            qualifiers={{ guiLabel: "Point group:" }}
          />
          <CCP4i2TaskElement
            {...taskProps}
            itemName="IGNORE_SYMMETRY"
            qualifiers={{
              guiLabel: "Ignore symmetry in model file",
            }}
          />
        </FieldRow>

        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Helical parameters" }}
          containerHint="BlockLevel"
        >
          <FieldRow>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="TWIST"
              qualifiers={{
                guiLabel: "Helical twist (degrees):",
              }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="RISE"
              qualifiers={{
                guiLabel: "Helical rise (\u00C5):",
              }}
            />
          </FieldRow>
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Symmetry origin" }}
          containerHint="BlockLevel"
        >
          <FieldRow>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="CENTER_X"
              qualifiers={{ guiLabel: "X:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="CENTER_Y"
              qualifiers={{ guiLabel: "Y:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="CENTER_Z"
              qualifiers={{ guiLabel: "Z:" }}
            />
          </FieldRow>
        </CCP4i2ContainerElement>

        <CCP4i2ContainerElement
          {...taskProps}
          itemName=""
          qualifiers={{ guiLabel: "Symmetry axes" }}
          containerHint="BlockLevel"
        >
          <FieldRow>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS1_X"
              qualifiers={{ guiLabel: "Axis 1 X:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS1_Y"
              qualifiers={{ guiLabel: "Axis 1 Y:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS1_Z"
              qualifiers={{ guiLabel: "Axis 1 Z:" }}
            />
          </FieldRow>
          <FieldRow>
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS2_X"
              qualifiers={{ guiLabel: "Axis 2 X:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS2_Y"
              qualifiers={{ guiLabel: "Axis 2 Y:" }}
            />
            <CCP4i2TaskElement
              {...taskProps}
              itemName="AXIS2_Z"
              qualifiers={{ guiLabel: "Axis 2 Z:" }}
            />
          </FieldRow>
        </CCP4i2ContainerElement>

        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            gap: 1,
            flexWrap: "wrap",
          }}
        >
          <CCP4i2TaskElement
            {...taskProps}
            itemName="BLURUSE"
            qualifiers={{
              guiLabel: "Blur map with B-value:",
            }}
            sx={{ width: "auto" }}
          />
          {isTruthy(blurUse) && (
            <Box sx={{ width: "8rem" }}>
              <CCP4i2TaskElement
                {...taskProps}
                itemName="BLUR"
                qualifiers={{ guiLabel: " " }}
              />
            </Box>
          )}
        </Box>

        <CCP4i2TaskElement
          {...taskProps}
          itemName="CROSS_VALIDATION"
          qualifiers={{
            guiLabel: "Run cross validation with half maps",
          }}
        />
      </CCP4i2ContainerElement>
    </>
  );
};
