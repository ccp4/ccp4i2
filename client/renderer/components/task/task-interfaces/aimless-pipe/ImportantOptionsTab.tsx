import React from "react";
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { InlineField } from "../../task-elements/inline-field";
import { BoolToggle } from "../../task-elements/shared-hooks";

interface ImportantOptionsTabProps extends CCP4i2TaskInterfaceProps {
  sdcorrOptions: string;
  analysisOverride: BoolToggle;
  intensitiesOverride: BoolToggle;
  sdcorrOverride: BoolToggle;
  sdcorrRefine: BoolToggle;
  scalingDetails: BoolToggle;
  outputUnmerged: BoolToggle;
  vis: {
    showSimilarity: () => boolean;
  };
}

export const ImportantOptionsTab: React.FC<ImportantOptionsTabProps> = (
  props
) => {
  const {
    sdcorrOptions,
    analysisOverride,
    intensitiesOverride,
    sdcorrOverride,
    sdcorrRefine,
    scalingDetails,
    outputUnmerged,
    ...taskProps
  } = props;

  return (
    <>
      {/* Pointless symmetry options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Options for symmetry determination in Pointless",
        }}
        containerHint="FolderLevel"
      >
        <InlineField
          label="Maximum resolution for scoring set by CC(1/2) in P1 >"
          hint="[usual method, default 0.6]"
        >
          <CCP4i2TaskElement
            itemName="CCHALFLIMIT"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <InlineField
          label="Maximum resolution for scoring set by I/sigma(I) >"
          hint="[fall-back method, default 6]"
        >
          <CCP4i2TaskElement
            itemName="ISIGLIMIT"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <InlineField label="Tolerance for comparing lattices (degrees or equivalent on lengths)">
          <CCP4i2TaskElement
            itemName="TOLERANCE"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
      </CCP4i2ContainerElement>

      {/* Aimless scaling and merging options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Options for scaling and merging in Aimless",
        }}
        containerHint="FolderLevel"
      >
        {/* Analysis override */}
        <CCP4i2TaskElement
          itemName="ANALYSIS_OVERRIDE"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "override default parameters for estimation of maximum resolution",
          }}
          onChange={analysisOverride.onChange}
        />
        {analysisOverride.value && (
          <CCP4i2ContainerElement
            {...taskProps}
            itemName=""
            qualifiers={{ guiLabel: " " }}
            containerHint="BlockLevel"
          >
            <InlineField label="Set minimum information content (average bits/reflection)">
              <CCP4i2TaskElement
                itemName="INFOCONTENTTHRESHOLD"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Set minimum CC(1/2)">
              <CCP4i2TaskElement
                itemName="CCMINIMUM"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Set minimum anomalous CC(1/2)">
              <CCP4i2TaskElement
                itemName="CCANOMMINIMUM"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Set minimum MnI/sigI">
              <CCP4i2TaskElement
                itemName="ISIGMINIMUM"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>
        )}

        {/* Intensities and partials override */}
        <CCP4i2TaskElement
          itemName="INTENSITIES_OVERRIDE"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "override default parameters for selection of intensities and treatment of partials",
          }}
          onChange={intensitiesOverride.onChange}
        />
        {intensitiesOverride.value && (
          <CCP4i2ContainerElement
            {...taskProps}
            itemName=""
            qualifiers={{ guiLabel: " " }}
            containerHint="BlockLevel"
          >
            <InlineField
              label="Use"
              width="16rem"
              hint="(profile-fitted for weak intensities, summation for strong)"
            >
              <CCP4i2TaskElement
                itemName="INTENSITIES_OPTIONS"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <CCP4i2TaskElement
                itemName="PARTIALS_TEST"
                {...taskProps}
                qualifiers={{
                  guiLabel:
                    "Only accept partials with total fraction between",
                }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="PARTIALS_FRACLOW"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">and</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="PARTIALS_FRACHIGH"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
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
                itemName="PARTIALS_SCALE"
                {...taskProps}
                qualifiers={{ guiLabel: "Scale partials in range" }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="PARTIALS_SCALE_MIN"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">to lower acceptance limit</Typography>
            </Box>
            <CCP4i2TaskElement
              itemName="ACCEPT_OVERLOADS"
              {...taskProps}
              qualifiers={{ guiLabel: "accept overloaded observations" }}
            />
            <CCP4i2TaskElement
              itemName="ACCEPT_EDGES"
              {...taskProps}
              qualifiers={{
                guiLabel:
                  "accept observations on edge of tile or detector",
              }}
            />
            <CCP4i2TaskElement
              itemName="ACCEPT_XDS_MISFITS"
              {...taskProps}
              qualifiers={{
                guiLabel:
                  "accept observations flagged by XDS as outliers (MISFITS)",
              }}
            />
          </CCP4i2ContainerElement>
        )}

        {/* SD correction override */}
        <CCP4i2TaskElement
          itemName="SDCORRECTION_OVERRIDE"
          {...taskProps}
          qualifiers={{
            guiLabel: "override default parameters for SD correction",
          }}
          onChange={sdcorrOverride.onChange}
        />
        {sdcorrOverride.value && (
          <CCP4i2ContainerElement
            {...taskProps}
            itemName=""
            qualifiers={{ guiLabel: " " }}
            containerHint="BlockLevel"
          >
            <CCP4i2TaskElement
              itemName="SDCORRECTION_REFINE"
              {...taskProps}
              qualifiers={{ guiLabel: "refine SDcorrection parameters" }}
              onChange={sdcorrRefine.onChange}
            />
            {sdcorrRefine.value && (
              <>
                <InlineField
                  width="12rem"
                  hint={
                    sdcorrOptions === "SAME"
                      ? "same parameters for each run"
                      : sdcorrOptions === "SIMILAR"
                        ? "similar parameters for each run"
                        : "different parameters for each run"
                  }
                >
                  <CCP4i2TaskElement
                    itemName="SDCORRECTION_OPTIONS"
                    {...taskProps}
                    qualifiers={{ guiLabel: " " }}
                  />
                </InlineField>
                {sdcorrOptions === "SIMILAR" && (
                  <Box sx={{ pl: 3 }}>
                    <Typography
                      variant="body2"
                      sx={{ fontStyle: "italic", mb: 0.5 }}
                    >
                      Similarity restraint SDs:
                    </Typography>
                    <FieldRow equalWidth={false} size="xs">
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SIMILARITY_SDFAC"
                        {...taskProps}
                        qualifiers={{ guiLabel: "sdFac" }}
                      />
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SIMILARITY_SDB"
                        {...taskProps}
                        qualifiers={{ guiLabel: "sdB" }}
                      />
                      <CCP4i2TaskElement
                        itemName="SDCORRECTION_SIMILARITY_SDADD"
                        {...taskProps}
                        qualifiers={{ guiLabel: "sdAdd" }}
                      />
                    </FieldRow>
                  </Box>
                )}
                <Box
                  sx={{
                    display: "flex",
                    alignItems: "center",
                    gap: 1,
                    flexWrap: "wrap",
                  }}
                >
                  <CCP4i2TaskElement
                    itemName="SDCORRECTION_FIXSDB"
                    {...taskProps}
                    qualifiers={{ guiLabel: "fix sdB parameters" }}
                  />
                  <Typography variant="body1">
                    SD to tie sdB to zero
                  </Typography>
                  <Box sx={{ width: "8rem" }}>
                    <CCP4i2TaskElement
                      itemName="SDCORRECTION_TIESDB_SD"
                      {...taskProps}
                      qualifiers={{ guiLabel: " " }}
                    />
                  </Box>
                </Box>
              </>
            )}
            {!sdcorrRefine.value && (
              <>
                <CCP4i2TaskElement
                  itemName="SDCORRECTION_SET"
                  {...taskProps}
                  qualifiers={{
                    guiLabel: "Set SD correction parameters",
                  }}
                />
                <FieldRow equalWidth={false} size="xs">
                  <CCP4i2TaskElement
                    itemName="SDCORRECTION_SDFAC"
                    {...taskProps}
                    qualifiers={{ guiLabel: "sdFac" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SDCORRECTION_SDB"
                    {...taskProps}
                    qualifiers={{ guiLabel: "sdB" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SDCORRECTION_SDADD"
                    {...taskProps}
                    qualifiers={{ guiLabel: "sdAdd" }}
                  />
                </FieldRow>
              </>
            )}
          </CCP4i2ContainerElement>
        )}

        {/* Scaling details override */}
        <CCP4i2TaskElement
          itemName="SCALING_DETAILS"
          {...taskProps}
          qualifiers={{
            guiLabel: "override default parameters for scaling details",
          }}
          onChange={scalingDetails.onChange}
        />
        {scalingDetails.value && (
          <CCP4i2ContainerElement
            {...taskProps}
            itemName=""
            qualifiers={{ guiLabel: " " }}
            containerHint="BlockLevel"
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
                itemName="CYCLES_FLAG"
                {...taskProps}
                qualifiers={{ guiLabel: "Refine scale factors for" }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="CYCLES_N"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">cycles</Typography>
            </Box>

            <InlineField
              label="use multiple processors"
              hint="determined from number of reflections"
              width="auto"
            >
              <CCP4i2TaskElement
                itemName="PARALLEL"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>

            {/* Selection thresholds */}
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <CCP4i2TaskElement
                itemName="SELECT1"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="SELECT_IOVSDMIN"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">
                minimum I/sd for 1st round scaling
              </Typography>
              <CCP4i2TaskElement
                itemName="SELECT2"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="SELECT_EMIN"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">
                minimum E for 2nd round scaling
              </Typography>
            </Box>

            {/* Tie restraints: [label] [checkbox] [SD input] */}
            <InlineField
              label="Restrain neighbouring scale factors on rotation axis with SD"
              after={
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="TIE_ROTATION_SD"
                    {...taskProps}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              }
              width="auto"
            >
              <CCP4i2TaskElement
                itemName="TIE_ROTATION"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField
              label="Restrain surface parameters to a sphere with SD"
              after={
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="TIE_SURFACE_SD"
                    {...taskProps}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              }
              width="auto"
            >
              <CCP4i2TaskElement
                itemName="TIE_SURFACE"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField
              label="Restrain neighbouring B-factors on rotation axis with SD"
              after={
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="TIE_BFACTOR_SD"
                    {...taskProps}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              }
              width="auto"
            >
              <CCP4i2TaskElement
                itemName="TIE_BFACTOR"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField
              label="Restrain B-factors to zero with SD"
              after={
                <Box sx={{ width: "8rem" }}>
                  <CCP4i2TaskElement
                    itemName="TIE_BZERO_SD"
                    {...taskProps}
                    qualifiers={{ guiLabel: " " }}
                  />
                </Box>
              }
              width="auto"
            >
              <CCP4i2TaskElement
                itemName="TIE_BZERO"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>
        )}

        {/* Output options */}
        <Box
          sx={{
            display: "flex",
            alignItems: "center",
            gap: 2,
            flexWrap: "wrap",
            mt: 1,
          }}
        >
          <CCP4i2TaskElement
            itemName="OUTPUT_UNMERGED"
            {...taskProps}
            qualifiers={{
              guiLabel: "output unmerged data as well as merged",
            }}
            onChange={outputUnmerged.onChange}
          />
          {outputUnmerged.value && (
            <CCP4i2TaskElement
              itemName="ORIGINAL_HKL"
              {...taskProps}
              qualifiers={{
                guiLabel: "output original hkl instead of reduced hkl",
              }}
            />
          )}
        </Box>
      </CCP4i2ContainerElement>

      {/* FreeR extension options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Options for FreeR set extension" }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" sx={{ fontStyle: "italic" }}>
          If you are extending an existing FreeR set, it must match the
          observed data in unit cell and Laue group
        </Typography>
        <Typography variant="body2" sx={{ fontStyle: "italic" }}>
          The unit cells should match to the lower resolution of the two
          datasets
        </Typography>
        <CCP4i2TaskElement
          itemName="OVERRIDE_CELL_DIFFERENCE"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "allow existing freeR set to have different unit cells",
          }}
        />
        <Typography
          variant="body2"
          sx={{ fontStyle: "italic", color: "error.main" }}
        >
          Be sure you know what you are doing: the cells must be very
          similar even if outside the test limits
        </Typography>
      </CCP4i2ContainerElement>
    </>
  );
};
