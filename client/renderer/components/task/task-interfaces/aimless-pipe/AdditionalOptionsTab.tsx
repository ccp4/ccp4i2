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
import { Box, Typography } from "@mui/material";
import { CCP4i2TaskInterfaceProps } from "../task-container";
import { CCP4i2TaskElement } from "../../task-elements/task-element";
import { CCP4i2ContainerElement } from "../../task-elements/ccontainer";
import { FieldRow } from "../../task-elements/field-row";
import { InlineField } from "../../task-elements/inline-field";
import { BoolToggle } from "../../task-elements/shared-hooks";

interface AdditionalOptionsTabProps extends CCP4i2TaskInterfaceProps {
  scalingProtocol: string;
  scalesRotType: string;
  scalesBrotType: string;
  removeLattice: BoolToggle;
  keepLattice: BoolToggle;
  outlierOverride: BoolToggle;
  expertOptions: BoolToggle;
  vis: {
    hasRotationScales: () => boolean;
    hasSecondaryBeam: () => boolean;
    showTileXY: () => boolean;
    rotSpacing: () => boolean;
    rotNbins: () => boolean;
    brotSpacing: () => boolean;
    brotNbins: () => boolean;
    showLatticeThreshold: () => boolean;
    showBatchList: () => boolean;
  };
}

export const AdditionalOptionsTab: React.FC<AdditionalOptionsTabProps> = (
  props
) => {
  const {
    scalingProtocol,
    scalesRotType,
    scalesBrotType,
    removeLattice,
    keepLattice,
    outlierOverride,
    expertOptions,
    vis,
    ...taskProps
  } = props;

  return (
    <>
      {/* Options for symmetry determination in Pointless */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Options for symmetry determination in Pointless",
        }}
        containerHint="FolderLevel"
      >
        <Typography
          variant="body2"
          sx={{ fontStyle: "italic", color: "primary.main", mb: 0.5 }}
        >
          Choice of cell setting conventions:
        </Typography>
        <CCP4i2TaskElement
          itemName="SET_SETTING"
          {...taskProps}
          qualifiers={{ guiLabel: " " }}
        />
      </CCP4i2ContainerElement>

      {/* Option to remove centred lattice absences in Pointless */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel:
            "Option to remove centred lattice absences in Pointless",
        }}
        containerHint="FolderLevel"
      >
        <Typography variant="body2" sx={{ fontStyle: "italic" }}>
          If the input data belong to a primitive lattice (P), the data are
          checked for additional centred lattice symmetry
        </Typography>
        <Typography variant="body2" sx={{ fontStyle: "italic" }}>
          By default, extra centred lattice reflections will be removed if
          the estimated &quot;probability&quot; of a centred lattice
        </Typography>
        <Typography variant="body2" sx={{ fontStyle: "italic", mb: 1 }}>
          is greater than the THRESHOLD set here or by default. KEEP and
          REMOVE options are unconditional
        </Typography>
        <InlineField label="Probability threshold for removing centred lattice reflections">
          <CCP4i2TaskElement
            itemName="LATTICE_CENTERING_THRESHOLD"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>
        <CCP4i2TaskElement
          itemName="KEEP_LATTICE_CENTERING"
          {...taskProps}
          qualifiers={{
            guiLabel: "Always keep centred lattice reflections (KEEP)",
          }}
          onChange={keepLattice.onChange}
        />
        <CCP4i2TaskElement
          itemName="REMOVE_LATTICE_CENTERING"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "Unconditionally remove centred lattice absences (REMOVE)",
          }}
          onChange={removeLattice.onChange}
        />
        {removeLattice.value && (
          <Box sx={{ pl: 3 }}>
            <CCP4i2TaskElement
              itemName="LATTICE_CENTERING"
              {...taskProps}
              qualifiers={{ guiLabel: "Lattice centering type" }}
            />
          </Box>
        )}
      </CCP4i2ContainerElement>

      {/* Options for scaling in Aimless */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{ guiLabel: "Options for scaling in Aimless" }}
        containerHint="FolderLevel"
      >
        {/* Scale protocol with relative B-factor */}
        <InlineField
          label="Scale"
          width="20rem"
          after={
            scalingProtocol !== "ONLYMERGE" ? (
              <CCP4i2TaskElement
                itemName="BFACTOR_SCALE"
                {...taskProps}
                qualifiers={{ guiLabel: "with relative B-factor" }}
              />
            ) : undefined
          }
        >
          <CCP4i2TaskElement
            itemName="SCALING_PROTOCOL"
            {...taskProps}
            qualifiers={{ guiLabel: " " }}
          />
        </InlineField>

        {/* Rotation scales (ROTATION or SECONDARY protocols) */}
        {vis.hasRotationScales() && (
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
              <Typography variant="body1">
                Define scale ranges along rotation axis by
              </Typography>
              <Box sx={{ width: "12rem" }}>
                <CCP4i2TaskElement
                  itemName="SCALES_ROTATION_TYPE"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <CCP4i2TaskElement
                itemName="SCALES_ROTATION_SPACING"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
                visibility={vis.rotSpacing}
              />
              <CCP4i2TaskElement
                itemName="SCALES_ROTATION_NBINS"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
                visibility={vis.rotNbins}
              />
              {scalesRotType === "SPACING" && (
                <Typography variant="body1">degrees</Typography>
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
              <Typography variant="body1">
                Define B-factor ranges along rotation axis by
              </Typography>
              <Box sx={{ width: "12rem" }}>
                <CCP4i2TaskElement
                  itemName="SCALES_BROTATION_TYPE"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <CCP4i2TaskElement
                itemName="SCALES_BROTATION_SPACING"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
                visibility={vis.brotSpacing}
              />
              <CCP4i2TaskElement
                itemName="SCALES_BROTATION_NBINS"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
                visibility={vis.brotNbins}
              />
              {scalesBrotType === "SPACING" && (
                <Typography variant="body1">degrees</Typography>
              )}
            </Box>
          </CCP4i2ContainerElement>
        )}

        {/* Secondary beam correction (SECONDARY protocol only) */}
        {vis.hasSecondaryBeam() && (
          <>
            <InlineField label="Maximum order of spherical harmonics for secondary beam correction (eg 4 or 6)">
              <CCP4i2TaskElement
                itemName="SCALES_SECONDARY_NSPHHARMONICS"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <CCP4i2ContainerElement
              {...taskProps}
              itemName=""
              qualifiers={{ guiLabel: " " }}
              containerHint="BlockLevel"
            >
              <InlineField label="Tile scaling for CCD detectors" width="12rem">
                <CCP4i2TaskElement
                  itemName="SCALES_TILETYPE"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </InlineField>
              {vis.showTileXY() && (
                <FieldRow equalWidth={false} size="xs">
                  <CCP4i2TaskElement
                    itemName="SCALES_NTILEX"
                    {...taskProps}
                    qualifiers={{ guiLabel: "Tiles on X" }}
                  />
                  <CCP4i2TaskElement
                    itemName="SCALES_NTILEY"
                    {...taskProps}
                    qualifiers={{ guiLabel: "Tiles on Y" }}
                  />
                </FieldRow>
              )}
            </CCP4i2ContainerElement>
          </>
        )}

        {/* Outlier rejection */}
        <CCP4i2TaskElement
          itemName="OUTLIER_OVERRIDE"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "override default parameters for outlier rejection",
          }}
          onChange={outlierOverride.onChange}
        />
        {outlierOverride.value && (
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
              <Typography variant="body1">
                Reject outliers if &gt;
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="OUTLIER_SDMAX"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">from mean, or &gt;</Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="OUTLIER_SDMAX2"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">if 2 observations</Typography>
            </Box>
            <Box
              sx={{
                display: "flex",
                alignItems: "center",
                gap: 1,
                flexWrap: "wrap",
              }}
            >
              <Typography variant="body1">
                Reject outliers between I+ and I- if &gt;
              </Typography>
              <Box sx={{ width: "8rem" }}>
                <CCP4i2TaskElement
                  itemName="OUTLIER_SDMAXALL"
                  {...taskProps}
                  qualifiers={{ guiLabel: " " }}
                />
              </Box>
              <Typography variant="body1">,</Typography>
              <CCP4i2TaskElement
                itemName="OUTLIER_SDMAXALL_ADJUST"
                {...taskProps}
                qualifiers={{
                  guiLabel:
                    "increase for large anomalous differences",
                }}
              />
            </Box>
            <CCP4i2TaskElement
              itemName="OUTLIER_COMBINE"
              {...taskProps}
              qualifiers={{
                guiLabel: "compare outliers within each datasets",
              }}
            />
            <InlineField label="Set maximum E to reject unreasonably large intensities">
              <CCP4i2TaskElement
                itemName="OUTLIER_EMAX"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </CCP4i2ContainerElement>
        )}
      </CCP4i2ContainerElement>

      {/* Options for both Pointless and Aimless */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Options for both Pointless and Aimless",
        }}
        containerHint="FolderLevel"
      >
        <Typography
          variant="body2"
          sx={{ fontStyle: "italic", color: "primary.main", mb: 0.5 }}
        >
          Override automatic definition of runs to mark discontinuities in
          data
        </Typography>
        <CCP4i2TaskElement
          itemName="RUN_MODE"
          {...taskProps}
          qualifiers={{
            guiLabel: "Run selection options",
            guiMode: "radio",
          }}
        />
        <CCP4i2TaskElement
          itemName="RUN_BATCHLIST"
          {...taskProps}
          qualifiers={{
            guiLabel:
              "Batch ranges to define runs (after any renumbering in Pointless), with optional high resolution limit",
          }}
          visibility={vis.showBatchList}
        />
      </CCP4i2ContainerElement>

      {/* Expert options */}
      <CCP4i2ContainerElement
        {...taskProps}
        itemName=""
        qualifiers={{
          guiLabel: "Expert options, not for normal use",
        }}
        containerHint="FolderLevel"
      >
        <CCP4i2TaskElement
          itemName="EXPERT_OPTIONS"
          {...taskProps}
          qualifiers={{ guiLabel: "show expert options" }}
          onChange={expertOptions.onChange}
        />
        {expertOptions.value && (
          <Box sx={{ pl: 3 }}>
            <CCP4i2TaskElement
              itemName="ALLOW_NONCHIRAL"
              {...taskProps}
              qualifiers={{
                guiLabel: "Allow non-chiral space groups",
              }}
            />
            <InlineField label="CC(1/2) inner shell disaster limit">
              <CCP4i2TaskElement
                itemName="CCHALFDISASTER"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
            <InlineField label="Minimum multiplicity, inner disaster limit">
              <CCP4i2TaskElement
                itemName="MULTIPLICITYDISASTER"
                {...taskProps}
                qualifiers={{ guiLabel: " " }}
              />
            </InlineField>
          </Box>
        )}
      </CCP4i2ContainerElement>
    </>
  );
};
