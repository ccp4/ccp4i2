import { CCP4i2TaskInterfaceProps } from "./task-container";
import {
  CCP4i2TaskElement,
  CCP4i2TaskElementProps,
} from "../task-elements/task-element";
import { CCP4i2Tab, CCP4i2Tabs } from "../task-elements/tabs";
import { doRetrieve, useApi } from "../../../api";
import { useJob, usePrevious } from "../../../utils";
import { CCP4i2ContainerElement } from "../task-elements/ccontainer";
import { useCallback, useEffect, useMemo } from "react";
import { showMtzColumnDialog, parseMtzColumns } from "../task-elements/mtz-column-dialog";
import { Job } from "../../../types/models";
import { useCCP4i2Window } from "../../../app-context";
import {
  Card,
  CardContent,
  CardHeader,
  Grid2,
  Paper,
  Stack,
  Typography,
} from "@mui/material";

const TaskInterface: React.FC<CCP4i2TaskInterfaceProps> = (props) => {
  const api = useApi();
  const { job } = props;
  const { useFileDigest, useTaskItem, mutateContainer, mutateValidation, uploadFileParam } =
    useJob(job.id);
  const { cootModule } = useCCP4i2Window();

  const { item: HKLINItem, value: HKLINValue } = useTaskItem("HKLIN");

  const oldHKLINValue = usePrevious(HKLINValue);

  const { data: HKLINDigest } = useFileDigest(HKLINItem?._objectPath);

  const { item: HKLIN_OBSItem } = useTaskItem("HKLIN_OBS");

  const { update: updateSPACEGROUP } = useTaskItem("SPACEGROUP");

  const { update: updateUNITCELL } = useTaskItem("UNITCELL");

  const { update: updateWAVELENGTH } = useTaskItem("WAVELENGTH");

  const { update: updateCRYSTALNAME } = useTaskItem("CRYSTALNAME");

  const { update: updateDATASETNAME } = useTaskItem("DATASETNAME");

  const { update: setHKLIN_OBS_COLUMNS } = useTaskItem("HKLIN_OBS_COLUMNS");

  const { update: setHKLIN_OBS_CONTENT_FLAG } = useTaskItem(
    "HKLIN_OBS_CONTENT_FLAG"
  );

  const { value: HKLIN_FORMATValue, update: updateHKLIN_FORMAT } =
    useTaskItem("HKLIN_FORMAT");

  // Define a useCallback which is what will be called when the file digest changes.
  //This is wrapped in a useCallback, since it will include calls to methods that will be defined dynamically
  // and we do not want the actual function to be called when those methods are updated (hence not symply put into the useEffect)
  const handleDigestChanged = useCallback(
    (digest: any) => {
      // API returns {success: true, data: {...}} - extract the data
      const digestData = digest?.data;
      if (
        !digestData ||
        Object.keys(digestData).length == 0 ||
        !updateSPACEGROUP ||
        !updateWAVELENGTH ||
        !updateHKLIN_FORMAT ||
        !updateUNITCELL ||
        job?.status != 1
      )
        return;

      const asyncFunc = async () => {
        let parametersChanged = false;
        {
          const result = await updateSPACEGROUP(
            digestData.spaceGroup.replace(/\s+/g, "")
          );
          parametersChanged = parametersChanged || Boolean(result);
        }
        if (digestData.wavelength) {
          //Note some ancient MTZ files lack a wavelength
          const result = await updateWAVELENGTH(digestData.wavelength);
          parametersChanged = parametersChanged || Boolean(result);
        }
        {
          const result = await updateHKLIN_FORMAT(digestData.format.toUpperCase());
          parametersChanged = parametersChanged || Boolean(result);
        }
        {
          const result = await updateUNITCELL(digestData.cell);
          parametersChanged = parametersChanged || Boolean(result);
        }
        if (parametersChanged) {
          await mutateContainer();
          await mutateValidation();
        }
      };
      return asyncFunc();
    },
    [updateSPACEGROUP, updateWAVELENGTH, updateHKLIN_FORMAT, updateUNITCELL]
  );

  //And here the simple useEffect
  useEffect(() => {
    if (!handleDigestChanged) return;
    const asyncFunc = async () => {
      await handleDigestChanged(HKLINDigest);
    };
    asyncFunc();
  }, [HKLINDigest]);

  // Process the selected columns after user picks from dialog
  const processColumnSelection = useCallback(
    async (columnPath: string, file: File) => {
      if (!setHKLIN_OBS_CONTENT_FLAG) return;
      if (!setHKLIN_OBS_COLUMNS) return;

      // API returns {success: true, data: {...}} - extract the data
      const digestData = HKLINDigest?.data;

      const match = columnPath.match(/\[([^\]]+)\]/);
      console.log({ match });
      if (match) {
        await setHKLIN_OBS_COLUMNS(match[1]);
        const columnNames = match[1].split(",").map((name) => name.trim());
        const columnTypes = columnNames.map(
          (name) =>
            digestData?.listOfColumns?.find(
              (col: { columnLabel: string }) => col.columnLabel === name
            )?.columnType
        );
        const datasetNames = columnNames.map(
          (name) =>
            digestData?.listOfColumns?.find(
              (col: { columnLabel: string }) => col.columnLabel === name
            )?.dataset
        );

        const signature = columnTypes.join("");

        let datasetIndex: number = -1;

        // Check if datasetNames conditions are met
        if (
          datasetNames.length === columnNames.length &&
          datasetNames.length > 0 &&
          datasetNames.every((name) => name === datasetNames[0])
        ) {
          await updateDATASETNAME(datasetNames[0]);
          const consensusDatasetName = datasetNames[0];

          // Try to find matching index in datasets array
          datasetIndex = digestData?.datasets?.indexOf(consensusDatasetName) ?? -1;
        }
        // Handle crystal names if available
        let selectedCrystalName = "Xtal1"; // Default fallback
        if (
          datasetIndex !== -1 &&
          Array.isArray(digestData?.crystalNames) &&
          datasetIndex < digestData.crystalNames.length
        ) {
          // Option a) Use crystal name at same index as consensus dataset
          selectedCrystalName = digestData.crystalNames[datasetIndex];
        } else if (Array.isArray(digestData?.crystalNames)) {
          // Option b) Find last entry that is not "HKL_base"
          const validCrystalNames = digestData.crystalNames.filter(
            (name: string) => name !== "HKL_base"
          );

          if (validCrystalNames.length > 0) {
            selectedCrystalName =
              validCrystalNames[validCrystalNames.length - 1];
          }
          // If no valid names found, keep default "Xtal1"
        }
        updateCRYSTALNAME(selectedCrystalName);

        // Handle wavelengths if available
        let selectedWavelength = "1.0"; // Default fallback
        if (
          datasetIndex !== -1 &&
          Array.isArray(digestData?.wavelengths) &&
          datasetIndex < digestData.wavelengths.length
        ) {
          // Option a) Use wavelength at same index as consensus dataset
          selectedWavelength = digestData.wavelengths[datasetIndex];
        } else if (Array.isArray(digestData?.wavelengths)) {
          // Option b) Find last entry that is not "HKL_base"
          const validWavelengths = digestData.wavelengths.filter(
            (name: string) => name !== "HKL_base"
          );

          if (validWavelengths.length > 0) {
            selectedWavelength = validWavelengths[validWavelengths.length - 1];
          } else if (digestData.wavelengths.length > 0) {
            selectedWavelength = digestData.wavelengths[0]; // Fallback to first wavelength
          }
          // If no valid wavelengths found, keep default "1.0"
        }
        updateWAVELENGTH(selectedWavelength);

        // Determine content flag based on signature
        const contentFlag = ["KMKM", "GLGL", "JQ", "FQ"].indexOf(signature);
        console.log({ signature, contentFlag, columnPath, columnNames });
        if (contentFlag > -1) {
          await setHKLIN_OBS_CONTENT_FLAG(contentFlag + 1);
        }
      }

      // Upload the file with selected columns using centralized uploadFileParam (with intent tracking)
      if (columnPath && columnPath.trim().length > 0) {
        await uploadFileParam({
          objectPath: HKLIN_OBSItem._objectPath,
          file: file,
          fileName: file.name,
          columnSelector: columnPath,
        });
        // Note: uploadFileParam already handles mutateContainer and mutateValidation
      }
    },
    [
      job,
      HKLIN_OBSItem,
      setHKLIN_OBS_COLUMNS,
      setHKLIN_OBS_CONTENT_FLAG,
      HKLINDigest,
      updateDATASETNAME,
      updateCRYSTALNAME,
      updateWAVELENGTH,
      uploadFileParam,
    ]
  );

  //Here handle a case where a new MTZ or mmcif file is uploaded as HKLIN.
  //Uses Promise-based dialog instead of component state
  const handleHKLINFileChange = useCallback(
    async (hklinValue: any) => {
      if (
        !hklinValue?.dbFileId ||
        !hklinValue?.baseName ||
        !oldHKLINValue ||
        !cootModule ||
        job?.status != 1
      )
        return;
      if (JSON.stringify(hklinValue) === JSON.stringify(oldHKLINValue)) return;

      // Download the file
      const downloadURL = `files/${hklinValue.dbFileId}/download_by_uuid`;
      const arrayBuffer = await doRetrieve(downloadURL, hklinValue.baseName);
      const blob = new Blob([arrayBuffer], {
        type: "application/CCP4-mtz-file",
      });
      const file = new File([blob], hklinValue.baseName, {
        type: "application/CCP4-mtz-file",
      });

      // Check if it's an MTZ file that needs column selection
      const isMtzFile = file.name.toLowerCase().endsWith(".mtz");
      if (!isMtzFile) {
        return;
      }

      // Parse the MTZ file columns
      const columnNames = await parseMtzColumns(file, cootModule);
      if (!columnNames) {
        return;
      }

      // Show the column selection dialog
      const columnPath = await showMtzColumnDialog(columnNames, HKLIN_OBSItem);
      if (!columnPath) {
        return; // User cancelled
      }

      // Process the selected columns
      await processColumnSelection(columnPath, file);
    },
    [oldHKLINValue, cootModule, job?.status, HKLIN_OBSItem, processColumnSelection]
  );

  //And then a simple useEffect
  useEffect(() => {
    const asyncFunc = async () => {
      if (HKLINValue) {
        await handleHKLINFileChange(HKLINValue);
      }
    };
    asyncFunc();
  }, [HKLINValue]);

  return (
    <>
      <CCP4i2Tabs {...props}>
        <CCP4i2Tab label="Main inputs" key="1">
          <CCP4i2ContainerElement
            {...props}
            itemName=""
            qualifiers={{ guiLabel: "Input data", initiallyOpen: true }}
            key="Input data"
            containerHint="FolderLevel"
          >
            <CCP4i2TaskElement
              {...props}
              key="HKLIN"
              itemName="HKLIN"
              qualifiers={{ guiLabel: "Reflections" }}
            />
            <Card>
              <CardHeader variant="primary" title="Additional information" />
              <CardContent>
                <Grid2 container direction="row" sx={{ mb: 2 }} spacing={2}>
                  <Grid2 size={{ xs: 6, md: 4 }}>
                    <CCP4i2TaskElement
                      {...props}
                      key="SPACEGROUP"
                      itemName="SPACEGROUP"
                      qualifiers={{ guiLabel: "Space group" }}
                    />
                  </Grid2>
                  <Grid2 size={{ xs: 6, md: 8 }}>
                    <CCP4i2TaskElement
                      {...props}
                      key="UNITCELL"
                      itemName="UNITCELL"
                    />
                  </Grid2>
                  <Stack direction="row" spacing={2}>
                    <CCP4i2TaskElement
                      {...props}
                      key="CRYSTALNAME"
                      itemName="CRYSTALNAME"
                      qualifiers={{ guiLabel: "Crystal name" }}
                      sx={{ maxWidth: "3rem" }}
                    />
                    <CCP4i2TaskElement
                      {...props}
                      key="DATASETNAME"
                      itemName="DATASETNAME"
                      qualifiers={{ guiLabel: "Dataset name" }}
                      sx={{ maxWidth: "3rem" }}
                    />
                  </Stack>
                </Grid2>
                <CCP4i2TaskElement
                  {...props}
                  key="WAVELENGTH"
                  itemName="WAVELENGTH"
                  qualifiers={{ guiLabel: "Wavelength" }}
                />
              </CardContent>
            </Card>

            {false && (
              <CCP4i2TaskElement
                {...props}
                key="HKLIN_FORMAT"
                itemName="HKLIN_FORMAT"
                qualifiers={{ guiLabel: "HKLIN format" }}
              />
            )}
            {HKLIN_FORMATValue === "MTZ" && (
              <MtzPanel {...props} itemName="" digest={HKLINDigest} />
            )}
            {HKLIN_FORMATValue === "MMCIF" && (
              <MmcifPanel {...props} itemName="" digest={HKLINDigest} />
            )}
            <CCP4i2TaskElement {...props} itemName="FREERFLAG" />
          </CCP4i2ContainerElement>
        </CCP4i2Tab>
      </CCP4i2Tabs>
    </>
  );
};

interface MmcifPanelProps extends CCP4i2TaskElementProps {
  digest: any;
}

const MmcifPanel: React.FC<MmcifPanelProps> = (props) => {
  const api = useApi();
  const { digest, job } = props;
  // API returns {success: true, data: {...}} - extract the data
  const digestData = digest?.data;
  const { useTaskItem, mutateContainer, mutateValidation } = useJob(job.id);
  const { value: MMCIF_SELECTED_BLOCKValue } = useTaskItem(
    "MMCIF_SELECTED_BLOCK"
  );
  const oldMMCIF_SELECTED_BLOCKValue = usePrevious(MMCIF_SELECTED_BLOCKValue);

  const { update: setMMCIF_SELECTED_ISINTENSITY } = useTaskItem(
    "MMCIF_SELECTED_ISINTENSITY"
  );

  const { update: setMMCIF_SELECTED_COLUMNS } = useTaskItem(
    "MMCIF_SELECTED_COLUMNS"
  );

  const { value: MMCIF_SELECTED_INFOValue, update: setMMCIF_SELECTED_INFO } =
    useTaskItem("MMCIF_SELECTED_INFO");

  const handleSelectedBlockChange = useCallback(
    async (mmcifSelectedBlockName) => {
      if (!digestData?.rblock_infos || job?.status != 1) return;
      if (
        !setMMCIF_SELECTED_COLUMNS ||
        !setMMCIF_SELECTED_ISINTENSITY ||
        !setMMCIF_SELECTED_INFO
      )
        return;
      if (mmcifSelectedBlockName === oldMMCIF_SELECTED_BLOCKValue) return;
      if (!digestData?.format || digestData?.format !== "mmcif") return;
      const asyncFunc = async () => {
        if (digestData?.rblock_infos) {
          const selectedBlock = digestData.rblock_infos.find(
            (info: { bname: string }) => info.bname === mmcifSelectedBlockName
          );
          if (selectedBlock) {
            const newColumns =
              selectedBlock.labelsets?.columnsets[0]?.columnlabels?.join(", ");
            await setMMCIF_SELECTED_COLUMNS(newColumns);
            const isIntensity = selectedBlock.info?.includes("I") ? 1 : -1;
            await setMMCIF_SELECTED_ISINTENSITY(isIntensity);
            const blockInfo = `${selectedBlock.info}\nhkl list type: ${selectedBlock.hklcheckformat}\nHighest resolution: ${selectedBlock.highres}`;
            await setMMCIF_SELECTED_INFO(blockInfo);
          }
        } else {
          await setMMCIF_SELECTED_COLUMNS("");
          await setMMCIF_SELECTED_ISINTENSITY("");
          await setMMCIF_SELECTED_INFO("");
        }
        mutateContainer();
        mutateValidation();
      };
      asyncFunc();
    },
    [
      oldMMCIF_SELECTED_BLOCKValue,
      digestData,
      setMMCIF_SELECTED_COLUMNS,
      setMMCIF_SELECTED_ISINTENSITY,
      setMMCIF_SELECTED_INFO,
    ]
  );

  useEffect(() => {
    if (!handleSelectedBlockChange) return;
    const asyncFunc = async () => {
      await handleSelectedBlockChange(MMCIF_SELECTED_BLOCKValue);
    };
    asyncFunc();
  }, [MMCIF_SELECTED_BLOCKValue]);

  return (
    digestData?.rblock_infos && (
      <>
        <Card>
          <CardHeader variant="primary" title="Mmcif file information" />
          <CardContent>
            <CCP4i2TaskElement
              {...props}
              key="MMCIF_SELECTED_BLOCK"
              itemName="MMCIF_SELECTED_BLOCK"
              qualifiers={{
                guiLabel: "Selected block",
                guiMode: "multiLineRadio",
                onlyEnumerators: true,
                enumerators: digestData.rblock_infos.map(
                  (info: { bname: string }) => info.bname
                ),
              }}
            />
            <pre>{MMCIF_SELECTED_INFOValue}</pre>
            <Stack direction={"row"} spacing={2} sx={{ mb: 2 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ width: "100%" }}
                key="MMCIF_SELECTED_COLUMNS"
                itemName="MMCIF_SELECTED_COLUMNS"
                qualifiers={{
                  guiLabel: "Selected columns",
                }}
                disabled={true}
              />
              <CCP4i2TaskElement
                {...props}
                sx={{ width: "100%" }}
                key="MMCIF_SELECTED_ISINTENSITY"
                itemName="MMCIF_SELECTED_ISINTENSITY"
                qualifiers={{
                  guiLabel: "Selected is intensity",
                }}
                disabled={true}
              />
            </Stack>
          </CardContent>
        </Card>
      </>
    )
  );
};

const MtzPanel: React.FC<MmcifPanelProps> = (props) => {
  const { digest, job } = props;
  const { useTaskItem } = useJob(job.id);
  return (
    <Paper sx={{ border: "1px solid black", px: 2, py: 1, mb: 2 }}>
      <Typography variant="h6" gutterBottom>
        MTZ file information
      </Typography>
      <CCP4i2TaskElement
        {...props}
        key="HKLIN_OBS"
        itemName="HKLIN_OBS"
        qualifiers={{ guiLabel: "Parsed MTZ file" }}
        disabled={() => true}
      />
      <Stack direction="row">
        <CCP4i2TaskElement
          {...props}
          key="HKLIN_OBS_CONTENT_FLAG"
          itemName="HKLIN_OBS_CONTENT_FLAG"
          qualifiers={{ guiLabel: "Content flag" }}
          disabled={() => false}
        />
        <CCP4i2TaskElement
          {...props}
          key="HKLIN_OBS_COLUMNS"
          itemName="HKLIN_OBS_COLUMNS"
          qualifiers={{ guiLabel: "Columns" }}
          disabled={() => false}
        />
      </Stack>
    </Paper>
  );
};

export default TaskInterface;
