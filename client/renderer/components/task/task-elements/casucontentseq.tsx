import { CCP4i2TaskElement, CCP4i2TaskElementProps } from "./task-element";
import { CCP4i2ContainerElement } from "./ccontainer";
import { Card, CardContent, CardHeader, Grid2 } from "@mui/material";
import { useApi } from "../../../api";
import { useJob, usePrevious, valueOfItem } from "../../../utils";
import { ErrorInfo } from "./error-info";
import { apiGet } from "../../../api-fetch";
import { useCallback, useEffect } from "react";
import { useParameterChangeIntent } from "../../../providers/parameter-change-intent-provider";

export const CAsuContentSeqElement: React.FC<CCP4i2TaskElementProps> = (
  props
) => {
  const api = useApi();
  const { itemName, job } = props;
  const { useTaskItem, getValidationColor, mutateContainer } = useJob(job.id);
  const { clearIntentForPath } = useParameterChangeIntent();

  const { item } = useTaskItem(itemName);
  const { update: setPolymerType } = useTaskItem(
    `${item._objectPath}.polymerType`
  );
  const { update: setName } = useTaskItem(`${item._objectPath}.name`);
  const { update: setSequence } = useTaskItem(`${item._objectPath}.sequence`);
  const { update: setDescription } = useTaskItem(
    `${item._objectPath}.description`
  );
  const { update: setNCopies } = useTaskItem(`${item._objectPath}.nCopies`);
  const setSEQUENCEFromSEQIN = useCallback(
    async (seqinDigestResponse: any, annotation: string) => {
      // API returns {success: true, data: {...}} - extract the data
      const seqinDigest = seqinDigestResponse?.data;
      console.log("Setting SEQUENCE from SEQIN digest", seqinDigest);
      if (
        !setSequence ||
        !setName ||
        !setPolymerType ||
        !setDescription ||
        !setNCopies ||
        !item ||
        !seqinDigest ||
        job?.status != 1
      ) {
        console.log("Cannot set SEQUENCE from SEQIN - missing data", seqinDigest);
        return;
      }
      if (seqinDigest?.moleculeType) {
        console.log("Seqin digest was a sequence file");
        const { name, moleculeType, sequence } = seqinDigest || {};
        const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
        // Each setParameter call triggers mutations, so container/validation
        // will be updated automatically after each call
        await setPolymerType(moleculeType);
        await setName(sanitizedName);
        await setSequence(sequence);
        await setDescription(annotation);
        await setNCopies(1);
        // Clear intents so child components will sync their local state
        // with the new container data (otherwise wasRecentlyChanged blocks the sync)
        clearIntentForPath(`${item._objectPath}.polymerType`);
        clearIntentForPath(`${item._objectPath}.name`);
        clearIntentForPath(`${item._objectPath}.sequence`);
        clearIntentForPath(`${item._objectPath}.description`);
        clearIntentForPath(`${item._objectPath}.nCopies`);
        // Force a final container refresh to ensure modal displays updated data
        await mutateContainer();
        props.onChange?.({ name, moleculeType, sequence });
      } else if (seqinDigest?.composition) {
        console.log("Seqin digest was a coordinate file");
        const { name, moleculeType, sequence } = {
          name: `Chain_${seqinDigest.composition.peptides[0]}`,
          moleculeType: "PROTEIN",
          sequence:
            seqinDigest.sequences[seqinDigest.composition.peptides[0]] || "",
        };
        const sanitizedName = name.replace(/[^a-zA-Z0-9]/g, "_");
        await setPolymerType(moleculeType);
        await setName(sanitizedName);
        await setSequence(sequence);
        await setDescription(annotation);
        await setNCopies(1);
        // Clear intents so child components will sync their local state
        // with the new container data (otherwise wasRecentlyChanged blocks the sync)
        clearIntentForPath(`${item._objectPath}.polymerType`);
        clearIntentForPath(`${item._objectPath}.name`);
        clearIntentForPath(`${item._objectPath}.sequence`);
        clearIntentForPath(`${item._objectPath}.description`);
        clearIntentForPath(`${item._objectPath}.nCopies`);
        // Force a final container refresh to ensure modal displays updated data
        await mutateContainer();
        props.onChange?.({ name, moleculeType, sequence });
      }
    },
    [
      setSequence,
      setName,
      setPolymerType,
      setDescription,
      setNCopies,
      job,
      item,
      mutateContainer,
      clearIntentForPath,
    ]
  );

  return (
    <Card sx={{ border: "3px solid", borderColor: getValidationColor(item) }}>
      <CardHeader
        title={item._qualifiers.guiLabel}
        action={<ErrorInfo {...props} />}
      />
      <CardContent sx={{ my: 0, py: 0, pt: 2 }}>
        <Grid2 container rowSpacing={0} sx={{ my: 0, py: 0 }}>
          {item && (
            <Grid2 size={{ xs: 4 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ my: 0, py: 0, minWidth: "10rem" }}
                itemName={`${item._objectPath}.nCopies`}
                qualifiers={{
                  guiLabel: "nCopies",
                  enumerators: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                }}
              />
            </Grid2>
          )}
          {item && (
            <Grid2 size={{ xs: 4 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ my: 0, py: 0, minWidth: "10rem" }}
                itemName={`${item._objectPath}.polymerType`}
                qualifiers={{
                  guiLabel: "polymerType",
                }}
              />
            </Grid2>
          )}
          {item && (
            <Grid2 size={{ xs: 4 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ my: 0, py: 0, minWidth: "10rem" }}
                itemName={`${item._objectPath}.name`}
                qualifiers={{
                  guiLabel: "name",
                }}
              />
            </Grid2>
          )}
          {["description", "sequence"].map((key) => (
            <Grid2 key={key} size={{ xs: 12 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ my: 0, py: 0, minWidth: "calc(100% - 4rem)", mr: 2 }}
                itemName={`${item._objectPath}.${key}`}
                qualifiers={{
                  guiLabel: key,
                  guiMode: "multiLine",
                }}
              />
            </Grid2>
          ))}
          {["source"].map((key) => (
            <Grid2 key={key} size={{ xs: 12 }}>
              <CCP4i2TaskElement
                {...props}
                sx={{ my: 0, py: 0 }}
                itemName={`${item._objectPath}.${key}`}
                qualifiers={{
                  guiLabel: key,
                  guiMode: "multiLine",
                  mimeTypeName: "application/CCP4-seq",
                  downloadModes: ["uniprotFasta", "ebiPdb"],
                }}
                onChange={async (updatedItem: any) => {
                  console.log("Fetch file for param", updatedItem);
                  const { dbFileId, annotation } = valueOfItem(updatedItem);
                  const digest = await apiGet(
                    `files/${dbFileId}/digest_by_uuid`
                  );
                  console.log({ digest, annotation });
                  setSEQUENCEFromSEQIN(digest, annotation);
                }}
                suppressMutations={true}
              />
            </Grid2>
          ))}
        </Grid2>
      </CardContent>
    </Card>
  );
};
