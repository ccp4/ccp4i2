/**
 * Does a project file satisfy a field's requiredContentFlag?
 *
 * Only observed data (CObsDataFile, server CCP4XtalData.CAN_CONVERT_TO)
 * defines conversions between content flags: IPAIR=1, FPAIR=2, IMEAN=3,
 * FMEAN=4, and an IPAIR file can serve a field that asks for FMEAN. Every
 * other class numbers its flags differently -- for coordinates 1 is PDB
 * format and 2 is mmCIF -- so for them the flag must match exactly. Applying
 * the observed-data table to every class let a PDB-format model through a
 * filter that asked for mmCIF.
 */
export const OBSERVED_DATA_TYPE = "application/CCP4-mtz-observed";

const OBSERVED_CAN_CONVERT_TO: Record<number, number[]> = {
  1: [1, 2, 3, 4], // IPAIR -> IPAIR, FPAIR, IMEAN, FMEAN
  2: [2, 4], // FPAIR -> FPAIR, FMEAN
  3: [3, 4], // IMEAN -> IMEAN, FMEAN
  4: [4], // FMEAN -> FMEAN
};

export const canConvertToRequired = (
  fileContent: number | null | undefined,
  requiredFlags: number[],
  fileType?: string | null
): boolean => {
  if (fileContent == null) return false;
  if (fileType !== OBSERVED_DATA_TYPE) return requiredFlags.includes(fileContent);
  const convertibleTo = OBSERVED_CAN_CONVERT_TO[fileContent];
  if (!convertibleTo) return false;
  return requiredFlags.some((required) => convertibleTo.includes(required));
};
