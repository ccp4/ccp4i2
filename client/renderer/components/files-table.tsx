import {
  Avatar,
  IconButton,
  LinearProgress,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  Tooltip,
} from "@mui/material";
import { Delete, Download } from "@mui/icons-material";
import { File } from "../types/models";
import { fileSize } from "../pipes";
import { useApi } from "../api";
import { useDeleteDialog } from "../providers/delete-dialog";
import { KeyedMutator } from "swr";
import { useMemo } from "react";

export const fileTypeMapping: any = {
  Unknown: "DataFile",
  "application/CCP4-seq": "SeqDataFile",
  "chemical/x-pdb": "PdbDataFile",
  MultiPDB: "",
  "application/CCP4-mtz": "MtzDataFile",
  "application/CCP4-unmerged-mtz": "MtzDataFile",
  "application/CCP4-unmerged-experimental": "UnmergedDataFile",
  "application/CCP4-map": "MapDataFile",
  "application/refmac-dictionary": "DictDataFile",
  "application/refmac-TLS": "TLSDataFile",
  "application/CCP4-mtz-freerflag": "FreeRDataFile",
  "application/CCP4-mtz-observed": "ObsDataFile",
  "application/CCP4-mtz-phases": "PhsDataFile",
  "application/CCP4-mtz-map": "MapCoeffsDataFile",
  Dummy: "",
  "application/CCP4-seqalign": "SeqAlignDataFile",
  "application/CCP4-mtz-mini": "MiniMtzDataFile",
  "application/coot-script": "CootHistoryDataFile",
  "application/refmac-external-restraints": "RefmacRestraintsDataFile",
  "application/CCP4-scene": "SceneDataFile",
  "application/CCP4-shelx-FA": "ShelxFADataFile",
  "application/phaser-sol": "PhaserSolDataFile",
  "chemical/x-mdl-molfile": "MDLMolDataFile",
  "application/iMosflm-xml": "ImosflmXmlDataFile",
  "application/CCP4-image": "ImageFile",
  "application/CCP4-generic-reflections": "GenericReflDataFile",
  "application/HHPred-alignments": "HhpredDataFile",
  "application/Blast-alignments": "BlastDataFile",
  "chemical/x-pdb-ensemble": "EnsemblePdbDataFile",
  "application/CCP4-asu-content": "AsuDataFile",
  "application/dials-jfile": "DialsJsonFile",
  "application/dials-pfile": "DialsPickleFile",
  "application/phaser-rfile": "PhaserRFileDataFile",
  "application/refmac-keywords": "RefmacKeywordFile",
  "application/x-pdf": "PDFDataFile",
  "application/postscript": "PostscriptDataFile",
  "application/EBI-validation-xml": "EBIValidationXMLDataFile",
  "chemical/x-cif": "MmcifReflDataFile",
};

export default function FilesTable({
  files,
  mutate,
}: {
  files: File[] | undefined;
  mutate: KeyedMutator<File[]>;
}) {
  const deleteDialog = useDeleteDialog();
  const api = useApi();
  const fileTypeIcon = (type: string) => {
    return Object.keys(fileTypeMapping).includes(type)
      ? fileTypeMapping[type]
      : "ccp4";
  };

  if (files === undefined) return <LinearProgress />;
  if (files.length === 0) return <></>;

  function handleDelete(file: File) {
    if (deleteDialog)
      deleteDialog({
        type: "show",
        what: file.name,
        onDelete: () => {
          api.delete(`files/${file.id}`).then(() => mutate());
        },
      });
  }

  return (
    <Table size="small">
      <TableHead>
        <TableRow>
          <TableCell>Type</TableCell>
          <TableCell>Name</TableCell>
          <TableCell>Size</TableCell>
          <TableCell></TableCell>
        </TableRow>
      </TableHead>
      <TableBody>
        {files.map((file: File) => (
          <TableRow key={file.id}>
            <TableCell title={file.type}>
              <Avatar
                src={`/api/proxy/djangostatic/qticons/${fileTypeIcon(
                  file.type
                )}.png`}
                sx={{ height: "1.5rem", width: "1.5rem" }}
              />
            </TableCell>
            <TableCell>{file.name}</TableCell>
            {/*<TableCell>{fileSize(file.size)}</TableCell>*/}
            <TableCell>
              {/*<Tooltip title="Export file">
                <IconButton href={file.file} download>
                  <Download />
                </IconButton>
              </Tooltip>*/}
              <Tooltip title="Delete file">
                <IconButton onClick={() => handleDelete(file)}>
                  <Delete />
                </IconButton>
              </Tooltip>
            </TableCell>
          </TableRow>
        ))}
      </TableBody>
    </Table>
  );
}
