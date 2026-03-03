import { useMemo } from "react";
import { useCCP4i2Window } from "../../../app-context";
import { Job } from "../../../types/models";
import { LinearProgress } from "@mui/material";
import AimlessPipeInterface from "./aimless-pipe";
import AuspexInterface from "./auspex";
import LidiaAcedrgNewInterface from "./LidiaAcedrgNew";
import ClustalwInterface from "./clustalw";
import Crank2Interface from "./crank2";
import GenericInterface from "./generic";
import ImportMergedInterface from "./import_merged";
import ModelcraftInterface from "./modelcraft";
import MolrepPipeInterface from "./molrep_pipe";
import MolrepSelfrot from "./molrep_selfrot";
import ProsmartRefmacInterface from "./prosmart-refmac";
import PhaserSimpleInterface from "./phaser_simple";
import PhaserRNPPipeline from "./phaser_rnp_pipeline";
import ServalcatPipeInterface from "./servalcat-pipe";
import SubstituteLigandInterface from "./SubstituteLigand";
import ProvideAsuContentsInterface from "./ProvideAsuContents";
import ProvideSequenceInterface from "./ProvideSequence";
import PhaserEPInterface from "./phaser_EP";
import PhaserEPLLGInterface from "./phaser_EP_LLG";
import PhaserPipelineInterface from "./phaser_pipeline";
import ParrotInterface from "./parrot";
import SHELXInterface from "./shelx";
import ShelxeMRInterface from "./shelxeMR";
import CSymmatchInterface from "./csymmatch";
import FreerFlagInterface from "./freerflag";
import GesamtInterface from "./gesamt";
import SplitMtzInterface from "./splitMtz";
import Xia2MultiplexInterface from "./xia2_multiplex";
import { useJob } from "../../../utils";

import AcornInterface from "./acorn";
import AMPLEInterface from "./ample";
import ArcimboldoInterface from "./arcimboldo";
import MakeLinkInterface from "./MakeLink";

// Auto-generated interfaces from legacy GUI files
import GeneratedAcedrgLinkInterface from "./generated/AcedrgLink";
import GeneratedAlternativeImportXIA2Interface from "./generated/AlternativeImportXIA2";
import GeneratedMakeProjectsAndDoLigandPipelineInterface from "./generated/MakeProjectsAndDoLigandPipeline";
import GeneratedPrepareDepositInterface from "./generated/PrepareDeposit";
import GeneratedProvideAlignmentInterface from "./generated/ProvideAlignment";
import SIMBADInterface from "./SIMBAD";
import GeneratedShelxCDInterface from "./generated/ShelxCD";
import GeneratedSubtractNativeInterface from "./generated/SubtractNative";
import GeneratedTestObsConversionsInterface from "./generated/TestObsConversions";
import GeneratedAddFractionalCoordsInterface from "./generated/add_fractional_coords";
import GeneratedAddingStatsToMmcifI2Interface from "./generated/adding_stats_to_mmcif_i2";
import GeneratedCcp4mgEditModelInterface from "./generated/ccp4mg_edit_model";
import GeneratedCcp4mgEditNomrbumpInterface from "./generated/ccp4mg_edit_nomrbump";
import GeneratedCcp4mgGeneralInterface from "./generated/ccp4mg_general";
import GeneratedChltofomInterface from "./generated/chltofom";
import GeneratedCmapcoeffInterface from "./generated/cmapcoeff";
import GeneratedComitInterface from "./generated/comit";
import GeneratedCoordinateSelectorInterface from "./generated/coordinate_selector";
import GeneratedCoot1Interface from "./generated/coot1";
import GeneratedCootRebuildInterface from "./generated/coot_rebuild";
import GeneratedCootRsrMorphInterface from "./generated/coot_rsr_morph";
import GeneratedCpattersonInterface from "./generated/cpatterson";
import GeneratedCtruncateInterface from "./generated/ctruncate";
import GeneratedDensityCalculatorInterface from "./generated/density_calculator";
import GeneratedDrMrModelbuildPipelineInterface from "./generated/dr_mr_modelbuild_pipeline";
import GeneratedEditbfacInterface from "./generated/editbfac";
import I2DimpleInterface from "./i2Dimple";
import GeneratedImportMosflmInterface from "./generated/import_mosflm";
import GeneratedLorestrI2Interface from "./generated/lorestr_i2";
import GeneratedMergeMtzInterface from "./generated/mergeMtz";
import GeneratedMetalCoordInterface from "./generated/metalCoord";
import GeneratedModelASUCheckInterface from "./generated/modelASUCheck";
import GeneratedMordaI2Interface from "./generated/morda_i2";
import MrBumpBasicInterface from "./mrbump_basic";
import GeneratedMosflmInterface from "./generated/mosflm";
import GeneratedPdbRedoApiInterface from "./generated/pdb_redo_api";
import GeneratedPdbviewEditInterface from "./generated/pdbview_edit";
import GeneratedPhaserEnsemblerInterface from "./generated/phaser_ensembler";
import GeneratedPhaserPhilInterface from "./generated/phaser_phil";
import GeneratedPisapipeInterface from "./generated/pisapipe";
import GeneratedPointlessReindexToMatchInterface from "./generated/pointless_reindexToMatch";
import GeneratedPrivateerInterface from "./generated/privateer";
import GeneratedQtpisaInterface from "./generated/qtpisa";
import GeneratedSheetbendInterface from "./generated/sheetbend";
import GeneratedValidateProteinInterface from "./generated/validate_protein";
import GeneratedXia2DialsInterface from "./generated/xia2_dials";

import GeneratedXia2SsxReduceInterface from "./generated/xia2_ssx_reduce";
import GeneratedZanudaInterface from "./generated/zanuda";

export interface CCP4i2TaskInterfaceProps {
  job: Job;
}

export interface TaskContainerProps {
  jobId?: number;
}

export const TaskContainer: React.FC<TaskContainerProps> = ({ jobId: propJobId }) => {
  const { jobId: contextJobId } = useCCP4i2Window();
  // Prefer prop jobId over context jobId for immediate rendering
  const jobId = propJobId ?? contextJobId;
  const { job, container } = useJob(jobId);

  const taskInterface = useMemo(() => {
    switch (job?.task_name) {
      case null:
        return <LinearProgress />;
      case "aimless_pipe":
        return <AimlessPipeInterface job={job} />;
      case "AUSPEX":
        return <AuspexInterface job={job} />;
      case "clustalw":
        return <ClustalwInterface job={job} />;
      case "crank2":
        return <Crank2Interface job={job} />;
      case "csymmatch":
        return <CSymmatchInterface job={job} />;
      case "freerflag":
        return <FreerFlagInterface job={job} />;
      case "gesamt":
        return <GesamtInterface job={job} />;
      case "import_merged":
        return <ImportMergedInterface job={job} />;
      case "LidiaAcedrgNew":
        return <LidiaAcedrgNewInterface job={job} />;
      case "modelcraft":
        return <ModelcraftInterface job={job} />;
      case "molrep_pipe":
        return <MolrepPipeInterface job={job} />;
      case "molrep_selfrot":
        return <MolrepSelfrot job={job} />;
      case "parrot":
        return <ParrotInterface job={job} />;
      case "phaser_EP":
        return <PhaserEPInterface job={job} />;
      case "phaser_EP_LLG":
        return <PhaserEPLLGInterface job={job} />;
      case "phaser_pipeline":
        return <PhaserPipelineInterface job={job} />;
      case "phaser_simple":
        return <PhaserSimpleInterface job={job} />;
      case "phaser_rnp_pipeline":
        return <PhaserRNPPipeline job={job} />;
      case "prosmart_refmac":
        return <ProsmartRefmacInterface job={job} />;
      case "ProvideAsuContents":
        return <ProvideAsuContentsInterface job={job} />;
      case "ProvideSequence":
        return <ProvideSequenceInterface job={job} />;
      case "servalcat_pipe":
        return <ServalcatPipeInterface job={job} />;
      case "SubstituteLigand":
        return <SubstituteLigandInterface job={job} />;
      case "shelx":
        return <SHELXInterface job={job} />;
      case "shelxeMR":
        return <ShelxeMRInterface job={job} />;
      case "splitMtz":
        return <SplitMtzInterface job={job} />;

      case "acorn":
        return <AcornInterface job={job} />;
      case "AMPLE":
        return <AMPLEInterface job={job} />;
      case "MakeLink":
        return <MakeLinkInterface job={job} />;

      // Auto-generated interfaces from legacy GUI files
      case "AcedrgLink":
        return <GeneratedAcedrgLinkInterface job={job} />;
      case "AlternativeImportXIA2":
        return <GeneratedAlternativeImportXIA2Interface job={job} />;
      case "MakeProjectsAndDoLigandPipeline":
        return <GeneratedMakeProjectsAndDoLigandPipelineInterface job={job} />;
      case "PrepareDeposit":
        return <GeneratedPrepareDepositInterface job={job} />;
      case "ProvideAlignment":
        return <GeneratedProvideAlignmentInterface job={job} />;
      case "SIMBAD":
        return <SIMBADInterface job={job} />;
      case "ShelxCD":
        return <GeneratedShelxCDInterface job={job} />;
      case "SubtractNative":
        return <GeneratedSubtractNativeInterface job={job} />;
      case "TestObsConversions":
        return <GeneratedTestObsConversionsInterface job={job} />;
      case "add_fractional_coords":
        return <GeneratedAddFractionalCoordsInterface job={job} />;
      case "adding_stats_to_mmcif_i2":
        return <GeneratedAddingStatsToMmcifI2Interface job={job} />;
      case "arcimboldo":
        return <ArcimboldoInterface job={job} />;
      case "ccp4mg_edit_model":
        return <GeneratedCcp4mgEditModelInterface job={job} />;
      case "ccp4mg_edit_nomrbump":
        return <GeneratedCcp4mgEditNomrbumpInterface job={job} />;
      case "ccp4mg_general":
        return <GeneratedCcp4mgGeneralInterface job={job} />;
      case "chltofom":
        return <GeneratedChltofomInterface job={job} />;
      case "cmapcoeff":
        return <GeneratedCmapcoeffInterface job={job} />;
      case "comit":
        return <GeneratedComitInterface job={job} />;
      case "coordinate_selector":
        return <GeneratedCoordinateSelectorInterface job={job} />;
      case "coot1":
        return <GeneratedCoot1Interface job={job} />;
      case "coot_rebuild":
        return <GeneratedCootRebuildInterface job={job} />;
      case "coot_rsr_morph":
        return <GeneratedCootRsrMorphInterface job={job} />;
      case "cpatterson":
        return <GeneratedCpattersonInterface job={job} />;
      case "ctruncate":
        return <GeneratedCtruncateInterface job={job} />;
      case "density_calculator":
        return <GeneratedDensityCalculatorInterface job={job} />;
      case "dr_mr_modelbuild_pipeline":
        return <GeneratedDrMrModelbuildPipelineInterface job={job} />;
      case "editbfac":
        return <GeneratedEditbfacInterface job={job} />;
      case "i2Dimple":
        return <I2DimpleInterface job={job} />;
      case "import_mosflm":
        return <GeneratedImportMosflmInterface job={job} />;
      case "lorestr_i2":
        return <GeneratedLorestrI2Interface job={job} />;
      case "mergeMtz":
        return <GeneratedMergeMtzInterface job={job} />;
      case "metalCoord":
        return <GeneratedMetalCoordInterface job={job} />;
      case "modelASUCheck":
        return <GeneratedModelASUCheckInterface job={job} />;
      case "morda_i2":
        return <GeneratedMordaI2Interface job={job} />;
      case "mrbump_basic":
        return <MrBumpBasicInterface job={job} />;
      case "mosflm":
        return <GeneratedMosflmInterface job={job} />;
      case "pdb_redo_api":
        return <GeneratedPdbRedoApiInterface job={job} />;
      case "pdbview_edit":
        return <GeneratedPdbviewEditInterface job={job} />;
      case "phaser_ensembler":
        return <GeneratedPhaserEnsemblerInterface job={job} />;
      case "phaser_phil":
        return <GeneratedPhaserPhilInterface job={job} />;
      case "pisapipe":
        return <GeneratedPisapipeInterface job={job} />;
      case "pointless_reindexToMatch":
        return <GeneratedPointlessReindexToMatchInterface job={job} />;
      case "privateer":
        return <GeneratedPrivateerInterface job={job} />;
      case "qtpisa":
        return <GeneratedQtpisaInterface job={job} />;
      case "sheetbend":
        return <GeneratedSheetbendInterface job={job} />;
      case "validate_protein":
        return <GeneratedValidateProteinInterface job={job} />;
      case "xia2_dials":
        return <GeneratedXia2DialsInterface job={job} />;
      case "xia2_multiplex":
        return <Xia2MultiplexInterface job={job} />;
      case "xia2_ssx_reduce":
        return <GeneratedXia2SsxReduceInterface job={job} />;
      case "zanuda":
        return <GeneratedZanudaInterface job={job} />;

      default:
        return job && <GenericInterface job={job} />;
    }
  }, [job, container]);

  if (!jobId) return <LinearProgress />;
  if (!container) return <LinearProgress />;
  if (!job) return <LinearProgress />;

  return <>{taskInterface}</>;
};
