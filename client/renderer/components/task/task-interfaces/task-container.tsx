import { useMemo } from "react";
import { useCCP4i2Window } from "../../../app-context";
import { Job } from "../../../types/models";
import { LinearProgress } from "@mui/material";
import AimlessPipeInterface from "./aimless_pipe";
import LidiaAcedrgNewInterface from "./LidiaAcedrgNew";
import ClustalwInterface from "./clustalw";
import Crank2Interface from "./crank2";
import GenericInterface from "./generic";
import ImportMergedInterface from "./import_merged";
import ModelcraftInterface from "./modelcraft";
import MolrepSelfrot from "./molrep_selfrot";
import ProsmartRefmacInterface from "./prosmart_refmac";
import PhaserSimpleInterface from "./phaser_simple";
import PhaserRNPPipeline from "./phaser_rnp_pipeline";
import ServalcatPipeInterface from "./servalcat_pipe";
import SubstituteLigandInterface from "./SubstituteLigand";
import ProvideAsuContentsInterface from "./ProvideAsuContents";
import ProvideSequenceInterface from "./ProvideSequence";
import PhaserEPInterface from "./phaser_EP";
import PhaserEPLLGInterface from "./phaser_EP_LLG";
import PhaserPipelineInterface from "./phaser_pipeline";
import ParrotInterface from "./parrot";
import SHELXInterface from "./shelx";
import CSymmatchInterface from "./csymmatch";
import { useJob } from "../../../utils";

// Auto-generated interfaces from legacy GUI files
import GeneratedAMPLEInterface from "./generated/AMPLE";
import GeneratedAcedrgLinkInterface from "./generated/AcedrgLink";
import GeneratedAlternativeImportXIA2Interface from "./generated/AlternativeImportXIA2";
import GeneratedLidiaAcedrgInterface from "./generated/LidiaAcedrg";
import GeneratedMakeLinkInterface from "./generated/MakeLink";
import GeneratedMakeProjectsAndDoLigandPipelineInterface from "./generated/MakeProjectsAndDoLigandPipeline";
import GeneratedPrepareDepositInterface from "./generated/PrepareDeposit";
import GeneratedProvideAlignmentInterface from "./generated/ProvideAlignment";
import GeneratedSIMBADInterface from "./generated/SIMBAD";
import GeneratedShelxCDInterface from "./generated/ShelxCD";
import GeneratedShelxCEInterface from "./generated/ShelxCE";
import GeneratedSubtractNativeInterface from "./generated/SubtractNative";
import GeneratedTestObsConversionsInterface from "./generated/TestObsConversions";
import GeneratedAddFractionalCoordsInterface from "./generated/add_fractional_coords";
import GeneratedAddingStatsToMmcifI2Interface from "./generated/adding_stats_to_mmcif_i2";
import GeneratedArcimboldo from "./generated/arcimboldo";
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
import GeneratedI2DimpleInterface from "./generated/i2Dimple";
import GeneratedImportMosflmInterface from "./generated/import_mosflm";
import GeneratedLorestrI2Interface from "./generated/lorestr_i2";
import GeneratedMergeMtzInterface from "./generated/mergeMtz";
import GeneratedMetalCoordInterface from "./generated/metalCoord";
import GeneratedModelASUCheckInterface from "./generated/modelASUCheck";
import GeneratedMordaI2Interface from "./generated/morda_i2";
import GeneratedMosflmInterface from "./generated/mosflm";
import GeneratedNautilusBuildRefineInterface from "./generated/nautilus_build_refine";
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
import GeneratedXia2MultiplexInterface from "./generated/xia2_multiplex";
import GeneratedXia2SsxReduceInterface from "./generated/xia2_ssx_reduce";
import GeneratedZanudaInterface from "./generated/zanuda";

export interface CCP4i2TaskInterfaceProps {
  job: Job;
}

export const TaskContainer = () => {
  const { jobId } = useCCP4i2Window();
  const { job, container } = useJob(jobId);

  const taskInterface = useMemo(() => {
    switch (job?.task_name) {
      case null:
        return <LinearProgress />;
      case "aimless_pipe":
        return <AimlessPipeInterface job={job} />;
      case "clustalw":
        return <ClustalwInterface job={job} />;
      case "crank2":
        return <Crank2Interface job={job} />;
      case "csymmatch":
        return <CSymmatchInterface job={job} />;
      case "import_merged":
        return <ImportMergedInterface job={job} />;
      case "LidiaAcedrgNew":
        return <LidiaAcedrgNewInterface job={job} />;
      case "modelcraft":
        return <ModelcraftInterface job={job} />;
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

      // Auto-generated interfaces from legacy GUI files
      case "AMPLE":
        return <GeneratedAMPLEInterface job={job} />;
      case "AcedrgLink":
        return <GeneratedAcedrgLinkInterface job={job} />;
      case "AlternativeImportXIA2":
        return <GeneratedAlternativeImportXIA2Interface job={job} />;
      case "LidiaAcedrg":
        return <GeneratedLidiaAcedrgInterface job={job} />;
      case "MakeLink":
        return <GeneratedMakeLinkInterface job={job} />;
      case "MakeProjectsAndDoLigandPipeline":
        return <GeneratedMakeProjectsAndDoLigandPipelineInterface job={job} />;
      case "PrepareDeposit":
        return <GeneratedPrepareDepositInterface job={job} />;
      case "ProvideAlignment":
        return <GeneratedProvideAlignmentInterface job={job} />;
      case "SIMBAD":
        return <GeneratedSIMBADInterface job={job} />;
      case "ShelxCD":
        return <GeneratedShelxCDInterface job={job} />;
      case "ShelxCE":
        return <GeneratedShelxCEInterface job={job} />;
      case "SubtractNative":
        return <GeneratedSubtractNativeInterface job={job} />;
      case "TestObsConversions":
        return <GeneratedTestObsConversionsInterface job={job} />;
      case "add_fractional_coords":
        return <GeneratedAddFractionalCoordsInterface job={job} />;
      case "adding_stats_to_mmcif_i2":
        return <GeneratedAddingStatsToMmcifI2Interface job={job} />;
      case "arcimboldo":
        return <GeneratedArcimboldo job={job} />;
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
        return <GeneratedI2DimpleInterface job={job} />;
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
      case "mosflm":
        return <GeneratedMosflmInterface job={job} />;
      case "nautilus_build_refine":
        return <GeneratedNautilusBuildRefineInterface job={job} />;
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
        return <GeneratedXia2MultiplexInterface job={job} />;
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
