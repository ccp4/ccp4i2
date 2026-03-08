import { useMemo } from "react";
import { useCCP4i2Window } from "../../../app-context";
import { Job } from "../../../types/models";
import { LinearProgress } from "@mui/material";
import { useJob } from "../../../utils";

import AcedrgLinkInterface from "./AcedrgLink";
import AcornInterface from "./acorn";
import AddFractionalCoordsInterface from "./add_fractional_coords";
import AddingStatsToMmcifI2Interface from "./adding_stats_to_mmcif_i2";
import AimlessPipeInterface from "./aimless-pipe";
import AlternativeImportXIA2Interface from "./AlternativeImportXIA2";
import AMPLEInterface from "./ample";
import ArcimboldoInterface from "./arcimboldo";
import ArpWarpClassicInterface from "./arp_warp_classic";
import AuspexInterface from "./auspex";
import BuccaneerBuildRefineMrInterface from "./buccaneer_build_refine_mr";
import BusterInterface from "./buster";
import Ccp4mgEditModelInterface from "./ccp4mg_edit_model";
import Ccp4mgEditNomrbumpInterface from "./ccp4mg_edit_nomrbump";
import Ccp4mgGeneralInterface from "./ccp4mg_general";
import ChainsawInterface from "./chainsaw";
import ChltofomInterface from "./chltofom";
import Cif2mtzInterface from "./cif2mtz";
import ClustalwInterface from "./clustalw";
import CmapcoeffInterface from "./cmapcoeff";
import ComitInterface from "./comit";
import Coot1Interface from "./coot1";
import CootRebuildInterface from "./coot_rebuild";
import CootRsrMorphInterface from "./coot_rsr_morph";
import CoordinateSelectorInterface from "./coordinate_selector";
import CpattersonInterface from "./cpatterson";
import Crank2Interface from "./crank2";
import CSymmatchInterface from "./csymmatch";
import CtruncateInterface from "./ctruncate";
import DensityCalculatorInterface from "./density_calculator";
import DialsImageInterface from "./dials_image";
import DialsRlatticeInterface from "./dials_rlattice";
import DrMrModelbuildPipelineInterface from "./dr_mr_modelbuild_pipeline";
import DuiInterface from "./dui";
import EditbfacInterface from "./editbfac";
import EdstatsInterface from "./edstats";
import FindmyseqInterface from "./findmyseq";
import FreerFlagInterface from "./freerflag";
import GenericInterface from "./generic";
import GesamtInterface from "./gesamt";
import I2DimpleInterface from "./i2Dimple";
import ImosflmInterface from "./imosflm";
import ImportMergedInterface from "./import_merged";
import ImportMosflmInterface from "./import_mosflm";
import ImportSerialInterface from "./import_serial";
import ImportSerialPipeInterface from "./import_serial_pipe";
import ImportXia2Interface from "./import_xia2";
import LidiaAcedrgNewInterface from "./LidiaAcedrgNew";
import LorestrI2Interface from "./lorestr_i2";
import MakeLinkInterface from "./MakeLink";
import MakeMonsterInterface from "./MakeMonster";
import MakeProjectsAndDoLigandPipelineInterface from "./MakeProjectsAndDoLigandPipeline";
import MatthewsInterface from "./matthews";
import MergeMtzInterface from "./mergeMtz";
import MetalCoordInterface from "./metalCoord";
import ModelASUCheckInterface from "./modelASUCheck";
import ModelcraftInterface from "./modelcraft";
import MolrepDenInterface from "./molrep_den";
import MolrepPipeInterface from "./molrep_pipe";
import MolrepSelfrot from "./molrep_selfrot";
import MordaI2Interface from "./morda_i2";
import MosflmInterface from "./mosflm";
import MrBumpBasicInterface from "./mrbump_basic";
import MrparseInterface from "./mrparse";
import MtzutilsInterface from "./mtzutils";
import NautilusBuildRefineInterface from "./nautilus_build_refine";
import NewProjectFromMergedInterface from "./newProject_fromMerged";
import PairefInterface from "./pairef";
import ParrotInterface from "./parrot";
import PdbRedoApiInterface from "./pdb_redo_api";
import PdbsetUiInterface from "./pdbset_ui";
import PdbviewEditInterface from "./pdbview_edit";
import PhaserEnsemblerInterface from "./phaser_ensembler";
import PhaserEPInterface from "./phaser_EP";
import PhaserEPLLGInterface from "./phaser_EP_LLG";
import PhaserMrInterface from "./phaser_mr";
import PhaserPhilInterface from "./phaser_phil";
import PhaserPipelineInterface from "./phaser_pipeline";
import PhaserRNPPipeline from "./phaser_rnp_pipeline";
import PhaserSimpleInterface from "./phaser_simple";
import PhaserSingleMRInterface from "./phaser_singleMR";
import PhasertngPicardInterface from "./phasertng_picard";
import PisapipeInterface from "./pisapipe";
import PointlessReindexToMatchInterface from "./pointless_reindexToMatch";
import PrepareDepositInterface from "./PrepareDeposit";
import PrivateerInterface from "./privateer";
import ProsmartInterface from "./prosmart";
import ProsmartRefmacInterface from "./prosmart-refmac";
import ProvideAlignmentInterface from "./ProvideAlignment";
import ProvideAsuContentsInterface from "./ProvideAsuContents";
import ProvideSequenceInterface from "./ProvideSequence";
import ProvideTLSInterface from "./ProvideTLS";
import PyphaserMrInterface from "./pyphaser_mr";
import QtpisaInterface from "./qtpisa";
import ScaleitInterface from "./scaleit";
import SculptorInterface from "./sculptor";
import ServalcatPipeInterface from "./servalcat-pipe";
import SheetbendInterface from "./sheetbend";
import ShelxCDInterface from "./ShelxCD";
import SHELXInterface from "./shelx";
import ShelxeMRInterface from "./shelxeMR";
import SIMBADInterface from "./SIMBAD";
import SlicendiceInterface from "./slicendice";
import SplitMtzInterface from "./splitMtz";
import SubstituteLigandInterface from "./SubstituteLigand";
import SubtractNativeInterface from "./SubtractNative";
import TableoneInterface from "./tableone";
import TestObsConversionsInterface from "./TestObsConversions";
import UniqueInterface from "./unique";
import ValidateProteinInterface from "./validate_protein";
import Xia2DialsInterface from "./xia2_dials";
import Xia2MultiplexInterface from "./xia2_multiplex";
import Xia2SsxReduceInterface from "./xia2_ssx_reduce";
import ZanudaInterface from "./zanuda";

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
      case "AcedrgLink":
        return <AcedrgLinkInterface job={job} />;
      case "acorn":
        return <AcornInterface job={job} />;
      case "add_fractional_coords":
        return <AddFractionalCoordsInterface job={job} />;
      case "adding_stats_to_mmcif_i2":
        return <AddingStatsToMmcifI2Interface job={job} />;
      case "aimless_pipe":
        return <AimlessPipeInterface job={job} />;
      case "AlternativeImportXIA2":
        return <AlternativeImportXIA2Interface job={job} />;
      case "AMPLE":
        return <AMPLEInterface job={job} />;
      case "arcimboldo":
        return <ArcimboldoInterface job={job} />;
      case "arp_warp_classic":
        return <ArpWarpClassicInterface job={job} />;
      case "AUSPEX":
        return <AuspexInterface job={job} />;
      case "buccaneer_build_refine_mr":
        return <BuccaneerBuildRefineMrInterface job={job} />;
      case "buster":
        return <BusterInterface job={job} />;
      case "ccp4mg_edit_model":
        return <Ccp4mgEditModelInterface job={job} />;
      case "ccp4mg_edit_nomrbump":
        return <Ccp4mgEditNomrbumpInterface job={job} />;
      case "ccp4mg_general":
        return <Ccp4mgGeneralInterface job={job} />;
      case "chainsaw":
        return <ChainsawInterface job={job} />;
      case "chltofom":
        return <ChltofomInterface job={job} />;
      case "cif2mtz":
        return <Cif2mtzInterface job={job} />;
      case "clustalw":
        return <ClustalwInterface job={job} />;
      case "cmapcoeff":
        return <CmapcoeffInterface job={job} />;
      case "comit":
        return <ComitInterface job={job} />;
      case "coot1":
        return <Coot1Interface job={job} />;
      case "coot_rebuild":
        return <CootRebuildInterface job={job} />;
      case "coot_rsr_morph":
        return <CootRsrMorphInterface job={job} />;
      case "coordinate_selector":
        return <CoordinateSelectorInterface job={job} />;
      case "cpatterson":
        return <CpattersonInterface job={job} />;
      case "crank2":
        return <Crank2Interface job={job} />;
      case "csymmatch":
        return <CSymmatchInterface job={job} />;
      case "ctruncate":
        return <CtruncateInterface job={job} />;
      case "density_calculator":
        return <DensityCalculatorInterface job={job} />;
      case "dials_image":
        return <DialsImageInterface job={job} />;
      case "dials_rlattice":
        return <DialsRlatticeInterface job={job} />;
      case "dr_mr_modelbuild_pipeline":
        return <DrMrModelbuildPipelineInterface job={job} />;
      case "dui":
        return <DuiInterface job={job} />;
      case "editbfac":
        return <EditbfacInterface job={job} />;
      case "edstats":
        return <EdstatsInterface job={job} />;
      case "findmyseq":
        return <FindmyseqInterface job={job} />;
      case "freerflag":
        return <FreerFlagInterface job={job} />;
      case "gesamt":
        return <GesamtInterface job={job} />;
      case "i2Dimple":
        return <I2DimpleInterface job={job} />;
      case "imosflm":
        return <ImosflmInterface job={job} />;
      case "import_merged":
        return <ImportMergedInterface job={job} />;
      case "import_mosflm":
        return <ImportMosflmInterface job={job} />;
      case "import_serial":
        return <ImportSerialInterface job={job} />;
      case "import_serial_pipe":
        return <ImportSerialPipeInterface job={job} />;
      case "import_xia2":
        return <ImportXia2Interface job={job} />;
      case "LidiaAcedrgNew":
        return <LidiaAcedrgNewInterface job={job} />;
      case "lorestr_i2":
        return <LorestrI2Interface job={job} />;
      case "MakeLink":
        return <MakeLinkInterface job={job} />;
      case "MakeMonster":
        return <MakeMonsterInterface job={job} />;
      case "MakeProjectsAndDoLigandPipeline":
        return <MakeProjectsAndDoLigandPipelineInterface job={job} />;
      case "matthews":
        return <MatthewsInterface job={job} />;
      case "mergeMtz":
        return <MergeMtzInterface job={job} />;
      case "metalCoord":
        return <MetalCoordInterface job={job} />;
      case "modelASUCheck":
        return <ModelASUCheckInterface job={job} />;
      case "modelcraft":
        return <ModelcraftInterface job={job} />;
      case "molrep_den":
        return <MolrepDenInterface job={job} />;
      case "molrep_pipe":
        return <MolrepPipeInterface job={job} />;
      case "molrep_selfrot":
        return <MolrepSelfrot job={job} />;
      case "morda_i2":
        return <MordaI2Interface job={job} />;
      case "mosflm":
        return <MosflmInterface job={job} />;
      case "mrbump_basic":
        return <MrBumpBasicInterface job={job} />;
      case "mrparse":
        return <MrparseInterface job={job} />;
      case "mtzutils":
        return <MtzutilsInterface job={job} />;
      case "nautilus_build_refine":
        return <NautilusBuildRefineInterface job={job} />;
      case "newProject_fromMerged":
        return <NewProjectFromMergedInterface job={job} />;
      case "pairef":
        return <PairefInterface job={job} />;
      case "parrot":
        return <ParrotInterface job={job} />;
      case "pdb_redo_api":
        return <PdbRedoApiInterface job={job} />;
      case "pdbset_ui":
        return <PdbsetUiInterface job={job} />;
      case "pdbview_edit":
        return <PdbviewEditInterface job={job} />;
      case "phaser_ensembler":
        return <PhaserEnsemblerInterface job={job} />;
      case "phaser_EP":
        return <PhaserEPInterface job={job} />;
      case "phaser_EP_LLG":
        return <PhaserEPLLGInterface job={job} />;
      case "phaser_mr":
        return <PhaserMrInterface job={job} />;
      case "phaser_phil":
        return <PhaserPhilInterface job={job} />;
      case "phaser_pipeline":
        return <PhaserPipelineInterface job={job} />;
      case "phaser_rnp_pipeline":
        return <PhaserRNPPipeline job={job} />;
      case "phaser_simple":
        return <PhaserSimpleInterface job={job} />;
      case "phaser_singleMR":
        return <PhaserSingleMRInterface job={job} />;
      case "phasertng_picard":
        return <PhasertngPicardInterface job={job} />;
      case "pisapipe":
        return <PisapipeInterface job={job} />;
      case "pointless_reindexToMatch":
        return <PointlessReindexToMatchInterface job={job} />;
      case "PrepareDeposit":
        return <PrepareDepositInterface job={job} />;
      case "privateer":
        return <PrivateerInterface job={job} />;
      case "prosmart":
        return <ProsmartInterface job={job} />;
      case "prosmart_refmac":
        return <ProsmartRefmacInterface job={job} />;
      case "ProvideAlignment":
        return <ProvideAlignmentInterface job={job} />;
      case "ProvideAsuContents":
        return <ProvideAsuContentsInterface job={job} />;
      case "ProvideSequence":
        return <ProvideSequenceInterface job={job} />;
      case "ProvideTLS":
        return <ProvideTLSInterface job={job} />;
      case "pyphaser_mr":
        return <PyphaserMrInterface job={job} />;
      case "qtpisa":
        return <QtpisaInterface job={job} />;
      case "scaleit":
        return <ScaleitInterface job={job} />;
      case "sculptor":
        return <SculptorInterface job={job} />;
      case "servalcat_pipe":
        return <ServalcatPipeInterface job={job} />;
      case "sheetbend":
        return <SheetbendInterface job={job} />;
      case "ShelxCD":
        return <ShelxCDInterface job={job} />;
      case "shelx":
        return <SHELXInterface job={job} />;
      case "shelxeMR":
        return <ShelxeMRInterface job={job} />;
      case "SIMBAD":
        return <SIMBADInterface job={job} />;
      case "slicendice":
        return <SlicendiceInterface job={job} />;
      case "splitMtz":
        return <SplitMtzInterface job={job} />;
      case "SubstituteLigand":
        return <SubstituteLigandInterface job={job} />;
      case "SubtractNative":
        return <SubtractNativeInterface job={job} />;
      case "tableone":
        return <TableoneInterface job={job} />;
      case "TestObsConversions":
        return <TestObsConversionsInterface job={job} />;
      case "unique":
        return <UniqueInterface job={job} />;
      case "validate_protein":
        return <ValidateProteinInterface job={job} />;
      case "xia2_dials":
        return <Xia2DialsInterface job={job} />;
      case "xia2_multiplex":
        return <Xia2MultiplexInterface job={job} />;
      case "xia2_ssx_reduce":
        return <Xia2SsxReduceInterface job={job} />;
      case "zanuda":
        return <ZanudaInterface job={job} />;

      default:
        return job && <GenericInterface job={job} />;
    }
  }, [job, container]);

  if (!jobId) return <LinearProgress />;
  if (!container) return <LinearProgress />;
  if (!job) return <LinearProgress />;

  return <>{taskInterface}</>;
};
