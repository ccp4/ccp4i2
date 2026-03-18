/**
 * Task interface lookup and container component.
 *
 * Maps task names to their custom React interface components.
 * Tasks not listed here fall back to GenericInterface.
 *
 * See also: lib/task-registry.ts for category/ordering in the task chooser.
 */
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
import PhaserEPAUTOInterface from "./phaser_EP_AUTO";
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

/**
 * Custom task interface components, keyed by task name.
 * Tasks not listed here will render with GenericInterface.
 *
 * To add a new custom interface:
 *   1. Create the component file in this directory
 *   2. Import it above
 *   3. Add one line to this record: "task_name": TaskNameInterface,
 */
const TASK_INTERFACES: Record<
  string,
  React.ComponentType<CCP4i2TaskInterfaceProps>
> = {
  AcedrgLink: AcedrgLinkInterface,
  acorn: AcornInterface,
  add_fractional_coords: AddFractionalCoordsInterface,
  adding_stats_to_mmcif_i2: AddingStatsToMmcifI2Interface,
  aimless_pipe: AimlessPipeInterface,
  AlternativeImportXIA2: AlternativeImportXIA2Interface,
  AMPLE: AMPLEInterface,
  arcimboldo: ArcimboldoInterface,
  arp_warp_classic: ArpWarpClassicInterface,
  AUSPEX: AuspexInterface,
  buccaneer_build_refine_mr: BuccaneerBuildRefineMrInterface,
  buster: BusterInterface,
  ccp4mg_edit_model: Ccp4mgEditModelInterface,
  ccp4mg_edit_nomrbump: Ccp4mgEditNomrbumpInterface,
  ccp4mg_general: Ccp4mgGeneralInterface,
  chainsaw: ChainsawInterface,
  chltofom: ChltofomInterface,
  cif2mtz: Cif2mtzInterface,
  clustalw: ClustalwInterface,
  cmapcoeff: CmapcoeffInterface,
  comit: ComitInterface,
  coot1: Coot1Interface,
  coot_rebuild: CootRebuildInterface,
  coot_rsr_morph: CootRsrMorphInterface,
  coordinate_selector: CoordinateSelectorInterface,
  cpatterson: CpattersonInterface,
  crank2: Crank2Interface,
  csymmatch: CSymmatchInterface,
  ctruncate: CtruncateInterface,
  density_calculator: DensityCalculatorInterface,
  dials_image: DialsImageInterface,
  dials_rlattice: DialsRlatticeInterface,
  dr_mr_modelbuild_pipeline: DrMrModelbuildPipelineInterface,
  dui: DuiInterface,
  editbfac: EditbfacInterface,
  edstats: EdstatsInterface,
  findmyseq: FindmyseqInterface,
  freerflag: FreerFlagInterface,
  gesamt: GesamtInterface,
  i2Dimple: I2DimpleInterface,
  imosflm: ImosflmInterface,
  import_merged: ImportMergedInterface,
  import_mosflm: ImportMosflmInterface,
  import_serial: ImportSerialInterface,
  import_serial_pipe: ImportSerialPipeInterface,
  import_xia2: ImportXia2Interface,
  LidiaAcedrgNew: LidiaAcedrgNewInterface,
  lorestr_i2: LorestrI2Interface,
  MakeLink: MakeLinkInterface,
  MakeMonster: MakeMonsterInterface,
  MakeProjectsAndDoLigandPipeline: MakeProjectsAndDoLigandPipelineInterface,
  matthews: MatthewsInterface,
  mergeMtz: MergeMtzInterface,
  metalCoord: MetalCoordInterface,
  modelASUCheck: ModelASUCheckInterface,
  modelcraft: ModelcraftInterface,
  molrep_den: MolrepDenInterface,
  molrep_pipe: MolrepPipeInterface,
  molrep_selfrot: MolrepSelfrot,
  morda_i2: MordaI2Interface,
  mosflm: MosflmInterface,
  mrbump_basic: MrBumpBasicInterface,
  mrparse: MrparseInterface,
  mtzutils: MtzutilsInterface,
  nautilus_build_refine: NautilusBuildRefineInterface,
  newProject_fromMerged: NewProjectFromMergedInterface,
  pairef: PairefInterface,
  parrot: ParrotInterface,
  pdb_redo_api: PdbRedoApiInterface,
  pdbset_ui: PdbsetUiInterface,
  pdbview_edit: PdbviewEditInterface,
  phaser_ensembler: PhaserEnsemblerInterface,
  phaser_EP: PhaserEPInterface,
  phaser_EP_AUTO: PhaserEPAUTOInterface,
  phaser_EP_LLG: PhaserEPLLGInterface,
  phaser_mr: PhaserMrInterface,
  phaser_phil: PhaserPhilInterface,
  phaser_pipeline: PhaserPipelineInterface,
  phaser_rnp_pipeline: PhaserRNPPipeline,
  phaser_simple: PhaserSimpleInterface,
  phaser_singleMR: PhaserSingleMRInterface,
  phasertng_picard: PhasertngPicardInterface,
  pisapipe: PisapipeInterface,
  pointless_reindexToMatch: PointlessReindexToMatchInterface,
  PrepareDeposit: PrepareDepositInterface,
  privateer: PrivateerInterface,
  prosmart: ProsmartInterface,
  prosmart_refmac: ProsmartRefmacInterface,
  ProvideAlignment: ProvideAlignmentInterface,
  ProvideAsuContents: ProvideAsuContentsInterface,
  ProvideSequence: ProvideSequenceInterface,
  ProvideTLS: ProvideTLSInterface,
  qtpisa: QtpisaInterface,
  scaleit: ScaleitInterface,
  sculptor: SculptorInterface,
  servalcat_pipe: ServalcatPipeInterface,
  sheetbend: SheetbendInterface,
  ShelxCD: ShelxCDInterface,
  shelx: SHELXInterface,
  shelxeMR: ShelxeMRInterface,
  SIMBAD: SIMBADInterface,
  slicendice: SlicendiceInterface,
  splitMtz: SplitMtzInterface,
  SubstituteLigand: SubstituteLigandInterface,
  SubtractNative: SubtractNativeInterface,
  tableone: TableoneInterface,
  TestObsConversions: TestObsConversionsInterface,
  unique: UniqueInterface,
  validate_protein: ValidateProteinInterface,
  xia2_dials: Xia2DialsInterface,
  xia2_multiplex: Xia2MultiplexInterface,
  xia2_ssx_reduce: Xia2SsxReduceInterface,
  zanuda: ZanudaInterface,
};

export interface TaskContainerProps {
  jobId?: number;
}

export const TaskContainer: React.FC<TaskContainerProps> = ({ jobId: propJobId }) => {
  const { jobId: contextJobId } = useCCP4i2Window();
  // Prefer prop jobId over context jobId for immediate rendering
  const jobId = propJobId ?? contextJobId;
  const { job, container } = useJob(jobId);

  const taskInterface = useMemo(() => {
    if (!job) return <LinearProgress />;
    const Component = TASK_INTERFACES[job.task_name] ?? GenericInterface;
    return <Component job={job} />;
  }, [job, container]);

  if (!jobId) return <LinearProgress />;
  if (!container) return <LinearProgress />;
  if (!job) return <LinearProgress />;

  return <>{taskInterface}</>;
};
