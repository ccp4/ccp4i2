
from ccp4i2.baselayer import QtCore
from qtgui import CCP4TaskWidget
import os, sys, re, subprocess
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

class CTaskArpWarpClassic(CCP4TaskWidget.CTaskWidget):

    TASKNAME = 'arp_warp_classic'
    TASKVERSION = 0.0
    TASKMODULE = 'test' if sys.platform.startswith('win') else 'model_building'
    TASKTITLE='ARP/wARP'
    DESCRIPTION='Build model (ARP/wARP classic)'
    WHATNEXT = ['coot_rebuild', 'prosmart_refmac', 'edstats']
    MGDISPLAYFILES = ['FPHIOUT']

    def drawContents(self):
        self.labels = dict()
        notoggle = True
        notoggle = False

        folder = self.openFolder(folderFunction='inputData', title='Input Data')

        self.createLine([
            'label', 'Run ARP/wARP for',
            'widget', 'AWA_ARP_MODE'])

        self.createLine([
            'subtitle', 'Select experimental data',
            'Selection of anomalous data provies additional option in refinement, see below'])

        self.openSubFrame(frame=[True])

        self.createLine([
            'message', 'XXXXX Enter input MTZ file name (HKLIN)',
            'widget','AWA_FOBS'])

        self.createLine([
            'message', 'XXXXX Enter input MTZ file name (HKLIN)',
            'widget', 'AWA_FREE'])

        self.closeSubFrame()

        self.createLine([
            'subtitle', 'Provide model or phases for calculation of starting maps',
            'Select model after molecular replacement or phases from density modification'])

        self.openSubFrame(frame=[True])

        self.createLine([
            'message', 'XXXXX Enter input MTZ file name (HKLIN)',
            'widget', 'AWA_PHINI'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['WARPNTRACEPHASES']])

        self.createLine([
            'message', 'XXXXX Enter name of file containing molecule coordinates (MODELIN)',
            'widget', 'AWA_MODELIN'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['WARPNTRACEMODEL', 'MOLREP']])

        self.closeSubFrame()

        tooltip = 'To build polyalanine model unselect sequence and define the numer of residues and molecules in the AU.'
        tooltip += ' Sequence is compulsory if you want to use NCS model extension and/or NCS restraints.'
        self.createLine([ 'subtitle', 'Define contents of the asymmetric unit', tooltip])

        self.openSubFrame(frame=[True])

        self.createLine([
            'message', 'Enter name of file containing the sequence (SEQIN)',
            'widget','AWA_SEQIN'])

        self.closeSubFrame()

        tooltip = '(1) Use of anomalous data is only active if the experimental data contain anomalous pairs.'
        tooltip += ' (2) The phases used in refinement can be the same as used for calculation of starting maps.'
        self.createLine(['subtitle', 'Define phase information for refinements', tooltip])

        self.openSubFrame(frame=[True])

        self.createLine([
            'widget', '-guiMode', 'radio', 'AWA_REF_MODE'])

        self.createLine([
            'message', 'XXXXX Enter input MTZ file name (HKLIN)',
            'widget', 'AWA_PHREF'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'open', ['AWA_HL', 'AWA_PHASED']])

        self.createLine([
            'message', 'XXXXX Enter name of file containing heavy atom coordinates (MODELIN)',
            'widget', 'AWA_HEAVYIN'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'open', ['AWA_SAD']])

        line = self.createLine([
            'label', 'Scattering',
            'widget', 'AWA_ANO_OPTION'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'open', ['AWA_SAD']])

        self.createLine([
            'widget', 'AWA_SCAT_LAMBDA_AWA_SAD',
            'label', 'A.'],
                appendLine = line,
                toggle = [] if notoggle else ['AWA_ANO_OPTION', 'open', ['LAMBDA']])

        self.createLine([
            'label', 'for atom',
            'widget', 'AWA_SCAT_ATOM',
            'label', "f' =",
            'widget', 'AWA_SCAT_FP_AWA_SAD',
            'label', "f'' =",
            'widget', 'AWA_SCAT_FDP_AWA_SAD'],
                appendLine = line,
                toggle = [] if notoggle else ['AWA_ANO_OPTION', 'open', ['MEASURED']])

        self.closeSubFrame()

        self.closeFolder()

        folder = self.openFolder(folderFunction='controlParameters', title='Parameters')

        self.createLine([
            'subtitle', 'ARP/wARP flow parameters'])

        line = self.createLine([
            'message', 'Each ARP-REFMAC cycle consists of refinement with Refmac and remodelling dummy-atom substructure with ARP',
            'label', 'Do',
            'widget', 'AWA_BIG_CYCLES',
            'label', 'macrocycles, each including',
            'widget', 'AWA_SMALL_CYCLES',
            'label', '.'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'hide', ['MOLREP']])

        self.label_init('AWA_TOTAL_CYCLES', line, 4, 'ARP-REFMAC cycles')

        self.createLine([
            'label', 'Do',
            'widget', 'AWA_SMALL_CYCLES',
            'label', 'ARP-REFMAC cycles (NB: no autotracing in this mode)'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['MOLREP']])

        self.openSubFrame(
            toggle = [] if notoggle else ['AWA_ARP_MODE', 'hide', ['MOLREP']])

        self.createLine([
            'widget', 'AWA_USE_COND',
            'label', 'Use Conditional Restraints for free atoms.'])

        self.createLine([
            'widget', 'AWA_FORCE_COND',
            'label', 'Impose the use of Conditional Restraints for very large structures, time consuming.'])

        self.createLine([
            'widget', 'AWA_FAKE_DATA',
            'label', 'Use Fake Data.'])

        self.createLine([
            'widget', 'AWA_NCS_RESTRAINTS',
            'label', 'Use Non-Crystallographic Symmetry Restraints.'])

        self.createLine([
            'widget', 'AWA_NCS_EXTENSION',
            'label', 'Use Non-Crystallographic Symmetry information to extend chains.'])

        self.createLine([
            'widget', 'AWA_LOOPS',
            'label', 'Use Loopy to build loops between chain fragments at the end of the run.'])

        self.createLine([
            'widget', 'AWA_BUILD_SIDE',
            'label', 'Start docking the autotraced chains to sequence after',
            'widget', 'AWA_SIDE_AFTER',
            'label', 'building cycles.'],
                toggleFunction = [self.AWA_SEQIN_isset, ['AWA_SEQIN']])

        self.createLine([
            'widget', 'AWA_IS_SEMET',
            'label', 'Sequence dock as a Se-Methionine substituted protein.'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'open', ['AWA_SAD']])

        self.createLine([
            'widget', 'AWA_ALBE',
            'label', 'Search for helices and strands before each building cycle.'])

        self.createLine([
            'widget', 'AWA_SKIP_BUILD',
            'label', 'Skip the chain tracing for the first',
            'widget', 'AWA_SKIP_CYCLES',
            'label', 'building cycles.'])

        self.closeSubFrame()

        self.createLine([
            'widget', 'AWA_FREEBUILD',
            'label', 'Before the start of autobuilding construct new free atoms model in map'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['WARPNTRACEMODEL']])

        self.createLine([
            'widget', 'AWA_FLATTEN',
            'label', 'Before the start of autobuilding perform density modification with DM'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['WARPNTRACEMODEL']])

        self.createLine([
            'label', 'Iterate the autotracing up to',
            'widget', 'AWA_MULTITRACE',
            'label', 'times.'],
                toggle = [] if notoggle else ['AWA_ARP_MODE', 'open', ['WARPNTRACEPHASES', 'WARPNTRACEMODEL']])

        self.createLine([
            'label', 'Add atoms in density above',
            'widget', 'AWA_ADDATOM_SIGMA',
            'label', 'sigma, and remove atoms in density below',
            'widget', 'AWA_REMATOM_SIGMA',
            'label', 'sigma.'])

        self.createLine([
            'label', 'Add and remove',
            'widget', 'AWA_UP_UPDATE',
            'label', 'times atoms than calculated automatically.'])

        self.createLine([
            'subtitle', 'Refmac parameters'])

        self.createLine([
            'widget', 'AWA_TWIN',
            'label', 'Attempt to correct for data collected from a twinned crystal.'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'hide', ['AWA_SAD']])

        self.createLine([
            'widget', 'AWA_NCYCLES',
            'label', 'cycles of refinement with REFMAC per ARP-REFMAC cycle.'])

        self.createLine([
            'label', 'Use phase restraints with blurring factor of',
            'widget', 'AWA_PHASE_BLUR'],
                toggle = [] if notoggle else ['AWA_REF_MODE', 'hide', ['AWA_REFML', 'AWA_SAD']])

        self.createLine([
            'widget', 'AWA_WEIGHT_MODE',
            'label', 'matrix weighting for Xray / Geometry'])

        self.createLine([
            'label', 'Matrix weight for Xray / Geometry',
            'widget', 'AWA_WMAT'],
                toggle = [] if notoggle else ['AWA_WEIGHT_MODE', 'hide', ['AUTO']])

        self.createLine([
            'widget', 'AWA_RIDGE_RESTRAINTS',
            'label', 'Jelly Body restraints'])

        self.createLine([
            'label', 'Use for scaling the',
            'widget', 'AWA_SCALE',
            'label', 'scaling model, with an',
            'widget', 'AWA_SCANIS',
            'label', 'scaling B factor.'])

        self.createLine([
            'label', 'Scaling and sigmaa calculation will be done with the',
            'widget', 'AWA_REFMAC_REF_SET',
            'label', 'set of reflections'],
                toggleFunction = [self.AWA_FREE_isset, ['AWA_FREE']])

        self.createLine([
            'widget', 'AWA_SOLVENT',
            'label', 'Solvent Mask correction'])

        self.closeFolder()

        if self.isEditable():
            self.handle_AWA_REF_MODE_changed()
            self.container.controlParameters.AWA_REF_MODE.dataChanged.connect( self.handle_AWA_REF_MODE_changed)
            self.handle_AWA_ARP_MODE_changed()
            self.container.controlParameters.AWA_ARP_MODE.dataChanged.connect( self.handle_AWA_ARP_MODE_changed)

            self.handle_AWA_FOBS_changed()
            self.container.inputData.AWA_FOBS.dataChanged.connect( self.handle_AWA_FOBS_changed)

            self.handle_AWA_TOTAL_CYCLES_changed()
            self.handle_AWA_TOTAL_CYCLES_changed()
            self.container.controlParameters.AWA_BIG_CYCLES.dataChanged.connect( self.handle_AWA_TOTAL_CYCLES_changed)
            self.container.controlParameters.AWA_SMALL_CYCLES.dataChanged.connect( self.handle_AWA_TOTAL_CYCLES_changed)

        b1 = bool(self.container.controlParameters.AWA_MOCKYES.get())
        if b1 or 'AWA_DEVEL' in os.environ:
            folder = self.openFolder(folderFunction='controlParameters', title='Developer')
            tooltip = 'There is a bug if you see this section'
            self.createLine(['subtitle', 'Developer options', tooltip])
            self.openSubFrame(frame=[True])

            self.createLine([
                'widget', 'AWA_MOCKYES',
                'label', 'mock run using precalculated results from PSP example'])

            self.createLine([
                'widget', 'AWA_JSRVIEW',
                'label', 'show report in jsrview'],
                    toggle = [] if notoggle else ['AWA_MOCKYES', 'open', [True]])

            self.createLine([
                'label', 'Pause for',
                'widget', 'AWA_MOCKPAUSE',
                'label', 'millisecond after each line read from the log file'],
                    toggle = [] if notoggle else ['AWA_MOCKYES', 'open', [True]])

            self.closeSubFrame()
            self.closeFolder()

    @QtCore.Slot()
    def handle_AWA_TOTAL_CYCLES_changed(self):
        small = self.container.controlParameters.AWA_SMALL_CYCLES.get()
        big = self.container.controlParameters.AWA_BIG_CYCLES.get()
        text = 'ARP-REFMAC cycles'
        if type(small) is int and type(big) is int:
            text = 'ARP-REFMAC cycles (%d cycles in total)' %(small* big)

        self.label_set('AWA_TOTAL_CYCLES', text)

    def label_init(self, varname, line, index, default):
        container = getattr(self.container.guiParameters, varname)
        widget = line.layout().itemAt(index).widget()
        setattr(self, 'label_' + varname, (container, widget))
        text = default
        if container.isSet():
            value = str(container.get())
            if value:
                text = value

        self.label_set(varname, text)

    def label_set(self, varname, text):
        container, widget = getattr(self, 'label_' + varname)
        container.set(text)
        widget.setText(text)

    @QtCore.Slot()
    def handle_AWA_ARP_MODE_changed(self):
        use_xyz = self.container.controlParameters.AWA_ARP_MODE == 'WARPNTRACEMODEL'
        self.update_requirements(('AWA_MODELIN', 'AWA_PHINI'), ('AWA_MODELIN',) if use_xyz else ('AWA_PHINI',))

    @QtCore.Slot()
    def handle_AWA_REF_MODE_changed(self):
        use_ano = self.container.controlParameters.AWA_REF_MODE == 'AWA_SAD'
        use_none = self.container.controlParameters.AWA_REF_MODE == 'AWA_REFML'
        self.update_requirements(('AWA_HEAVYIN', 'AWA_PHREF'), ('AWA_HEAVYIN',) if use_ano else ('AWA_PHREF',) if not use_none else ())

        qlt = 'As anomalous Fs' if use_ano else 'As mean Fs'
        self.container.inputData.AWA_FOBS.setQualifier('guiLabel', qlt)
        self.getWidget('AWA_FOBS').layout().itemAt(1).widget().setText(qlt)

    def update_requirements(self, list_all, list_required):
        for key in list_all:
            qf = self.getWidget(key)
            qcb = getattr(qf, 'jobCombo', None)
            if qcb:
                cdf_is_required = key in list_required
                cdf_was_required = qf.model.qualifiers('mustExist')
                if cdf_is_required is not cdf_was_required:
                    cdf = getattr(self.container.inputData, key)
                    cdf.setQualifier('allowUndefined', not cdf_is_required)
                    cdf.setQualifier('mustExist', cdf_is_required)
                    qcb.setItemText(qcb.findData(-1), '.. must be selected' if cdf_is_required else '.. is not used')
                    qf.validate()

    @QtCore.Slot()
    def handle_AWA_FOBS_changed(self):
        type_dict = getattr(self, 'type_dict', dict(GLGL=1, KMKM=2, JQ=3, FQ=4))
        obs_cdf = self.container.inputData.AWA_FOBS
        obs_cdf.loadFile()

        old_type = getattr(self, 'type', None)
        self.type = type_dict.get(obs_cdf.fileContent.columnSignature(), None)
        type_changed = self.type != old_type

        old_wavelength = getattr(self, 'wavelength', None)
        wavelength_list = obs_cdf.fileContent._value['wavelengths']
        self.wavelength = wavelength_list[-1] if wavelength_list else None
        wavelength_changed = self.wavelength != old_wavelength

        if type_changed:
            obs_bg_method = getattr(self.getWidget('AWA_REF_MODE').widget, 'buttons', None)
            if obs_bg_method is not None:
                for obs_button in self.getWidget('AWA_REF_MODE').widget.buttons():
                    value = str(obs_button.objectName())
                    if value == 'AWA_REFML':
                        obs_button_none = obs_button

                    elif value == 'AWA_SAD':
                        obs_button_ano = obs_button

                obs_was_ano = obs_button_ano.isEnabled()
                obs_is_ano = self.type in (None, 1, 2)

                if not obs_was_ano and obs_is_ano:
                    obs_button_ano.setEnabled(True)

                elif obs_was_ano and not obs_is_ano:
                    if obs_button_ano.isChecked():
                        obs_button_none.click()

                    obs_button_ano.setEnabled(False)

                #self.container.controlParameters.AWA_SCAT_LAMBDA_AWA_SAD.setQualifier('allowUndefined', not obs_is_ano)
                self.getWidget('AWA_SCAT_LAMBDA_AWA_SAD').widget.editingFinished.emit()

        if wavelength_changed:
            wl_widget = self.getWidget('AWA_SCAT_LAMBDA_AWA_SAD').widget
            wl_widget.setText(str(round(self.wavelength, 6)) if self.wavelength != None else '1.0')
            wl_widget.editingFinished.emit()

    def AWA_SEQIN_isset(self):
        return self.container.inputData.AWA_SEQIN.isSet()

    def AWA_FREE_isset(self):
        return self.container.inputData.AWA_FREE.isSet()

