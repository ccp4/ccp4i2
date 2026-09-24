"""CPanddaDataset - one dataset as the PanDDA orchestrator's DATASETS list
declares it, and the orchestrator's run KPIs.

The orchestrator (``pandda_campaign``) takes a declared list, never a campaign
(design note, docs/pandda-campaign-design.md, section 3.1): validity stays
pure container logic, the task runs under i2run against a hand-built list,
and the list in params.xml is the record of exactly what was submitted. Each
item is one crystal: a label, the refined model and the reflections PanDDA
will analyse, and the ligand dictionary when there is one.

Resolvable by the def.xml class-name lookup via ``ccp4i2.core.CPanddaDataset``.
"""
from ccp4i2.core.base_object.class_metadata import content
from ccp4i2.core.CCP4Data import CData
from ccp4i2.core.CCP4PerformanceData import CPerformanceIndicator


class CPanddaDataset(CData):
    """One crystal for a PanDDA run: its label, model, reflections and
    (optionally) ligand dictionary."""

    class Meta:
        contents_order = ['DTAG', 'XYZIN', 'HKLIN', 'DICT', 'PROJECT_UUID', 'SOURCE_JOB_UUID']
        qualifiers = {"allowUndefined": False, "guiLabel": "Dataset"}

    DTAG = content(
        "CString", guiLabel='Label',
        toolTip="The name this dataset is known by (usually its project); "
                "written to Projects.csv. PanDDA itself sees a clean xtal-NNNN name")
    XYZIN = content(
        "CPdbDataFile", allowUndefined=False, mustExist=True, fromPreviousJob=True,
        guiLabel='Refined model',
        toolTip='The apo model PanDDA analyses, as dimple left it')
    HKLIN = content(
        "CMtzDataFile", allowUndefined=False, mustExist=True, fromPreviousJob=True,
        guiLabel='Reflections',
        toolTip="The MTZ dimple wrote: structure factors, map coefficients and the "
                "free-R set. The free-R label is normalised at staging")
    DICT = content(
        "CDictDataFile", allowUndefined=True, fromPreviousJob=True,
        guiLabel='Ligand dictionary',
        toolTip='Restraint dictionary for the ligand soaked into this crystal; '
                'bond orders are normalised at staging')
    # Ownership, carried explicitly: the files above are imported into the
    # *orchestrator's* project when the job runs and take its uuid, so the
    # dataset's own project -- where fan-out must put the receipt -- would be
    # lost. The manifest keys on this (design note 3.3), never on a name.
    PROJECT_UUID = content(
        "CString", allowUndefined=True, guiLabel='Project',
        toolTip="UUID of the project this dataset belongs to; where its receipt will go")
    SOURCE_JOB_UUID = content(
        "CString", allowUndefined=True, guiLabel='Source job',
        toolTip='UUID of the job whose outputs these files are (the dimple run)')


class CPanddaRunPerformance(CPerformanceIndicator):
    """What the run did, as KPIs."""

    class Meta:
        contents_order = ['nDatasets', 'nDatasetsProcessed', 'nDatasetsAnalysed', 'nEvents', 'wallSeconds']
        qualifiers = {"allowUndefined": True}

    nDatasets = content("CInt", guiLabel='Datasets staged')
    nDatasetsProcessed = content("CInt", guiLabel='Datasets processed',
                                 toolTip='processed_datasets/ directories PanDDA wrote')
    nDatasetsAnalysed = content("CInt", guiLabel='Datasets analysed',
                                toolTip='Datasets PanDDA characterised and searched (those with a Z-map)')
    nEvents = content("CInt", guiLabel='Events',
                      toolTip='Rows in the run-level events table')
    wallSeconds = content("CFloat", guiLabel='Wall time (s)')
