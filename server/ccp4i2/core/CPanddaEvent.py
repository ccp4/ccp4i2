"""CPanddaEvent - one PanDDA event, with its files as members.

The receipt task (``pandda_events``) declares ``outputData.EVENTS`` as a
``CList`` of this class, so each event carries its own event map and
candidate pose *inside* the event record. That is what makes the event<->file
association structural: it lives in ``params.xml`` and needs no manifest
(design note, docs/pandda-campaign-design.md, sections 7.1 and 10.1), and it is
what the glean walks to give every nested file its own ``File`` row
(``EVENTS[3].EVENT_MAP``).

Numbers are what PanDDA wrote, in PanDDA's units: ``OPTIMAL_CONTOUR`` is an
absolute map value, not sigma (section 7.3), and ``EVENT_IDX`` is an ordinal
within one run only (section 7.4). ``SITE_IDX`` is PanDDA's own site number
from its events table; the campaign-site reference of section 9 is v2 and is
added alongside it, not in its place.

Resolvable by the def.xml class-name lookup via ``ccp4i2.core.CPanddaEvent``,
like ``CDmDomain``. Core deliberately holds only the data; reading the PanDDA
tree is the wrapper's job.
"""
from ccp4i2.core.base_object.class_metadata import content
from ccp4i2.core.CCP4Data import CData
from ccp4i2.core import CCP4MathsData  # noqa: F401  (CXyz must be registered)
from ccp4i2.core.CCP4PerformanceData import CPerformanceIndicator


class CPanddaEvent(CData):
    """One event from one PanDDA run over one dataset: its scores, where it
    is, the BDC-corrected event map and the autobuilt candidate pose."""

    class Meta:
        contents_order = [
            'EVENT_IDX', 'SITE_IDX', 'BDC', 'SCORE', 'BUILD_SCORE', 'RSCC',
            'HIT_PROBABILITY', 'OPTIMAL_CONTOUR', 'DISPLAY_CONTOUR', 'CENTROID', 'LIGAND_ID',
            'EVENT_MAP', 'POSE', 'SCENE',
        ]
        qualifiers = {"allowUndefined": True}

    EVENT_IDX = content(
        "CInt", guiLabel='Event',
        toolTip="PanDDA's event number within this run; not stable across runs")
    SITE_IDX = content(
        "CInt", guiLabel='Site',
        toolTip="PanDDA's site number from its events table, if the table was written")
    BDC = content(
        "CFloat", guiLabel='BDC',
        toolTip='Background density correction factor for the event map')
    SCORE = content(
        "CFloat", guiLabel='Event score',
        toolTip="PanDDA's score for the event density")
    BUILD_SCORE = content(
        "CFloat", guiLabel='Build score',
        toolTip='Score of the selected autobuilt pose')
    RSCC = content(
        "CFloat", guiLabel='RSCC',
        toolTip='Real-space correlation of the selected pose with the event map')
    HIT_PROBABILITY = content(
        "CFloat", guiLabel='Hit probability',
        toolTip="PanDDA's hit-in-site probability from its events table, if written")
    OPTIMAL_CONTOUR = content(
        "CFloat", guiLabel='Optimal contour',
        toolTip='Contour level PanDDA chose for the pose, in absolute map units (not sigma)')
    DISPLAY_CONTOUR = content(
        "CFloat", guiLabel='Display contour',
        toolTip="Where to open the event map, in absolute map units: 1.5 times the spread of the "
                "map over its non-zero region, capped by the optimal contour. The optimal contour is "
                "where the build scored best and on a poorly characterised run can sit above the map's peak")
    CENTROID = content(
        "CXyz", guiLabel='Centroid',
        toolTip='Centroid of the event density, in orthogonal Angstroms')
    LIGAND_ID = content(
        "CString", guiLabel='Ligand',
        toolTip="The ligand's component code from the dictionary PanDDA used; the pose "
                "copy carries this name, where PanDDA itself writes LIG")
    EVENT_MAP = content(
        "CMapDataFile", guiLabel='Event map',
        toolTip='BDC-corrected event map for this event')
    POSE = content(
        "CPdbDataFile", guiLabel='Candidate pose',
        toolTip='Autobuilt ligand pose: a candidate to be judged, not the model of record')
    SCENE = content(
        "CMoorhenSceneDataFile", guiLabel='Scene',
        toolTip='A Moorhen scene of this event: apo model, Z-map at z=3, the event map at its '
                'display contour, the pose with its dictionary, centred on the event')

    def has_build(self):
        return self.BUILD_SCORE.isSet() or self.POSE.isSet()


class CPanddaReceiptPerformance(CPerformanceIndicator):
    """What a receipt expected against what it delivered, as KPIs, so a
    short receipt is visible in the job list without opening it."""

    class Meta:
        contents_order = ['nEventsExpected', 'nEventsDelivered',
                          'nPosesExpected', 'nPosesDelivered', 'bestBuildScore']
        qualifiers = {"allowUndefined": True}

    nEventsExpected = content("CInt", guiLabel='Events expected',
                              toolTip='Events PanDDA recorded for this dataset')
    nEventsDelivered = content("CInt", guiLabel='Event maps delivered',
                               toolTip='Event maps found and copied')
    nPosesExpected = content("CInt", guiLabel='Poses expected',
                             toolTip='Events PanDDA recorded a build for')
    nPosesDelivered = content("CInt", guiLabel='Poses delivered',
                              toolTip='Candidate poses found and copied')
    bestBuildScore = content("CFloat", guiLabel='Best build score',
                             toolTip='Highest build score over the delivered poses')
