#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

from ccp4i2.wrappers.xia2_dials.script import xia2_dials


class xia2_xds(xia2_dials.xia2_dials):
    """xia2 driving XDS.

    Identical to the DIALS task but for which pipeline runs, so it inherits
    everything and states only the differences.
    """

    TASKNAME = "xia2_xds"
    XML_ROOT_TAG = "Xia2Xds"

    # The DIALS scope goes; the XDS ones stay. Unlike xia2_dials this task
    # leaves xia2.settings.pipeline in the interface, because XDS has several
    # pipelines worth choosing between.
    PHIL_EXCLUDE_SCOPES = [
        "dials",
        "strategy",
        "xia2.settings.input.image",
    ]

    #: Set from the user's choice in get_command_target(), not fixed.
    PIPELINE_ARGS = []

    #: XDS variants only — the full upstream list includes the DIALS ones.
    XDS_PIPELINES = ["3d", "3dd", "3di", "3dii"]
    DEFAULT_PIPELINE = "3dii"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._restrict_pipeline_choice()

    def _restrict_pipeline_choice(self):
        """Offer only the XDS pipelines, defaulting to 3dii.

        xia2's master_phil lists the DIALS pipelines alongside the XDS ones
        and marks dials as the default. Offering those here would let the user
        select DIALS from the XDS task.
        """
        try:
            settings = self.container.controlParameters.xia2.xia2__settings
            pipeline = settings.xia2__settings__pipeline
        except AttributeError:
            return
        pipeline.set_qualifier("enumerators", list(self.XDS_PIPELINES))
        pipeline.set_qualifier("onlyEnumerators", True)
        pipeline.set_qualifier("default", self.DEFAULT_PIPELINE)
        # xia2's default is dials, which is not one of the choices offered
        # here — leaving it would both mis-run the task and fail validation
        # against the narrowed enumerators.
        if str(pipeline) not in self.XDS_PIPELINES:
            pipeline.set(self.DEFAULT_PIPELINE)

    def get_command_target(self):
        """The pipeline the user chose, named after the phil file so it wins.

        The phil the base class writes carries xia2's own default of dials, so
        the choice has to be asserted on the command line rather than left to
        the file.
        """
        try:
            pipeline = str(
                self.container.controlParameters.xia2.xia2__settings
                .xia2__settings__pipeline
            )
        except AttributeError:
            pipeline = self.DEFAULT_PIPELINE
        if pipeline not in self.XDS_PIPELINES:
            pipeline = self.DEFAULT_PIPELINE
        return [f"pipeline={pipeline}"]

    def _collect_pickles_and_jsons(self):
        # XDS produces no DIALS pickles or JSONs to collect.
        pass

    @staticmethod
    def _get_annotation(prefix, suffix):
        """Form suitable annotation strings"""
        return prefix + " from XDS integration of " + suffix
