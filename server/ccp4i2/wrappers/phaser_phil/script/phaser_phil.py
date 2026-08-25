from ccp4i2.core.PhilPluginScript import PhilPluginScript


class phaser_phil(PhilPluginScript):
    """Phaser, driven through its PHIL interface.

    The parameters come from the phaser installation's own
    phenix_interface/__init__.params at runtime, rather than from a def.xml
    generated from a snapshot of it.

    This also settles the FIXME the previous implementation carried: it wrote
    a file of bare definition=value lines and noted that fetching them against
    Phaser's master phil would be better. build_working_phil() does exactly
    that, so the job gets a validated, fully-formed phil.
    """

    TASKNAME = "phaser_phil"
    TASKCOMMAND = "ccp4-python"

    #: Phaser ships its parameters as a file inside the package.
    PHIL_PARAMS_FILE = "phaser:phenix_interface/__init__.params"

    def get_command_target(self):
        """Phaser's entry point, which ccp4-python runs with the phil file."""
        from phaser.command_line import main

        return [main.__file__]
