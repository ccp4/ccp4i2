from __future__ import print_function
import traceback
import json
import logging
import traceback

from ccp4i2.core import CCP4File
from ccp4i2.core import CCP4Data
from ccp4i2.core import CCP4Container
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.db import models
from ..plugins.get_plugin import get_job_plugin
from ..containers.json_encoder import CCP4i2JsonEncoder


logger = logging.getLogger(f"ccp4i2:{__name__}")


def get_job_container(job):
    # Placeholder implementation
    return CCP4Container.CContainer()


def json_for_job_container(job: models.Job):
    plugin = get_job_plugin(job)
    container = plugin.container
    return json.dumps(container, cls=CCP4i2JsonEncoder)
