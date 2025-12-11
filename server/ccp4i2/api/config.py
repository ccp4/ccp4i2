import sys
import pathlib
from django.apps import AppConfig
# DISABLED: Old ccp4i2 imports
# from ccp4i2.googlecode import diff_match_patch_py3
# sys.path.append(str(pathlib.Path(diff_match_patch_py3.__file__).parent.parent))


class ApiConfig(AppConfig):
    name = "ccp4i2.api"
