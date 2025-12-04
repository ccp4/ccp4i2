import sys
import pathlib
from django.core.management.base import BaseCommand
# DISABLED: Old ccp4i2 imports
# from ccp4i2.googlecode import diff_match_patch_py3
# sys.path.append(str(pathlib.Path(diff_match_patch_py3.__file__).parent.parent))
# from core import CCP4DataManager
# from core.CCP4DataManager import CDataManager


class Command(BaseCommand):

    help = "Generate start point for template bindings"
    requires_system_checks = []

    def handle(self, **kwargs):
        data_manager: CDataManager = CCP4DataManager.CDataManager()
        data_manager.buildClassLookup()
        skip_classes = [
            "CFloat",
            "CInt",
            "CBoolean",
            "CString",
            "CDict",
            "CList",
            "CCollection",
        ]
        skip_fields = [
            "project",
            "baseName",
            "relPath",
            "annotation",
            "dbFileId",
            "subType",
            "contentFlag",
        ]

        for key, the_class in sorted(
            data_manager.clsLookup.items(), key=lambda a: a[0]
        ):
            the_class = data_manager.getClass(key)
            the_instance = the_class()
            parent_classes = [cls.__name__ for cls in the_class.__mro__]
            parent_class = the_class.__bases__[-1].__name__

            if the_class.__name__ in skip_classes:
                continue

            if parent_class in skip_classes and parent_class != "CList":
                print(f"export type {key} = {parent_class};")
                continue

            if parent_class == "CDataFile":
                extra_fields = [
                    f for f in the_instance.CONTENTS if f not in skip_fields
                ]
                if not extra_fields:
                    print(f"type {key} = CDataFile;")
                    continue

            if parent_class == "CList":
                print(
                    f"export type {key} = {the_instance.SUBITEM['class'].__name__}[];"
                )
                continue

            print(
                f"export interface {key} {f'extends {parent_class}' if parent_class != 'CData' else ''} {{"
            )
            for content_name, content_info in the_instance.CONTENTS.items():
                if content_info["class"].__name__ == "CList":
                    sub_item = content_info.get(
                        "subItem", {"class": {"__name__": any}}
                    ).get("class", "any")
                    try:
                        sub_item_class_name = sub_item.__name__
                    except Exception as err:
                        print(err)
                        sub_item_class_name = "any"
                    print(
                        f"\t{content_name}: {sub_item_class_name}[]"
                        # f"\t{content_name}: {the_instance['subItem'].__class__.__name__}[];"
                    )
                elif (key in "CDataFile") or not (
                    ("CDataFile" in parent_classes) and content_name in skip_fields
                ):
                    print(f"\t{content_name}: {content_info['class'].__name__};")
            print("};")
