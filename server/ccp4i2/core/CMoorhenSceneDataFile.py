"""CMoorhenSceneDataFile - a Moorhen scene a job authored about its outputs.

A scene is a portable YAML document (client/renderer/lib/scene, mirrored for
the server in ``ccp4i2/scene_contracts``): which files to load, how to draw
them, which maps to contour and where, and where the camera sits. A task
that knows what its outputs should look like -- a PanDDA receipt knows the
event centroid, the map's contour level and which pose goes with which
dictionary -- writes one as an output, and the Moorhen preview of the job
applies it instead of guessing from raw files.

Resolvable by the def.xml class-name lookup via
``ccp4i2.core.CMoorhenSceneDataFile``.
"""
from ccp4i2.core.base_object.cdata_file import CDataFile


class CMoorhenSceneDataFile(CDataFile):
    class Meta:
        qualifiers = {
            "mimeTypeName": "application/moorhen-scene",
            "mimeTypeDescription": "Moorhen scene",
            "fileExtensions": ["scene.yaml", "yaml", "yml"],
            "guiLabel": "Moorhen scene",
            "toolTip": "A scene the job authored: what to load, how to draw it, where to look",
        }
