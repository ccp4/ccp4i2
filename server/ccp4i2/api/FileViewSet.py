import logging
from django.http import FileResponse, JsonResponse
from rest_framework.response import Response
from rest_framework import status
from rest_framework.parsers import MultiPartParser, JSONParser
from rest_framework.permissions import IsAuthenticated

# Modern utilities
from ..lib.utils.files.digest import digest_file
from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder

# Modern utilities
from ..lib.utils.files.preview import preview_file
from xml.etree import ElementTree as ET
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.viewsets import ModelViewSet, ReadOnlyModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from . import serializers
from ..db import models
from ..lib.response import api_success, api_error

logger = logging.getLogger(f"ccp4i2:{__name__}")


class FileViewSet(ModelViewSet):
    queryset = models.File.objects.all()
    serializer_class = serializers.FileSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]
    filterset_fields = ["job"]

    def get_queryset(self):
        queryset = models.File.objects.all()
        job_id = self.request.query_params.get("job")
        if job_id is not None:
            queryset = queryset.filter(job_id=job_id)
        file_type = self.request.query_params.get("type")
        if file_type is not None:
            queryset = queryset.filter(type=file_type)
        directory = self.request.query_params.get("directory")
        if directory is not None:
            queryset = queryset.filter(directory=directory)
        # `?project={id}` filters files transitively through Job→Project.
        # Replaces the deleted ProjectViewSet.files @action; the legacy
        # `?job__project={id}` form is retained as an alias for back-compat.
        project_id = self.request.query_params.get("project")
        if project_id is not None:
            queryset = queryset.filter(job__project_id=project_id)
        job_project = self.request.query_params.get("job__project")
        if job_project is not None:
            queryset = queryset.filter(job__project_id=job_project)
        parent_isnull = self.request.query_params.get("job__parent__isnull")
        if parent_isnull is not None:
            queryset = queryset.filter(job__parent__isnull=parent_isnull.lower() == "true")
        return queryset

    # by_uuid() / download_by_uuid() / digest_by_uuid() @actions removed —
    # consumers should hit /files_by_uuid/{uuid}/ (FileByUuidViewSet) instead.
    # The old actions abused the {pk} URL slot to carry a uuid; the new
    # viewset uses lookup_field="uuid" so the URL pattern is honest.

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def download(self, request, pk=None):
        try:
            the_file = models.File.objects.get(id=pk)
            return FileResponse(open(the_file.path, "rb"), filename=the_file.name)
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def digest(self, request, pk=None):
        """
        Get digest (summary) of file contents.

        Uses modern digest utility that:
        - Loads CDataFileContent subclass based on file type
        - Returns JSON dictionary representation
        - Tries to infer file type if not explicit

        Returns dict with file metadata like sequences, composition, cell, etc.
        """
        try:
            the_file = models.File.objects.get(id=pk)
            result = digest_file(the_file)

            # Check if digest returned an error
            if isinstance(result, dict) and result.get("status") == "Failed":
                return api_error(result.get("reason", "Digest failed"), status=400)

            return api_success(result)

        except models.File.DoesNotExist:
            logger.exception("File %s not found", pk)
            return api_error("File not found", status=404)
        except Exception as err:
            logger.exception("Failed to digest file %s", pk, exc_info=err)
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.FileSerializer,
    )
    def preview(self, request, pk=None):
        try:
            the_file = models.File.objects.get(id=pk)
            the_viewer = request.data.get("viewer")
            preview_file(the_viewer, str(the_file.path))
            return api_success({"previewed": True})
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except ValueError as err:
            logger.warning("Unsupported viewer requested: %s", err)
            return api_error(str(err), status=400)
        except Exception as err:
            logger.exception("Failed to preview file %s", pk, exc_info=err)
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def file_path(self, request, pk=None):
        """Return the absolute filesystem path for a file (Electron desktop only)."""
        try:
            the_file = models.File.objects.get(id=pk)
            return api_success({"path": str(the_file.path)})
        except models.File.DoesNotExist:
            return api_error("File not found", status=404)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def molblock(self, request, pk=None):
        """
        Convert a ligand dictionary CIF file to MolBlock format.

        This endpoint takes a refmac dictionary CIF file (type application/refmac-dictionary)
        and converts it to MolBlock format suitable for 2D rendering with RDKit.

        Returns:
            JSON with 'molblock' key containing the MolBlock string,
            'ligand_code' with the 3-letter ligand code from the CIF,
            or error if conversion fails.
        """
        try:
            the_file = models.File.objects.get(id=pk)

            # Verify file type is a ligand dictionary
            if the_file.type.name != "application/refmac-dictionary":
                return api_error(
                    f"File type must be application/refmac-dictionary, got {the_file.type.name}",
                    status=400
                )

            # Check file exists
            if not the_file.path.exists():
                return api_error("File not found on disk", status=404)

            file_path = str(the_file.path)

            # Extract ligand code from CIF using gemmi
            ligand_code = None
            try:
                import gemmi
                doc = gemmi.cif.read_file(file_path)
                for block in doc:
                    # Try to get comp_id from atom loop
                    for comp_id in block.find_loop("_chem_comp_atom.comp_id"):
                        ligand_code = comp_id
                        break
                    if ligand_code:
                        break
                    # Fallback: try block name (often the ligand code)
                    if block.name and len(block.name) <= 4:
                        ligand_code = block.name
                        break
            except Exception as e:
                logger.warning("Could not extract ligand code from CIF: %s", e)

            # Import the conversion function
            from ..wrappers.acedrgNew.script.cifToMolBlock import cifFileToMolBlock

            # Convert CIF to MolBlock
            molblock = cifFileToMolBlock(file_path)

            if not molblock:
                return api_error("Failed to convert CIF to MolBlock - check server logs for details", status=400)

            return api_success({
                "molblock": molblock,
                "ligand_code": ligand_code,
            })

        except models.File.DoesNotExist:
            logger.exception("File %s not found", pk)
            return api_error("File not found", status=404)
        except ImportError as err:
            logger.exception("Failed to import cifToMolBlock: %s", err)
            return api_error("MolBlock conversion not available (missing RDKit/gemmi)", status=500)
        except Exception as err:
            logger.exception("Failed to convert file %s to molblock", pk, exc_info=err)
            return api_error(str(err), status=500)


class FileByUuidViewSet(ReadOnlyModelViewSet):
    """
    Read-only view of files addressed by UUID rather than Django PK.

    CCP4i2 parameter files (def.xml, fileUse refs, exported job bundles)
    persist file references as UUIDs — they have no knowledge of the
    deployment-local Django integer ID. This viewset gives consumers a
    canonical way to resolve those refs:

        GET /files_by_uuid/{uuid}/           — file metadata
        GET /files_by_uuid/{uuid}/download/  — file content
        GET /files_by_uuid/{uuid}/digest/    — file content digest

    Mutations remain on the id-addressed FileViewSet; uuid is a lookup
    surface for read paths only.
    """
    queryset = models.File.objects.all()
    serializer_class = serializers.FileSerializer
    lookup_field = "uuid"
    permission_classes = [IsAuthenticated]

    @action(detail=True, methods=["get"])
    def download(self, request, uuid=None):
        try:
            the_file = self.get_object()
            return FileResponse(open(the_file.path, "rb"), filename=the_file.name)
        except models.File.DoesNotExist as err:
            logging.exception("File with uuid %s not found", uuid, exc_info=err)
            return api_error(str(err), status=404)

    @action(detail=True, methods=["get"])
    def digest(self, request, uuid=None):
        """
        Get digest (summary) of file contents.

        Uses modern digest utility that:
        - Loads CDataFileContent subclass based on file type
        - Returns JSON dictionary representation
        - Tries to infer file type if not explicit

        Returns dict with file metadata like sequences, composition, cell, etc.
        """
        try:
            the_file = self.get_object()
            result = digest_file(the_file)
            if isinstance(result, dict) and result.get("status") == "Failed":
                return api_error(result.get("reason", "Digest failed"), status=400)
            return api_success(result)
        except models.File.DoesNotExist:
            logger.exception("File with uuid %s not found", uuid)
            return api_error("File not found", status=404)
        except Exception as err:
            logger.exception("Failed to digest file with uuid %s", uuid, exc_info=err)
            return api_error(str(err), status=500)
