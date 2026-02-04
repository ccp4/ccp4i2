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
from rest_framework.viewsets import ModelViewSet
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
        """
        Optionally filter files by job ID.

        Query parameters:
            job: Filter files to only those created by the specified job ID
        """
        queryset = models.File.objects.all()
        job_id = self.request.query_params.get("job")
        if job_id is not None:
            queryset = queryset.filter(job_id=job_id)
        return queryset

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def by_uuid(self, request, pk=None):
        try:
            the_file = models.File.objects.get(uuid=pk)
            serializer = serializers.FileSerializer(the_file, many=False)
            return Response(serializer.data)
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

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
    def download_by_uuid(self, request, pk=None):
        try:
            the_file = models.File.objects.get(uuid=pk)
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
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def digest_by_uuid(self, request, pk=None):
        """
        Get digest (summary) of file contents by UUID.

        Uses modern digest utility that:
        - Loads CDataFileContent subclass based on file type
        - Returns JSON dictionary representation
        - Tries to infer file type if not explicit

        Returns dict with file metadata like sequences, composition, cell, etc.
        """
        try:
            the_file = models.File.objects.get(uuid=pk)
            result = digest_file(the_file)

            # Check if digest returned an error
            if isinstance(result, dict) and result.get("status") == "Failed":
                return api_error(result.get("reason", "Digest failed"), status=400)

            return api_success(result)

        except models.File.DoesNotExist:
            logger.exception("File with UUID %s not found", pk)
            return api_error("File not found", status=404)
        except Exception as err:
            logger.exception("Failed to digest file with UUID %s", pk, exc_info=err)
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
