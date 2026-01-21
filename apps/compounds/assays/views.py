"""
Assays ViewSets

DRF ViewSets for assay models with reversion support.

Permission model:
- Read access: Anyone (including 'user' operating level)
- Write access: Requires 'contributor' or 'admin' operating level
"""

import reversion
from django_filters.rest_framework import DjangoFilterBackend
from rest_framework import viewsets, filters, status
from rest_framework.decorators import action
from rest_framework.response import Response
from reversion.models import Version

from users.permissions import IsContributorOrReadOnly
from .analysis import analyse_assay, analyse_data_series

from .models import (
    DilutionSeries,
    Protocol,
    ProtocolDocument,
    Assay,
    DataSeries,
    AnalysisResult,
    Hypothesis,
    FittingMethod,
)
from .serializers import (
    FittingMethodSerializer,
    DilutionSeriesSerializer,
    ProtocolSerializer,
    ProtocolDetailSerializer,
    ProtocolDocumentSerializer,
    ProtocolDocumentCreateSerializer,
    AssayListSerializer,
    AssayDetailSerializer,
    AssayCreateSerializer,
    DataSeriesListSerializer,
    DataSeriesDetailSerializer,
    DataSeriesCreateSerializer,
    AnalysisResultSerializer,
    HypothesisListSerializer,
    HypothesisDetailSerializer,
    TableOfValuesImportSerializer,
)


class ReversionMixin:
    """Mixin to add reversion support to ViewSets."""

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Created via API")
            return instance

    def perform_update(self, serializer):
        with reversion.create_revision():
            instance = serializer.save()
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Updated via API")
            return instance

    def perform_destroy(self, instance):
        with reversion.create_revision():
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment(f"Deleted via API")
            instance.delete()

    @action(detail=True, methods=['get'])
    def history(self, request, pk=None):
        """Get version history for this object."""
        obj = self.get_object()
        versions = Version.objects.get_for_object(obj)
        return Response([{
            'id': v.id,
            'date': v.revision.date_created,
            'user': v.revision.user.email if v.revision.user else None,
            'comment': v.revision.comment,
        } for v in versions[:50]])


class FittingMethodViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Fitting Methods."""
    queryset = FittingMethod.objects.filter(is_active=True)
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.SearchFilter, filters.OrderingFilter]
    search_fields = ['name', 'description']
    ordering = ['name']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            from .serializers import FittingMethodDetailSerializer
            return FittingMethodDetailSerializer
        return FittingMethodSerializer


class DilutionSeriesViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Dilution Series."""
    queryset = DilutionSeries.objects.all()
    serializer_class = DilutionSeriesSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [filters.OrderingFilter]
    ordering = ['unit']


class ProtocolViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Protocols."""
    queryset = Protocol.objects.select_related('preferred_dilutions', 'created_by')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['analysis_method']
    search_fields = ['name', 'comments']
    ordering_fields = ['name', 'created_at']
    ordering = ['name']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            return ProtocolDetailSerializer
        return ProtocolSerializer

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save(created_by=self.request.user if self.request.user.is_authenticated else None)
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment("Created via API")
            return instance

    @action(detail=True, methods=['get'])
    def assays(self, request, pk=None):
        """List assays using this protocol."""
        protocol = self.get_object()
        assays = protocol.assays.all()
        serializer = AssayListSerializer(assays, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['post'])
    def propagate_target(self, request, pk=None):
        """
        Update all assays using this protocol to use the protocol's target.

        This is useful when the protocol's default target is changed and
        the user wants to apply that change to existing assays.

        Returns the count of updated assays.
        """
        protocol = self.get_object()

        if not protocol.target:
            return Response({
                'status': 'error',
                'error': 'Protocol has no target set'
            }, status=status.HTTP_400_BAD_REQUEST)

        # Count assays that would be updated (different target or null)
        assays_to_update = protocol.assays.exclude(target=protocol.target)
        count = assays_to_update.count()

        if count == 0:
            return Response({
                'status': 'success',
                'updated': 0,
                'message': 'All assays already have this target'
            })

        # Update all assays to use the protocol's target
        with reversion.create_revision():
            assays_to_update.update(target=protocol.target)
            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment(f"Propagated target from protocol: {protocol.name}")

        return Response({
            'status': 'success',
            'updated': count,
            'target_id': str(protocol.target.id),
            'target_name': protocol.target.name
        })


class AssayViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Assays."""
    queryset = Assay.objects.select_related('protocol', 'target', 'created_by')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['protocol', 'target']
    ordering_fields = ['created_at']
    ordering = ['-created_at']

    def perform_destroy(self, instance):
        """Delete assay and clean up associated analysis results."""
        # Delete analysis results first (they're not CASCADE from DataSeries)
        analysis_ids = list(
            instance.data_series.filter(analysis__isnull=False)
            .values_list('analysis_id', flat=True)
        )
        if analysis_ids:
            AnalysisResult.objects.filter(id__in=analysis_ids).delete()

        # Call parent to delete assay (cascades to data series)
        super().perform_destroy(instance)

    def get_serializer_class(self):
        if self.action == 'create':
            return AssayCreateSerializer
        if self.action == 'retrieve':
            return AssayDetailSerializer
        return AssayListSerializer

    @action(detail=True, methods=['get'])
    def data_series(self, request, pk=None):
        """List all data series for this assay."""
        assay = self.get_object()
        series = assay.data_series.select_related('compound', 'analysis')
        serializer = DataSeriesListSerializer(series, many=True)
        return Response(serializer.data)

    @action(detail=True, methods=['post'])
    def discover_series(self, request, pk=None):
        """
        Parse the data file and discover data series.

        This is a placeholder - actual implementation would parse
        the Excel/CSV file and create DataSeries objects.
        """
        assay = self.get_object()
        # TODO: Implement actual parsing logic from legacy code
        return Response({
            'status': 'not_implemented',
            'message': 'Data series discovery not yet implemented'
        })

    @action(detail=True, methods=['post'])
    def analyse_all(self, request, pk=None):
        """
        Run curve fitting analysis on all data series in this assay.

        Returns summary of analysis results.
        """
        import logging
        logger = logging.getLogger(__name__)

        try:
            assay = self.get_object()
            logger.info(f"Starting analysis for assay {assay.id} with {assay.data_series.count()} data series")
            results = analyse_assay(assay)
            logger.info(f"Analysis completed: {results['successful']} successful, {results['failed']} failed")
            return Response({
                'status': 'completed',
                'total': results['total'],
                'successful': results['successful'],
                'failed': results['failed'],
                'valid': results['valid'],
                'invalid': results['invalid'],
            })
        except Exception as e:
            logger.exception(f"Analysis failed for assay {pk}: {e}")
            return Response({
                'status': 'error',
                'error': str(e),
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

    @action(detail=True, methods=['post'])
    def upload_images(self, request, pk=None):
        """
        Batch upload images for Table-Of-Values analysis.

        For protocols using 'table_of_values' analysis method, images are
        pre-analyzed externally. This endpoint matches uploaded filenames
        to the 'Image File' value in each DataSeries analysis results.

        Expects multipart/form-data with one or more files.

        Returns summary of matched/unmatched files.
        """
        assay = self.get_object()
        files = request.FILES.getlist('files')

        if not files:
            return Response({
                'status': 'error',
                'error': 'No files provided'
            }, status=status.HTTP_400_BAD_REQUEST)

        # Get all data series for this assay that have analysis results
        data_series_qs = assay.data_series.select_related('analysis').filter(
            analysis__isnull=False
        )

        matched = []
        unmatched = []

        for f in files:
            found = False
            for series in data_series_qs:
                # Check if this file matches the 'Image File' in analysis results
                image_filename = series.analysis.results.get('Image File')
                if image_filename and image_filename == f.name:
                    # Save image to plot_image field (preserves original filename)
                    series.plot_image.save(f.name, f, save=True)
                    matched.append({
                        'filename': f.name,
                        'data_series_id': str(series.id),
                        'compound_name': series.compound_name,
                    })
                    found = True
                    break

            if not found:
                unmatched.append(f.name)

        return Response({
            'status': 'completed',
            'matched': len(matched),
            'unmatched': len(unmatched),
            'matched_files': matched,
            'unmatched_files': unmatched,
        })

    @action(detail=True, methods=['post'])
    def import_table_of_values(self, request, pk=None):
        """
        Bulk import pre-analyzed data for Table-Of-Values protocols.

        For protocols using 'table_of_values' analysis method, data comes
        pre-analyzed externally. This endpoint creates DataSeries and
        AnalysisResult records from the spreadsheet data.

        Expects JSON body:
        {
            "compound_column": "Compound",  // Column with compound names
            "kpi_column": "KPI",            // Column containing the KPI name
            "image_column": "Image File",   // Optional: column with image filenames
            "data": [
                {"Compound": "NCL-00042", "KPI": "EC50", "EC50": 234.5, ...},
                ...
            ]
        }

        The KPI column value (e.g., "EC50") must match a column name in the data.
        All rows must have the same KPI value.

        Returns summary of created records.
        """
        import re
        import logging
        logger = logging.getLogger(__name__)

        assay = self.get_object()

        # Validate protocol is table_of_values
        if assay.protocol.analysis_method != 'table_of_values':
            return Response({
                'status': 'error',
                'error': f"Protocol analysis method is '{assay.protocol.analysis_method}', not 'table_of_values'"
            }, status=status.HTTP_400_BAD_REQUEST)

        # Validate request data
        serializer = TableOfValuesImportSerializer(data=request.data)
        if not serializer.is_valid():
            return Response({
                'status': 'error',
                'errors': serializer.errors
            }, status=status.HTTP_400_BAD_REQUEST)

        validated = serializer.validated_data
        compound_column = validated['compound_column']
        kpi_column = validated['kpi_column']
        image_column = validated.get('image_column', '')
        data_rows = validated['data']

        created_series = []
        errors = []

        with reversion.create_revision():
            for idx, row in enumerate(data_rows):
                try:
                    compound_name = row.get(compound_column)
                    if not compound_name:
                        errors.append({
                            'row': idx,
                            'error': f"Missing compound name in column '{compound_column}'"
                        })
                        continue

                    # Get KPI value (the column name that holds the primary metric)
                    kpi_value = row.get(kpi_column)

                    # Build the results dictionary from all columns
                    results = {}
                    for col_name, col_value in row.items():
                        # Skip the compound column from results
                        if col_name == compound_column:
                            continue
                        results[col_name] = col_value

                    # Ensure KPI is properly set
                    results['KPI'] = kpi_value
                    results['fit_successful'] = True

                    # If there's an image column, ensure it's stored as 'Image File'
                    if image_column and image_column in row:
                        results['Image File'] = row[image_column]

                    # Try to match compound by name
                    compound = None
                    from compounds.registry.models import Compound
                    match = re.match(r'^NCL-0*(\d+)$', str(compound_name), re.IGNORECASE)
                    if match:
                        try:
                            reg_number = int(match.group(1))
                            compound = Compound.objects.filter(reg_number=reg_number).first()
                        except ValueError:
                            pass

                    # Create AnalysisResult
                    analysis = AnalysisResult.objects.create(
                        status='valid',
                        results=results
                    )

                    # Create DataSeries
                    series = DataSeries.objects.create(
                        assay=assay,
                        compound=compound,
                        compound_name=str(compound_name),
                        row=idx,
                        start_column=0,
                        end_column=0,
                        extracted_data={},
                        analysis=analysis
                    )

                    created_series.append({
                        'id': str(series.id),
                        'compound_name': compound_name,
                        'compound_matched': compound is not None,
                        'kpi': kpi_value,
                        'kpi_value': results.get(kpi_value) if kpi_value else None,
                    })

                except Exception as e:
                    logger.exception(f"Error creating series for row {idx}: {e}")
                    errors.append({
                        'row': idx,
                        'compound_name': row.get(compound_column, 'unknown'),
                        'error': str(e)
                    })

            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment(f"Imported {len(created_series)} data series from Table of Values")

        return Response({
            'status': 'completed',
            'created': len(created_series),
            'errors_count': len(errors),
            'series': created_series,
            'errors': errors if errors else None,
        })

    @action(detail=False, methods=['post'])
    def parse_adme_preview(self, request):
        """
        Parse an ADME file and return preview data for the frontend.

        Expects multipart/form-data with a single 'file' field.

        Returns:
        {
            "status": "success",
            "parser": {
                "vendor": "NCU",
                "assay_type": "liver_microsome_stability",
                "protocol_slug": "ncu-lm",
                "kpi_field": "t1_2_min"
            },
            "metadata": {...},
            "results": [
                {
                    "compound_id": "NCL-00042",
                    "compound_matched": true,
                    "compound_reg_number": "NCL-00000042",
                    "species": "Human",
                    "is_control": false,
                    "results": {...},
                    "kpi_value": 45.2
                },
                ...
            ],
            "errors": [...],
            "summary": {
                "total": 10,
                "matched": 8,
                "unmatched": 2,
                "controls": 1
            }
        }
        """
        import os
        import re
        import tempfile
        import logging
        from pathlib import Path

        logger = logging.getLogger(__name__)

        # Get uploaded file
        uploaded_file = request.FILES.get('file')
        if not uploaded_file:
            return Response({
                'status': 'error',
                'error': 'No file provided'
            }, status=status.HTTP_400_BAD_REQUEST)

        # Save to temp file for parsing - preserve original filename for parser detection
        # Parser detection relies on filename pattern (e.g., ADME-NCU-LM-YYYYMMDD.xlsx)
        import shutil
        tmp_dir = tempfile.mkdtemp()
        tmp_path = Path(tmp_dir) / uploaded_file.name
        with open(tmp_path, 'wb') as tmp:
            for chunk in uploaded_file.chunks():
                tmp.write(chunk)

        try:
            # Import and detect parser
            from compounds.assays.importers import detect_parser, ParseResult

            parser = detect_parser(tmp_path)
            if not parser:
                return Response({
                    'status': 'error',
                    'error': f'Could not detect ADME file type for: {uploaded_file.name}',
                    'hint': 'File should match pattern like ADME-NCU-LM-YYYYMMDD.xlsx'
                }, status=status.HTTP_400_BAD_REQUEST)

            # Parse the file
            parse_result: ParseResult = parser.parse(tmp_path)

            # Match compounds
            from compounds.registry.models import Compound

            results_preview = []
            matched_count = 0
            unmatched_count = 0
            control_count = 0

            for parsed in parse_result.results:
                if parsed.is_control:
                    control_count += 1
                    continue

                # Try to match compound
                compound = self._match_compound(parsed.compound_id)
                compound_matched = compound is not None

                if compound_matched:
                    matched_count += 1
                else:
                    unmatched_count += 1

                # Get KPI value
                kpi_value = parsed.results.get(parser.kpi_field)

                results_preview.append({
                    'compound_id': parsed.compound_id,
                    'compound_matched': compound_matched,
                    'compound_reg_number': compound.reg_number if compound else None,
                    'compound_db_id': str(compound.id) if compound else None,
                    'species': parsed.species,
                    'assay_format': parsed.assay_format,
                    'is_control': parsed.is_control,
                    'results': parsed.results,
                    'kpi_value': kpi_value,
                    'flags': parsed.results.get('flags', []),
                })

            # Format errors for response
            errors_list = [{
                'message': e.message,
                'severity': e.severity,
                'row': e.row,
                'column': e.column,
            } for e in parse_result.errors]

            return Response({
                'status': 'success',
                'filename': uploaded_file.name,
                'parser': {
                    'vendor': parser.vendor,
                    'assay_type': parser.assay_type,
                    'protocol_slug': parser.protocol_slug,
                    'kpi_field': parser.kpi_field,
                },
                'metadata': parse_result.metadata,
                'results': results_preview,
                'errors': errors_list,
                'summary': {
                    'total': len(parse_result.results),
                    'matched': matched_count,
                    'unmatched': unmatched_count,
                    'controls': control_count,
                }
            })

        except Exception as e:
            logger.exception(f"ADME parse failed: {e}")
            return Response({
                'status': 'error',
                'error': str(e)
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)

        finally:
            # Clean up temp directory and file
            try:
                shutil.rmtree(tmp_dir)
            except Exception:
                pass

    @action(detail=False, methods=['post'])
    def import_adme(self, request):
        """
        Import ADME data after preview confirmation.

        Expects JSON body:
        {
            "filename": "ADME-NCU-LM-20231218.xlsx",
            "parser_slug": "ncu-lm",
            "results": [
                {
                    "compound_id": "NCL-00042",
                    "compound_db_id": "uuid-or-null",
                    "species": "Human",
                    "results": {...}
                },
                ...
            ],
            "skip_unmatched": false,
            "comments": "Optional import notes"
        }

        Returns summary of created records.
        """
        import logging
        logger = logging.getLogger(__name__)

        # Validate request data
        filename = request.data.get('filename', 'ADME Import')
        parser_slug = request.data.get('parser_slug')
        results_data = request.data.get('results', [])
        skip_unmatched = request.data.get('skip_unmatched', False)
        comments = request.data.get('comments', '')
        target_id = request.data.get('target')  # Optional target UUID

        if not parser_slug:
            return Response({
                'status': 'error',
                'error': 'parser_slug is required'
            }, status=status.HTTP_400_BAD_REQUEST)

        if not results_data:
            return Response({
                'status': 'error',
                'error': 'No results to import'
            }, status=status.HTTP_400_BAD_REQUEST)

        # Get or create protocol
        from compounds.assays.importers import get_parser_by_slug
        parser = get_parser_by_slug(parser_slug)
        if not parser:
            return Response({
                'status': 'error',
                'error': f'Unknown parser: {parser_slug}'
            }, status=status.HTTP_400_BAD_REQUEST)

        protocol, _ = Protocol.objects.get_or_create(
            name=f"{parser.vendor} {parser.assay_type}",
            defaults={
                'analysis_method': 'pharmaron_adme',
                'comments': f'Auto-created for {parser.vendor} {parser.assay_type} imports',
            }
        )

        # Resolve target: use provided target_id, or fall back to protocol's default target
        from compounds.registry.models import Target
        target = None
        if target_id:
            try:
                target = Target.objects.get(id=target_id)
            except Target.DoesNotExist:
                logger.warning(f"Target {target_id} not found, ignoring")
        elif protocol.target:
            target = protocol.target

        created_series = []
        errors = []
        skipped = 0

        with reversion.create_revision():
            # Create Assay for this import
            assay = Assay.objects.create(
                protocol=protocol,
                target=target,
                created_by=request.user if request.user.is_authenticated else None,
                comments=f"Imported from {filename}" + (f"\n{comments}" if comments else ""),
            )

            from compounds.registry.models import Compound

            for idx, item in enumerate(results_data):
                try:
                    compound_id = item.get('compound_id')
                    compound_db_id = item.get('compound_db_id')
                    species = item.get('species')
                    results = item.get('results', {})

                    # Skip controls
                    if item.get('is_control'):
                        continue

                    # Get compound if matched
                    compound = None
                    if compound_db_id:
                        try:
                            compound = Compound.objects.get(id=compound_db_id)
                        except Compound.DoesNotExist:
                            pass

                    # Skip unmatched if requested
                    if not compound and skip_unmatched:
                        skipped += 1
                        continue

                    # Create AnalysisResult
                    analysis = AnalysisResult.objects.create(
                        status='valid',
                        results={
                            **results,
                            'KPI': parser.kpi_field,
                        }
                    )

                    # Create DataSeries
                    series = DataSeries.objects.create(
                        assay=assay,
                        compound=compound,
                        compound_name=str(compound_id),
                        row=idx,
                        start_column=0,
                        end_column=0,
                        extracted_data=results.get('time_course', {}),
                        analysis=analysis,
                    )

                    created_series.append({
                        'id': str(series.id),
                        'compound_id': compound_id,
                        'compound_matched': compound is not None,
                        'species': species,
                        'kpi_value': results.get(parser.kpi_field),
                    })

                except Exception as e:
                    logger.exception(f"Error creating series for {item.get('compound_id')}: {e}")
                    errors.append({
                        'compound_id': item.get('compound_id'),
                        'error': str(e)
                    })

            if request.user.is_authenticated:
                reversion.set_user(request.user)
            reversion.set_comment(f"Imported {len(created_series)} ADME data series from {filename}")

        return Response({
            'status': 'completed',
            'assay_id': str(assay.id),
            'created': len(created_series),
            'skipped': skipped,
            'errors_count': len(errors),
            'series': created_series,
            'errors': errors if errors else None,
        })

    def _match_compound(self, compound_id: str):
        """
        Try to match a compound ID to a registered compound.

        Supports:
        - NCL-XXXXXXXX format (direct reg_number match)
        - Plain numeric IDs (tries NCL-{padded})
        """
        import re
        from compounds.registry.models import Compound

        if not compound_id:
            return None

        compound_id = str(compound_id).strip()

        # Direct match on reg_number
        try:
            return Compound.objects.get(reg_number__iexact=compound_id)
        except Compound.DoesNotExist:
            pass

        # Try NCL format if numeric
        if compound_id.isdigit():
            ncl_id = f"NCL-{compound_id.zfill(8)}"
            try:
                return Compound.objects.get(reg_number__iexact=ncl_id)
            except Compound.DoesNotExist:
                pass

        # Try extracting number from NCL-XXXXX pattern
        match = re.match(r'NCL-?(\d+)', compound_id, re.IGNORECASE)
        if match:
            number = match.group(1)
            ncl_id = f"NCL-{number.zfill(8)}"
            try:
                return Compound.objects.get(reg_number__iexact=ncl_id)
            except Compound.DoesNotExist:
                pass

        return None


class DataSeriesViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Data Series."""
    queryset = DataSeries.objects.select_related('assay', 'compound', 'analysis', 'dilution_series')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['assay', 'compound']
    search_fields = ['compound_name']
    ordering_fields = ['compound_name']
    ordering = ['compound_name']

    def get_serializer_class(self):
        if self.action == 'create':
            return DataSeriesCreateSerializer
        if self.action == 'retrieve':
            return DataSeriesDetailSerializer
        return DataSeriesListSerializer

    @action(detail=True, methods=['post'])
    def analyse(self, request, pk=None):
        """
        Run curve fitting analysis on this data series.

        Returns the analysis result.
        """
        series = self.get_object()
        fitting_method = series.assay.protocol.get_effective_fitting_method()

        try:
            analysis = analyse_data_series(series, fitting_method)
            return Response({
                'status': analysis.status,
                'results': analysis.results,
                'kpi_value': analysis.kpi_value,
            })
        except Exception as e:
            return Response({
                'status': 'error',
                'error': str(e),
            }, status=status.HTTP_500_INTERNAL_SERVER_ERROR)


class AnalysisResultViewSet(viewsets.ReadOnlyModelViewSet):
    """Read-only operations for Analysis Results."""
    queryset = AnalysisResult.objects.all()
    serializer_class = AnalysisResultSerializer
    filter_backends = [DjangoFilterBackend]
    filterset_fields = ['status']


class ProtocolDocumentViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Protocol Documents."""
    queryset = ProtocolDocument.objects.select_related('protocol', 'created_by')
    serializer_class = ProtocolDocumentSerializer
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.OrderingFilter]
    filterset_fields = ['protocol']
    ordering = ['-created_at']

    def get_serializer_class(self):
        """Use create serializer for POST (file upload), read serializer otherwise."""
        if self.action == 'create':
            return ProtocolDocumentCreateSerializer
        return ProtocolDocumentSerializer

    def perform_create(self, serializer):
        with reversion.create_revision():
            instance = serializer.save(
                created_by=self.request.user if self.request.user.is_authenticated else None
            )
            if self.request.user.is_authenticated:
                reversion.set_user(self.request.user)
            reversion.set_comment("Uploaded via API")
            return instance


class HypothesisViewSet(ReversionMixin, viewsets.ModelViewSet):
    """CRUD operations for Hypotheses."""
    queryset = Hypothesis.objects.select_related('target', 'parent', 'product_compound')
    permission_classes = [IsContributorOrReadOnly]
    filter_backends = [DjangoFilterBackend, filters.SearchFilter, filters.OrderingFilter]
    filterset_fields = ['target', 'status']
    search_fields = ['smiles', 'rationale']
    ordering_fields = ['created_at', 'updated_at', 'status']
    ordering = ['-created_at']

    def get_serializer_class(self):
        if self.action == 'retrieve':
            return HypothesisDetailSerializer
        return HypothesisListSerializer

    @action(detail=False, methods=['get'])
    def by_status(self, request):
        """Group hypotheses by status."""
        from django.db.models import Count
        status_counts = (
            Hypothesis.objects
            .values('status')
            .annotate(count=Count('id'))
            .order_by('status')
        )
        return Response(list(status_counts))
