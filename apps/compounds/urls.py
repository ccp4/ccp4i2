"""
Compounds App URL Configuration

Combines routes from registry, assays, and constructs sub-apps.
These URLs are included by the main ccp4i2 urls.py when compounds is enabled.
"""

from django.urls import path, include
from rest_framework.routers import DefaultRouter

from compounds.admin_views import import_legacy_fixtures, import_status
from compounds.registry.views import (
    SupplierViewSet,
    TargetViewSet,
    CompoundViewSet,
    BatchViewSet,
    BatchQCFileViewSet,
    CompoundTemplateViewSet,
)
from compounds.assays.views import (
    FittingMethodViewSet,
    DilutionSeriesViewSet,
    ProtocolViewSet,
    ProtocolDocumentViewSet,
    AssayViewSet,
    DataSeriesViewSet,
    AnalysisResultViewSet,
    HypothesisViewSet,
)
from compounds.assays.aggregation_views import AggregationViewSet
from compounds.constructs.views import (
    ConstructProjectViewSet,
    PlasmidViewSet,
    ProteinViewSet,
    ProteinSynonymViewSet,
    ProteinUseViewSet,
    CassetteViewSet,
    CassetteUseViewSet,
    SequencingResultViewSet,
    ExpressionTagTypeViewSet,
    ProteaseViewSet,
    ExpressionTagLocationViewSet,
    ExpressionTagViewSet,
)

router = DefaultRouter()

# Registry routes
router.register(r'suppliers', SupplierViewSet, basename='supplier')
router.register(r'targets', TargetViewSet, basename='target')
router.register(r'compounds', CompoundViewSet, basename='compound')
router.register(r'batches', BatchViewSet, basename='batch')
router.register(r'batch-qc-files', BatchQCFileViewSet, basename='batch-qc-file')
router.register(r'compound-templates', CompoundTemplateViewSet, basename='compound-template')

# Assay routes
router.register(r'fitting-methods', FittingMethodViewSet, basename='fitting-method')
router.register(r'dilution-series', DilutionSeriesViewSet, basename='dilution-series')
router.register(r'protocols', ProtocolViewSet, basename='protocol')
router.register(r'protocol-documents', ProtocolDocumentViewSet, basename='protocol-document')
router.register(r'assays', AssayViewSet, basename='assay')
router.register(r'data-series', DataSeriesViewSet, basename='data-series')
router.register(r'analysis-results', AnalysisResultViewSet, basename='analysis-result')
router.register(r'hypotheses', HypothesisViewSet, basename='hypothesis')

# Aggregation routes
router.register(r'aggregations', AggregationViewSet, basename='aggregation')

# Construct Database routes
router.register(r'construct-projects', ConstructProjectViewSet, basename='construct-project')
router.register(r'plasmids', PlasmidViewSet, basename='plasmid')
router.register(r'proteins', ProteinViewSet, basename='protein')
router.register(r'protein-synonyms', ProteinSynonymViewSet, basename='protein-synonym')
router.register(r'protein-uses', ProteinUseViewSet, basename='protein-use')
router.register(r'cassettes', CassetteViewSet, basename='cassette')
router.register(r'cassette-uses', CassetteUseViewSet, basename='cassette-use')
router.register(r'sequencing-results', SequencingResultViewSet, basename='sequencing-result')
router.register(r'expression-tag-types', ExpressionTagTypeViewSet, basename='expression-tag-type')
router.register(r'proteases', ProteaseViewSet, basename='protease')
router.register(r'expression-tag-locations', ExpressionTagLocationViewSet, basename='expression-tag-location')
router.register(r'expression-tags', ExpressionTagViewSet, basename='expression-tag')

urlpatterns = [
    # Admin endpoints for legacy data import
    path('admin/import-legacy/', import_legacy_fixtures, name='import-legacy-fixtures'),
    path('admin/import-status/', import_status, name='import-status'),

    # Router-based endpoints
    path('', include(router.urls)),
]
