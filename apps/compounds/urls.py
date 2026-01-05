"""
Compounds App URL Configuration

Combines routes from registry and assays sub-apps.
These URLs are included by the main ccp4i2 urls.py when compounds is enabled.
"""

from django.urls import path, include
from rest_framework.routers import DefaultRouter

from compounds.registry.views import (
    SupplierViewSet,
    TargetViewSet,
    CompoundViewSet,
    BatchViewSet,
    BatchQCFileViewSet,
    CompoundTemplateViewSet,
)
from compounds.assays.views import (
    DilutionSeriesViewSet,
    ProtocolViewSet,
    AssayViewSet,
    DataSeriesViewSet,
    AnalysisResultViewSet,
    HypothesisViewSet,
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
router.register(r'dilution-series', DilutionSeriesViewSet, basename='dilution-series')
router.register(r'protocols', ProtocolViewSet, basename='protocol')
router.register(r'assays', AssayViewSet, basename='assay')
router.register(r'data-series', DataSeriesViewSet, basename='data-series')
router.register(r'analysis-results', AnalysisResultViewSet, basename='analysis-result')
router.register(r'hypotheses', HypothesisViewSet, basename='hypothesis')

urlpatterns = [
    path('', include(router.urls)),
]
