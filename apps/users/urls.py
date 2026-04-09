"""URL routing for users app."""

from django.urls import include, path
from rest_framework.routers import DefaultRouter

from .views import AdminListView, CurrentUserView, OperatingLevelView, UserViewSet

router = DefaultRouter()
router.register('users', UserViewSet, basename='user')

app_name = 'users'

urlpatterns = [
    path('me/', CurrentUserView.as_view(), name='current-user'),
    path('me/operating-level/', OperatingLevelView.as_view(), name='operating-level'),
    path('admins/', AdminListView.as_view(), name='admin-list'),
    path('', include(router.urls)),
]
