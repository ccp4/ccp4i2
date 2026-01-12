"""URL routing for users app."""

from django.urls import include, path
from rest_framework.routers import DefaultRouter

from .views import AdminListView, CurrentUserView, UserViewSet

router = DefaultRouter()
router.register('users', UserViewSet, basename='user')

app_name = 'users'

urlpatterns = [
    path('me/', CurrentUserView.as_view(), name='current-user'),
    path('admins/', AdminListView.as_view(), name='admin-list'),
    path('', include(router.urls)),
]
