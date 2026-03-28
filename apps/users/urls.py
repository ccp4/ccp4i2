# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
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
