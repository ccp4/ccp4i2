"""Admin configuration for Azure Extensions models."""

from django.contrib import admin

from .models import StagedUpload


@admin.register(StagedUpload)
class StagedUploadAdmin(admin.ModelAdmin):
    list_display = ["uuid", "upload_type", "status", "original_filename", "created_at"]
    list_filter = ["status", "upload_type"]
    search_fields = ["uuid", "original_filename", "requested_by"]
    readonly_fields = ["uuid", "created_at"]
    ordering = ["-created_at"]
