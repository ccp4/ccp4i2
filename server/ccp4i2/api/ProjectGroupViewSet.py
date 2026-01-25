"""
ViewSet for managing project groups (campaigns).

Project groups allow organizing multiple projects together, particularly useful
for fragment screening campaigns where a parent project provides reference
coordinates and FreeR flags, and member projects each represent a dataset
soaked with a different compound.
"""
import logging
import os
import zipfile
import tempfile
from pathlib import Path
from django.http import FileResponse
from django.conf import settings
from rest_framework.viewsets import ModelViewSet
from rest_framework.parsers import FormParser, MultiPartParser, JSONParser
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from rest_framework import status
from . import serializers
from ..db import models
from ..lib.response import api_success, api_error

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectGroupViewSet(ModelViewSet):
    """
    ViewSet for managing project groups.

    Provides CRUD operations for ProjectGroup model plus custom actions
    for managing group memberships and retrieving related data.

    Actions:
        - list: List all project groups (supports ?type=fragment_set filter)
        - create: Create a new project group
        - retrieve: Get a project group with its memberships
        - update/partial_update: Update a project group
        - destroy: Delete a project group
        - parent_project: Get the parent project for a fragment set
        - member_projects: Get all member projects with job summaries
        - add_member: Add a project to the group
        - remove_member: Remove a project from the group
    """

    queryset = models.ProjectGroup.objects.prefetch_related("memberships").all()
    serializer_class = serializers.ProjectGroupSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]

    def get_serializer_class(self):
        """Use detailed serializer for retrieve action."""
        if self.action == "retrieve":
            return serializers.ProjectGroupDetailSerializer
        return serializers.ProjectGroupSerializer

    def get_queryset(self):
        """Filter by type if provided in query params."""
        queryset = super().get_queryset()
        group_type = self.request.query_params.get("type")
        if group_type:
            queryset = queryset.filter(type=group_type)
        return queryset

    @action(detail=True, methods=["get"], )
    def parent_project(self, request, pk=None):
        """
        Get the parent project for this group.

        For fragment sets, the parent project contains reference coordinates
        and FreeR flags that are used for all member datasets.

        Returns:
            Response: Serialized parent project data, or null if none set.
        """
        try:
            group = self.get_object()
            parent_membership = group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.PARENT
            ).first()

            if not parent_membership:
                return Response(None)

            serializer = serializers.ProjectSerializer(parent_membership.project)
            return Response(serializer.data)

        except Exception as e:
            logger.exception(
                "Failed to get parent project for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["get"], )
    def member_projects(self, request, pk=None):
        """
        Get all member projects with their job summaries and full job list.

        Returns member projects (not parent) along with:
        - Aggregated job status counts
        - Full list of jobs with status for icon display
        - KPIs from the most recent finished job

        Returns:
            Response: List of member projects with job data.
        """
        try:
            group = self.get_object()
            member_memberships = group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.MEMBER
            ).select_related("project")

            result = []
            for membership in member_memberships:
                project = membership.project
                project_data = serializers.ProjectSerializer(project).data

                # Get all jobs for this project with prefetched key values
                jobs = (
                    models.Job.objects.filter(project=project)
                    .prefetch_related("float_values", "char_values")
                    .order_by("number")
                )

                # Job summary counts
                job_summary = {
                    "total": jobs.count(),
                    "finished": sum(1 for j in jobs if j.status == models.Job.Status.FINISHED),
                    "failed": sum(1 for j in jobs if j.status == models.Job.Status.FAILED),
                    "running": sum(1 for j in jobs if j.status == models.Job.Status.RUNNING),
                    "pending": sum(1 for j in jobs if j.status == models.Job.Status.PENDING),
                }

                # Full job list for icon display
                jobs_list = []
                for job in jobs:
                    jobs_list.append({
                        "id": job.id,
                        "uuid": str(job.uuid),
                        "number": job.number,
                        "task_name": job.task_name,
                        "status": job.status,
                        "title": job.title,
                    })

                # Get KPIs - look through all jobs for the relevant key values
                # (matching legacy behavior: find last job with each key type)
                kpis = {}
                for job in reversed(list(jobs)):  # Most recent first
                    for fv in job.float_values.all():
                        if fv.key_id not in kpis:
                            kpis[fv.key_id] = fv.value
                    for cv in job.char_values.all():
                        if cv.key_id not in kpis:
                            kpis[cv.key_id] = cv.value

                project_data["job_summary"] = job_summary
                project_data["jobs"] = jobs_list
                project_data["kpis"] = kpis
                result.append(project_data)

            return Response(result)

        except Exception as e:
            logger.exception(
                "Failed to get member projects for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["get"], )
    def parent_files(self, request, pk=None):
        """
        Get reference files (coordinates and FreeR) from the parent project.

        Returns files with specific job_param_name values that indicate
        reference data for fragment screening:
        - XYZOUT: Output coordinates (PDB files)
        - FREEROUT: FreeR flag files

        Returns:
            Response: Dict with 'coordinates' and 'freer' file lists.
        """
        try:
            group = self.get_object()
            parent_membership = group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.PARENT
            ).first()

            if not parent_membership:
                return Response({"coordinates": [], "freer": []})

            parent_project = parent_membership.project

            # Get coordinate files (XYZOUT from coordinate_selector jobs)
            coord_files = models.File.objects.filter(
                job__project=parent_project,
                job_param_name="XYZOUT",
                type__name="chemical/x-pdb",
            )

            # Get FreeR files (FREEROUT from freerflag jobs)
            freer_files = models.File.objects.filter(
                job__project=parent_project,
                job_param_name="FREEROUT",
                type__name="application/CCP4-mtz-freerflag",
            )

            return Response(
                {
                    "coordinates": serializers.FileSerializer(
                        coord_files, many=True
                    ).data,
                    "freer": serializers.FileSerializer(freer_files, many=True).data,
                }
            )

        except Exception as e:
            logger.exception(
                "Failed to get parent files for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["post"], )
    def add_member(self, request, pk=None):
        """
        Add a project to the group.

        Request body:
            - project_id: ID of the project to add (required)
            - type: Membership type - 'parent' or 'member' (default: 'member')

        Note: A group can only have one parent project. Adding a new parent
        will fail if one already exists.

        Returns:
            Response: The created membership data.
        """
        try:
            group = self.get_object()
            project_id = request.data.get("project_id")
            membership_type = request.data.get(
                "type", models.ProjectGroupMembership.MembershipType.MEMBER
            )

            if not project_id:
                return api_error("project_id is required", status=400)

            # Validate membership type
            valid_types = [
                models.ProjectGroupMembership.MembershipType.PARENT,
                models.ProjectGroupMembership.MembershipType.MEMBER,
            ]
            if membership_type not in valid_types:
                return api_error(
                    f"Invalid membership type. Must be one of: {valid_types}",
                    status=400,
                )

            # Check if trying to add parent when one already exists
            if membership_type == models.ProjectGroupMembership.MembershipType.PARENT:
                existing_parent = group.memberships.filter(
                    type=models.ProjectGroupMembership.MembershipType.PARENT
                ).exists()
                if existing_parent:
                    return api_error(
                        "Group already has a parent project. Remove it first.",
                        status=400,
                    )

            # Get the project
            try:
                project = models.Project.objects.get(pk=project_id)
            except models.Project.DoesNotExist:
                return api_error("Project not found", status=404)

            # Check if project is already in group
            if group.memberships.filter(project=project).exists():
                return api_error("Project is already a member of this group", status=400)

            # Create membership
            membership = models.ProjectGroupMembership.objects.create(
                group=group, project=project, type=membership_type
            )

            serializer = serializers.ProjectGroupMembershipSerializer(membership)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        except Exception as e:
            logger.exception("Failed to add member to group %s", pk, exc_info=e)
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["delete"],
        ,
        url_path=r"members/(?P<project_id>\d+)",
    )
    def remove_member(self, request, pk=None, project_id=None):
        """
        Remove a project from the group.

        URL: DELETE /api/ccp4i2/projectgroups/{group_id}/members/{project_id}/

        Returns:
            Response: Success message on deletion.
        """
        try:
            group = self.get_object()

            try:
                membership = group.memberships.get(project_id=project_id)
            except models.ProjectGroupMembership.DoesNotExist:
                return api_error("Project is not a member of this group", status=404)

            project_name = membership.project.name
            membership.delete()

            return Response(
                {
                    "status": "success",
                    "message": f"Project '{project_name}' removed from group",
                }
            )

        except Exception as e:
            logger.exception(
                "Failed to remove member %s from group %s", project_id, pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=False, methods=["post"], , url_path="create_with_parent")
    def create_with_parent(self, request):
        """
        Create a new campaign with an auto-created parent project.

        This is the preferred way to create a fragment screening campaign.
        It creates both the ProjectGroup and a new Project to serve as the
        parent, linking them with a PARENT membership.

        Request body:
            - name: Name for both the campaign and parent project (required)
            - type: Group type, defaults to 'fragment_set'

        Returns:
            Response: The created ProjectGroup data.
        """
        try:
            name = request.data.get("name")
            group_type = request.data.get("type", "fragment_set")

            if not name:
                return api_error("name is required", status=400)

            # Check if a project with this name already exists (case-insensitive)
            if models.Project.objects.filter(name__iexact=name).exists():
                return api_error(
                    f"A project named '{name}' already exists. Please choose a different campaign name.",
                    status=400
                )

            # Check if a campaign with this name already exists
            if models.ProjectGroup.objects.filter(name__iexact=name).exists():
                return api_error(
                    f"A campaign named '{name}' already exists. Please choose a different name.",
                    status=400
                )

            # Create the parent project using the serializer
            # This ensures proper directory setup (CCP4_JOBS, CCP4_IMPORTED_FILES, etc.)
            project_serializer = serializers.ProjectSerializer(data={
                "name": name,
                "description": f"Parent project for {name} campaign",
                "directory": "__default__",  # Serializer will generate proper path
            })
            project_serializer.is_valid(raise_exception=True)
            project = project_serializer.save()
            logger.info("Created parent project %s for campaign %s", project.id, name)

            # Create the campaign (ProjectGroup)
            group = models.ProjectGroup.objects.create(
                name=name,
                type=group_type,
            )
            logger.info("Created campaign %s", group.id)

            # Link them with a PARENT membership
            models.ProjectGroupMembership.objects.create(
                group=group,
                project=project,
                type=models.ProjectGroupMembership.MembershipType.PARENT,
            )
            logger.info("Linked project %s as parent of campaign %s", project.id, group.id)

            serializer = serializers.ProjectGroupSerializer(group)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        except Exception as e:
            logger.exception("Failed to create campaign with parent", exc_info=e)
            return api_error(str(e), status=500)

    @action(detail=True, methods=["post"], )
    def set_parent(self, request, pk=None):
        """
        Set or change the parent project for this group.

        This is a convenience method that removes any existing parent
        and sets the specified project as the new parent.

        Request body:
            - project_id: ID of the project to set as parent (required)

        Returns:
            Response: The created parent membership data.
        """
        try:
            group = self.get_object()
            project_id = request.data.get("project_id")

            if not project_id:
                return api_error("project_id is required", status=400)

            # Get the project
            try:
                project = models.Project.objects.get(pk=project_id)
            except models.Project.DoesNotExist:
                return api_error("Project not found", status=404)

            # Remove existing parent if any
            group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.PARENT
            ).delete()

            # Remove project from members if it's already a member
            group.memberships.filter(project=project).delete()

            # Create new parent membership
            membership = models.ProjectGroupMembership.objects.create(
                group=group,
                project=project,
                type=models.ProjectGroupMembership.MembershipType.PARENT,
            )

            serializer = serializers.ProjectGroupMembershipSerializer(membership)
            return Response(serializer.data, status=status.HTTP_201_CREATED)

        except Exception as e:
            logger.exception("Failed to set parent for group %s", pk, exc_info=e)
            return api_error(str(e), status=500)

    @action(detail=False, methods=["get", "post"], )
    def project_campaigns(self, request):
        """
        Get campaign membership info for multiple projects.

        For GET requests (legacy, limited by URL length):
            Query params:
                - project_ids: Comma-separated list of project IDs
                - include_members: If 'true', also include campaigns where project is a member

        For POST requests (preferred for large numbers of projects):
            Request body (JSON):
                - project_ids: List of project IDs
                - include_members: Boolean, whether to include member campaigns

        Returns:
            Response: Dict mapping project_id to campaign info.
            Each entry includes:
                - campaign_id: The campaign ID
                - campaign_name: The campaign name
                - member_count: Number of member projects (for parent entries)
                - membership_type: 'parent' or 'member'
        """
        try:
            # Handle both GET (query params) and POST (request body)
            if request.method == "POST":
                project_ids = request.data.get("project_ids", [])
                include_members = request.data.get("include_members", False)
                # Handle case where project_ids might be passed as comma-separated string in POST too
                if isinstance(project_ids, str):
                    project_ids = [int(pid) for pid in project_ids.split(",") if pid]
            else:
                project_ids_param = request.query_params.get("project_ids", "")
                include_members = request.query_params.get("include_members", "false").lower() == "true"
                if not project_ids_param:
                    return Response({})
                project_ids = [int(pid) for pid in project_ids_param.split(",") if pid]

            if not project_ids:
                return Response({})

            result = {}

            # Find all campaigns where these projects are parents
            parent_memberships = models.ProjectGroupMembership.objects.filter(
                project_id__in=project_ids,
                type=models.ProjectGroupMembership.MembershipType.PARENT,
            ).select_related("group")

            for membership in parent_memberships:
                member_count = membership.group.memberships.filter(
                    type=models.ProjectGroupMembership.MembershipType.MEMBER
                ).count()
                result[str(membership.project_id)] = {
                    "campaign_id": membership.group.id,
                    "campaign_name": membership.group.name,
                    "member_count": member_count,
                    "membership_type": "parent",
                }

            # Optionally find campaigns where these projects are members
            if include_members:
                member_memberships = models.ProjectGroupMembership.objects.filter(
                    project_id__in=project_ids,
                    type=models.ProjectGroupMembership.MembershipType.MEMBER,
                ).select_related("group")

                for membership in member_memberships:
                    project_id_str = str(membership.project_id)
                    # Don't overwrite parent info if project is both parent and member
                    if project_id_str not in result:
                        result[project_id_str] = {
                            "campaign_id": membership.group.id,
                            "campaign_name": membership.group.name,
                            "membership_type": "member",
                        }

            return Response(result)

        except Exception as e:
            logger.exception("Failed to get project campaigns", exc_info=e)
            return api_error(str(e), status=500)

    @action(detail=True, methods=["get"], )
    def pandda_data(self, request, pk=None):
        """
        Get PANDDA-ready data from member projects.

        Finds i2Dimple and LidiaAcedrgNew jobs in member projects and
        returns information needed to build a PANDDA export.

        Returns:
            Response: Dict with 'datasets' list containing project name,
                      dimple job info, and acedrg job info for each member.
        """
        try:
            group = self.get_object()
            member_memberships = group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.MEMBER
            ).select_related("project")

            datasets = []
            for membership in member_memberships:
                project = membership.project

                # Find the most recent finished i2Dimple job
                dimple_job = (
                    models.Job.objects.filter(
                        project=project,
                        task_name__in=["i2Dimple", "dimple"],
                        status=models.Job.Status.FINISHED,
                    )
                    .order_by("-id")
                    .first()
                )

                # Find the most recent finished LidiaAcedrgNew job
                acedrg_job = (
                    models.Job.objects.filter(
                        project=project,
                        task_name__in=["LidiaAcedrgNew", "acedrg"],
                        status=models.Job.Status.FINISHED,
                    )
                    .order_by("-id")
                    .first()
                )

                if dimple_job:
                    datasets.append({
                        "project_name": project.name,
                        "project_id": project.id,
                        "dimple_job_id": dimple_job.id if dimple_job else None,
                        "acedrg_job_id": acedrg_job.id if acedrg_job else None,
                        "has_dimple": dimple_job is not None,
                        "has_acedrg": acedrg_job is not None,
                    })

            return Response({
                "datasets": datasets,
                "total_ready": len([d for d in datasets if d["has_dimple"]]),
                "total_with_dict": len([d for d in datasets if d["has_acedrg"]]),
            })

        except Exception as e:
            logger.exception(
                "Failed to get PANDDA data for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["post"], )
    def export_pandda(self, request, pk=None):
        """
        Build and return a PANDDA-ready ZIP file.

        Collects final.pdb, final.mtz from i2Dimple jobs and dict.cif
        from LidiaAcedrgNew jobs for all member projects with finished
        dimple jobs.

        Returns:
            FileResponse: ZIP file containing PANDDA-ready data structure:
                datasets/
                    <project_name>/
                        final.pdb
                        final.mtz
                        dict.cif (if available)
        """
        try:
            group = self.get_object()
            member_memberships = group.memberships.filter(
                type=models.ProjectGroupMembership.MembershipType.MEMBER
            ).select_related("project")

            # Create temporary directory for ZIP
            temp_dir = tempfile.mkdtemp(prefix="pandda_export_")
            zip_path = os.path.join(temp_dir, f"pandda_{group.name}.zip")

            included_count = 0

            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for membership in member_memberships:
                    project = membership.project

                    # Find the most recent finished i2Dimple job
                    dimple_job = (
                        models.Job.objects.filter(
                            project=project,
                            task_name__in=["i2Dimple", "dimple"],
                            status=models.Job.Status.FINISHED,
                        )
                        .order_by("-id")
                        .first()
                    )

                    if not dimple_job:
                        continue

                    # Dimple job directory
                    dimple_dir = Path(dimple_job.job_directory)
                    final_pdb = dimple_dir / "final.pdb"
                    final_mtz = dimple_dir / "final.mtz"

                    if not final_pdb.exists() or not final_mtz.exists():
                        logger.warning(
                            "Dimple outputs not found for project %s", project.name
                        )
                        continue

                    # Add dimple outputs
                    dataset_path = f"datasets/{project.name}"
                    zf.write(str(final_pdb), f"{dataset_path}/final.pdb")
                    zf.write(str(final_mtz), f"{dataset_path}/final.mtz")
                    included_count += 1

                    # Find LidiaAcedrgNew job for dictionary
                    acedrg_job = (
                        models.Job.objects.filter(
                            project=project,
                            task_name__in=["LidiaAcedrgNew", "acedrg"],
                            status=models.Job.Status.FINISHED,
                        )
                        .order_by("-id")
                        .first()
                    )

                    if acedrg_job:
                        acedrg_dir = Path(acedrg_job.job_directory)
                        # Check for both possible dictionary names
                        dict_cif = acedrg_dir / "LIG.cif"
                        if not dict_cif.exists():
                            dict_cif = acedrg_dir / "DRG.cif"

                        if dict_cif.exists():
                            zf.write(str(dict_cif), f"{dataset_path}/dict.cif")

            if included_count == 0:
                return api_error(
                    "No datasets with finished Dimple jobs found", status=400
                )

            # Return the ZIP file
            response = FileResponse(
                open(zip_path, "rb"),
                content_type="application/zip",
                as_attachment=True,
                filename=f"pandda_{group.name}.zip",
            )
            return response

        except Exception as e:
            logger.exception(
                "Failed to export PANDDA for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)
