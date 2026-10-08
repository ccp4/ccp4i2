"""
ViewSet for managing project groups (campaigns).

Project groups allow organizing multiple projects together, particularly useful
for fragment screening campaigns where a parent project provides reference
coordinates and FreeR flags, and member projects each represent a dataset
soaked with a different compound.
"""
import logging
import math
from django.db.models import Count, Q
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
from ..lib.kpi_values import kpi_map
from ..lib.response import api_success, api_error
from ..lib import pandda_export
from ..wrappers.pandda_campaign.script import pandda_invocation
from ..lib import campaign_scene
from ..lib import campaign_events
from ..lib.utils.jobs import pandda_site_index

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

    queryset = models.ProjectGroup.objects.all()
    serializer_class = serializers.ProjectGroupSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]

    def get_serializer_class(self):
        """Use detailed serializer for retrieve action."""
        if self.action == "retrieve":
            return serializers.ProjectGroupDetailSerializer
        return serializers.ProjectGroupSerializer

    def get_queryset(self):
        """Filter by type, and count members in the same query as the rows.

        Two costs per campaign, both paid on every listing, neither visible
        from the page that pays them.

        ``member_count`` was a ``SerializerMethodField`` calling
        ``obj.memberships.filter(...).count()``. A ``.filter()`` on a related
        manager cannot use ``prefetch_related`` -- the prefetch caches the
        whole set, and filtering it goes back to the database -- so the old
        ``prefetch_related("memberships")`` was loaded and then ignored, and a
        COUNT ran per campaign anyway. An annotation gets it in the same query
        as the rows.

        ``ProjectGroupSerializer`` declares ``fields = "__all__"``, which
        includes the ``projects`` M2M, so DRF fetched every member project id
        of every campaign -- another query each, and a payload nothing on the
        landing page reads. Prefetching makes that one query for the whole
        page; dropping the field would be better still, but it is part of the
        response shape and not this change's to remove.
        """
        queryset = super().get_queryset().annotate(
            member_total=Count(
                "memberships",
                filter=Q(memberships__type=models.ProjectGroupMembership
                         .MembershipType.MEMBER),
                distinct=True,
            )
        )
        queryset = queryset.prefetch_related(
            "projects" if self.action == "list" else "memberships")
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
        - The dataset x site matrix: ``site_cells`` (one entry per campaign
          site: the nearest PanDDA event within the site's radius, and the
          verdict), ``frame_mismatch`` and ``current_model_job``

        A fixed number of queries for the whole page, whatever the number of
        members: every member's jobs, KPI values, verdicts and events are
        fetched once each and grouped here, and no file is opened per row
        (the dataset's cell rides on its CampaignEvent rows; the parent's is
        one header read, cached).

        Returns:
            Response: List of member projects with job data.
        """
        try:
            group = self.get_object()
            member_memberships = list(
                group.memberships.filter(
                    type=models.ProjectGroupMembership.MembershipType.MEMBER
                ).select_related("project").prefetch_related("project__tags")
            )
            project_ids = [m.project_id for m in member_memberships]

            # Every site of this campaign, so a row can say how much of its
            # evaluation is outstanding and carry one cell per site. The
            # denominator is the campaign's current site count, which keeps
            # the fraction comparable down the column rather than varying per
            # dataset.
            sites = list(group.site_set.all())
            site_total = len(sites)

            # Every verdict at this campaign's sites for these datasets: one
            # query, restricted to this campaign, because a project can belong
            # to more than one.
            evaluations_by_project = {}
            verdicts = {}
            for evaluation in models.SiteEvaluation.objects.filter(
                site__group=group, project_id__in=project_ids
            ).select_related("site"):
                evaluations_by_project.setdefault(
                    evaluation.project_id, []).append(evaluation)
                verdicts[(evaluation.project_id, evaluation.site_id)] = (
                    evaluation.verdict)

            # Every job of every member, and their KPI values, in three
            # queries. The values are selected by project rather than
            # prefetched by job id, so the parameter count is bounded by the
            # number of members rather than of jobs.
            jobs_by_project = {pid: [] for pid in project_ids}
            for job in models.Job.objects.filter(
                project_id__in=project_ids
            ).order_by("number"):
                jobs_by_project[job.project_id].append(job)
            floats_by_job = {}
            for value in models.JobFloatValue.objects.filter(
                job__project_id__in=project_ids
            ):
                floats_by_job.setdefault(value.job_id, []).append(value)
            chars_by_job = {}
            for value in models.JobCharValue.objects.filter(
                job__project_id__in=project_ids
            ):
                chars_by_job.setdefault(value.job_id, []).append(value)

            matrix = campaign_events.site_matrix(
                group, sites, jobs_by_project, verdicts
            )

            result = []
            for membership in member_memberships:
                project = membership.project
                project_data = serializers.ProjectSerializer(project).data
                jobs = jobs_by_project[project.id]

                # Job summary counts
                job_summary = {
                    "total": len(jobs),
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
                # (matching legacy behavior: find last job with each key type).
                # kpi_map drops values JSON cannot carry, so one NaN R-factor
                # cannot take the whole campaign overview down with it.
                kpis = {}
                for job in reversed(jobs):  # Most recent first
                    context = f"job {job.number} of {project.name}"
                    for key, value in kpi_map(floats_by_job.get(job.id, ()), context).items():
                        kpis.setdefault(key, value)
                    for key, value in kpi_map(chars_by_job.get(job.id, ()), context).items():
                        kpis.setdefault(key, value)

                # What was found at each site, and how much is still unlooked
                # at, as chips. Kept for compatibility alongside site_cells.
                #
                # Only hit and unclear are listed. A rich campaign has 30-40
                # sites and most are empty for most datasets, so sending every
                # verdict would put a 40-entry list on every row to render a
                # cell that shows two chips. "Empty" is recoverable from the
                # evaluated count, and the per-dataset view has the detail.
                evaluations = []
                evaluated = 0
                for evaluation in evaluations_by_project.get(project.id, ()):
                    evaluated += 1
                    if evaluation.verdict in ("hit", "unclear"):
                        evaluations.append({
                            "site_id": evaluation.site_id,
                            "site_name": evaluation.site.name,
                            "verdict": evaluation.verdict,
                        })

                # Hits first, so a dataset with many unclears never has its
                # hits pushed out of a capped chip list.
                evaluations.sort(
                    key=lambda e: (e["verdict"] != "hit", e["site_name"])
                )

                project_data["job_summary"] = job_summary
                project_data["jobs"] = jobs_list
                project_data["kpis"] = kpis
                project_data["site_evaluations"] = evaluations
                project_data["sites_evaluated"] = evaluated
                project_data["sites_total"] = site_total
                # The dataset x site matrix (lib/campaign_events.site_matrix):
                # one cell per site, keyed by site uuid -- the nearest PanDDA
                # event within the site's radius, or null, and the verdict, or
                # null for "nobody has looked" -- plus frame_mismatch and
                # current_model_job.
                project_data.update(matrix[project.id])
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

    @action(detail=False, methods=["post"], url_path="create_with_parent")
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
            datasets = []
            for project, dimple_job, acedrg_job in pandda_export.collect_pandda_datasets(group):
                if not dimple_job:
                    continue
                datasets.append({
                    "project_name": project.name,
                    "project_id": project.id,
                    "dimple_job_id": dimple_job.id,
                    "acedrg_job_id": acedrg_job.id if acedrg_job else None,
                    "has_dimple": True,
                    "has_acedrg": acedrg_job is not None,
                    "resolution": pandda_export.dataset_resolution(dimple_job),
                })

            # What PanDDA would actually process this set at, and whether the
            # data or the comparator floor decided that. A campaign smaller
            # than max_shell_datasets is processed at the resolution of its
            # worst crystal however good the rest are, so the client can warn
            # before anyone pays for a run.
            resolutions = [d["resolution"] for d in datasets if d["resolution"]]
            shell, dragged = pandda_invocation.shell_resolution(resolutions)
            return Response({
                "datasets": datasets,
                "total_ready": len([d for d in datasets if d["has_dimple"]]),
                "total_with_dict": len([d for d in datasets if d["has_acedrg"]]),
                "resolution": {
                    "best": min(resolutions) if resolutions else None,
                    "worst": max(resolutions) if resolutions else None,
                    "processing": shell,
                    "dragged_by_comparator_floor": dragged,
                    "max_shell_datasets": pandda_invocation.MAX_SHELL_DATASETS,
                },
            })

        except Exception as e:
            logger.exception(
                "Failed to get PANDDA data for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["get"], )
    def summary_scene(self, request, pk=None):
        """
        Build a Moorhen summary scene for a fragment campaign.

        Returns a scene that overlays every discovered fragment hit on the
        parent reference structure: the parent as a ribbon, and each hit
        dataset's ligand as sticks scoped to its own restraint dictionary.
        A dataset counts as a hit when its refined coordinates contain the
        ligand its dictionary describes (see ``lib.campaign_scene``).

        Returns:
            Response: ``{"scene": <MoorhenScene>, "stats": {...}}`` where
                ``scene`` is a JSON-serialisable scene the frontend applies
                directly via the scene resolver.
        """
        try:
            group = self.get_object()
            return Response(campaign_scene.build_summary_scene(group))
        except Exception as e:
            logger.exception(
                "Failed to build summary scene for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["get"],
        url_path=r"sites/(?P<site_id>[0-9]+)/scene",
    )
    def site_scene(self, request, pk=None, site_id=None):
        """
        Build a Moorhen scene for one binding site of a fragment campaign.

        The exemplar as a ribbon with the pocket residues as sticks, and the
        ligand of every dataset judged a hit at this site drawn on top, each
        fitted onto the exemplar locally on the site. Membership is decided
        by the recorded verdicts alone (see ``lib.campaign_scene``).

        Query parameters:
            include=unclear: also draw datasets with an ``unclear`` verdict,
                in a muted colour.
            superpose=none: draw every dataset in its own frame. Superposition
                is the default, so the parameter turns it off -- for seeing the
                frames as deposited, or diagnosing a fit that went wrong.

        Returns:
            Response: ``{"scene": <MoorhenScene>, "stats": {...}}``, the
                same envelope as ``summary_scene``.
        """
        try:
            group = self.get_object()
            try:
                site = group.site_set.get(id=site_id)
            except models.CampaignSite.DoesNotExist:
                return api_error("Site not found in this campaign", status=404)

            include = {
                token.strip()
                for token in request.query_params.get("include", "").split(",")
            }
            return Response(
                campaign_scene.build_site_scene(
                    group,
                    site,
                    include_unclear="unclear" in include,
                    superpose=request.query_params.get("superpose") != "none",
                )
            )
        except Exception as e:
            logger.exception(
                "Failed to build scene for site %s of group %s", site_id, pk,
                exc_info=e,
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
            result = pandda_export.build_pandda_zip(group)

            if result.included_count == 0:
                return api_error(
                    "No datasets with finished Dimple jobs found", status=400
                )

            return FileResponse(
                open(result.zip_path, "rb"),
                content_type="application/zip",
                as_attachment=True,
                filename=result.zip_path.name,
            )

        except Exception as e:
            logger.exception(
                "Failed to export PANDDA for group %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @staticmethod
    def _site_payload(site):
        """The wire shape of a site.

        `id` is new in migration 0024 and is what everything else refers to: a
        site used to be identified by its name and its index in a JSON list, so
        renaming one silently orphaned every reference to it. The rest of the
        shape is unchanged, so existing readers keep working.
        """
        payload = {
            "id": site.id,
            "uuid": str(site.uuid),
            "name": site.name,
            "origin": [site.origin_x, site.origin_y, site.origin_z],
            "order": site.order,
            # How close (A) a PanDDA event must be to count as at this site.
            "radius": site.radius,
        }
        # Present only where the caller annotated it (the list view). Deleting
        # a site deletes every verdict recorded there, across every dataset,
        # and the difference between losing nothing and losing forty datasets'
        # work is the whole content of that warning -- without it a
        # confirmation is just a button to click through.
        count = getattr(site, "evaluation_count", None)
        if count is not None:
            payload["evaluation_count"] = count
        if site.quat:
            payload["quat"] = site.quat
        if site.zoom is not None:
            payload["zoom"] = site.zoom
        return payload

    @staticmethod
    def _parse_site(data, index=0):
        """Validate one incoming site. Returns (fields, error_message)."""
        if not isinstance(data, dict):
            return None, f"Site {index} must be an object"
        if "name" not in data:
            return None, f"Site {index} missing required field 'name'"
        origin = data.get("origin")
        if not isinstance(origin, list) or len(origin) != 3:
            return None, f"Site {index} 'origin' must be an array of 3 numbers"
        try:
            x, y, z = (float(v) for v in origin)
        except (TypeError, ValueError):
            return None, f"Site {index} 'origin' must be an array of 3 numbers"

        quat = data.get("quat")
        if quat is not None and (not isinstance(quat, list) or len(quat) != 4):
            return None, f"Site {index} 'quat' must be an array of 4 numbers"

        zoom = data.get("zoom")
        if zoom is not None:
            try:
                zoom = float(zoom)
            except (TypeError, ValueError):
                return None, f"Site {index} 'zoom' must be a number"

        fields = {
            "name": str(data["name"])[:100],
            "origin_x": x,
            "origin_y": y,
            "origin_z": z,
            "quat": quat,
            "zoom": zoom,
        }
        if data.get("radius") is not None:
            radius, error = ProjectGroupViewSet._parse_radius(data["radius"])
            if error:
                return None, f"Site {index} {error}"
            fields["radius"] = radius
        return fields, None

    @staticmethod
    def _parse_radius(value):
        """A site radius: a positive, finite number of Angstroms. (value, error)."""
        try:
            radius = float(value)
        except (TypeError, ValueError):
            return None, "'radius' must be a number"
        if not math.isfinite(radius) or radius <= 0:
            return None, "'radius' must be a positive number of Angstroms"
        return radius, None

    @action(detail=True, methods=["get", "post"])
    def sites(self, request, pk=None):
        """
        List this campaign's binding sites, or add one.

        GET returns the list, ordered. POST adds a single site and returns it.

        Sites were a JSON list on the group until migration 0024, replaced
        wholesale on every write. Adding one at a time means two people editing
        a campaign no longer overwrite each other: the old PUT sent the entire
        array, so the last writer silently discarded the other's new sites.
        Use PATCH/DELETE on `sites/<id>/` to change or remove one.
        """
        try:
            group = self.get_object()

            if request.method == "GET":
                sites = group.site_set.annotate(
                    evaluation_count=Count("evaluations")
                )
                return Response([self._site_payload(s) for s in sites])

            fields, error = self._parse_site(request.data)
            if error:
                return api_error(error, status=400)

            if group.site_set.filter(name=fields["name"]).exists():
                return api_error(
                    f"This campaign already has a site called '{fields['name']}'",
                    status=409,
                )

            last = group.site_set.order_by("-order").first()
            site = models.CampaignSite.objects.create(
                group=group, order=(last.order + 1) if last else 0, **fields
            )
            logger.info("Added site '%s' to campaign %s", site.name, pk)
            return Response(self._site_payload(site), status=201)

        except Exception as e:
            logger.exception("Failed to manage sites for group %s", pk, exc_info=e)
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["get"],
        url_path=r"evaluations/(?P<project_id>[0-9]+)",
    )
    def project_evaluations(self, request, pk=None, project_id=None):
        """Every verdict recorded for one dataset, including the empties.

        The overview payload in `member_projects` carries only hits and
        unclears, because it renders a row per dataset and most sites are
        empty for most datasets. A view of a single dataset can afford the
        whole picture and needs it: without the empties it cannot tell
        "somebody looked and found nothing" from "nobody has looked yet",
        which is the distinction these rows exist to keep. A control that
        confused the two would write the wrong one back.

        Ordered by site, so a caller can zip this against the site list.
        """
        try:
            group = self.get_object()

            if not group.memberships.filter(project_id=project_id).exists():
                return api_error(
                    "That project is not part of this campaign", status=404
                )

            evaluations = (
                models.SiteEvaluation.objects.filter(
                    project_id=project_id, site__group=group
                )
                .select_related("site")
                .order_by("site__order", "site__id")
            )

            return Response([
                {
                    "site_id": evaluation.site_id,
                    "site_name": evaluation.site.name,
                    "project_id": int(project_id),
                    "verdict": evaluation.verdict,
                    "evaluator": evaluation.evaluator,
                    "note": evaluation.note,
                }
                for evaluation in evaluations
            ])

        except Exception as e:
            logger.exception(
                "Failed to list evaluations of project %s in group %s",
                project_id, pk, exc_info=e,
            )
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["put", "delete"],
        url_path=r"sites/(?P<site_id>[0-9]+)/evaluation/(?P<project_id>[0-9]+)",
    )
    def site_evaluation(self, request, pk=None, site_id=None, project_id=None):
        """Record, change or withdraw what was found at one site in one dataset.

        PUT with {"verdict": "hit"|"empty"|"unclear"} and an optional
        "evaluator" and "note". DELETE withdraws the verdict.

        Withdrawing is not the same as recording "empty": no row means nobody
        has looked, while "empty" asserts that somebody looked and found
        nothing. That distinction is the reason these are rows rather than
        tags, so DELETE really does remove the row rather than writing an
        "empty" over it.
        """
        try:
            group = self.get_object()

            try:
                site = group.site_set.get(id=site_id)
            except models.CampaignSite.DoesNotExist:
                return api_error("Site not found in this campaign", status=404)

            if not group.memberships.filter(project_id=project_id).exists():
                return api_error(
                    "That project is not part of this campaign", status=404
                )

            if request.method == "DELETE":
                deleted, _ = models.SiteEvaluation.objects.filter(
                    project_id=project_id, site=site
                ).delete()
                return Response(status=204 if deleted else 404)

            verdict = (request.data or {}).get("verdict")
            valid = [choice[0] for choice in models.SiteEvaluation.Verdict.choices]
            if verdict not in valid:
                return api_error(
                    f"verdict must be one of {', '.join(valid)}", status=400
                )

            evaluation, created = models.SiteEvaluation.objects.update_or_create(
                project_id=project_id,
                site=site,
                defaults={
                    "verdict": verdict,
                    "evaluator": str(request.data.get("evaluator", ""))[:150],
                    "note": str(request.data.get("note", "")),
                },
            )
            logger.info(
                "%s %s at site '%s' for project %s",
                "Recorded" if created else "Updated", verdict, site.name, project_id,
            )
            return Response(
                {
                    "site_id": site.id,
                    "site_name": site.name,
                    "project_id": int(project_id),
                    "verdict": evaluation.verdict,
                    "evaluator": evaluation.evaluator,
                    "note": evaluation.note,
                },
                status=201 if created else 200,
            )

        except Exception as e:
            logger.exception(
                "Failed to set evaluation for site %s of group %s",
                site_id, pk, exc_info=e,
            )
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["patch", "delete"],
        url_path=r"sites/(?P<site_id>[0-9]+)",
    )
    def site_detail(self, request, pk=None, site_id=None):
        """Update or delete one site, addressed by its stable id.

        Renaming through here keeps every evaluation attached, which is the
        point of sites having ids at all.
        """
        try:
            group = self.get_object()
            try:
                site = group.site_set.get(id=site_id)
            except models.CampaignSite.DoesNotExist:
                return api_error("Site not found in this campaign", status=404)

            if request.method == "DELETE":
                name = site.name
                site.delete()  # cascades to this site's evaluations
                logger.info("Deleted site '%s' from campaign %s", name, pk)
                return Response(status=204)

            data = request.data
            if not isinstance(data, dict):
                return api_error("Expected an object", status=400)

            if "name" in data:
                name = str(data["name"])[:100]
                if group.site_set.filter(name=name).exclude(id=site.id).exists():
                    return api_error(
                        f"This campaign already has a site called '{name}'",
                        status=409,
                    )
                site.name = name

            if "origin" in data:
                origin = data["origin"]
                if not isinstance(origin, list) or len(origin) != 3:
                    return api_error(
                        "'origin' must be an array of 3 numbers", status=400
                    )
                try:
                    site.origin_x, site.origin_y, site.origin_z = (
                        float(v) for v in origin
                    )
                except (TypeError, ValueError):
                    return api_error(
                        "'origin' must be an array of 3 numbers", status=400
                    )

            if "quat" in data:
                quat = data["quat"]
                if quat is not None and (
                    not isinstance(quat, list) or len(quat) != 4
                ):
                    return api_error(
                        "'quat' must be an array of 4 numbers", status=400
                    )
                site.quat = quat

            if "zoom" in data:
                zoom = data["zoom"]
                if zoom is not None:
                    try:
                        zoom = float(zoom)
                    except (TypeError, ValueError):
                        return api_error("'zoom' must be a number", status=400)
                site.zoom = zoom

            if "order" in data:
                try:
                    site.order = int(data["order"])
                except (TypeError, ValueError):
                    return api_error("'order' must be an integer", status=400)

            if "radius" in data:
                radius, error = self._parse_radius(data["radius"])
                if error:
                    return api_error(error, status=400)
                site.radius = radius

            site.save()
            return Response(self._site_payload(site))

        except Exception as e:
            logger.exception(
                "Failed to update site %s of group %s", site_id, pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(detail=True, methods=["get"], url_path="pandda-sites")
    def pandda_sites(self, request, pk=None):
        """The sites a PanDDA run found over this campaign, and what sits at each.

        The run's own site axis, not the curated one: membership is PanDDA's
        ``site_idx`` and no verdict is needed, which is what makes it usable
        the moment a run has been fanned out. ``sites`` is a navigable index
        -- per site its centroid, its rollups, and every event with the
        receipt job to open it in -- and ``runs`` lists every run this
        campaign holds receipts from, so a caller can offer the choice.

        Query parameters:
            run=<uuid>: index this run's receipts. Default: the run whose
                newest receipt is newest.

        Returns:
            Response: ``{"run", "runs", "tree", "sites", "unsited", "stats"}``
                (see ``lib.utils.jobs.pandda_site_index``).
        """
        try:
            group = self.get_object()
            return Response(pandda_site_index.build_site_index(
                group, run_job_uuid=request.query_params.get("run") or None))
        except Exception as e:
            logger.exception("Failed to index PanDDA sites for group %s", pk, exc_info=e)
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["get"],
        url_path=r"pandda-sites/(?P<site_idx>[0-9]+)/scene",
    )
    def pandda_site_scene(self, request, pk=None, site_idx=None):
        """A Moorhen scene of every ligand PanDDA built at one site of one run.

        The exemplar drawn once, and the autobuilt pose of every event at this
        site on top of it, each fitted onto the exemplar. No verdict is
        required and none is implied: these are candidates to be judged
        together, which is the view the curated site scene cannot give before
        anybody has judged them.

        Query parameters:
            run=<uuid>: which run's sites these are. Default: the latest.
            superpose=none: draw every dataset in its own frame.

        Returns:
            Response: ``{"scene": <MoorhenScene>, "stats": {...}}``, the same
                envelope as ``site_scene``.
        """
        try:
            group = self.get_object()
            index = pandda_site_index.build_site_index(
                group, run_job_uuid=request.query_params.get("run") or None)
            wanted = int(site_idx)
            site = next((s for s in index["sites"] if s["site_idx"] == wanted), None)
            if site is None:
                return api_error(
                    f"Site {wanted} is not a site of this run", status=404)
            return Response(campaign_scene.build_run_site_scene(
                group, site, run=index["run"],
                superpose=request.query_params.get("superpose") != "none",
            ))
        except Exception as e:
            logger.exception(
                "Failed to build PanDDA site scene %s for group %s", site_idx, pk,
                exc_info=e)
            return api_error(str(e), status=500)
