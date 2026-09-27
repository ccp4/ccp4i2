"""``member_count`` on the campaigns list: right, and flat in the number of
campaigns.

The count itself is not new. What it cost was: a ``SerializerMethodField``
calling ``obj.memberships.filter(...).count()`` cannot use the viewset's
``prefetch_related`` -- filtering a related manager goes back to the database
-- so the list ran one extra query per campaign while also prefetching every
membership row it then ignored. On an instance with fifty campaigns that is
fifty avoidable round trips on a page whose whole job is to render fifty
names.

So the assertion that matters here is not a query number, which would be
brittle, but that the number does not grow when the campaigns do.

Uses the api/ conftest, which auto-applies django_db(transaction=True) and
sets AllowAny on the viewsets -- do NOT add @pytest.mark.django_db here.
"""

from django.db import connection
from django.test.utils import CaptureQueriesContext
from rest_framework.test import APIClient

from ccp4i2.db import models

MEMBER = models.ProjectGroupMembership.MembershipType.MEMBER
PARENT = models.ProjectGroupMembership.MembershipType.PARENT


def _campaign(root, name, n_members, with_parent=True):
    group = models.ProjectGroup.objects.create(
        name=name, type=models.ProjectGroup.GroupType.FRAGMENT_SET)
    if with_parent:
        parent = models.Project.objects.create(
            name=f"{name}_parent", directory=str(root / f"{name}_parent"))
        models.ProjectGroupMembership.objects.create(
            group=group, project=parent, type=PARENT)
    for i in range(n_members):
        project = models.Project.objects.create(
            name=f"{name}_x{i:04d}", directory=str(root / f"{name}_x{i:04d}"))
        models.ProjectGroupMembership.objects.create(
            group=group, project=project, type=MEMBER)
    return group


def _list():
    response = APIClient().get("/api/ccp4i2/projectgroups/")
    assert response.status_code == 200, response.content
    body = response.json()
    rows = body.get("results", body) if isinstance(body, dict) else body
    return {row["name"]: row for row in rows}


def test_member_count_counts_members_and_not_the_parent(
        bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    _campaign(test_project_path, "three", 3)
    _campaign(test_project_path, "none", 0)

    rows = _list()
    # The parent is a membership too, and counting it would overstate every
    # campaign by one.
    assert rows["three"]["member_count"] == 3
    assert rows["none"]["member_count"] == 0


def test_the_count_survives_a_type_filter(bypass_api_permissions, test_project_path):
    test_project_path.mkdir(parents=True, exist_ok=True)
    _campaign(test_project_path, "frags", 2)
    response = APIClient().get("/api/ccp4i2/projectgroups/?type=fragment_set")
    body = response.json()
    rows = body.get("results", body) if isinstance(body, dict) else body
    assert [r["member_count"] for r in rows if r["name"] == "frags"] == [2]


def test_listing_more_campaigns_does_not_cost_more_queries(
        bypass_api_permissions, test_project_path):
    """The point of the annotation. One campaign and then five, and the query
    count must be the same -- a per-row count would make it grow by four.

    Measured rather than asserted against a literal: the absolute number is
    an implementation detail and would break on any unrelated change, while
    its *growth* is the defect this guards against.
    """
    test_project_path.mkdir(parents=True, exist_ok=True)
    _campaign(test_project_path, "solo", 2)

    with CaptureQueriesContext(connection) as first:
        _list()
    baseline = len(first.captured_queries)

    for i in range(4):
        _campaign(test_project_path, f"more{i}", 3)

    with CaptureQueriesContext(connection) as second:
        rows = _list()
    assert len(second.captured_queries) == baseline, (
        f"{baseline} queries for 1 campaign, "
        f"{len(second.captured_queries)} for 5")
    assert len(rows) == 5
    assert rows["more0"]["member_count"] == 3
