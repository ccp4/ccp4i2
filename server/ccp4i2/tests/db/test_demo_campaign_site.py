"""The demo campaign's site sits on the reference fragment, not its mirror.

make_demo_campaign places the site at the centroid of the reference
structure's bound fragment. CampaignSite.origin is Moorhen's view origin, the
NEGATED centre, and .position (what events are matched against) negates it
back; the command stored the centroid itself as the origin, so the site landed
at the mirror image through zero, 50 A from the pocket, and none of the
campaign's PanDDA events ever showed in its column.
"""
from django.test import TestCase

from ...db import models
from ...db.management.commands.make_demo_campaign import create_site_at


class DemoCampaignSiteTest(TestCase):
    def test_the_site_position_is_the_centre_it_was_given(self):
        group = models.ProjectGroup.objects.create(name="demo")
        site = create_site_at(group, "Acetyl-lysine pocket", [-24.2, -7.8, 3.6])
        assert site.position == [-24.2, -7.8, 3.6]
        assert site.origin == [24.2, 7.8, -3.6]
