"""Tests for the ProjectTag forest.

Tags are the mechanism for organising the project list, and once a nested
browser reads them the tree has to hold up: no duplicate nodes, no cycles,
paths that survive a rename or a re-parent, and roll-up that never writes
implied tags back into the database.
"""

from django.core.exceptions import ValidationError
from django.db.utils import IntegrityError
from django.test import TestCase

from ...db.models import Project, ProjectTag

SEP = ProjectTag.PATH_SEPARATOR


class ProjectTagPathTests(TestCase):
    def test_root_path_is_its_text(self):
        tag = ProjectTag.objects.create(text="SARS-CoV-2")
        self.assertEqual(tag.path, "SARS-CoV-2")
        self.assertEqual(tag.depth, 0)

    def test_child_path_includes_ancestry(self):
        root = ProjectTag.objects.create(text="SARS-CoV-2")
        child = ProjectTag.objects.create(text="Mpro", parent=root)
        grandchild = ProjectTag.objects.create(text="soaks", parent=child)

        self.assertEqual(child.path, f"SARS-CoV-2{SEP}Mpro")
        self.assertEqual(grandchild.path, f"SARS-CoV-2{SEP}Mpro{SEP}soaks")
        self.assertEqual(grandchild.display_path, "SARS-CoV-2/Mpro/soaks")
        self.assertEqual(grandchild.depth, 2)

    def test_duplicate_root_text_is_rejected(self):
        """The hole the materialised path exists to close."""
        ProjectTag.objects.create(text="SARS")
        with self.assertRaises(IntegrityError):
            ProjectTag.objects.create(text="SARS")

    def test_same_text_under_different_parents_is_fine(self):
        a = ProjectTag.objects.create(text="2024")
        b = ProjectTag.objects.create(text="2025")
        ProjectTag.objects.create(text="soaks", parent=a)
        ProjectTag.objects.create(text="soaks", parent=b)
        self.assertEqual(ProjectTag.objects.filter(text="soaks").count(), 2)

    def test_slash_in_text_does_not_collide_with_nesting(self):
        """Why the separator is not "/": these are distinct nodes."""
        flat = ProjectTag.objects.create(text="A/B")
        a = ProjectTag.objects.create(text="A")
        nested = ProjectTag.objects.create(text="B", parent=a)
        self.assertNotEqual(flat.path, nested.path)


class ProjectTagRewriteTests(TestCase):
    def setUp(self):
        self.root = ProjectTag.objects.create(text="root")
        self.child = ProjectTag.objects.create(text="child", parent=self.root)
        self.grandchild = ProjectTag.objects.create(text="grandchild", parent=self.child)

    def test_rename_rewrites_descendants(self):
        self.root.text = "renamed"
        self.root.save()

        self.child.refresh_from_db()
        self.grandchild.refresh_from_db()
        self.assertEqual(self.child.path, f"renamed{SEP}child")
        self.assertEqual(self.grandchild.path, f"renamed{SEP}child{SEP}grandchild")

    def test_reparent_rewrites_descendants(self):
        other = ProjectTag.objects.create(text="other")
        self.child.parent = other
        self.child.save()

        self.child.refresh_from_db()
        self.grandchild.refresh_from_db()
        self.assertEqual(self.child.path, f"other{SEP}child")
        self.assertEqual(self.grandchild.path, f"other{SEP}child{SEP}grandchild")

    def test_promote_to_root_rewrites_descendants(self):
        self.child.parent = None
        self.child.save()

        self.grandchild.refresh_from_db()
        self.assertEqual(self.child.path, "child")
        self.assertEqual(self.grandchild.path, f"child{SEP}grandchild")


class ProjectTagCycleTests(TestCase):
    def test_tag_cannot_be_its_own_parent(self):
        tag = ProjectTag.objects.create(text="loop")
        tag.parent = tag
        with self.assertRaises(ValidationError):
            tag.save()

    def test_tag_cannot_be_reparented_under_its_own_descendant(self):
        root = ProjectTag.objects.create(text="root")
        child = ProjectTag.objects.create(text="child", parent=root)
        root.parent = child
        with self.assertRaises(ValidationError):
            root.save()


class ProjectTagQueryTests(TestCase):
    def setUp(self):
        self.root = ProjectTag.objects.create(text="SARS-CoV-2")
        self.mpro = ProjectTag.objects.create(text="Mpro", parent=self.root)
        self.soaks = ProjectTag.objects.create(text="soaks", parent=self.mpro)
        self.unrelated = ProjectTag.objects.create(text="SARS-CoV-2 archive")

        self.top = Project.objects.create(name="top", directory="/tmp/top")
        self.deep = Project.objects.create(name="deep", directory="/tmp/deep")
        self.root.projects.add(self.top)
        self.soaks.projects.add(self.deep)

    def test_descendants_excludes_self_and_prefix_lookalikes(self):
        descendants = set(self.root.descendants())
        self.assertEqual(descendants, {self.mpro, self.soaks})
        self.assertNotIn(self.unrelated, descendants)

    def test_self_and_descendants_includes_self(self):
        self.assertEqual(
            set(self.root.self_and_descendants()), {self.root, self.mpro, self.soaks}
        )

    def test_roll_up_collects_projects_from_below(self):
        self.assertEqual(
            set(self.root.tagged_projects()), {self.top, self.deep}
        )

    def test_roll_up_does_not_write_implied_tags(self):
        self.root.tagged_projects()
        self.assertEqual(set(self.root.projects.all()), {self.top})
        self.assertEqual(set(self.deep.tags.all()), {self.soaks})

    def test_direct_only_ignores_descendants(self):
        self.assertEqual(
            set(self.root.tagged_projects(include_descendants=False)), {self.top}
        )
