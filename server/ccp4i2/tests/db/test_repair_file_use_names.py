"""The repair leaves a row alone whenever it cannot be certain.

`manage.py repair_file_use_names` gives a FileUse row that records a bare list
index (`"[0]"`) the list's name back. The valuable half is what it REFUSES to
do: a wrong provenance label is worse than a bare one, and old databases are
where the odd cases live. A healthy database exercises none of that, so these
seed the damage deliberately.
"""
from io import StringIO
from pathlib import Path
from shutil import rmtree

from django.core.management import call_command
from django.test import TestCase, override_settings

from ...db.models import File, FileType, FileUse, Job, Project

PROJECTS_DIR = Path(__file__).parent.parent / "CCP4I2_REPAIR_FILEUSE_TEST_DIR"
PDB = "chemical/x-pdb"
MAP = "application/CCP4-mtz-map"


def run(**kwargs):
    out = StringIO()
    call_command("repair_file_use_names", stdout=out, **kwargs)
    return out.getvalue()


@override_settings(CCP4I2_PROJECTS_DIR=PROJECTS_DIR)
class RepairBase(TestCase):
    """A coot1 job, because its declarations are the awkward real case: one pdb
    list in inputData, a pdb list in outputData, and three map lists of which
    two share a mimeType and differ only by requiredSubType."""

    task = "coot1"

    def setUp(self):
        # Take the directory away again, as the sibling snapshot tests do.
        # Without this the tier leaves CCP4I2_REPAIR_FILEUSE_TEST_DIR behind in
        # the tree, untracked, for everyone who runs it.
        self.addCleanup(rmtree, PROJECTS_DIR, ignore_errors=True)
        self.project_dir = PROJECTS_DIR / "repair"
        self.job_dir = self.project_dir / "CCP4_JOBS" / "job_1"
        self.job_dir.mkdir(parents=True, exist_ok=True)
        self.project = Project.objects.create(
            name="repair", directory=str(self.project_dir)
        )
        self.job = Job.objects.create(
            project=self.project, number="1", task_name=self.task, status=6
        )

    def make_file(self, name, mime, sub_type, on_disk=True):
        file_type, _ = FileType.objects.get_or_create(name=mime)
        if on_disk:
            (self.job_dir / name).write_text("x")
        return File.objects.create(
            name=name,
            directory=File.Directory.JOB_DIR,
            type=file_type,
            sub_type=sub_type,
            job=self.job,
            job_param_name=name,
        )

    def use(self, file_obj, index=0, role=FileUse.Role.IN):
        return FileUse.objects.create(
            file=file_obj, job=self.job, role=role, job_param_name=f"[{index}]"
        )

    def write_params(self, mapping):
        """A params.xml putting each dbFileId inside a named list."""
        blocks = "".join(
            f"<{name}>"
            + "".join(f"<x><dbFileId>{i}</dbFileId></x>" for i in ids)
            + f"</{name}>"
            for name, ids in mapping.items()
        )
        (self.job_dir / "params.xml").write_text(
            f"<{self.task}><container><inputData>{blocks}"
            f"</inputData></container></{self.task}>"
        )

    @staticmethod
    def dbid(file_obj):
        return str(file_obj.uuid).replace("-", "")


class TestItRepairsWhatItCan(RepairBase):
    def test_a_map_is_named_by_its_required_sub_type(self):
        """Two of coot1's map lists share a mimeType; only requiredSubType
        tells FPHIIN_LIST from DELFPHIIN_LIST."""
        fphi = self.make_file("MAPOUT_1.mtz", MAP, 1)
        row = self.use(fphi)
        self.write_params({"FPHIIN_LIST": [self.dbid(fphi)]})

        run(apply=True)

        row.refresh_from_db()
        self.assertEqual(row.job_param_name, "FPHIIN_LIST[0]")

    def test_the_index_is_preserved(self):
        fphi = self.make_file("MAPOUT_2.mtz", MAP, 1)
        row = self.use(fphi, index=3)
        self.write_params({"FPHIIN_LIST": [self.dbid(fphi)]})

        run(apply=True)

        row.refresh_from_db()
        self.assertEqual(row.job_param_name, "FPHIIN_LIST[3]")

    def test_a_dry_run_changes_nothing(self):
        fphi = self.make_file("MAPOUT_3.mtz", MAP, 1)
        row = self.use(fphi)
        self.write_params({"FPHIIN_LIST": [self.dbid(fphi)]})

        output = run()

        row.refresh_from_db()
        self.assertEqual(row.job_param_name, "[0]")
        self.assertIn("Would repair: 1", output)

    def test_running_twice_changes_nothing_the_second_time(self):
        fphi = self.make_file("MAPOUT_4.mtz", MAP, 1)
        self.use(fphi)
        self.write_params({"FPHIIN_LIST": [self.dbid(fphi)]})

        run(apply=True)
        second = run(apply=True)

        self.assertIn("Nothing to do", second)

    def test_the_role_picks_the_container_half(self):
        """A pdb file is ambiguous between coot1's XYZIN_LIST and its XYZOUT
        list until the row's role says which half to look in."""
        pdb = self.make_file("out.pdb", PDB, 0)
        row = self.use(pdb, role=FileUse.Role.OUT)
        self.write_params({"XYZOUT": [self.dbid(pdb)]})

        run(apply=True)

        row.refresh_from_db()
        self.assertEqual(row.job_param_name, "XYZOUT[0]")


class TestItRefusesWhenUnsure(RepairBase):
    def _unchanged(self, row, output, reason):
        row.refresh_from_db()
        self.assertEqual(row.job_param_name, "[0]")
        self.assertIn("Would repair: 0", output)
        self.assertIn(reason, output)

    def test_a_file_missing_from_disk_is_left_alone(self):
        """Provenance that cannot be confirmed is not relabelled."""
        fphi = self.make_file("gone.mtz", MAP, 1, on_disk=False)
        row = self.use(fphi)
        self.write_params({"FPHIIN_LIST": [self.dbid(fphi)]})

        self._unchanged(row, run(), "missing from disk")

    def test_a_params_file_that_disagrees_wins(self):
        """The mapping says FPHIIN_LIST; params.xml says the file is in
        DELFPHIIN_LIST. Contradiction means hands off, not a guess."""
        fphi = self.make_file("odd.mtz", MAP, 1)
        row = self.use(fphi)
        self.write_params({"DELFPHIIN_LIST": [self.dbid(fphi)]})

        self._unchanged(row, run(), "params.xml puts this file in")

    def test_a_file_absent_from_params_is_left_alone(self):
        other = self.make_file("elsewhere.mtz", MAP, 1)
        row = self.use(other)
        self.write_params({"FPHIIN_LIST": ["0" * 32]})

        self._unchanged(row, run(), "does not list this file")

    def test_a_file_with_no_sub_type_is_left_alone(self):
        """File.sub_type is nullable (File.type is not), and old data may never
        have recorded one. With no sub_type no map list can be chosen, because
        they differ only by requiredSubType."""
        untyped = self.make_file("nosub.mtz", MAP, None)
        row = self.use(untyped)
        self.write_params({"FPHIIN_LIST": [self.dbid(untyped)]})

        self._unchanged(row, run(), "candidate lists")

    def test_an_ambiguous_sub_type_is_left_alone(self):
        """sub_type 9 matches no map list's requiredSubType, so there is no
        candidate at all -- reported, not forced."""
        odd = self.make_file("weird.mtz", MAP, 9)
        row = self.use(odd)
        self.write_params({"FPHIIN_LIST": [self.dbid(odd)]})

        self._unchanged(row, run(), "candidate lists")

    def test_an_unreadable_task_is_left_alone(self):
        """A task whose plugin will not import -- removed, or needing CCP4 that
        is not installed -- must not abort the run."""
        Job.objects.filter(pk=self.job.pk).update(task_name="no_such_task_at_all")
        fphi = self.make_file("MAPOUT_5.mtz", MAP, 1)
        row = self.use(fphi)

        self._unchanged(row, run(), "cannot read the declarations")


class TestItLeavesLegacyShapesAlone(RepairBase):
    def test_a_plain_list_name_is_not_touched(self):
        """Qt-era ccp4i2 had no separate name for a CList member, so imported
        rows carry the plain list name repeated. They resolve correctly already
        (name plus an index over the ordered files), and inventing indices for
        them would be worse than leaving them."""
        first = self.make_file("a.pdb", PDB, 0)
        second = self.make_file("b.pdb", PDB, 0)
        for file_obj in (first, second):
            FileUse.objects.create(
                file=file_obj,
                job=self.job,
                role=FileUse.Role.IN,
                job_param_name="XYZIN_LIST",
            )

        output = run(apply=True)

        self.assertIn("Nothing to do", output)
        self.assertEqual(
            list(
                FileUse.objects.filter(job=self.job)
                .values_list("job_param_name", flat=True)
                .distinct()
            ),
            ["XYZIN_LIST"],
        )

    def test_an_already_qualified_name_is_not_touched(self):
        fphi = self.make_file("done.mtz", MAP, 1)
        FileUse.objects.create(
            file=fphi,
            job=self.job,
            role=FileUse.Role.IN,
            job_param_name="FPHIIN_LIST[0]",
        )

        self.assertIn("Nothing to do", run(apply=True))
