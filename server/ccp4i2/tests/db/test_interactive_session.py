"""The interactive session: open, drop, finish.

An interactive task's "program" is a window in the app, so its job is driven
by a session row rather than a child process (docs/moorhen-task-design.md).
This covers that mechanism directly -- open a session, drop files into it,
finish it -- against the session row and the drop directory.

It replaces ``tests/i2run/test_moorhen.py``, which drove the real moorhen
plugin and impersonated the window with a ``daemon=True`` thread. That test
could not fail cleanly: the thread gave up after 120 s and died silently,
while the plugin waits on the session row with no timeout -- correct for a
window a human closes when ready, fatal in a harness. A lost race hung the
whole i2run suite indefinitely instead of failing (1.63 s alone; an infinite
hang 89 job-creating tests deep). The i2run tier should not exercise
interactive tools at all.

What is gained by dropping down a level: no thread, no race, no CCP4, and
every *disposition* gets covered rather than just the one the thread happened
to reach. The drop-directory layout asserted here (``output<N>.<ext>`` plus
``output<N>.meta.json``) is the contract a task's harvest reads, so the
handover is pinned even though no plugin runs.
"""

import json
from pathlib import Path

from django.test import TestCase

from ccp4i2.db import models
from ccp4i2.lib.utils.jobs import interactive


class InteractiveSessionTests(TestCase):
    """The session lifecycle, without a window and without a job process."""

    def setUp(self):
        # A real directory: drop_file writes into it.
        self.project_dir = Path(self.mkdtemp())
        self.project = models.Project.objects.create(
            name="interactive_test", directory=str(self.project_dir)
        )
        self.job = models.Job.objects.create(
            project=self.project,
            number="1",
            title="Moorhen",
            task_name="moorhen",
            status=models.Job.Status.PENDING,
        )
        self.job.directory.mkdir(parents=True, exist_ok=True)

        self.model = self.project_dir / "gamma_model.pdb"
        self.model.write_text("CRYST1\nATOM      1  N   ALA A   1\nEND\n")

    def mkdtemp(self):
        import tempfile

        directory = tempfile.mkdtemp(prefix="ccp4i2_interactive_")
        self.addCleanup(self._rmtree, directory)
        return directory

    @staticmethod
    def _rmtree(directory):
        import shutil

        shutil.rmtree(directory, ignore_errors=True)

    def refetched_job(self):
        """The job as a *fresh* instance, the way the next request sees it.

        Django caches the reverse one-to-one, so a Job held across calls keeps
        the JobInteractiveSession it first read. In production the plugin's
        ensure_session_for_runner and the API's finish_session are different
        processes, so each reads the row afresh; re-fetching here models that
        rather than testing Django's instance cache. Without it, finish_session
        sees a stale dispatched=False and dispatches the job for real -- which
        in this test meant actually spawning a job-runner subprocess.
        """
        return models.Job.objects.get(pk=self.job.pk)

    # -- opening ---------------------------------------------------------

    def test_open_session_sets_running_without_a_process(self):
        """Run on an interactive task opens a session instead of dispatching."""
        interactive.open_session(self.job)
        self.job.refresh_from_db()

        self.assertEqual(self.job.status, models.Job.Status.RUNNING)
        self.assertIsNone(self.job.process_id, "no child process should be claimed")
        self.assertTrue(interactive.drop_dir_for(self.job).is_dir())

    def test_open_session_is_idempotent_for_a_reopened_window(self):
        interactive.open_session(self.job)
        interactive.open_session(self.job)  # window reopened

        self.assertEqual(
            models.JobInteractiveSession.objects.filter(job=self.job).count(), 1
        )

    def test_open_session_refuses_a_job_already_running_a_process(self):
        self.job.status = models.Job.Status.QUEUED
        self.job.save()

        with self.assertRaises(interactive.SessionError):
            interactive.open_session(self.job)

    # -- dropping --------------------------------------------------------

    def test_drop_file_writes_the_harvest_contract(self):
        """output<N>.<ext> plus output<N>.meta.json -- what a harvest reads."""
        interactive.open_session(self.job)

        dropped = interactive.drop_file(
            self.job, self.model, kind="model", annotation="from the window"
        )

        drop_dir = interactive.drop_dir_for(self.job)
        target = drop_dir / dropped["name"]
        self.assertTrue(target.is_file())
        self.assertEqual(target.suffix, ".pdb")
        self.assertEqual(dropped["number"], 1)

        meta = json.loads(target.with_suffix(".meta.json").read_text())
        self.assertEqual(meta["kind"], "model")
        self.assertEqual(meta["annotation"], "from the window")
        self.assertEqual(meta["original_name"], "gamma_model.pdb")

    def test_drops_are_numbered_in_save_order(self):
        interactive.open_session(self.job)

        first = interactive.drop_file(self.job, self.model, annotation="first")
        second = interactive.drop_file(self.job, self.model, annotation="second")

        self.assertEqual((first["number"], second["number"]), (1, 2))
        outputs = interactive.session_outputs(self.job)
        self.assertEqual([o["annotation"] for o in outputs], ["first", "second"])

    def test_drop_file_rejects_an_unknown_kind(self):
        interactive.open_session(self.job)

        with self.assertRaises(interactive.SessionError):
            interactive.drop_file(self.job, self.model, kind="banana")

    def test_drop_file_refuses_a_job_with_no_session(self):
        with self.assertRaises(interactive.SessionError):
            interactive.drop_file(self.job, self.model)

    # -- finishing -------------------------------------------------------

    def test_finish_hands_over_to_a_waiting_plugin(self):
        """The case the old i2run test raced a thread to reach.

        ensure_session_for_runner is what the waiting plugin calls to say it
        owns this job; finishing must then let it harvest rather than
        dispatching the job a second time.
        """
        interactive.open_session(self.job)
        interactive.ensure_session_for_runner(self.job)
        interactive.drop_file(self.job, self.model, annotation="from the window")

        result = interactive.finish_session(self.refetched_job())

        self.assertEqual(result["disposition"], "harvesting")
        self.assertTrue(result["session"]["finished"])
        self.assertEqual(len(result["outputs"]), 1)

    def test_finish_with_nothing_saved_marks_the_job_for_deletion(self):
        interactive.open_session(self.job)

        result = interactive.finish_session(self.job)

        self.assertEqual(result["disposition"], "deleted")
        self.job.refresh_from_db()
        self.assertEqual(self.job.status, models.Job.Status.TO_DELETE)

    def test_a_closed_window_keeps_the_session_open_for_reconnect(self):
        interactive.open_session(self.job)
        interactive.drop_file(self.job, self.model)

        result = interactive.finish_session(self.job, finished=False)

        self.assertEqual(result["disposition"], "kept_open")
        self.assertFalse(result["session"]["finished"])

    def test_finishing_twice_is_reported_not_repeated(self):
        interactive.open_session(self.job)
        interactive.ensure_session_for_runner(self.job)
        interactive.drop_file(self.job, self.model)
        interactive.finish_session(self.refetched_job())

        again = interactive.finish_session(self.refetched_job())

        self.assertEqual(again["disposition"], "already_finished")

    def test_drop_after_finish_is_refused(self):
        interactive.open_session(self.job)
        interactive.ensure_session_for_runner(self.job)
        interactive.drop_file(self.job, self.model)
        interactive.finish_session(self.refetched_job())

        with self.assertRaises(interactive.SessionError):
            interactive.drop_file(self.refetched_job(), self.model)

    def test_cancel_ends_the_session_and_harvests_nothing(self):
        interactive.open_session(self.job)
        interactive.drop_file(self.job, self.model)

        interactive.cancel_session(self.job)

        state = interactive.session_state(self.job)
        self.assertTrue(state["session"]["finished"])
        with self.assertRaises(interactive.SessionError):
            interactive.drop_file(self.job, self.model)
