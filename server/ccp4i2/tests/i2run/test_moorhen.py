"""i2run moorhen: the plugin waits on the session row until a window
finishes it. Here a thread plays the window: it waits for the session,
drops a model through the session library, and finishes."""

import threading
import time

from django.db import connection

from .utils import demoData, i2run


def _play_the_window(saved, model_path):
    from ccp4i2.db import models
    from ccp4i2.lib.utils.jobs import interactive

    try:
        deadline = time.time() + 120
        session = None
        while time.time() < deadline:
            session = models.JobInteractiveSession.objects.filter(
                job__task_name="moorhen", dispatched=True).first()
            if session is not None:
                break
            time.sleep(0.5)
        assert session is not None, "the waiting plugin never opened a session"
        job = session.job
        saved.append(interactive.drop_file(job, model_path, annotation="from the window"))
        saved.append(interactive.finish_session(job)["disposition"])
    finally:
        connection.close()


def test_moorhen_waits_for_the_window_and_harvests_its_save():
    model = demoData("gamma", "gamma_model.pdb")
    saved = []
    window = threading.Thread(target=_play_the_window, args=(saved, model), daemon=True)
    window.start()

    with i2run(["moorhen", "--XYZIN_LIST", model]) as job:
        window.join(timeout=10)
        assert saved and saved[-1] == "harvesting"
        assert (job / "XYZOUT_0.pdb").exists()
        assert (job / "program.xml").read_text().count("<number_output_files>1<") == 1
