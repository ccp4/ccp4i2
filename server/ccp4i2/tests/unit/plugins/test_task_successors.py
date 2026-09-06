"""A superseded task names a registered successor, and the classic Phaser
tasks are superseded by their PHIL counterparts. No CCP4 needed."""
from ccp4i2.core.tasks import TASKS, get_successor


def test_every_successor_is_a_registered_task():
    for name, task in TASKS.items():
        if task.successor:
            assert task.successor in TASKS, f"{name} names a successor that is not registered"
            assert task.successor != name


def test_a_successor_is_not_itself_superseded():
    for name, task in TASKS.items():
        if task.successor:
            assert not TASKS[task.successor].successor, f"{name}'s successor is itself superseded"


def test_classic_phaser_tasks_are_superseded():
    expected = {
        "phaser_MR_AUTO": "phaser_mr_auto_phil", "phaser_MR_RNP": "phaser_mr_rnp_phil",
        "phaser_MR_FRF": "phaser_mr_frf_phil", "phaser_MR_FTF": "phaser_mr_ftf_phil",
        "phaser_MR_PAK": "phaser_mr_pak_phil", "phaser_EP_AUTO": "phaser_ep_auto_phil",
        "phaser_simple": "phaser_simple_phil", "phaser_pipeline": "phaser_pipeline_phil",
        "phaser_EP": "phaser_ep_phil", "phaser_rnp_pipeline": "phaser_rnp_pipeline_phil",
    }
    for classic, successor in expected.items():
        assert get_successor(classic) == successor
    assert get_successor("phaser_mr_auto_phil") is None
    assert get_successor("no_such_task") is None
