"""usable_gib is the cgroup limit inside a container, the host total outside."""
from ccp4i2.lib.utils import system_memory as sm


def test_cgroup_v2_limit_wins_over_the_host(tmp_path, monkeypatch):
    v2 = tmp_path / "memory.max"; v2.write_text("17179869184\n")        # 16 GiB
    monkeypatch.setattr(sm, "host_total_bytes", lambda: 256 * 2 ** 30)
    assert sm.cgroup_limit_bytes([v2]) == 16 * 2 ** 30
    assert sm.usable_gib([v2]) == 16.0


def test_unlimited_or_absent_cgroup_falls_back_to_the_host(tmp_path, monkeypatch):
    monkeypatch.setattr(sm, "host_total_bytes", lambda: 64 * 2 ** 30)
    v2 = tmp_path / "memory.max"; v2.write_text("max\n")
    assert sm.cgroup_limit_bytes([v2]) is None and sm.usable_gib([v2]) == 64.0
    v1 = tmp_path / "memory.limit_in_bytes"; v1.write_text("9223372036854771712\n")   # v1 "unlimited"
    assert sm.cgroup_limit_bytes([v1]) is None
    assert sm.usable_gib([tmp_path / "missing"]) == 64.0
    monkeypatch.setattr(sm, "host_total_bytes", lambda: None)
    assert sm.usable_gib([tmp_path / "missing"]) is None


def test_v1_limit_is_read_when_v2_is_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(sm, "host_total_bytes", lambda: 256 * 2 ** 30)
    v1 = tmp_path / "memory.limit_in_bytes"; v1.write_text("34359738368")              # 32 GiB
    assert sm.usable_gib([tmp_path / "memory.max", v1]) == 32.0
