from ...lib.utils import puid


def test_puid():
    assert len(puid(length=10)) == 10
    assert len(puid(length=20)) == 20
    assert isinstance(puid(), str)
