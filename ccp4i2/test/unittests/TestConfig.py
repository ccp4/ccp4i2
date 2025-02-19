from ccp4i2.core.CCP4Config import CConfig, CONFIG


def test_default():
    CConfig.insts = None
    CONFIG()
    assert CONFIG().developer
    assert not CONFIG().graphical
    assert CONFIG().dbFile is None
    assert CONFIG().dbUser is None
    assert CONFIG().maxRunningProcesses == 4


def test_custom():
    CConfig.insts = None
    CONFIG(developer=False, graphical=True)
    assert not CONFIG().developer
    assert CONFIG().graphical
