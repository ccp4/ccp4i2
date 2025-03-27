from tarfile import open as taropen
from tempfile import TemporaryDirectory
from pytest import fixture, mark
from .utils import download, i2run


@fixture(scope="module", name="image_dir")
def image_dir_fixture():
    url = "https://zenodo.org/record/10271/files/th_8_2.tar.bz2"
    with download(url) as tarPath:
        with taropen(tarPath) as tar:
            with TemporaryDirectory() as tmpDir:
                tar.extractall(tmpDir, tar.getmembers()[:20])
                yield tmpDir


@mark.parametrize("task", ["xia2_dials", "xia2_xds"])
def test_xia2(task, image_dir):
    args = [task]
    args += ["--IMAGE_FILE"]
    args += ["imageFile/baseName=th_8_2_0001.cbf"]
    args += ["imageFile/relPath=" + image_dir]
    args += ["imageStart=1"]
    args += ["imageEnd=20"]
    with i2run(args) as job:
        assert False, "Check output files"
