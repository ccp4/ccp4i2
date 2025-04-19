from sys import platform
from tarfile import open as taropen
from tempfile import TemporaryDirectory
from gemmi import read_mtz_file
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


def run_test(task, image_dir):
    args = [task]
    args += ["--IMAGE_FILE"]
    args += ["imageFile/baseName=th_8_2_0001.cbf"]
    args += ["imageFile/relPath=" + image_dir]
    args += ["imageStart=1"]
    args += ["imageEnd=20"]
    with i2run(args) as job:
        for name in ("freer", "NATIVE_SWEEP1_INTEGRATE", "obs"):
            read_mtz_file(str(job / f"AUTOMATIC_DEFAULT_{name}.mtz"))


def test_xia2_dials(image_dir):
    run_test("xia2_dials", image_dir)


@mark.skipif(platform == "win32", reason="Not supported on Windows")
def test_xia2_xds(image_dir):
    run_test("xia2_xds", image_dir)
