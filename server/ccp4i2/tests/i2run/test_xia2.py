from sys import platform
from tarfile import open as taropen
from tempfile import TemporaryDirectory
import re
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


def run_test_from_image_file(task, image_dir):
    args = [task]
    args += ["--IMAGE_FILE"]
    args += ["imageFile/baseName=th_8_2_0001.cbf"]
    args += ["imageFile/relPath=" + image_dir]
    args += ["imageStart=1"]
    args += ["imageEnd=20"]
    run_test(args)


def run_test_from_image_directory(task, image_dir):
    args = [task]
    args += ["--IMAGE_DIRECTORY", str(image_dir)]
    run_test(args)


def run_test(args):
    with i2run(args) as job:
        for name in ("freer", "NATIVE_SWEEP1_INTEGRATE", "obs"):
            read_mtz_file(str(job / f"AUTOMATIC_DEFAULT_{name}.mtz"))
        log = (job / "log.txt").read_text()
        assert re.search(r"spacegroup: (.*)\n", log).group(1) == "P 41 21 2"
        assert float(re.search(r"High reso[^\d]+([\d\.]+)", log).group(1)) < 1.3
        assert float(re.search(r"Low reso[^\d]+([\d\.]+)", log).group(1)) > 35
        assert int(re.search(r"Total unique +(\d+)", log).group(1)) > 11000


@mark.skipif(platform == "win32", reason="Not supported on Windows with Python 3.9 xia2")
def test_xia2_dials_file(image_dir):
    run_test_from_image_file("xia2_dials", image_dir)


def test_xia2_dials_directory(image_dir):
    run_test_from_image_directory("xia2_dials", image_dir)


@mark.skipif(platform == "win32", reason="Not supported on Windows")
def test_xia2_xds_file(image_dir):
    run_test_from_image_file("xia2_xds", image_dir)
