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
        log = (job / "log.txt").read_text()
        assert re.search(r"spacegroup: (.*)\n", log).group(1) == "P 41 21 2"
        assert float(re.search(r"High reso[^\d]+([\d\.]+)", log).group(1)) < 1.3
        assert float(re.search(r"Low reso[^\d]+([\d\.]+)", log).group(1)) > 35
        assert int(re.search(r"Total unique +(\d+)", log).group(1)) > 11000


@mark.skip(reason="Skipping temporarily for comprehensive test")
def test_xia2_dials(image_dir):
    run_test("xia2_dials", image_dir)


@mark.skip(reason="Skipping temporarily for comprehensive test")
def test_xia2_xds(image_dir):
    run_test("xia2_xds", image_dir)
