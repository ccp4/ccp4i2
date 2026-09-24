"""The invocation contract as code (design note 4.1-4.5, 6.1, 6.4; assertions
in 15.2): both defensive arguments in the argv, RAY_TMPDIR set and no PanDDA
switch we did not mean to set in the environment, the failure catalogue,
the progress symbol, and the sizing estimate. CCP4-free, pure."""
import pytest

from ccp4i2.wrappers.pandda_campaign.script import pandda_invocation as c


def test_argv_carries_both_defensive_literals():
    argv = c.build_argv("/s/datasets", "/o/pandda2_out", 6)
    assert argv[argv.index("--dataset_range") + 1] == "0-999999999"
    assert argv[argv.index("--ligand_pdb_regex") + 1] == "ligand.pdb"
    assert argv[argv.index("--pdb_regex") + 1] == "final.pdb"
    assert argv[argv.index("--mtz_regex") + 1] == "final.mtz"
    assert argv[argv.index("--ligand_cif_regex") + 1] == "dict.cif"
    assert argv[argv.index("--ligand_dir_regex") + 1] == "compound"
    assert argv[argv.index("--local_cpus") + 1] == "6"
    assert argv[argv.index("--data_dirs") + 1] == "/s/datasets"
    assert argv[argv.index("--min_characterisation_datasets") + 1] == "25", "PanDDA's default, explicit"
    assert c.build_argv("/s", "/o", 1, 4)[-1] == "4"
    assert argv[argv.index("--out_dir") + 1] == "/o/pandda2_out"


def test_environment_sets_scratch_and_strips_stray_pandda_switches():
    base = {"PATH": "/usr/bin", "CCP4": "/ccp4", "PANDDA_USE_CROWTHER": "0",
            "PANDDA_ANYTHING": "1", "RAY_TMPDIR": "/tmp/ray"}
    env = c.build_env(base, "/job/ray_scratch")
    assert env["RAY_TMPDIR"] == "/job/ray_scratch"
    assert not any(k.startswith("PANDDA_") for k in env), "presence-checked switches: a stray 0 enables"
    assert env["PATH"] == "/usr/bin" and env["CCP4"] == "/ccp4"
    assert "PANDDA_USE_CROWTHER" in base, "the caller's environment is untouched"


def test_progress_is_the_last_line_or_unknown():
    assert c.parse_progress("") is None
    assert c.parse_progress("some ray noise\n") is None
    text = "x\nPANDDA_PROGRESS: dataset 1/201\nray\nPANDDA_PROGRESS: dataset 17/201\n"
    assert c.parse_progress(text) == (17, 201)


@pytest.mark.parametrize("text,expected", [
    ("Traceback ...\nMemoryError\n", "oom"),
    ("... Killed: 9", "oom"),
    ("No RFree Flag found!", "free_r_label"),
    ("KeyError: \"block 'comp_LIG'\"", "ligand_block"),
    ("0/201 datasets passed range filter", "dataset_range_zeroed"),
    ("OSError: [Errno 28] No space left on device", "ray_scratch_full"),
    ("sh: refmac5: command not found", "ccp4_missing"),
    ("OSError: validate_socket_filename failed: AF_UNIX path length cannot exceed 103 bytes: /Users/x/...",
     "socket_path_too_long"),
    ("something else entirely", "unclassified_crash"),
])
def test_failure_classification_is_by_pattern(text, expected):
    name, code, prompt = c.classify_failure(text)
    assert name == expected
    assert 210 <= code <= 217 and prompt


def test_every_catalogue_code_is_distinct():
    codes = [code for _n, _p, code, _t in c.FAILURE_CATALOGUE] + [c.UNCLASSIFIED[2]]
    assert len(set(codes)) == len(codes)


def test_cell_classes_and_sizing():
    baz2b = (82.5, 97.0, 58.1, 90, 90, 90)          # small bromodomain
    cdk4 = (58.0, 64.0, 186.0, 90, 90, 90)          # long axis
    assert c.cell_volume_class(baz2b) == "small"
    assert c.cell_volume_class(cdk4) == "large"
    assert c.sizing_hint(3, [baz2b, cdk4]) == {"datasets": 3, "cell_volume_class": "large"}
    assert c.sizing_hint(0, []) == {"datasets": 0, "cell_volume_class": "small"}
    small = c.estimate_peak_gib(60, "small", 4)
    large = c.estimate_peak_gib(120, "large", 6)
    assert 8 < small < 32, small                    # 6.1: a small cell with <=60 fits 16-32 GB
    assert 30 < large < 60, large                   # 6.1: a single large-cell shell is 30-40 GB+
    assert c.estimate_peak_gib(120, "large", 2) < large, "workers multiply it"
    assert c.estimate_peak_gib(300, "large", 6) == large, "comparators cap at 60"


def test_probe_reads_the_ccp4_launcher(tmp_path):
    root = tmp_path / "share" / "mamba"
    site = root / "envs" / "pandda2" / "lib" / "python3.9" / "site-packages"
    (site / "pandda_2_gemmi-0.0.1.dist-info").mkdir(parents=True)   # the CCP4 bundle's name
    (site / "pandda_2_gemmi-0.0.1.dist-info" / "direct_url.json").write_text('{"url": "file:///jenkins/build"}')
    (site / "pandda_gemmi" / "pandda").mkdir(parents=True)
    (site / "pandda_gemmi" / "pandda" / "pandda.py").write_text("print('PANDDA_PROGRESS: dataset', flush=True)\n")
    launcher = tmp_path / "bin" / "pandda2.analyse"
    launcher.parent.mkdir()
    launcher.write_text(f"#!/bin/sh\n\nexec {tmp_path}/micromamba/bin/micromamba run -r {root} -n pandda2 pandda2.analyse \"$@\"\n")
    probe = c.probe_executable(launcher)
    assert probe["version"] == "0.0.1"
    assert probe["distribution"] == "pandda_2_gemmi"
    assert probe["origin"] == "file:///jenkins/build"
    assert probe["progress_signal"] is True
    assert probe["site_packages"] == str(site)
    assert "micromamba run" in probe["launcher"]
    assert c.probe_executable(tmp_path / "missing")["version"] is None


def test_scratch_defaults_to_a_path_ray_can_open_a_socket_under(tmp_path):
    short = "/tmp/j"
    assert c.scratch_fits(short)
    assert c.default_scratch_dir(short, "abcdef0123") == c.Path(short) / "ray_scratch"
    # a project store under a home directory: what a laptop actually has
    long = "/Users/someone/.ccp4i2-pandda/projects/baz2b_demo_campaign_5e9i/CCP4_JOBS/job_3"
    assert not c.scratch_fits(c.Path(long) / "ray_scratch")
    fallback = c.default_scratch_dir(long, "abcdef0123456789")
    assert str(fallback) == "/tmp/ccp4i2-ray-abcdef01"
    assert c.scratch_fits(fallback)


def test_a_hollow_run_is_told_from_a_real_one(tmp_path):
    from .synthetic_tree import event_record, make_tree
    tree = make_tree(tmp_path / "out", {"xtal-0000": [event_record(1)], "xtal-0001": [], "xtal-0002": []})
    (tree / "processed_datasets" / "xtal-0002" / "xtal-0002-z_map.native.ccp4").unlink()
    log = ("    xtal-0002 : Filtered because no ligand data!                     \n"
           "comparators : before : 5\n"
           "    NOT ENOUGH COMPARATOR DATASETS: 5! SKIPPING!                     \n"
           "    NOT ENOUGH COMPARATOR DATASETS: 5! SKIPPING!                     \n")
    summary = c.summarise_output_tree(tree, log)
    assert (summary["processed"], summary["analysed"], summary["events"], summary["complete"]) == (3, 2, 1, True)
    assert summary["unanalysed"] == ["xtal-0002"]
    assert summary["reasons"] == ["xtal-0002: Filtered because no ligand data!",
                                  "NOT ENOUGH COMPARATOR DATASETS: 5! SKIPPING!"]
    for d in (tree / "processed_datasets").iterdir():
        zmap = d / f"{d.name}-z_map.native.ccp4"
        if zmap.exists():
            zmap.unlink()
    assert c.summarise_output_tree(tree, "")["analysed"] == 0
