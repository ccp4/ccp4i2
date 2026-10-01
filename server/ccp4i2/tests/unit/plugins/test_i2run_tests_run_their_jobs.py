"""Every i2run test runs its job.

i2run() is a context manager: the job runs on entering it. test_ample called
it as a plain statement, so it built the manager and never entered it, and
"passed" in a tenth of a second for as long as it has existed without AMPLE
ever running.
"""
import ast
from pathlib import Path

TESTS = Path(__file__).resolve().parents[2] / "i2run"


def bare_calls(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return [node.lineno for node in ast.walk(tree)
            if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call)
            and getattr(node.value.func, "id", None) == "i2run"]


def test_no_i2run_call_is_left_unentered():
    found = {p.name: bare_calls(p) for p in sorted(TESTS.glob("test_*.py"))}
    assert {k: v for k, v in found.items() if v} == {}, \
        "i2run(...) as a statement never runs the job: use 'with i2run(args) as job:'"
