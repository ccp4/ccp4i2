"""i2run can build its arguments for tasks with program column groups.

A CProgramColumnGroup answers any attribute name, with None when it names no
column, so hasattr(node, "get_merged_metadata") was true and calling it raised
TypeError: i2run could not so much as list the arguments of ctruncate (or of
cad_copy_column and mtzutils), let alone run one.
"""
import pytest

from ccp4i2.cli.i2run.i2run_components import KeywordExtractor


@pytest.mark.parametrize("task", ["ctruncate", "cad_copy_column", "mtzutils"])
def test_arguments_of_a_task_with_column_groups(task):
    keywords = KeywordExtractor.extract_from_task_name(task)
    assert keywords, f"no arguments for {task}"
