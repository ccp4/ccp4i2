import pytest

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from ccp4i2.core.base_object.fundamental_types import CInt, CList

class TestExample:
    @classmethod
    def setup_class(cls):
        # Class-level setup (runs once before all tests)
        cls.shared_resource = "initialized"

    def setup_method(self, method):
        # Method-level setup (runs before each test)
        self.data = [1, 2, 3]

    def test_one(self):
        a = CInt(1)
        assert a == 1
        assert a == CInt(1)
        b = CInt(2)
        assert b == 2
        assert a + b == 3
        assert a * 2 == 2
        assert a - a == 0
        assert a - b == -1
        assert a / b == 0.5
        assert b ** 2 == 4

    def test_arrays(self):
        arr = CList([CInt(i) for i in range(5)], name="arr")
        assert arr.name == "arr"
        assert sum(arr) == CInt(10)
        assert arr[1].parent == arr
