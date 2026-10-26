"""Shared pytest setup for the RMG input converter tests."""

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

DATA = os.path.join(HERE, "data")
PSEUDO = os.path.join(DATA, "pseudo")


@pytest.fixture
def data_dir():
    return DATA


def read_data(*parts):
    with open(os.path.join(DATA, *parts)) as f:
        return f.read()
