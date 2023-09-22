import pytest

from revonto.multiple_testing import CorrectionBase, MultiCorrectionFactory


@pytest.fixture
def bonferroni() -> CorrectionBase:
    return MultiCorrectionFactory("bonferroni").corr_obj


@pytest.fixture
def sm_bonferroni() -> CorrectionBase:
    return MultiCorrectionFactory("sm_bonferroni").corr_obj


def test_bonferroni(bonferroni: CorrectionBase):
    assert bonferroni.set_correction([0.1, 0.2, 0.5, 0.05, 0.01]).tolist() == [
        0.5,
        1,
        1,
        0.25,
        0.05,
    ]


def test_sm_bonferroni(sm_bonferroni: CorrectionBase):
    assert sm_bonferroni.set_correction([0.1, 0.2, 0.5, 0.05, 0.01]).tolist() == [
        0.5,
        1,
        1,
        0.25,
        0.05,
    ]
