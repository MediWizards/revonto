import pytest

from revonto.pvalcalc import pvalue_calculate


def test_available_pvalue_calculate():
    assert (
        pytest.approx(pvalue_calculate(1, 2, 3, 40, "fisher_scipy_stats"))
        == 0.1461538461538462
    )
    assert (
        pytest.approx(pvalue_calculate(1, 2, 3, 40, "binomial_scipy_stats"))
        == 0.1443750
    )


def test_pvalue_calculate_exception():
    with pytest.raises(ValueError):
        pvalue_calculate(1, 1, 1, 1, "notamethod")
