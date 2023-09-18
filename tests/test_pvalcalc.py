import pytest

from revonto.pvalcalc import PValueFactory

fisherobj = PValueFactory("fisher_scipy_stats")
binomialobj = PValueFactory("binomial_scipy_stats")


def test_fisher_example():
    assert (
        pytest.approx(fisherobj.pval_obj.calc_pvalue(2, 12, 16, 29), rel=1e-6)
        == 0.0007206922
    )  # round the fisher exact score to a reasonable amount


def test_fisher_exception():
    with pytest.raises(ValueError):
        fisherobj.pval_obj.calc_pvalue(2, 10, 1, 20)  # c will be -1


def test_binomial_example():
    assert (
        pytest.approx(binomialobj.pval_obj.calc_pvalue(2, 12, 16, 29), rel=1e-6)
        == 0.00841718
    )
