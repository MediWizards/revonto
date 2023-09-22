import pytest

from revonto.reverse_lookup import GOReverseLookupStudy


def test_reverse_lookup_study(annotations_test, godag_test):
    studyset = ["GO:0000002", "GO:0005829"]

    study = GOReverseLookupStudy(annotations_test, godag_test)

    results = study.run_study(studyset)

    assert results[0].object_id == "UniProtKB:A0A024RBG1"
    assert pytest.approx(results[0].pvals["uncorrected"]) == 0.14285714285714288
    assert (
        pytest.approx(results[0].pvals["bonferroni"]) == 0.14285714285714288
    )  # only one test was done
