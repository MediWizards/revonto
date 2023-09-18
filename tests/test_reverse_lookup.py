import pytest

from revonto.reverse_lookup import GOReverseLookupStudy


def test_reverse_lookup_study(annotations_test, godag_test):
    studyset = ["GO:0000002", "GO:0005829"]

    study = GOReverseLookupStudy(annotations_test, godag_test)

    results = study.run_study(studyset)

    assert results[0].object_id == "UniProtKB:A0A024RBG1"
    assert pytest.approx(results[0].p_uncorrected) == 0.40
    assert pytest.approx(results[0].p_bonferroni) == 0.40  # only one test was done
