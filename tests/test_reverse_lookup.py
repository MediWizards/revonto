import pytest

from revonto.reverse_lookup import GOReverseLookupStudy

import os

def test_reverse_lookup_study():
    from revonto.associations import Annotations
    anno = Annotations(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf"))

    from revonto.ontology import GODag
    godag = GODag(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo"))

    studyset = ["GO:0000002", "GO:0005829"]

    study = GOReverseLookupStudy(anno, godag)

    results = study.run_study(studyset)

    assert results[0].object_id == "UniProtKB:A0A024RBG1"
    assert pytest.approx(results[0].p_uncorrected) == 0.40
