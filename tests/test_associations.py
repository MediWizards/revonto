import pytest

import os
from revonto.associations import Annotations

anno = Annotations(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf"))

def test_header():
    assert anno.version == "2.2"
    assert anno.date == "2023-07-29T02:43"

def test_UniProtKBA0A024RBG1_assoc():
    assert "UniProtKB:A0A024RBG1" in anno
    assert len(anno["UniProtKB:A0A024RBG1"]) == 3
    assert next(obj.relationship for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == "enables"
    assert next(obj.NOTrelation for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == False
    assert next(obj.reference for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == "GO_REF:0000043"
    assert next(obj.evidence_code for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == "IEA"
    assert next(obj.taxon for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == "taxon:9606"
    assert next(obj.date for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0046872") == "20230703"

def test_UniProtKBA0A024RBG1_cardinality_0_fields_assoc():
    assert next(obj.relationship for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == "located_in"
    assert next(obj.NOTrelation for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == False
    assert next(obj.reference for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == "GO_REF:0000052"
    assert next(obj.evidence_code for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == "IDA"
    assert next(obj.taxon for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == "taxon:9606"
    assert next(obj.date for obj in anno["UniProtKB:A0A024RBG1"] if obj.term_id == "GO:0005829") == "20230619"
