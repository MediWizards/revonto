import pytest

import os
from revonto.associations import Annotations, propagate_associations
from revonto.ontology import GODag

anno = Annotations(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/human_test.gaf"))

def test_header():
    assert anno.version == "2.2"
    assert anno.date == "2023-07-29T02:43"

def test_UniProtKBA0A024RBG1_assoc():
    assert "GO:0002250" in anno
    assert len(anno["GO:0002250"]) == 2
    assert next(obj.relationship for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == "involved_in"
    assert next(obj.NOTrelation for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == False
    assert next(obj.reference for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == "GO_REF:0000043"
    assert next(obj.evidence_code for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == "IEA"
    assert next(obj.taxon for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == "taxon:9606"
    assert next(obj.date for obj in anno["GO:0002250"] if obj.object_id == "UniProtKB:A0A075B6H7") == "20230703"

def test_UniProtKBA0A024RBG1_cardinality_0_fields_assoc():
    assert next(obj.relationship for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == "located_in"
    assert next(obj.NOTrelation for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == False
    assert next(obj.reference for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == "GO_REF:0000052"
    assert next(obj.evidence_code for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == "IDA"
    assert next(obj.taxon for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == "taxon:9606"
    assert next(obj.date for obj in anno["GO:0005829"] if obj.object_id == "UniProtKB:A0A024RBG1") == "20230619"

godag = GODag(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/go1.obo"))

def test_propagate_associations():
    propagate_associations(godag, anno)
    assert "GO:0000001" in anno
    assert "GO:0000003" not in anno
    assert len(anno["GO:0000001"]) == 2