import pytest

import os
import re
from revonto.associations import Annotations, Annotation, propagate_associations
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

def test_enforce_data_structure():
    annodict = Annotations()
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    with pytest.raises(ValueError, match=re.escape("All elements in the set for key GO:1234 must be Annotation objects.")):
        annodict["GO:1234"] = {"A"}
    with pytest.raises(ValueError, match=re.escape("Value for key GO:1234 must be a set of Annotation objects.")):
        annodict["GO:1234"] = anno1
    with pytest.raises(ValueError, match=re.escape("All Annotation objects must have the same term_id as the key (GO:5678).")):
        annodict["GO:5678"] = {anno1}

def test_add_two_anotations_dicts():
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    anno2 = Annotation(object_id="ABC2", term_id="GO:1234")
    anno3 = Annotation(object_id="ABC2", term_id="GO:5678")

    annodict1 = Annotations()
    annodict1["GO:1234"] = {anno1}
    annodict2 = Annotations()
    annodict2["GO:1234"] = {anno2}
    annodict3 = Annotations()
    annodict3["GO:5678"] = {anno3}

    assert annodict1+annodict2 == annodict2+annodict1 # addition should be the same regardles of order

    annodict12 = annodict1 + annodict2
    assert len(annodict12) and "GO:1234" in annodict12 # only one key should exist
    assert len(annodict12["GO:1234"]) == 2 # two Annotation s should be present under this key

    annodict11 = annodict1 + annodict1
    assert len(annodict11) == 1 and len(annodict11["GO:1234"]) == 1 # adding same Annotation to the set should not change it

    annodict13 = annodict1 + annodict3
    assert len(annodict13) == 2 # two keys should be present
    
