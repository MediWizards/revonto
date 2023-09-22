import pytest

from revonto.associations import (
    Annotation,
    Annotations,
    match_annotations_to_godag,
    propagate_associations,
)
from revonto.ontology import GODag, GOTerm


def test_header(annotations_test):
    assert annotations_test.version == "2.2"
    assert annotations_test.date == "2023-07-29T02:43"


@pytest.mark.skip
def test_UniProtKBA0A024RBG1_assoc(annotations_test):
    assert "GO:0002250" in annotations_test
    assert len(annotations_test["GO:0002250"]) == 2
    assert (
        next(
            obj.relationship
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "involved_in"
    )
    assert (
        next(
            obj.NOTrelation
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        is False
    )
    assert (
        next(
            obj.reference
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "GO_REF:0000043"
    )
    assert (
        next(
            obj.evidence_code
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "IEA"
    )
    assert (
        next(
            obj.taxon
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "taxon:9606"
    )
    assert (
        next(
            obj.date
            for obj in annotations_test["GO:0002250"]
            if obj.object_id == "UniProtKB:A0A075B6H7"
        )
        == "20230703"
    )


@pytest.mark.skip
def test_UniProtKBA0A024RBG1_cardinality_0_fields_assoc(annotations_test):
    assert (
        next(
            obj.relationship
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "located_in"
    )
    assert (
        next(
            obj.NOTrelation
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        is False
    )
    assert (
        next(
            obj.reference
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "GO_REF:0000052"
    )
    assert (
        next(
            obj.evidence_code
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "IDA"
    )
    assert (
        next(
            obj.taxon
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "taxon:9606"
    )
    assert (
        next(
            obj.date
            for obj in annotations_test["GO:0005829"]
            if obj.object_id == "UniProtKB:A0A024RBG1"
        )
        == "20230619"
    )


def test_propagate_associations(annotations_test, godag_test):
    propagated_annotations = propagate_associations(annotations_test, godag_test)
    assert any(anno.term_id == "GO:0000001" for anno in propagated_annotations)
    assert not any(anno.term_id == "GO:0000003" for anno in propagated_annotations)
    assert (
        sum(1 for anno in propagated_annotations if anno.term_id == "GO:0000001") == 2
    )


def test_dict_from_attr():
    annoset = Annotations()
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    anno2 = Annotation(object_id="ABC2", term_id="GO:1234")
    anno3 = Annotation(object_id="ABC2", term_id="GO:5678")
    annoset.update([anno1, anno2, anno3])

    dict_by_term_id = annoset.dict_from_attr("term_id")
    assert len(dict_by_term_id) == 2
    assert len(dict_by_term_id["GO:1234"]) == 2

    dict_by_object_id = annoset.dict_from_attr("object_id")
    assert len(dict_by_object_id) == 2
    assert len(dict_by_object_id["ABC2"]) == 2


def test_annotation_set_operations():
    anno1 = Annotation(object_id="ABC1", term_id="GO:1234")
    anno2 = Annotation(object_id="ABC2", term_id="GO:1234")
    anno3 = Annotation(object_id="ABC2", term_id="GO:5678")

    annoset1 = Annotations()
    annoset1.add(anno1)
    annoset2 = Annotations()
    annoset2.add(anno2)
    annoset3 = Annotations()
    annoset3.add(anno3)

    assert annoset1.union(annoset2) == annoset2.union(
        annoset1
    )  # union should be the same regardles of order

    assert len(annoset1.union(annoset2, annoset3)) == 3

    assert annoset1.intersection(annoset2) == set()

    assert annoset1.intersection(annoset1.union(annoset2)) == annoset1


def test_match_annotations_to_godag():
    annoset = Annotations()
    annoset.update(
        [
            Annotation(object_id="ABC1", term_id="GO:1234"),
            Annotation(object_id="ABC2", term_id="GO:1234"),
            Annotation(object_id="ABC2", term_id="GO:5678"),
        ]
    )

    godag = GODag()
    godag["GO:1234"] = GOTerm("GO:1234")

    matched_anno = match_annotations_to_godag(annoset, godag)

    assert len(annoset) == 3
    assert len(matched_anno) == 2
