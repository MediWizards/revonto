import pytest

from revonto.ortholog import find_orthologs, gOrth


def test_gOrth():
    result = gOrth(["ZDB-GENE-040912-6", "ZDB-GENE-170217-1"], "drerio", "hsapiens")
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000005421", "ENSG00000105852", "ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
    }


def test_missing_gOrth():
    result = gOrth(["A0A087WV62"], "drerio", "hsapiens")
    assert result == {
        "A0A087WV62": [],
    }


def test_find_orthologs():
    result = find_orthologs(
        ["ZDB-GENE-040912-6", "ZDB-GENE-170217-1", "ZDB-GENE-021119-1"],
        "7955",
        "9606",
        database="gOrth",
    )
    assert result == {
        "ZDB-GENE-040912-6": ["ENSG00000005421", "ENSG00000105852", "ENSG00000105854"],
        "ZDB-GENE-170217-1": ["ENSG00000168938"],
        "ZDB-GENE-021119-1": ["ENSG00000151577"],
    }


@pytest.mark.parametrize("db", ["gOrth"])
def test_multiple_same_incoming(db):
    result = find_orthologs(
        ["ZDB-GENE-170217-1", "ZDB-GENE-170217-1"], "7955", "9606", database=db
    )
    assert len(result) == 1
