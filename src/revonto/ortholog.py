from collections import defaultdict
from typing import Union

import requests


def NCBITaxon_to_gProfiler(taxon):
    """_summary_

    Args:
        taxon (_type_): _description_

    Returns:
        _type_: _description_
    """
    taxon_equivalents = {
        9606: "hsapiens",
        7955: "drerio",
    }
    return taxon_equivalents[taxon]


def gOrth(
    source_ids: list, source_taxon: str, target_taxon: str
) -> dict[str, list[str]]:
    """_summary_

    Args:
        source_ids (list): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_
    """
    r = requests.post(
        url="https://biit.cs.ut.ee/gprofiler_archive3/e108_eg55_p17/api/orth/orth/",
        json={
            "organism": source_taxon,
            "target": target_taxon,
            "query": source_ids,
        },
    )

    target_ids = defaultdict(list)
    result: list[dict] = r.json()["result"]
    for entry in result:
        entry_source_id = entry["incoming"]
        if entry["ortholog_ensg"] != "N/A":
            target_ids[entry_source_id].append(entry["ortholog_ensg"])
        else:
            target_ids[entry_source_id] = []
    return target_ids


def find_orthologs(
    source_ids: Union[str, list[str]],
    source_taxon: str,
    target_taxon: str,
    database: str,
) -> dict[str, list[str]]:
    """_summary_

    Args:
        source_ids (Union[str, list[str]]): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_
        database (str): _description_

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: _description_
    """
    if isinstance(source_ids, str):
        source_ids = [source_ids]

    if database == "gOrth":
        source_taxon = NCBITaxon_to_gProfiler(source_taxon)
        target_taxon = NCBITaxon_to_gProfiler(target_taxon)
        target_ids = gOrth(source_ids, source_taxon, target_taxon)
    elif database == "local_files":
        raise NotImplementedError
    else:
        ValueError(
            f"database {database} is not available as a source of ortholog information"
        )

    return target_ids
