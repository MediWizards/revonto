from collections import defaultdict
from typing import Union

import requests

from .utils import NCBITaxon_to_gProfiler


def gOrth(
    source_ids: list[str], source_taxon: str, target_taxon: str
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

    target_ids = defaultdict(list, {k: [] for k in source_ids})  # initialise with keys
    result: list[dict] = r.json()["result"]
    for entry in result:
        entry_source_id = entry["incoming"]
        if entry["ortholog_ensg"] not in ["N/A", "None", None]:
            target_ids[entry_source_id].append(entry["ortholog_ensg"])
    return target_ids


def find_orthologs(
    source_ids: Union[str, list[str], set[str]],
    source_taxon: str,
    target_taxon: str = "9606",
    database: str = "gOrth",
) -> dict[str, list[str]]:
    """_summary_

    Args:
        source_ids (Union[str, list[str], set[str]]): _description_
        source_taxon (str): _description_
        target_taxon (str): _description_
        database (str): _description_

    Raises:
        NotImplementedError: _description_

    Returns:
        _type_: _description_
    """
    if not isinstance(source_taxon, str) and not isinstance(target_taxon, str):
        raise TypeError("taxons must be str")

    if isinstance(source_ids, str):
        source_ids_list = [source_ids]
    if isinstance(source_ids, set):
        source_ids_list = list(source_ids)
    if isinstance(source_ids, list):
        source_ids_list = source_ids

    if database == "gOrth":
        source_taxon = NCBITaxon_to_gProfiler(source_taxon)
        target_taxon = NCBITaxon_to_gProfiler(target_taxon)
        if not source_taxon or not target_taxon:
            return {}
        target_ids_dict = gOrth(source_ids_list, source_taxon, target_taxon)
    elif database == "local_files":
        raise NotImplementedError
    else:
        ValueError(
            f"database {database} is not available as a source of ortholog information"
        )

    return target_ids_dict
