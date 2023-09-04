"""
Read and store Gene Ontology's GAF (GO Annotation File).
"""
from typing import Set

class Annotation():
    """
    Each annotation holds the following variables:
    (GO) term_id
    relationship (beaware of NOT)
    reference
    evidence_code (object)
    taxon
    date
    """

class Annotations(dict[str, Set[Annotation]]):
    """
    Store Annotations as a Dict with key "ID" and value a set of all Annotation objects connected to that product.
    This is how the classes preforming the study expect the data to be formated.

    ID can be anythig (genename, DB:ID) depending on your pre and post analysis.
    """

class AnnoParserBase():
    """
    There is more than one type of annotation file. 
    Therefore we will use a base class to standardize the data and the methods.
    
    Currently we only support GAF, beacuse we need 
    """
    
    # Expected values for a Qualifier
    exp_qualifiers = set([
        # Seen in both GAF and gene2go
        'not', 'contributes_to', 'colocalizes_with',
    ])

    exp_nss = set(['BP', 'MF', 'CC'])

    def __init__(self) -> None:
        pass


class GafParser(AnnoParserBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""
