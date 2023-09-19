"""
Read and store Gene Ontology's GAF (GO Annotation File).
"""
from __future__ import annotations as an

from typing import TYPE_CHECKING, Any, Generator

if TYPE_CHECKING:
    from .ontology import GODag

import copy
import os


class Annotation:
    """
    Each annotation holds the following variables:
    object_id (unique identifier of the product) - can be genename, DB:ID, ...
    (GO) term_id
    relationship (beaware of NOT)
    reference
    evidence_code (object)
    taxon
    date
    """

    def __init__(
        self,
        object_id=None,
        term_id="",
        relationship=None,
        NOTrelation=False,
        reference=None,
        evidence_code=None,
        taxon=None,
        date=None,
        **kwargs,
    ) -> None:
        # mandatory - this makes an annotation "unique", rest is just metadata
        self.object_id = object_id
        self.term_id = term_id
        # optional but recommended
        self.relationship = relationship
        self.NOTrelation = NOTrelation
        self.reference = reference
        self.evidence_code = evidence_code
        self.taxon = taxon
        self.date = date
        # you can add any number of others TODO: Maybe optional object class like goatools

    def copy(self) -> Annotation:
        return copy.deepcopy(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Annotation):
            return NotImplemented
        if self.object_id == other.object_id and self.term_id == other.term_id:
            return True
        else:
            return False

    def __hash__(self):
        return hash((self.object_id, self.term_id))


class Annotations(set[Annotation]):
    """Store annotations as a set of Annotation objects"""

    def __init__(self, file: str = ""):
        super().__init__()
        if file != "":
            self.version, self.date = self.load_assoc_file(file)

    def load_assoc_file(self, assoc_file):
        """read association file"""

        extension = os.path.splitext(assoc_file)[1]
        if extension == ".gaf":
            reader = GafParser(assoc_file)
        elif extension == ".gpad":
            raise NotImplementedError("GPAD files are not yet supported")
        else:
            raise NotImplementedError(f"{extension} files are not yet supported")

        for rec in reader:
            self.add(rec)

        return reader.version, reader.date

    def dict_from_attr(self, attribute: str) -> dict[str, set[Annotation]]:
        """groups annotations by attribute and store it in dictionary

        Args:
            attribute (str): which Annotation attribute to group by

        Raises:
            ValueError: if attribute is not in Annotation class

        Returns:
            _type_: dictionary of sets of Annotation objects, grouped by attribute
        """
        if not hasattr(Annotation(), attribute):
            raise ValueError(f"Attribute {attribute} not in Annotation class.")

        grouped_dict: dict[str, set[Annotation]] = {}
        for anno in self:
            attribute_value = getattr(anno, attribute)
            if attribute_value == "" or attribute_value is None:
                attribute_value = "None"
            grouped_dict.setdefault(attribute_value, set()).add(anno)

        return grouped_dict


class AnnoParserBase:
    """
    There is more than one type of annotation file.
    Therefore we will use a base class to standardize the data and the methods.

    Currently we only support GAF, beacuse we need
    """

    def __init__(self, assoc_file) -> None:
        if os.path.isfile(assoc_file):
            self.assoc_file = assoc_file
        else:
            raise FileNotFoundError(f"{assoc_file} not found")
        self.version = None
        self.date = None

    def __iter__(self):
        raise NotImplementedError("Call derivative class!")


class GafParser(AnnoParserBase):
    """Reads a Gene Annotation File (GAF). Returns an iterable. One association at a time."""

    def __init__(self, assoc_file) -> None:
        super().__init__(assoc_file)

    def __iter__(self) -> Generator[Annotation, Any, Any]:
        with open(self.assoc_file) as fstream:
            hdr = True

            for line in fstream:
                line = line.rstrip()
                if hdr:
                    if not self._init_hdr(line):
                        hdr = False
                if not hdr and line:
                    values = line.split("\t")
                    rec_curr = Annotation()
                    self._add_to_ref(rec_curr, values)
                    yield rec_curr

    def _init_hdr(self, line: str):
        """save gaf version and date"""
        if line[:14] == "!gaf-version: ":
            self.version = line[14:]
            return True
        if line[:17] == "!date-generated: ":
            self.date = line[17:]
            return True
        if line[0] != "!":
            return False
        return True

    def _add_to_ref(self, rec_curr: Annotation, values):
        """populate Annotation object with values from line"""
        rec_curr.object_id = values[0] + ":" + values[1]
        rec_curr.term_id = values[4]
        rec_curr.relationship = values[3]  # change to object
        if "NOT" in values[3]:
            rec_curr.NOTrelation = True
        rec_curr.reference = values[5]
        rec_curr.evidence_code = values[6]  # change to object
        rec_curr.taxon = values[12]  # change to object
        rec_curr.date = values[13]


class EvidenceCodes:
    """
    class which holds information about evidence codes.
    upon creation the fields are populated accordint to the evicence code in __init__
    currently not used
    """

    codes = {}

    def __init__(self, code) -> None:
        if code not in self.codes:
            pass


# in future maybe move it to update_associations.py
def propagate_associations(anno: Annotations, godag: GODag):
    """
    Iterate through the ontology and assign all childrens' annotations to each term.
    """
    anno_term_dict = anno.dict_from_attr(
        "term_id"
    )  # create a dictionary with annotations grouped by term_id
    propagated_anno = Annotations()
    propagated_anno.update(anno)
    for term_id, term in godag.items():
        annotations_to_append = anno_term_dict.get(term_id, set())
        for parent in term.get_all_parents():
            for entry in annotations_to_append:
                entry_to_append = (
                    entry.copy()
                )  # make a copy, since we need to change the term_id
                entry_to_append.term_id = parent
                # TODO: change evidence code or something to mark the propagated associations
                propagated_anno.add(entry_to_append)
    return propagated_anno


def match_annotations_to_godag(anno: Annotations, godag: GODag):
    """match that all goterms in Annotations are also in GODag.

    Args:
        anno (Annotations): _description_
        godag (GODag): _description_
    """
    all_goterms_in_godag = godag.keys()
    matched_annoobj = Annotations()
    for annoobj in anno:
        if annoobj.term_id in all_goterms_in_godag:
            matched_annoobj.add(annoobj)
    return matched_annoobj
