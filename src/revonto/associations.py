"""
Read and store Gene Ontology's GAF (GO Annotation File).
"""
from __future__ import annotations as an
from typing import Set, Generator, TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:
    from .ontology import GODag
import os
import copy


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


class Annotations(dict[str, Set[Annotation]]):
    """
    Store Annotations as a Dict with key "term_id" and value a set of all Annotation objects connected to that go term.
    This is how the classes preforming the study expect the data to be formated.
    """

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
            self.setdefault(rec.term_id, set()).add(rec)

        return reader.version, reader.date

    def __setitem__(self, key, value):
        if not isinstance(value, set):
            raise ValueError(f"Value for key {key} must be a set of Annotation objects.")
        if not all(isinstance(annotation, Annotation) for annotation in value):
            raise ValueError(f"All elements in the set for key {key} must be Annotation objects.")
        if not all(annotation.term_id == key for annotation in value):
            raise ValueError(f"All Annotation objects must have the same term_id as the key ({key}).")
        super().__setitem__(key, value)

    def __add__(self, other):
        """
        Combine two Annotations objects using the + operator.
        :param other: Another Annotations object to be merged with this one.
        :return: A new Annotations object containing the combined data.
        """
        combined_annotations = Annotations()

        # Merge the current object into the new one
        for term_id, annotations in self.items():
            combined_annotations.setdefault(term_id, set()).update(annotations)

        # Merge the other object into the new one
        for term_id, annotations in other.items():
            combined_annotations.setdefault(term_id, set()).update(annotations)

        return combined_annotations


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
def propagate_associations(godag: GODag, anno: Annotations):
    """
    Iterate through the ontology and assign all childrens' annotations to each term.
    """

    for term_id, term in godag.items():
        annotations_to_append = anno.get(term_id, {})
        for parent in term.get_all_parents():
            for entry in annotations_to_append:
                entry_to_append = (
                    entry.copy()
                )  # make a copy, since we need to change the term_id
                entry_to_append.term_id = parent
                anno.setdefault(parent, set()).add(entry_to_append)


def anno2objkey(anno: Annotations) -> Dict[str, Set[Annotation]]:
    """Change Annotations dict to have object_id as keys"""
    # Should it be moved to Annotations class?
    new_anno = {}
    for goid, goassocset in anno.items():
        for assoc in goassocset:
            new_anno.setdefault(assoc.object_id, set()).add(assoc)
    return new_anno
