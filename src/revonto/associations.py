"""
Read and store Gene Ontology's GAF (GO Annotation File).
"""
from __future__ import annotations as an

from typing import TYPE_CHECKING, Any, Generator, Optional

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
        self.object_id = object_id  # genename is nice, but not unique trans-species! the use of genename as an object_id is highly discouraged. If you use genename, at least add some kind of species identifier to it. You can use add taxon_to_object_id()
        self.term_id = term_id
        self.taxon = (
            taxon  # just in case you decide to use genename, the taxon is needed
        )
        # optional but recommended
        self.relationship = relationship
        self.NOTrelation = NOTrelation
        self.reference = reference
        self.evidence_code = evidence_code
        self.date = date
        # you can add any number of others TODO: Maybe optional object class like goatools

    def copy(self) -> Annotation:
        return copy.deepcopy(self)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Annotation):
            return NotImplemented
        if (
            self.object_id == other.object_id
            and self.term_id == other.term_id
            and self.taxon == other.taxon
        ):
            return True
        else:
            return False

    def __hash__(self):
        return hash((self.object_id, self.term_id, self.taxon))


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
        """groups annotations by attribute and store it in dictionary.
        MOVE TO READ ME: If you use gene name as an identifier, some of the transpecies annotations might be grouped together

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

    def propagate_associations(self, godag: GODag) -> None:
        """
        Iterate through the ontology and assign all childrens' annotations to each term.
        """
        anno_term_dict = self.dict_from_attr(
            "term_id"
        )  # create a dictionary with annotations grouped by term_id

        for term_id, term in godag.items():
            annotations_to_append = anno_term_dict.get(term_id, set())
            for parent in term.get_all_parents():
                for entry in annotations_to_append:
                    entry_to_append = (
                        entry.copy()
                    )  # make a copy, since we need to change the term_id
                    entry_to_append.term_id = parent
                    # TODO: change evidence code or something to mark the propagated associations
                    self.add(entry_to_append)

    def match_annotations_to_godag(self, godag: GODag) -> None:
        """match that all goterms in Annotations are also in GODag.

        Args:
            anno (Annotations): _description_
            godag (GODag): _description_
        """
        all_goterms_in_godag = godag.keys()
        items_to_remove = set()
        for annoobj in self:
            if annoobj.term_id not in all_goterms_in_godag:
                items_to_remove.add(annoobj)
        self.difference_update(items_to_remove)

    def add_taxon_to_object_id(self) -> None:
        """Append taxon to object_id. Especially useful if you decide to use gene names for object_id.

        Args:
            anno (Annotations): _description_
        """

        for annoobj in self:
            if annoobj.taxon:
                annoobj.object_id = annoobj.object_id + "-" + annoobj.taxon

    def find_orthologs(self, taxon, prune=False) -> None:
        """_summary_

        Args:
            taxon (_type_): _description_
            prune (bool, optional): _description_. Defaults to False.
        """        
        all_object_ids = set()
        for anno in self:
            all_object_ids.add(anno.object_id)
        orthologs_dict = 




class AnnoParserBase:
    """
    There is more than one type of annotation file.
    Therefore we will use a base class to standardize the data and the methods.
    """

    def __init__(self, assoc_file) -> None:
        if os.path.isfile(assoc_file):
            self.assoc_file = assoc_file
        else:
            raise FileNotFoundError(f"{assoc_file} not found")
        self.version: Optional[str] = None
        self.date: Optional[str] = None

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
        rec_curr.relationship = values[3]  # TODO:change to object
        if "NOT" in values[3]:
            rec_curr.NOTrelation = True
        rec_curr.reference = values[5]
        rec_curr.evidence_code = values[6]  # TODO:change to object
        rec_curr.taxon = values[12][6:]  # remove "taxon" TODO:change to object
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
