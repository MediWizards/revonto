"""
Read and store Gene Ontology's GAF (GO Annotation File).
"""
from typing import Set
import os

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
    def __init__(self) -> None:
        self.term_id = None
        self.relationship = None
        self.NOTrelation = False
        self.reference = []
        self.evidence_code = None
        self.taxon = None
        self.date = None

class Annotations(dict[str, Set[Annotation]]):
    """
    Store Annotations as a Dict with key "ID" and value a set of all Annotation objects connected to that product.
    This is how the classes preforming the study expect the data to be formated.

    ID can be anythig (genename, DB:ID) depending on your pre and post analysis.
    """
    def __init__(self, assoc_file="go-basic.obo"):
        super().__init__()
        self.version, self.date = self.load_assoc_file(assoc_file)

    def load_assoc_file(self, assoc_file):
        """read association file"""

        extension = os.path.splitext(assoc_file)[1]
        if extension == ".gaf":
            reader = GafParser(assoc_file)
        if extension == ".gpad":
            raise NotImplementedError("GPAD files are not yet supported")
        
        for obj_id, rec in reader:
            self.setdefault(obj_id, set()).add(rec)

        return reader.version, reader.date

        

class AnnoParserBase():
    """
    There is more than one type of annotation file. 
    Therefore we will use a base class to standardize the data and the methods.
    
    Currently we only support GAF, beacuse we need 
    """
    
    def __init__(self,assoc_file) -> None:
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

    def __iter__(self) -> (str, Annotation):
        with open(self.assoc_file) as fstream:
            hdr = True

            for line in fstream:
                line = line.rstrip()
                if hdr:
                    if not self._init_hdr(line):
                        hdr = False
                if not hdr and line:
                    values = line.split("\t")
                    object_id = values[0] + ":" + values[1]
                    rec_curr = Annotation()
                    self._add_to_ref(rec_curr, values)
                    yield object_id, rec_curr

    def _init_hdr(self, line:str):
        """save gaf version and date"""
        if line[:14] == "!gaf-version: ":
            self.version = line[14:]
            return True
        if line[:17] == "!date-generated: ":
            self.date = line [17:]
            return True
        if line[0] != "!":
            return False
        return True
    
    def _add_to_ref(self, rec_curr:Annotation, values):
        """populate Annotation object with values from line"""
        rec_curr.term_id = values[4]
        rec_curr.relationship = values[3] #change to object
        if "NOT" in values[3]: rec_curr.NOTrelation = True
        rec_curr.reference = values[5]
        rec_curr.evidence_code = values[6] #change to object
        rec_curr.taxon = values[12] #change to object
        rec_curr.date = values[13]

class EvidenceCodes():
    """
    class which holds information about evidence codes.
    upon creation the fields are populated accordint to the evicence code in __init__
    currently not used
    """
    codes = {}
    def __init__(self, code) -> None:
        if code not in self.codes:
            pass
