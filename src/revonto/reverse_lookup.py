from __future__ import annotations
from typing import TYPE_CHECKING, List
if TYPE_CHECKING:
    from .associations import Annotations
    from .ontology import GODag

from revonto.pvalcalc import PValueFactory
from revonto.associations import anno2objkey
class ReverseLookupRecord(object):
    """Represents one result (from a single product) in the ReverseLookupStudy"""
    def __init__(self, objid, **kwargs):
        self.object_id = objid
        self.name = "n.a"
        self.method_flds = []
        self.kws = kwargs
        # Ex: ratio_in_pop ratio_in_study study_items p_uncorrected pop_items
        for key, val in kwargs.items():
            setattr(self, key, val)
        _stucnt = kwargs.get("ratio_in_study", (0, 0))
        self.study_count = _stucnt[0]
        self.study_n = _stucnt[1]
        _popcnt = kwargs.get("ratio_in_pop", (0, 0))
        self.pop_count = _popcnt[0]
        self.pop_n = _popcnt[1]

class GOReverseLookupStudy():
    """Runs pvalue test, as well as multiple corrections"""
    def __init__(
        self,
        anno: Annotations,  #this is annotation object. This is the population. (preprocess it to add orthologs or to propagate associations to parents). NOTE: species you add to the association object affect the result; only include the target species and the ones ortologs were searched for.
        obo_dag: GODag, #check if it is needed?
        alpha=0.05,
        pvalcalc="fisher_scipy_stats",
        methods=None,
        ):
        self.anno = anno
        self.obo_dag = obo_dag
        self.alpha = alpha
        if methods is None:
            self.methods = ["bonferroni"] #add statsmodel multipletest
        self.pval_obj = PValueFactory(pvalcalc).pval_obj

    def run_study(self, study, **kws) -> List[ReverseLookupRecord]:
        """Run Gene Ontology Reverse Lookup Study"""

        if len(study) == 0:
            return []
        
        #process kwargs
        methods = kws.get("methods", self.methods)
        alpha = kws.get("alpha", self.alpha)

        #calculate the uncorrected pvalues using the pvalcalc of choice
        results = self.get_pval_uncorr(study) #results is a list of ReverseLookupRecord objects
        if not results:
            return []
        
        #do multipletest corrections on uncorrected pvalues, add to ReverseLookupRecord objects
        #self._run_multitest_corr(results, methods, alpha, study)

        #results.sort(key=lambda r: [r.enrichment, r.NS, r.p_uncorrected])
        return results #list of ReverseLookupRecord objects
    
    def get_pval_uncorr(self, study) -> List[ReverseLookupRecord]:
        """Calculate the uncorrected pvalues for study items."""
        results = []

        anno2objkeydict = anno2objkey(self.anno)  #dictionary with all annotations with object (product) id as keys instead of goterms
        study2annoobjid = set() #list of all annotation objects id from goterms in study
        for goid in study:
            for annoobj in self.anno.get(goid, set()):
                study2annoobjid.add(annoobj.object_id)
        

        for objid in study2annoobjid:
            #for each object id (product id) calculate pvalue
            study_items = set(termanno for termanno in anno2objkeydict[objid] if termanno.term_id in study)
            study_count = len(study_items) #for each object id (product id) check how many goterms in study are associated to it
            study_n = len(study) # N of study set

            pop_count = len(anno2objkeydict[objid]) # total number of goterms an objectid (product id) is associated in the whole population set
            pop_n = len(self.anno) # total number of goterms in population set

            one_record = ReverseLookupRecord(
                objid,
                p_uncorrected=self.pval_obj.calc_pvalue(study_count, study_n, pop_count, pop_n),
                study_items=study_items,
                population_items=anno2objkeydict[objid],
                ratio_in_study=(study_count, study_n),
                ratio_in_pop=(pop_count, pop_n)
            )

            results.append(one_record)

        return results






