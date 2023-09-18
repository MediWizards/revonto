from __future__ import annotations
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from .associations import Annotations
    from .ontology import GODag

from .pvalcalc import PValueFactory
from .multiple_testing import MultiCorrectionFactory


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

    def set_corrected_pval(self, method, pvalue):
        """Add object attribute based on method name."""
        method = method.replace("-", "_")
        self.method_flds.append(method)
        fieldname = "".join(["p_", method])
        setattr(self, fieldname, pvalue)


class GOReverseLookupStudy:
    """Runs pvalue test, as well as multiple corrections"""

    def __init__(
        self,
        anno: Annotations,  # this is annotation object. This is the population. (preprocess it to add orthologs or to propagate associations to parents). NOTE: species you add to the association object affect the result; only include the target species and the ones ortologs were searched for.
        obo_dag: GODag,  # check if it is needed?
        alpha=0.05,
        pvalcalc="fisher_scipy_stats",
        methods=None,
    ):
        self.anno = anno
        self.obo_dag = obo_dag
        self.alpha = alpha
        if methods is None:
            self.methods = ["bonferroni"]  # add statsmodel multipletest
        self.pval_obj = PValueFactory(pvalcalc).pval_obj

    def run_study(self, study, **kws) -> List[ReverseLookupRecord]:
        """_summary_

        Args:
            study (_type_): list of all goterms (term_id) for a process-

        Returns:
            List[ReverseLookupRecord]: _description_
        """
        """Run Gene Ontology Reverse Lookup Study"""

        if len(study) == 0:
            return []

        # process kwargs
        methods = kws.get("methods", self.methods)
        alpha = kws.get("alpha", self.alpha)

        # calculate the uncorrected pvalues using the pvalcalc of choice
        results = self.get_pval_uncorr(
            study
        )  # results is a list of ReverseLookupRecord objects
        if not results:
            return []

        # do multipletest corrections on uncorrected pvalues, add to ReverseLookupRecord objects
        self._run_multitest_corr(results, methods, alpha)

        # 'keep_if' can be used to keep only significant GO terms. Example:
        #     >>> keep_if = lambda nt: nt.p_fdr_bh < 0.05 # if results are significant
        #     >>> goea_results = goeaobj.run_study(geneids_study, keep_if=keep_if)
        if "keep_if" in kws:
            keep_if = kws["keep_if"]
            results = [r for r in results if keep_if(r)]

        # Default sort order:
        #results.sort(key=lambda r: [r.p_uncorrected])

        return results  # list of ReverseLookupRecord objects

    def get_pval_uncorr(self, study) -> List[ReverseLookupRecord]:
        """Calculate the uncorrected pvalues for study items."""
        results = []

        dict_by_object_id = self.anno.dict_from_attr("object_id")
        dict_by_term_id = self.anno.dict_from_attr("term_id")

        study2annoobjid = (
            set()
        )  # list of all annotation objects id from goterms in study
        for term_id in study:
            for annoobj in dict_by_term_id.get(term_id, set()):
                study2annoobjid.add(annoobj.object_id)

        for object_id in study2annoobjid:
            # for each object id (product id) calculate pvalue
            study_items = set(
                anno_obj.term_id
                for anno_obj in dict_by_object_id[object_id]
                if anno_obj.term_id in study
            )
            study_count = len(
                study_items
            )  # for each object id (product id) check how many goterms in study are associated to it
            study_n = len(study)  # N of study set

            population_items = set(
                anno_obj.term_id for anno_obj in dict_by_object_id[object_id]
            )

            pop_count = len(
                population_items
            )  # total number of goterms an objectid (product id) is associated in the whole population set
            pop_n = len(dict_by_term_id)  # total number of goterms in population set

            one_record = ReverseLookupRecord(
                object_id,
                p_uncorrected=self.pval_obj.calc_pvalue(
                    study_count, study_n, pop_count, pop_n
                ),
                study_items=study_items,
                population_items=population_items,
                ratio_in_study=(study_count, study_n),
                ratio_in_pop=(pop_count, pop_n),
            )

            results.append(one_record)

        return results

    def _run_multitest_corr(self, results, methods, a):
        pvals = [r.p_uncorrected for r in results]
        for method in methods:
            corrected_pvals = MultiCorrectionFactory(method).corr_obj.set_correction(
                pvals, a
            )
            self._update_pvalcorr(results, method, corrected_pvals)

    @staticmethod
    def _update_pvalcorr(results, method, corrected_pvals):
        """Add data members to store multiple test corrections."""
        if corrected_pvals is None:
            return
        for rec, val in zip(results, corrected_pvals):
            rec.set_corrected_pval(method, val)
