class ReverseLookupRecord(object):
    """Represents one result (from a single product) in the ReverseLookupStudy"""
    def __init__(self, name, **kwargs):
        # Methods seen in current enrichment result
        # pylint: disable=invalid-name
        self.name = name
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
        self.enrichment = self._init_enrichment()

    def set_corrected_pval(self, nt_method, pvalue):
        """Add object attribute based on method name."""
        self.method_flds.append(nt_method)
        fieldname = "".join(["p_", nt_method.fieldname])
        setattr(self, fieldname, pvalue)