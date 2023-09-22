"""Options for calculating uncorrected p-values."""


class PvalueCalcBase:
    """Base class for calculating p-value"""

    def __init__(self, pval_fnc) -> None:
        self.pval_fnc = pval_fnc

    def calc_pvalue(
        self, study_count: int, study_n: int, pop_count: int, pop_n: int
    ) -> float:
        """pvalues are calculated in derived classes"""
        raise NotImplementedError("Must be called from derived classes")


class FisherScipyStats(PvalueCalcBase):
    """Fisher exact test"""

    def __init__(self) -> None:
        from scipy import stats

        super().__init__(stats.fisher_exact)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected pvalues."""
        #               YES       NO
        # study_genes    a scnt   b    | a+b = study_n
        # not s_genes    c        d    |  c+d
        #             --------   -----
        #   pop_genes  a+c pcnt   b+d  a+b+c+d = pop_n

        avar = study_count
        bvar = study_n - study_count
        cvar = pop_count - study_count
        dvar = pop_n - pop_count - bvar
        if cvar < 0:
            raise ValueError(
                f"STUDY={avar}/{bvar} POP={cvar}/{dvar} scnt({study_count}) stot({study_n}) pcnt({pop_count}) ptot({pop_n})"
            )
        _, p_uncorrected = self.pval_fnc([[avar, bvar], [cvar, dvar]])
        return p_uncorrected


class BinomialScipyStats(PvalueCalcBase):
    """Binomial test"""

    def __init__(self) -> None:
        from scipy import stats

        super().__init__(stats.binomtest)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected pvalues."""
        k = study_count
        n = study_n
        p = pop_count / pop_n
        if pop_n == 0:
            raise ZeroDivisionError("pop_n is equal to zero")
        p_uncorrected = self.pval_fnc(k, n, p).pvalue
        return p_uncorrected


class PValueFactory:
    """Factory for choosing an algorithm"""

    options = {
        "fisher_scipy_stats": FisherScipyStats,
        "binomial_scipy_stats": BinomialScipyStats,
    }

    def __init__(self, pvalcalc="fisher_scipy_stats") -> None:
        if pvalcalc not in self.options:
            raise ValueError(f"pvalcalc must be one of {self.options.keys()}")
        self.pval_fnc_name = pvalcalc
        self.pval_obj: PvalueCalcBase = self._init_pval_obj()

    def _init_pval_obj(self):
        """Returns pvalue object based on the input"""
        return self.options[self.pval_fnc_name]()
