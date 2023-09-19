from typing import Any

import numpy as np
from numpy.typing import NDArray
from statsmodels.stats.multitest import (  # TODO: in future only import when needed and only once
    multipletests,
)


class CorrectionBase:
    """Base class for local multiple test correction calculations."""

    def __init__(self):
        pass

    def set_correction(self, pvals, a=0.05) -> NDArray[Any]:
        # the purpose of multiple correction is to lower the alpha
        # instead of the canonical value (like .05)
        return NotImplemented


class Bonferroni(CorrectionBase):
    """
    >>> Bonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals
    array([ 0.05 ,  0.05 ,  0.15 ,  0.25 ,  0.025])
    """

    def set_correction(self, pvals: list[float], a=0.05) -> NDArray[Any]:
        """Do Bonferroni multiple test correction on original p-values."""
        p_uncorrected = np.array(pvals)
        n = len(p_uncorrected)

        corrected_pvals = n * p_uncorrected
        corrected_pvals[corrected_pvals > 1] = 1
        return (
            corrected_pvals  # return NDarray, since statsmodels also return this type
        )


class BH_FDR(CorrectionBase):
    """Benjamini–Hochberg False Discovery Rate"""


class Statsmodels(CorrectionBase):
    def __init__(self, method_name):
        self.statsmodels_multitest = multipletests
        self.method = method_name

    def set_correction(self, pvals, a=0.05) -> NDArray[Any]:
        results = self.statsmodels_multitest(pvals, a, self.method)
        corrected_pvals = results[
            1
        ]  # reject_lst, pvals_corrected, alphacSidak, alphacBonf

        return corrected_pvals


class MultiCorrectionFactory:
    """Factory for choosing a multiple testing correction function."""

    options = {
        "bonferroni": Bonferroni(),
        "sidak": None,
        "holm": None,
        "fdr": BH_FDR(),
        "sm_bonferroni": Statsmodels(
            "bonferroni"
        ),  #  0) Bonferroni one-step correction
        "sm_sidak": Statsmodels("sidak"),  #  1) Sidak one-step correction
        "sm_holm-sidak": Statsmodels(
            "holm-sidak"
        ),  #  2) Holm-Sidak step-down method using Sidak adjustments
        "sm_holm": Statsmodels(
            "holm"
        ),  #  3) Holm step-down method using Bonferroni adjustments
        "sm_simes-hochberg": Statsmodels(
            "simes-hochberg"
        ),  #  4) Simes-Hochberg step-up method  (independent)
        "sm_hommel": Statsmodels(
            "hommel"
        ),  #  5) Hommel closed method based on Simes tests (non-negative)
        "sm_fdr_bh": Statsmodels(
            "fdr_bh"
        ),  #  6) FDR Benjamini/Hochberg  (non-negative)
        "sm_fdr_by": Statsmodels("fdr_by"),  #  7) FDR Benjamini/Yekutieli (negative)
        "sm_fdr_tsbh": Statsmodels(
            "frd_tsbh"
        ),  #  8) FDR 2-stage Benjamini-Hochberg (non-negative)
        "sm_fdr_tsbky": Statsmodels(
            "fdr_tsbky"
        ),  #  9) FDR 2-stage Benjamini-Krieger-Yekutieli (non-negative)
        "sm_fdr_gbs": Statsmodels(
            "fdr_gbs"
        ),  # 10) FDR adaptive Gavrilov-Benjamini-Sarkar
    }

    def __init__(self, correction="bonferroni") -> None:
        if correction not in self.options:
            raise ValueError(f"correction must be one of {self.options.keys()}")
        self.corr_fnc_name = correction
        self.corr_obj: CorrectionBase = self._init_corr_obj()

    def _init_corr_obj(self):
        """Returns pvalue object based on the input"""
        return self.options[self.corr_fnc_name]
