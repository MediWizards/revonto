import numpy as np

class MultiCorrectionFactory():
    """Factory for choosing a multiple testing correction function."""
    
class CorrectionBase():
    """Base class for local multiple test correction calculations."""
    def __init__(self, pvals, a=.05):
        self.pvals = self.corrected_pvals = np.array(pvals)
        self.n = len(self.pvals)    # number of multiple tests
        self.a = a                  # type-1 error cutoff for each test

        self.set_correction()
        # Reset all pvals > 1 to 1
        self.corrected_pvals[self.corrected_pvals > 1] = 1

    def set_correction(self):
        # the purpose of multiple correction is to lower the alpha
        # instead of the canonical value (like .05)
        pass
    
class Bonferroni(CorrectionBase):
    """
    >>> Bonferroni([0.01, 0.01, 0.03, 0.05, 0.005], a=0.05).corrected_pvals
    array([ 0.05 ,  0.05 ,  0.15 ,  0.25 ,  0.025])
    """
    def set_correction(self):
        """Do Bonferroni multiple test correction on original p-values."""
        self.corrected_pvals *= self.n
        
class BH_FDR(CorrectionBase):
    """Benjamini–Hochberg False Discovery Rate"""
