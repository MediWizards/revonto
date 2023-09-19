from revonto.associations import Annotations, propagate_associations, match_annotations_to_godag
from revonto.ontology import GODag
from revonto.reverse_lookup import GOReverseLookupStudy, results_intersection

# import GO term lists. They are in separate file to keep this file cleaner.
from studysets_for_cancer_inflamation import studyset_cancer, studyset_infla
# There is no need to store GO term lists like that. You can store it in any format you want.
# You only need to write a custom parser to read the file and produce a list of GO terms

godag = GODag("go.obo")

anno = Annotations( "goa_human.gaf")

# godag defines the population of GOTerms, if you don't intend to use the whole go.obo,
# make sure to match the Annotations object to it. You only need to match it once after the import.
# Latter modifications of Annotations object (such as propagation) dont "unmatch" it, 
# therefore you can run match after propagation and the result will be the same.
# If you combine (union, add...) multiple Annotations objects, make sure that each object was matched
# or you need to match the resulting combined Annotations object.

#from revonto.associations import match_annotations_to_godag
#matched_anno = match_annotations_to_godag(anno, godag)

# If you would like to include indirect annotaions (from children) propagate them!
propagated_anno = propagate_associations(anno, godag)

# you now have two Annotations objects: anno and propagated_anno
# be sure to use the one you want in code below

# ortholog function will come here and will modify Annotations object

# setup the study. Please note that anno Annotations object is used. Change it so propagated_anno if you want.
study = GOReverseLookupStudy(anno, godag, alpha=0.05, pvalcalc="fisher_scipy_stats", methods=["bonferroni"])

# run the study with the studysets
results_infla = study.run_study(studyset_infla)
results_cancer = study.run_study(studyset_cancer)

# the results_xxx are lists of ReverseLookupRecords. Here are some important attributes:
# object_id - unique product identifier, p_{uncorrected|method name} - these atributes hold p-values
# p_uncorrected is p value from pvalcalc, p_bonferroni is corrected by Bonferroni, p_fdr, p_sm_bonferroni...
# study_count is the count of the intersection of all goterms associated with object_id in Annotations object and the studyset
# study_n is the count of all goterms in the 
# tudyset
# pop_count is the the count of all the goterms associated with object_id in Annotataions object
# pop_n is the count of all the goterms in the GODag object (obo) - please note that depending on your analysis
# Note: GODag doesn't necessarly include all the GOTerms in GO. Perhaps you built the study subset only from Molecular Function part of GO. Or with a subset.

# check which were the significant products from each subset
significant_infla = [r for r in results_infla if r.p_bonferroni < 0.05]
significant_cancer = [r for r in results_cancer if r.p_bonferroni < 0.05]

print([r.object_id for r in significant_infla])
print([r.object_id for r in significant_cancer])

# intersect the results list. The results_intersection can recieve any number of lists
# and gives back a dict of all object_id s which are present in all the lists.
# key = object_id, value = list of ReverseLookupRecords - one matched from each list
significant_intersection = results_intersection(significant_infla, significant_cancer)

print(significant_intersection.keys())