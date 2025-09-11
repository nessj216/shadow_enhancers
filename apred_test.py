# import the module
from adpred import ADpred as adp

# with own sequence:
myProtein = adp.protein("myProteinName", "ATREFERTATREFERTAADDWLNDCWATREFERTA")

# with uniprot identifier (example gcn4)
myProtein = adp.protein('GCN4_YEAST')

# if secondary structure is not known:
myProtein.predict()   # This will predict secondary structure

# If you wish to first get the secondary structure
myProtein.predict_second_struct()

# do saturated mutagenesis analysis to 30mer between residue positions 108 and 138.
myProtein.saturated_mutagenesis(start=108, end=138)

# By default the WT structure is used for all mutants, however,
# If you wish to recalculate the secondary structure for each mutant
myProtein.saturated_mutagenesis(start=108, end=138, 'second_struct_on_each_mutant')