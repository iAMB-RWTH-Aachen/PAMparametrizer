from PAModelpy.utils.pam_generation import merge_enzyme_complexes,get_protein_gene_mapping
from cobra.io import read_sbml_model
import pandas as pd

file = 'Results/1_preprocessing/proteinAllocationModel_iML1515_EnzymaticData_250203.xlsx'

db = pd.read_excel(file, sheet_name='ActiveEnzymes')

model = read_sbml_model('Models/iML1515.xml')

protein2gene, gene2protein = get_protein_gene_mapping(db, model)


# Ensure the enzyme complexes are merged on one row
eco_enzymes_mapped = merge_enzyme_complexes(db, gene2protein)[['rxn_id', 'enzyme_id']]

print(eco_enzymes_mapped)