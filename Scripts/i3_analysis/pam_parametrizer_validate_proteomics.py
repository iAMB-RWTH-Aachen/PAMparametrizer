import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from PAModelpy.utils.pam_generation import set_up_pam
from Modules.utils.pamparametrizer_analysis import get_results_from_simulations

PAM_KCAT_FILES = [os.path.join('Results', '3_analysis', 'parameter_files',
                                   f'proteinAllocationModel_EnzymaticData_iML1515_{file_nmbr}.xlsx') for file_nmbr in
                      range(1, 7)]
ORI_PAM_KCAT_FILE = os.path.join('Results', '1_preprocessing',
                                 'proteinAllocationModel_iML1515_EnzymaticData_241209.xlsx')
ECOLI_PROTEOME_DATA_PATH = os.path.join('Data', 'Ecoli_phenotypes', 'proteome_data_extract_schmidt2016.xlsx')
ECOLI_MODEL_FILE_PATH = os.path.join('Models', 'iML1515.xml')
UNIPROT_INFO_FILE = os.path.join('Data', 'Databases', 'uniprotkb_ecolik12_240726.xlsx')

locustag_regex =r'\b([b|s]\d{4})\b'

def get_simulated_protein_conc(pam_param_files: list[str] = [ORI_PAM_KCAT_FILE] + PAM_KCAT_FILES,
                               model_file: str = ECOLI_MODEL_FILE_PATH,
                               substrate_rates:list[str]=[-10, -6, -4.5])-> dict:
    proteomics_results = dict()

    for file in pam_param_files:
        model_name = file.split('.')[0].split('_')[-1]
        model = set_up_pam(file, model_file, sensitivity=False)
        enzymes = [enz.id for enz in model.enzyme_variables if enz._model is not None]
        proteomics_results[model_name] = get_results_from_simulations(pamodel=model,
                                                                      substrate_rates=substrate_rates,
                                                                      proteins_to_save=enzymes)['proteins']

    return proteomics_results

def get_proteomics_literature(proteomics_file_path: str = ECOLI_PROTEOME_DATA_PATH) -> pd.DataFrame:
    proteome_df = pd.read_excel(proteomics_file_path,
                                sheet_name='ProteinMasses',
                                engine='openpyxl',
                                index_col=0)

    # only get the chemostat cultivations on glucose
    proteome_glc = proteome_df[['Glucose', 'Chemostat µ=0.5', 'Chemostat µ=0.35']]  # unit: fg/cell
    proteome_glc.columns = ['Batch', 'mu_5', 'mu_35']
    return proteome_glc

def map_genes_to_uniprot(proteome_glc: pd.DataFrame,
                         simulation_results: pd.DataFrame,
                         uniprot_file: str =UNIPROT_INFO_FILE,
                         locustag_regex: str=locustag_regex) -> pd.DataFrame:
    # use information from uniprot to map the locus tag ids to protein identifiers
    uniprot_info_df = pd.read_excel(uniprot_file)
    # only keep those who are in the model
    uniprot_in_model = uniprot_info_df.loc[
        [prot in simulation_results.enzyme_id.to_list() for prot in uniprot_info_df.Entry]]

    # get the gene id from the gene names
    uniprot_in_model['b_number'] = uniprot_in_model['Gene Names'].str.extract(locustag_regex)
    uniprot_df = uniprot_in_model[['b_number', 'Entry']]
    proteome_glc_mapped = pd.merge(proteome_glc, uniprot_df, how='right',
                                   left_on='Bnumber', right_on='b_number')

    #normalize the protein concentrations
    proteome_glc_normalized = proteome_glc_mapped.copy()
    # proteomics data
    for col in proteome_glc.columns:
        proteome_glc_normalized[col] = proteome_glc_mapped[col].div(proteome_glc_mapped[col].sum())
    # convert to long_format
    proteome_glc_long = pd.melt(proteome_glc_normalized, value_vars=['Batch', 'mu_5', 'mu_35'],
                                id_vars=['b_number', 'Entry'],
                                var_name='experiment', value_name='fraction')
    return proteome_glc_long

def normalize_simulated_protein_concentrations(proteomics_results:pd.DataFrame) -> pd.DataFrame:
    prot_predicted_long = pd.DataFrame(columns=['method', 'enzyme_id'])

    for model, protein_df in proteomics_results.items():
        protein_df = parse_protein_df(protein_df)
        prot_predicted_long = pd.merge(prot_predicted_long, protein_df[['method', 'enzyme_id', 'fraction']],
                                       on=['method', 'enzyme_id'], how='outer', copy=False,
                                       suffixes=('', f'_{model}'))

    return prot_predicted_long

def parse_protein_df(df):
    substrate2method = {-10: 'Batch', -6: 'mu_5', -4.5: 'mu_35'}
    df = df.dropna(how='any')
    df.enzyme_id = df.enzyme_id.str.split('_')
    # df['fraction_per_peptide'] = df.apply(lambda x: x.fraction / len(x.enzyme_id), axis=1)
    # df = df.explode('enzyme_id')
    #
    # df['method'] = df.apply(lambda x: substrate2method[x['substrate_uptake']], axis=1)
    # # Calculate sum of fraction grouped by substrate_uptake
    # sum_by_substrate = df.groupby('substrate_uptake')['fraction_per_peptide'].transform('sum')
    # # Divide fraction values by the sum
    # df['fraction'] = df['fraction'] / sum_by_substrate
    return df

if __name__ == '__main__':
    simulation_proteins_results = get_simulated_protein_conc()
    proteomics_glc = get_proteomics_literature()
    proteomics_mapped = map_genes_to_uniprot(proteomics_glc, simulation_proteins_results['1'])
    simulated_proteins_long = normalize_simulated_protein_concentrations(simulation_proteins_results)
    print(simulated_proteins_long.to_markdown())
