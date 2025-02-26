import pandas as pd
import re
import os
from typing import Union, Tuple, Literal
from cobra.io.sbml import read_sbml_model

# load PAMpy modules
from PAModelpy import PAModel
from PAModelpy.utils.pam_generation import set_up_pam, parse_reaction2protein, _order_enzyme_complex_id
from PAModelpy.EnzymeSectors import ActiveEnzymeSector, UnusedEnzymeSector, TransEnzymeSector, CustomSector
from PAModelpy.configuration import Config


DEFAULT_MOLMASS = 39959.4825 #kDa
DEFAULT_KCAT = 11 #s-1

def create_pamodel_from_diagnostics_file(file_path:str,
                                         model: PAModel,
                                         sheet_name: str = 'Best_Individuals')-> PAModel:
    """
    Modifies a Protein Allocation Model using information about turnover numbers from a diagnostics file (result from PAMparametrizer)

    Args:
        file_path (str): path to the diagnostics xlsx file. The file should at least have the following columns:
            - run_id: the iteration of the PAMparametrizer
            - rxn_id: the id of the reaction to modify. This can be the catalytic reaction id (CE_<rxn_id>_<enzyme_id>)
            - enzyme_id: the id of the enzyme to modify
            - direction: 'f' or 'b', determines the directionality of the reaction
            - kcat[s-1]: the new kcat value in 1/s
        model (PAModel): the PAM to adjust
        sheet_name (str): name of the sheet with the information about the modifications

    Returns:
        PAModel: with adjusted parameters

    """
    best_individual_df = pd.read_excel(file_path, sheet_name=sheet_name)
    for _, group in best_individual_df.groupby('run_id'):
        for _, row in group.iterrows():
            rxn_id = _extract_reaction_id_from_catalytic_reaction_id(row['rxn_id'])
            enzyme_id = _order_enzyme_complex_id(row['enzyme_id'])
            kcat_dict = {rxn_id: {row['direction']: row['kcat[s-1]']}}
            model.change_kcat_value(enzyme_id=enzyme_id, kcats=kcat_dict)
    return model

def get_rxn2kcat_protein2gene_dict(param_file_path:str,
                                   model_file_path: str
                                   ) -> Tuple[
    dict[str, dict[str,dict[
        Literal['f', 'b', 'molmass', 'protein_reaction_association'
        ], float]
    ]],
    dict[str,str]]:
    # create enzyme objects for each gene-associated reaction
    pam = set_up_pam(param_file_path, model_file_path, sensitivity=False)
    enzyme_db = pd.read_excel(param_file_path, sheet_name='ActiveEnzymes')
    rxn2protein, protein2gene = parse_reaction2protein(enzyme_db, pam)

    ae_sector = pam.sectors.ActiveEnzymeSector
    new_rxn2prot = ae_sector.rxn2protein.copy()
    for rxn, enz_dict in ae_sector.rxn2protein.items():
        if rxn[:2]=='CE': continue

        for enzyme_id, enzyme_dict in enz_dict.items():
            protein_reaction = enzyme_dict['protein_reaction_association']
            if not ae_sector._enzyme_is_enzyme_complex(protein_reaction, enzyme_id): continue

            for pr in protein_reaction:
                if not len(pr) > 1: continue

                enzyme_complex_id = '_'.join(pr)
                new_rxn2prot[rxn] = {**new_rxn2prot[rxn],
                                     **{enzyme_complex_id: enzyme_dict}}
    return new_rxn2prot, protein2gene

def _extract_reaction_id_from_catalytic_reaction_id(input_str: str,
                                                    default_enzyme_id_pattern: str = r'E[0-9][0-9]*|Enzyme_*') -> str:
    # Define the regex pattern for protein IDs, obtained from UniProtKB, 2024-08-07
    # https://www.uniprot.org/help/accession_numbers
    protein_id_pattern = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'

    # Remove the 'CE_' prefix if it exists
    if input_str.startswith('CE_'):
        input_str = input_str[3:]

    # Define the regex pattern to match protein IDs
    protein_id_regex = re.compile(r'_' + protein_id_pattern + r'|' + r'_' + default_enzyme_id_pattern)
     # split off all protein ids from the reaction
    reaction_id = protein_id_regex.split(input_str)[0]
    # Remove any trailing or leading underscores that might remain
    return reaction_id.strip('_')

def _get_rxn2kcat_as_series(rxn2kcat: dict[str, dict],
                                   name: str):

    kcats = {}
    for rxn,enz_dict in rxn2kcat.items():
        for enz, kcat_dict in enz_dict.items():
            for direction, kcat in kcat_dict.items():
                if len(direction) == 1:
                    kcats[f"{rxn}_{enz}_{direction}"] = kcat
    return pd.Series(kcats, name = name)


def search_index_in_parameter_file(df:pd.DataFrame,
                                   protein:str,
                                   reaction:str,
                                   direction:str) -> pd.Index:
    """
    Get the index of an enzyme-reaction-direction relationship in a dataframe containing information of on the parameters
    """
    all_protein_rows = df[df.enzyme_id == protein]
    all_rxn_and_protein_rows = all_protein_rows[all_protein_rows.rxn_id == reaction]
    the_row = all_rxn_and_protein_rows[all_rxn_and_protein_rows.direction == direction]
    return the_row.index

def create_new_aes_parameter_file(old_param_file:str,
                                  diagnostics_file_path:str,
                                  new_aes_file: str,
                                  diagnostics_sheet_name: str = 'Best_Individuals') -> None:
    """
    Updates the 'ActiveEnzymes' sheet in a parameter file with new kcat values
    based on the results from the PAMparametrizer and saves the modified file.

    This function reads an existing parameter file and updates the kcat values
    in the 'ActiveEnzymes' sheet based on the best parameterization results
    from a diagnostics file. The updated file is then saved to a new location, as defined in the new_aes_file.

    Args:
        old_param_file (str): Path to the existing parameter file (Excel format).
        diagnostics_file_path (str): Path to the diagnostics file containing
            the best kcat parameterizations (Excel format).
        new_aes_file (str): Path where the updated parameter file will be saved.
        diagnostics_sheet_name (str): name of the sheet where the information on the new parameters is stored.
            Implemented to enable easy testing.

    Returns:
        None: The function writes the updated parameter file to the specified location.

    Raises:
        FileNotFoundError: If `old_param_file` or `diagnostics_file_path` do not exist.
    """

    # Check if input files exist
    if not os.path.isfile(old_param_file):
        raise FileNotFoundError(f"Parameter file not found: {old_param_file}")
    if not os.path.isfile(diagnostics_file_path):
        raise FileNotFoundError(f"Diagnostics file not found: {diagnostics_file_path}")

    parameter_files = pd.read_excel(old_param_file, sheet_name=None)
    aes_parameter_file = parameter_files['ActiveEnzymes']

    parametrization_results = pd.read_excel(
        diagnostics_file_path, sheet_name=diagnostics_sheet_name).drop_duplicates(
        ['rxn_id', 'direction', 'enzyme_id'], keep='last') #the last one refers to the last modification of the kcat value

    # extract reaction id from catalytic event id
    parametrization_results['rxn_id'] = [_extract_reaction_id_from_catalytic_reaction_id(rid) for rid in
                                         parametrization_results['rxn_id']]

    for index, row in parametrization_results.iterrows():
        aes_index = search_index_in_parameter_file(aes_parameter_file, row.enzyme_id, row.rxn_id, row.direction)
        aes_parameter_file.loc[aes_index, 'kcat_values'] = row['kcat[s-1]']

    write_mode = 'w'
    kwargs = {}

    #make sure the file is either replaced or created
    if os.path.isfile(new_aes_file):
        write_mode = 'a'
        kwargs = {'if_sheet_exists': 'replace'}

    with pd.ExcelWriter(new_aes_file,
            mode=write_mode, engine='openpyxl', **kwargs) as writer:
        aes_parameter_file.to_excel(writer, sheet_name='ActiveEnzymes', index=False)
        for sheet, df in parameter_files.items():
            if sheet != 'ActiveEnzymes':
                df.to_excel(writer, sheet_name = sheet, index=False)



def setup_yeast_pam(pam_info_file:str= os.path.join('Data', 'proteinAllocationModel_yeast9_EnzymaticData_TurnUp.xlsx'),
                     model:str = 'Models/yeast9.xml', config:Config = None,
                     total_protein: Union[bool, float] = 0.28697423725932236, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True) -> PAModel:
    if config is None:
        config = set_up_yeast_config()
    yeast_pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    if total_protein:
        #the data we used to calibrate the model includes also housekeeping proteins in the UP section. Thus we
        #can assume include this section in the total protein constraint (no correction needed)
        protein_growth_relation = CustomSector(
            id_list= ['r_2111'], #biomass formation
            name='total_protein_growth_relation',
            cps_0=[0],
            cps_s=[-0.45431074007034344],
            mol_mass=1e6
        )
        yeast_pam.add_sector(protein_growth_relation)
    return yeast_pam


def set_up_yeast_config():
    config = Config()
    config.TOTAL_PROTEIN_CONSTRAINT_ID = "TotalProteinConstraint"
    config.P_TOT_DEFAULT = 0.388  # g_protein/g_cdw
    config.CO2_EXHANGE_RXNID = "r_1672"
    config.GLUCOSE_EXCHANGE_RXNID = "r_1714"
    config.BIOMASS_REACTION = "r_2111"
    config.OXYGEN_UPTAKE_RXNID = "r_1992"
    config.ACETATE_EXCRETION_RXNID = "r_1634"
    config.PHYS_RXN_IDS = [
    config.BIOMASS_REACTION,
    config.GLUCOSE_EXCHANGE_RXNID,
    config.ACETATE_EXCRETION_RXNID,
    config.CO2_EXHANGE_RXNID,
    config.OXYGEN_UPTAKE_RXNID]
    config.ENZYME_ID_REGEX = r'(Y[A-P][LR][0-9]{3}[CW])'
    return config

def setup_pputida_pam(pam_info_file:str= os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iJN1463_EnzymaticData_250117.xlsx'),
                     model:str = 'Models/iJN1463.xml',
                     total_protein: Union[bool, float] = 0.3, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'BIOMASS_KT2440_WT3'
    #make sure default enzyme ids for locus tags without enzyme name are parsed correctly
    config.ENZYME_ID_REGEX += r'|Enzyme_*|Enzyme_PP_[0-9]+'
    pputida_pam = set_up_pam(pam_info_file, model, config,
                     total_protein, active_enzymes, translational_enzymes,
                     unused_enzymes, sensitivity = sensitivity)
    return pputida_pam


def setup_cglutanicum_pam(pam_info_file:str= os.path.join(
                                             'Results', '1_preprocessing',
                                             'proteinAllocationModel_iCGB21FR_EnzymaticData_250214.xlsx'),
                     model:str = 'Models/iCGB21FR_annotated_copyable.xml',
                     total_protein: Union[bool, float] = 0.3, active_enzymes: bool = True,
                    translational_enzymes: bool = True, unused_enzymes: bool = True, sensitivity = True):
    config = Config()
    config.reset()
    config.BIOMASS_REACTION = 'Growth'

    model = read_sbml_model(model)
    #change medium to CGXII by removing 3,4-Dihydroxybenzoate
    model.reactions.get_by_id('EX_34dhbz_e').lower_bound = 0

    cglutanicum_pam = set_up_pam(pam_info_file,
                                 model = model,
                                 config = config,
                                 total_protein = total_protein,
                                 active_enzymes = active_enzymes,
                                 translational_enzymes = translational_enzymes,
                                 unused_enzymes = unused_enzymes,
                                 sensitivity = sensitivity)
    return cglutanicum_pam
