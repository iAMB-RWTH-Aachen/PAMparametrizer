import os
import pandas as pd
import numpy as np
from typing import Literal, Tuple, Iterable, Union, List, Dict
import datetime
import ast

from Modules.utils.sector_config_functions import SectorParameterDict
from Modules.PAM_parametrizer import SectorConfig


def save_sector_information_to_excel(
        param_vs_lin_rxn: SectorParameterDict,
        param_vs_growth: SectorParameterDict,
        biomass_rxn:str,
        lin_rxn_id:str,
        sector_id:Literal['Translational', 'UnusedEnzyme'],
        pam_data_file:str = None,
        substrate_range: Iterable[Union[float, int]] = np.arange(-4,0,1)
)-> None:
    """Saves parameter information for a given protein allocation sector to an Excel file.

    This function writes parameter values related to both growth and linear fluxes
    into a structured Excel file with pre-existing PAM sector data (can be dummy values). Providing the relation to
    growth rate helps to translate the relation to other conditions

    Args:
        param_vs_lin_rxn (SectorParameterDict): Parameters for the linear reaction (slope and intercept).
        param_vs_growth (SectorParameterDict): Parameters for the biomass growth (slope and intercept).
        biomass_rxn (str): Reaction ID for the biomass reaction.
        lin_rxn_id (str): Reaction ID for the linear reaction to which the sector is related.
        sector_id (Literal['Translational', 'UnusedEnzyme']): Sector identifier to update, the identifier relates to the
            sheet name as defined in the Excel file which contains the sector information.
        pam_data_file (str, optional): Path to the input Excel file with PAM data. If None, attempts to load default.
        substrate_range (Iterable[Union[float, int]], optional): Substrate uptake values which are related to linear
            growth (no overflow-like metabolic phenotype). Defaults to np.arange(-4, 0, 1).

    Raises:
        KeyError: If no file path is provided and default file is not found.

    Returns:
        None: Writes the updated PAM data to an Excel file the Result/1_preprocessing directory.

    Notes:
        - ParameterDict = dict[Literal['slope', 'intercept'], float]
        - Output file name is in the following format:
         'Results / 1_preprocessing / proteinAllocationModel_iML1515_EnzymaticData_<yymmdd>.xlsx')
    """
    pam_parameter_information, output_file_path = _get_pam_parameter_information_from_excel(pam_data_file)

    slope_id, intercept_id = ('tps_mu', 'tps_0') if sector_id == 'Translational' else ('ups_mu', 'ups_0')


    sector_vs_growthrate = pd.Series({
        'id_list':biomass_rxn,
        slope_id:param_vs_growth['slope'],
        intercept_id:param_vs_growth['intercept'],
        'mol_mass':pam_parameter_information[sector_id].set_index('Parameter').loc['mol_mass', 'Value'],
        'substrate_range': list(substrate_range)
    })

    sector_vs_glucose = pd.Series({
        'id_list':lin_rxn_id,
        slope_id:param_vs_lin_rxn['slope'],
        intercept_id:param_vs_lin_rxn['intercept'],
        'mol_mass':pam_parameter_information['Translational'].set_index('Parameter').loc['mol_mass', 'Value'],
                })

    substrate_range_row = {
        'Parameter':'substrate_range',
        'Value':[],
        'Unit':'mmol/gDW/h',
        'Description': 'Range of values of susbstrate uptake which are associated with linear growth regime (no fermentation)'
    }

    with pd.ExcelWriter(output_file_path) as writer:
        for sheet, df in pam_parameter_information.items():
            if sheet == sector_id:
                df = pd.concat([df,pd.DataFrame([substrate_range_row])])
                df = df.set_index('Parameter')
                df['Value_for_growth'] = sector_vs_growthrate
                df['Value'] = sector_vs_glucose
                #change order for better readability
                df = df.reset_index()[['Parameter','Value','Value_for_growth','Unit','Description']]
            df.to_excel(writer, sheet_name=sheet, index = False)


def _get_pam_parameter_information_from_excel(pam_data_file: str)-> Tuple[dict[str, pd.DataFrame], str]:
    """Loads PAM parameter information from an Excel file.

       If a file matching the current date exists in the default Results directory,
       it will be used instead of the provided `pam_data_file`.

       Args:
           pam_data_file (str): Path to the PAM Excel file. If None and no default file exists, raises an error.

       Raises:
           KeyError: If no valid Excel file path is available.

       Returns:
           Tuple[dict[str, pd.DataFrame], str]:
               - Dictionary with sheet names as keys and DataFrames as values.
               - Path to the Excel file used.
       """
    current_date = datetime.datetime.now().strftime('%y%m%d')
    output_file_path = os.path.join('Results', '1_preprocessing',
                                    f'proteinAllocationModel_iML1515_EnzymaticData_{current_date}.xlsx')


    if os.path.isfile(output_file_path): pam_data_file = output_file_path
    elif pam_data_file is None: raise KeyError('Please provide a path to an excel file with information about the sector parameters')
    pam_parameter_information = pd.read_excel(pam_data_file, sheet_name = None)

    return pam_parameter_information, output_file_path


def set_up_sector_config(pam_info_file: str,
                         sectors_not_related_to_growth: List[str]
                         ) -> Dict[str,SectorConfig]:
    """Initializes sector configuration objects from PAM parameter Excel file.

    This function reads the PAM parameter Excel file and sets up configuration objects
    for each sector not linearly related to the growth rate, but to another linear reaction
    (such as substrate uptake rate). This allows the relation to be extrapolated to other
    conditions/substrates It extracts the slope, intercept, and substrate range associated
    with growth for each specified sector.

    Args:
        pam_info_file (str): Path to the Excel file containing PAM parameter data.
        sectors_not_related_to_growth (List[str]): List of sector identifiers
            (e.g., 'TranslationalProteinSector', 'UnusedEnzymeSector') that are not
            linearly related to the growth rate, but to for example the substrate uptake rate

    Returns:
        Dict[str, SectorConfig]: Dict of the parsed SectorConfig typed dicts with the corresponding model sector name
            as keys
    """

    sectors_to_sheets = {
        'TranslationalProteinSector': 'Translational',
        'UnusedEnzymeSector': 'UnusedEnzyme'
    }

    sector_configs = {}
    for sector_id in sectors_not_related_to_growth:
        if sector_id not in sectors_to_sheets:
            raise ValueError(f"Unknown sector: {sector_id}. Must be one of {list(sectors_to_sheets.keys())}.")

        sector_to_growth = pd.read_excel(
            pam_info_file,
            sheet_name=sectors_to_sheets[sector_id]
        ).set_index('Parameter')[['Value_for_growth']]

        sector_configs[sector_id] = SectorConfig(
            sectorname=sector_id,
            slope=sector_to_growth.filter(like='_mu', axis=0)['Value_for_growth'].iloc[0],
            intercept=sector_to_growth.filter(like='_0', axis=0)['Value_for_growth'].iloc[0],
            #convert string to list with ast.literal_eval
            substrate_range= ast.literal_eval(sector_to_growth.loc['substrate_range', 'Value_for_growth'])
        )
    return sector_configs