import os
import pandas as pd
import numpy as np
from typing import Literal, Tuple, Iterable, Union
import datetime

from Modules.utils.sector_config_functions import ParameterDict


#save to excel
def save_sector_information_to_excel(
        param_vs_lin_rxn: ParameterDict,
        param_vs_growth: ParameterDict,
        biomass_rxn:str,
        lin_rxn_id:str,
        sector_id:Literal['Translational', 'UnusedEnzyme'],
        pam_data_file:str = None,
        substrate_range: Iterable[Union[float, int]] = np.arange(-4,0,1)
)-> None:
    pam_parameter_information, output_file_path = _get_pam_parameter_information_from_excel(pam_data_file)
    slope_id, intercept_id = 'tps_mu','tps_0'
    if sector_id != 'Translational': slope_id, intercept_id = 'ups_mu','ups_0'

    sector_vs_growthrate = pd.Series({
        'id_list':biomass_rxn,
        slope_id:param_vs_growth['slope'],
        intercept_id:param_vs_growth['intercept'],
        'mol_mass':pam_parameter_information[sector_id].set_index('Parameter').loc['mol_mass', 'Value'],
        'substrate_range': substrate_range
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
            if sheet == 'Translational':
                df = pd.concat([df,pd.Series(substrate_range_row)], ignore_index =True)
                df = df.set_index('Parameter')
                df['Value_for_growth'] = sector_vs_growthrate
                df['Value'] = sector_vs_glucose
                #change order for better readability
                df = df.reset_index()[['Parameter','Value','Value_for_growth','Unit','Description']]
            df.to_excel(writer, sheet_name=sheet, index = False)


def _get_pam_parameter_information_from_excel(pam_data_file: str)-> Tuple[dict[str, pd.DataFrame], str]:
    current_date = datetime.datetime.now()
    # Format the date as yymmdd
    formatted_date = current_date.strftime('%y%m%d')
    output_file_path = os.path.join('Results',
                                '1_preprocessing',
                                f'proteinAllocationModel_iML1515_EnzymaticData_{formatted_date}.xlsx')
    if os.path.isfile(output_file_path): pam_data_file = output_file_path
    elif pam_data_file is None: raise KeyError('Please provide a path to an excel file with information about the sector parameters')
    pam_parameter_information = pd.read_excel(pam_data_file, sheet_name = None)

    return pam_parameter_information, output_file_path