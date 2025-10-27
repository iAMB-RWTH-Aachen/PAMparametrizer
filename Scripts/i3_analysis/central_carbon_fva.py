from cobra.flux_analysis import flux_variability_analysis
import os
import pandas as pd

from cobra.io.sbml import read_sbml_model
from PAModelpy.utils import set_up_pam

if __name__ == '__main__':
    NUM_MODELS = 10
    FVA_RESULT_FILE = os.path.join('Results', '3_analysis','fva_central_carbon_metabolism.xlsx')
    models_to_check = {i: os.path.join('Results', '3_analysis', 'parameter_files',
                                    f'proteinAllocationModel_EnzymaticData_iML1515_{i}.xlsx') for i in range(1, NUM_MODELS+1)}
    models_to_check['GotEnzymes'] = (os.path.join('Results','2_parametrization', 'proteinAllocationModel_iML1515_EnzymaticData_multi.xlsx'))
    models_to_check = {name: set_up_pam(file, sensitivity=False) for name, file in models_to_check.items()}
    models_to_check['iML1515'] = read_sbml_model(os.path.join('Models', 'iML1515.xml'))

    for name, model in models_to_check.items():
        result = flux_variability_analysis(model = model,loopless=True)
        if not os.path.exists(FVA_RESULT_FILE): kwargs = {'mode':'w'}
        else:kwargs = {'mode':'a', 'if_sheet_exists':'replace'}
        with pd.ExcelWriter(FVA_RESULT_FILE, engine='openpyxl', **kwargs) as writer:
            result.to_excel(writer, sheet_name=name)


