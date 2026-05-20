import pandas as pd
import numpy as np
from pathlib import Path

BRENDA_PARAMETER_FILE_NAME = Path('Data/Databases/BRENDA_Kegg_Bareven.xlsx')
PARAMETER_FILE_NAME = Path('Results/1_preprocessing/proteinAllocationModel_iML1515_EnzymaticData_250912.xlsx')
RANDOM_PARAMETERS_FILE_NAME = Path('Results/1_preprocessing/proteinAllocationModel_iML1515_EnzymaticData_random.xlsx')
DIFFUSION_LIMIT=1e6

if __name__ == '__main__':
    np.random.seed(0)
    brenda_kegg_parameters = pd.read_excel(BRENDA_PARAMETER_FILE_NAME, sheet_name='KineticTable')


    param_files = pd.read_excel(PARAMETER_FILE_NAME, sheet_name=None)
    kcat_data = param_files['ActiveEnzymes']
    kcat_data['kcat_values'] = [val if val<DIFFUSION_LIMIT else DIFFUSION_LIMIT
                                for val in np.random.lognormal(
            np.log10(brenda_kegg_parameters['kcat (1/sec)'].mean(skipna=True)),
            sigma = np.log10(brenda_kegg_parameters['kcat (1/sec)'].std(skipna=True)),
            size = kcat_data.shape[0]
        )
                                ]
    param_files['ActiveEnzymes'] = kcat_data

    with pd.ExcelWriter(RANDOM_PARAMETERS_FILE_NAME) as writer:
        for sheetname, df in param_files.items():
            df.to_excel(writer, sheet_name=sheetname, index=False)
