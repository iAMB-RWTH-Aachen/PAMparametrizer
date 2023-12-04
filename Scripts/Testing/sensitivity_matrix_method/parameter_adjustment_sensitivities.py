import os
import pandas as pd

from PAM_Parametrization.Scripts.pam_generation import setup_toy_pam

from PAModelpy.PAModel import PAModel
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.pamodel_coeff_data import PAModelCoeffData
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.kcat_linear_model import KcatLinearModel


# Reference data from toy model simulations:
DATA_DIR = os.path.join(os.path.split(os.getcwd())[0], 'Data')
RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
LMBDA = 0

def calculate_rhs(pamodeldata:PAModelCoeffData , growth_rate: float):
    error = pamodeldata.pamodel.objective.value-growth_rate
    right_hand_side = error*(1-pamodeldata.sum_flux_csc_normalized_proteome_csc)
    return right_hand_side

def set_up_kcat_linear_model(pamodel:PAModel,valid_data: pd.DataFrame, substrate_rates:list, growth_rates: list):
    esc_matrix: [list(), list()] = list()
    right_hand_side = list()
    for substrate, growth_rate in zip(substrate_rates, growth_rates):
        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', substrate, 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
            pamodel_data = PAModelCoeffData(pamodel=pamodel)
            right_hand_side += [calculate_rhs(pamodel_data, growth_rate)]
            esc_matrix += [pamodel_data.enzyme_sensitivity_coefficients]

    return KcatLinearModel(substrate_rates, esc_matrix, right_hand_side, pamodel_data.rxn2enzyme_from_esc, valid_data, LMBDA).gurobi_model


if __name__ == '__main__':
    valid_data_df = pd.read_csv(RESULT_DF_FILE)
    #start kcat: [1, 0.5, 1, 0.5 ,0.45, 1.5]
    #aimed end kcat: [1, 0.5, 5, 0.1, 0.25, 1.5]

    toy_pamodel = setup_toy_pam()

    substr_uptake_rates = valid_data_df['R1_ub'].to_list()
    growth_rates = valid_data_df['R7'].to_list()
    colnames = list()
    for colname in valid_data_df.columns:
        if colname == 'R7':
            colnames.append('R6')
        elif colname == 'R8':
            colnames.append('R4')
        elif colname == 'R9':
            colnames.append('R3')
        else:
            colnames.append(colname)

    valid_data_df.columns = colnames

    linear_model = set_up_kcat_linear_model(toy_pamodel, valid_data_df,substr_uptake_rates, growth_rates)
    linear_model.optimize()
    linear_model.write('kcat_linear_model.lp')
    for var in linear_model.getVars():
        print(var.VarName, var.X)
