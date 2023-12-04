import gurobipy as gp
from gurobipy import GRB
import pandas as pd
import numpy as np

from PAModelpy.PAModel import PAModel
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.pamodel_coeff_data import PAModelCoeffData
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.kcat_linear_model import KcatLinearModel

def lasso_regression(X, y, lambda_value):
    m, n = X.shape

    # Create a Gurobi model
    model = gp.Model()

    # Add variables for the regression coefficients and their absolute values
    beta = model.addVars(n, name="beta")
    abs_beta = model.addVars(n, name="abs_beta")

    # Add constraints to link beta and abs_beta
    model.addConstrs((abs_beta[j] >= beta[j] for j in range(n)))
    model.addConstrs((abs_beta[j] >= -beta[j] for j in range(n)))

    # Add objective function
    objective = 0.5 * gp.quicksum((y[i] - gp.quicksum(X[i, j] * beta[j] for j in range(n)))**2 for i in range(m))
    objective += lambda_value * gp.quicksum(abs_beta[j] for j in range(n))

    model.setObjective(objective, GRB.MINIMIZE)

    # Optimize the model
    model.optimize()

    if model.status == GRB.OPTIMAL:
        # Get the optimal values of the regression coefficients
        beta_values = [beta[j].x for j in range(n)]
        return beta_values
    else:
        print("Optimization failed. Status:", model.status)
        return None
def calculate_rhs(pamodeldata:PAModelCoeffData, predicted_growth: float , growth_rate: float):
    error = predicted_growth-growth_rate
    right_hand_side = error*(1-pamodeldata.sum_flux_csc_normalized_proteome_csc)
    return right_hand_side
def set_up_kcat_lasso_lp(pamodel:PAModel,valid_data: pd.DataFrame, substrate_id:str, biomass_id: str ):
    esc_matrix: [np.array, np.array] = list()
    right_hand_side = list()

    for index, row in valid_data.iterrows():
        pamodel.objective = {pamodel.reactions.get_by_id(biomass_id): 1}
        substrate = row[substrate_id]
        growth_rate = row[biomass_id]
        conditions = {rxn: value for rxn, value in row[1:].items() if (rxn != substrate_id and rxn!= substrate_id)}

        pamodel.change_reaction_bounds(rxn_id='R1',
                                       lower_bound=0, upper_bound=substrate)
        print('Running simulations with ', round(substrate, 3), 'mmol/g_cdw/h of substrate going into the system')
        pamodel.optimize()
        if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
            pamodel_data = PAModelCoeffData(pamodel=pamodel)
            right_hand_side += [calculate_rhs(pamodel_data, pamodel_data.pamodel.objective.value, growth_rate)]
            esc_matrix += [np.array(pamodel_data.enzyme_sensitivity_coefficients)]#[np.array([pamodel_data.enzyme_sensitivity_coefficients[i] for i in pamodel_data.get_top_n_esc(n=3).index])]

        # pamodel.change_reaction_bounds(biomass_id, growth_rate, growth_rate)
        for rxn, value in conditions.items():
            pamodel.objective = {pamodel.reactions.get_by_id(rxn):-1}
            model_rxn = pamodel.reactions.get_by_id(rxn)
            prev_lb, prev_ub = model_rxn.lower_bound, model_rxn.upper_bound
            # pamodel.change_reaction_bounds(rxn, value, value)

            pamodel.optimize()
            if pamodel.solver.status == 'optimal' and pamodel.objective.value > 0:
                pamodel_data = PAModelCoeffData(pamodel=pamodel)
                right_hand_side += [calculate_rhs(pamodel_data,pamodel.reactions.get_by_id(biomass_id).flux, growth_rate)]
                esc_matrix += [np.array(pamodel_data.enzyme_sensitivity_coefficients)]#[np.array([pamodel_data.enzyme_sensitivity_coefficients[i] for i in pamodel_data.get_top_n_esc(n=3).index])]
            pamodel.change_reaction_bounds(rxn, prev_lb, prev_ub)
        # pamodel.change_reaction_bounds(biomass_id, 0, 100)

    return lasso_regression(X = np.array(esc_matrix), y=right_hand_side, lambda_value=0.1)
    # return KcatLinearModel(substrate_rates, esc_matrix, right_hand_side, pamodel_data.rxn2enzyme_from_esc, valid_data, LMBDA).gurobi_model


# Example usage:
import numpy as np
import os

from PAM_Parametrization.Scripts.pam_generation import setup_toy_pam

from PAModelpy.PAModel import PAModel
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.pamodel_coeff_data import PAModelCoeffData
from PAM_Parametrization.Scripts.Testing.sensitivity_matrix_method.Objects.kcat_linear_model import KcatLinearModel

# Reference data from toy model simulations:
DATA_DIR = os.path.join(os.path.split(os.getcwd())[0], 'Data')
RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
LMBDA = 0.5

valid_data_df = pd.read_csv(RESULT_DF_FILE)
#start kcat: [1, 0.5, 1, 0.5 ,0.45, 1.5]
#aimed end kcat: [1, 0.5, 5, 0.1, 0.25, 1.5]
kcat_start = [i/5 for i in [1, 0.5, 1, 0.5 ,0.45, 1.5]]

toy_pamodel = setup_toy_pam(kcat_fwd=kcat_start)

substr_uptake_rates = valid_data_df['R1_ub'].to_list()
growth_rates = valid_data_df['R7'].to_list()
# colnames = list()
# for colname in valid_data_df.columns:
#     if colname == 'R7':
#         colnames.append('R6')
#     elif colname == 'R8':
#         colnames.append('R4')
#     elif colname == 'R9':
#         colnames.append('R3')
#     else:
#         colnames.append(colname)

# valid_data_df.columns = colnames

result = set_up_kcat_lasso_lp(toy_pamodel, valid_data_df,'R1_ub', 'R7')
print("Estimated kcat difference:", result)

# # Generate some random data for demonstration purposes
# np.random.seed(42)
# X = np.random.rand(100, 5)
# true_beta = np.array([2, -3, 1.5, 0, 0])
# y = X.dot(true_beta) + np.random.normal(scale=0.1, size=100)
#
# # Set the regularization parameter lambda
# lambda_value = 0.1
#
# # Solve the Lasso regression problem
# result = lasso_regression(X, y, lambda_value)
#
# # Print the results
# print("True beta:", true_beta)
# print("Estimated beta:", result)