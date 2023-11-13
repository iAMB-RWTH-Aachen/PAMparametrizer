import os
import pandas as pd
from dataclasses import dataclass
import warnings
import gurobipy as gp

from PAM_Parametrization.Scripts.pam_generation import setup_toy_pam

from PAModelpy.PAModel import PAModel

# Reference data from toy model simulations:
DATA_DIR = os.path.join(os.getcwd(), 'Data')
RESULT_DF_FILE = os.path.join(DATA_DIR, 'toy_model_simulations_ga.csv')
LMBDA = 0.5

@dataclass
class PAModelCoeffData:
    pamodel : PAModel

    @property
    def max_protein_fraction(self) -> float:
        return self.pamodel.constraints[self.pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub
    @property
    def sum_of_flux_csc(self)-> float:
        if self.pamodel!= 'optimal':
            warnings.warn('PAModel is not optimal, cannot return sum of flux CSCs')
            return
        return sum(self.pamodel.capacity_sensitivity_coefficients[
                       self.pamodel.capacity_sensitivity_coefficients['constraint'].isin(['flux_ub', 'flux_lb'])
                   ].coefficient.to_list())
    @property
    def proteome_csc(self) -> float:
        return sum(self.pamodel.capacity_sensitivity_coefficients[
                       self.pamodel.capacity_sensitivity_coefficients['constraint'] == 'proteome'
                   ].coefficient.to_list())
    @property
    def sum_of_enzymes(self)->float:
        """
        calculates sum of enzymes in mg/gcdw/h
        :return: sum of enzymes
        """
        return self.pamodel.calculate_sum_of_enzymes()
    @property
    def sum_flux_csc_normalized_proteome_csc(self) -> float:
        return self.sum_FCSC - (1-self.max_protein_fraction)*self.proteome_csc

    @property
    def enzyme_sensitivity_coefficients(self) -> list:
        if self.pamodel!= 'optimal':
            warnings.warn('PAModel is not optimal, cannot return ESCs')
            return
        return self.pamodel.enzyme_sensitivity_coefficients.coefficients.list()

    @property
    def rxn2enzyme_from_esc(self):
        if self.pamodel!= 'optimal':
            warnings.warn('PAModel is not optimal, cannot return ESCs related rxn2enzymes')
            return
        enzymes = self.pamodel.enzyme_sensitivity_coefficients.enzyme_id.list()
        rxns = self.pamodel.enzyme_sensitivity_coefficients.rxn_id.list()
        rxns2enzymes = [f'{rxn}_{enz}' for rxn, enz in zip(rxns, enzymes)]
        return rxns2enzymes

@dataclass
class KcatLinearModel:
    esc_matrix: [list, list]
    right_hand_side: list
    identifiers: [str, str]
    lmbda: float

    @property
    def gurobi_model(self):
        model = gp.Model()
        objective_expression = gp.LinExpr()
        #loop over all simulations
        for esc_simulation in self.esc_matrix:
        #add variables and constraints
            for id, esc, rhs in zip(self.identifiers, esc_simulation, self.right_hand_side):
                var = model.addVar(name= id)
                model.addConstr(esc*var == rhs, name=id)
                #make objective expression absolute by adding to other constraintsm, see https://lpsolve.sourceforge.net/5.1/absolute.htm
                abs_var = self.create_absolute_variable(model,gp.LinExpr(var), id+'_var')
                abs_error = self.create_absolute_variable(model, gp.LinExpr(esc*var - rhs), id+'_error')
                #lasso regression term
                objective_expression.addTerms(self.lmbda, abs_var)
                #error termo
                objective_expression.addTerms(1.0, abs_error)

        model.setObjective(objective_expression, gp.GRB.MINIMIZE)
        return model

    def create_absolute_variable(model:gp.Model, expression: gp.LinExpr, name:str)-> gp.Var:
        # make objective expression absolute by adding to other constraintsm, see https://lpsolve.sourceforge.net/5.1/absolute.htm
        abs_var = model.addVar(name)
        model.addConstr(expression <= abs_var, name= id+'_pos')
        model.addConstr(-expression <= abs_var, name= id+'_neg')
        return abs_var

def calculate_rhs(pamodeldata:PAModelCoeffData , growth_rate: float):
    error = pamodeldata.pamodel.objective.value-growth_rate
    right_hand_side = [error(1-pamodeldata.sum_flux_csc_normalized_proteome_csc)]
    return right_hand_side

def set_up_kcat_linear_model(pamodel:PAModel, substrate_rates:list, growth_rates: list):
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

    return KcatLinearModel(esc_matrix, right_hand_side, pamodel_data.rxn2enzyme_from_esc, LMBDA).gurobi_model


if __name__ == '__main__':
    valid_data_df = pd.read_csv(RESULT_DF_FILE)
    #start kcat: [1, 0.5, 1, 0.5 ,0.45, 1.5]
    #aimed end kcat: [1, 0.5, 5, 0.1, 0.25, 1.5]

    toy_pamodel = setup_toy_pam()

    substr_uptake_rates = valid_data_df['R1_ub'].list()
    growth_rates = valid_data_df['R7'].list()

    set_up_kcat_linear_model(toy_pamodel, substr_uptake_rates, growth_rates)
