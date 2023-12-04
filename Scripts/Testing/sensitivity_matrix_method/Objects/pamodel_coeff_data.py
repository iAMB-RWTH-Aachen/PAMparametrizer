from dataclasses import dataclass
import warnings
import gurobipy as gp
import pandas as pd

from PAModelpy.PAModel import PAModel

@dataclass
class PAModelCoeffData:
    pamodel : PAModel

    @property
    def max_protein_fraction(self) -> float:
        return self.pamodel.constraints[self.pamodel.TOTAL_PROTEIN_CONSTRAINT_ID].ub
    @property
    def sum_of_flux_csc(self)-> float:
        if self.pamodel.solver.status!= 'optimal':
            warnings.warn('PAModel is not optimal, cannot return sum of flux CSCs')
            return 0.0
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
        calculates sum of enzymes in mmol/gcdw/h
        :return: sum of enzymes
        """
        # self.pamodel.calculate_sum_of_enzymes() cannot be used because this returns the sum of enzymes in mg/gcdw/h
        # coefficients needed for calculating fraction of utilized enzymes should be in mmol/gcdw/h
        return sum([enz.concentration for enz in self.pamodel.enzyme_variables])
    @property
    def sum_flux_csc_normalized_proteome_csc(self) -> float:
        return self.sum_of_flux_csc + (1-self.sum_of_enzymes/self.max_protein_fraction)*self.proteome_csc

    @property
    def enzyme_sensitivity_coefficients(self) -> list:
        if self.pamodel.solver.status != 'optimal':
            warnings.warn('PAModel is not optimal, cannot return ESCs')
            return []
        return self.pamodel.enzyme_sensitivity_coefficients.coefficient.to_list()

    @property
    def rxn2enzyme_from_esc(self) -> list[str,str]:
        if self.pamodel.solver.status!= 'optimal':
            warnings.warn('PAModel is not optimal, cannot return ESCs related rxn2enzymes')
            return []
        enzymes = self.pamodel.enzyme_sensitivity_coefficients.enzyme_id.to_list()
        rxns = self.pamodel.enzyme_sensitivity_coefficients.rxn_id.to_list()
        rxns2enzymes = [f'{rxn}_{enz}' for rxn, enz in zip(rxns, enzymes)]
        return rxns2enzymes

    def get_top_n_esc(self, n:int = 3)-> pd.DataFrame:
        esc_df = self.pamodel.enzyme_sensitivity_coefficients.copy()
        # ordered_esc = esc_df.sort_values(by='coefficient', ascending = True)
        ordered_esc = esc_df[esc_df.rxn_id.isin(['R3', 'R4', 'R5'])]
        return ordered_esc#.iloc[:n+1,:]


