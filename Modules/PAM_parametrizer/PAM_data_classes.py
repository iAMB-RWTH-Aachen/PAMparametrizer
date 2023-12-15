from dataclasses import dataclass
import pandas as pd
from PAModelpy.configuration import Config

@dataclass
class ValidationData:
    valid_data_df: pd.DataFrame
    reactions_to_validate : [str]
    biomass_reaction_extension : str = 'BIOMASS'
    exchange_reaction_extension: str = 'EX'
    reactions_to_plot = [Config.ACETATE_EXCRETION_RXNID, Config.CO2_EXHANGE_RXNID, Config.OXYGEN_UPTAKE_RXNID, Config.BIOMASS_REACTION]


    @property
    def growth_rates(self) -> list:
        growth_rates = [data for data in self.valid_data_df.columns
                        if data.split('_')[0] == self.biomass_reaction_extension]
        return growth_rates

    @property
    def reactions_to_validate(self) -> list:
        reactions_to_validate = [data for data in self.valid_data_df.columns
                                 if data.split('_')[0] == self.exchange_reaction_extension]
        return reactions_to_validate

@dataclass
class HyperParameters:
    threshold_error = 0.1
    threshold_iteration = 100
    number_of_kcats_to_mutate = 5
    number_of_bins = 5
    bin_resolution = 5

@dataclass
class ParametrizationResults:
    error_df: pd.DataFrame
    fac_df: pd.DataFrame
    bins_to_change: pd.DataFrame
    sensitive_enzymes: pd.DataFrame = pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id'])
    fluxes = []

    def initiate_result_dfs(self, enzyme_ids, reactions_to_validate, growth_rate) -> None:
        self.error_df = pd.DataFrame(columns=['bin', 'substrate'] + reactions_to_validate + growth_rate)
        self.fac_df = pd.DataFrame(columns=['bin', 'substrate'] + enzyme_ids)
        self.initiate_bins_to_change()

    def initiate_bins_to_change(self) -> None:
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])