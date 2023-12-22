from dataclasses import dataclass, field
import pandas as pd
from PAModelpy.configuration import Config
from typing import Union, Callable
from pathlib import Path

from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from Modules.genetic_algorithm_parametrization import GAPOUniform

@dataclass
class ValidationData:
    valid_data_df: pd.DataFrame
    _reactions_to_validate : [str] = field(default_factory=list)
    biomass_reaction_extension : str = 'BIOMASS'
    exchange_reaction_extension: str = 'EX'
    _reactions_to_plot = [Config.ACETATE_EXCRETION_RXNID, Config.CO2_EXHANGE_RXNID, Config.OXYGEN_UPTAKE_RXNID, Config.BIOMASS_REACTION]


    def _get_biomass_reactions(self) -> list:
        return [data for data in self.valid_data_df.columns
                        if data.split('_')[0] == self.biomass_reaction_extension]
    @property
    def biomass_reactions(self) -> list:
        return self._get_biomass_reactions()

    def _get_reactions_to_validate(self) -> list:
        return [data for data in self.valid_data_df.columns
                                 if data.split('_')[0] == self.exchange_reaction_extension]

    @property
    def reactions_to_validate(self) -> list:
        return self._get_reactions_to_validate()

    @reactions_to_validate.setter
    def reactions_to_validate(self, reactions_to_validate:list):
        self._reactions_to_validate = reactions_to_validate



@dataclass
class HyperParameters:
    threshold_error = 0.1
    threshold_iteration = 100
    number_of_kcats_to_mutate = 5
    number_of_bins = 5
    bin_resolution = 5
    bin_split_deviation_threshold = 0.20
    genetic_algorithm_filename_base = 'genetic_algorithm_run_'
    genetic_algorithm: Callable = GAPOUniform
    genetic_algorithm_hyperparams = {'mutation_probability':0.5, # probability with which an individual (solution) is mutated in a generation
                                'mutation_rate':1, # probability with which an attribute (e.g. gene) of an individual is mutated
                                'population_size':10, # number of individuals (solution) per population
                                'crossover_probability':0.8, # probability with which two indivduals/offsprings are crossed over
                                'number_generations':20, # number of consecutive generations per gene flow event
                                'number_gene_flow_events':2, # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
                                'init_attribute_probability':0.,
                                'fitness_class': 'Fitfun_params_uniform',
                                'processes': 2,  # number of parallel workers
                                'time_limit': 600,  # time limit in seconds
                                'folderpath_save': Path(r"Results"),  # path for saving results
                                'overwrite_intermediate_results': True  # if true, saved intermediate results are overwritten
                                     }
    #change fitness function if Gaussian sampling is applied
    if isinstance(genetic_algorithm, GAPOGauss):
        genetic_algorithm_hyperparams['fitness_class'] = "Fitfun_params_gaussian"

@dataclass
class ParametrizationResults:
    error_df: pd.DataFrame = None
    esc_df: pd.DataFrame = None
    bins_to_change: pd.DataFrame = None
    sensitive_enzymes: pd.DataFrame = pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id'])
    fluxes_df = pd.DataFrame()
    substrate_range = []
    _color = 440154

    def initiate_result_dfs(self, reactions_to_validate, biomass_reaction) -> None:
        self.error_df = pd.DataFrame(columns=['bin']+ reactions_to_validate)
        self.esc_df = pd.DataFrame(columns=['bin', 'substrate', 'enzyme_id', 'rxn_id'])
        self.fluxes_df = pd.DataFrame(columns= ['bin','substrate'] + reactions_to_validate)
        self.initiate_bins_to_change()

    def initiate_bins_to_change(self) -> None:
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])

    def add_fluxes(self, pamodel, bin_id: Union[float, int], substrate_uptake_rate: Union[float, int]):
        fluxes_to_save = []
        for rxn in pamodel.reactions:
            if rxn.id in self.fluxes_df.columns:
                fluxes_to_save += [rxn.flux]
        self.fluxes_df.loc[len(self.fluxes_df)] = [bin_id, substrate_uptake_rate] + fluxes_to_save

    def add_enzyme_sensitivity_coefficients(self, enzyme_sensitivity_coeff: pd.DataFrame,
                                            bin_id:Union[float, int], substrate_uptake_rate: Union[float, int]):
        # we are only interested in the enzymes
        enzyme_sensitivity_coeff['bin'] = [float(bin_id)]*len(enzyme_sensitivity_coeff)
        enzyme_sensitivity_coeff['substrate'] = [substrate_uptake_rate]*len(enzyme_sensitivity_coeff)
        self.esc_df = pd.concat([self.esc_df, enzyme_sensitivity_coeff])
        self.esc_df.reset_index(drop = True)