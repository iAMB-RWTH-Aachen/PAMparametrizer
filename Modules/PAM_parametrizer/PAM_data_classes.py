from dataclasses import dataclass, field
import pandas as pd
from PAModelpy.configuration import Config
from typing import Union, Callable
from pathlib import Path
import random

from Modules.genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from Modules.genetic_algorithm_parametrization import GAPOUniform

@dataclass
class ValidationData:
    valid_data:dict
    sampled_valid_data: dict = None
    _reactions_to_validate : {str:str} = field(default_factory=dict)
    biomass_reaction_extension : str = 'BIOMASS'
    exchange_reaction_extension: str = 'EX'
    _reactions_to_plot = [Config.ACETATE_EXCRETION_RXNID, Config.CO2_EXHANGE_RXNID, Config.OXYGEN_UPTAKE_RXNID, Config.BIOMASS_REACTION]


    def _get_biomass_reactions(self) -> list:
        valid_data_df = random.choice(list(self.valid_data.values()))
        return [data for data in valid_data_df.columns
                        if data.split('_')[0] == self.biomass_reaction_extension]
    @property
    def biomass_reactions(self) -> list:
        return self._get_biomass_reactions()

    def _get_reactions_to_validate(self) -> dict:
        reactions_to_validate_dict = {}
        for substr_upt_rxn, valid_data_df in self.valid_data.items():
            reactions_to_validate_dict[substr_upt_rxn] = [data for data in self.valid_data_df.columns
                                 if data.split('_')[0] == self.exchange_reaction_extension]
        return reactions_to_validate_dict

    @property
    def reactions_to_validate(self) -> dict:
        return self._get_reactions_to_validate()

    @reactions_to_validate.setter
    def reactions_to_validate(self, reactions_to_validate:dict):
        self._reactions_to_validate = reactions_to_validate


@dataclass
class HyperParameters:
    threshold_error = 0.9
    threshold_iteration = 100
    number_of_kcats_to_mutate = 5
    number_of_bins = 5
    bin_resolution = 5
    bin_split_deviation_threshold = 0.20
    genetic_algorithm_filename_base = 'genetic_algorithm_run_'
    filename_extension = ''
    genetic_algorithm: Callable = GAPOUniform
    genetic_algorithm_hyperparams = {'mutation_probability':0.5, # probability with which an individual (solution) is mutated in a generation
                                'mutation_rate':0.5, # probability with which an attribute (e.g. gene) of an individual is mutated
                                'population_size':10, # number of individuals (solution) per population
                                'crossover_probability':0.8, # probability with which two indivduals/offsprings are crossed over
                                'number_generations':20, # number of consecutive generations per gene flow event
                                'number_gene_flow_events':2, # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
                                'init_attribute_probability':0.,
                                'fitness_class': 'Fitfun_params_uniform',
                                'processes': 2,  # number of parallel workers
                                'time_limit': 600,  # time limit in seconds
                                # 'error_weights': {'EX_ac_e':5}, # reaction which should have a different impact on the error calculation than other reactions
                                'folderpath_save': Path(r"Results"),  # path for saving results
                                'overwrite_intermediate_results': True,  # if true, saved intermediate results are overwritten
                                'print_progress': True # if True, progress of the genetic algorithm is printed
                                     }
    #change fitness function if Gaussian sampling is applied
    if isinstance(genetic_algorithm, GAPOGauss):
        genetic_algorithm_hyperparams['fitness_class'] = "Fitfun_params_gaussian"

@dataclass
class ParametrizationResults:
    substrate_uptake_reactions: list
    error: dict[str:pd.DataFrame] = field(default_factory=dict)
    esc_df: pd.DataFrame = None
    bins_to_change: pd.DataFrame = pd.DataFrame(columns=['bin', 'split', 'merge'])
    sensitive_enzymes: pd.DataFrame = pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id'])
    fluxes: dict[str:pd.DataFrame] = field(default_factory=dict)
    substrate_range: dict[str:list] = field(default_factory=dict)
    _color = 440154
    best_individuals = pd.DataFrame(columns = ['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error'])
    computational_time = pd.DataFrame(columns = ['run_id', 'time_s', 'time_h'])
    final_errors = pd.DataFrame(columns=['run_id', 'r_squared'])

    def initiate_result_dfs(self, reactions_to_validate, biomass_reaction) -> None:
        self.esc_df = pd.DataFrame(columns=['bin', 'substrate', 'enzyme_id', 'rxn_id'])
        for rxn in self.substrate_uptake_reactions:
            self.error[rxn] = pd.DataFrame(columns=['bin'] + reactions_to_validate[rxn])
            self.fluxes[rxn] = pd.DataFrame(columns= ['bin','substrate'] + reactions_to_validate[rxn])
        self.initiate_bins_to_change()

    def initiate_bins_to_change(self) -> None:
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])

    def add_fluxes(self, pamodel, bin_id: Union[float, int, str],
                   substrate_uptake_rate: Union[float, int], substrate_uptake_reaction:str,
                   fluxes_abs:bool = True):
        new_row = [bin_id, substrate_uptake_rate]+[0]*len(self.fluxes[substrate_uptake_reaction].columns[2:])
        self.fluxes[substrate_uptake_reaction].loc[len(self.fluxes[substrate_uptake_reaction])] = new_row
        for rxn in pamodel.reactions:
            if rxn.id in self.fluxes[substrate_uptake_reaction].columns:
                if fluxes_abs: flux = abs(rxn.flux)
                else: flux = rxn.flux
                self.fluxes[substrate_uptake_reaction].iloc[-1, self.fluxes[substrate_uptake_reaction].columns.get_loc(rxn.id)] = flux

    def add_fluxes_from_fluxdict(self, flux_dict:dict, bin_id: Union[float, int, str],
                                 substrate_uptake_rate: Union[float, int], substrate_uptake_reaction:str,
                                 fluxes_abs:bool = True):
        new_row = [bin_id, substrate_uptake_rate]+[0]*len(self.fluxes[substrate_uptake_reaction].columns[2:])
        self.fluxes[substrate_uptake_reaction].loc[len(self.fluxes[substrate_uptake_reaction])] = new_row
        for rxn_id, flux in flux_dict.items():
            if rxn_id in self.fluxes[substrate_uptake_reaction].columns:
                if fluxes_abs: flux = abs(flux)
                else: flux = flux
                self.fluxes[substrate_uptake_reaction].iloc[-1, self.fluxes[substrate_uptake_reaction].columns.get_loc(rxn_id)] = flux

    def add_enzyme_sensitivity_coefficients(self, enzyme_sensitivity_coeff: pd.DataFrame,
                                            bin_id:Union[float, int, str], substrate_uptake_rate: Union[float, int]):
        # we are only interested in the enzymes
        if not isinstance(bin_id, str): bin_id = float(bin_id)
        enzyme_sensitivity_coeff['bin'] = [bin_id]*len(enzyme_sensitivity_coeff)
        enzyme_sensitivity_coeff['substrate'] = [substrate_uptake_rate]*len(enzyme_sensitivity_coeff)
        self.esc_df = pd.concat([self.esc_df, enzyme_sensitivity_coeff])
        self.esc_df.reset_index(drop = True)

    def add_best_individuals(self, run_id: str, best_indiv_enz_rxn_kcat: list, ga_error:float) -> None:
        self.best_individuals.loc[len(self.best_individuals)] = [run_id] + best_indiv_enz_rxn_kcat + [ga_error]
    def add_computational_time(self, run_id: str, time_in_sec: float):
        self.computational_time.loc[len(self.computational_time)] = [run_id, time_in_sec, time_in_sec/3600]

    def add_final_error(self, run_id: str, final_error: float):
        self.final_errors.loc[len(self.final_errors)] = [run_id, final_error]