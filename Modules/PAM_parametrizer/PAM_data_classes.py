from dataclasses import dataclass, field
import pandas as pd
from PAModelpy.configuration import Config
from PAModelpy import TransEnzymeSector, UnusedEnzymeSector, CustomSector
from typing import Union, Callable, TypedDict, Iterable, List
import os
from cobra import DictList

from ..genetic_algorithm_parametrization import GAPOGaussian as GAPOGauss
from ..genetic_algorithm_parametrization import GAPOUniform

class SectorConfig(TypedDict):
    """
    Dictionary to store the information on how to configure a protein sector for new carbon sources.
    Can be used as a 'bride' between different carbon sources by relating the sector to the growth rate

    Args:
        sectorname (str): the sector identifier for which the relation is defined as present in the PAModel
        slope (float): slope between protein fraction allocated to the sector and the growth rate [g/gCDW/h]
        intercept (float): the protein fraction allocated to the sector at zero growth [g/gCDW]
        substrate_range (iterable of floats or ints): substrate range which for all the substrates is not related to fermentative phenotypes
    """
    sectorname:str
    slope:float
    intercept:float
    substrate_range:Iterable[Union[float, int]]


@dataclass
class ValidationData:
    """
        Container for experimental validation data used in PAM parametrization.

        This class organizes experimental measurements (e.g., uptake rates, growth,
        exchange fluxes) into a standardized format that can be compared to model outputs.
        It provides helpers for extracting biomass and exchange reactions and allows
        explicit control of which reactions are validated.

        Attributes:
            valid_data (pd.DataFrame): Experimental validation data, where columns
                correspond to reactions and rows to conditions or timepoints.
            id (str): Substrate uptake reaction id.
            validation_range (list): Range of substrate uptake or conditions tested.
            sampled_valid_data (pd.DataFrame, optional): Subset of sampled data used in parametrization,
                automatically generated during the optimization procedure.
            inactive_exchanges (list): List of exchange reactions assumed inactive
                (e.g. exchanges which should not carry flux).
            sector_configs (dict): Optional sector configuration used for this substrate.
            _reactions_to_validate (list): Internal storage of reaction IDs to validate.
            biomass_reaction_extension (str): Prefix for biomass reactions (default "BIOMASS").
            exchange_reaction_extension (str): Prefix for exchange reactions (default "EX").
            _reactions_to_plot (list): Predefined reactions to include in diagnostic plots.

        Properties:
            biomass_reactions (list): Extracted biomass reactions from the validation data.
            reactions_to_validate (list): List of reactions considered for validation.

        """

    valid_data:pd.DataFrame
    id: str
    validation_range: list
    sampled_valid_data: pd.DataFrame = None
    inactive_exchanges:Union[List, None] = field(default_factory=list)
    sector_configs: dict = field(default_factory=dict)
    _reactions_to_validate : list = field(default_factory=list)
    biomass_reaction_extension : str = 'BIOMASS'
    exchange_reaction_extension: str = 'EX'
    _reactions_to_plot = [Config.ACETATE_EXCRETION_RXNID, Config.CO2_EXHANGE_RXNID, Config.OXYGEN_UPTAKE_RXNID, Config.BIOMASS_REACTION]

    def __post_init__(self):
        self._add_inactive_exchanges_to_validation_df()

    def _add_inactive_exchanges_to_validation_df(self) -> None:
        """Adds inactive exchange reactions to the validation dataframe with zero values."""
        for ex_rxn in self.inactive_exchanges:
            self.valid_data[ex_rxn] = 0

    def _get_biomass_reactions(self) -> list:
        """Extracts biomass reactions from the validation dataframe."""
        return [data for data in self.valid_data.columns
                        if data.split('_')[0] == self.biomass_reaction_extension]
    @property
    def biomass_reactions(self) -> list:
        return self._get_biomass_reactions()

    def _get_reactions_to_validate(self) -> list:
        """Returns the list of reactions to validate, either predefined or inferred."""
        if len(self._reactions_to_validate)==0:
            rxns2validate = [data for data in self.valid_data.columns if data[-3:]!="_ub"]#data.split('_')[0] == self.exchange_reaction_extension]
        else:
            rxns2validate = self._reactions_to_validate
        return rxns2validate

    @property
    def reactions_to_validate(self) -> list:
        return self._reactions_to_validate()

    @reactions_to_validate.setter
    def reactions_to_validate(self, reactions_to_validate:list):
        self._reactions_to_validate = reactions_to_validate


@dataclass
class HyperParameters:
    """
       Hyperparameters controlling the parametrization and genetic algorithm process.

       Defines convergence thresholds, number of bins, kcats to mutate, and GA-specific
       settings such as mutation rate, population size, and number of generations.

       Attributes:
           threshold_error (float): Allowed r² threshold for terminating the algorithm.
           threshold_iteration (int): Maximum number of iterations before stopping.
           threshold_convergence (float): Difference in r² tolerated for convergence.
           threshold_nmbr_convergence (int): Number of consecutive converged iterations required
            before increasing number of kcats to mutate.
           number_of_kcats_to_mutate (int): Number of kcat values mutated per GA iteration.
           number_of_bins (int): Number of bins used for discretization (if applicable).
           bin_resolution (int): Resolution of each bin (if applicable).
           bin_split_deviation_threshold (float): Threshold for splitting bins during optimization  (if applicable).
           genetic_algorithm_filename_base (str): Base filename for GA results.
           filename_extension (str): Extension added to output filenames.
           genetic_algorithm (Callable): Genetic algorithm implementation to use.
           genetic_algorithm_hyperparams (dict): Dictionary of GA hyperparameters.
       """

    threshold_error = 0.9
    threshold_iteration = 100
    threshold_convergence = 0.09 #what is the difference in r_squared allowed to call it converging?
    threshold_nmbr_convergence = 2 #how many iterations have to be converged?
    number_of_kcats_to_mutate = 5
    number_of_bins = 5
    bin_resolution = 5
    bin_split_deviation_threshold = 0.20
    genetic_algorithm_filename_base = 'genetic_algorithm_run_'
    filename_extension = ''
    genetic_algorithm: Callable = GAPOUniform
    genetic_algorithm_hyperparams = {'mutation_probability':0.9, # probability with which an individual (solution) is mutated in a generation
                                'mutation_rate':0.9, # probability with which an attribute (e.g. gene) of an individual is mutated
                                'population_size':5, # number of individuals (solution) per population
                                'crossover_probability':0.8, # probability with which two indivduals/offsprings are crossed over
                                'number_generations':10, # number of consecutive generations per gene flow event
                                'number_gene_flow_events':2, # number of gene flow events, i.e. merging of multiple
                                 # populations independently evolved on parallel workers
                                'init_attribute_probability':0.,
                                'fitness_class': 'Fitfun_params_uniform',
                                'processes': 2,  # number of parallel workers
                                'time_limit': 60000,  # time limit in seconds
                                'error_weights': {}, # reaction which should have a different impact on the error calculation than other reactions
                                'folderpath_save':os.path.join('Results', '2_parametrization'),  # path for saving results
                                'overwrite_intermediate_results': True,  # if true, saved intermediate results are overwritten
                                'print_progress': True # if True, progress of the genetic algorithm is printed
                                     }
    #change fitness function if Gaussian sampling is applied
    if isinstance(genetic_algorithm, GAPOGauss):
        genetic_algorithm_hyperparams['fitness_class'] = "Fitfun_params_gaussian"

@dataclass
class ParametrizationResults:
    """
        Stores and manages results from PAM parametrization runs.

        This class collects flux predictions, enzyme sensitivity coefficients, errors,
        best individuals, and computational times during parametrization. It provides
        helpers to add new results, track errors, and manage multiple substrates.

        Attributes:
            substrate_uptake_reactions (list): List of substrate uptake reaction IDs.
            esc_df (pd.DataFrame): Enzyme sensitivity coefficients across bins.
            sensitive_enzymes (pd.DataFrame): Aggregated sensitivities per enzyme.
            flux_results (DictList): FluxResults objects for each substrate.
            _color (int): Default color code for plotting.
            best_individuals (pd.DataFrame): Records of best individuals from GA runs.
            computational_time (pd.DataFrame): Timing of runs in seconds and hours.
            final_errors (pd.DataFrame): Final error values (e.g., r²) per run.

        """
    substrate_uptake_reactions:list
    esc_df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=['bin', 'split', 'merge']))
    sensitive_enzymes: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=['bin', 'mean_sensitivity', 'enzyme_id']))
    flux_results: DictList = None
    _color: int = 440154
    best_individuals: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=['run_id', 'enzyme_id', 'direction', 'rxn_id', 'kcat[s-1]', 'ga_error']))
    computational_time: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=['run_id', 'time_s', 'time_h']))
    final_errors: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=['run_id', 'r_squared']))

    def initiate_flux_results(self):
        """Initializes FluxResults for each substrate uptake reaction."""
        self.flux_results = DictList([FluxResults(substr_rxn) for substr_rxn in self.substrate_uptake_reactions])

    def initiate_result_dfs(self, reactions_to_validate: dict) -> None:
        """Creates fresh results DataFrames for fluxes and sensitivities."""
        self.esc_df = pd.DataFrame(columns=['bin', 'substrate', 'enzyme_id', 'rxn_id'])
        self.flux_results = DictList([FluxResults(substr_rxn) for substr_rxn in self.substrate_uptake_reactions])
        for flux_result in self.flux_results:
            flux_result.initiate_result_dfs(reactions_to_validate[flux_result.id])
        self.initiate_bins_to_change()

    def initiate_bins_to_change(self) -> None:
        """Initializes a DataFrame for bin split/merge operations."""
        self.bins_to_change = pd.DataFrame(columns=['bin', 'split', 'merge'])

    def add_new_substrate_source(self, substrate_uptake_id:str, reactions_to_validate:list) -> None:
        """Adds a new substrate and initializes its results container."""
        self.substrate_uptake_reactions.append(substrate_uptake_id)
        flux_results = FluxResults(substrate_uptake_id)
        flux_results.initiate_result_dfs(reactions_to_validate)
        self.flux_results.append(flux_results)

    def add_fluxes(self, pamodel, bin_id: Union[float, int, str], substrate_reaction_id: str,
                   substrate_uptake_rate: Union[float, int], fluxes_abs:bool = True):
        """Adds fluxes from a PAM run into the results."""
        self.flux_results.get_by_id(substrate_reaction_id).add_fluxes(pamodel, bin_id,
                   substrate_uptake_rate, fluxes_abs)

    def add_fluxes_from_fluxdict(self, flux_dict:dict, bin_id: Union[float, int, str], substrate_reaction_id: str,
                                 substrate_uptake_rate: Union[float, int], fluxes_abs:bool = True):
        """Adds fluxes directly from a flux dictionary."""
        self.flux_results.get_by_id(substrate_reaction_id).add_fluxes_from_fluxdict(flux_dict, bin_id,
                                                                                    substrate_uptake_rate, fluxes_abs)

    def add_enzyme_sensitivity_coefficients(self, enzyme_sensitivity_coeff: pd.DataFrame,
                                            bin_id:Union[float, int, str], substrate_uptake_rate: Union[float, int]):
        """Appends enzyme sensitivity coefficients to the results DataFrame."""
        # we are only interested in the enzymes
        if not isinstance(bin_id, str): bin_id = float(bin_id)
        enzyme_sensitivity_coeff['bin'] = [bin_id]*len(enzyme_sensitivity_coeff)
        enzyme_sensitivity_coeff['substrate'] = [substrate_uptake_rate]*len(enzyme_sensitivity_coeff)
        self.esc_df = pd.concat([self.esc_df, enzyme_sensitivity_coeff])
        self.esc_df.reset_index(drop = True)

    def add_error_to_error_df(self, substrate_uptake_id: str,bin_id:Union[int, str], error: float):
        """ Records error for a given substrate and bin."""
        self.flux_results.get_by_id(substrate_uptake_id).add_error_to_error_df(bin_id, error)

    def add_best_individuals(self, run_id: str, best_indiv_enz_rxn_kcat: list, ga_error:float) -> None:
        """Logs the best-performing GA individuals for a run."""
        self.best_individuals.loc[len(self.best_individuals)] = [run_id] + best_indiv_enz_rxn_kcat + [ga_error]
        self.best_individuals = self.best_individuals.drop_duplicates(['run_id', 'enzyme_id', 'direction', 'rxn_id'],
                                                                      keep = 'last')

    def add_computational_time(self, run_id: str, time_in_sec: float):
        """Tracks computational time for a run."""
        self.computational_time.loc[len(self.computational_time)] = [run_id, time_in_sec, time_in_sec/3600]

    def add_final_error(self, run_id: str, final_error: float):
        """Records the final error for a run."""
        self.final_errors.loc[len(self.final_errors)] = [run_id, final_error]

    def remove_simulations_from_flux_df(self, substrate_uptake_id, bin_id):
        """Removes specific simulations from stored flux results."""
        self.flux_results.get_by_id(substrate_uptake_id).remove_simulations_from_flux_df(bin_id)

@dataclass
class FluxResults:
    """
    Stores flux results and errors for a single substrate uptake reaction.

    Each FluxResults object corresponds to one substrate uptake ID and
    contains DataFrames of predicted fluxes and errors across bins.

    Attributes:
        id (str): Substrate uptake reaction ID.
        error_df (pd.DataFrame): DataFrame of errors per bin.
        fluxes_df (pd.DataFrame): DataFrame of fluxes per bin and reaction.
        substrate_range (list): Range of substrate uptake values tested.
    """

    id: str
    error_df: pd.DataFrame = None
    fluxes_df: pd.DataFrame = None
    substrate_range: list = field(default_factory=list)

    def initiate_result_dfs(self, reactions_to_validate: list) -> None:
        """Initializes DataFrames for fluxes and errors."""
        self.error_df = pd.DataFrame(columns=['bin'] + reactions_to_validate)
        self.fluxes_df = pd.DataFrame(columns=['bin', 'substrate'] + reactions_to_validate)

    def add_fluxes(self, pamodel, bin_id: Union[float, int, str],
                   substrate_uptake_rate: Union[float, int], fluxes_abs:bool = True):
        """Adds fluxes from a PAM run."""
        new_row = [bin_id, substrate_uptake_rate]+[0]*len(self.fluxes_df.columns[2:])
        self.fluxes_df.loc[len(self.fluxes_df)] = new_row
        for rxn in pamodel.reactions:
            if rxn.id in self.fluxes_df.columns:
                if fluxes_abs: flux = abs(rxn.flux)
                else: flux = rxn.flux
                self.fluxes_df.iloc[-1, self.fluxes_df.columns.get_loc(rxn.id)] = flux

    def add_fluxes_from_fluxdict(self, flux_dict:dict, bin_id: Union[float, int, str],
                                 substrate_uptake_rate: Union[float, int], fluxes_abs:bool = True):
        """Adds fluxes directly from a dictionary of reaction fluxes."""
        new_row = [bin_id, substrate_uptake_rate]+[0]*len(self.fluxes_df.columns[2:])
        self.fluxes_df.loc[len(self.fluxes_df)] = new_row
        for rxn_id, flux in flux_dict.items():
            if rxn_id in self.fluxes_df.columns:
                if fluxes_abs: flux = abs(flux)
                else: flux = flux
                self.fluxes_df.iloc[-1, self.fluxes_df.columns.get_loc(rxn_id)] = flux

    def add_error_to_error_df(self, bin_id:Union[int, str], error: float):
        """Records error values for a given bin."""
        self.error_df.loc[len(self.error_df)] = [bin_id] + error

    def remove_simulations_from_flux_df(self, bin_id)-> None:
        """Removes simulations corresponding to a bin."""
        self.fluxes_df = self.fluxes_df[self.fluxes_df['bin'] != bin_id]