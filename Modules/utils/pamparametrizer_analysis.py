import pandas as pd
import numpy as np
import cobra
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Iterable, Union, Literal, Dict, List, Tuple, Callable
from scipy.stats import entropy
from scipy.cluster.hierarchy import fcluster
from sklearn.decomposition import PCA

import PAModelpy
from PAModelpy import PAModel

from .sector_config_functions import change_sector_parameters_with_config_dict, SectorParameterDict

#######
#SETUP METHODS
#######

def set_up_pam_parametrizer_and_get_substrate_uptake_rates(set_up_parametrizer: Callable,
                                                           parametrizer_kwargs:Dict = {'max_substrate_uptake_rate': -0.1},
                                                           substrate_uptake_id: str = 'EX_glc__D_e') -> Tuple:
    parametrizer = set_up_parametrizer(**parametrizer_kwargs)
    parametrizer._init_results_objects()
    substrate_rates = parametrizer._init_validation_df([parametrizer.min_substrate_uptake_rate,
                                                        parametrizer.max_substrate_uptake_rate])[substrate_uptake_id]
    substrate_rates = list(set(substrate_rates)) #make sure there are no duplicates
    substrate_rates = sorted(substrate_rates)
    return parametrizer, substrate_rates




#######
#SIMULATION TOOLS
#######
def get_results_from_simulations(model: Union[PAModel, 'Model'],
                                 substrate_rates: Union[Iterable[float], Iterable[Iterable[float]]],
                                 substrate_ids: Union[str, list[str]] = ['EX_glc__D_e'],
                                 fluxes_to_save: list[str] = None,
                                 proteins_to_save:list[str] = None,
                                 transl_sector_config: Union[dict[str,SectorParameterDict], bool]=True) -> dict[str, pd.DataFrame]:

    if isinstance(substrate_ids, str):
        substrate_ids = [substrate_ids]
    if not isinstance(substrate_rates[0], Iterable):
        substrate_rates = [substrate_rates]

    solution_information = _set_up_solution_info(fluxes_to_save, proteins_to_save)

    for substrate_list, substrate_id in zip(substrate_rates, substrate_ids):
        transl_sector = transl_sector_config
        if isinstance(transl_sector_config, dict):
            if substrate_id not in transl_sector_config: continue
            transl_sector = transl_sector_config[substrate_id]
        model.optimize()
        if isinstance(model, PAModel): _set_up_pamodel_for_simulations(model, substrate_id, transl_sector)

        for substrate in substrate_list:
            try:
                model.change_reaction_bounds(rxn_id=substrate_id,
                                           lower_bound=substrate, upper_bound=0)
            except:
                model.reactions.get_by_id(substrate_id).lower_bound = substrate
                model.reactions.get_by_id(substrate_id).upper_bound = 0

            model.optimize()
            print(f'Running simulations with {round(substrate, 2)} mmol/g_cdw/h of substrate ({substrate_id}) going into the system')
            sol_pam = model.optimize()

            if model.solver.status == 'optimal' and model.objective.value > 0:
                if fluxes_to_save is not None:
                    solution_information['fluxes'] = save_fluxes(sol_pam,
                                                                 model,
                                                                 fluxes_to_save,
                                                                 substrate,
                                                                 solution_information['fluxes'],
                                                                 substrate_id)
                if proteins_to_save is not None:
                    solution_information['proteins'] = save_proteins(model,
                                                                     proteins_to_save,
                                                                     substrate,
                                                                     solution_information['proteins'],
                                                                     substrate_id)
                #reset model
                try:
                    model.change_reaction_bounds(rxn_id=substrate_id,
                                               lower_bound=0, upper_bound=1e3)
                except:
                    model.reactions.get_by_id(substrate_id).lower_bound = 0
                    model.reactions.get_by_id(substrate_id).upper_bound = 1e3


    return solution_information #TODO seems not to save the fluxes correctly


def get_results_from_simulations_fixed_mu(pamodel: PAModel,
                                          growth_rates: Iterable[float],
                                          substrate_id: str = 'EX_glc__D_e',
                                          fluxes_to_save: list[str] = None,
                                          proteins_to_save:list[str] = None,
                                          transl_sector_config=True,
                                          method_ids: list[str]= None) -> dict[str, pd.DataFrame]:
    _set_up_pamodel_for_simulations(pamodel, substrate_id, transl_sector_config)
    solution_information = _set_up_solution_info(fluxes_to_save, proteins_to_save, method_ids)

    #leave the substrate rate open until maximal physiological accurate uptake rate
    pamodel.change_reaction_bounds(substrate_id, -15, 0)

    for method, mu in zip(method_ids, growth_rates):
        # first get the max growth rate
        pamodel.objective = pamodel.BIOMASS_REACTION
        pamodel.change_reaction_bounds(rxn_id=pamodel.BIOMASS_REACTION,
                                       lower_bound=0, upper_bound=mu)
        print('Running simulations with a growth rate of', mu, '1/h')
        pamodel.optimize()
        growth_rate = pamodel.objective.value

        #do another simulation to minimize substrate uptake rate to remove any futile cycles
        pamodel.objective = substrate_id
        pamodel.change_reaction_bounds(rxn_id=pamodel.BIOMASS_REACTION,
                                       lower_bound=growth_rate, upper_bound=growth_rate)


        sol_pam = pamodel.optimize()
        substrate_rate = pamodel.reactions.get_by_id(substrate_id).flux
        print(f'The simulated growth rate is {growth_rate} 1/h with a substrate uptake rate of {substrate_rate} mmol/g_cdw/h')
        if pamodel.solver.status == 'optimal' and growth_rate > 0:
            if fluxes_to_save is not None:
                solution_information['fluxes'] = save_fluxes(sol_pam, pamodel, fluxes_to_save, substrate_rate,
                                                             solution_information['fluxes'], substrate_id, method)
            if proteins_to_save is not None:
                solution_information['proteins'] = save_proteins(pamodel, proteins_to_save, substrate_rate,
                                                                 solution_information['proteins'], substrate_id,method)

    return solution_information

def _set_up_pamodel_for_simulations(pamodel:PAModel,
                                   substrate_id: str,
                                   transl_sector_config:Union[bool, SectorParameterDict]) -> None:
    if not isinstance(transl_sector_config, dict) and transl_sector_config:
        transl_sector_config = {'slope': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_mu[0],
                                'intercept': pamodel.sectors.get_by_id('TranslationalProteinSector').tps_0[0]}

    if transl_sector_config is not False:
        change_sector_parameters_with_config_dict(pamodel=pamodel,
                                                     sector_config = transl_sector_config,
                                                     substrate_uptake_id = substrate_id,
                                                  sector_id = 'TranslationalProteinSector')

def _set_up_solution_info(fluxes_to_save: list[str],
                          proteins_to_save:list[str],
                          method_ids:list[str] = None
                          ) -> dict[str, pd.DataFrame]:
    solution_information = {}
    if fluxes_to_save is not None:
        solution_information['fluxes'] = pd.DataFrame(columns=[ 'substrate_id','substrate'] + fluxes_to_save)
    if proteins_to_save is not None:
        solution_information['proteins'] = pd.DataFrame(
            columns=['substrate_id','enzyme_id','fraction', 'growth_rate', 'substrate_uptake'])
    if method_ids is not None:
        for key, df in solution_information.items():
            df['method'] = []
            solution_information[key]=df
    return solution_information


def save_fluxes(solution: cobra.Solution,
                pamodel:PAModel,
                fluxes_to_save: list[str],
                substrate_rate: float,
                flux_df: pd.DataFrame,
                substrate_id:str,
                method:str = None) -> pd.DataFrame:
    solution_flux = [substrate_id, substrate_rate] + [solution[rxn] if rxn in pamodel.reactions else np.NaN for rxn in
                                   fluxes_to_save]
    if method is not None: solution_flux.append(method)
    flux_df.loc[len(flux_df)] = solution_flux
    flux_df.sort_values('substrate', ascending=False)
    return flux_df

def save_proteins(pamodel:PAModel,
                proteins_to_save: list[str],
                substrate_rate: float,
                protein_df: pd.DataFrame,
                substrate_id:str,
                method: str = None) -> pd.DataFrame:
    total_conc = 0
    for enzid in proteins_to_save:
        enzyme = pamodel.enzymes.get_by_id(enzid)
        conc = enzyme.concentration
        total_conc += conc
        new_row = [substrate_id,enzid, conc, pamodel.reactions.get_by_id(pamodel.BIOMASS_REACTION).flux, substrate_rate]
        if method is not None: new_row.append(method)
        protein_df.loc[len(protein_df)] = new_row
    return protein_df

def calculate_error_for_reactions(validation_df: pd.DataFrame,
                                   flux_df: pd.DataFrame,
                                  rxns_to_validate: list) -> float:
    # calculate error for different exchange rates
    error = []
    for rxn, substrate_id in zip(rxns_to_validate, flux_df.substrate_id):
        # only select the rows which are filled with data
        validation_data = validation_df.dropna(axis=0, subset=rxn)
        try: validation_data = validation_data.loc[substrate_id]
        except: pass
        # if there are no reference data points, continue to the next reaction
        if len(validation_data) == 0:
            continue
        r_squared = calculate_r_squared_for_reaction(rxn, validation_df, substrate_id,
                                                           flux_df[flux_df.substrate_id == substrate_id])
        error += [r_squared]
    return error

def calculate_r_squared_for_reaction(reaction_id: str, validation_data: pd.DataFrame,
                                     substrate_uptake_id: str,
                                      fluxes: pd.DataFrame) -> float:
    substr_rxn = substrate_uptake_id
    # Take the absolute value of substrate uptake to avoid issues with reaction directionality
    validation_data[substr_rxn] = [round(abs(flux),4) for flux in validation_data[substr_rxn]]
    simulated_data = pd.DataFrame({substr_rxn: [round(abs(flux),4) for flux in fluxes['substrate']],
                                   'simulation': fluxes[reaction_id]})
    ref_data_rxn = pd.merge(validation_data,simulated_data,on=substr_rxn, how='inner')
    # error: squared difference
    ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[reaction_id] - x['simulation']) ** 2)

    # calculate R^2:
    data_average = np.nanmean(validation_data[reaction_id])
    residual_ss = np.nansum(ref_data_rxn.error)
    total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[reaction_id]])
    # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
    if (residual_ss == 0) | (total_ss == 0):
        r_squared = 0
    else:
        r_squared = 1 - residual_ss / total_ss
    return r_squared

def calculate_difference_simulation_experiment(validation_df:pd.DataFrame,
                                               flux_df:pd.DataFrame,
                                               rxns_to_validate:list[str],
                                               substr_rxn:str):
    differences = []
    for rxn in rxns_to_validate:
        # only select the rows which are filled with data
        validation_data = validation_df.dropna(axis=0, subset=rxn)
        # if there are no reference data points, continue to the next reaction
        if len(validation_data) == 0:
            continue
        validation_data[substr_rxn] = [round(abs(flux),4) for flux in validation_data[substr_rxn]]
        simulated_data = pd.DataFrame({substr_rxn: [round(abs(flux_df['substrate']),4)],
                                   'simulation': flux_df[rxn]})
        ref_data_rxn = pd.merge(validation_data,simulated_data,on=substr_rxn, how='inner')
        # error: squared difference
        differences += [row[rxn] - row['simulation'] if not np.isnan(row.simulation) else row[rxn] for i,row in ref_data_rxn.iterrows()]
    return differences


########
#ANALYSIS OF KCATS AND ENZYMES
########


def calculate_kcat_differences(df_grouped: pd.DataFrame,
                               rxn2kcat: dict[str, dict]) -> pd.DataFrame:
    """
    Calculate the changes in the kcat values for all simulation. Designed to be used in an apply method on a dataframe
    grouped by the different alternative models of a single condition.

    Args:
        df_grouped: group of a grouped dataframe (by model), containing at leat the columns: 'run_id', 'enzyme_id', and 'kcat[s-1]'

    Return:
        pd.DataFrame with changed kcat values (additional columns: 'previous kcat', 'kcat_change', 'absolute_kcat_change',
        'relative_change', and 'absolute_relative_change'
    """
    group_sorted = df_grouped.sort_values('run_id')
    group_sorted = group_sorted.groupby(
        'enzyme_id', group_keys=False).apply(calculate_relative_change, rxn2kcat).reset_index(drop=True)
    return group_sorted


def calculate_relative_change(group:pd.DataFrame,
                              rxn2kcat: dict[str,dict]) -> pd.DataFrame:
    previous_kcat = get_previous_kcat_values(group, rxn2kcat)

    group['previous_kcat'] = previous_kcat
    group['kcat_change'] = group['kcat[s-1]'] - previous_kcat  # .abs()
    group['absolute_kcat_change'] = group['kcat_change'].abs()

    group['relative_change'] = group['kcat_change'] / previous_kcat  # .abs()
    group['absolute_relative_change'] = group['relative_change'].abs()
    group.replace([np.inf, -np.inf], 0, inplace=True)

    return group


def get_previous_kcat_values(group: pd.DataFrame,
                             rxn2kcat:dict[str, dict]) -> float:
    kcat_default = 33
    previous_kcat = group['kcat[s-1]'].shift(1)
    # if there is no previous kcat change, get it from original kcat dataset
    for i in group.index:
        if pd.isna(previous_kcat.loc[i]):
            entry = group.loc[i]
            if entry['enzyme_id'] not in rxn2kcat[entry['rxn_id']]:
                original_kcat_value = kcat_default
            else:
                original_kcat_value = rxn2kcat[entry['rxn_id']][entry['enzyme_id']][entry['direction']]
            previous_kcat.loc[i] = original_kcat_value
    return previous_kcat


def convert_peptide_to_enzyme_concentrations(peptide_df: pd.DataFrame,
                                             enzyme_db:pd.DataFrame,
                                             concentration_columns: List[str]) -> pd.DataFrame:
    """
    Computes the enzyme concentration for each reaction-enzyme relation, considering enzyme complexes and isozymes.

    The function classifies each enzyme-reaction pair into:
    - **Homomeric Enzyme**: A single enzyme with no isozymes.
    - **Enzyme Complex**: A functional unit composed of multiple enzyme subunits.
    - **Isozyme**: Different enzymes catalyzing the same reaction.

    Concentration Calculation Rules:
    - Homomeric enzymes retain their original concentration.
    - Enzyme complexes take the minimum concentration of all participating peptides.
    - Isozymes take the sum of the concentration of all enzymes contributing to the reaction.

    Args:
        peptide_df (pd.DataFrame): DataFrame containing enzyme-specific peptide concentrations.
            Must include:
            - `enzyme_id` (str): The unique identifier for the enzyme.
            - At least one of the columns in `concentration_columns` containing enzyme concentrations.

        enzyme_db (pd.DataFrame): DataFrame containing enzyme-reaction associations.
            Must include:
            - `rxn_id` (str): The reaction identifier.
            - `enzyme_id` (str): The enzyme(s) catalyzing the reaction (single ID or complex in "_"-joined format).

        concentration_columns (list of str): List of column names in `peptide_df` for which concentrations should be computed.

    Returns:
        pd.DataFrame: A DataFrame containing:
            - `rxn_id` (str): Reaction identifier.
            - `enzyme_id` (str): Enzyme ID(s) associated with the reaction.
            - Computed concentration columns based on `concentration_columns`.
            - `enzyme_type` (str): Classification as "Homomer", "Complex", or "Isozyme".
    """

    # Create a dictionary for quick lookups of concentrations for each enzyme_id
    peptide_concentration = {col: peptide_df.set_index("enzyme_id")[col].to_dict() for col in concentration_columns}

    # Store results
    results = []

    # Process each reaction as a group
    for rxn_id, group in enzyme_db.groupby("rxn_id"):
        enzyme_entries = group["enzyme_id"].tolist()
        enzyme_sets = [eid.split("_") for eid in enzyme_entries]  # Handle enzyme complexes

        # Track classification and concentration storage
        enzyme_type = "Homomer"  # Default assumption
        # Store concentration data per enzyme
        concentration_data = {col: {} for col in concentration_columns}

        if len(enzyme_entries) > 1:  # More than one enzyme for the reaction → Potential Isozyme
            enzyme_type = "Isozyme"

            # Sum concentrations for all isozymes
            for col in concentration_columns:
                total_conc = sum(
                    min(peptide_concentration[col].get(e, 0) for e in enzymes) for enzymes in enzyme_sets
                )
                # Assign the same summed concentration to all participating enzymes
                for enzyme in enzyme_entries:
                    concentration_data[col][enzyme] = total_conc


        else:  # Single enzyme (homomer or complex)
            for enzyme, enzyme_ids in zip(enzyme_entries, enzyme_sets):
                enzyme_type = "Homomer" if len(enzyme_ids) == 1 else "Complex"

                for col in concentration_columns:
                    concentration_data[col][enzyme] = (
                        min(peptide_concentration[col].get(e, 0) for e in enzyme_ids)
                        if len(enzyme_ids) > 1 else peptide_concentration[col].get(enzyme, 0)
                    )
        # Store result for this reaction
        # Append row for each enzyme
        for enzyme in enzyme_entries:
            results.append({
                "rxn_id": rxn_id,
                "enzyme_id": enzyme,
                **{col: concentration_data[col][enzyme] for col in concentration_columns},
                "enzyme_type": enzyme_type
            })

    return pd.DataFrame(results)


def parse_enzyme_complex_id(enzyme_id: str):
    """Convert an enzyme complex id to a list of peptide identifiers while ignoring the default enzyme ids"""
    if 'Enzyme_' in enzyme_id:
        sorted_enzyme_id = [enzyme_id.replace('Enzyme_', '')]
    else:
        # make sure the enzyme complex is in the right order, otherwise match will fail
        sorted_enzyme_id = sorted(enzyme_id.split('_'))
    return sorted_enzyme_id


def normalize_simulated_protein_concentrations(df: pd.DataFrame,
                                enzyme_db: pd.DataFrame,
                                ue_sector: PAModelpy.UnusedEnzymeSector):
    """Normalizes protein fractions in a simulation dataset.

    This function processes simulation results by converting predicted enzyme concentrations form mmol/gCDW to
    g/gCDW and by normalizing protein concentrations.
    It accounts for unused enzyme fractions and converts enzyme-level simulated concentration into
    enzyme-level fractions. The final concentration values are normalized based
    on the total protein concentration in the system, taking unused proteins into account.

    Args:
        df (pd.DataFrame):
            The simulation dataframe containing columns:
            - 'method' (str): The simulation method identifier.
            - 'fraction' (float): Protein fraction values in mmol/g_cdw.
            - 'growth_rate' (float): Growth rate for each condition.
            - 'molMass' (int): molar mass of the enzyme in kDa.
        enzyme_db (pd.DataFrame):
            The enzyme database mapping peptides to their corresponding enzymes.
        ue_sector (UnusedEnzymeSector):
            An object representing the unused enzyme sector, containing attributes:
            - `ups_mu` (float): Scaling factor based on growth rate.
            - `ups_0` (list[float]): Baseline unused enzyme concentration.

    Returns:
        pd.DataFrame: A dataframe with normalized enzyme concentrations, containing:
        - 'enzyme_id' (str): The enzyme identifier.
        - 'normalized_fraction' (float): The normalized protein fraction.
        - 'method' (str): The simulation method.
        - 'enzyme_type' (str): if the enzyme is an isozyme, homomer or complex
    """

    df = df.dropna(how='any')

    #convert mmol to grams of protein
    #mmol*1e-3 -> mol
    #mol*g/mol -> mol
    df['fraction'] = df.fraction *1e-3 * df.molMass

    peptide_simulated_df_list = []  # Store results in a list

    for method, simulation_df in df.groupby('method'):
        # unused enzymes should be added to the total sum to be comparable to 'in vivo' data
        unused_enzyme_sum = simulation_df.growth_rate.iloc[0] * ue_sector.ups_mu + ue_sector.ups_0[0]
        total_conc = simulation_df.fraction.sum() * 1e6 + unused_enzyme_sum

        # sum the concentrations of isozymes
        simulation_df = convert_peptide_to_enzyme_concentrations(simulation_df, enzyme_db, ['fraction'])

        simulation_df['normalized_fraction'] = simulation_df['fraction'] * 1e6 / total_conc
        simulation_df['method'] = method

        peptide_simulated_df_list.append(simulation_df[['enzyme_id', 'rxn_id','normalized_fraction', 'method', 'enzyme_type']])

    return pd.concat(peptide_simulated_df_list, ignore_index=True)

#########
#STATISTICS AND CLUSTERING
########


def get_clusters_from_clustermap(clustermap: sns.matrix.ClusterGrid,
                                 df: pd.DataFrame,
                                 nrow_clusters: int = 10,
                                 ncol_clusters: int = 2):
    df_clustered = df.copy()

    # Extract row and column linkage
    row_linkage = clustermap.dendrogram_row.linkage
    col_linkage = clustermap.dendrogram_col.linkage

    # Determine cluster membership using scipy's fcluster function
    row_clusters = fcluster(row_linkage, t=nrow_clusters, criterion='maxclust')
    col_clusters = fcluster(col_linkage, t=ncol_clusters, criterion='maxclust')

    # Add cluster labels to your data
    df_clustered['Row_Cluster'] = row_clusters

    #     # Example output
    #     print("Row Clusters:", max(row_clusters))
    #     print(df_clustered[['Row_Cluster']])
    #     for cluster in range(1,nrow_clusters+1):
    #         clustered_kcats = df_clustered[df_clustered['Row_Cluster']==cluster]
    #         print(len(clustered_kcats.Row_Cluster))
    #         print('median', clustered_kcats.median())
    #         print('mean', clustered_kcats.mean())
    #         print('stdev', clustered_kcats.std())

    #     print("\nColumn Clusters:")
    col_clusters_dict = {clst: [] for clst in range(1, ncol_clusters + 1)}
    for col, cluster in zip(df_clustered.columns, col_clusters):
        col_clusters_dict[cluster].append(col)
    #         print(f"{col}: Cluster {cluster}")
    return df_clustered, col_clusters_dict


def select_clustered_rows_by_variation(
    clustered_data: pd.DataFrame,
    column_clusters: dict[str,list],
    per_cluster: bool = True,
    select_highest: bool = True,
    num_rows: int = 10,
    metric: Literal['MAD', 'ENT', 'STD', 'CV'] = 'MAD'
) -> pd.DataFrame:
    """
    Select rows from clustered data based on variation metrics across column clusters.

    This function calculates variation metrics (e.g., MAD, entropy, CV, or standard deviation)
    for rows across defined column clusters and selects rows with the highest or lowest variation
    globally or within each row cluster.

    Args:
    -----------
    clustered_data : pd.DataFrame
        DataFrame containing the data with row clusters labeled under a 'Row_Cluster' column.
        Other columns correspond to the features for which variation is computed.
    column_clusters : dict
        A dictionary where keys are column cluster identifiers and values are lists of column names
        belonging to each column cluster.
    per_cluster : bool, default=True
        If `True`, selects rows with the highest or lowest variation within each row cluster.
        If `False`, selects rows globally based on variation.
    select_highest : bool, default=True
        If `True`, selects rows with the highest variation. If `False`, selects rows with the lowest variation.
    num_rows : int, default=10
        Number of rows to select if `per_cluster` is `False`.
    metric : {'MAD', 'ENT', 'STD', 'CV'}, default='MAD'
        The metric to calculate variation:
        - 'MAD': Median Absolute Deviation.
        - 'ENT': Shannon Entropy.
        - 'STD': Standard Deviation.
        - 'CV': Coefficient of Variation.

    Returns:
    --------
    pd.DataFrame
        A DataFrame containing the selected rows with additional columns for variation metrics.

    Notes:
    ------
    - Rows are normalized using z-score normalization before calculating variation metrics.
    - If `per_cluster` is `True`, the function identifies the row with the highest/lowest variation
      within each row cluster.

    """

    # Normalize rows using z-score normalization
    normalized_data = row_wise_zscore_normalization(
        clustered_data.drop(columns=['Row_Cluster'])
    ).merge(
        clustered_data[['Row_Cluster']],
        left_index=True,
        right_index=True
    )

    variation_df = _compute_variation_metrics_per_row(normalized_data, column_clusters)

    # Combine variation metrics with the normalized data
    data_with_variation = normalized_data.merge(
        variation_df, left_index=True, right_on='Row_Index'
    )

    # Select rows based on the specified conditions
    if per_cluster:
        if select_highest:
            selected_rows = data_with_variation.loc[
                data_with_variation.groupby('Row_Cluster')[metric].idxmax()
            ]
        else:
            selected_rows = data_with_variation.loc[
                data_with_variation.groupby('Row_Cluster')[metric].idxmin()
            ]
    else:
        selected_rows = data_with_variation.sort_values(
            metric, ascending=not select_highest
        ).head(num_rows)

    return selected_rows

def row_wise_zscore_normalization(df: pd.DataFrame):
    return df.apply(lambda row: (row - row.mean()) / row.std(), axis=1)

def _compute_variation_metrics_per_row(normalized_data: pd.DataFrame,
                                       column_clusters: dict[str, list]) -> pd.DataFrame:
    # Compute variation metrics
    variation_metrics = []
    for row_index, row in normalized_data.iterrows():
        cluster_means = [
            row[cols].mean() for cols in column_clusters.values()
        ]
        variation_metrics.append({
            'Row_Index': row_index,
            'MAD': np.median(np.abs(cluster_means - np.median(cluster_means))),
            'ENT': entropy(cluster_means),
            'CV': np.std(cluster_means) / np.mean(cluster_means) if np.mean(cluster_means) != 0 else np.nan,
            'STD': np.std(cluster_means)
        })

    return pd.DataFrame(variation_metrics)

#########
#PLOTTING
#########
def plot_histogram_logspace(ax: plt.Axes, data:pd.DataFrame,
                            color:Union[str,float],
                            label: str, x_label:str, n_bins:int=50,
                            relative:bool=False, ymax: float = None) -> None:
    if not relative:
        kwargs = {}
        ax.set_yscale('log')
        if not ymax: ymax = 1e3
    else:
        kwargs = {'weights': np.ones_like(data) / len(data)}
        if not ymax: ymax = 0.35

    hist, bins = np.histogram(data, bins=n_bins, **kwargs)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))

    bin_heights, bin_borders, _ = ax.hist(data, bins=logbins, histtype='step',
                                          stacked=True, fill=False, label=label, color=color, **kwargs)

    ax.vlines(13.7, 0, ymax, linestyle='dotted')
    ax.set_ylim([0, ymax])

    ax.set_xlabel(x_label)
    ax.set_xscale('log')

    ax.set_xlim([1e-3, 1e9])

def plot_PCA_graph(input_df, columns_to_analyse: list[str], hue:str, values:str):
    # Pivot the data to have alternatives as rows and enzymes as columns
    # This creates a wide matrix required for PCA
    data_matrix = input_df.pivot_table(
        index=columns_to_analyse,
        columns=hue,
        values=values,
        aggfunc='mean',
        fill_value=0
    )

    data_matrix.reset_index(names=columns_to_analyse, inplace=True)
    mapping_to_analyse = dict()
    for col in columns_to_analyse:
        mapping_to_analyse[col] = data_matrix[col]

    data_for_pca = data_matrix.drop(columns=columns_to_analyse)

    # Apply PCA
    pca = PCA(n_components=2)
    reduced_data = pca.fit_transform(data_for_pca)
    explained_variance = pca.explained_variance_ratio_
    print(f"PCA explained variance: {explained_variance}")

    # Create a DataFrame for plotting
    reduced_df = pd.DataFrame(reduced_data, columns=['PC1', 'PC2'])

    for col, col_series in mapping_to_analyse.items():
        reduced_df[col] = col_series

    # Merge back descriptions for coloring
    mapping = input_df[columns_to_analyse+[hue]].drop_duplicates()
    reduced_df = pd.merge(reduced_df, mapping, on=columns_to_analyse, how='left')

    # Plot the reduced data
    plt.figure(figsize=(12, 8))
    sns.scatterplot(
        x='PC1', y='PC2', hue='COG description', data=reduced_df, palette='tab20', s=100, alpha=0.8
    )

    # Add titles and labels
    plt.title('PCA of Enzyme Changes per Alternative', fontsize=16)
    plt.xlabel(f'Principal Component 1 ({round(explained_variance[0], 3) * 100}%)', fontsize=12)
    plt.ylabel(f'Principal Component 2 ({round(explained_variance[1], 3) * 100}%)', fontsize=12)

    # Position the legend outside the plot
    plt.legend(
        bbox_to_anchor=(1.05, 0.5),  # Position the legend to the right-hand side, centered vertically
        loc='center left',  # Align it to the left of the bounding box
        title=hue
    )
    plt.tight_layout()
    plt.show()