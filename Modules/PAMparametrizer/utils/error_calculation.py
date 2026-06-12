import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from typing import Union, Iterable, Optional, List, Tuple
from PAModelpy import PAModel


def calculate_r_squared_for_reaction(reaction_id: str, validation_data: pd.DataFrame,
                                     substrate_uptake_id: str,
                                      fluxes: pd.DataFrame) -> float:
    validation_df = validation_data.copy()
    substr_rxn = substrate_uptake_id + '_ub' if substrate_uptake_id+ '_ub' in validation_df.columns else substrate_uptake_id

    # Take the absolute value of substrate uptake to avoid issues with reaction directionality
    validation_df[substr_rxn] = [round(abs(flux),4) for flux in validation_df[substr_rxn]]
    simulated_data = pd.DataFrame({substr_rxn: [round(abs(flux),4) for flux in fluxes['substrate']],
                                   'simulation': fluxes[reaction_id]})
    ref_data_rxn = pd.merge(validation_df,simulated_data,on=substr_rxn, how='inner')
    # error: squared difference
    ref_data_rxn = ref_data_rxn.assign(error=lambda x: (x[reaction_id] - x['simulation']) ** 2)

    # calculate R^2:
    data_average = np.nanmean(validation_df[reaction_id])
    residual_ss = np.nansum(ref_data_rxn.error)

    if len(ref_data_rxn[reaction_id])==1:
        return custom_error(ref_data_rxn['simulation'].iloc[0],ref_data_rxn[reaction_id].iloc[0])
    else:
        total_ss = np.nansum([(data - data_average) ** 2 for data in ref_data_rxn[reaction_id]])
    # calculating r_squared is only feasible of the numerator and the denomenator are both nonzero
    if residual_ss == 0:
        r_squared = 1
    elif total_ss == 0:
        r_squared = 0
    else:
        r_squared = 1 - residual_ss / total_ss
    return r_squared


def calculate_smape_for_reaction(reaction_id: str, validation_data: pd.DataFrame,
                                     substrate_uptake_id: str,
                                      fluxes: pd.DataFrame) -> float:
    substr_rxn = substrate_uptake_id + '_ub'
    # Take the absolute value of substrate uptake to avoid issues with reaction directionality
    validation_df = validation_data.copy()
    validation_df[substr_rxn] = [round(abs(flux), 4) for flux in validation_df[substr_rxn]]
    simulated_data = pd.DataFrame({substr_rxn: [round(abs(flux), 4) for flux in fluxes['substrate']],
                                   'simulation': fluxes[reaction_id]})
    ref_data_rxn = pd.merge(validation_df, simulated_data, on=substr_rxn, how='inner')
    return calculate_symmetric_mean_absolute_percentage_error(ref_data_rxn[reaction_id], ref_data_rxn['simulation'])



def custom_error(observed, simulated, lambda_factor=1.0):
    """
    Calculate custom error where error is 1 if distance between observed and simulated is 0,
    and approaches 0 as distance increases using an exponential decay function.

    Args:
    observed (float): The observed datapoint value.
    simulated (float): The simulated value.
    lambda_factor (float): The scaling factor for the exponential decay. Higher values lead to faster decay.

    Returns:
    float: The calculated error.
    """
    distance = abs(observed - simulated)
    error = np.exp(-lambda_factor * distance)
    return error

def nanaverage(data:Union[list],weights:dict = None,axis:int = None) -> Iterable:
    masked_data = np.ma.masked_array(data, np.isnan(data))
    average = np.ma.average(masked_data, axis=axis, weights=weights)
    return average


def calculate_symmetric_mean_absolute_percentage_error(y_true:Iterable[float],
                                                       y_pred:Iterable[float]) -> float:
    """
    Compute Symmetric Mean Absolute Percentage Error (sMAPE)

    Parameters:
        y_true (array-like): True values
        y_pred (array-like): Predicted values

    Returns:
        float: sMAPE value
    """
    y_true = np.nan_to_num(np.array(y_true))
    y_pred = np.nan_to_num(np.array(y_pred))

    numerator = np.abs(y_true - y_pred)
    denominator = (np.abs(y_true) + np.abs(y_pred)) / 2

    # Avoid division by zero by replacing zero denominators with a small constant
    denominator = np.where(denominator == 0, 1e-10, denominator)

    smape_value = np.mean(numerator / denominator) * 100
    return smape_value


def _align_frames(
    validation_df: pd.DataFrame,
    flux_df: pd.DataFrame,
    rxns_to_validate,
    substr_rxn: str,
    substrate_sim: str = "substrate",
) -> pd.DataFrame:
    """
    Merge experimental and simulated data on the substrate identifier.
    The function also makes sure that all reactions we are interested in exist
    in both tables (otherwise they are dropped automatically).

    Returns a DataFrame with columns:
        - substrate (the joining key)
        - <exp_rxn_1>, <exp_rxn_2>, …            (experimental values)
        - <sim_rxn_1>, <sim_rxn_2>, …            (simulation values)
    """

    for df,rxn in zip([flux_df, validation_df],[substrate_sim, substr_rxn]):
        df[rxn] = df[rxn].abs()
    exp_sub = validation_df[[substr_rxn] + rxns_to_validate].copy().sort_values(substr_rxn)
    sim_sub = flux_df[[substrate_sim] + rxns_to_validate].copy().sort_values(substrate_sim)
    sim_sub = sim_sub.rename(columns={substrate_sim: substr_rxn})
    merged = pd.merge_asof(exp_sub, sim_sub,
                      on=substr_rxn,
                      suffixes=('_exp', '_sim'),
                           direction = 'nearest',
                           tolerance = 1e-5)
    return merged

def calulate_pearson_correlation_simulation_vs_experiment(
    validation_df: pd.DataFrame,
    flux_df: pd.DataFrame,
    rxns_to_validate: List[str],
    substr_rxn: str,
    substrate_sim: str = "substrate",
    absolute_rxns: Optional[List[str]] = ['PGI', 'EX_ac_e', 'EX_glc__D_e', 'substrate']
) -> Tuple[float, float]:
    """
    1️⃣  Align experimental and simulated data on the substrate column.
    2️⃣  (Optional) take absolute values for the reactions listed in `absolute_rxns`.
    3️⃣  Stack all reactions into two long 1‑D arrays.
    4️⃣  Compute a *single* Pearson r and its two‑sided p‑value.
    """
    if absolute_rxns:
        for rxn in absolute_rxns:
            if rxn in flux_df.columns:
                flux_df[rxn] = flux_df[rxn].abs()
    merged = _align_frames(validation_df, flux_df, rxns_to_validate, substr_rxn, substrate_sim)
    merged_clean = merged.dropna(axis =1)
    rxns_to_validate = [rxn for rxn in rxns_to_validate if (f"{rxn}_exp" in merged_clean.columns)and(f"{rxn}_sim" in merged_clean.columns)]


    exp_cols = [f"{rxn}_exp" for rxn in rxns_to_validate]
    sim_cols = [f"{rxn}_sim" for rxn in rxns_to_validate]

    merged_clean = merged.dropna(subset = exp_cols+sim_cols ,axis =0)

    exp_long = merged_clean[exp_cols].values.ravel()
    sim_long = merged_clean[sim_cols].values.ravel()

    if np.std(exp_long) == 0 or np.std(sim_long) == 0:
        raise ValueError("One of the aggregated vectors is constant → Pearson undefined.")
    r, p = pearsonr(exp_long, sim_long)
    return r, p
