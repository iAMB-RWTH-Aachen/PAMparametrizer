import pandas as pd
import numpy as np
from typing import Union, Iterable
import traceback

from sympy.logic.inference import valid


def calculate_r_squared_for_reaction(reaction_id: str, validation_data: pd.DataFrame,
                                     substrate_uptake_id: str,
                                      fluxes: pd.DataFrame) -> float:
    substr_rxn = substrate_uptake_id + '_ub'
    # Take the absolute value of substrate uptake to avoid issues with reaction directionality
    validation_df = validation_data.copy()
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

