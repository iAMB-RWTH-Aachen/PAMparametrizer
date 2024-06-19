import pandas as pd
import numpy as np

def adaptive_sampling(exp_data: pd.DataFrame, num_samples=10):
    """
    Perform adaptive sampling based on data point distribution and variability.

    Parameters:
    - x_exp (array-like): X-coordinates of experimental data.
    - y_exp (array-like): Y-coordinates of experimental data.
    - num_samples (int): Number of samples to generate.
    - min_density (int): Minimum density of samples per unit distance.
    - max_density (int): Maximum density of samples per unit distance.

    Returns:
    - sampled_x (array): Sampled x-coordinates.
    """
    sampled_indices = []
    exp_data = exp_data.reset_index(drop=True)
    remaining_df = exp_data.copy()

    for _ in range(num_samples-2):
        # Compute NaN counts for each column
        nan_counts = remaining_df.isnull().sum()

        # Determine the column with the maximum NaN count
        max_nan_column = nan_counts.idxmax()  # Compute sampling probabilities based on NaN counts in the max_nan_column
        probabilities = np.where(remaining_df[max_nan_column].isnull(), 0, 1).astype(float)
        # Normalize probabilities to sum up to 1
        if np.nansum(probabilities) != 0:
            probabilities /= np.nansum(probabilities)
        else:
            probabilities = [1/len(probabilities) for i in probabilities]

        # Sample index with replacement based on probabilities
        sampled_index = np.random.choice(remaining_df.index, p=probabilities)

        # Remove sampled row from remaining DataFrame
        remaining_df = remaining_df.drop(index=sampled_index)

        # Add sampled index to list
        sampled_indices.append(sampled_index)

        # Break loop if remaining DataFrame is empty
        if remaining_df.empty:
            break

    # Create sampled DataFrame
    #always add first and last datapoint to dataset
    if 0 not in sampled_indices: sampled_indices += [0]
    if len(exp_data)-1 not in sampled_indices: sampled_indices += [len(exp_data)-1]
    sampled_df = exp_data.loc[sampled_indices]
    return sampled_df