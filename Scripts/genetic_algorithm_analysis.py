import pandas as pd
from genetic_algorithm_parametrization.core_parametrization_gaussian import GAPO

"""
function library related to running the genetic algorithm
"""

def get_kcat_error_from_ga(genetic_algorithm: GAPO):
    result_path = str(genetic_algorithm.folderpath_save.joinpath(genetic_algorithm.filename_save)) + '.xlsx'
    results = pd.read_excel(result_path, sheet_name='final_population').sort_values(by=['fitness_weighted_sum'],
                                                                                    ascending=True)

    best_indiv_kcats = [float(kcat) for kcat in results.iloc[0][['attributes']].iloc[0].split(',')]
    lowest_error = results.iloc[0][['fitness_weighted_sum']].iloc[0]

    return {'kcats': best_indiv_kcats, 'error': lowest_error}


def reverse_nested_dictionary(dict_to_reverse:dict, new_key:str='parameter')->dict:
    """
    must make sure all the sub dictionaries have the same layout
    :param dict_to_reverse: nested dictionary to reverse in the form of:
           `{key:{subkey:value, subkey2:value2, subkey3:value3}}`
    :param new_key: str: new key to save the old keys of the upperlevel dictionaries
    :return: reversed_dict: unnested dictionary where the keys of the subdictionary are the keys of upper dictionary
            `{new_key: list_of_keys, subkey:list_of_value, subkey2:list_of_value2, subkey3:list_of_value3}`
    """
    reversed_dict = {new_key: []}
    for key, sub_dict in dict_to_reverse.items():
        reversed_dict[new_key] += [key]
        for subkey,value in sub_dict.items():
            if subkey not in reversed_dict.keys():
                reversed_dict[subkey] = [value]
            else:
                reversed_dict[subkey] += [value]
    return reversed_dict

