import itertools
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import zscore
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt

"""
How to use:
expression_matrix, location_matrix, axis_matrix, size = read_data(pd.read_csv(network_0.csv), 5, "IMC")
--
expression_matrix: all the marker expressions (num_channels cols)
location_matrix: all x, y coordinates (2 cols)
axis_matrix: all major and minor axes (2 cols)
size: number of cells in graph (int)
"""


def read_data(matrix, num_channels, data_type):
    if data_type == "IMC":
        return get_imc_matrix(matrix, num_channels)
    elif data_type == "IHC":
        return get_ihc_matrix(matrix, num_channels)
    else:
        return Exception("Did not specify data type from {IMC, IHC}")


def get_imc_matrix(matrix, num_channels):
    """
    When we want to generate a single value for the cell, we normally sum the flux and f_buffer column.
    :param given pd read matrix
    :return: expression and location matrix in Dataframe
    """
    read_matrix = matrix
    expression_matrix, location_matrix, axis_matrix = {}, {}, {}
    # reading in expressions
    for channel in np.arange(num_channels):
        if np.sum(read_matrix["flux_" + str(channel + 1).zfill(2)]) == 0:
            break
        expression_matrix["channel_" + str(channel + 1).zfill(2)] = \
            read_matrix["flux_" + str(channel + 1).zfill(2)] \
            + read_matrix["f_buffer_" + str(channel + 1).zfill(2)]
    expression_matrix = pd.DataFrame.from_dict(expression_matrix)

    # reading locations
    location_matrix["x"] = read_matrix["X_image"]
    location_matrix["y"] = read_matrix["Y_image"]
    location_matrix = pd.DataFrame.from_dict(location_matrix)

    # reading axis
    axis_matrix["Major"] = read_matrix["ell_smaj"]
    axis_matrix["Minor"] = read_matrix["ell_smin"]
    axis_matrix = pd.DataFrame.from_dict(axis_matrix)
    axis_matrix = axis_matrix.to_numpy()

    expression_matrix = expression_matrix.to_numpy()
    # TODO: standardize or not?
    # for index in np.arange(expression_matrix.shape[0]):
    #     # logged_value = np.log(expression_matrix[index, :] + 1)
    #     expression_matrix[index, :] = zscore(expression_matrix[index, :])
    location_matrix = location_matrix.to_numpy()

    size = expression_matrix.shape[0]
    return expression_matrix, location_matrix, axis_matrix, size


def get_ihc_matrix(matrix, num_channels):
    """
    :param given pd read matrix
    :return: expression and location matrix in Dataframe
    """
    read_matrix = matrix
    expression_matrix, location_matrix, axis_matrix = {}, {}, {}
    # reading in expressions
    expression_matrix = read_matrix.iloc[:, 5:5 + num_channels]

    # reading locations
    location_matrix["x"] = read_matrix["Cell_X_Position"]
    location_matrix["y"] = read_matrix["Cell_Y_Position"]
    location_matrix = pd.DataFrame.from_dict(location_matrix)

    # reading axis
    axis_matrix["Major"] = read_matrix["Membrane_Major_Axis"]
    axis_matrix["Minor"] = read_matrix["Membrane_Minor_Axis"]
    axis_matrix = pd.DataFrame.from_dict(axis_matrix)
    axis_matrix = axis_matrix.to_numpy()

    expression_matrix = expression_matrix.to_numpy()
    for index in np.arange(expression_matrix.shape[0]):
        # logged and z-scored
        # logged_value = np.log(expression_matrix[index, :] + 1)
        expression_matrix[index, :] = zscore(expression_matrix[index, :])
    location_matrix = location_matrix.to_numpy()

    size = expression_matrix.shape[0]
    return expression_matrix, location_matrix, axis_matrix, size


print("extracting data...")
FILENAME = snakemake.input[0]
matrix = pd.read_csv(FILENAME)
expression_matrix, location_matrix, axis_matrix, size = read_data(matrix, 36, "IMC")
np.save(snakemake.output[0], expression_matrix)
np.save(snakemake.output[1], location_matrix)
np.save(snakemake.output[2], axis_matrix)
