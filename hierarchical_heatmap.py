import numpy as np
import pandas as pd
import phenograph  #https://github.com/dpeerlab/PhenoGraph
# import scipy.sparse as ss
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
import seaborn as sns


def get_location_matrix(sample_cube):
    # partition the cube to have a dataframe with x, y coordinates
    location_matrix = {}
    location_matrix["x"] = sample_cube["Cell_X_Position"]
    location_matrix["y"] = sample_cube["Cell_Y_Position"]
    location_matrix = pd.DataFrame.from_dict(location_matrix)
    return location_matrix


def to_99_percentile(expression):
    #https://rdrr.io/github/BodenmillerGroup/bbRtools/man/censor_dat.html
    data = np.copy(expression)
    for i in range(data.shape[1]):
        data[:, i] = data[:, i] / np.percentile(data[:, i], 99)
        data[:, i] = np.where(data[:, i] > 1, 1, data[:, i])
    return data


def to_cluster_mean(data, labels):
    cluster_mean = np.zeros((np.max(labels)+1, data.shape[1]))
    for i in range(np.max(labels)+1):
        # list of indexes with certain label
        indexes = np.where(labels == i)[0]
        array_with_same_label = []
        for j in indexes:
            array_with_same_label.append(data[j])
        array_with_same_label = np.array(array_with_same_label)
        # find mean of each cluster
        cluster_mean[i] = array_with_same_label.mean(axis=0)
    return cluster_mean


# sample_names = file.Sample_Name.unique()
# for sample in sample_names:
#     sample_cube = file.loc[file['Sample_Name'] == sample]
#     data = sample_cube.iloc[:, 5:12].to_numpy()
#     labels = np.load("communities_{}.npy".format(sample))
#     data = to_99_percentile(data)
#     cluster_mean = to_cluster_mean(data, labels)
#     z_score_cluster_mean = zscore(cluster_mean)
#     # plot to heatmap
#     col_names = ["Membrane_CD8", "Membrane_LAG3", "Membrane_CD4", "Membrane_CD20", "Membrane_PD1", "Membrane_FoxP3", "Membrane_CD10"]
#     z_score_cluster_mean_df = pd.DataFrame(data=z_score_cluster_mean, columns=col_names)
#     heatmap = sns.clustermap(z_score_cluster_mean_df)
#     heatmap.savefig("heatmap_{}.png".format(sample))

cube = pd.read_csv(snakemake.input[0])
data = cube.to_numpy()
labels = np.load(snakemake.input[1])
print("processing heatmap...")
data = to_99_percentile(data)
cluster_mean = to_cluster_mean(data, labels)
z_score_cluster_mean = zscore(cluster_mean)
# plot to heatmap
z_score_cluster_mean_df = pd.DataFrame(data=z_score_cluster_mean, columns=cube.columns)
heatmap = sns.clustermap(z_score_cluster_mean_df)
plt.title("Hierarchical Heatmap")
heatmap.savefig(snakemake.output[0])

