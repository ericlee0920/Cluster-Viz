import sys
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import itertools
from numba import jit
from numba.typed import List
import seaborn as sns


def get_location_matrix(sample_cube):
    # partition the cube to have a dataframe with x, y coordinates
    location_matrix = {}
    # TODO: mind the names
    location_matrix["x"] = sample_cube["Cell_X_Position"]
    location_matrix["y"] = sample_cube["Cell_Y_Position"]
    location_matrix = pd.DataFrame.from_dict(location_matrix)
    return location_matrix


def get_labels(sample_cube, label_name):
    # label according to label_name
    label_with = sample_cube[label_name]
    label_names = label_with.unique()
    label_dict = {label_names[i]: i for i in range(len(label_names))}
    label_done = [label_dict[label_with.to_numpy()[i]] for i in range(len(label_with))]
    return label_done


def get_neighbor_graph(location_matrix, labels, threshold=0.85):
    # calculate pairwise distance to form similarity matrix
    pairwise_distance = squareform(pdist(location_matrix))
    pairwise_distance /= np.max(pairwise_distance)
    distance_upperbound = np.max(pairwise_distance)
    similarity_matrix = distance_upperbound - pairwise_distance

    # set user threshold, default at 0.98, make adjacency matrix
    threshold = threshold
    adjacency_matrix = np.where(similarity_matrix < threshold, 0, similarity_matrix)
    adjacency_matrix = np.where(similarity_matrix == distance_upperbound, 0, adjacency_matrix)
    neighbour_graph = nx.from_numpy_matrix(adjacency_matrix)

    # set labels to neighbor graph
    nx.set_node_attributes(neighbour_graph, dict(enumerate(labels)), 'label')
    return neighbour_graph


def get_union_graph(location_matrix, labels, axis_matrix):
    # get upper triangle
    pairwise_distance = squareform(pdist(location_matrix))
    pairwise_distance = np.triu(pairwise_distance, 0)
    major = axis_matrix[:, 0]
    minor = axis_matrix[:, 1]
    pairwise_indexes = list(itertools.combinations(range(len(labels)), 2))
    typed_indexes = List()
    [typed_indexes.append(i) for i in pairwise_indexes]
    # if dist[A, B] < Ar + Br, remain, else set to 0
    pairwise_distance = pairwise_distance_helper(pairwise_distance, typed_indexes, major, minor)
    neighbour_graph = nx.from_numpy_matrix(pairwise_distance)
    # set labels to neighbor graph
    nx.set_node_attributes(neighbour_graph, dict(enumerate(labels)), 'label')
    # nx.write_gpickle(neighbour_graph, "graph.gpickle")
    return neighbour_graph


@jit(nopython=True)
def pairwise_distance_helper(pairwise_dist, pairwise_indexes, major, minor):
    pairwise_distance = np.copy(pairwise_dist)
    for x, y in pairwise_indexes:
        # if pairwise_distance[x, y] > major[x] + major[y]:
        # if pairwise_distance[x, y] > minor[x] + minor[y]:
        if pairwise_distance[x, y] > ((major[x]+minor[x])/2 + (major[y]+minor[y])/2):
            pairwise_distance[x, y] = 0
    return pairwise_distance

def plot_neighbor_graph(neighbor_graph, location_matrix):
    # plot the neighbor graph
    pos = dict([i for i in enumerate(location_matrix)])
    node_color = list(nx.get_node_attributes(neighbor_graph, "label").values())
    # plt.figure(figsize=(30, 28))
    # Change the node_size here
    nodeTypes = list(range(max(node_color)+1))
    nodeTypeDict = {i: [] for i in nodeTypes}
    np.random.seed(100)
    colors = []
    for i in nodeTypes:
        colors.append('#%06X' % np.random.randint(0, 0xFFFFFF))
    nodeColorDict = dict(zip(nodeTypes, colors))

    for i in range(len(node_color)):
        label_here = node_color[i]
        nodeTypeDict[label_here].append(i)
    nodePos = {}
    for i in range(len(pos)):
        nodePos[i] = (pos[i][0], pos[i][1])

    fig, ax = plt.subplots(1, figsize=(16, 16))
    for nt in nodeTypes:
        nlist = nodeTypeDict[nt]
        ncolor = nodeColorDict[nt]
        nx.draw_networkx_nodes(neighbor_graph,
                               pos=nodePos,
                               ax=ax,
                               node_color=ncolor,
                               nodelist=nlist,
                               label=nt,
                               node_size=40
                               )
        nx.draw_networkx_edges(neighbor_graph, pos, width=1.5, alpha=0.5)
    ax.legend(scatterpoints=1)
    # nx.draw_networkx(neighbor_graph, node_color=node_color, pos=pos, with_labels=False, ax=ax)


print("plotting spatial graph...")
print("wait around 5 minutes...")
data = np.load(snakemake.input[0])
location_matrix = np.load(snakemake.input[1])
axis_matrix = np.load(snakemake.input[2])
labels = np.load(snakemake.input[3])
# # neighbor_graph = get_neighbor_graph(location_matrix, labels)
neighbor_graph = get_union_graph(location_matrix, labels, axis_matrix)
plot_neighbor_graph(neighbor_graph, location_matrix)
plt.legend()
plt.title("Neighbor Graph of IMC")
plt.savefig(snakemake.output[0])
plt.show()


pair_number = np.arange(neighbor_graph.number_of_edges())
label_dict = dict(zip(pair_number, labels))
edges = np.array(neighbor_graph.edges)
edges[:, 0] = [label_dict[edges[:, 0][i]] for i in pair_number]
edges[:, 1] = [label_dict[edges[:, 1][i]] for i in pair_number]
for i in pair_number:
    if edges[i, 0] > edges[i, 1]:
        edges[i, 0], edges[i, 1] = edges[i, 1], edges[i, 0]
edges = pd.DataFrame(edges, columns=["node1", "node2"])
edges = edges.sort_values(by=['node1', 'node2'])
edges["type"] = [str(edges.iloc[i, 0])+"-"+str(edges.iloc[i, 1]) for i in pair_number]
print("Types of edges in your graph: {}".format(edges["type"].unique().shape[0]))

print("graphing edge distributions...")
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(80, 30))
sns.countplot(x="type", palette="ch:.25", ax=ax, data=edges)
plt.title("Edge Distributions")
plt.savefig(snakemake.output[1])
plt.show()

col_names = list(pd.read_csv(snakemake.input[4])["col_name"])
cube = pd.DataFrame(data, columns=col_names)
cube.to_csv(snakemake.output[2], index=False)

#t-SNE
print("forming tSNE...")
from sklearn.manifold import TSNE
from matplotlib import colors
tsne_obj = TSNE(n_components=2, random_state=1).fit_transform(data)
np.random.seed(100)
nodeTypes = list(range(max(labels)+1))
color = []
for i in nodeTypes:
    color.append('#%06X' % np.random.randint(0, 0xFFFFFF))
fig, ax = plt.subplots(1, figsize=(16, 16))
scatter = ax.scatter(tsne_obj[:, 0], tsne_obj[:, 1], c=labels, cmap=colors.ListedColormap(color))
legend = ax.legend(*scatter.legend_elements(num=max(labels)+1), loc="lower left", title="Classes")
ax.add_artist(legend)
plt.savefig(snakemake.output[3])

