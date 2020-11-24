import numpy as np
import scipy.sparse as ss
import phenograph
# """
# For a dataset of N rows, communities will be a length N vector of integers specifying a community assignment for each row in the data.
# Any rows assigned -1 were identified as outliers and should not be considered as a member of any $community$.
# $graph$ is a N x N scipy.sparse matrix representing the weighted graph used for community detection.
# Q is the modularity score for communities as applied to graph.
# """

print("forming clusters...")
FILENAME = snakemake.input[0]
data = np.load(FILENAME)
communities, graph, Q = phenograph.cluster(data)
np.save(snakemake.output[0], communities)
# ss.save_npz('graph.npz', graph)
print("\nThe modularity score of your data after clustering: {}\n".format(Q))
