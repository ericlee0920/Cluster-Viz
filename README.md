# Expression Cluster Visualization
## Updated: November 26, 2020.
### Background and Rationale
**Cluster-Viz** is a workflow for the visualization of single cell expression data, specifically from [Imaging Mass Cytometry™](https://www.fluidigm.com/applications/imaging-mass-cytometry) (IMC). IMC is an extension of mass cytometry that preserves spatial context, rendering not only protein expression data but also spatial coordinates. High dimensional expression data when coupled with spatial information allows the intricate analysis of cell to cell interactions. It is our belief that plotting by spatial coordinates improve beyond visualizations using dimension reduction such as [t-SNE](https://lvdmaaten.github.io/tsne/). The objective of Cluster-Viz addresses the need to visualize high dimensional data and extract the underlying topological structure for a single patient sample.
This workflow can also be generalized to other single cell expression data with changes to the input data format which is described below in detail. 

The underlying approach to this workflow has four modules:
- Data Extraction: pre-processing the data for downstream analysis.
- Community Clustering: using Phenograph to cluster single cells to clusters of similar cell types based on expression values.
- Spatial Graph and Edge Distribution: produce cell type labelled spatial graphs, a bar graph of edge distributions, and a t-SNE plot for comparison.
- Hierarchical Heatmap: produce a z-scored mean marker hierarchical heat map.

<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/DAG.png?raw=true" width="300" height="300">

Package Dependencies:
  - snakemake-minimal = 5.29.0
  - jinja2 = 2.11.2
  - networkx >= 2.5
  - matplotlib = 3.3.3
  - graphviz > 2.40
  - numpy = 1.19.4
  - pandas = 1.1.4
  - scipy = 1.5.3
  - scikit-learn > 0.23
  - seaborn = 0.11.0
  - numba = 0.51.2
  - pip = 20.2.4
  - pip:
    - Phenograph (>=1.5.7) should be most recent

### Usage
The execution of the following files requires a Linux/Unix-like OS. If you are using Windows, please make sure you are setting up a VM.

1. Clone this file in a proper directory. This will download all the files and datasets necessary for execution. The output will be a folder named "Cluster-Viz".
```
git clone https://github.com/ericlee0920/Cluster-Viz.git
```
2. Create an environment with all the dependencies. This will manage all the packages under the same collection for execution.
```
conda env create --name cluster --file environment.yaml
```
3. Activate the environment.
```
conda activate cluster
```
4. Now run Cluster-Viz by executing the following code. This will create three files in the results folder: `neighbor_graph.png`, `edge_distribution.png`, `heatmap.png`.
```
snakemake --cores 1 results/heatmap.png
```
5. Generate a DAG from the workflow:
```
snakemake --dag results/heatmap.png | dot -Tsvg > dag.svg
```

### Input
Cluster-Viz takes in two csv files: a data cube for a single patient sample and a marker list. All the files should be placed in the folder `data/`. 

The sample dataset here is a processed subset of a breast cancer IMC dataset for a single patient. All protein expressions and relevant metadata is in the data cube called `out.csv` and a its marker list in `column_names.csv`. The metadata is described below.

Data cubes should have the following 13 properties as columns:
*If you are using expression data from other technologies and platforms, the bolded features below are required in the csv. If the data you are working with do not differentiate the source of intensity, use the columns flux_01-50 and set all of f_buffer_01-50 to 0.*

- X, Y: coordinates with the origin at the bottom-left of the IMC image.
- X_image, Y_image: coordinates with the origin at the top-left of the IMC image.
- area_cnt, area_minCircle, area_ellipse: morphological parameters, such as count of pixels in cell, minimal area with cell as a circle, area with cell as an ellipse.
- ell_angle, ell_smaj, ell_smin, ell_e: morphological parameters with cell as an ellipse, including the angle, the major axis, the minor axis, and the eccentricity.
- flux_01-50: intensity of different IMC markers in the nuclear area of the cell. The number of channels refer to the column names in rows in `column_names.csv`.
- f_buffer_01-50: intensity of different IMC markers outside of the nuclear area of the cell. The number of channels refer to the column names in rows in `column_names.csv`.

Marker list should have a single column of channel names in the dataset corresponding to the number of flux or flux buffer columns in the data cube. Row 1 of this file refers to flux_01 and f_buffer_01. 

### Output
Cluster-Viz produces three types of visualizations in png files for interpretation of expression data: labelled spatial graphs, edge type distribution plots, and z-scored mean marker hierarchical heat maps. All the files should be placed in the folder `result/`. Intermediate files are temporary and are removed.

- Labelled spatial graphs: This graph shows the distribution of cell types in their exact spatial coordinates. Nodes refer to cells, and edges connect cells that are overlapping. The colors refer to the cell type cluster determined by Phenograph. We can see same cell types likely to clump together, which is a common behavior in cancer images.

**NOTE**: Phenograph is a stochastic method using the Louvain modularity algorithm to construct clusters. The labels are not deterministic in nature, thus we see label switching. Also, due to multiple subclones present in this image and also being a subset of a larger graph, there are possibly a few nodes that may be clustered differently in each run. Phenograph also has a tendency to overcluster, the overclustering will be resolved by looking at the heatmap produced.

Run #1  | Run #2
------------- | -------------
<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run1_neighbor_graph.png?raw=true" width="500" height="500"> | <img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run2_neighbor_graph.png?raw=true" width="500" height="500">

- t-SNE plot: This is a visualization for high dimensional data based on Stochastic Neighbor Embedding. Axes do not refer to spatial coordinates. A t-SNE plot is included to show that comparing to the labelled spatial graphs, the labelled spatial graphs can allow users to visualize the inherent spatial structure better than a visualization that only considers expression values. We cannot identify that these two runs show the same expression values, yet in the labelled spatial graphs we can see the input is the same.

Run #1  | Run #2
------------- | -------------
<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run1_tSNE.png?raw=true" width="500" height="500"> | <img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run2_tSNE.png?raw=true" width="500" height="500">

- Edge type distribution plots: This graph shows the distributions of all the type of edges or interactions in the graph. The x axis refers to the type of edges in *_* format (e.g. 1_4) where numbers refer to cell type. In this plot, we can see that due to label switching, the edge interaction distribution changes. This graph indicates that the peaks show abundance of edges with the same cell type, such as 0_0 and 1_1. This is the same observation seen above in the labelled spatial graph.

Run #1  | Run #2
------------- | -------------
<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run1_edge_distribution.png?raw=true" width="500" height="500"> | <img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run2_edge_distribution.png?raw=true" width="500" height="500">

- Z-scored mean marker hierarchical heat maps: This shows the z-scored mean marker expression of all type of single cell phenotypic clusters identified in clustering. The numbers refer to cell type or cluster number. The bottom row refers to the marker names. Colors on the color bar refer to the measurement mean of each marker in a specific cell type. Here, we can see minor changes in marker intensity profile for each cluster. This is due to a few nodes that may be clustered differently in each run. This problem will be resolved with larger datasets. As for resolving overclustering from Phenograph, usually a manual merge by pathologists is needed to merge clusters that seem to be the same cell type together. This visualization supports pathologists to merge clusters.

Run #1  | Run #2
------------- | -------------
<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run1_heatmap.png?raw=true" width="450" height="450"> | <img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/sample_runs/run2_heatmap.png?raw=true" width="450" height="450">
