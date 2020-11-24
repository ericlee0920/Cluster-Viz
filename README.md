# Expression Cluster Visualization
## Updated: November 24, 2020.
This pipeline performs clustering on high dimensional expression single-cell data in a non-spatial context as well as the visualization of the clusters in network graphs and hierarchical clustering heatmaps. 

### Background and Rationale
•	include the what’s and why’s – also your aims
•	include any package dependencies that are required (bullet points are ok for this)

### DAG
<img src="https://github.com/ericlee0920/Cluster-Viz/blob/main/DAG.png?raw=true" width="300" height="300">

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
snakemake -cores 1 results/heatmap.png
```
5. Generate a DAG from the workflow:
```
snakemake --dag results/heatmap.png | dot -Tsvg > dag.svg
```

### Input
Describe the format of the input data, explaining all fields.


### Output
Describe the format of the output including files and visualizations. Treat this section like the results of a paper. You can look at readthedocs pages of popular bioinformatics tools to get inspired for this.
