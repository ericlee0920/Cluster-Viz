rule data_extraction:
    input:
        "data/out.csv",
    output:
        temp("results/expression_matrix.npy"),
        temp("results/location_matrix.npy"),
        temp("results/axis_matrix.npy")
    script:
        "extraction.py"

rule community_clustering:
    input:
        "results/expression_matrix.npy",
    output:
        temp("results/communities.npy")
    script:
        "community_analysis.py"

rule spatial_graph_and_edge_distribution:
    input:
        "results/expression_matrix.npy",
        "results/location_matrix.npy",
        "results/axis_matrix.npy",
        "results/communities.npy",
        "data/column_names.csv"
    output:
        "results/neighbor_graph.png",
        "results/edge_distribution.png",
        temp("results/cube.csv"),
        "results/tSNE.png"
    script:
        "spatial_graph.py"

rule hierarchical_heatmap:
    input:
        "results/cube.csv",
        "results/communities.npy"
    output:
        "results/heatmap.png"
    script:
        "hierarchical_heatmap.py"
