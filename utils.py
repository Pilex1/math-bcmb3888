import networkx as nx
import numpy as np
import scipy as sp
from networkx.algorithms import community
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community import k_clique_communities
from community import community_louvain
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
from sklearn import metrics

def load_network():

    G0 = nx.read_weighted_edgelist("4932.protein.links.v11.5.txt",comments="#",nodetype=str)
    # print(f"number of nodes in original dataset: ", len(G0.nodes))

    #removing the prefix in proteins
    protein_info = pd.read_csv("Protein_info.txt", sep='\t')
    map_dic = protein_info.set_index('#string_protein_id').to_dict()['preferred_name']

    G = nx.relabel_nodes(G0, map_dic)

    # remove essential proteins
    essential_proteins = pd.read_csv("yeast essential proteins.csv", header=None)[1]
    # print()
    # print(essential_proteins)
    G.remove_nodes_from(essential_proteins)
    # print(f"number of nodes after removing essential proteins: ", len(G.nodes))  

    # delete those edges with a combined score of <= threshold_score (small confidence)
    threshold_score = 500
    for edge in G.edges: 
        weight = list(G.get_edge_data(edge[0],edge[1]).values())
        if(weight[0] <= threshold_score):
            G.remove_edge(edge[0],edge[1])

    return G
        