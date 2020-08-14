#! /usr/bin/env python

"""
# ----------------------------------------------------------------------------------------------
# encoding: utf-8
# ReadData.py
# Susan D. Ghiassian
# Last Modified: 2020-08-14
# This code read the input data to be used in the network-based ranking of drugs as described in:
# 
# A systems-based approach to identify potential antivirals with a COVID-19 showcase
#
# by Mengran Wang, Johanna Withers, Piero Ricchiuto, Michael McAnally, 
# Helia Sanchez, Alif Saleh, Slava Akmaev and Dina Ghiassian 
# 
# ----------------------------------------------------------------------------------------------
"""

import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import pickle as pcl
import time
import scipy

def load_S2E():
# This function converts gene symbol to Entrez ID
    S2S = pcl.load(open('%s/Symbol2Synonyms.pcl'%data_dir,'rb'))
    E2S = pcl.load(open('%s/Entrez2Symbol.pcl'%data_dir,'rb'))
    S2E = {}
    for e,s in E2S.items():
        S2E[s] = e
        syns = S2S[s]
        if syns == '-':
            continue
        for syn in syns:
            S2E[syn] = e
    return S2E

def load_HI():
    # Read and generate the network

    file=open('./%s/Human_Interactome.txt'%data_dir,'r') #Edge list containing at least two columns representing two interacting partners
    initial_data=file.read().splitlines()[1::]
    file.close()
    G_o = nx.Graph()
    for row in initial_data:
        n=row.strip().split('\t')

        G_o.add_node(int(n[0]))
        G_o.add_node(int(n[1]))
        G_o.add_edge(int(n[0]), int(n[1]))
    
    # Isolate the largest connected component for follow up analysis
    Clusters=sorted(list(nx.connected_component_subgraphs(G_o)), key=len, reverse=True)
    G=Clusters[0]
    G_nodes=G.nodes()
    
    # Remove self loops
    self_loops = G.selfloop_edges()
    G.remove_edges_from(self_loops)
    return HI

def load_drugbank(target_type = 0):
    #target type
    # 0: include all drugs all targets
    # 1: include only drugs with at least one human polypeptide and all targets
    # 2: include only drugs with al least one human polypeptide and only polypeptide targets

    from collections import defaultdict
    import pandas as pd

    df = pd.read_csv('%s/Data/DrugTarget/Drugbank.csv'%curr_dir)
    data= zip(df['ID'].values, df['entrez_id'].values,df['Type'].values)
    drug2target2type=defaultdict(dict)
    for drug,entrez,t_type in data:
        try:
            #drug,entrez,t_type = line[0],line[1],line[2]
            drug2target2type[drug][int(entrez)] = t_type
        except ValueError: continue

    if target_type == 0 :
        drug2target = dict([(d,set(t2t.keys())) for d,t2t in drug2target2type.items()])
    if target_type == 1:

        drug2target = dict([(d,set(t2t.keys())) for d,t2t in drug2target2type.items() if 'Polypeptide' in set(t2t.values())])
    if target_type == 2:
        drug2target = dict([(d,set([tar for tar,tar_type in t2t.items() if tar_type == 'Polypeptide' ])) for d,t2t in drug2target2type.items() if 'Polypeptide' in set(t2t.values())])
    return drug2target

