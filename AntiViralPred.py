#! /usr/bin/env python

"""
# ----------------------------------------------------------------------------------
# encoding: utf-8
# AntiViralPred.py
# Susan D. Ghiassian
# Last Modified: 2020-08-14
# This code runs the network-based ranking of drugs as described in
# 
# A systems-based approach to identify potential antivirals with a COVID-19 showcase
#
# by Mengran Wang, Johanna Withers, Piero Ricchiuto, Michael McAnally, 
# Helia Sanchez, Alif Saleh, Slava Akmaev and Dina Ghiassian 
# 
# ------------------------------------------------------------------------------------
"""


import pandas as pd
import numpy as np
import networkx as nx
from collections import defaultdict
import pickle as pcl
import time
import scipy
import ReadData


def get_virus_host_graph(virus):
    #Specify virus as either 'HIV' or 'SARS-CoV-2'
    S2E = ReadData.load_S2E()
    VHG = nx.Graph() #virus-host graph
    if virus == 'SARS-CoV-2':
        f = '%s/SARS-COV2/bait_prey_high_confidence.xlsx'%data_dir
        Gordon_table = pd.read_excel(f)
        n = Gordon_table.shape[0] - 1 #first line is header
        for i in range(n): 
            Bait = list(Gordon_table.loc[i+1])[0] #skipping first line
            Prey = list(Gordon_table.loc[i+1])[2]
            try:
                Prey = S2E[Prey] 
            except KeyError:
                continue
            VHG.add_edge(Bait,Prey) 
    if virus == 'HIV':
        f = '%s/HIV/NCBI/HIV-1_physical_interactions.tsv'%data_dir
        NCBI_table = pd.read_table(f,delimiter='\t')  
        n = NCBI_table.shape[0] - 1 #first line is header
        for i in range(n):
            if list(NCBI_table.loc[i+1])[3] >= 2: #include only the ones with at leasst two interactions
                v = list(NCBI_table.loc[i+1])[0]
                h = list(NCBI_table.loc[i+1])[1]
                VHG.add_edge(v,str(h))
    
    return(VHG)


def build_three_layer_network(drug_source,virus,target_type,HI=None,VHG=None):
    S2E = ReadData.load_S2E()    

    if HI==None:
        HI = ReadData.load_HI()
    multi_graph = HI.copy()
    if VHG == None:
        VHG = get_virus_host_graph(virus)
    for edge in VHG.edges():
        multi_graph.add_edge(edge[0],edge[1])
    D2T = ReadData.load_drugbank(target_type = target_type)
    for drug,targets in D2T.items():
        for target in targets:    
            multi_graph.add_edge(drug,str(target))
    return(multi_graph)        


def dense(double_dict):
    from collections import defaultdict
    double_dict_dense = defaultdict(dict)
    for k,v2v in double_dict.items(): 
        for v1,v2 in v2v.items(): 
            if v2!=0: 
                double_dict_dense[k][v1]=v2 
    return(double_dict_dense)

def L3(G):
    #Clusters=sorted(list(nx.connected_components(G)), key=len, reverse=True)
    #G_nodes = list(Clusters[0])
    #G = nx.subgraph(G,G_nodes)
    
    # This is to trim the network and speed up the calculations
    #G2 = nx.Graph(G) #make an unfrozen copy of G so that links can be removed from it
    #self_loops = nx.selfloop_edges(G) 
    #G2.remove_edges_from(self_loops)
    
    import numpy as np
    #import time
    #start_time = time.time()
    G_nodes = G.nodes()
    A = nx.adjacency_matrix(G, nodelist=G_nodes)
    D = np.diag(1./np.sqrt(np.asarray(A.todense()).sum(axis=0)))
    D_sparse = scipy.sparse.csr_matrix(D)
    
    T = A.dot(D_sparse).dot(A).dot(D_sparse).dot(A)
    #print("--- %s seconds ---" % (time.time() - start_time))

    #file=open('similarities.txt','w')
    #np.savetxt('similarities.txt', T.todense())
    return(T)





if __name__ == "__main__":
    data_dit = './Data'
    virus = 'SARS-CoV-2' #'HIV'
    drug_source = 'DrugBank' 
    ranking_version = 0 #0,1,2    
    
    #Build network
    HI = ReadData.load_HI()
    HI_nodes = HI.nodes()
    
    VHD = build_three_layer_network(drug_source,virus,target_type = ranking_version,HI = HI,VHG=None)
    VHD_nodes = list(VHD.nodes())
    
    if virus == 'SARS-CoV-2':
        f = '%s/SARS-COV2/bait_prey_high_confidence.xlsx'%data_dir
        Gordon_table = pd.read_excel(f)
        viral_proteins = set(list(Gordon_table.iloc[:,0])[1:])
        
    if virus == 'HIV':
        f = '%s/HIV/NCBI/HIV-1_physical_interactions.tsv'%data_dir
        NCBI_table = pd.read_table(f,delimiter = '\t')
        viral_proteins = set(list(NCBI_table.iloc[:,0])[1:])
        viral_proteins = list(viral_proteins & set(VHD_nodes))        

    

    D2T = ReadData.load_drugbank(target_type = ranking_version)
    drugs = list(D2T.keys())


    T = L3(VHD) #Build score matrix


    drug_ind,viral_ind = [],[]
    for d in drugs:
        if d in VHD_nodes:
            drug_ind.append(VHD_nodes.index(d))
    for v in set(viral_proteins) & set(VHD_nodes):
        viral_ind.append(VHD_nodes.index(v))

            
    D2V = defaultdict(dict)
    for d in drugs:
        if d not in VHD_nodes:
            continue
        i = VHD_nodes.index(d)
        for v in viral_proteins:
            j = VHD_nodes.index(v)
            D2V[d][v] = T[i,j]

#    
    with open('%s_virus_drug_norml3_matrix_%s_hc_2evidence_%d_ranking.txt'%(virus,drug_source,ranking_version),'w') as fout: 
        fout.write('\t'.join(['DrugID']+list(viral_proteins))+'\n') 
        for drug,v2s in D2V.items():
            s_list = []
            for virus,norml3 in v2s.items():
                s_list.append(norml3) 
            fout.write('\t'.join([drug]+list(map(str,s_list)))+'\n') 
            


