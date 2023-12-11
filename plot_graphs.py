import pickle
import os
import networkx as nx
import pandas as pd
from pathlib import Path
import json

def main():
    results_directory = r'ParamOptimizationResults\12_10_2023_06_57'
    res = None
    for subdir in os.listdir(results_directory):
        if str(subdir).startswith("gene_"):
            pkl_path = os.path.join(results_directory, subdir,'new_graph.pkl')
            print(subdir)
            print("Loading graph...")
            with open(pkl_path, 'rb') as f:
                G = pickle.load(f)
    
            if res is None:
                res = set([x for x in G.neighbors('APC')])
            else:
                print(res == set([x for x in G.neighbors('APC')]))

def calc_global_rank(data_dir):
    gene_scores = {}
    gene_patients = {}
    APC_alternatives = []
    for file in Path(data_dir).glob('*.csv'):
        df = pd.read_csv(file, index_col=0)
        for pathway in df.index:
            snv, weight = max(df.loc[pathway].items(), key=lambda item: item[1])
            gene_scores[snv] = gene_scores.get(snv, 0) + weight
            gene_patients[snv] = gene_patients.get(snv, set()).union({file.stem})
            if 'APC' in df.loc[pathway].keys() and df.loc[pathway]['APC'] > 0 and 'APC' != snv:
                APC_alternatives.append((snv, weight, df.loc[pathway]['APC']))
                 
    return (gene_scores, gene_patients, APC_alternatives)

if __name__ == '__main__':
    (gene_scores, gene_patients, APC_alternatives) = calc_global_rank(r'Data\DriverMaxSetApproximation\BaseData')
    # with open(r'ParamOptimizationResults\12_10_2023_06_57\gene_param=0.1_gene_penalty_patient_discount_param=0.1\sorted_gene_names_by_weight.json', 'r') as f:
    #     res = json.load(f)
    # print((sorted(gene_scores.items(), key=lambda x:x[1],reverse=True))[:10])
    for snv,w1,w2 in APC_alternatives:
        print(snv,w1,w2, gene_scores[snv], len(gene_patients[snv]))