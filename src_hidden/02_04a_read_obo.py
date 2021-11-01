#!/usr/bin/env python
import networkx as nx
import re
import os
import argparse
import pickle

'''
This was mdified from embed_ontologies to only be for doid
'''

def get_terms(OBO_file):
    with open(OBO_file, 'r') as f:
        text = f.read()

    terms = [term for term in text.split('\n\n')
             if term.startswith('[Term]')]             
    return terms
    

def get_graph(terms):

    # ID terms have at least 1 capital letter, a colon, and at least 1 digit
    ID_pattern = re.compile('[A-Za-z]+:\d+')

    graph = nx.DiGraph()
    ID_name_map = {}

    for idx, term in enumerate(terms):
        if 'is_obsolete: true' in term:
            continue  # skip obsolete terms
        
        lines = term.split('\n')
        for line in lines:
            if line.startswith('id:'):
                ID = ID_pattern.search(line).group(0)

            elif line.startswith('name:'):
                ID_name_map[ID] = line.split(': ')[1]

            elif line.startswith('is_a:'):
                parent = ID_pattern.search(line).group(0)
                graph.add_edge(parent, ID)
                        
            #Ignoring 'develops from' and 'related to'
            elif (line.startswith('relationship: part_of')):
                parent = ID_pattern.search(line).group(0)
                graph.add_edge(parent, ID)

    simple_cycles = list(nx.simple_cycles(graph))
    print('The simpe cycles are',simple_cycles)
    if not nx.is_directed_acyclic_graph(graph):
        raise nx.NetworkXException('graph is not acyclic')

    return graph, ID_name_map


if __name__ == '__main__':

    fp = '/mnt/research/compbio/krishnanlab/projects/GenePlexus/repos/GenePlexusBackend/data_hidden2/GSCs/'
    FN = fp + 'doid.obo'

        
    terms = get_terms(FN)
    graph, ID_name_map = get_graph(terms)
    

    print(len(list(graph.nodes())))
    print(len(ID_name_map))
    print(list(graph.edges())[0:10])
    
    with open(fp+'doid_terms.tsv', 'w') as f:
        for ID in ID_name_map:
            f.write('{}\t{}\n'.format(ID, ID_name_map[ID]))

    with open(fp+'doid_ancestors.tsv', 'w') as f:
        for node in graph.nodes():
            ancestors = nx.ancestors(graph, node)
            if (len(ancestors) > 0):
                f.write('{}\t{}\n'.format(node, ', '.join(ancestors)))

    with open(fp+'doid_descendants.tsv', 'w') as f:
        for node in graph.nodes():
            descendants = nx.descendants(graph, node)
            if (len(descendants) > 0):
                f.write('{}\t{}\n'.format(node, ', '.join(descendants)))
