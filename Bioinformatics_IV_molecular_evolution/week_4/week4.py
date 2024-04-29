from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
from coursera_bioinformatics.utils import protein_utils
import pyperclip
import re
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_1.week1 import *
import copy

class Node: 
    def __init__(self, node_str: str):
        self.node_str = node_str
    
    def set_node_value(self, value): 
        self.node_value = value

class Edge: 
    def __init__(self, start_node: Node, end_node: Node, weight: int=0):
        self.start_node = start_node
        self.end_node = end_node
        self.weight = weight
    
    def set_weight(self, weight): 
        self.weight = weight

class Graph:
    def __init__(self, nodes: list[Node], edges: list[Edge]): 
        self.nodes = nodes
        self.edges = edges
    
    def get_and_set_adjacency_list(self): 
        adjacency_list = {}
        edges = copy.deepcopy(self.edges)
        for edge in edges:
            start_node = edge.start_node.node_str
            end_node = edge.end_node.node_str
            weight = edge.weight

            if not(start_node in adjacency_list): 
                adjacency_list[start_node] = {}

            adjacency_list[start_node][end_node] = weight
        
        self.adjacency_list = adjacency_list
        return adjacency_list

    def show_graph_in_adjacency_list(self, if_print=True): 
        to_print = ""
        if self.adjacency_list == None:
            pass
        for node in self.adjacency_list.keys(): 
            for next_node in self.adjacency_list[node].keys(): 
                temp = node + "->" + next_node + ": " + str(self.adjacency_list[node][next_node]) + "\n"
                to_print += temp
        
        to_print = to_print.strip()

        if if_print: 
            print(to_print)
        
        pyperclip.copy(to_print)
        return to_print
    
    # DFS find path
    def find_all_paths(self, start_node: str, end_node: str, path=[]): 
        path = path + [start_node]

        if start_node == end_node:
            return [path]

        paths = []

        if not(start_node in self.adjacency_list.keys()): 
            return []

        for node in self.adjacency_list[start_node].keys(): 
            if not(node in path): 
                new_paths = self.find_all_paths(node, end_node, path)
                for new_path in new_paths: 
                    paths.append(new_path)
        
        return paths
    
# the same algorithm in Bioinformatics III week 1
# the path length is refered to the total of all node_vals
# node value equals the value of the node-th value in the spectrum
# and the values are already stored in graph.node.node_value
def find_longest_path(graph: Graph, source: Node, sink: Node): 
    dp_table = {source.node_str: 0}
    backtrack_table = {source.node_str: None}

    for node in [item for item in graph.nodes if item.node_str != source.node_str]: 
        dp_table[node.node_str] = - float("inf")
        backtrack_table[node.node_str] = None

        for edge in [item for item in graph.edges if item.end_node.node_str == node.node_str]: 
            prev_node = edge.start_node
            new_length = dp_table[prev_node.node_str] + node.node_value
            if new_length > dp_table[node.node_str]: 
                dp_table[node.node_str] = new_length
                backtrack_table[node.node_str] = prev_node.node_str
    
    path = do_backtrack(backtrack_table, source.node_str, sink.node_str)
    
    return path

def do_backtrack(backtrack_table: dict, source: str, sink: str): 
    path = []
    node = sink
    path.append(node)
    while True and node != None: 
        node = backtrack_table[node]
        path.append(node)
        if node == source: 
            break
    
    path = [item for item in reversed(path)]
    return path

# Code Challenge: Construct the graph of a spectrum.
# Input: A space-delimited list of integers Spectrum.
# Output: Graph(Spectrum).
def construct_graph_of_ideal_spectrum(spectrum_str: str) -> Graph:
    spectrum = [int(item) for item in universal_utils.parse_arr(spectrum_str)]
    spectrum.append(0)
    spectrum = sorted(spectrum)
    spectrum_len = len(spectrum)

    edges = []
    nodes = [Node(str(node)) for node in spectrum]
    for i in range(spectrum_len): 
        for j in range(i + 1, spectrum_len): 
            weight = spectrum[j] - spectrum[i]
            if not(weight in protein_utils.ALL_AA_WEIGHTS): 
                continue
            
            start_node_value = str(spectrum[i])
            start_node = Node(str(start_node_value))
            start_node.set_node_value(start_node_value)
            end_node_value = str(spectrum[j])
            end_node = Node(str(end_node_value))
            end_node.set_node_value(end_node_value)

            edge = Edge(start_node, end_node, weight)
            edges.append(edge)

    graph = Graph(nodes, edges)
    adjacency_list = graph.get_and_set_adjacency_list()

    for node in adjacency_list.keys(): 
        for next_node in adjacency_list[node].keys(): 
            weight = adjacency_list[node][next_node]
            AA = protein_utils.get_AA_by_weight(weight)
            adjacency_list[node][next_node] = AA
    
    return graph

# return the ideal spectrum of a given peptide
# all the AAs in the peptide are expressed in a single letter
# an ideal spectrum is the collection of integer masses of all the prefixes and suffixes of a peptide
# i.e. assume every whole peptide only breaks once in the process of ionization, and subpeptides will not further break into smaller parts
def get_ideal_spectrum(peptide: str): 
    spectrum = []
    peptide_len = len(peptide)
    spectrum.append(protein_utils.get_weight_by_peptide(peptide))
    
    for i in range(1, peptide_len): 
        prefix_peptide = peptide[0 : i]
        suffix_peptide = peptide[i : peptide_len + 1]
        prefix_mass = protein_utils.get_weight_by_peptide(prefix_peptide)
        suffix_mass = protein_utils.get_weight_by_peptide(suffix_peptide)
        spectrum.append(prefix_mass)
        spectrum.append(suffix_mass)
    
    return sorted(spectrum)

# Input: A space-delimited list of integers Spectrum.
# Output: An amino acid string that explains Spectrum.
def construct_peptide_from_ideal_spectrum(spectrum_str): 
    graph = construct_graph_of_ideal_spectrum(spectrum_str)

    nodes = sorted([node.node_str for node in graph.nodes], key=lambda x: int(x))
    source_node = nodes[0]
    sink_node = nodes[-1]

    paths = graph.find_all_paths(source_node, sink_node, [])

    # universal_utils.print_arr(paths)

    peptides = []
    for path in paths: 
        peptide = ""
        for i in range(1, len(path)): 
            weight = int(path[i]) - int(path[i - 1])
            AA = protein_utils.get_AA_by_weight(weight)
            peptide += AA
        peptides.append(peptide)
    
    # peptides_with_correct_ideal_spectrum = []
    for peptide in peptides: 
        spectrum = [str(item) for item in get_ideal_spectrum(peptide)]
        if spectrum_str == " ".join(spectrum): 
            # peptides_with_correct_ideal_spectrum.append(peptide)
            return peptide
    
    # return peptides

# convert a peptide into a mass vector
# Input: An amino acid string Peptide.
# Output: The peptide vector Peptide'.
# a peptide vector refers to all the prefix masses of a peptide
def get_peptide_vector(peptide: str): 
    peptide_weight = protein_utils.get_weight_by_peptide(peptide)
    prefix_masses = []
    for i in range(1, len(peptide) + 1): 
        prefix_peptide = peptide[0 : i]
        prefix_mass = protein_utils.get_weight_by_peptide(prefix_peptide)
        prefix_masses.append(prefix_mass)

    peptide_vector = []
    for mass in range(1, peptide_weight + 1): 
        if mass in prefix_masses: 
            peptide_vector.append(1)
        else:
            peptide_vector.append(0)
    
    return peptide_vector

# Input: A binary vector P.
# Output: A peptide whose peptide vector is equal to P (if such a peptide exists).
def construct_peptide_from_vector(peptide_vector: list): 
    prefix_masses = []
    for i in range(len(peptide_vector)): 
        if peptide_vector[i] == 1: 
            prefix_masses.append(i + 1)
    
    peptide = ""
    for prefix_mass in prefix_masses: 
        current_AA_weight = prefix_mass - protein_utils.get_weight_by_peptide(peptide)
        peptide += protein_utils.get_AA_by_weight(current_AA_weight)

    return peptide

def peptide_sequencing_by_amplitude_spectral_vector(amplitude_spectral_vector: list[int]): 
    peptide_weight = len(amplitude_spectral_vector)
    vector = [0] + amplitude_spectral_vector

    min_AA_weight = min(protein_utils.ALL_AA_WEIGHTS)
    # connect node i to node j by a directed edge if j − i is equal to the mass of an amino acid
    edges = []
    nodes = []
    # possible_start_nodes = [0]
    for i in range(peptide_weight + 1 - min_AA_weight): 
        # if not(i in possible_start_nodes): 
        #     continue
        for j in range(i + min_AA_weight, peptide_weight + 1): 
            weight = j - i
            if weight in protein_utils.ALL_AA_WEIGHTS: 
                start_node = Node(str(i))
                end_node = Node(str(j))

                edge = Edge(start_node, end_node, weight)
                edges.append(edge)
                # possible_start_nodes.append(j)
    for i in range(peptide_weight + 1): 
        node = Node(str(i))
        node.set_node_value(vector[i])
        nodes.append(node)
    
    # print(edges)
    graph = Graph(nodes, edges)
    # graph = Graph(possible_start_nodes, edges)
    graph.get_and_set_adjacency_list()

    # there are too many paths from "0" to str(peptide_weight) in the DAG!
    # paths = graph.find_all_paths("0", str(peptide_weight), [])
    # peptides = []
    # for path in paths: 
    #     peptide = ""
    #     for i in range(1, len(path)): 
    #         weight = int(path[i]) - int(path[i - 1])
    #         AA = protein_utils.get_AA_by_weight(weight)
    #         peptide += AA
    #     peptides.append(peptide) 
    
    # # universal_utils.print_arr(peptides)
    # scores_peptides_mapping = {}
    # for peptide in peptides: 
    #     score = score_peptide_to_amplitude_spectral_vector(peptide, amplitude_spectral_vector)
    #     if not(score in scores_peptides_mapping): 
    #         scores_peptides_mapping[score] = []
    #     scores_peptides_mapping[score].append(peptide)
    
    # max_score = max(scores_peptides_mapping.keys())
    # return scores_peptides_mapping[max_score][0]

    # use the same algorithm in Bioinformatics III week 1
    path = find_longest_path(graph, nodes[0], nodes[-1])
    peptide = ""
    for i in range(1, len(path)): 
        weight = int(path[i]) - int(path[i - 1])
        AA = protein_utils.get_AA_by_weight(weight)
        peptide += AA
    return peptide

# score a peptide against a spectrum
def score_peptide_to_amplitude_spectral_vector(peptide: str, amplitude_spectral_vector: list[int]): 
    score = 0
    if protein_utils.get_weight_by_peptide(peptide) != len(amplitude_spectral_vector): 
        return - float("inf")
    
    peptide_vector = get_peptide_vector(peptide)

    for i in range(len(peptide_vector)): 
        if peptide_vector[i] == 1: 
            score += amplitude_spectral_vector[i]

    return score


if __name__ == "__main__": 
    # 1.3 step 5
    # spectrum_str = "57 71 154 185 301 332 415 429 486"
    # graph = construct_graph_of_ideal_spectrum(spectrum_str)
    # graph.show_graph_in_adjacency_list()

    # with open("./datasets/dataset_30262_5.txt", "r") as f: 
    #     spectrum_str = f.read().strip()
    #     graph = construct_graph_of_ideal_spectrum(spectrum_str)
    #     graph.show_graph_in_adjacency_list()

    # 1.3 step 7
    # peptide = "REDCA"
    # peptide = "GGDTN"
    # universal_utils.print_arr(get_ideal_spectrum(peptide))

    # spectrum_str = "57 71 154 185 301 332 415 429 486"
    # print(construct_peptide_from_ideal_spectrum(spectrum_str))

    # extra dataset
    # spectrum_str = "103 131 259 287 387 390 489 490 577 636 690 693 761 840 892 941 1020 1070 1176 1198 1247 1295 1334 1462 1481 1580 1599 1743 1762 1842 1861 2005 2024 2123 2142 2270 2309 2357 2406 2428 2534 2584 2663 2712 2764 2843 2911 2914 2968 3027 3114 3115 3214 3217 3317 3345 3473 3501 3604"
    # print(construct_peptide_from_ideal_spectrum(spectrum_str))

    # with open("./datasets/dataset_30262_8.txt", "r") as f: 
    #     spectrum_str = f.read().strip()
    #     peptide = construct_peptide_from_ideal_spectrum(spectrum_str)
    #     print(peptide)

    # 1.5 step 5
    # peptide = "XZZXX"
    # universal_utils.print_arr(get_peptide_vector(peptide))

    # with open("./datasets/dataset_30264_5.txt", "r") as f: 
    #     peptide = f.read().strip()
    #     pyperclip.copy(universal_utils.print_arr(get_peptide_vector(peptide)))

    # 1.5 step 6
    # peptide_vector = [int(item) for item in universal_utils.parse_arr("0 0 0 1 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 0 1")]
    # peptide = construct_peptide_from_vector(peptide_vector)
    # print(peptide)
    # pyperclip.copy(peptide)

    # with open("./datasets/dataset_30264_6.txt", "r") as f: 
    #     peptide_vector = [int(item) for item in universal_utils.parse_arr(f.read().strip())]
    #     peptide = construct_peptide_from_vector(peptide_vector)
    #     print(peptide)
    #     pyperclip.copy(peptide)

    # 1.5 step 13
    # amplitude_spectral_vector = [int(item) for item in universal_utils.parse_arr("0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8")]
    # peptide = peptide_sequencing_by_amplitude_spectral_vector(amplitude_spectral_vector)
    # print(peptide)

    with open("./datasets/dataset_30264_6.txt", "r") as f: 
        amplitude_spectral_vector = [int(item) for item in universal_utils.parse_arr(f.read().strip())]
        peptide = peptide_sequencing_by_amplitude_spectral_vector(amplitude_spectral_vector)
        print(peptide)
