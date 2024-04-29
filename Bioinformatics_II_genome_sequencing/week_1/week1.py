import pyperclip
import random
import time
from itertools import product

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

# get all possible composition of a strand
def get_kmer_composition(strand, k):
    results = []
    for i in range(len(strand) - k + 1):
        current_kmer = strand[i : i + k]
        results.append(current_kmer)
    
    return results


# reconstruct a string using every k-mer element in a composition list
# the first k-mer overlaps k - 1 nucleotides with the second, and the second k-mer overlaps k - 1 nucleotides with the third, ... (if in the correct order)
# but the order of the k-mers is not what they appear to construct the string
def string_construction(composition, k):
    # first, find the first kmer to construct the string
    first_kmer_list = []
    for kmer in composition: 
        other_kmers = [item for item in composition if item != kmer]
        is_first = True
        # 当前遍历到的kmer，如果没有其他kmer以它的前k-1位结束，就说明当前kmer可能是第一个
        for other_kmer in other_kmers:
            if kmer[0 : k - 1] == other_kmer[1 : k]:
                is_first = False
        
        if is_first:
            first_kmer_list.append(kmer)
    
    # print(first_kmer_list)
    # TAA
    # but the first kmer can be one of all the kmers

    result_list = []

    ideal_len = len(composition) + k - 1

    # second, start from each first_kmer, try to construct the entire string
    for first_kmer in first_kmer_list:
        # 暂存result
        result = first_kmer
        other_kmers = [item for item in composition if item != first_kmer]

        prev_result_len = 0
        # 如果走到最后发现result的长度不是正确的长度，就打乱重来
        # 感觉不是dp啊。。。但是有用！
        while len(result) != ideal_len:
            current_result_len = len(result)
            for other_kmer in other_kmers:
                if result[current_result_len - (k - 1) :] == other_kmer[0 : k - 1]:
                    result += other_kmer[-1]
                    other_kmers.remove(other_kmer)
            
            # print(result)
            # print(other_kmers)
            # print()

            prev_result_len = current_result_len
            current_result_len = len(result)
            # 如果长度不再能增加，就重置other-kmers数组并打乱
            if prev_result_len == current_result_len:
                other_kmers = [item for item in composition if item != first_kmer]
                result = first_kmer
                random.shuffle(other_kmers)

        result_list.append(result)

    return result_list


# input: a collection Patterns of k-mers
# output: the overlap graph Overlap(Patterns), in the form of an adjacency list
# the nodes and their edges can be returned in any order
def get_overlap_graph(kmers):
    overlap_graph = {}
    k = len(kmers[0])

    for kmer in kmers:
        # traverse all the kmers
        # add current kmer as a key, and the value is a set
        overlap_graph[kmer] = set()
        suffix = kmer[1 :]

        # traverse all the other kmers
        for other_kmer in [item for item in kmers if item != kmer]:
            prefix = other_kmer[: k - 1]

            # add current other_kmer to the set if the suffix overlaps with prefix
            if prefix == suffix:
                overlap_graph[kmer].add(other_kmer)
    
    return overlap_graph

# print dict-like graph
def print_graph_dict(graph_dict):
    if graph_dict is None:
        return
    
    result = ""
    for key in graph_dict.keys():
        if len(graph_dict[key]) == 0:
            continue

        temp = key + ": "

        for kmer in graph_dict[key]:
            temp += kmer + " "
        
        temp = temp.strip()
        print(temp)
        result += temp + "\n"
    
    result = result.strip()
    pyperclip.copy(result)

# input: an int k and a strand
# output: DeBruijnk(strand), in the form of an adjacency list
# len(node) should be k - 1 in de Bruijn graph
def get_de_Bruijn_graph_by_strand(strand, k):
    if k > len(strand):
        print("invalid k bigger than len(strand)")
        return None

    de_Bruijn_graph = {}
    
    for i in range(len(strand) - (k - 1) + 1):
        current_node = strand[i : i + k - 1]

        # the node which current node has a path to it
        next_node = ""
        if i != (len(strand) - (k - 1)):
            next_node = strand[i + 1 : i + k]

        # current node should be unique in the graph, so it is added into dict.keys()
        if current_node not in de_Bruijn_graph.keys():
            de_Bruijn_graph[current_node] = []

        # but the node which current node has a path to it can be duplicated, so next node is added into a list
        if (next_node != ""):
            de_Bruijn_graph[current_node].append(next_node)

    return de_Bruijn_graph

# input: a collection of k-mers composition
# output: the adjacency list of the de Bruijn graph DeBruijn(composition)
def get_de_Bruijn_graph_by_composition(composition):
    k = len(composition[0])

    # first, generate all the prefix nodes of the de Bruijn graph
    # a node should be a (k-1)-mer
    nodes = set()
    for kmer in composition:
        nodes.add(kmer[: k - 1])
        nodes.add(kmer[1 :])
    
    # second, generate the de Bruijn graph
    de_Bruijn_graph = {}
    for node in nodes:
        # traverse all the prefix nodes
        de_Bruijn_graph[node] = []

        for kmer in composition:
            # connect the prefix node to its suffix node by a directed edge
            prefix = kmer[: k - 1]
            suffix = kmer[1 :]

            if prefix == node:
                de_Bruijn_graph[node].append(suffix)
    
    return de_Bruijn_graph



if __name__ == "__main__":
    # 1.2 step 3
    # print_arr(get_kmer_composition("CAATCCAAC", 5))
    # with open("dataset_197_3.txt", "r") as f:
    #     k = int(f.readline().strip())
    #     strand = f.readline().strip()
    #     results = get_kmer_composition(strand, k)
    #     print_arr(get_kmer_composition(strand, k))
    #     result = ""
    #     for kmer in results:
    #         result += kmer + " "
    #     result.strip()
    #     pyperclip.copy(result)

    # 1.2 step 6
    # composition = "AAT  ATG  ATG  ATG  CAT  CCA  GAT  GCC  GGA  GGG  GTT  TAA  TGC  TGG  TGT".split("  ")
    # print(string_construction(composition, 3))

    # 1.3 step 3
    # with open("dataset_198_3.txt", "r") as f:
    #     composition = f.read().strip().split(" ")
        # print(string_construction(composition, len(composition[0])))
        # result = composition[0]
        # k = len(result)
        # for item in composition[1 :]:
        #     result += item[-1]
        
        # print(result)

        # print(string_construction(composition, len(composition[0])))


    # 1.3 step 10
    # kmers = "ATGCG GCATG CATGC AGGCA GGCAT GGCAC".split(" ")
    # print_graph_dict(get_graph_dict(kmers))
    # with open("dataset_198_10.txt", "r") as f:
    #     kmers = f.read().strip().split(" ")
    #     print_graph_dict(get_graph_dict(kmers))

    # 1.3 step 13 Construct a 4-universal string.
    # kmers = ["".join(p) for p in product("01", repeat=4)]
    # k = len(kmers[0])
    # ideal_len = len(kmers) + k - 1
    # result_list = set()
    # # 现在任何一个kmer都可以作为开头
    # for kmer in kmers:
    #     result = kmer
    #     other_kmers = [item for item in kmer if item != kmer]

    #     prev_result_len = 0
    #     while len(result) != ideal_len:
    #         current_result_len = len(result)
    #         for other_kmer in other_kmers:
    #             if result[current_result_len - (k - 1) :] == other_kmer[0 : k - 1]:
    #                 result += other_kmer[-1]
    #                 other_kmers.remove(other_kmer)
            
    #         prev_result_len = current_result_len
    #         current_result_len = len(result)
    #         # 如果长度不再能增加，就重置other-kmers数组并打乱
    #         if prev_result_len == current_result_len:
    #             other_kmers = [item for item in kmers if item != kmer]
    #             result = kmer
    #             random.shuffle(other_kmers)
    #     result_list.add(result)
    # # 竟然每个元素作为开头都能构建出一个4-universal string
    # print(result_list)

    # 1.4 step 6 
    # strand = "AAGATTCTCTAAGA"
    # k = 4
    # print_graph_dict(get_de_Bruijn_graph(strand, k))
    # with open("dataset_199_6.txt", "r") as f:
    #     k = int(f.readline().strip())
    #     strand = f.readline().strip()

    #     print_graph_dict(get_de_Bruijn_graph(strand, k))

    # 1.4 step 7
    # strand = "TAATGCCATGGGATGTT"
    # print("de Bruijn k = 2")
    # print_graph_dict(get_de_Bruijn_graph(strand, 2))
    # print("de Bruijn k = 3")
    # print_graph_dict(get_de_Bruijn_graph(strand, 3))
    # print("de Bruijn k = 4")
    # print_graph_dict(get_de_Bruijn_graph(strand, 4))
    # print("de Bruijn k = 5")
    # print_graph_dict(get_de_Bruijn_graph(strand, 5))

    # strand1 = "TAATGCCATGGGATGTT"
    # strand2 = "TAATGGGATGCCATGTT"
    # print("de Bruijn graph of TAATGCCATGGGATGTT:")
    # print_graph_dict(get_de_Bruijn_graph_by_strand(strand1, 3))
    # print("de Bruijn graph of TAATGGGATGCCATGTT:")
    # print_graph_dict(get_de_Bruijn_graph_by_strand(strand2, 3))

    # 1.5 step 10
    # composition = "GAGG CAGG GGGG GGGA CAGG AGGG GGAG".split(" ")
    # print_graph_dict(get_de_Bruijn_graph_by_composition(composition))
    # with open("dataset_200_8.txt", "r") as f:
    #     composition = f.read().strip().split(" ")

    #     print_graph_dict(get_de_Bruijn_graph_by_composition(composition))

    # week 1 Quiz
    # problem 4
    print("problem 4")
    k = 3
    binary_strings = "1001101100 0101010100 1011100010 0111010010 0111010001 0011101000".split(" ")
    for binary_string in binary_strings:
        composition = set(get_kmer_composition(binary_string, k))
        ideal_length = len(binary_string) - k + 1
        if len(composition) == ideal_length:
            print(binary_string)

    # problem 5
    print("\nproblem 5")
    kmers = \
"""GCGA
CAAG
AAGA
GCCG
ACAA
AGTA
TAGG
AGTA
ACGT
AGCC
TTCG
AGTT
AGTA
CGTA
GCGC
GCGA
GGTC
GCAT
AAGC
TAGA
ACAG 
TAGA
TCCT
CCCC
GCGC
ATCC
AGTA
AAGA
GCGA
CGTA""".split("\n")
    print(kmers)

    count = 0
    for kmer in kmers:
        if kmer[: 4 - 1] == "AAG":
            count += 1
        
    print("outdegree of AAG:", count)
    