from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
import pyperclip
import re
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_1.week1 import *
import copy

pattern_num = r'^[0-9]+$'

# return an unweighed, directed tree, i.e. from root to leaves
def parse_directed_tree_with_mapping(edges_str): 
    edges = edges_str.split("\n")
    tree = Tree("")
    leaves_node_count = 0
    mapping = {}
    for edge in edges: 
        nodes = edge.split("->")
        node_0 = nodes[0]
        node_1 = nodes[1]

        # 如果这个node不是数字，则是需要通过mapping进行index和具体名称的映射
        if not(re.match(pattern_num, node_0)): 
            temp = node_0
            node_0 = str(leaves_node_count)
            leaves_node_count += 1
            mapping[node_0] = temp
        if not(re.match(pattern_num, node_1)): 
            temp = node_1
            node_1 = str(leaves_node_count)
            leaves_node_count += 1
            mapping[node_1] = temp

        if not(node_0 in tree.adjacency_list.keys()): 
            tree.adjacency_list[node_0] = {}
        if not(node_1 in tree.adjacency_list.keys()): 
            tree.adjacency_list[node_1] = {}
        
        # add unweighed, directed edge
        tree.adjacency_list[node_0][node_1] = 0
        # tree.adjacency_list[node_1][node_0] = 0

    # tree.refresh_tree()
    
    return tree, mapping

# transfer a undirected tree (i.e. two edges node_0->node_1 and node_1->node_0 both exist)
# into a directed tree, i.e. if node_0 is a leaf, node_1 is an internal node, then only the edge node_1->node_0 exists
def transfer_tree_to_directed(tree: Tree): 
    most_inner_node = sorted(tree.nodes, key=lambda x: int(x))[-1]
    # edge: (start, end, weight)
    directed_edges = []
    traversed_nodes = []

    traversed_nodes.append(most_inner_node)
    for node in tree.adjacency_list[most_inner_node].keys(): 
        directed_edges.append((most_inner_node, node, tree.adjacency_list[most_inner_node][node]))
        traversed_nodes.append(node)
        traverse_tree_return_directed_edges(tree, node, directed_edges, traversed_nodes)
    
    # now all directed edges have been added
    directed_tree = Tree("")
    for edge in directed_edges: 
        start_node = edge[0]
        end_node = edge[1]
        weight = edge[2]
        if start_node not in directed_tree.adjacency_list.keys(): 
            directed_tree.adjacency_list[start_node] = {}
        directed_tree.adjacency_list[start_node][end_node] = weight
    for leaf in tree.leaves: 
        directed_tree.adjacency_list[leaf] = {}
    
    return directed_tree


def traverse_tree_return_directed_edges(tree: Tree, node, directed_edges: list, traversed_nodes: list):
    for next_node in [item for item in tree.adjacency_list[node].keys() if not(item in traversed_nodes)]: 
        directed_edges.append((node, next_node, tree.adjacency_list[node][next_node]))
        traversed_nodes.append(next_node)
        traverse_tree_return_directed_edges(tree, next_node, directed_edges, traversed_nodes)

def parse_undirected_tree_with_mapping(edges_str): 
    edges = edges_str.split("\n")
    tree = Tree("")
    leaves_node_count = 0
    leaves_str_list = []
    mapping = {}
    for edge in edges: 
        nodes = edge.split("->")
        node_0 = nodes[0]
        node_1 = nodes[1]

        # 如果这个node不是数字，则是需要通过mapping进行index和具体名称的映射
        if not(re.match(pattern_num, node_0)) and not(node_0 in leaves_str_list): 
            temp = node_0
            node_0 = str(leaves_node_count)
            leaves_node_count += 1
            mapping[node_0] = temp
            leaves_str_list.append(temp)
        else:
            for node_str in mapping: 
                if mapping[node_str] == node_0:
                    node_0 = node_str
        if not(re.match(pattern_num, node_1)) and not(node_1 in leaves_str_list): 
            temp = node_1
            node_1 = str(leaves_node_count)
            leaves_node_count += 1
            mapping[node_1] = temp
            leaves_str_list.append(temp)
        else:
            for node_str in mapping: 
                if mapping[node_str] == node_1: 
                    node_1 = node_str

        if not(node_0 in tree.adjacency_list.keys()): 
            tree.adjacency_list[node_0] = {}
        if not(node_1 in tree.adjacency_list.keys()): 
            tree.adjacency_list[node_1] = {}
        
        # add unweighed, directed edge
        tree.adjacency_list[node_0][node_1] = 0

    tree.refresh_tree()
    
    return tree, mapping


# SmallParsimony(T, Character)
#     for each node v in tree T
#         Tag(v) ← 0
#         if v is a leaf
#             Tag(v) ← 1
#             for each symbol k in the alphabet
#                 if Character(v) = k
#                     sk(v) ← 0
#                 else
#                     sk(v) ← ∞
#     while there exist ripe nodes in T
#         v ← a ripe node in T
#         Tag(v) ← 1
#         for each symbol k in the alphabet
#             sk(v) ← minimumall symbols i {si(Daughter(v))+αi,k} + minimumall symbols j {sj(Son(v))+αj,k}
#     return minimum over all symbols k {sk(v)}
def small_parsimony_labeling(tree, characters, mapping): 
    m = len(mapping[list(mapping.keys())[0]])
    leaves = set(mapping.keys())
    min_parsimony_score = 0
    for index in range(m):
        tag = {}
        dp_table = {}
        # backtrack table {node: {char: [left_char, right_char]}}
        backtrack_table = {}
        for node in tree.adjacency_list.keys(): 
            tag[node] = 0
            dp_table[node] = {}
            # if node is a leaf (otherwise it will not appear in the keys of dict mapping)
            if node in leaves: 
                tag[node] = 1

                for char in characters: 
                    if mapping[node][index] == char: 
                        dp_table[node][char] = 0
                    else:
                        dp_table[node][char] = float("inf")
        
        # while there still exist ripe nodes
        last_ripe_node = ""
        while True: 
            # We call an internal node of T ripe if 
            # its tag is 0 but its children’s tags are both 1
            ripe_nodes = get_ripe_nodes(tree, tag)
            if len(ripe_nodes) == 0: 
                temp_min_score = float("inf")
                for char in characters: 
                    temp = dp_table[last_ripe_node][char]
                    if temp < temp_min_score: 
                        temp_min_score = temp
                
                min_parsimony_score += temp_min_score
                break

            # choose a ripe node
            ripe_node = ripe_nodes.pop()
            tag[ripe_node] = 1
            backtrack_table[ripe_node] = {}
            for char in characters: 
                child_nodes = list(tree.adjacency_list[ripe_node].keys())
                left_node = child_nodes[0]
                right_node = child_nodes[1]
                left_min = get_min_extending_val_and_char(dp_table, characters, left_node, char)
                right_min = get_min_extending_val_and_char(dp_table, characters, right_node, char)
                dp_table[ripe_node][char] = left_min[0] + right_min[0]
                backtrack_table[ripe_node][char] = [left_min[1], right_min[1]]

            # the last ripe node, i.e. root
            last_ripe_node = ripe_node
        
        # perform backtrack using values recorded in backtrack_table
        temp_mapping = {}
        last_ripe_node_char = min(characters, key=dp_table[last_ripe_node].get)
        temp_mapping[last_ripe_node] = last_ripe_node_char
        do_backtrack(tree, last_ripe_node, last_ripe_node_char, dp_table, backtrack_table, characters, temp_mapping)
        if index == 0: 
            for node in [item for item in temp_mapping.keys() if item not in leaves]:
                mapping[node] = temp_mapping[node]
        else:
            for node in [item for item in temp_mapping.keys() if item not in leaves]:
                mapping[node] += temp_mapping[node]
        
        # print("###### index: %d ######" %index)
        # for node in [item for item in temp_mapping.keys() if item not in leaves]: 
        #     print("%s = %s, val(%d %d %d %d)" %(node, temp_mapping[node], 
        #                                         dp_table[node]["A"], 
        #                                         dp_table[node]["T"], 
        #                                         dp_table[node]["C"], 
        #                                         dp_table[node]["G"]))
        # print()

    # for node in mapping: 
    #     print("%s: %s" %(node, mapping[node]))
    
    # construct directed tree to a complete weighed tree
    for node in tree.adjacency_list.keys(): 
        for child in tree.adjacency_list[node].keys(): 
            s1 = mapping[node]
            s2 = mapping[child]
            weight = universal_utils.get_hamming_distance(s1, s2)

            tree.adjacency_list[node][child] = weight
            tree.adjacency_list[child][node] = weight
    
    tree.refresh_tree()

    return min_parsimony_score, tree, mapping

# backtrack function recursive
def do_backtrack(tree, current_node, current_node_char, dp_table, backtrack_table, characters, temp_mapping): 
    left = list(tree.adjacency_list[current_node].keys())[0]
    right = list(tree.adjacency_list[current_node].keys())[1]

    left_char = backtrack_table[current_node][current_node_char][0]
    right_char = backtrack_table[current_node][current_node_char][1]
    temp_mapping[left] = left_char
    temp_mapping[right] = right_char
    if len(tree.adjacency_list[left].keys()) != 0: 
        do_backtrack(tree, left, left_char, dp_table, backtrack_table, characters, temp_mapping)
    if len(tree.adjacency_list[right].keys()) != 0: 
        do_backtrack(tree, right, right_char, dp_table, backtrack_table, characters, temp_mapping)

# We call an internal node of T ripe if 
# its tag is 0 but its children’s tags are both 1
def get_ripe_nodes(tree, tag): 
    potential_ripe_nodes = []
    for node in tag.keys(): 
        if tag[node] == 0: 
            potential_ripe_nodes.append(node)
    ripe_nodes = potential_ripe_nodes[:]
    for node in potential_ripe_nodes: 
        child_nodes = tree.adjacency_list[node].keys()
        for child in child_nodes: 
            if tag[child] == 0:
                ripe_nodes.remove(node)
                break
    
    return ripe_nodes

# according to current symbol being traversed, 
# as well as value recorded in the dp_table dictionary, 
# return the min value of extending from the child to itself
# if different chars can be taken from characters list
def get_min_extending_val_and_char(dp_table, characters, node, current_symbol): 
    min_extending_val = float("inf")
    min_extending_char = ""
    for char in characters: 
        temp = dp_table[node][char]
        if char != current_symbol: 
            temp += 1
        
        if temp < min_extending_val: 
            min_extending_val = temp
            min_extending_char = char
        elif temp == min_extending_val and char == current_symbol:
            min_extending_char = char
    
    return min_extending_val, min_extending_char
            
# wrapping function for generate rooted tree using small parsimony algorithm
def small_parsimony_generate_rooted_tree(edges_str): 
    parsed_edges = parse_directed_tree_with_mapping(edges_str)
    unweighed_tree = parsed_edges[0]
    mapping = parsed_edges[1]
    min_parsimony_score, tree, mapping = small_parsimony_labeling(unweighed_tree, DNA_RNA_utils.DNA_BASES, mapping)

    to_print = ""
    to_print += str(min_parsimony_score) + "\n"
    to_print += Tree.show_tree_in_adjacent_list(tree, mapping=mapping)
    pyperclip.copy(to_print)

    return min_parsimony_score, tree, mapping

# wrapping function for generate unrooted tree using small parsimony algorithm
def small_parsimony_generate_unrooted_tree(edges_str): 
    # When the position of the root in the tree is unknown, 
    # we can simply assign the root to any edge that we like, 
    # apply SmallParsimony to the resulting rooted tree, and then remove the root. 
    
    # edges_str should be processed as an directed tree (from root to leaves)
    # edges = edges_str.split("\n")
    # edges = [edges[i] for i in range(len(edges)) if i % 2 == 1]
    # edges_str = "\n".join(edges)
    # parsed_edges = parse_directed_tree_with_mapping(edges_str)
    # unweighed_tree = parsed_edges[0]
    # mapping = parsed_edges[1]
    parsed = parse_undirected_tree_with_mapping(edges_str)
    tree = parsed[0]
    mapping = parsed[1]
    unweighed_tree = transfer_tree_to_directed(tree)

    # select an edge, and add the root node to it
    start_node = sorted(list(unweighed_tree.adjacency_list.keys()), key=lambda x: int(x))[-1]
    end_node = list(unweighed_tree.adjacency_list[start_node].keys())[0]
    new_node = str(int(start_node) + 1)
    del unweighed_tree.adjacency_list[start_node][end_node]
    unweighed_tree.adjacency_list[new_node] = {}
    unweighed_tree.adjacency_list[new_node][start_node] = 0
    unweighed_tree.adjacency_list[new_node][end_node] = 0

    # perform the small parsimony algorithm
    min_parsimony_score, tree, mapping = small_parsimony_labeling(unweighed_tree, DNA_RNA_utils.DNA_BASES, mapping)

    # delete the root node added before
    # and restore the splitted edge
    weight_new_start = tree.adjacency_list[new_node][start_node]
    weight_new_end = tree.adjacency_list[new_node][end_node]
    weight_start_end = weight_new_start + weight_new_end
    tree.adjacency_list[start_node][end_node] = weight_start_end
    tree.adjacency_list[end_node][start_node] = weight_start_end
    del tree.adjacency_list[new_node]
    del tree.adjacency_list[start_node][new_node]
    del tree.adjacency_list[end_node][new_node]

    tree.refresh_tree()

    to_print = ""
    to_print += str(min_parsimony_score) + "\n"
    to_print += Tree.show_tree_in_adjacent_list(tree, mapping=mapping)
    pyperclip.copy(to_print)

    return min_parsimony_score, tree, mapping

# do small parsimony to given tree
def small_parsimony_generate_unrooted_tree_by_tree(tree: Tree, mapping): 
    # delete all key-val pairs in mapping other than leaves nodes
    for node in [item for item in mapping.keys() if not(item) in tree.leaves]: 
        del mapping[node]

    unweighed_tree = transfer_tree_to_directed(tree)

    # select an edge, and add the root node to it
    start_node = sorted(list(unweighed_tree.adjacency_list.keys()), key=lambda x: int(x))[-1]
    end_node = list(unweighed_tree.adjacency_list[start_node].keys())[0]
    new_node = str(int(start_node) + 1)
    del unweighed_tree.adjacency_list[start_node][end_node]
    unweighed_tree.adjacency_list[new_node] = {}
    unweighed_tree.adjacency_list[new_node][start_node] = 0
    unweighed_tree.adjacency_list[new_node][end_node] = 0

    # perform the small parsimony algorithm
    min_parsimony_score, tree, mapping = small_parsimony_labeling(unweighed_tree, DNA_RNA_utils.DNA_BASES, mapping)

    # delete the root node added before
    # and restore the splitted edge
    weight_new_start = tree.adjacency_list[new_node][start_node]
    weight_new_end = tree.adjacency_list[new_node][end_node]
    weight_start_end = weight_new_start + weight_new_end
    tree.adjacency_list[start_node][end_node] = weight_start_end
    tree.adjacency_list[end_node][start_node] = weight_start_end
    del tree.adjacency_list[new_node]
    del tree.adjacency_list[start_node][new_node]
    del tree.adjacency_list[end_node][new_node]

    tree.refresh_tree()

    return min_parsimony_score, tree, mapping

# Input: Two internal nodes a and b specifying an edge e, followed by an adjacency list of an unrooted binary tree.
# Output: Two adjacency lists representing the nearest neighbors of the tree with respect to e. Separate the adjacency lists with a blank line.
def get_nearest_neighbors(edge : tuple[str, str], tree : Tree):
    node_0 = edge[0]
    node_1 = edge[1]

    # select 4 subtrees, 2 from node_0, 2 from node_1, because tree should be binary tree
    node_0_0 = sorted([node for node in list(tree.adjacency_list[node_0].keys()) if node != node_1], key=lambda x: int(x))[0]
    node_0_1 = sorted([node for node in list(tree.adjacency_list[node_0].keys()) if node != node_1], key=lambda x: int(x))[1]
    node_1_0 = sorted([node for node in list(tree.adjacency_list[node_1].keys()) if node != node_0], key=lambda x: int(x))[0]
    node_1_1 = sorted([node for node in list(tree.adjacency_list[node_1].keys()) if node != node_0], key=lambda x: int(x))[1]

    # the first neighbor
    neighbor_1 = Tree.copy_tree(tree)
    weight_0_01 = neighbor_1.adjacency_list[node_0][node_0_1]
    weight_1_10 = neighbor_1.adjacency_list[node_1][node_1_0]

    del neighbor_1.adjacency_list[node_0][node_0_1]
    del neighbor_1.adjacency_list[node_0_1][node_0]
    del neighbor_1.adjacency_list[node_1][node_1_0]
    del neighbor_1.adjacency_list[node_1_0][node_1]

    neighbor_1.adjacency_list[node_0][node_1_0] = weight_1_10
    neighbor_1.adjacency_list[node_1_0][node_0] = weight_1_10
    neighbor_1.adjacency_list[node_1][node_0_1] = weight_0_01
    neighbor_1.adjacency_list[node_0_1][node_1] = weight_0_01

    neighbor_1.refresh_tree()

    # the second neighbor
    neighbor_2 = Tree.copy_tree(tree)
    weight_0_01 = neighbor_2.adjacency_list[node_0][node_0_1]
    weight_1_11 = neighbor_2.adjacency_list[node_1][node_1_1]

    del neighbor_2.adjacency_list[node_0][node_0_1]
    del neighbor_2.adjacency_list[node_0_1][node_0]
    del neighbor_2.adjacency_list[node_1][node_1_1]
    del neighbor_2.adjacency_list[node_1_1][node_1]

    neighbor_2.adjacency_list[node_0][node_1_1] = weight_1_11
    neighbor_2.adjacency_list[node_1_1][node_0] = weight_1_11
    neighbor_2.adjacency_list[node_1][node_0_1] = weight_0_01
    neighbor_2.adjacency_list[node_0_1][node_1] = weight_0_01

    neighbor_2.refresh_tree()

    neighbor_1_str = neighbor_1.show_unweighed_tree(if_show=False)
    # print()
    neighbor_2_str = neighbor_2.show_unweighed_tree(if_show=False)
    to_copy = neighbor_1_str + "\n\n" + neighbor_2_str
    pyperclip.copy(to_copy)

    return (neighbor_1, neighbor_2)

# NearestNeighborInterchange(Strings)
#      score ← ∞
#      generate an arbitrary unrooted binary tree Tree with |Strings| leaves
#      label the leaves of Tree by arbitrary strings from Strings
#      solve  the  Small Parsimony in an Unrooted Tree Problem for Tree
#      label the internal nodes of Tree according to a most parsimonious labeling
#      newScore ← the parsimony score of Tree
#      newTree ← Tree
#      while newScore < score
#           score ← newScore
#           Tree ← newTree
#           for each internal edge e in Tree
#                for each nearest neighbor NeighborTree of Tree with respect to the edge e
#                     solve the Small Parsimony in an Unrooted Tree Problem for NeighborTree
#                     neighborScore ← the minimum parsimony score of NeighborTree
#                     if neighborScore < newScore
#                          newScore ← neighborScore
#                          newTree ← NeighborTree
#      return newTree
def nearest_neighbor_interchange_algorithm(edges_str): 
    score = float("inf")
    parsed = parse_undirected_tree_with_mapping(edges_str)
    tree = parsed[0]
    mapping = parsed[1]
    new_score, tree, mapping = small_parsimony_generate_unrooted_tree_by_tree(tree, mapping)
    to_print = ""
    new_tree = tree
    step_count = 0
    while new_score < score: 
        score = new_score
        tree = copy.deepcopy(new_tree)

        internal_edges = get_internal_edges(tree)
        for internal_edge in internal_edges:
            neighbors = get_nearest_neighbors(internal_edge, tree)
            # i.e. change the shape of the previous tree
            for neighbor in neighbors: 
                score_temp, tree_temp, mapping_temp = small_parsimony_generate_unrooted_tree_by_tree(neighbor, mapping)
                if score_temp < new_score: 
                    new_score = score_temp
                    new_tree = tree_temp
                    mapping = mapping_temp
                    print(new_score)
                    tree_print_temp = new_tree.show_tree_in_adjacent_list(mapping=mapping)
                    print()
                    to_print += str(new_score) + "\n" + tree_print_temp + "\n\n"
                    step_count += 1
    
    to_print.strip()
    pyperclip.copy(to_print)
    print("step count: %d" %step_count)
    return step_count, to_print




def get_internal_edges(tree: Tree): 
    internal_nodes = [node for node in tree.nodes if not(node in tree.leaves)]
    internal_edges = []
    traversed_nodes = []
    for node in internal_nodes: 
        traversed_nodes.append(node)
        next_nodes = [item for item in tree.adjacency_list[node].keys() \
                          if item in internal_nodes and not(item in traversed_nodes)]
        for next_node in next_nodes: 
            internal_edges.append((node, next_node))
    
    return internal_edges

if __name__ == "__main__": 
    # 1.2 step 9
#     leaves_node_num = 4
#     edges_str = """4->CAAATCCC
# 4->ATTGCGAC
# 5->CTGCGCTG
# 5->ATGGACGA
# 6->4
# 6->5"""
#     unweighed_tree = parse_directed_tree_with_mapping(edges_str)[0]
#     mapping = parse_directed_tree_with_mapping(edges_str)[1]
#     small_parsimony_generate_rooted_tree(leaves_node_num, edges_str)

    # extra dataset
    # with open("./datasets/Small_Parsimony.txt", "r") as f: 
    #     # input
    #     f.readline()
    #     leaves_node_num = int(f.readline().strip())
    #     edges_str = ""
    #     temp = ""
    #     while temp != "Output\n": 
    #         edges_str += temp
    #         temp = f.readline()
    #     edges_str = edges_str.strip()
    #     unweighed_tree = parse_directed_tree_with_mapping(edges_str)[0]
    #     mapping = parse_directed_tree_with_mapping(edges_str)[1]
    #     min_parsimony_score, tree, mapping = small_parsimony_labeling(unweighed_tree, DNA_RNA_utils.DNA_BASES, mapping)
    #     Tree.show_tree_in_adjacent_list(unweighed_tree, mapping=mapping)
    #     print(min_parsimony_score)
        
    # with open("./datasets/dataset_30291_9.txt", "r") as f: 
    #     leaves_node_num = int(f.readline().strip())
    #     edges_str = f.read().strip()
    #     small_parsimony_generate_rooted_tree(edges_str)

    # 1.2 step 12
#     leaves_node_num = 4
#     edges_str = """TCGGCCAA->4
# 4->TCGGCCAA
# CCTGGCTG->4
# 4->CCTGGCTG
# CACAGGAT->5
# 5->CACAGGAT
# TGAGTACC->5
# 5->TGAGTACC
# 4->5
# 5->4"""
#     small_parsimony_generate_unrooted_tree(edges_str)

    # with open("./datasets/dataset_30291_11.txt", "r") as f: 
    #     leaves_node_num = int(f.readline().strip())
    #     edges_str = f.read().strip()
    #     small_parsimony_generate_unrooted_tree(edges_str)

    # 1.3 step 6
#     edge = ("5", "4")
#     edges_str = """0->4
# 4->0
# 1->4
# 4->1
# 2->5
# 5->2
# 3->5
# 5->3
# 4->5
# 5->4"""
#     tree = Tree(edges_str)
#     neighbors = get_nearest_neighbors(("5", "4"), tree)

    # with open("./datasets/dataset_30292_6.txt", "r") as f: 
    #     edge = f.readline().strip().split(" ")
    #     edges_str = f.read().strip()
    #     tree = Tree(edges_str)
    #     get_nearest_neighbors(edge, tree)


    # 1.3 step 9
#     leaves_node_num = 5
#     edges_str = """GCAGGGTA->5
# TTTACGCG->5
# CGACCTGA->6
# GATTCCAC->6
# 5->TTTACGCG
# 5->GCAGGGTA
# 5->7
# TCCGTAGT->7
# 7->5
# 7->6
# 7->TCCGTAGT
# 6->GATTCCAC
# 6->CGACCTGA
# 6->7"""
#     nearest_neighbor_interchange_algorithm(edges_str)

    # extra dataset
#     leaves_node_num = 8
#     edges_str = """AGCTCAGCGCCCCGGAGCACCCTCCTGAAGTATGCACATT->11
# CCGCCGTACACACAGTCTTGAACCATTTACCGCAGTTTCC->10
# GTCCAAGAGTATGTGAAACCTGCAGTGACGAAGGCGAGAT->8
# ACTGAGGTACGGGTTATACCGCGCATCTGCGAGTAAAACA->11
# 10->CCGCCGTACACACAGTCTTGAACCATTTACCGCAGTTTCC
# 10->CGGTCCTCTAGGAGCTTGTCTTTATCTGCCGCCGACATGC
# 10->12
# ATTATGGGGCACTGAGCATACGCAAACGACTATGCTTTCC->9
# 9->ATTATGGGGCACTGAGCATACGCAAACGACTATGCTTTCC
# 9->TGCGGCGGGCGCGCCCAAACAGCGTGACCAAGTCGATGCA
# 9->13
# 8->GTCCAAGAGTATGTGAAACCTGCAGTGACGAAGGCGAGAT
# 8->CCACACGTGGCTGTTATATGATATTAATATATTTAATCTT
# 8->13
# 11->AGCTCAGCGCCCCGGAGCACCCTCCTGAAGTATGCACATT
# 11->ACTGAGGTACGGGTTATACCGCGCATCTGCGAGTAAAACA
# 11->12
# CCACACGTGGCTGTTATATGATATTAATATATTTAATCTT->8
# 13->8
# 13->9
# 13->12
# TGCGGCGGGCGCGCCCAAACAGCGTGACCAAGTCGATGCA->9
# 12->11
# 12->10
# 12->13
# CGGTCCTCTAGGAGCTTGTCTTTATCTGCCGCCGACATGC->10"""
#     nearest_neighbor_interchange_algorithm(edges_str)

    with open("./datasets/dataset_30292_8.txt", "r") as f: 
        leaves_node_num = int(f.readline().strip())
        edges_str = f.read().strip()
        nearest_neighbor_interchange_algorithm(edges_str)

    # week 3 Quiz
    # no coding quiz though