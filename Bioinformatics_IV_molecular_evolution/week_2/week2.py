from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
import pyperclip
import re
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_1.week1 import *
import copy

# Discrepancy(T, D) = Σ(di,j(T) - Di,j) ^ 2, and 0 <= j <= i <= n - 1
def discrepancy_between_tree_and_distance_matrix(tree, distance_matrix): 
    additive_distance_matrix = tree.leaves_distance_matrix[:]
    
    discrepancy = 0
    leaves_node_num = len(distance_matrix)
    # print(leaves_node_num)
    for i in range(leaves_node_num): 
        for j in range(i + 1, leaves_node_num): 
            temp = (additive_distance_matrix[i][j] - distance_matrix[i][j]) ** 2
            discrepancy += temp
    
    return discrepancy

# UPGMA(D, n) 
#     Clusters ← n single-element clusters labeled 1, ... , n
#     construct a graph T with n isolated nodes labeled by single elements 1, ... , n
#     for every node v in T 
#         Age(v) ← 0
#     while there is more than one cluster
#         find the two closest clusters Ci and Cj 
#         merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
#         add a new node labeled by cluster Cnew to T
#         connect node Cnew to Ci and Cj by directed edges 
#         Age(Cnew) ← DCi, Cj / 2
#         remove the rows and columns of D corresponding to Ci and Cj
#         remove Ci and Cj from Clusters
#         add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
#         add Cnew to Clusters
#     root ← the node in T corresponding to the remaining cluster
#     for each edge (v, w) in T
#         length of (v, w) ← Age(v) - Age(w)
#     return T
def UPGMA_generate_tree(leaves_distance_matrix, total_node_num): 
    distance_matrix = copy.deepcopy(leaves_distance_matrix)
    clusters = {}
    tree = Tree("")
    ages = {}
    for i in range(total_node_num): 
        current_node = str(i)
        clusters[current_node] = [current_node]
        tree.adjacency_list[current_node] = {}
        ages[current_node] = 0
    
    current_node_num = total_node_num
    while len(clusters) != 1: 
        # find the two closest clusters Ci and Cj
        min_distance = float("inf")
        c_i = ""
        c_j = ""
        # if we can itinerate the distance matrix with the items in the set of cluster roots
        # then we don't have to remove columns and rows of Ci and Cj in the distance matrix
        for node_i in clusters.keys(): 
            for node_j in clusters.keys():
                if node_j == node_i:
                    continue
                i_index = int(node_i)
                j_index = int(node_j)
                distance = distance_matrix[i_index][j_index]
                if distance < min_distance:
                    c_i = node_i
                    c_j = node_j
                    min_distance = distance
        
        c_new = str(current_node_num)
        ages[c_new] = distance_matrix[int(c_i)][int(c_j)] / 2
        current_node_num += 1

        # merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements
        # add Cnew to Clusters
        clusters[c_new] = []
        # a cluster does not involve internal nodes
        for node in clusters[c_i] + clusters[c_j]:
            clusters[c_new].append(node)
        len_c_i = len(clusters[c_i])
        len_c_j = len(clusters[c_j])
        del clusters[c_i]
        del clusters[c_j]

        # add a new node labeled by cluster Cnew to T
        # connect node Cnew to Ci and Cj by directed edges 
        tree.adjacency_list[c_new] = {}
        # at the same time, update the edge length using the ages array
        tree.adjacency_list[c_new][c_i] = ages[c_new] - ages[c_i]
        tree.adjacency_list[c_i][c_new] = ages[c_new] - ages[c_i] 
        tree.adjacency_list[c_new][c_j] = ages[c_new] - ages[c_j]
        tree.adjacency_list[c_j][c_new] = ages[c_new] - ages[c_j]

        # add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        distance_matrix.append([])
        for i in range(current_node_num): 
            temp = (distance_matrix[int(c_i)][i] * len_c_i + distance_matrix[int(c_j)][i] * len_c_j) / (len_c_i + len_c_j)
            if i == int(c_new):
                temp = 0
                distance_matrix[i].append(temp)
            else:
                distance_matrix[i].append(temp)
                distance_matrix[int(c_new)].append(temp)

    tree.refresh_tree()
    return tree

def get_total_distance_arr(distance_matrix): 
    n = len(distance_matrix)
    total_distance_arr = [0 for _ in range(n)]
    for i in range(n): 
        total_distance_arr[i] = sum(distance_matrix[i])
    
    return total_distance_arr

def get_neighbor_joining_matrix(distance_matrix):
    n = len(distance_matrix)
    neighbor_joining_matrix = [
        [0 for _ in range(n)] for _ in range(n)
    ]

    # define TotalDistanceD(i) as the sum ∑1≤k≤n Di,k of distances from leaf i to all other leaves
    total_distance_arr = get_total_distance_arr(distance_matrix)

    # D*i,i = 0
    # D*i,j = (n - 2) · Di,j - TotalDistanceD(i) - TotalDistanceD(j)
    for i in range(n): 
        for j in range(n): 
            if i == j: 
                neighbor_joining_matrix[i][j] = 0
            else:
                neighbor_joining_matrix[i][j] = \
                (n - 2) * distance_matrix[i][j] - total_distance_arr[i] - total_distance_arr[j]
    
    return neighbor_joining_matrix

def find_min_in_matrix(matrix): 
    min_val = float("inf")
    min_i, min_j = 0, 0
    n = len(matrix)
    for i in range(n): 
        for j in range(n): 
            if i == j: 
                continue
            val = matrix[i][j]
            if val < min_val: 
                min_val = val
                min_i = i
                min_j = j
    
    return min_i, min_j

# a wrapper function to call the recursion neighbor joining algorithm
def neighbor_joining_generate_tree(leaves_distance_matrix): 
    total_node_num = len(leaves_distance_matrix)
    mapping = None
    tree = neighbor_joining_algorithm(leaves_distance_matrix, total_node_num, mapping)
    tree.refresh_tree()
    return tree

# NeighborJoining(D)
#     n ← number of rows in D
#     if n = 2
#         T ← tree consisting of a single edge of length D1,2
#         return T
#     D* ← neighbor-joining matrix constructed from the distance matrix D
#     find elements i and j such that D*i,j is a minimum non-diagonal element of D*
#     Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
#     limbLengthi ← (1/2)(Di,j + Δ)
#     limbLengthj ← (1/2)(Di,j - Δ)
#     add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
#     D ← D with rows i and j removed
#     D ← D with columns i and j removed
#     T ← NeighborJoining(D)
#     add two new limbs (connecting node m with leaves i and j) to the tree T
#     assign length limbLengthi to Limb(i)
#     assign length limbLengthj to Limb(j)
#     return T
def neighbor_joining_algorithm(leaves_distance_matrix, total_node_num, mapping): 
    distance_matrix = copy.deepcopy(leaves_distance_matrix)
    n = len(distance_matrix)
    # recursion termination
    if n == 2: 
        weight = leaves_distance_matrix[0][1]
        node_0 = list(mapping.keys())[0]
        node_1 = list(mapping.keys())[1]
        # adjacency_list_str = "%s->%s:%d\n%s->%s:%d" %(
        #     node_0, node_1, weight, node_1, node_0, weight
        # )
        # tree = Tree(adjacency_list_str)
        tree = Tree("")
        tree.adjacency_list[node_0] = {node_1: weight}
        tree.adjacency_list[node_1] = {node_0: weight}

        return tree

    # record the relationship between node and index in the matrix
    if mapping == None: 
        # 记录字符串形式的node名称和distance matrix中角标的映射关系
        mapping = {}
        for k in range(n): 
            current_node = str(k)
            mapping[current_node] = k

    total_distance_arr = get_total_distance_arr(distance_matrix)
    neighbor_joining_matrix = get_neighbor_joining_matrix(distance_matrix)
    # find elements i and j such that D*i,j is a minimum non-diagonal element of D*
    min_i_j = find_min_in_matrix(neighbor_joining_matrix)
    index_i = min_i_j[0]
    index_j = min_i_j[1]
    node_i = universal_utils.get_key_by_distincet_val(mapping, index_i)
    node_j = universal_utils.get_key_by_distincet_val(mapping, index_j)

    delta_distance = (total_distance_arr[index_i] - total_distance_arr[index_j]) / (n - 2)
    limb_i = (distance_matrix[index_i][index_j] + delta_distance) / 2
    limb_j = (distance_matrix[index_i][index_j] - delta_distance) / 2

    # for line in neighbor_joining_matrix: 
    #     universal_utils.print_arr(line)

    node_m = str(total_node_num)
    total_node_num += 1
    distance_matrix.append([0 for _ in range(n + 1)])
    for k in range(n): 
        temp = (distance_matrix[k][index_i] + distance_matrix[k][index_j] - distance_matrix[index_i][index_j]) / 2
        distance_matrix[k].append(temp)
        distance_matrix[n][k] = temp
    
    # remove cols and rows with index_i and index_j
    if index_i > index_j: 
        index_i = index_j
    distance_matrix.pop(index_j)
    distance_matrix.pop(index_i)
    for k in range(len(distance_matrix)): 
        distance_matrix[k].pop(index_j)
        distance_matrix[k].pop(index_i)

    # udpate mapping dictionary!!
    del mapping[node_i]
    del mapping[node_j]
    for k in range(n): 
        node = universal_utils.get_key_by_distincet_val(mapping, k)
        if node: 
            if k > index_i and k < index_j: 
                mapping[node] -= 1
            elif k > index_j: 
                mapping[node] -= 2

    mapping[node_m] = len(distance_matrix) - 1

    # for line in distance_matrix: 
    #     universal_utils.print_arr(line)
    # print(mapping)

    tree = neighbor_joining_algorithm(distance_matrix, total_node_num, mapping)
    tree.adjacency_list[node_m][node_i] = limb_i
    tree.adjacency_list[node_m][node_j] = limb_j
    tree.adjacency_list[node_i] = {}
    tree.adjacency_list[node_j] = {}
    tree.adjacency_list[node_i][node_m] = limb_i
    tree.adjacency_list[node_j][node_m] = limb_j

    return tree

if __name__ == "__main__": 
    # 1.1 step 2
#     adjacency_list_str = """0->5:3
# 1->5:4
# 2->4:1
# 3->4:2
# 5->0:3
# 5->1:4
# 5->4:5
# 4->2:1
# 4->3:2
# 4->5:5"""
#     distance_matrix_str = """0\t3\t4\t3\n3\t0\t4\t5\n4\t4\t0\t2\n3\t5\t2\t0"""

#     tree = Tree(adjacency_list_str)
#     distance_matrix = universal_utils.parse_arr_2_dimension(distance_matrix_str)

#     universal_utils.print_arr(tree.leaves_distance_matrix)
#     universal_utils.print_arr(distance_matrix)
    
#     print("discrepancy between tree and non-additive distance matrix: %d" %discrepancy_between_tree_and_distance_matrix(tree, distance_matrix))

    # 1.2 step 8
#     total_node_num = 4
#     distance_matrix = universal_utils.parse_arr_2_dimension("""0	20	17	11
# 20	0	20	13
# 17	20	0	10
# 11	13	10	0""")
#     tree = UPGMA_generate_tree(distance_matrix, total_node_num)
#     Tree.show_tree_in_adjacent_list(tree)

        #        6
        #      /   \
        #     5     1
        #   /   \     
        # 0       4
        #       /   \
        #     2       3

    # extra dataset
    # with open("./UPGMA.txt", "r") as f:
    #     f.readline()
    #     total_node_num = int(f.readline().strip())
    #     leaves_distance_matrix_str = ""
    #     temp = f.readline()
    #     while temp != "Output\n": 
    #         leaves_distance_matrix_str += temp
    #         temp = f.readline()
    #     leaves_distance_matrix_str = leaves_distance_matrix_str.strip()
    #     leaves_distance_matrix = universal_utils.parse_arr_2_dimension(leaves_distance_matrix_str)
    #     tree = UPGMA_generate_tree(leaves_distance_matrix, total_node_num)
    #     Tree.show_tree_in_adjacent_list(tree, "f")

    # with open("./dataset_30288_8.txt", "r") as f:
    #     total_node_num = int(f.readline().strip())
    #     leaves_distance_matrix = universal_utils.parse_arr_2_dimension(f.read().strip())
    #     tree = UPGMA_generate_tree(leaves_distance_matrix, total_node_num)
    #     Tree.show_tree_in_adjacent_list(tree)

    # 1.2 step 10
    # with open("./coronavirus_distance_matrix_nonadditive.txt", "r") as f:
    #     species = universal_utils.parss_arr(f.readline().strip())
    #     # universal_utils.print_arr(species)
    #     leaves_distance_matrix = []
    #     for line in f.read().strip().split("\n"): 
    #         leaves_distance_matrix.append(
    #             [int(item) for item in line.split("\t")[1: ]]
    #         )
    #     tree = UPGMA_generate_tree(leaves_distance_matrix, len(species))
    #     tree.refresh_tree()
    #     Tree.show_tree_in_adjacent_list(tree, "f")

    # 1.3 step 6
#     distance_matrix = universal_utils.parse_arr_2_dimension("""0\t13\t21\t22
# 13\t0\t12\t13
# 21\t12\t0\t13
# 22\t13\t13\t0""")
#     distance_matrix = universal_utils.parse_arr_2_dimension("""0\t3\t4\t3
# 3\t0\t4\t5
# 4\t4\t0\t2
# 3\t5\t2\t0""")
    distance_matrix = universal_utils.parse_arr_2_dimension("""0\t23\t27\t20
23\t0\t30\t28
27\t30\t0\t30
20\t28\t30\t0""")
    # neighbor_joining_matrix = get_neighbor_joining_matrix(distance_matrix)
    # universal_utils.print_arr(neighbor_joining_matrix)
    tree = neighbor_joining_generate_tree(distance_matrix)
    Tree.show_tree_in_adjacent_list(tree, "f")

    # 1.3 step 7
    # with open("./dataset_30289_7.txt", "r") as f: 
    #     total_node_num = int(f.readline().strip())
    #     distance_matrix = universal_utils.parse_arr_2_dimension(f.read().strip())
    #     tree = neighbor_joining_generate_tree(distance_matrix)
    #     Tree.show_tree_in_adjacent_list(tree, "f")

    # 1.3 step 8
    # with open("./coronavirus_distance_matrix_nonadditive.txt", "r") as f:
    #     species = universal_utils.parss_arr(f.readline().strip())
    #     # universal_utils.print_arr(species)
    #     leaves_distance_matrix = []
    #     for line in f.read().strip().split("\n"): 
    #         leaves_distance_matrix.append(
    #             [int(item) for item in line.split("\t")[1: ]]
    #         )
    #     tree = neighbor_joining_generate_tree(leaves_distance_matrix)
    #     Tree.show_tree_in_adjacent_list(tree, "f")

    # week 2 Quiz
    # problem 2
#     print("###### problem 2 ######")
# #     distance_matrix = universal_utils.parse_arr_2_dimension("0\t13\t16\t10\n13\t0\t21\t15\n16\t21\t0\t18\n10\t15\t18\t0")
# #     adjacency_list_str = """0->4:5 
# # 1->4:9
# # 2->4:12
# # 3->4:6
# # 4->0:5
# # 4->1:9
# # 4->2:12
# # 4->3:6"""
#     distance_matrix = universal_utils.parse_arr_2_dimension("0 14 17 17\n14 0 13 7\n17 13 0 20\n17 7 20 0")
#     adjacency_list_str = """0->4:9
# 4->0:9
# 2->4:8
# 4->2:8
# 1->5:2
# 5->1:2
# 3->5:6
# 5->3:6
# 4->5:4
# 5->4:4"""
#     tree = Tree(adjacency_list_str)
#     for line in tree.leaves_distance_matrix: 
#         universal_utils.print_arr(line)
#     discrepancy = discrepancy_between_tree_and_distance_matrix(tree, distance_matrix)
#     print("discrepancy between tree and matrix: %d" %discrepancy)
#     print()

#     # problem 3
#     # distance_matrix = universal_utils.parse_arr_2_dimension("0\t20\t9\t11\n20\t0\t17\t11\n9\t17\t0\t8\n11\t11\t8\t0")
#     # tree = UPGMA_generate_tree(distance_matrix, len(distance_matrix))
#     # tree.refresh_tree()
#     # Tree.show_tree_in_adjacent_list(tree, "f")

#     # problem 4
#     print("###### problem 4 ######")
#     distance_matrix = universal_utils.parse_arr_2_dimension("0 14 17 17\n14 0 7 13\n17 7 0 16\n17 13 16 0")
#     neighbor_joining_matrix = get_neighbor_joining_matrix(distance_matrix)
#     print("D*[2][3] = %d" %neighbor_joining_matrix[2][3])
#     print()

#     # problem 5
#     print("###### problem 5 ######")
#     distance_matrix = universal_utils.parse_arr_2_dimension("0\t20\t9\t11\n20\t0\t17\t11\n9\t17\t0\t8\n11\t11\t8\t0")
#     tree = neighbor_joining_generate_tree(distance_matrix)
#     Tree.show_tree_in_adjacent_list(tree, "f")
    
