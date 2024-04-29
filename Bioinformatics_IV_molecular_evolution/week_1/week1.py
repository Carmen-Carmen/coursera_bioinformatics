from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
import pyperclip
import re
import copy

class Graph:
    def __init__(self, adjacency_list_str):
        self.adjacency_list = Graph.parse_adjacency_list(adjacency_list_str)
        self.in_and_out_degrees = Graph.get_in_and_out_degrees(self.adjacency_list)
        self.nodes = list(self.in_and_out_degrees.keys())
        self.distance_matrix = Graph.get_distance_matrix(self.adjacency_list, self.nodes)

    def parse_adjacency_list(str):
        adjacency_list = {}
        edges = str.split("\n")
        for edge in edges:
            if ":" in edge: 
                weight = int(edge.split(":")[1])
                from_node = edge.split(":")[0].split("->")[0]
                to_node = edge.split(":")[0].split("->")[1]

                if not(from_node in adjacency_list):
                    adjacency_list[from_node] = {}
                adjacency_list[from_node][to_node] = weight
            else:
                from_node = edge.split("->")[0]
                to_node = edge.split("->")[1]
        
                if not(from_node in adjacency_list):
                    adjacency_list[from_node] = {}
                adjacency_list[from_node][to_node] = 0

        return adjacency_list
    
    def get_in_and_out_degrees(adjacency_list):
        in_and_out_degrees = {}
        # get all nodes
        nodes = set()
        for from_node in adjacency_list.keys():
            nodes.add(from_node)
            for to_node in adjacency_list[from_node].keys():
                nodes.add(to_node)
        
        # initialize the in_and_out_degrees dict
        for node in nodes:
            in_and_out_degrees[node] = {
                "in": 0, 
                "out": 0, 
            }

        # calculate in and out degrees
        for from_node in adjacency_list.keys():
            in_and_out_degrees[from_node]["out"] = len(adjacency_list[from_node].keys())
            for to_node in adjacency_list[from_node].keys():
                in_and_out_degrees[to_node]["in"] += 1
        
        return in_and_out_degrees

    # get the distance between all nodes in the graph
    # implementing the Floyd-Marshall algorithm
    def get_distance_matrix(adjacency_list, nodes):
        n = len(nodes)
        distance_matrix = [
            [float("inf")] * n for _ in range(n)
        ]

        # initialize the distance matrix with direct edge weights
        for i in range(n):
            # the shortest distance from one node to itself is 0
            distance_matrix[i][i] = 0
            for j in range(n):
                node_i = str(i)
                node_j = str(j)
                if i != j and \
                node_j in adjacency_list[node_i].keys(): 
                # which means there is an edge directly connecting node_i and node_j
                    distance_matrix[i][j] = adjacency_list[node_i][node_j]
        
        # update the distance matrix using intermediate nodes
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if distance_matrix[i][k] != float("inf") and distance_matrix[k][j] != float("inf"):
                        distance_matrix[i][j] = min(
                            distance_matrix[i][j], distance_matrix[i][k] + distance_matrix[k][j]
                        )

        return distance_matrix
        
    def show_distance_matrix(self):
        to_show = ""
        for row in self.distance_matrix:
            to_show += universal_utils.print_arr(row)
        
        to_show = to_show.strip()
        pyperclip.copy(to_show)
        return to_show

    def refresh_graph(self):
        self.in_and_out_degrees = Graph.get_in_and_out_degrees(self.adjacency_list)
        self.nodes = list(self.in_and_out_degrees.keys())
        self.distance_matrix = Graph.get_distance_matrix(self.adjacency_list, self.nodes)

class Tree(Graph):
    def __init__(self, adjacency_list_str):
        if adjacency_list_str != "": 
            # call the constructor of the parent calss
            super().__init__(adjacency_list_str)
            self.leaves = Tree.get_leaves(self.in_and_out_degrees)
            self.leaves_distance_matrix = Tree.get_distance_matrix_between_leaves(self.distance_matrix, self.leaves)
        else:
            self.adjacency_list = {}
    
    def get_leaves(in_and_out_degrees):
        leaves = []
        for node in in_and_out_degrees.keys():
            if in_and_out_degrees[node]["in"] == 1 and \
                in_and_out_degrees[node]["out"] == 1:
                leaves.append(node)
        
        leaves = sorted(leaves)

        return leaves
    
    def get_distance_matrix_between_leaves(distance_matrix, leaves):
        leaves_distance_matrix = []
        # remove the rows and cols which do not belong to a leaf node
        leaves = [int(node) for node in leaves]
        for i in leaves:
            leaves_distance_matrix.append(distance_matrix[i])
        
        # print(leaves_distance_matrix)

        for i in range(len(leaves_distance_matrix)):
            temp = []
            for j in leaves:
                temp.append(
                    leaves_distance_matrix[i][j]
                )
            leaves_distance_matrix[i] = temp
        
        return leaves_distance_matrix
    
    def show_distance_matrix_between_leaves(self):
        to_show = ""
        for row in self.leaves_distance_matrix:
            to_show += universal_utils.print_arr(row) + "\n"
        
        to_show = to_show.strip()
        pyperclip.copy(to_show)
        return to_show

    # For each j, we can compute LimbLength(j) by 
    # finding the minimum value of ( Di,j + Dj,k − Di,k ) / 2
    # over all pairs of leaves i and k (where i ≠ j and k ≠ j).
    # leave_index: a 0-based index
    def get_limb_length(leaves_distance_matrix, leave_index):
        n = len(leaves_distance_matrix)
        
        limb_length = float('inf')
        for i in range(0, n):
            if i == leave_index:
                continue
            for j in range(0, n):
                if j == leave_index:
                    continue
                if j == i:
                    continue

                temp = (
                    leaves_distance_matrix[i][leave_index] + \
                    leaves_distance_matrix[leave_index][j] - \
                    leaves_distance_matrix[i][j]
                ) // 2

                if temp < limb_length:
                    limb_length = temp
        
        return limb_length

    # default show edge length in integer
    def show_tree_in_adjacent_list(self, format="d", mapping = None):
        adjacency_list = self.adjacency_list
        to_print = ""
        # 是否存在node的index和具体名称之间的关系
        if mapping: 
            for start_node in sorted(list(adjacency_list.keys()), key=lambda x: int(x)):
                for end_node in sorted(list(adjacency_list[start_node].keys()), key=lambda x: int(x)):
                    start_node_str = str(start_node)
                    end_node_str = str(end_node)
                    if start_node in mapping.keys(): 
                        start_node_str = mapping[start_node]
                    if end_node in mapping.keys(): 
                        end_node_str = mapping[end_node]
                    if format == "d":
                        to_print += "%s->%s:%d" %(start_node_str, end_node_str, adjacency_list[start_node][end_node]) + "\n"
                    elif format == "f":
                        to_print += "%s->%s:%.4f" %(start_node_str, end_node_str, adjacency_list[start_node][end_node]) + "\n"
        else:
            for start_node in sorted(list(adjacency_list.keys()), key=lambda x: int(x)):
                for end_node in sorted(list(adjacency_list[start_node].keys()), key=lambda x: int(x)):
                    if format == "d":
                        to_print += "%s->%s:%d" %(start_node, end_node, adjacency_list[start_node][end_node]) + "\n"
                    elif format == "f":
                        to_print += "%s->%s:%.4f" %(start_node, end_node, adjacency_list[start_node][end_node]) + "\n"
        
        to_print = to_print.strip()
        pyperclip.copy(to_print)
        print(to_print)
        return to_print
    
    def show_unweighed_tree(self, mapping = None, if_show=True): 
        adjacency_list = self.adjacency_list
        to_print = ""
        # 是否存在node的index和具体名称之间的关系
        if mapping: 
            for start_node in sorted(list(adjacency_list.keys()), key=lambda x: int(x)):
                for end_node in sorted(list(adjacency_list[start_node].keys()), key=lambda x: int(x)):
                    if start_node in mapping.keys(): 
                        start_node_str = mapping[start_node]
                    if end_node in mapping.keys(): 
                        end_node_str = mapping[end_node]
                    to_print += "%s->%s" %(start_node_str, end_node_str)+ "\n"
        else:
            for start_node in sorted(list(adjacency_list.keys()), key=lambda x: int(x)):
                for end_node in sorted(list(adjacency_list[start_node].keys()), key=lambda x: int(x)):
                    to_print += "%s->%s" %(start_node, end_node) + "\n"
        
        to_print = to_print.strip()
        pyperclip.copy(to_print)
        if if_show:
            print(to_print)
        return to_print

    
    def refresh_tree(self):
        self.refresh_graph()
        self.leaves = Tree.get_leaves(self.in_and_out_degrees)
        self.leaves_distance_matrix = Tree.get_distance_matrix_between_leaves(self.distance_matrix, self.leaves)

    # DFS path finding!
    def find_path(tree, start, end, path=None):
        if path is None:
            path = []
        path.append(start)

        if start == end: 
            return path
        
        for child in tree.adjacency_list[start].keys(): 
            if child not in path: 
                new_path = Tree.find_path(tree, child, end, path[:])
                if new_path: 
                    return new_path
        
        return None
    
    def copy_tree(tree): 
        new_tree = Tree("")
        new_tree.adjacency_list = copy.deepcopy(tree.adjacency_list)

        new_tree.refresh_tree()
        return new_tree

# AdditivePhylogeny(D)
#     n ← number of rows in D
#     if n = 2
#         return the tree consisting of a single edge of length D1,2
#     limbLength ← Limb(D, n)
#     for j ← 1 to n - 1
#         Dj,n ← Dj,n - limbLength
#         Dn,j ← Dj,n
#     (i, k) ← two leaves such that Di,k = Di,n + Dn,k
#     x ← Di,n
#     D ← D with row n and column n removed
#     T ← AdditivePhylogeny(D)
#     v ← the (potentially new) node in T at distance x from i on the path between i and k
#     add leaf n back to T by creating a limb (v, n) of length limbLength
#     return T
total_node_num_global = 0   # starting from the deepest recursive call, total_node_num_global += 1
def additive_phelogeny_generate_tree(leaves_distance_matrix, total_node_num):
    n = len(leaves_distance_matrix)

    # end of recursion: return the tree consisting of a single edge of length D0,1
    if n == 2:
        node_0 = "0"
        node_1 = "1"
        weight = leaves_distance_matrix[0][1]
        adjacency_list_str = node_0 + "->" + node_1 + ": " + str(weight) + "\n" + \
                            node_1 + "->" + node_0 + ": " + str(weight)
        tree = Tree(adjacency_list_str)

        global total_node_num_global
        total_node_num_global = total_node_num

        return tree

    # get the limb length of the last leaf node
    limb_length = Tree.get_limb_length(leaves_distance_matrix, n - 1)
    # create bold distance matrix
    for j in range(n):
        if j == n - 1:
            continue
        leaves_distance_matrix[j][n - 1] = leaves_distance_matrix[j][n - 1] - limb_length
        leaves_distance_matrix[n - 1][j] = leaves_distance_matrix[j][n - 1] 
    
    # two leaves make Di,k = Di,n + Dn,k
    selected_i_k = (0, 0)
    to_break = False
    for i in range(0, n - 1):
        for k in range(0, n - 1):
            if leaves_distance_matrix[i][k] == leaves_distance_matrix[i][n - 1] + leaves_distance_matrix[n - 1][k]: 
                selected_i_k = (i, k)
                to_break = True
                break
        if to_break:
            break

    # 要插入的位置出现在了i->k的路径上
    node_i = str(selected_i_k[0])
    node_k = str(selected_i_k[1])
    length_i_n = leaves_distance_matrix[selected_i_k[0]][n - 1]

    # remove the row and column of index = n - 1
    leaves_distance_matrix = [item[: n - 1] for item in leaves_distance_matrix[: n - 1]]

    # recursive call of this function
    tree = additive_phelogeny_generate_tree(leaves_distance_matrix, total_node_num)
    
    # the potentially new node in tree at distance x from i on the path between i and k
    new_node = str(total_node_num_global)
    node_n = str(n - 1)
    # find the path from i to k
    path = Tree.find_path(tree, node_i, node_k)
    path_len = 0
    current_node = path.pop(0)
    next_node = ""
    flag = -1    # 0: 加在已存在的internal node; 1: 新加一个internal node
    # find the (current_node, next_node) limb to add node_n
    while len(path) != 0:
        next_node = path.pop(0)
        path_len += tree.adjacency_list[current_node][next_node]

        if path_len == length_i_n:
            flag = 0
            break
        elif path_len > length_i_n: 
            flag = 1
            break

        current_node = next_node
    
    if flag == 0: 
    # add node_n to an existing internal node
        internal_node_to_add = path.pop(0)
        tree.adjacency_list[internal_node_to_add][node_n] = limb_length
        tree.adjacency_list[node_n] = {}
        tree.adjacency_list[node_n][internal_node_to_add] = limb_length
    elif flag == 1:
    # add node_n to new_node, and new_node should be added on the limb current_node - next_node
        internal_node_to_add = new_node
        
        tree.adjacency_list[internal_node_to_add] = {}
        tree.adjacency_list[node_n] = {}

        # current_node <-> internal_node_to_add
        if node_i == current_node: 
            tree.adjacency_list[current_node][internal_node_to_add] = length_i_n
            tree.adjacency_list[internal_node_to_add][current_node] = length_i_n
        else:
            tree.adjacency_list[current_node][internal_node_to_add] = length_i_n - (path_len - tree.adjacency_list[current_node][next_node])
            tree.adjacency_list[internal_node_to_add][current_node] = length_i_n - (path_len - tree.adjacency_list[current_node][next_node])
        # next_node <-> internal_node_to_add
        tree.adjacency_list[internal_node_to_add][next_node] = path_len - length_i_n
        tree.adjacency_list[next_node][internal_node_to_add] = path_len - length_i_n

        # internal_node_to_add <-> node_n
        tree.adjacency_list[internal_node_to_add][node_n] = limb_length
        tree.adjacency_list[node_n][internal_node_to_add] = limb_length

        del tree.adjacency_list[current_node][next_node]
        del tree.adjacency_list[next_node][current_node]

        total_node_num_global += 1
    
    return tree

if __name__ == "__main__":
    # testing the DNA_RNA_utils module
    # e.g. to find the origin of the DNA strand which an mRNA is transcribed from
    # mRNA_strand = DNA_RNA_utils.generate_random_strand(50, "RNA")
    # print(mRNA_strand)
    # transcription_origin_DNA_strand = DNA_RNA_utils.get_reverse_complement_strand(
    #     mRNA_strand, DNA_RNA_utils.complement_direction.RNA_to_DNA
    # )
    # print(transcription_origin_DNA_strand)

    # 1.2 step 11
#     node_num = 4
#     adjacency_list_str = """0->4:11
# 1->4:2
# 2->5:6
# 3->5:7
# 4->0:11
# 4->1:2
# 4->5:4
# 5->4:4
# 5->3:7
# 5->2:6"""
#     tree = Tree(adjacency_list_str)
#     tree.show_distance_matrix()
#     tree.show_distance_matrix_between_leaves()

    # with open("./dataset_30284_12.txt", "r") as f:
    #     f.readline()
    #     adjacency_list_str = f.read().strip()
    #     tree = Tree(adjacency_list_str)
    #     tree.show_distance_matrix_between_leaves()

    # 1.3 step 10:
#     arr_2_dim = """0	13	21	22
# 13	0	12	13
# 21	12	0	13
# 22	13	13	0"""
#     leaves_distance_matrix = universal_utils.parse_arr_2_dimension(arr_2_dim)
#     print(Tree.get_limb_length(leaves_distance_matrix, 1))

    # with open("./datasets/Limb_Length.txt", "r") as f:
    #     f.readline()
    #     n = int(f.readline().strip())
    #     leave_index = int(f.readline().strip())
    #     arr_2_dim = ""
    #     for _ in range(n):
    #         arr_2_dim += f.readline()
    #     arr_2_dim = arr_2_dim.strip()
    #     leaves_distance_matrix = universal_utils.parse_arr_2_dimension(arr_2_dim)
    #     result = Tree.get_limb_length(leaves_distance_matrix, leave_index)
    #     f.readline()
    #     correct = int(f.readline().strip())
    #     if correct == result:
    #         print("result = %d, correct!" %(result))

    # with open("./datasets/dataset_30285_11.txt", "r") as f:
    #     f.readline()
    #     leave_index = int(f.readline().strip())
    #     arr_2_dim = f.read().strip()
    #     leaves_distance_matrix = universal_utils.parse_arr_2_dimension(arr_2_dim)
    #     result = Tree.get_limb_length(leaves_distance_matrix, leave_index)
    #     print("result: %d" %result)

    # 1.4 step 6
    leaves_distance_matrix = universal_utils.parse_arr_2_dimension(
        """0	13	21	22
13	0	12	13
21	12	0	13
22	13	13	0"""
    )
    tree = additive_phelogeny_generate_tree(leaves_distance_matrix, 4)
    tree.refresh_tree()
    Tree.show_tree_in_adjacent_list(tree)

    # extra dataset
    with open("./datasets/Additive_Phylogeny.txt", "r") as f:
        f.readline()
        total_node_num = int(f.readline().strip())
        leaves_distance_matrix_str = ""
        temp = f.readline()
        while temp != "Output\n": 
            leaves_distance_matrix_str += temp
            temp = f.readline()
        leaves_distance_matrix_str = leaves_distance_matrix_str.strip()
        leaves_distance_matrix = universal_utils.parse_arr_2_dimension(leaves_distance_matrix_str)
        tree = additive_phelogeny_generate_tree(leaves_distance_matrix, total_node_num)
        Tree.show_tree_in_adjacent_list(tree)

    with open("./datasets/dataset_30286_6.txt", "r") as f:
        total_node_num = int(f.readline().strip())
        leaves_distance_matrix = universal_utils.parse_arr_2_dimension(f.read().strip())
        tree = additive_phelogeny_generate_tree(leaves_distance_matrix, total_node_num)
        tree.refresh_tree()
        Tree.show_tree_in_adjacent_list(tree)

    # 1.4 step 11
    # SARS-like coronaviruse distance matrix
    # with open("./datasets/coronavirus_distance_matrix_additive.txt", "r") as f:
    #     species = universal_utils.parse_arr(f.readline().strip())
    #     # universal_utils.print_arr(species)
    #     leaves_distance_matrix = []
    #     for line in f.read().strip().split("\n"): 
    #         leaves_distance_matrix.append(
    #             [int(item) for item in line.split("\t")[1: ]]
    #         )
    #     tree = additive_phelogeny_generate_tree(leaves_distance_matrix, len(species))
    #     tree.refresh_tree()
    #     Tree.show_tree_in_adjacent_list(tree)

    # Quiz
    # problem 5
    # leaves_distance_matrix = universal_utils.parse_arr_2_dimension(
    #     "0\t13\t16\t10\n13\t0\t21\t15\n16\t21\t0\t18\n10\t15\t18\t0"
    # )
    # limb_length = Tree.get_limb_length(leaves_distance_matrix, 3)
    # print("limb length of node \"l\": %d" %limb_length)