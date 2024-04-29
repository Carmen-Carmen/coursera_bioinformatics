from itertools import product
import time
import pyperclip
import random
import timeit
import os
import copy

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

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

def find_Eulerian_cycle_by_graph(graph):
    # method 1 不太对好像，无法避免角标越界
    # EulerianCycle(Graph)
    # form a cycle Cycle by randomly walking in Graph (don't visit the same edge twice!)
    # while there are unexplored edges in Graph
    #     select a node newStart in Cycle with still unexplored edges
    #     form Cycle’ by traversing Cycle (starting at newStart) and then randomly walking 
    #     Cycle ← Cycle’
    # return Cycle
    # graph_copy = graph.copy()

    # cycle = []
    # start_node = random.choice(list(graph.keys()))
    # next_node = graph[start_node][0]
    # # the edge of start_node -> next_node should be marked as traversed (remove)
    # graph[start_node].remove(next_node)
    # # add the two nodes to cycle list
    # cycle.append(start_node)
    # cycle.append(next_node)

    # # prev_node records "next_node" when each time enters the while loop
    # prev_node = next_node
    # while next_node != start_node:
    #     print(next_node)
    #     prev_node = next_node
    #     next_node = graph[next_node][0]
    #     graph[prev_node].remove(next_node)
    #     cycle.append(next_node)
    # # when it gets out of the while loop, next_node == start_node
    # graph[prev_node].remove(next_node)
    # cycle.append(next_node)

    # while not(is_all_traversed(graph)):
    #     # select a node newStart in Cycle with still unexplored edges
    #     for node in graph.keys():
    #         if len(graph[node]) != 0:
    #             start_node = node
    #             break
    #     graph = graph_copy.copy()
    #     next_node = graph[start_node][0]

    #     # the edge of start_node -> next_node should be marked as traversed (remove)
    #     graph[start_node].remove(next_node)
    #     # add the two nodes to cycle list
    #     cycle.append(start_node)
    #     cycle.append(next_node)

    #     # prev_node records "next_node" when each time enters the while loop
    #     prev_node = next_node
    #     while next_node != start_node:
    #         prev_node = next_node
    #         next_node = graph[next_node][0]
    #         graph[prev_node].remove(next_node)
    #         cycle.append(next_node)
    #     # when it gets out of the while loop, next_node == start_node
    #     graph[prev_node].remove(next_node)
    #     cycle.append(next_node)
    
    # return cycle

    # method 2
    # http://www.graph-magics.com/articles/euler.php
    # Algorithm for directed graphs:
    # 1. Start with an empty stack and an empty circuit (eulerian path).
    # - If all vertices have same out-degrees as in-degrees - choose any of them.
    # - If all but 2 vertices have same out-degree as in-degree, and one of those 2 vertices has 
    #   out-degree with one greater than its in-degree, and the other has in-degree with one greater 
    #   than its out-degree - then choose the vertex that has its out-degree with one greater than its in-degree.
    # - Otherwise no euler circuit or path exists.
    in_and_out_degree = get_in_and_out_degree(graph)
    is_in_and_out_same = True
    start_node = ""
    for node in in_and_out_degree.keys():
        in_degree = in_and_out_degree[node]["in-degree"]
        out_degree = in_and_out_degree[node]["out-degree"]
        if in_degree != out_degree:
            is_in_and_out_same = False
            if out_degree == in_degree + 1:
                start_node = node
                break
    
    if is_in_and_out_same:
        start_node = random.choice(list(graph.keys()))
    elif start_node == "":
        # 说明经过遍历graph的所有key，都没有找到out-degree = in-degree + 1, 
        # 且还出现了in和out不等的情况，那就是不能形成eulerian circuit
        print("graph unable to form Eulerian cycle!")
        return None
    
    # 2. - If current vertex has no out-going edges (i.e. neighbors) - add it to circuit, 
    #       remove the last vertex from the stack and set it as the current one. 
    # - Otherwise (in case it has out-going edges, i.e. neighbors) - add the vertex to the stack, 
    #   take any of its neighbors, remove the edge between that vertex and selected neighbor, 
    #   and set that neighbor as the current vertex.
    # 3. Repeat step 2 until the current vertex has no more out-going edges (neighbors) and the stack is empty.
    stack = []
    current_node = start_node
    cycle = []

    if len(graph[current_node]) != 0: 
        stack.append(current_node)
        next_node = random.choice(graph[current_node])
        graph[current_node].remove(next_node)
        current_node = next_node

    while len(stack) != 0:
        if len(graph[current_node]) != 0:
            stack.append(current_node)
            next_node = random.choice(graph[current_node])
            graph[current_node].remove(next_node)

            current_node = next_node
        else:
            cycle.append(current_node)
            current_node = stack.pop()
    
    if cycle[0] != current_node:
        print("graph unable to form Eulerian cycle!") 
        return None
    
    # beware that the tour (cycle) is given in revers order
    cycle.append(current_node)
    cycle.reverse()

    return cycle

def find_Eulerian_path_by_graph(graph):
    # complete_graph_node(graph)
    # graph_copy = graph.copy()
    # in_and_out_degree = get_in_and_out_degree(graph)
    # start_node = ""
    # end_node = ""
    # for node in in_and_out_degree.keys():
    #     in_degree = in_and_out_degree[node]["in-degree"]
    #     out_degree = in_and_out_degree[node]["out-degree"]
    #     if in_degree != out_degree:
    #         if out_degree == in_degree + 1:
    #             start_node = node
    #         elif out_degree + 1 == in_degree:
    #             end_node = node

    # if start_node == "":
    #     # 说明经过遍历graph的所有key，都没有找到out-degree = in-degree + 1, 
    #     # 说明不能形成Eulerian path
    #     print("graph unable to form Eulerian path!")
    #     return None
    
    # # 2. - If current vertex has no out-going edges (i.e. neighbors) - add it to circuit, 
    # #       remove the last vertex from the stack and set it as the current one. 
    # # - Otherwise (in case it has out-going edges, i.e. neighbors) - add the vertex to the stack, 
    # #   take any of its neighbors, remove the edge between that vertex and selected neighbor, 
    # #   and set that neighbor as the current vertex.
    # # 3. Repeat step 2 until the current vertex has no more out-going edges (neighbors) and the stack is empty.
    # while True:
    #     stack = []
    #     current_node = start_node
    #     path = []

    #     if len(graph[current_node]) != 0: 
    #         stack.append(current_node)
    #         next_node = random.choice(graph[current_node])
    #         graph[current_node].remove(next_node)
    #         current_node = next_node

    #     while len(stack) != 0:
    #         if len(graph[current_node]) != 0:
    #             stack.append(current_node)
    #             next_node = random.choice(graph[current_node])
    #             graph[current_node].remove(next_node)

    #             current_node = next_node
    #         else:
    #             path.append(current_node)
    #             current_node = stack.pop()
        
    #     if is_all_traversed(graph) and path[0] == end_node:
    #         # beware that the tour (path) is given in revers order
    #         path.append(current_node)
    #         path.reverse()

    #         return path
    #     else:
    #         graph = graph_copy.copy()

    # Another method from the comments
    # 1. find the start_node (out-degree = in-degree + 1) and end_node (in-degree = out-degree + 1)
    complete_graph_node(graph)

    in_and_out_degree = get_in_and_out_degree(graph)
    start_node = ""
    end_node = ""
    for node in in_and_out_degree.keys():
        if start_node != "" and end_node != "":
            break
        else:
            in_degee = in_and_out_degree[node]["in-degree"]
            out_degee = in_and_out_degree[node]["out-degree"]
            if in_degee == out_degee + 1:
                end_node = node
            elif in_degee + 1 == out_degee:
                start_node = node
    
    path = []
    if start_node == "" or end_node == "":
        # print("graph unable to form Eulerian path!")
        # 只是说这个graph是个cycle，但是也不能说不存在Eulerian path呀！
        cycle = find_Eulerian_cycle_by_graph(graph)
        # print(cycle)
        # 此时这些连接是真实存在的，因此直接返回cycle即可
        path = cycle
    else:
    # 2. add a path from end_node to start_node
        graph[end_node].append(start_node)


    # 3. now that there must be an Eulerian cycle in the graph
        # perform the find_Eulerian_cycle algorithm, and adjust the cycle to begin with start_node
        # remove the added edge
        cycle = find_Eulerian_cycle_by_graph(graph)
        cycle.pop()
        # print(cycle)
        path = []
        if cycle[0] == start_node and cycle[-1] == end_node:
            path = cycle
        else:
            start_index = cycle.index(start_node)
            end_index = cycle.index(end_node)
            # 存在1个以上的start_node，则end_node必定只有一个，可通过end_node在cycle数组中的index来确定start_index
            # vice versa
            if cycle.count(start_node) != 1:
                end_index = cycle.index(end_node)
                start_index = end_index + 1
            elif cycle.count(end_node) != 1:
                start_index = cycle.index(start_node)
                end_index = start_index - 1

            # construct the path
            for i in range(start_index, len(cycle)):
                path.append(cycle[i])
            for i in range(0, start_index):
                path.append(cycle[i])
                
    # print(path)
    return path
        
# get the in-degree and out-degree of each node from a dict-like graph
def get_in_and_out_degree(graph):
    in_and_out_degree = {}
    for key in graph.keys():
        in_and_out_degree[key] = {
            "in-degree" : 0, 
            "out-degree" : len(graph[key])
        }

    for key in graph.keys():
        for node in graph[key]:
            # this is a list
            in_and_out_degree[node]["in-degree"] += 1
    
    # result = ""
    # with open("degree.log", "a") as f:
    #     f.write("\n")
    #     f.write("######## %s ########\n" %time.time())
    #     for key in in_and_out_degree: 
    #         result += "%s: in-degree (%d), out-degree (%d)\n" \
    #             %(key, in_and_out_degree[key]["in-degree"], in_and_out_degree[key]["out-degree"])
        
    #     f.write(result)

    return in_and_out_degree
    
def complete_graph_node(graph):
    nodes = set()
    for key in graph.keys():
        for node in graph[key]:
            nodes.add(node)
    
    for node in nodes:
        if not(node in graph.keys()):
            graph[node] = []
    
    return

def is_all_traversed(graph):
    for val in graph.values():
        if len(val) != 0:
            return False
    
    return True

# the brute method for string construction problem
def string_construction_brute(composition, k):
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

    # print_graph_dict(de_Bruijn_graph)
    return de_Bruijn_graph

# a new method for string construction problem
# 1. form a de Bruijn graph from given composition
# 2. find a Eulerian path in the de Bruijn graph
# 3. construct the origin string by the path
def string_construction_new(composition):
    de_Bruijn_graph = get_de_Bruijn_graph_by_composition(composition)
    return string_construction_by_de_Bruijn_graph(de_Bruijn_graph)

def string_construction_by_de_Bruijn_graph(de_Bruijn_graph):
    path = find_Eulerian_path_by_graph(de_Bruijn_graph)
    if path is None:
        # print("err in generating path")
        return None
    # print(path)
    ori_str = ""
    ori_str += path[0][:len(path[0]) - 1]
    for node in path:
        ori_str += node[-1]
    
    return ori_str

def get_k_universal_circular_string(k):
    composition = ["".join(p) for p in product("01", repeat=(k))]
    de_Bruijn_graph = get_de_Bruijn_graph_by_composition(composition)
    cycle = find_Eulerian_cycle_by_graph(de_Bruijn_graph)
    ori_str = ""
    ori_str += cycle[0][:len(cycle[0]) - 1]
    for node in cycle:
        ori_str += node[-1]
    
    return ori_str[: len(ori_str) - (k - 1)]

# generate the (k, d)-mer composition of a given strand
def get_paired_kmers_composition(strand, k, d):
    pairs = []
    for i in range(len(strand) - (2 * k + d) + 1):
        temp = strand[i : i + (2 * k + d)]
        read1 = temp[: k]
        read2 = temp[-k :]
        pairs.append("(%s|%s)" %(read1, read2))
    
    pairs.sort()
    result = ""
    for pair in pairs:
        result += pair + " "
    result.strip()
    return result

# StringSpelledByGappedPatterns problem
# composition should be given in the form below: 
# GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA
def string_construction_by_paired_kmers_composition(paired_composition, k, d):
    # # 错误的本质：
    # # 本质上错在了没有利用paired的特性，把pair好的两个edge又拆开来分别建立de Bruijn图了
    # # 应当先基于paired_composition (即左右都要符合kmer[1 :] = kmer'[: k - 1]) 建立de Bruijn图
    # paired_composition = paired_composition.split(" ")
    # prefix_composition = []
    # suffix_composition = []
    # for pair in paired_composition:
    #     pair = pair.split("|")
    #     prefix_composition.append(pair[0])
    #     suffix_composition.append(pair[1])
    

    # # 现在的问题：
    # # 目前的算法会根据prefix_composition和suffix_composition各自分别生成一条strand
    # # 在生成过程中可能会分别随机从两个composition中取node作为start和end (如果所有node都是balanced)
    # # 然而实际上应当根据同一个pair作为开头来生成
    # prefix_strand = string_construction_new(prefix_composition)
    # # print(prefix_strand + "\n")
    # suffix_strand = string_construction_new(suffix_composition)
    # # print(suffix_strand + "\n")
    # if prefix_strand is None or suffix_strand is None:
    #     print("err in generating strand!")
    #     return
    # # 有可能产生其他路线，这个时候就需要通过d来判断是否正确重合
    # while prefix_strand[k + d :] != suffix_strand[: - (k + d)]:
    #     # print(prefix_strand[k + d :])
    #     # print(suffix_strand[: - (k + d)])
    #     prefix_strand = string_construction_new(prefix_composition)
    #     # print(prefix_strand + "\n")
    #     suffix_strand = string_construction_new(suffix_composition)
    #     # print(suffix_strand + "\n")
    
    # return prefix_strand[: k + d] + suffix_strand

    # 重新写一遍
    # paired_de_Bruijn_graph = get_de_Bruijn_graph_by_paired_composition(paired_composition, k)
    # # 将这个paired图拆解为2个de_Bruijn_graph
    # de_Bruijn_graph1 = split_paired_de_Bruijn_graph(paired_de_Bruijn_graph)[0]
    # de_Bruijn_graph2 = split_paired_de_Bruijn_graph(paired_de_Bruijn_graph)[1]
    
    # print("paired:")
    # print_graph_dict(paired_de_Bruijn_graph)
    # print("graph1:")
    # print_graph_dict(de_Bruijn_graph1)
    # print("graph2:")
    # print_graph_dict(de_Bruijn_graph2)

    # # 深拷贝太慢了，还是直接生成吧
    # # de_Bruijn_graph1_copy = copy.deepcopy(de_Bruijn_graph1)
    # # de_Bruijn_graph2_copy = copy.deepcopy(de_Bruijn_graph2)

    # prefix_strand = string_construction_by_de_Bruijn_graph(de_Bruijn_graph1)
    # suffix_strand = string_construction_by_de_Bruijn_graph(de_Bruijn_graph2)

    # print("prefix:", prefix_strand)
    # print("suffix:", suffix_strand)
    # print()

    # # 有可能产生其他路线，这个时候就需要通过d来判断是否正确重合
    # while prefix_strand[k + d :] != suffix_strand[: - (k + d)]:
    #     # 在find path的过程中graph里的edge都被毁掉了
    #     # 注意一定要用深拷贝啊，因为我的val是list对象，如果用dict.copy()的浅拷贝，依然会对原来dict对象所指向的list的内存空间进行修改
    #     # 但是深拷贝太慢了
    #     # de_Bruijn_graph1 = copy.deepcopy(de_Bruijn_graph1_copy)
    #     # de_Bruijn_graph2 = copy.deepcopy(de_Bruijn_graph2_copy)

    #     de_Bruijn_graph1 = split_paired_de_Bruijn_graph(paired_de_Bruijn_graph)[0]
    #     de_Bruijn_graph2 = split_paired_de_Bruijn_graph(paired_de_Bruijn_graph)[1]

    #     print("paired:")
    #     print_graph_dict(paired_de_Bruijn_graph)
    #     print("graph1:")
    #     print_graph_dict(de_Bruijn_graph1)
    #     print("graph2:")
    #     print_graph_dict(de_Bruijn_graph2)

    #     # print(prefix_strand[k + d :])
    #     # print(suffix_strand[: - (k + d)])
    #     prefix_strand = string_construction_by_de_Bruijn_graph(de_Bruijn_graph1)
    #     suffix_strand = string_construction_by_de_Bruijn_graph(de_Bruijn_graph2)
    #     print("prefix:", prefix_strand)
    #     print("suffix:", suffix_strand)
    #     print()
    
    # return prefix_strand[: k + d] + suffix_strand

    # 其实没有这么复杂，只需要用paired_de_Bruijn_graph去生成一个path，再拆成俩字符串即可
    paired_de_Bruijn_graph = get_de_Bruijn_graph_by_paired_composition(paired_composition, k)
    # print("paired graph:")
    # print_graph_dict(paired_de_Bruijn_graph)
    path = find_Eulerian_path_by_graph(paired_de_Bruijn_graph)

    prefix_strand = path[0].split("|")[0]
    suffix_strand = path[0].split("|")[1]
    for i in range(1, len(path)):
        prefix_strand += path[i].split("|")[0][-1]
        suffix_strand += path[i].split("|")[1][-1]
    
    # print(prefix_strand, suffix_strand)
    
    while prefix_strand[k + d :] != suffix_strand[: - (k + d)]:
        paired_de_Bruijn_graph = get_de_Bruijn_graph_by_paired_composition(paired_composition, k)
        # print("paired graph:")
        # print_graph_dict(paired_de_Bruijn_graph)
        path = find_Eulerian_path_by_graph(paired_de_Bruijn_graph)

        prefix_strand = path[0].split("|")[0]
        suffix_strand = path[0].split("|")[1]
        for i in range(1, len(path)):
            prefix_strand += path[i].split("|")[0][-1]
            suffix_strand += path[i].split("|")[1][-1]
    
        # print(prefix_strand, suffix_strand)

    return prefix_strand[: k + d] + suffix_strand
    
def get_de_Bruijn_graph_by_paired_composition(paired_composition, k):
    de_Bruijn_graph = {}
    pairs = paired_composition.split(" ")
    nodes = set()
    for pair in pairs:
        kmer1 = pair.split("|")[0]
        kmer2 = pair.split("|")[1]

        prefix = "%s|%s" %(kmer1[: k - 1], kmer2[: k - 1])
        suffix = "%s|%s" %(kmer1[1 :], kmer2[1 :])

        nodes.add(prefix)
        nodes.add(suffix)
    
    # print(nodes)

    for node in nodes:
        de_Bruijn_graph[node] = []

        for pair in pairs:
            kmer1 = pair.split("|")[0]
            kmer2 = pair.split("|")[1]

            prefix = "%s|%s" %(kmer1[: k - 1], kmer2[: k - 1])
            suffix = "%s|%s" %(kmer1[1 :], kmer2[1 :])

            if  node.split("|")[0] == prefix.split("|")[0] and \
                node.split("|")[1] == prefix.split("|")[1]:
                de_Bruijn_graph[node].append(suffix)
            
    return de_Bruijn_graph

def split_paired_de_Bruijn_graph(paired_de_Bruijn_graph):
    # print("spliting...")

    de_Bruijn_graph1 = {}
    de_Bruijn_graph2 = {}
    for node in paired_de_Bruijn_graph.keys():
        upper_node = node.split("|")[0]
        lower_node = node.split("|")[1]

        de_Bruijn_graph1[upper_node] = []
        de_Bruijn_graph2[lower_node] = []
    
    for node in paired_de_Bruijn_graph.keys():
        upper_node = node.split("|")[0]
        lower_node = node.split("|")[1]

        for next_node in paired_de_Bruijn_graph[node]:
            # print(next_node)
            upper_next_node = next_node.split("|")[0]
            lower_next_node = next_node.split("|")[1]

            de_Bruijn_graph1[upper_node].append(upper_next_node)
            de_Bruijn_graph2[lower_node].append(lower_next_node)
    
    return [de_Bruijn_graph1, de_Bruijn_graph2]

# Contig Generation Problem: Generate the contigs from a collection of reads (with imperfect coverage).
# i.e. reads compose an imperfect coverage of the genome
def get_contigs_by_reads(reads):
    # get de Bruijn graph
    de_Bruijn_graph = get_de_Bruijn_graph_by_composition(reads)

    maximal_non_branching_paths = get_maximal_non_branching_paths_by_de_Bruijn_graph(de_Bruijn_graph)

    contigs = []
    for path in maximal_non_branching_paths:
        contig = path[0]
        path = path[1 :]
        for nodes in path:
            contig += nodes[-1]

        contigs.append(contig)
    
    return contigs

    return get_contigs_by_de_Bruijn_graph(de_Bruijn_graph)

def get_maximal_non_branching_paths_by_de_Bruijn_graph(de_Bruijn_graph):
    complete_graph_node(de_Bruijn_graph)
    # print_graph_dict(de_Bruijn_graph)
    # for node in in_and_out_degree.keys():
    #     in_degree = in_and_out_degree[node]["in-degree"]
    #     out_degree = in_and_out_degree[node]["out-degree"]
    #     print("%s: in-degree(%d), out-degree(%d)" %(node, in_degree, out_degree))

    maximal_non_branching_paths = []
    # definition of maximal non-branching path:
    # a non-branching path that cannot be extended into a longer non-branching path
    # i.e. start from a node v which one of in(v) and out(v) does not equals 1
    # end in a node w which one of in(w) and out(w) does not equals 1
    for node in de_Bruijn_graph.keys():
        if len(de_Bruijn_graph[node]) == 0:
            continue

        if not(is_non_branching_node(de_Bruijn_graph, node)):
            # this node can be a starting node of a maximal non-branching path
            for out_node in de_Bruijn_graph[node]:
                next_node = out_node
                path = []
                path.append(node)

                path.append(next_node)
                while is_non_branching_node(de_Bruijn_graph, next_node):
                    out_nodes = de_Bruijn_graph[next_node]
                    if len(out_nodes) == 0:
                        break
                    else:
                        # if this node is a non-branching node, then there must be only 1 node in out_nodes
                        next_node = out_nodes[0]
                    path.append(next_node)
                        
                if len(path) != 1:
                    maximal_non_branching_paths.append(path)
    
    # print(maximal_non_branching_paths)
    # 再单独处理环结构之前，去掉所有之前已经遍历过的edge
    for path in maximal_non_branching_paths:
        for i in range(len(path) - 1):
            node = path[i]
            next_node = path[i + 1]

            de_Bruijn_graph[node].remove(next_node)

    # 单独处理in和out都是1的环结构
    cycles = []
    for node in [item for item in de_Bruijn_graph.keys() if is_non_branching_node(de_Bruijn_graph, item)]:
        # 必须是单独的环，所以加入这一步判断，这样不会把同一个环以不同起点开始的结果作为新的环
        is_in_cycle = False
        for cycle in cycles:
            if node in cycle:
                is_in_cycle = True
        if is_in_cycle:
            continue

        for out_node in de_Bruijn_graph[node]:
            cycle = []
            next_node = out_node

            cycle.append(node)
            cycle.append(next_node)
            while is_non_branching_node(de_Bruijn_graph, next_node):
                next_node = de_Bruijn_graph[next_node][0]

                if next_node == node:
                    cycle.append(next_node)
                    cycles.append(cycle)
                    break
                cycle.append(next_node)
    
    for cycle in cycles:
        maximal_non_branching_paths.append(cycle)
    
    return maximal_non_branching_paths

def is_non_branching_node(graph, node):
    in_and_out_degree = get_in_and_out_degree(graph)
    in_degree = in_and_out_degree[node]["in-degree"]
    out_degree = in_and_out_degree[node]["out-degree"]

    return in_degree == 1 and out_degree == 1

# split original reads to smaller length in order to get larger coverage of the genome
def split_reads_to_k_len(reads, k):
    ori_read_len = len(reads[0])
    if k > ori_read_len: 
        print("error k to split reads to! k should be smaller than the origin length of each read")
    
    splitted_reads = []
    for read in reads:
        for i in range(ori_read_len - k + 1):
            splitted_reads.append(read[i : i + k])

    return splitted_reads

def split_paired_reads_to_k_len(paired_reads, k):
    ori_read_len = len(paired_reads[0].split("|")[0])
    if k > ori_read_len: 
        print("error k to split reads to! k should be smaller than the origin length of each read")
    
    splitted_paired_reads = []
    for paired_read in paired_reads:
        prefix = paired_read.split("|")[0]
        suffix = paired_read.split("|")[1]
        for i in range(ori_read_len - k + 1):
            new_paired_read = prefix[i : i + k] + "|" + suffix[i : i + k]
            splitted_paired_reads.append(new_paired_read)
    
    return splitted_paired_reads

if __name__ == "__main__":
    # 1.2 step 2
    # graph = {
    #     "0": ["3"],
    #     "1": ["0"],
    #     "2": ["1", "6"],
    #     "3": ["2"],
    #     "4": ["2"],
    #     "5": ["4"],
    #     "6": ["5", "8"],
    #     "7": ["9"],
    #     "8": ["7"],
    #     "9": ["6"],
    # }

    # # print_graph_dict(graph)

    # # print(get_in_and_out_degree(graph))
    # print_arr(find_Eulerian_cycle_by_graph(graph))
    # with open("dataset_203_2.txt", "r") as f:
    #     graph = {}
    #     line = f.readline().strip()
    #     while line != "":
    #         parsing = line.split(": ")
    #         graph[parsing[0]] = parsing[1].split(" ")
    #         line = f.readline().strip()
        
    #     print_arr(find_Eulerian_cycle_by_graph(graph))

    # 1.2 step 6
#     lines =  \
# """0: 2
# 1: 3
# 2: 1
# 3: 0 4
# 6: 3 7
# 7: 8
# 8: 9
# 9: 6""".split("\n")
#     graph = {}
#     for line in lines:
#         parsing = line.split(": ")
#         key = parsing[0]
#         val = parsing[1].split(" ")
#         graph[key] = val
    
#     print_arr(find_Eulerian_path_by_graph(graph))
    # with open("dataset_203_6.txt", "r") as f:
    #     graph = {}
    #     line = f.readline().strip()
    #     while line != "":
    #         parsing = line.split(": ")
    #         graph[parsing[0]] = parsing[1].split(" ")
    #         line = f.readline().strip()
        
    #     print_arr(find_Eulerian_path_by_graph(graph))

    # 1.2 step 7
    # composition = "CTTA ACCA TACC GGCT GCTT TTAC".split(" ")
    # k = 4
    # print(string_construction_brute(composition, 4))
    # print(string_construction_new(composition))
    # with open("dataset_203_7.txt", "r") as f:
    #     k = int(f.readline().strip())
    #     composition = f.readline().strip().split(" ")

    #     print("###### new method ######")
    #     start = timeit.default_timer()
    #     str_new = string_construction_new(composition)
    #     end = timeit.default_timer()
    #     print(str_new)
    #     print("new method runtime:", str((end - start) * 1000), "ms")
    #     print("\n###### brute method ######")
    #     start = timeit.default_timer()
    #     str_brute = string_construction_brute(composition, k)[0]
    #     end = timeit.default_timer()
    #     print(str_brute)
    #     print("brute method runtime:", str((end - start) * 1000), "ms")
    #     print("\n" + str(str_new == str_brute))

    # 1.2 step 10
    # kmers = ["".join(p) for p in product("01", repeat=20)]
    # count = 0
    # for kmer in kmers:
    #     for other_kmer in kmers:
    #         if kmer[1 :] == other_kmer[: 19]:
    #             count += 1

    # print("num of edges:", count)

    # 1.2 step 11
    # print(get_k_universal_circular_string(3))
    # with open("dataset_203_11.txt", "r") as f:
    #     k = int(f.read().strip())
    #     print(get_k_universal_circular_string(k))

    # 1.3 step 6
    # strand = "TAATGCCATGGGATGTT"
    # k = 3
    # d = 2
    # print(get_paired_kmers_composition(strand, k, d))

    # 1.5 step 4
    # print(string_construction_by_paired_kmers_composition("GACC|GCGC ACCG|CGCC CCGA|GCCG CGAG|CCGG GAGC|CGGA", 4, 2))
    # print(string_construction_by_paired_kmers_composition("TCA|GCA TTC|TGC AAT|CAT ATT|ATG", 3, 1))
    # print(string_construction_by_paired_kmers_composition("GG|GA GT|AT TG|TA GA|AC AT|CT", 2, 1))
    # print(string_construction_by_paired_kmers_composition("GGG|GGG AGG|GGG GGG|GGT GGG|GGG GGG|GGG", 3, 2))
    # debuging
    # os.chdir("./StringReconstructionReadPairs")
    # inputs = os.listdir("./inputs")
    # outputs = os.listdir("./outputs")
    # try:
    #     inputs.remove(".DS_Store")
    #     outputs.remove(".DS_Store")
    # except Exception:
    #     pass
    # inputs.sort()
    # outputs.sort()
    # # print(inputs)
    # # print(outputs)
    # counter = 0
    # for input in inputs:
    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         params = f_input.readline().strip().split(" ")
    #         k = int(params[0])
    #         d = int(params[1])
    #         paired_composition = f_input.readline().strip()

    #         result = string_construction_by_paired_kmers_composition(paired_composition, k, d)

    #         # print(paired_composition)
        
    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.readline().strip()
    #         print("result:", result)
    #         print("correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")
        
    #     counter += 1
    #     print()

    # with open("dataset_6206_4.txt", "r") as f:
    #     params = f.readline().strip().split(" ")
    #     k = int(params[0])
    #     d = int(params[1])
    #     paired_composition = f.readline().strip()

    #     print(string_construction_by_paired_kmers_composition(paired_composition, k, d))

    # with open("dataset_204_16.txt", "r") as f:
    #     params = f.readline().strip().split(" ")
    #     k = int(params[0])
    #     d = int(params[1])
    #     paired_composition = f.readline().strip()

    #     print(string_construction_by_paired_kmers_composition(paired_composition, k, d))

    # 1.4 step 5
    # reads = "ATG ATG TGT TGG CAT GGA GAT AGA".split(" ")
    # print(get_contigs_by_reads(reads))

    # debuging
    # os.chdir("./ContigGeneration")
    # inputs = os.listdir("./inputs")
    # outputs = os.listdir("./outputs")
    # try:
    #     inputs.remove(".DS_Store")
    #     outputs.remove(".DS_Store")
    # except Exception:
    #     pass
    # inputs.sort()
    # outputs.sort()
    # # print(inputs)
    # # print(outputs)
    # counter = 0
    # for input in inputs:
    #     counter += 1
    #     # if "4" not in input:
    #     #     continue

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         reads = f_input.readline().strip().split(" ")

    #         result = ""
    #         for contig in get_contigs_by_reads(reads):
    #             result += contig + " "
    #         result = result.strip()

    #         # print(paired_composition)
        
    #     with open("./outputs/%s" %outputs[counter - 1], "r") as f_output:
    #         output = f_output.readline().strip()
    #         print("result:", result)
    #         print("correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter), "correct!")
    #         else:
    #             print("debug_%d" %(counter), "wrong!")
        
    #     print()

    # with open("dataset_205_5.txt", "r") as f:
    #     reads = f.readline().strip().split(" ")
    #     print_arr(get_contigs_by_reads(reads))

    # with open("dataset_6207_2.txt", "r") as f:
    #     graph = {}
    #     lines = f.read().split("\n")
    #     lines.remove("")
    #     # print(lines)
    #     for line in lines:
    #         parsing = line.split(": ")

    #         graph[parsing[0]] = parsing[1].split(" ")
        
    #     # print_graph_dict(graph)

    #     for path in get_maximal_non_branching_paths_by_de_Bruijn_graph(graph):
    #         print_arr(path)

    # 1.4 step 10 ### NOT SOLVED YET ###
    with open("./paired_reads_Carsonella_ruddii.txt", "r") as f:
        paired_composition_list = f.read().strip().split("\n")

        splitted_paired_reads = split_paired_reads_to_k_len(paired_composition_list, 35)

        paired_composition = ""
        for item in splitted_paired_reads:
            paired_composition += item + " "
        paired_composition = paired_composition.strip()


        k = 35
        d = 120 + (100 - k)

        print(string_construction_by_paired_kmers_composition(paired_composition, k, d))

        # treat the reads as kmers instead of kd-mers
        # reads = []
        # for pair in paired_composition_list:
        #     reads.append(pair.split("|")[0])
        #     reads.append(pair.split("|")[1])
        
        # splitted_reads = split_reads_to_k_len(reads, 50)
        # print_arr(splitted_reads)
        # print_graph_dict(get_de_Bruijn_graph_by_composition(splitted_reads))
        
        # print(string_construction_new(splitted_reads))
        # print(get_contigs_by_reads(splitted_reads))


    # Week 2 Quiz
    # print("###### problem 1 ######")
    # composition = "AAAT  AATG  ACCC  ACGC  ATAC  ATCA  ATGC  CAAA  CACC  CATA  CATC  CCAG  CCCA  CGCT  CTCA  GCAT  GCTC  TACG  TCAC  TCAT  TGCA".split("  ")
    # print(string_construction_new(composition))
    # print(string_construction_brute(composition, 4))

    # print("\n###### problem 2 ######")
    # raw = "1 -> 2,3,5  2 -> 1,4  3 -> 2,5  4 -> 1,2,5  5 -> 3,4".split("  ")
    # de_Bruijn_graph = {}
    # for item in raw:
    #     key = item.split(" -> ")[0]
    #     val = item.split(" -> ")[1].split(",")
    #     de_Bruijn_graph[key] = val
    
    # # print_graph_dict(de_Bruijn_graph)
    # in_and_out_degree = get_in_and_out_degree(de_Bruijn_graph)
    # for key in in_and_out_degree.keys():
    #     print("%s: in-degree (%d), out-degree (%d)" \
    #             %(key, in_and_out_degree[key]["in-degree"], in_and_out_degree[key]["out-degree"]))
    
    # # 1: in-degree (2), out-degree (3)
    # # 2: in-degree (3), out-degree (2)
    # # 3: in-degree (2), out-degree (2)
    # # 4: in-degree (2), out-degree (3)
    # # 5: in-degree (3), out-degree (2)
    # # so the minimum number of edges to be added is 2 (2 -> 1 and 5 -> 4)

    # print("\n###### problem 3 ######")
    # paired_composition = "(ACC|ATA)  (ACT|ATT)  (ATA|TGA)  (ATT|TGA)  (CAC|GAT)  (CCG|TAC)  (CGA|ACT)  (CTG|AGC)  (CTG|TTC)  (GAA|CTT)  (GAT|CTG)  (GAT|CTG)  (TAC|GAT)  (TCT|AAG)  (TGA|GCT)  (TGA|TCT)  (TTC|GAA)" \
    #                     .replace("  ", " ") \
    #                     .replace("(", "") \
    #                     .replace(")", "")
    # k = 3
    # d = 1
    # print(string_construction_by_paired_kmers_composition(paired_composition, k, d))
