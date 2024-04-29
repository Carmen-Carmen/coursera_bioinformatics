import random
import time
import pyperclip
import matplotlib.pyplot as plt

def print_arr(arr, is_signed = False):
    result = ""
    if is_signed:
        for item in arr:
            if item > 0:
                result += "+" + str(item) + " "
            else:
                result += str(item) + " "
    else:
        for item in arr:
            result += str(item) + " "
    result = result.strip()
    pyperclip.copy(result)
    
    print(result)

# Input: A chromosome Chromosome containing n synteny blocks.
# Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.
def transfer_chromosome_to_cycle(chromosome):
    nodes = [0 for _ in range(2 * len(chromosome))]
    for i in range(len(chromosome)):
        val = chromosome[i]
        # if i = 0, then the adjacent index in list nodes should be 0 and 1
        if val > 0:
            nodes[2 * i] = 2 * val - 1
            nodes[2 * i + 1] = 2 * val
        if val < 0:
            nodes[2 * i] = -(2 * val)
            nodes[2 * i + 1] = -(2 * val) - 1

    return nodes

# Input: A sequence Nodes of integers between 1 and 2n.
# Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.
def transfer_cycle_to_chromosome(nodes):
    chromosome = [0 for _ in range(len(nodes) // 2)]
    for i in range(len(nodes) // 2):
        if nodes[2 * i] < nodes[2 * i + 1]:
            chromosome[i] = nodes[2 * i + 1] // 2
        else:
            chromosome[i] = - nodes[2 * i] // 2
    
    return chromosome

# Input: A genome P. a genome may consist of multiple chromosomes
# Output: The collection of colored edges in the genome graph of P in the form (x, y).
def get_colored_edges(genome):
    edges = []
    for chromosome in genome:
        nodes = transfer_chromosome_to_cycle(chromosome)
        # assume that an n-element array (a1, . . . , an) has an 
        # invisible (n + 1)-th element that is equal to its first element
        nodes.append(nodes[0])
        for i in range(len(chromosome)):
            starting_node = nodes[2 * i + 1]
            ending_node = nodes[2 * i + 2]
            edges.append((starting_node, ending_node))
    
    return edges

# Input: The colored edges ColoredEdges of a genome graph.
# Output: The genome P corresponding to this genome graph.
def transfer_graph_to_genome(graph_edges):
    # this algorithm below should be applied to sorted edges
    genome = []
    cycle_starting_node = graph_edges[0][0]
    cycle_nodes = [graph_edges[0][0], graph_edges[0][1]]
    for edge in graph_edges[1 :]:
        # check if a cycle has been formed in cycle nodes
        if abs(cycle_nodes[-1] - cycle_starting_node) == 1 and \
        abs(edge[0] - cycle_nodes[-1]) != 1:
            # rearrange the order of nodes
            cycle_nodes = [cycle_nodes[-1]] + cycle_nodes[: -1]

            # generate the cyclic chromosome
            chromosome = transfer_cycle_to_chromosome(cycle_nodes)
            genome.append(chromosome)
            
            # update the starting node of a new cycle
            cycle_starting_node = edge[0]
            cycle_nodes = [edge[0], edge[1]]
            continue
        
        cycle_nodes.append(edge[0])
        cycle_nodes.append(edge[1])
    
    cycle_nodes = [cycle_nodes[-1]] + cycle_nodes[: -1]
    chromosome = transfer_cycle_to_chromosome(cycle_nodes)
    genome.append(chromosome)

    return genome

def find_edge_by_node(edges, node):
    for edge in edges:
        if node in edge:
            return edge
    
    return (-1, -1)

# given an edge and a node, return the node on the other side
def get_node_on_the_other_side(edge, node):
    if not (node in edge):
        return -1

    if node == edge[0]:
        return edge[1]
    else:
        return edge[0]

def parse_arr(str_list):
    return [int(item) for item in str_list.split(" ")]

def parse_genome(genome_str):
    genome = genome_str[1 : -1].split(")(")
    for i in range(len(genome)):
        genome[i] = parse_arr(genome[i])
    
    return genome

def parse_graph_edges(edges_str): 
    edges = edges_str[1 : -1].split("), (")
    for i in range(len(edges)):
        edge = [int(item) for item in edges[i].split(", ")]
        edges[i] = (edge[0], edge[1])
    
    return edges

def print_genome(genome):
    output = ""
    for chromosome in genome:
        temp = "("
        for synteny_block in chromosome:
            if synteny_block > 0:
                temp += "+" + str(synteny_block) + " "
            else:
                temp += str(synteny_block) + " "
        temp = temp.strip()
        temp += ")"
        output += temp
    
    pyperclip.copy(output)
    print(output)
    return output

# Code Challenge: Solve the 2-Break Distance Problem.
# Input: Genomes P and Q.
# Output: The 2-break distance d(P, Q).
def get_2_break_distance(genome_P, genome_Q):
    graph_P = get_colored_edges(genome_P)
    graph_Q = get_colored_edges(genome_Q)
    # print_arr(graph_P)
    # print_arr(graph_Q)
    graph = graph_P + graph_Q
    num_block = len(graph_P)

    num_cycle = 0
    cycle_start_node = graph[0][0]
    prev_node = graph[0][1]
    graph.remove(graph[0])
    to_update_cycle_start = False
    while len(graph) != 0:
        # update the starting node of a new cycle
        if to_update_cycle_start:
            cycle_start_node = graph[0][0]
            # also remember to update the prev_node cursor
            prev_node = graph[0][1]
            graph.remove(graph[0])
            to_update_cycle_start = False
            continue

        for edge in graph:
            # the edge is undirected
            if prev_node in edge:
                graph.remove(edge)
                if prev_node == edge[0]:
                    prev_node = edge[1]
                else:
                    prev_node = edge[0]
                break
        
        # check that if a cycle is formed
        if prev_node == cycle_start_node:
            num_cycle += 1
            to_update_cycle_start = True

    # 2-Break Distance Theorem: 
    # The 2-break distance between genomes P and Q is equal to Blocks(P, Q)− Cycles(P, Q)
    return num_block - num_cycle

# Input: The colored edges of a genome graph GenomeGraph, 
#       followed by indices i1 , i2 , i3 , and i4 .
# Output: The colored edges of the genome graph 
#       resulting from applying the 2-break operation 
#       2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 ).
def do_2_break_on_graph(graph_edges, start1, end1, start2, end2):
    # i.e. break edges: (start1, end1), (start2, end2)
    # and form edges: (start1, start2), (end1, end2)
    graph_edges = [
        edge for edge in graph_edges \
        if not((start1 in edge) and (end1 in edge)) and \
        not((start2 in edge) and (end2 in edge))        
    ]

    graph_edges.append((start1, start2))
    graph_edges.append((end1, end2))

    return graph_edges

# Input: A genome P, followed by indices i1 , i2 , i3 , and i4 .
# Output: The genome P' resulting from applying the 2-break operation 2-BreakOnGenome(GenomeGraph i1 , i2 , i3 , i4 ).
def do_2_break_on_genome(genome, start1, end1, start2, end2):
    graph_edges = get_colored_edges(genome)
    graph_edges = do_2_break_on_graph(graph_edges, start1, end1, start2, end2)

    # generate black edges as a reference
    black_edges = []
    block_num = 0
    for chromosome in genome:
        for block in chromosome:
            block_num += 1
            if block > 0:
                edge = (
                    2 * block - 1, 
                    2 * block
                )
                black_edges.append(edge)
            else:
                edge = (
                    2 * abs(block), 
                    2 * abs(block) - 1
                )
                black_edges.append(edge)
    
    # sort colored edges according to the black edges
    sorted_edges = []
    edge = graph_edges.pop(0)
    black_edge = find_edge_by_node(black_edges, edge[1])
    black_edges.remove(black_edge)
    if edge[1] != black_edge[0]:
        edge = (
                edge[1], 
                edge[0]
            )
    sorted_edges.append(edge)
    while len(sorted_edges) != block_num:
        edge = find_edge_by_node(graph_edges, black_edge[1])
        if edge == (-1, -1):
        # end of a cycle
            edge = graph_edges.pop(0)
            black_edge = find_edge_by_node(black_edges, edge[1])
            black_edges.remove(black_edge)
            if edge[1] != black_edge[0]:
                edge = (
                        edge[1], 
                        edge[0]
                    )
            sorted_edges.append(edge)
        else:
        # extension of the current cycle
            graph_edges.remove(edge)
            if edge[0] != black_edge[1]:
                edge = (
                    edge[1], 
                    edge[0]
                )
            sorted_edges.append(edge)
            black_edge = find_edge_by_node(black_edges, edge[1])
            black_edges.remove(black_edge)

            if black_edge[0] != edge[1]:
                black_edge = (
                    black_edge[1], 
                    black_edge[0]
                )

    return transfer_graph_to_genome(sorted_edges)

# two identical cyclic genome can be presented in different ways
def compare_cyclic_genome(target_genome, ref_genome):
    if len(target_genome) != len(ref_genome):
        return False
    
    for i in range(len(target_genome)):
        target_chm = target_genome[i]
        ref_chm = ref_genome[i]
        if len(target_chm) != len(ref_chm) :
            return False
        
        target_first_block = target_chm[0]
        tfb_index_in_ref = ref_chm.index(target_first_block)
        ref_chm = ref_chm[tfb_index_in_ref :] + ref_chm[: tfb_index_in_ref]

        if target_chm != ref_chm:
            return False
    
    return True

# ShortestRearrangementScenario(P, Q)
#      output P
#      RedEdges ← ColoredEdges(P)
#      BlueEdges ← ColoredEdges(Q)
#      BreakpointGraph ← the graph formed by RedEdges and BlueEdges
#      while BreakpointGraph has a non-trivial cycle Cycle
#           (i2,i3)<-An arbitrary edge from BlueEdges in a non trivial red-blue cycle
#           (i1,i2)<-An edge from RedEdges originating at node i1
#           (i3,i4)<-an edge from RedEdges originating at node i3
#           RedEdges ← RedEdges with edges (i1, i2) and (i3, i4) removed
#           RedEdges ← RedEdges with edges (i2, i3) and (i4, i1) added
#           BreakpointGraph ← the graph formed by RedEdges and BlueEdges
#           P ← 2-BreakOnGenome(P, i1 , i3 , i2 , i4 )
#           output P
# 2-Break Sorting Problem: Find a shortest transformation of one genome into another by 2-breaks.
# Input: Two genomes with circular chromosomes on the same set of synteny blocks.
# Output: The sequence of genomes resulting from applying a shortest 
#         sequence of 2-breaks transforming one genome into the other.
def do_2_break_sorting(genome_P, genome_Q):
    print_genome(genome_P)
    # print_genome(genome_Q)
    # print()

    distance = get_2_break_distance(genome_P, genome_Q)
    while distance != 0:
        graph_P = get_colored_edges(genome_P)
        graph_Q = get_colored_edges(genome_Q)

        # print(graph_P)
        # print(graph_Q)

        # select an arbitrary edge from BlueEdges in a non trivial red-blue cycle
        # i.e. an edge that appears in Q but not in P
        selected_edge = random.choice(
            [edge for edge in graph_Q if not(edge in graph_P)]
        )

        end1 = selected_edge[0]
        end2 = selected_edge[1]
        start1 = get_node_on_the_other_side(find_edge_by_node(graph_P, end1), end1)
        start2 = get_node_on_the_other_side(find_edge_by_node(graph_P, end2), end2)
        
        # print(selected_edge)
        # print((i1, i2))
        # print((i3, i4))
        # print()

        # update the smallest 2-break distance 
        # if the 2-break performed this time reduced the 2-distance
        prev_genome_P = genome_P.copy()
        genome_P = do_2_break_on_genome(genome_P, start1, end1, start2, end2)
        if get_2_break_distance(genome_P, genome_Q) < distance:
            print_genome(genome_P)
            # print_genome(genome_Q)
            # print()
            distance = get_2_break_distance(genome_P, genome_Q)
        else:
            # if 2-break distanced increased, 
            # genome P should be returned to the last state
            genome_P = prev_genome_P
            continue

def get_shared_k_mers(strand1, strand2, k):
    min_len = min(
        len(strand1), 
        len(strand2)
    )
    if k > min_len:
        return []
    
    results = []

    # O(n2), too slow
    # for i in range(0, len(strand1) - k + 1):
    #     k_mer1 = strand1[i : i + k]
    #     for j in range(0, len(strand2) - k + 1):
    #         k_mer2 = strand2[j : j + k]

    #         if k_mer1 == k_mer2 or \
    #         k_mer1 == get_reverse_complement_DNA_strand(k_mer2):
    #             results.append(
    #                 (i, j)
    #             )

    # try use hashing to store k-mers in strand2
    # O(n) time complexity
    k_mer2_dict = {}
    for j in range(0, len(strand2) - k + 1):
        k_mer2 = strand2[j : j + k]
        if k_mer2 in k_mer2_dict.keys():
            k_mer2_dict[k_mer2].append(j)
        else:
            k_mer2_dict[k_mer2] = [j]
    for i in range(0, len(strand1) - k + 1):
        k_mer1 = strand1[i : i + k]
        k_mer1_rev = get_reverse_complement_DNA_strand(k_mer1)
        if k_mer1 in k_mer2_dict.keys():
            for j_pos in k_mer2_dict[k_mer1]:
                results.append((i, j_pos))
        if k_mer1_rev in k_mer2_dict.keys():
            for j_pos in k_mer2_dict[k_mer1_rev]:
                results.append((i, j_pos))
    
    return results

def get_reverse_complement_DNA_strand(s):
    l = len(s)
    result = ""
    for i in range(0, l):
        base = s[i: i + 1]
        if base == "A":
            result += "T"
        elif base == "T":
            result += "A"
        elif base == "C":
            result += "G"
        elif base == "G":
            result += "C"
    
    reversed_result = ""
    for c in reversed(result):
        reversed_result += c
    
    return reversed_result

if __name__ == "__main__":
    # 1.5 step 4:
    # chromosome = parse_arr("+1 -2 -3 +4")
    # print_arr(transfer_chromosome_to_cycle(chromosome))
    # with open("./dataset_8222_4.txt", "r") as f:
    #     chromosome = parse_arr(f.read().strip()[1 : -1])
    #     print_arr(transfer_chromosome_to_cycle(chromosome))

    # 1.5 step 5:
    # nodes = parse_arr("12 11 10 9")
    # print_arr(transfer_cycle_to_chromosome(nodes), is_signed=True)
    # with open("./dataset_8222_5.txt", "r") as f:
    #     nodes = parse_arr(f.read().strip()[1 : -1])
    #     print_arr(transfer_cycle_to_chromosome(nodes), is_signed=True)

    # 1.5 step 7:
    # genome_str = "(+1 -2 -3)(+4 +5 -6)"
    # print(parse_genome(genome_str))
    # print(get_colored_edges(parse_genome(genome_str)))
    # with open("./dataset_8222_7.txt", "r") as f:
    #     genome_str = f.read().strip()
    #     print(get_colored_edges(parse_genome(genome_str)))

    # 1.5 step 8:
    # edges_str = "(2, 4), (3, 6), (5, 1), (7, 9), (10, 12), (11, 8)"
    # graph_edges = parse_graph_edges(edges_str)
    # print(transfer_graph_to_genome(graph_edges))

    # with open("./GraphToGenome.txt", "r") as f:
    #     f.readline()
    #     edges_str = f.readline().strip()
    #     graph_edges = parse_graph_edges(edges_str)
    #     genome = print_genome(transfer_graph_to_genome(graph_edges))
    #     print(genome)
    #     f.readline()
    #     if f.readline().strip() == genome:
    #         print("correct!")
    #     else:
    #         print("wrong!")

    # with open("./dataset_8222_8.txt", "r") as f:
    #     edges_str = f.readline().strip()
    #     graph_edges = parse_graph_edges(edges_str)
    #     genome = transfer_graph_to_genome(graph_edges)
    #     print(genome)

    # 1.2 step 4:
    # genome_P = parse_genome("(+1 +2 +3 +4 +5 +6)")
    # genome_Q = parse_genome("(+1 -3 -6 -5)(+2 -4)")
    # print(get_2_break_distance(genome_P, genome_Q))

    # with open("./2BreakDistance.txt", "r") as f:
    #     f.readline()
    #     genome_P = parse_genome(f.readline().strip())
    #     genome_Q = parse_genome(f.readline().strip())
    #     distance_2_break = get_2_break_distance(genome_P, genome_Q)
    #     print(distance_2_break)
    #     f.readline()
    #     correct_answer = int(f.readline().strip())
    #     if distance_2_break == correct_answer:
    #         print("correct!")
    #     else:
    #         print("wrong!")

    # with open("./dataset_288_4.txt", "r") as f:
    #     genome_P = parse_genome(f.readline().strip())
    #     genome_Q = parse_genome(f.readline().strip())
    #     print(get_2_break_distance(genome_P, genome_Q))

    # 1.6 step 2:
    # graph_edges = parse_graph_edges("(2, 4), (3, 8), (7, 5), (6, 1)")
    # indices = [1, 6, 3, 8]
    # start1 = indices[0]
    # end1 = indices[1]
    # start2 = indices[2]
    # end2 = indices[3]
    # print(do_2_break_on_graph(graph_edges, start1, end1, start2, end2))
    # with open("./dataset_8224_2.txt", "r") as f:
    #     graph_edges = parse_graph_edges(f.readline().strip())
    #     indices = [int(item) for item in f.readline().strip().split(", ")]
    #     start1 = indices[0]
    #     end1 = indices[1]
    #     start2 = indices[2]
    #     end2 = indices[3]
    #     print(do_2_break_on_graph(graph_edges, start1, end1, start2, end2))

    # 1.6 step 3:
    # genome = parse_genome("(+1 -3 -6 -5)(+2 -4)")
    # indices = [5, 12, 3, 7]
    # genome = parse_genome("(+1 -2 -4 +3)")
    # indices = [1, 6, 3, 8]
    # start1 = indices[0]
    # end1 = indices[1]
    # start2 = indices[2]
    # end2 = indices[3]
    # print_genome(do_2_break_on_genome(genome, start1, end1, start2, end2))

    # with open("./2BreakOnGenome.txt", "r") as f:
    #     f.readline()
    #     genome = parse_genome(f.readline().strip())
    #     indices = [int(item) for item in f.readline().strip().split(", ")]
    #     start1 = indices[0]
    #     end1 = indices[1]
    #     start2 = indices[2]
    #     end2 = indices[3]
    #     genome = print_genome(do_2_break_on_genome(genome, start1, end1, start2, end2))
    #     f.readline()
    #     correct_answer = f.readline().strip()
    #     target_genome = parse_genome(genome)
    #     ref_genome = parse_genome(correct_answer)
    #     if compare_cyclic_genome(target_genome, ref_genome):
    #         print("correct!")
    #     else:
    #         print("wrong!")

    # with open("./dataset_8224_3.txt", "r") as f:
    #     genome = parse_genome(f.readline().strip())
    #     indices = [int(item) for item in f.readline().strip().split(", ")]
    #     start1 = indices[0]
    #     end1 = indices[1]
    #     start2 = indices[2]
    #     end2 = indices[3]
    #     print_genome(do_2_break_on_genome(genome, start1, end1, start2, end2))

    # 1.2 step 5:
    # genome_P = parse_genome("(+1 -2 -3 +4)")
    # genome_Q = parse_genome("(+1 +2 -4 -3)")
    # genome_P = parse_genome("(+1 -7 +6 -10 +9 -8 +2 -11 -3 +5 +4)")
    # genome_Q = parse_genome("(+1 +2 +3 +4 +5 +6 +7 +8 +9 +10 +11)")
    # do_2_break_sorting(genome_P, genome_Q)

    # with open("./dataset_288_5.txt", "r") as f:
    #     genome_P = parse_genome(f.readline().strip())
    #     genome_Q = parse_genome(f.readline().strip())
    #     do_2_break_sorting(genome_P, genome_Q)

    # 1.4 step 3:
    # strand1 = "AAACTCATC"
    # strand2 = "TTTCAAATC"
    # k = 2
    # print("num of shared k-mers:", len(get_shared_k_mers(strand1, strand2, k)))

    # 1.4 step 5:
    # with open("./dataset_289_5.txt", "r") as f:
    #     k = int(f.readline().strip())
    #     strand1 = f.readline().strip()
    #     strand2 = f.readline().strip()
    #     shared_k_mers = get_shared_k_mers(strand1, strand2, k)
    #     to_copy = ""
    #     for pos in shared_k_mers:
    #         to_copy += str(pos) + "\n"
    #         print(pos)
    #     to_copy = to_copy.strip()
    #     pyperclip.copy(to_copy)
    #     x = [i[0] for i in shared_k_mers]
    #     y = [i[1] for i in shared_k_mers]
    #     plt.scatter(x, y, s = 5)
    #     plt.show()

    # 1.4 step 7: 
    # f1 = open("./E_coli.txt", "r")
    # f2 = open("./Salmonella_enterica.txt", "r")
    # strand1 = f1.read().strip()
    # strand2 = f2.read().strip()
    # k = 30
    # shared_k_mers = get_shared_k_mers(strand1, strand2, k)
    # print("shared 30-mers from E. coli and S. enterica genomes: %d" %(len(shared_k_mers)))
    # x = [i[0] for i in shared_k_mers]
    # y = [i[1] for i in shared_k_mers]
    # plt.scatter(x, y, s = 5)
    # plt.show()

    # try:
    #     f1.close()
    #     f2.close()
    # except Exception:
    #     pass

    # week 5 Quiz
    # problem 4
    print("\n### problem 4 ###")
    strand1 = "TCTTGCAGCTCGTCA"
    strand2 = "GTACTTTCAGAATCA"
    k = 3
    print(len(get_shared_k_mers(strand1, strand2, k)))

