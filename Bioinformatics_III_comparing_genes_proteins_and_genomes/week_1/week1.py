import math
import pyperclip
import os

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    result = result.strip()
    pyperclip.copy(result)
    
    print(result)
    

def combination(n, r):
    numerator = math.factorial(n)
    denominator = math.factorial(r) * math.factorial(n - r)
    combination = numerator // denominator
    return combination

###### example of recursive ######
def Hanoi_towers(n, start_peg, dest_peg):
    if n == 1:
        print("%d --> %d" %(start_peg, dest_peg))
        return
    
    transit_peg = 6 - start_peg - dest_peg

    Hanoi_towers(n - 1, start_peg, transit_peg)
    print("%d --> %d" %(start_peg, dest_peg))
    Hanoi_towers(n - 1, transit_peg, dest_peg)

    return

###### dynamic programming learning ######
# take the Change Problem as an example
# Change Problem: Find the minimum number of coins needed to make change.
# Input: An integer money and an array Coins of d positive integers.
# Output: The minimum number of coins with denominations Coins that changes money.

# implement recursive in solving this problem
# the method will be called too many times, thus impractical
def recursive_change_problem(money, coins):
    if money == 0:
        return 0
    
    # min_coin_num = infinite big
    min_coin_num = money // min(coins)
    # traverse all the coin values
    for coin_val in coins:
        if money >= coin_val:
            # min_coin_num = min(
            #                       min_coin_num(money - coins[1])) + 1, 
            #                       min_coin_num(money - coins[2])) + 1, 
            #                       ...
            #                       min_coin_num(money - coins[d])) + 1
            #                   )
            coin_num = recursive_change_problem(money - coin_val, coins)
            if coin_num + 1 < min_coin_num:
                min_coin_num = coin_num + 1
    
    return min_coin_num

def dp_change_problem(money, coins):
    # step 1: return only the minimum number of coins
    # # initialization dynamic programming table
    # max_num = money // min(coins)
    # table = [max_num for _ in range(money + 1)]
    # # table[money] = min_coin_num
    # table[0] = 0

    # for m in range(1, money + 1):
    #     for coin_val in coins:
    #         if m >= coin_val:
    #             if table[m - coin_val] + 1 < table[m]:
    #                 table[m] = table[m - coin_val] + 1
    
    # return table[money]

    # step 2: also return the coins needed for change
    # let the length of dp_table not exceed the value of the largest coin
    max_num = money // min(coins)
    max_coin_val = max(coins)
    dp_table = {}
    dp_table[0] = (0, [])

    for m in range(1, money + 1):
        # print("len of the dp table:", len(dp_table))
        # print(dp_table)

        # the value to the key is a list
        # the 1st element is the minimum number of coins
        # the 2nd element is a array of the coins needed
        dp_table[m] = [max_num, []]

        for coin_val in coins:
            if m >= coin_val:
                min_num_and_coins_prev = dp_table[m - coin_val]
                min_num_and_coins_current = dp_table[m]
                if min_num_and_coins_prev[0] + 1 < min_num_and_coins_current[0]:
                    # refresh min_coin_num
                    min_num_and_coins_current[0] = min_num_and_coins_prev[0] + 1
                    min_num_and_coins_current[1] = []
                    # append all the coins in the previous list and the current coin_val to the current list
                    for coin in min_num_and_coins_prev[1]:
                        min_num_and_coins_current[1].append(coin)
                    min_num_and_coins_current[1].append(coin_val)
            
        # items in the dp_table with keys smaller than (m - max_coin_val) are no longer needed in further dynamic programming
        # if (m - max_coin_val) in dp_table.keys():
        if m >= max_coin_val:
            # print("m:", m)
            # e.g. max_coin_val = 50, m = 50
            # then in the next cycle, m = 51, the smallest value of (m - coin_val) = 51 - 50 = 1
            # so it is okay to remove 0 (i.e. 50 - 50) in this cycle since dp_table[0] is no longer needed
            dp_table.pop(m - max_coin_val)
    
    # the return value of this method is a list
    return dp_table[money]

###### Manhattan Tourist problem ######
# ManhattanTourist(n, m, Down, Right)
#     s0, 0 ← 0
#     for i ← 1 to n
#         si, 0 ← si-1, 0 + downi-1, 0
#     for j ← 1 to m
#         s0, j ← s0, j−1 + right0, j-1
#     for i ← 1 to n
#         for j ← 1 to m
#             si, j ← max{si - 1, j + downi-1, j, si, j - 1 + righti, j-1}
#     return sn, m
# input: integers n and m, 
#       n * (m + 1) matrix Down
#       (n + 1) * m matrix Right
# output: the length of a longest path from (0, 0) to (n, m
def dp_manhattan_tourist(n, m, Down_weights, Right_weights):
    dp_table = [
        [
            0 for _ in range(m + 1)
        ] for _ in range(n + 1)
    ]
    # print(dp_table)
    dp_table[0][0] = 0

    # calculate path length for vertices in the first row
    for j in range(1, m + 1):
        dp_table[0][j] = dp_table[0][j - 1] + Right_weights[0][j - 1]
    # calculate path length for vertices in the first column
    for i in range(1, n + 1):
        dp_table[i][0] = dp_table[i - 1][0] + Down_weights[i - 1][0]
    # calculate path length for other vertices
    for j in range(1, m + 1):
        for i in range(1, n + 1):
            dp_table[i][j] = max(
                dp_table[i - 1][j] + Down_weights[i - 1][j], 
                dp_table[i][j - 1] + Right_weights[i][j - 1]
            )
    
    return dp_table[n][m]

# input: two strings
# output: the backtracking matrix containing pointers indicating 
# deletions(↓), insertions(→) or matches/mismatches(↘︎)
def get_LCS_backtrack_matrix(str1, str2):
    # row_num = len(str1), col_num = len(str2)
    backtrack_matrix = [
        [
            "" for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    # dp_table stores the length of an LCS to each node
    dp_table = [
        [
            0 for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            is_matched = 0
            if str1[i - 1] == str2[j - 1]:
                is_matched = 1
            dp_table[i][j] = max(
                dp_table[i - 1][j], 
                dp_table[i][j - 1], 
                dp_table[i - 1][j - 1] + is_matched
            )

            if dp_table[i][j] == dp_table[i - 1][j]:
                # deletion
                backtrack_matrix[i][j] = "↓"
            elif dp_table[i][j] == dp_table[i][j - 1]:
                # insertion
                backtrack_matrix[i][j] = "→"
            elif dp_table[i][j] == dp_table[i - 1][j - 1] + is_matched:
                # match/mismatch
                backtrack_matrix[i][j] = "↘︎"
            
    return backtrack_matrix 

# input: 
#       - backtrack matrix generated from str1 and str2
#       - str1, used to generate the path
#       - i, j: refer to the row_num and col_num in the matrix, respectively
# recursively output each node in the path
def get_LCS_path(backtrack_matrix, str1, i, j):
    # method 1: recursive
    # if i == 0 or j == 0:
    #     # when i == 0 or j == 0, backtrack_matrix[i][j] is ""
    #     return ""
    # if backtrack_matrix[i][j] == "↓":
    #     return get_LCS_path(backtrack_matrix, str1, i - 1, j)
    # elif backtrack_matrix[i][j] == "→":
    #     return get_LCS_path(backtrack_matrix, str1, i, j - 1)
    # else:
    #     return get_LCS_path(backtrack_matrix, str1, i - 1, j - 1) + str1[i - 1]

    # method 2: iterative, avoiding the problem of insufficient recursion depth
    LCS_path = ""
    while i > 0 and j > 0:
        if backtrack_matrix[i][j] == "↓":
            i -= 1
        elif backtrack_matrix[i][j] == "→":
            j -= 1
        else:
            LCS_path += str1[i - 1]
            i -= 1
            j -= 1

    LCS_path = "".join(reversed(LCS_path))
    return LCS_path

def visulization_backtrack_matrix(str1, str2, backtrack_matrix):
    for_print = "  "
    for j in range(len(str2)):
        for_print += str2[j] + " "
    print(for_print)
    for i in range(len(str1)):
        line = str1[i]
        for j in range(len(str2)):
            # backtrack_matrix[0][j] = backtrack_matrix[i][0] = "", should not be shown
            line += " " + backtrack_matrix[i + 1][j + 1]
        print(line)

def parse_digital_two_dimension_array(text):
    arr = text.split("\n")
    two_dimension_arr = []
    for line in arr:
        two_dimension_arr.append([int(digit) for digit in line.split(" ")])
    
    return two_dimension_arr

# Input: An integer representing the starting node to consider in a graph, 
#       followed by an integer representing the ending node to consider, 
#       followed by a list of edges in the graph. 
#       The edge notation "0 1 7" indicates that an edge connects node 0 to node 1 with weight 7.  
#       You may assume a given topological order corresponding to nodes in increasing order.
# Output: The length of a longest path in the graph, followed by a longest path as 
#       a sequence of space-separated node labels. (If multiple longest paths exist, you may return any one.)
def get_longest_path_from_DAG(topological_ordered_DAG):
    # construct a list of edges
    # each edge is represented by ("start->end", weight)
    # edges are stored in the given topological order
    in_and_out_degrees = {}
    edges = []
    for edge in topological_ordered_DAG.split("\n")[1: ]:
        splitted = edge.split(" ")
        start_node = splitted[0]
        end_node = splitted[1]
        weight = 0
        if len(splitted) != 2:
            weight = int(splitted[2])
        # edge
        edges.append(("%s->%s" %(start_node, end_node), weight))
        # in-degree and out-degree
        if not(start_node in in_and_out_degrees.keys()):
            in_and_out_degrees[start_node] = {
                "in-degree": 0, 
                "out-degree": 0, 
            }
        in_and_out_degrees[start_node]["out-degree"] += 1
        if not(end_node in in_and_out_degrees.keys()):
            in_and_out_degrees[end_node] = {
                "in-degree": 0, 
                "out-degree": 0, 
            }
        in_and_out_degrees[end_node]["in-degree"] += 1

    # source被规定为第一行的第一个node, sink被规定为第二个
    source = topological_ordered_DAG.split("\n")[0].split(" ")[0]
    sink = topological_ordered_DAG.split("\n")[0].split(" ")[1]
    to_remove = []
    for node in in_and_out_degrees.keys():
        if in_and_out_degrees[node]["in-degree"] == 0 and node != source:
            to_remove.append(node)
        if in_and_out_degrees[node]["out-degree"] == 0 and node != sink:
            to_remove.append(node)
    for node in to_remove:
        in_and_out_degrees.pop(node)
    
    # print("starting: %s, ending: %s" %(source, sink))

    nodes = list(in_and_out_degrees.keys())
    nodes = sorted(nodes, key=lambda x: int(x))
    # nodes = sorted(nodes)

    # dp_table = {source: [0, [source]]}
    # # topological order corresponding to nodes in increasing order
    # for node in [item for item in nodes if item != source]:
    #     # val: (path_len, path)
    #     dp_table[node] = [0, []]
    #     # traverse all edges with pointing to current node
    #     for edge in [item for item in edges if item[0].split("->")[1] == node]:
    #         prev_node = edge[0].split("->")[0]

    #         #  1. avoid key-error         2. make sure source node is in the path
    #         if not(prev_node in nodes) or not(source in dp_table[prev_node][1]):
    #             continue

    #         weight = edge[1]
    #         new_length = dp_table[prev_node][0] + weight
    #         if new_length >= dp_table[node][0]:
    #             dp_table[node][0] = new_length
    #             dp_table[node][1] = dp_table[prev_node][1].copy()
    #             dp_table[node][1].append(node)

    # return dp_table[sink]
    dp_table = {source: 0}
    backtrack_table = {source: None}
    for node in [item for item in nodes if item != source]:
        dp_table[node] = 0
        backtrack_table[node] = None
        for edge in [item for item in edges if item[0].split("->")[1] == node]:
            prev_node = edge[0].split("->")[0]
            if not(prev_node in nodes) or \
                (prev_node != source and not(source in backtrack_to_path(backtrack_table, source, prev_node))):
                continue
            weight = edge[1]
            new_length = dp_table[prev_node] + weight
            if new_length >= dp_table[node]:
                dp_table[node] = new_length
                backtrack_table[node] = prev_node
    
    return dp_table[sink], backtrack_to_path(backtrack_table, source, sink)

# backtrack to get the longest path
def backtrack_to_path(backtrack_table, source, sink):
    longest_path = []
    node = sink
    longest_path.append(node)
    while True and node != None:
        node = backtrack_table[node]
        longest_path.append(node)
        if node == source:
            break
    
    longest_path = [item for item in reversed(longest_path)]
    return longest_path


if __name__ == "__main__":
    # 1.3 step 4
    # n_col = 16
    # n_row = 12
    # print("number of paths:", combination(n_col + n_row, n_row))

    # 1.9 step 4
    # Hanoi_towers(10, 1, 3)

    # 1.5 step 6
    # coins = [5, 4, 1]
    # print("number of coins: %d" %recursive_change_problem(10, coins))

    # 1.5 step 8
    # coins = [5, 4, 1]
    # table = []
    # for money in range(13, 23):
    #     table.append(dp_change_problem(money, coins))
    
    # print_arr(table)

    # 1.5 step 10
    # sample
    # money = 40
    # coins = [int(coin_val) for coin_val in "50 25 20 10 5 1".split(" ")]
    # print("sample:", dp_change_problem(money, coins)[0])
    # randomized dataset
    # with open("./dataset_243_10.txt", "r") as f:
    #     money = int(f.readline().strip())
    #     coins = [int(coin_val) for coin_val in f.readline().strip().split(" ")]
    #     print("randomized dataset:", dp_change_problem(money, coins)[0])
    #     print_arr(dp_change_problem(money, coins)[1])

    # 1.6 step 10: 
    # Down_weights = parse_digital_two_dimension_array(
    #     "1 0 2 4 3\n4 6 5 2 1\n4 4 5 2 1\n5 6 8 5 3"
    # )
    # Right_weights = parse_digital_two_dimension_array(
    #     "3 2 4 0\n3 2 4 2\n0 7 3 3\n3 3 0 2\n1 3 2 2"
    # )
    # print("sample:", dp_manhattan_tourist(4, 4, Down_weights, Right_weights))

    # with open("./dataset_261_10.txt", "r") as f:
    #     args = f.readline().strip().split(" ")
    #     n = int(args[0])
    #     m = int(args[1])
    #     matrices = f.read().strip().split("-")
    #     Down_weights = parse_digital_two_dimension_array(matrices[0].strip())
    #     Right_weights = parse_digital_two_dimension_array(matrices[1].strip())
    #     print("randomized dataset:", dp_manhattan_tourist(n, m, Down_weights, Right_weights))

    # 1.8 step 3 & 4
    # backtrack_matrix = get_LCS_backtrack_matrix("HUMAN", "CHIMPANZEE")
    # for line in backtrack_matrix:
    #     print_arr(line)
    # print(get_LCS_path(backtrack_matrix, "HUMAN", len("HUMAN") - 1, len("CHIMPANZEE") - 1))

    # 1.8 step 5
    # str1 = "AACCTTGG"
    # str2 = "ACACTGTGA"
    # backtrack_matrix = get_LCS_backtrack_matrix(str1, str2)
    # visulization_backtrack_matrix(str1, str2, backtrack_matrix)
    # print("sample:", get_LCS_path(backtrack_matrix, str1, len(str1), len(str2)))

    # debug datasets
    # os.chdir("./LongestCommonSubsequence")
    # inputs = os.listdir("./inputs")
    # outputs = os.listdir("./outputs")
    # try:
    #     inputs.remove(".DS_Store")
    #     outputs.remove(".DS_Store")
    # except Exception:
    #     pass
    # inputs.sort()
    # outputs.sort()
    # counter = 0
    # for input in inputs:
    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()
    #         len1 = len(str1)
    #         len2 = len(str2)
    #         backtrack_matrix = get_LCS_backtrack_matrix(str1, str2)
    #         result = get_LCS_path(backtrack_matrix, str1, len1, len2)

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.readline().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")
        
    #     counter += 1
    #     print()

    # with open("./dataset_245_5.txt", "r") as f:
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     len1 = len(str1)
    #     len2 = len(str2)
    #     backtrack_matrix = get_LCS_backtrack_matrix(str1, str2)
    #     print(get_LCS_path(backtrack_matrix, str1, len1, len2))

    # 1.8 step 7
#     topological_ordered_DAG = """0 4
# 0 1 7
# 0 2 4
# 2 3 2
# 1 4 1
# 3 4 3"""
#     path_len = get_longest_path_from_DAG(topological_ordered_DAG)[0]
#     path = get_longest_path_from_DAG(topological_ordered_DAG)[1]
#     print("sample: ")
#     print(path_len)
#     print_arr(path)

    # debug dataset
    os.chdir("./LongestPath")
    inputs = os.listdir("./inputs")
    outputs = os.listdir("./outputs")
    try:
        inputs.remove(".DS_Store")
        outputs.remove(".DS_Store")
    except Exception:
        pass
    inputs.sort()
    outputs.sort()
    counter = 0
    for input in inputs:
        print("###### %s ######" %input)
        result = ""
        with open("./inputs/%s" %input, "r") as f_input:
            topological_ordered_DAG = f_input.read().strip()
            path_len = get_longest_path_from_DAG(topological_ordered_DAG)[0]
            path = get_longest_path_from_DAG(topological_ordered_DAG)[1]
            result += str(path_len) + "\n"
            for node in path:
                result += node + " "
            result = result[:-1]

        with open("./outputs/%s" %outputs[counter], "r") as f_output:
            output = f_output.read().strip()
            print("Result:", result)
            print("Correct:", output)
            if result == output:
                print("debug_%d" %(counter + 1), "correct!")
            else:
                print("debug_%d" %(counter + 1), "wrong!")
        
        counter += 1
        print()

    # with open("./dataset_245_7.txt", "r") as f:
    #     topological_ordered_DAG = f.read().strip()
    #     path_len = get_longest_path_from_DAG(topological_ordered_DAG)[0]
    #     path = get_longest_path_from_DAG(topological_ordered_DAG)[1]
    #     print("randomized dataset:")
    #     print(path_len)
    #     print_arr(path)
    
    # week1 Quiz
    # problem 1
#     print("###### problem 1 ######")
#     str1 = "GCGATC"
#     str2 = "CTGACG"
#     print(get_LCS_path(get_LCS_backtrack_matrix(str1, str2), str1, len(str1), len(str2)))

#     print("\n###### problem 3 ######")
#     # dynamic programming provided by ChatGPT
#     # def count_combinations(target_sum):
#     #     # Create a list to store the count of combinations for each sum up to the target sum
#     #     combinations = [0] * (target_sum + 1)

#     #     # Set the base cases
#     #     combinations[0] = 1

#     #     # Iterate through the numbers 2 and 3
#     #     for num in [2, 3]:
#     #         # Update the count of combinations for each sum
#     #         for i in range(num, target_sum + 1):
#     #             combinations[i] += combinations[i - num]

#     #     return combinations[target_sum]

#     # target_sum = 25
#     # num_combinations = count_combinations(target_sum)
#     # print(num_combinations)

#     def backtrack(remaining_sum, current_combination, combinations):
#         if remaining_sum == 0:
#             combinations.append(current_combination)
#         elif remaining_sum > 0:
#             # Try adding a 2 to the current combination
#             backtrack(remaining_sum - 2, current_combination + [2], combinations)
#             # Try adding a 3 to the current combination
#             backtrack(remaining_sum - 3, current_combination + [3], combinations)

#     def find_combinations(target_sum):
#         combinations = []
#         backtrack(target_sum, [], combinations)
#         return combinations

#     target_sum = 22
#     combinations = find_combinations(target_sum)
#     # for combination in combinations:
#     #     print((combination))
#     print(len(combinations))

#     print("\n###### problem 5 ######")
#     DAG = """a g
# a b 3
# a c 6
# a d 5
# b c 2
# b f 4
# c e 4
# c f 3
# c g 7
# d e 4
# d f 5
# e g 2
# f g 1"""
#     print(get_longest_path_from_DAG(DAG)[0])
#     print_arr(get_longest_path_from_DAG(DAG)[1])