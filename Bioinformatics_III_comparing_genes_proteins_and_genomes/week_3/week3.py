import math
import pyperclip
import os
import sys

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + "\t"
    result = result.strip()
    pyperclip.copy(result)
    
    print(result)

def print_2D_arr(arr):
    for line in arr:
        print_arr(line)
    print()

# Input: A match reward, 
#       a mismatch penalty, 
#       a gap opening penalty, 
#       a gap extension penalty, 
#       and two nucleotide strings.
# Output: The maximum alignment score between v and w, 
#       followed by an alignment of v and w achieving this maximum score.
def affine_gap_penalty_alignment(score_matrix, str1, str2):
    match = score_matrix[0]
    mismatch = score_matrix[1]
    gap_opening = score_matrix[2]
    gap_extension = score_matrix[3]
    print("match = %d\nmismatch = %d\ngap_opening = %d\ngap_extension = %d" %(match, mismatch, gap_opening, gap_extension))

    n_row = len(str1) + 1
    n_col = len(str2) + 1

    dp_table_middle = [
        [
            0 for _ in range(n_col)
        ] for _ in range(n_row)
    ]
    dp_table_lower = [
        [
            0 for _ in range(n_col)
        ] for _ in range(n_row)
    ]
    dp_table_upper = [
        [
            0 for _ in range(n_col)
        ] for _ in range(n_row)
    ]

    backtrack_matrix_lower = [
        [
            "↓" for _ in range(n_col)
        ] for _ in range(n_row)
    ]
    backtrack_matrix_middle = [
        [
            "↘︎" for _ in range(n_col)
        ] for _ in range(n_row)
    ]
    backtrack_matrix_upper = [
        [
            "→" for _ in range(n_col)
        ] for _ in range(n_row)
    ]

    # a) For the middle matrix, the left and top edges are [0,-sigma,-sigma-epsilon,-sigma-2*epsilon,-sigma-3*epsilon,...]
    if n_row > 1:
        dp_table_middle[1][0] = -gap_opening
    if n_col > 1:
        dp_table_middle[0][1] = -gap_opening
    for i in range(2, n_row):
        dp_table_middle[i][0] = dp_table_middle[i - 1][0] - gap_extension
    for j in range(2, n_col):
        dp_table_middle[0][j] = dp_table_middle[0][j - 1] - gap_extension
    # b) For the lower matrix, the top edge is negative infinity, the left edge is the same as the corresponding Middle value
    # c) For the upper matrix, the left edge is negative infinity, the top edge is the same as the corresponding Middle value
    for i in range(1, n_row):
        dp_table_lower[i][0] = dp_table_middle[i][0]
        dp_table_upper[i][0] = -999
    for j in range(1, n_col):
        dp_table_lower[0][j] = -999
        dp_table_upper[0][j] = dp_table_middle[0][j]
    dp_table_lower[0][0] = -999
    dp_table_upper[0][0] = -999

    max_score = 0
    for i in range(1, n_row):
        for j in range(1, n_col):
            # is_matched = match
            # if str1[i - 1] != str2[j - 1]:
            #     is_matched = -mismatch
            
            # dp_table_lower[i][j] = max(
            #     dp_table_lower[i - 1][j] - gap_extension, 
            #     dp_table_middle[i - 1][j] - gap_opening
            # )
            # dp_table_upper[i][j] = max(
            #     dp_table_upper[i][j - 1] - gap_extension, 
            #     dp_table_middle[i][j - 1] - gap_opening
            # )
            # dp_table_middle[i][j] = max(
            #     dp_table_lower[i][j], 
            #     dp_table_upper[i][j], 
            #     dp_table_middle[i - 1][j - 1] + is_matched
            # )

            # max_score = max(
            #     dp_table_lower[i][j], 
            #     dp_table_upper[i][j], 
            #     dp_table_middle[i][j], 
            # )

            # # gaps & extensions are more tolerated
            # if dp_table_lower[i][j] == max_score:
            #     backtrack_matrix[i][j] = "↓"
            # elif dp_table_upper[i][j] == max_score:
            #     backtrack_matrix[i][j] = "→"
            # elif dp_table_middle[i][j] == max_score:
            #     backtrack_matrix[i][j] = "↘︎"

            # source code refered to https://github.com/xuwd11/Coursera-Bioinformatics/blob/master/33_01_GlobalAlignmentGapPenalties.py
            score1 = dp_table_lower[i - 1][j] - gap_extension
            score2 = dp_table_middle[i - 1][j] - gap_opening
            lower_score = max(score1, score2)
            dp_table_lower[i][j] = lower_score
            if lower_score == score2:
                # from middle
                backtrack_matrix_lower[i][j] = "M"
            
            score1 = dp_table_upper[i][j - 1] - gap_extension
            score2 = dp_table_middle[i][j - 1] - gap_opening
            upper_score = max(score1, score2)
            dp_table_upper[i][j] = upper_score
            if upper_score == score2:
                # from middle
                backtrack_matrix_upper[i][j] = "M"
            
            is_matched = match
            if str1[i - 1] != str2[j - 1]:
                is_matched = - mismatch
            score3 = dp_table_middle[i - 1][j - 1] + is_matched
            middle_score = max(
                lower_score, 
                score3, 
                upper_score
            )
            dp_table_middle[i][j] = middle_score
            if middle_score == lower_score:
                # from lower
                backtrack_matrix_middle[i][j] = "L"
            elif middle_score == upper_score:
                # from upper
                backtrack_matrix_middle[i][j] = "U"
    max_score = max(
        dp_table_lower[n_row - 1][n_col - 1], 
        dp_table_middle[n_row - 1][n_col - 1], 
        dp_table_upper[n_row - 1][n_col - 1]
    )
    
    # print_2D_arr(dp_table_lower)
    # print_2D_arr(dp_table_middle)
    # print_2D_arr(dp_table_upper)
    # print_2D_arr(backtrack_matrix_lower)
    # print_2D_arr(backtrack_matrix_middle)
    # print_2D_arr(backtrack_matrix_upper)

    backtrack_matrices = [
        backtrack_matrix_lower, 
        backtrack_matrix_middle, 
        backtrack_matrix_upper
    ]
    alignments = get_affine_gap_penalty_alignments(backtrack_matrices, str1, str2)
    # print(alignments[0])
    # print(alignments[1])

    return max_score, alignments[0], alignments[1]
    
def get_affine_gap_penalty_alignments(backtrack_matrices, str1, str2):
    i = len(str1)
    j = len(str2)
    alignment1 = []
    alignment2 = []

    backtrack_matrix_lower = backtrack_matrices[0]
    backtrack_matrix_middle = backtrack_matrices[1]
    backtrack_matrix_upper = backtrack_matrices[2]

    # start from middle level
    level = "M"
    while i > 0 and j > 0:
        # print(alignment1)
        # print(alignment2)
        # print()

        if level == "L":
            if backtrack_matrix_lower[i][j] == "M":
                level = "M"
            alignment1.append(str1[i - 1])
            alignment2.append("-")
            i -= 1
            continue
        elif level == "U":
            if backtrack_matrix_upper[i][j] == "M":
                level = "M"
            alignment1.append("-")
            alignment2.append(str2[j - 1])
            j -= 1
            continue
        elif level == "M":
            if backtrack_matrix_middle[i][j] == "↘︎":
                alignment1.append(str1[i - 1])
                alignment2.append(str2[j - 1])
                i -= 1
                j -= 1
                continue
            else:
                level = backtrack_matrix_middle[i][j]        
    
    # fill up the remaining chars in str1 or str2
    while i > 0:
        alignment1.append(str1[i - 1])
        alignment2.append("-")
        i -= 1

    while j > 0:
        alignment1.append("-")
        alignment2.append(str2[j - 1])
        j -= 1

    # the alignment is generated in the reversed order
    alignment1 = [item for item in reversed(alignment1)]
    alignment2 = [item for item in reversed(alignment2)]

    return "".join(alignment1), "".join(alignment2)

def get_max_and_index(arr):
    max_val = max(arr)
    return max_val, arr.index(max_val)

# Input: Two strings and a matrix score.
# Output: A middle edge in the alignment graph of 
#       these strings (where the edge lengths are defined by score).
def get_middle_edge_in_alignment_graph(score_matrix, str1, str2):
    match = score_matrix[0]
    mismatch = score_matrix[1]
    indel = score_matrix[2]

    n_row = len(str1) + 1
    n_col = len(str2) + 1
    middle_col = (n_col - 1) // 2
    ceil_col = middle_col + 1
    # print("middle_col: %d\nceil_col: %d" %(middle_col, ceil_col))

    dp_table = [0 for _ in range(n_row)]
    for i in range(1, n_row):
        dp_table[i] = dp_table[i - 1] - indel

    middle_nodes_scores = [0 for _ in range(n_row)]
    if n_col == 2:
        for _ in range(n_row):
            middle_nodes_scores[_] += dp_table[_]
    ceil_nodes_scores = [0 for _ in range(n_row)]

    # from source
    # print_arr(dp_table)
    for j in range(1, ceil_col + 1):
        prev_i_minus_1 = dp_table[0]
        dp_table[0] -= indel
        for i in range(1, n_row):
            temp = dp_table[i]

            is_match = match
            if str1[i - 1] != str2[j - 1]:
                is_match = -mismatch

            dp_table[i] = max(
                prev_i_minus_1 + is_match, 
                dp_table[i] - indel,        # insertion
                dp_table[i - 1] - indel     # deletion
            )

            prev_i_minus_1 = temp
        
        if j == middle_col:
            for i in range(n_row):
                middle_nodes_scores[i] += dp_table[i]
        if j == ceil_col:
            for i in range(n_row):
                ceil_nodes_scores[i] += dp_table[i]

    # to sink
    dp_table = [0 for _ in range(n_row)]
    for i in range(1, n_row):
        dp_table[i] = dp_table[i - 1] - indel
    if n_col == 2:
        for i in range(n_row):
            ceil_nodes_scores[i] += [_ for _ in reversed(dp_table)][i]

    reversed_str1 = "".join(reversed(str1))
    reversed_str2 = "".join(reversed(str2))
    for j in range(1, n_col - middle_col):
        prev_i_minus_1 = dp_table[0]
        dp_table[0] -= indel
        for i in range(1, n_row):
            temp = dp_table[i]

            is_match = match
            if reversed_str1[i - 1] != reversed_str2[j - 1]:
                is_match = -mismatch

            dp_table[i] = max(
                prev_i_minus_1 + is_match, 
                dp_table[i] - indel,        # insertion
                dp_table[i - 1] - indel     # deletion
            )

            prev_i_minus_1 = temp

        if j == n_col - middle_col - 1:
            for i in range(n_row):
                middle_nodes_scores[i] += [_ for _ in reversed(dp_table)][i]
        if j == n_col - ceil_col - 1:
            for i in range(n_row):
                ceil_nodes_scores[i] += [_ for _ in reversed(dp_table)][i]
    
    starting_score, middle_row = get_max_and_index(middle_nodes_scores)

    if n_row == middle_row + 1:
        return [[middle_row, middle_col], [middle_row, middle_col + 1]]

    deletion = middle_nodes_scores[middle_row + 1]
    insertion = ceil_nodes_scores[middle_row]
    match_or_mismatch = ceil_nodes_scores[middle_row + 1]
    max_score = max(deletion, insertion, match_or_mismatch)
    if insertion == max_score:
        ceil_row = middle_row
        ceil_col = middle_col + 1
        return [middle_row, middle_col], [ceil_row, ceil_col]
    if deletion == max_score:
        ceil_row = middle_row + 1
        ceil_col = middle_col
        return [middle_row, middle_col], [ceil_row, ceil_col]
    if match_or_mismatch == max_score:
        ceil_row = middle_row + 1
        ceil_col = middle_col + 1
        return [middle_row, middle_col], [ceil_row, ceil_col]


linear_space_alignment_path = []    
# a recursive function implementing divide-and-conquer algorithm
# only linear space will be used, while runtime is still O(m*n)
def do_linear_space_alignment(score_matrix, str1, str2, top, bottom, left, right):
    # end of recursive
    if left == right:
        for _ in range(bottom - top):
            linear_space_alignment_path.append("↓")
    elif bottom == top:
        for _ in range(right - left):
            linear_space_alignment_path.append("→")
    else:
        mid_edge = get_mid_edge(score_matrix, str1, str2, top, bottom, left, right)
        mid_node = mid_edge[0]
        mid_node_next = mid_edge[1]
        # recursive to left top rectangle
        do_linear_space_alignment(score_matrix, str1, str2, top, mid_node[0], left, mid_node[1])

        # print(mid_node)
        # print(interpret_edge(mid_edge))
        # print(mid_node_next)
        # print()

        linear_space_alignment_path.append(interpret_edge(mid_edge))

        # recursive to right bottom rectangle
        do_linear_space_alignment(score_matrix, str1, str2, mid_node_next[0], bottom, mid_node_next[1], right)

# wrapper function of get_middle_edge_in_alignment
def get_mid_edge(score_matrix, str1, str2, top, bottom, left, right):
    mid_edge = get_middle_edge_in_alignment_graph(score_matrix, str1[top : bottom], str2[left : right])
    mid_edge[0][0] += top
    mid_edge[0][1] += left
    mid_edge[1][0] += top
    mid_edge[1][1] += left
    return mid_edge

# return the direction of the edge according to the coordinates
# of the starting and ending nodes of this edge
def interpret_edge(edge):
    starting_node = edge[0]
    ending_node = edge[1]

    if (starting_node[0] + 1 == ending_node[0]) and \
        (starting_node[1] + 1 == ending_node[1]):
        return "↘︎"
    elif (starting_node[0] + 1 == ending_node[0]) and \
        (starting_node[1] == ending_node[1]):
        return "↓"
    elif (starting_node[0] == ending_node[0]) and \
        (starting_node[1] + 1 == ending_node[1]):
        return "→"

# wrapper of the recursive function do_linear_space_alignment
def linear_space_alignment(score_matrix, str1, str2):
    top = 0
    bottom = len(str1)
    left = 0
    right = len(str2)

    do_linear_space_alignment(score_matrix, str1, str2, top, bottom, left, right)

    # backtrack using the linear space alignment path
    score = 0
    i, j = len(str1), len(str2)
    align1 = []
    align2 = []
    while len(linear_space_alignment_path) != 0:
        current_direction = linear_space_alignment_path.pop()
        # print(current_direction)
        if current_direction == "↘︎":
            align1.append(str1[i - 1])
            align2.append(str2[j - 1])
            if str1[i - 1] != str2[j - 1]:
                score -= score_matrix[1] 
            else:
                score += score_matrix[0]
            i -= 1
            j -= 1
        elif current_direction == "→":
            align1.append("-")
            align2.append(str2[j - 1])
            j -= 1
            score -= score_matrix[2]
        elif current_direction == "↓":
            align1.append(str1[i - 1])
            align2.append("-")
            i -= 1
            score -= score_matrix[2]
    
    aligned_str1 = "".join(reversed(align1))
    aligned_str2 = "".join(reversed(align2))
    # print("%d\n%s\n%s" %(score, aligned_str1, aligned_str2))
    return score, aligned_str1, aligned_str2

# Input: Three DNA strings of length at most 10.
# Output: The length of a longest common subsequence of these three strings, 
#       followed by a multiple alignment of the three strings corresponding to such an alignment.
# the score of a column of the alignment matrix is equal to 1 
# if all of the column's symbols are identical, and 0 if even one symbol disagrees.
def multiple_LCS_alignment(str1, str2, str3):
    len1 = len(str1) + 1
    len2 = len(str2) + 1
    len3 = len(str3) + 1

    dp_table = [
        [
            [
                0 for _ in range(len3)
            ] for _ in range(len2)
        ] for _ in range(len1)
    ]

    backtrack_table = [
        [
            [
                0 for _ in range(len3)
            ] for _ in range(len2)
        ] for _ in range(len1)
    ]
    for i in range(1, len1):
        backtrack_table[i][0][0] = 1
    for j in range(1, len2):
        backtrack_table[0][j][0] = 2
    for k in range(1, len3):
        backtrack_table[0][0][k] = 3
    for i in range(1, len1):
        for j in range(1, len2):
            backtrack_table[i][j][0] = 4
    for i in range(1, len1):
        for k in range(1, len3):
            backtrack_table[i][0][k] = 5
    for j in range(1, len2):
        for k in range(1, len3):
            backtrack_table[0][j][k] = 6

    for i in range(1, len1):
        for j in range(1, len2):
            for k in range(1, len3):
                is_match = 0
                if str1[i - 1] == str2[j - 1] and \
                str2[j - 1] == str3[k - 1]:
                    is_match = 1
                dp_table[i][j][k] = max(
                    dp_table[i - 1][j][k], 
                    dp_table[i][j - 1][k], 
                    dp_table[i][j][k - 1], 
                    dp_table[i - 1][j - 1][k], 
                    dp_table[i - 1][j][k - 1], 
                    dp_table[i][j - 1][k - 1], 
                    dp_table[i - 1][j - 1][k - 1] + is_match 
                )

                if dp_table[i][j][k] == dp_table[i - 1][j][k]:
                    backtrack_table[i][j][k] = 1
                    continue
                if dp_table[i][j][k] == dp_table[i][j - 1][k]:
                    backtrack_table[i][j][k] = 2
                    continue
                if dp_table[i][j][k] == dp_table[i][j][k - 1]:
                    backtrack_table[i][j][k] = 3
                    continue
                if dp_table[i][j][k] == dp_table[i - 1][j - 1][k]:
                    backtrack_table[i][j][k] = 4
                    continue
                if dp_table[i][j][k] == dp_table[i - 1][j][k - 1]:
                    backtrack_table[i][j][k] = 5
                    continue
                if dp_table[i][j][k] == dp_table[i][j - 1][k - 1]:
                    backtrack_table[i][j][k] = 6
                    continue
                if dp_table[i][j][k] == dp_table[i - 1][j - 1][k - 1] + is_match:
                    backtrack_table[i][j][k] = 7
                    continue
    
    align1 = []
    align2 = []
    align3 = []
    i = len1 - 1
    j = len2 - 1
    k = len3 - 1
    while i > 0 and j > 0 and k > 0:
        direction = backtrack_table[i][j][k]

        if direction == 1:
            align1.append(str1[i - 1])
            align2.append("-")
            align3.append("-")
            i -= 1
            continue
        if direction == 2:
            align1.append("-")
            align2.append(str2[j - 1])
            align3.append("-")
            j -= 1
            continue
        if direction == 3:
            align1.append("-")
            align2.append("-")
            align3.append(str3[k - 1])
            k -= 1
            continue
        if direction == 4:
            align1.append(str1[i - 1])
            align2.append(str2[j - 1])
            align3.append("-")
            i -= 1
            j -= 1
            continue
        if direction == 5:
            align1.append(str1[i - 1])
            align2.append("-")
            align3.append(str3[k - 1])
            i -= 1
            k -= 1
            continue
        if direction == 6:
            align1.append("-")
            align2.append(str2[j - 1])
            align3.append(str3[k - 1])
            j -= 1
            k -= 1
            continue
        if direction == 7:
            align1.append(str1[i - 1])
            align2.append(str2[j - 1])
            align3.append(str3[k - 1])
            i -= 1
            j -= 1
            k -= 1
            continue
    while i > 0 and j > 0:
        align1.append(str1[i - 1])
        align2.append(str2[j - 1])
        align3.append("-")
        i -= 1
        j -= 1
    while i > 0 and k > 0:
        align1.append(str1[i - 1])
        align2.append("-")
        align3.append(str3[k - 1])
        i -= 1
        k -= 1
    while j > 0 and k > 0:
        align1.append("-")
        align2.append(str2[j - 1])
        align3.append(str3[k - 1])
        j -= 1
        k -= 1
    while i > 0:
        align1.append(str1[i - 1])
        align2.append("-")
        align3.append("-")
        i -= 1
    while j > 0:
        align1.append("-")
        align2.append(str2[j - 1])
        align3.append("-")
        j -= 1
    while k > 0:
        align1.append("-")
        align2.append("-")
        align3.append(str3[k - 1])
        k -= 1
    
    aligned_str1 = "".join(reversed(align1))
    aligned_str2 = "".join(reversed(align2))
    aligned_str3 = "".join(reversed(align3))

    # print(dp_table[-1][-1][-1])
    # print(aligned_str1)
    # print(aligned_str2)
    # print(aligned_str3)

    return dp_table[-1][-1][-1], aligned_str1, aligned_str2, aligned_str3

if __name__ == "__main__":
    # 1.1 step 8:
    # score_matrix = [1, 3, 2, 1]
    # str1 = "GA"
    # str2 = "GTTA"
    # results = affine_gap_penalty_alignment(score_matrix, str1, str2)
    # print(results[0])
    # print(results[1])
    # print(results[2])

    # debug dataset
    # os.chdir("./AffineAlignment")
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
    #     # if not "9" in input:
    #     #     continue

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         score_matrix = [int(score) for score in f_input.readline().strip().split(" ")]
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()
    #         results = affine_gap_penalty_alignment(score_matrix, str1, str2)
    #         result += str(results[0]) + "\n" + results[1] + "\n" + results[2]
    #         print("str1: %s\nstr2: %s" %(str1, str2))

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter), "correct!")
    #         else:
    #             print("debug_%d" %(counter), "wrong!")

    #     counter += 1
        
    #     print()

    # with open("./dataset_249_8.txt", "r") as f:
    #     result = ""
    #     score_matrix = [int(score) for score in f.readline().strip().split(" ")]
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     results = affine_gap_penalty_alignment(score_matrix, str1, str2)
    #     result += str(results[0]) + "\n" + results[1] + "\n" + results[2]
    #     print("str1: %s\nstr2: %s" %(str1, str2))
    #     print(result)

    # 1.2 step 12:
    # score_matrix = [1, 1, 2]
    # str1 = "GAGA"
    # str2 = "GAT"
    # get_middle_edge_in_alignment_graph(score_matrix, str1, str2)

    # debug dataset
    # os.chdir("./MiddleEdge")
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
    #     # if not("7" in input):
    #     #     continue

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         score_matrix = [int(score) for score in f_input.readline().strip().split(" ")]
    #         str2 = f_input.readline().strip()
    #         str1 = f_input.readline().strip()
    #         print("str1: %s\nstr2: %s" %(str1, str2))
    #         results = get_middle_edge_in_alignment_graph(score_matrix, str1, str2)
    #         result += "%d %d\n%d %d" %(
    #             results[0][0], 
    #             results[0][1], 
    #             results[1][0], 
    #             results[1][1], 
    #         )

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result: \n%s" %result)
    #         print("Answer: \n%s" %output)
    #         if result == output:
    #             print("debug_%d" %(counter), "correct!")
    #         else:
    #             print("debug_%d" %(counter), "wrong!")

    #     counter += 1
        
    #     print()

    # with open("./dataset_250_12.txt", "r") as f:
    #     result = ""
    #     score_matrix = [int(score) for score in f.readline().strip().split(" ")]
    #     str2 = f.readline().strip()
    #     str1 = f.readline().strip()
    #     print("str1: %s\nstr2: %s" %(str1, str2))
    #     results = get_middle_edge_in_alignment_graph(score_matrix, str1, str2)
    #     result += "%d %d\n%d %d" %(
    #         results[0][0], 
    #         results[0][1], 
    #         results[1][0], 
    #         results[1][1], 
    #     )
    #     print(result)

    # 1.2 step 14, not completely solved :(
    # score_matrix = [1, 1, 2]
    # str1 = "GAGA"
    # str2 = "GAT"
    # linear_space_alignment(score_matrix, str1, str2)

    # debug dataset
    # os.chdir("./LinearSpaceAlignment")
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
    #     # counter = 4 - 1
    #     # if not(str(counter + 1) in input):
    #     #     continue

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         score_matrix = [int(score) for score in f_input.readline().strip().split(" ")]
    #         str2 = f_input.readline().strip()
    #         str1 = f_input.readline().strip()
    #         print("str1: %s\nstr2: %s" %(str1, str2))
    #         results = linear_space_alignment(score_matrix, str1, str2)
    #         result += str(results[0]) + "\n" + results[2] + "\n" + results[1]

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result: \n%s" %result)
    #         print("Answer: \n%s" %output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")

    #     counter += 1
        
    #     print()

    with open("./dataset_250_14.txt", "r") as f:
        result = ""
        score_matrix = [int(score) for score in f.readline().strip().split(" ")]
        str2 = f.readline().strip()
        str1 = f.readline().strip()
        results = linear_space_alignment(score_matrix, str1, str2)
        result += str(results[0]) + "\n" + results[1] + "\n" + results[2]
        print("str1: %s\nstr2: %s" %(str1, str2))
        print(result)

    # 1.3 step 5:
    # str1 = "ATATCCG"
    # str2 = "TCCGA"
    # str3 = "ATGTACTG"
    # multiple_LCS_alignment(str1, str2, str3)

    # debug dataset
    # os.chdir("./MultipleAlignment")
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
    #     # if not("7" in input):
    #     #     continue

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()
    #         str3 = f_input.readline().strip()
    #         print("str1: %s\nstr2: %s\nstr3: %s" %(str1, str2, str3))
    #         results = multiple_LCS_alignment(str1, str2, str3)
    #         result += "%d\n%s\n%s\n%s" %(
    #             results[0], 
    #             results[1], 
    #             results[2], 
    #             results[3]
    #         )

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result: \n%s" %result)
    #         print("Answer: \n%s" %output)
    #         if result == output:
    #             print("debug_%d" %(counter), "correct!")
    #         else:
    #             print("debug_%d" %(counter), "wrong!")

    #     counter += 1
        
    #     print()

    # with open("./dataset_251_5.txt", "r") as f:
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     str3 = f.readline().strip()
    #     print("str1: %s\nstr2: %s\nstr3: %s" %(str1, str2, str3))
    #     results = multiple_LCS_alignment(str1, str2, str3)
    #     result = "%d\n%s\n%s\n%s" %(
    #         results[0], 
    #         results[1], 
    #         results[2], 
    #         results[3]
    #     )
    #     print(result)

    # week 3 Quiz
    # problem 1
    # print("\n### problem 1 ###")
    # str1 = "A-C--GTTAC"
    # str2 = "ATGCAG---T"
    # score_matrix = [1, 1, 4, 1]
    # score = 0
    # for i in range(len(str1)):
    #     # match
    #     if str1[i] == str2[i]:
    #         score += score_matrix[0]
    #     # mismatch
    #     elif str1[i] != str2[i] and \
    #         str1[i] != "-" and str2[i] != "-":
    #         score -= score_matrix[1]
    #     # gap
    #     elif str1[i] == "-":
    #         if i > 0:
    #             if str1[i - 1] == "-":
    #                 score -= score_matrix[3] 
    #             else:
    #                 score -= score_matrix[2]
    #         else:
    #             score -= score_matrix[2]
    #     elif str2[i] == "-":
    #         if i > 0:
    #             if str2[i - 1] == "-":
    #                 score -= score_matrix[3] 
    #             else:
    #                 score -= score_matrix[2]
    #         else:
    #             score -= score_matrix[2]
    
    # print("score: %d" %score)

    # print("\n### problem 4 ###")
    # str1 = "CCAATACGAC"
    # str2 = "GCCTTACGCT"
    # str3 = "CCCTAGCGGC"
    # results = multiple_LCS_alignment(str1, str2, str3)
    # result = "%d\n%s\n%s\n%s" %(
    #     results[0], 
    #     results[1], 
    #     results[2], 
    #     results[3]
    # )
    # print(result)
