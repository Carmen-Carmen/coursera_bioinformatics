import math
import pyperclip
import os
import sys

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    result = result.strip()
    pyperclip.copy(result)
    
    print(result)

# the global alignment problem
# Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
# Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
def global_alignment(score_matrix, str1, str2):
    # score_matrix = [match_reward, mismatch_penalty, indel_penalty]

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
        dp_table[i][0] = dp_table[i - 1][0] - score_matrix[2]
    for j in range(1, len(str2) + 1):
        dp_table[0][j] = dp_table[0][j - 1] - score_matrix[2]

    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            is_matched = score_matrix[0]
            if str1[i - 1] != str2[j - 1]:
                is_matched = -score_matrix[1]

            score_match_or_mismatch = dp_table[i - 1][j - 1] + is_matched
            score_deletion = dp_table[i - 1][j] - score_matrix[2]
            score_insertion = dp_table[i][j - 1] - score_matrix[2]
            dp_table[i][j] = max(score_match_or_mismatch, score_deletion, score_insertion)

            if dp_table[i][j] == dp_table[i - 1][j] - score_matrix[2]:
                # deletion
                backtrack_matrix[i][j] = "↓"
            elif dp_table[i][j] == dp_table[i][j - 1] - score_matrix[2]:
                # insertion
                backtrack_matrix[i][j] = "→"
            elif dp_table[i][j] == dp_table[i - 1][j - 1] + is_matched:
                # match/mismatch
                backtrack_matrix[i][j] = "↘︎"

    alignments = get_global_alignments(backtrack_matrix, str1, str2)
    return dp_table[i][j], alignments[0], alignments[1]

# align str1 and str2 with the backtrack matrix
# output: aligned str1 and str2
def get_global_alignments(backtrack_matrix, str1, str2):
    i = len(str1)
    j = len(str2)
    alignment1 = []
    alignment2 = []

    while i > 0 and j > 0:
        # print(alignment1)
        # print(alignment2)
        # print()
        
        if backtrack_matrix[i][j] == "↓":
            alignment1.append(str1[i - 1])
            alignment2.append("-")
            i -= 1
        elif backtrack_matrix[i][j] == "→":
            alignment1.append("-")
            alignment2.append(str2[j - 1])
            j -= 1
        else:
            alignment1.append(str1[i - 1])
            alignment2.append(str2[j - 1])
            i -= 1
            j -= 1
    
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

# the global alignment problem
# Input: A match reward, a mismatch penalty, an indel penalty, and two nucleotide strings.
# Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score.
def local_alignment(score_matrix, str1, str2):
    indel_penalty = 5

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

    largest_score = -sys.maxsize
    largest_index = [0, 0]
    # init of 1st row and col
    dp_table[0][0] = score_matrix[str1[0]][str2[0]]
    for i in range(1, len(str1) + 1):
        AA1 = str1[i - 1]
        AA2 = str2[0]
        dp_table[i][0] = score_matrix[AA1][AA2]
        current_score = score_matrix[AA1][AA2]
        if current_score > largest_score:
                largest_score = current_score
                largest_index = [i, 0]
    for j in range(1, len(str2) + 1):
        AA1 = str1[0]
        AA2 = str2[j - 1]
        dp_table[0][j] = score_matrix[AA1][AA2]
        current_score = score_matrix[AA1][AA2]
        if current_score > largest_score:
                largest_score = current_score
                largest_index = [0, j]

    for i in range(1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            AA1 = str1[i - 1]
            AA2 = str2[j - 1]

            score_match_or_mismatch = dp_table[i - 1][j - 1] + score_matrix[AA1][AA2]
            score_deletion = dp_table[i - 1][j] - indel_penalty
            score_insertion = dp_table[i][j - 1] - indel_penalty
            current_score = max(0, 
                                 score_match_or_mismatch, 
                                 score_deletion, 
                                 score_insertion)
            dp_table[i][j] = current_score

            if current_score > largest_score:
                largest_score = current_score
                largest_index = [i, j]

            if dp_table[i][j] == dp_table[i - 1][j] - indel_penalty:
                # deletion
                backtrack_matrix[i][j] = "↓"
            elif dp_table[i][j] == dp_table[i][j - 1] - indel_penalty:
                # insertion
                backtrack_matrix[i][j] = "→"
            elif dp_table[i][j] == dp_table[i - 1][j - 1] + score_matrix[AA1][AA2]:
                # match/mismatch
                backtrack_matrix[i][j] = "↘︎"
            elif dp_table[i][j] == 0:
                # directly from source
                backtrack_matrix[i][j] = "S"

    # print(largest_score)
    # print_arr(largest_index)

    # position of sink is the node with largest score
    alignments = get_local_alignments(backtrack_matrix, str1, str2, 
                                      largest_index[0], largest_index[1])
    return largest_score, alignments[0], alignments[1]

# align str1 and str2 with the backtrack matrix, starting from sink(i, j)
# output: aligned str1 and str2
def get_local_alignments(backtrack_matrix, str1, str2, sink_i, sink_j):
    i = sink_i
    j = sink_j
    alignment1 = []
    alignment2 = []

    while i > 0 and j > 0:
        # print(alignment1)
        # print(alignment2)
        # print()
        
        if backtrack_matrix[i][j] == "↓":
            alignment1.append(str1[i - 1])
            alignment2.append("-")
            i -= 1
        elif backtrack_matrix[i][j] == "→":
            alignment1.append("-")
            alignment2.append(str2[j - 1])
            j -= 1
        elif backtrack_matrix[i][j] == "↘︎":
            alignment1.append(str1[i - 1])
            alignment2.append(str2[j - 1])
            i -= 1
            j -= 1
        else:
            # end when score = 0, i.e. linked directly to the source node
            # backtrack_matrix[i][j] == "S"
            break
    
    # the alignment is generated in the reversed order
    alignment1 = [item for item in reversed(alignment1)]
    alignment2 = [item for item in reversed(alignment2)]

    return "".join(alignment1), "".join(alignment2)

# return the matrix by a 2-dimensional dict, 
# so that the matching score of AA1 and AA2 can be get by matrix[AA1][AA2]
def get_PAM_250_matrix():
    with open("./PAM250.txt", "r") as f:
        AAs = f.readline().strip().split("  ")
        PAM_250_matrix = {}
        # print_arr(AAs)
        # print(len(AAs))
        while True:
            line = f.readline().strip()
            if line == "":
                break
            current_AA = line.split(" ")[0]
            PAM_250_matrix[current_AA] = {}
            # print(current_AA)
            scores = line.split(" ")[1:]
            scores = [int(score) for score in scores if score != " " and score != ""]
            # print_arr(scores)
            # print(len(scores))
            for i in range(len(scores)):
                pairing_AA = AAs[i]
                PAM_250_matrix[current_AA][pairing_AA] = scores[i]
        
        return PAM_250_matrix

# find the edit distance between 2 strings
# input: 2 strings
# output: the edit distance
def get_edit_distance(str1, str2):
    # 1. perform the global alignment
    # 0 for matches, because the 2 strings should be aligned in the most efficient way, i.e. less steps
    # -1 for mismatches or indels
    score_matrix = [0, 1, 1]
    alignment_results = global_alignment(score_matrix, str1, str2)
    aligned_str1 = alignment_results[1]
    aligned_str2 = alignment_results[2]
    # print(aligned_str1)
    # print(aligned_str2)

    # 2. calculate the hamming distance of the two aligned strs
    return get_hamming_distance(aligned_str1, aligned_str2)

# input: two strs with the same length
def get_hamming_distance(str1, str2):
    l = len(str1)
    hamming_distance = 0
    for i in range(l):
        if str1[i] != str2[i]:
            hamming_distance += 1
    
    return hamming_distance

def get_BLOSUM_62_matrix():
    with open("./BLOSUM62.txt", "r") as f:
        AAs = f.readline().strip().split("  ")
        BLOSUM_62_matrix = {}
        # print_arr(AAs)
        # print(len(AAs))
        while True:
            line = f.readline().strip()
            if line == "":
                break
            current_AA = line.split(" ")[0]
            BLOSUM_62_matrix[current_AA] = {}
            # print(current_AA)
            scores = line.split(" ")[1:]
            scores = [int(score) for score in scores if score != " " and score != ""]
            # print_arr(scores)
            # print(len(scores))
            for i in range(len(scores)):
                pairing_AA = AAs[i]
                BLOSUM_62_matrix[current_AA][pairing_AA] = scores[i]
        
        return BLOSUM_62_matrix
    
# Input: Two amino acid strings. str1 is longer than str2
# Output: A highest-scoring fitting alignment between v and w. Use the BLOSUM62 scoring table.
def fitting_alignment(score_matrix, str1, str2):
    indel_penalty = 1

    is_gene_strand = False
    if isinstance(score_matrix, list):
        indel_penalty = score_matrix[2]
        is_gene_strand = True

    # row_num = len(str1), col_num = len(str2)
    backtrack_matrix = [
        [
            "O" for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    # dp_table stores the length of an LCS to each node
    dp_table = [
        [
            0 for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    # init the first col, set all to 0
    for i in range(1, len(str1) + 1):
        dp_table[i][0] = 0

    # find the first occurence of str2[0] in str1
    starting_row = 0
    for i in range(len(str1) + 1):
        if str2[0] == str1[i]:
            starting_row = i
            break
    
    # init the starting row (first occurance of str2[0] in str2)
    # set dp_table[starting_row][j] to (-j * indel_penalty)
    for j in range(1, len(str2) + 1):
        dp_table[starting_row][j] = dp_table[starting_row][j - 1] - indel_penalty

    # starting from the first occurence of str2[0] in str1
    # or 1st row if str[0] not in str1
    for i in range(starting_row + 1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            AA1 = str1[i - 1]
            AA2 = str2[j - 1]
            score_match_or_mismatch = 0
            if is_gene_strand:
                if AA1 == AA2:
                    score_match_or_mismatch = dp_table[i - 1][j - 1] + score_matrix[0]
                else:
                    score_match_or_mismatch = dp_table[i - 1][j - 1] - score_matrix[1]
            else:
                score_match_or_mismatch = dp_table[i - 1][j - 1] + score_matrix[AA1][AA2]
            score_deletion = dp_table[i - 1][j] - indel_penalty
            score_insertion = dp_table[i][j - 1] - indel_penalty
            current_score = max( score_match_or_mismatch, 
                                 score_deletion, 
                                 score_insertion)
            dp_table[i][j] = current_score

            if dp_table[i][j] == score_deletion:
                # deletion
                backtrack_matrix[i][j] = "↓"
            elif dp_table[i][j] == score_insertion:
                # insertion
                backtrack_matrix[i][j] = "→"
            elif dp_table[i][j] == score_match_or_mismatch:
                # match/mismatch
                backtrack_matrix[i][j] = "↘︎"

    # for line in dp_table:
    #     print_arr(line)
    # for line in backtrack_matrix:
    #     print_arr(line)

    # find the node with the highest score in the last col
    # the i position of sink node
    sink_i = len(str1)
    highest_score_in_last_col = dp_table[sink_i][len(str2)]
    for i in range(len(str1) + 1):
        if dp_table[i][len(str2)] >= highest_score_in_last_col:
            sink_i = i
            highest_score_in_last_col = dp_table[i][len(str2)]
    
    # perform global_alignment between str1[source_i : sink_i] and str2
    alignments = get_fitting_alignments(backtrack_matrix, str1, str2, starting_row, sink_i)
    return dp_table[sink_i][len(str2)], alignments[0], alignments[1]

def get_fitting_alignments(backtrack_matrix, str1, str2, source_i, sink_i):
    i = sink_i
    j = len(str2)
    alignment1 = []
    alignment2 = []

    while i > source_i and j > 0:
        # print(alignment1)
        # print(alignment2)
        # print()
        
        if backtrack_matrix[i][j] == "↓":
            alignment1.append(str1[i - 1])
            alignment2.append("-")
            i -= 1
        elif backtrack_matrix[i][j] == "→":
            alignment1.append("-")
            alignment2.append(str2[j - 1])
            j -= 1
        else:
            alignment1.append(str1[i - 1])
            alignment2.append(str2[j - 1])
            i -= 1
            j -= 1
    
    # fill up the remaining chars in str2
    while j > 0:
        alignment1.append("-")
        alignment2.append(str2[j - 1])
        j -= 1

    # the alignment is generated in the reversed order
    alignment1 = [item for item in reversed(alignment1)]
    alignment2 = [item for item in reversed(alignment2)]

    return "".join(alignment1), "".join(alignment2)

# Input: 
#   - A match reward, a mismatch penalty, an indel penalty, 
#   - and two nucleotide strings v and w.
# Output: 
#   - The score of an optimal overlap alignment of v and w, 
#   - followed by an alignment of a suffix v' of v and a prefix w' of w achieving this maximum score.
def overlap_alignment(score_matrix, str1, str2):
    # score_matrix = [match_reward, mismatch_penalty, indel_penalty]

    # row_num = len(str1), col_num = len(str2)
    backtrack_matrix = [
        [
            "O" for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    # dp_table stores the length of an LCS to each node
    dp_table = [
        [
            0 for _ in range(len(str2) + 1)
        ] for _ in range(len(str1) + 1)
    ]

    for i in range(1, len(str1) + 1):
        dp_table[i][0] = 0

    # find the first occurence of str2[0] in str1
    starting_row = 0
    for i in range(len(str1)):
        if str2[0] == str1[i]:
            starting_row = i
            break

    for j in range(1, len(str2) + 1):
        dp_table[starting_row][j] = dp_table[starting_row][j - 1] - score_matrix[2]

    for i in range(starting_row + 1, len(str1) + 1):
        for j in range(1, len(str2) + 1):
            is_matched = score_matrix[0]
            if str1[i - 1] != str2[j - 1]:
                is_matched = -score_matrix[1]

            score_match_or_mismatch = dp_table[i - 1][j - 1] + is_matched
            score_deletion = dp_table[i - 1][j] - score_matrix[2]
            score_insertion = dp_table[i][j - 1] - score_matrix[2]
            dp_table[i][j] = max(score_match_or_mismatch, score_deletion, score_insertion)

            # for overlap alignment, matches and mismatches should be considered first
            # then, insertion, finally deletion
            # because the shortest LCS is expected
            if dp_table[i][j] == score_match_or_mismatch:
                # match/mismatch
                backtrack_matrix[i][j] = "↘︎"
            elif dp_table[i][j] == score_insertion:
                # insertion
                backtrack_matrix[i][j] = "→"
            elif dp_table[i][j] == score_deletion:
                # deletion
                backtrack_matrix[i][j] = "↓"
    
    # for line in dp_table:
    #     print_arr(line)
    # print()
    # for line in backtrack_matrix:
    #     print_arr(line)
    
    highest_score_in_last_row = dp_table[-1][0]
    for j in range(len(str2) + 1):
        if dp_table[-1][j] >= highest_score_in_last_row:
            highest_score_in_last_row = dp_table[-1][j]
    # find the first index of highest score in last row
    sink_j = dp_table[-1][1 :].index(highest_score_in_last_row) + 1

    alignments = get_overlap_alignments(backtrack_matrix, str1, str2, starting_row, sink_j)
    return highest_score_in_last_row, alignments[0], alignments[1]

# alignment between str1[source_i :] and str2[: sink_j]
def get_overlap_alignments(backtrack_matrix, str1, str2, source_i, sink_j):
    i = len(str1)
    j = sink_j
    alignment1 = []
    alignment2 = []

    while i > source_i and j > 0:
        if backtrack_matrix[i][j] == "↓":
            alignment1.append(str1[i - 1])
            alignment2.append("-")
            i -= 1
        elif backtrack_matrix[i][j] == "→":
            alignment1.append("-")
            alignment2.append(str2[j - 1])
            j -= 1
        else:
            alignment1.append(str1[i - 1])
            alignment2.append(str2[j - 1])
            i -= 1
            j -= 1

    # Start backtracking from this node and stop when you reach Column_0
    # fill up the remaining chars in str2
    while j > 0:
        alignment1.append("-")
        alignment2.append(str2[j - 1])
        j -= 1

    # the alignment is generated in the reversed order
    alignment1 = [item for item in reversed(alignment1)]
    alignment2 = [item for item in reversed(alignment2)]

    return "".join(alignment1), "".join(alignment2)

if __name__ == "__main__":
    # 1.2 step 3
    # str1 = "G"
    # str2 = "ACATACGATAG"
    # score_matrix = [int(score) for score in "3 1 2".split(" ")]

    # print(global_alignment(score_matrix, str1, str2)[0])
    # print(global_alignment(score_matrix, str1, str2)[1])
    # print(global_alignment(score_matrix, str1, str2)[2])

    # debug dataset
    # os.chdir("./GlobalAlignment")
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
    #         score_matrix = [int(score) for score in f_input.readline().strip().split(" ")]
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()
    #         results = global_alignment(score_matrix, str1, str2)
    #         result += str(results[0]) + "\n" + results[1] + "\n" + results[2]

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")
        
    #     counter += 1
    #     print()

    # with open("./dataset_247_3.txt", "r") as f:
    #     score_matrix = [int(score) for score in f.readline().strip().split(" ")]
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     results = global_alignment(score_matrix, str1, str2)
    #     result = str(results[0]) + "\n" + results[1] + "\n" + results[2]
    #     print(result)

    # 1.2 step 10:
    # PAM_250_matrix = get_PAM_250_matrix()
    # str1 = "MEANLY"
    # str2 = "PENALTY"
    # local_alignment(PAM_250_matrix, str1, str2)

    # with open("./dataset_247_10.txt", "r") as f:
    #     PAM_250_matrix = get_PAM_250_matrix()
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     results = local_alignment(PAM_250_matrix, str1, str2)
    #     result = str(results[0]) + "\n" + results[1] + "\n" + results[2]
    #     print(result)

    # 1.3 step 3:
    # get_edit_distance("GAGA", "GAT")
    # debug dataset
    # os.chdir("./EditDistance")
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
    #         result += str(get_edit_distance(str1, str2))

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")
        
    #     counter += 1
    #     print()

    # with open("./dataset_248_3.txt", "r") as f:
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     print(get_edit_distance(str1, str2))

    # 1.3 step 5:
    # BLOSUM_62_matrix = get_BLOSUM_62_matrix()
    # str1 = "DISCREPANTLY"
    # str2 = "PATENT"
    # results = fitting_alignment(BLOSUM_62_matrix, str1, str2)
    # print(results[0])
    # print(results[1])
    # print(results[2])

    # debug dataset
    # os.chdir("./FittingAlignment")
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
    #         score_matrix = get_BLOSUM_62_matrix()
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()
    #         results = fitting_alignment(score_matrix, str1, str2)
    #         result += str(results[0]) + "\n" + results[1] + "\n" + results[2]

    #     with open("./outputs/%s" %outputs[counter], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter + 1), "correct!")
    #         else:
    #             print("debug_%d" %(counter + 1), "wrong!")
        
    #     counter += 1
    #     print()

    # with open("./dataset_248_5.txt", "r") as f:
    #     BLOSUM_62_matrix = get_BLOSUM_62_matrix()
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()
    #     results = fitting_alignment(BLOSUM_62_matrix, str1, str2)
    #     print(results[0])
    #     print(results[1])
    #     print(results[2])

    # 1.3 step 7:
    # score_matrix = [int(score) for score in "1 1 2".split(" ")]
    # str1 = "GAGA"
    # str2 = "GAT"
    # overlap_alignment(score_matrix, str1, str2)

    # debug dataset
    # os.chdir("./OverlapAlignment")
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
    #     # if not "6" in input:
    #     #     continue
    #     counter = int(input[-5])

    #     print("###### %s ######" %input)
    #     result = ""
    #     with open("./inputs/%s" %input, "r") as f_input:
    #         score_matrix = [int(score) for score in f_input.readline().strip().split(" ")]
    #         str1 = f_input.readline().strip()
    #         str2 = f_input.readline().strip()

    #         print("str1: %s\nstr2: %s" %(str1, str2))

    #         results = overlap_alignment(score_matrix, str1, str2)
    #         result += str(results[0]) + "\n" + results[1] + "\n" + results[2]

    #     with open("./outputs/%s" %outputs[counter - 1], "r") as f_output:
    #         output = f_output.read().strip()
    #         print("Result:", result)
    #         print("Correct:", output)
    #         if result == output:
    #             print("debug_%d" %(counter), "correct!")
    #         else:
    #             print("debug_%d" %(counter), "wrong!")
        
    #     print()

    # with open("./dataset_248_7.txt", "r") as f:
    #     score_matrix = [int(score) for score in f.readline().strip().split(" ")]
    #     str1 = f.readline().strip()
    #     str2 = f.readline().strip()

    #     print("str1: %s\nstr2: %s" %(str1, str2))

    #     results = overlap_alignment(score_matrix, str1, str2)
    #     print(results[0])
    #     print(results[1])
    #     print(results[2])

    # week2 Quiz
    # # problem 1
    # print("###### problem 1 ######")
    # score_matrix = [1, 1, 2]
    # str1 = "A-C--GTTAC"
    # str2 = "ATGCAG---T"
    # score = 0
    # for i in range(len(str1)):
    #     if str1[i] == str2[i]:
    #         score += score_matrix[0]
    #     elif str1[i] == "-" or str2[i] == "-":
    #         score -= score_matrix[2]
    #     elif str1[i] != str2[i]:
    #         score -= score_matrix[1]
    # print("score: %d" %score)

    # # problem 2
    # print("\n###### problem 2 ######")
    # score_matrix = [1, 3, 1]
    # str1 = "ACAGTAGACAC"
    # str2 = "ATAC-AGATAC"
    # score = 0
    # for i in range(len(str1)):
    #     if str1[i] == str2[i]:
    #         score += score_matrix[0]
    #     elif str1[i] == "-" or str2[i] == "-":
    #         score -= score_matrix[2]
    #     elif str1[i] != str2[i]:
    #         score -= score_matrix[1]
    # print("score: %d" %score)

    # # problem 3
    # print("\n###### problem 3 ######")
    # score_matrix = [1, 1, 1]
    # str1 = "GATACACT"
    # str2 = "ACGACCACAGATACCGCTATTCACTATATCGTT"
    # results = fitting_alignment(score_matrix, str2, str1)
    # print(results[0])
    # print(results[1])
    # print(results[2])

    # # problem 4
    # print("\n###### problem 4 ######")
    # score_matrix = [1, 0, 2]
    # str1 = "GATACAGCACACTAGTACTACTTGAC"
    # str2 = "CTATAGTCTTAACGATATACGAC"
    # results = overlap_alignment(score_matrix, str1, str2)
    # print("score by my algorithm: %d" %results[0])
    # print(results[1])
    # print(results[2])
    
    # str1 = "CTAGTACTACTTGAC"
    # str2 = "CTA-TAGT-CTTAAC"
    # score = 0
    # for i in range(len(str1)):
    #     if str1[i] == str2[i]:
    #         score += score_matrix[0]
    #     elif str1[i] == "-" or str2[i] == "-":
    #         score -= score_matrix[2]
    #     elif str1[i] != str2[i]:
    #         score -= score_matrix[1]            
    # print("score of alignment provided: %d" %score)
    # print(str1)
    # print(str2)

    # 1.2 step 14 
    # only using global alignment
    with open("./dataset_250_14.txt", "r") as f:
        score_matrix = [int(score) for score in f.readline().strip().split(" ")]
        str1 = f.readline().strip()
        str2 = f.readline().strip()
        results = global_alignment(score_matrix, str1, str2)
        result = str(results[0]) + "\n" + results[1] + "\n" + results[2]
        print(result)
