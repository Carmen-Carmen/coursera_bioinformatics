from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
from coursera_bioinformatics.utils import protein_utils
import pyperclip
import re
import copy
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_4.week4 import *

# Input: A space-delimited spectral vector Spectrum' and an amino acid string Proteome.
# Output: A substring of Proteome with maximum score against Spectrum'.
def identify_peptide_from_proteome(spectral_vector: list[int], proteome: str): 
    max_score = - float("inf")
    max_score_peptides = []

    peptide_weight = len(spectral_vector)
    min_peptide_len = peptide_weight // max(protein_utils.ALL_AA_WEIGHTS)
    max_peptide_len = peptide_weight // min(protein_utils.ALL_AA_WEIGHTS)

    # print("max peptide length: %d" %max_peptide_len)
    # print("min peptide length: %d" %min_peptide_len)
    # traverse as few substrings as possible!
    for i in range(len(proteome) - max_peptide_len + 1): 
        for j in range(min_peptide_len, max_peptide_len + 1): 
            peptide = proteome[i : i + j]
            # print(peptide)
            score = score_peptide_to_amplitude_spectral_vector(peptide, spectral_vector)
            if score > max_score: 
                max_score_peptides = []
                max_score = score
            if score == max_score: 
                max_score_peptides.append(peptide)
    
    return max_score_peptides[0], max_score

# Input: A set of space-delimited spectral vectors SpectralVectors, 
#       an amino acid string Proteome, 
#       and an integer threshold.
# Output: The set PSMthreshold(Proteome, SpectralVectors).
def search_peptide_spectrum_matches(spectral_vector_list: list[list[int]], proteome: str, score_threshold: int): 
    peptide_spectrum_matches = []
    for vector in spectral_vector_list: 
        results = identify_peptide_from_proteome(vector, proteome)
        peptide = results[0]
        score = results[1]
        if score >= score_threshold: 
            peptide_spectrum_matches.append((peptide, vector))
    
    return peptide_spectrum_matches

# use dynamic programming
def calculate_size_of_spectrum_dictionary(spectral_vector: list[int], threshold: int, max_score: int): 
    vector = [0] + spectral_vector

    # method 1: use 2-d array as dp_table
    dp_table = [[0 for j in range(len(vector))] for i in range(max_score + 1)]

    # because a blank string scores 0
    dp_table[0][0] = 1

    for j in range(0, len(vector)): 
        for i in range(0, max_score + 1): 
            if i == 0 and j == 0: 
                continue

            peptide_count = 0
            for weight in protein_utils.ALL_AA_WEIGHTS: 
                if j - weight < 0 or i - vector[j] < 0 or i - vector[j] > max_score: 
                    continue  
                else:
                    peptide_count += dp_table[i - vector[j]][j - weight]
            
            dp_table[i][j] = peptide_count
    
    # for line in dp_table: 
    #     universal_utils.print_arr(line)

    # print("last column")
    # universal_utils.print_arr([dp_table[i][-1] for i in range(max_score + 1)])

    peptide_count = 0
    for score in range(threshold, max_score + 1): 
        peptide_count += dp_table[score][-1]
    
    return peptide_count

    # method 2: use a dict as the dp_table
    # recursively call
    # dp_table = {(0, 0): 1}
    # return sum([get_size(score, len(vector) - 1, vector, dp_table) for score in range(threshold, max_score + 1)])

def get_size(score: int, mass: int, vector: list[int], dp_table: dict): 
    if (score, mass) in dp_table.keys(): 
        return dp_table[(score, mass)]
    if score < 0 or mass < 0: 
        dp_table[(score, mass)] = 0
        return 0
    
    size = sum([get_size(score - vector[mass], mass - AA_weight, vector, dp_table) for AA_weight in protein_utils.ALL_AA_WEIGHTS])
    dp_table[(score, mass)] = size
    return size

def calculate_probability_of_spectrum_dictionary(spectral_vector: list[int], threshold: int, max_score: int): 
    vector = [0] + spectral_vector

    dp_table = [[0 for j in range(len(vector))] for i in range(max_score + 1)]
    dp_table[0][0] = 1

    for j in range(0, len(vector)): 
        for i in range(0, max_score + 1): 
            if i == 0 and j == 0: 
                continue

            probability = 0
            for weight in protein_utils.ALL_AA_WEIGHTS: 
                if j - weight < 0 or i - vector[j] < 0 or i - vector[j] > max_score: 
                    continue
                else:
                    probability += dp_table[i - vector[j]][j - weight] / 20
            
            dp_table[i][j] = probability
    
    # for line in dp_table: 
    #     universal_utils.print_arr(line)

    # print("last column")
    # universal_utils.print_arr([dp_table[i][-1] for i in range(max_score + 1)])

    probability = 0
    for score in range(threshold, max_score + 1): 
        probability += dp_table[score][-1]
    
    return probability 

def spectral_alignment_problem(peptide: str, spectral_vector: list, modification_num: int): 
    peptide_weight = protein_utils.get_weight_by_peptide(peptide)
    vector_len = len(spectral_vector)
    dp_table = [
        [
            [
                0 for k in range(modification_num + 1)
            ] for j in range(vector_len + 1)
        ] for i in range(peptide_weight + 1)
    ]
    for k in range(0, modification_num + 1): 
        dp_table[0][0][k] = - float("inf")
        for j in range(1, vector_len + 1): 
            dp_table[0][j][k] = spectral_vector[j - 1]
        for i in range(1, peptide_weight + 1):
            dp_table[i][0][k] = - float("inf")
    dp_table[0][0][0] = 0

    for k in range(0, modification_num + 1): 
        for i in range(peptide_weight + 1): 
            temp = ""
            for j in range(vector_len + 1): 
                temp += str(dp_table[i][j][k]) + " "
            print(temp)
        print()
    
    prefix_weights = []
    for i in range(1, len(peptide) + 1): 
        prefix = peptide[0 : i]
        prefix_weights.append(protein_utils.get_weight_by_peptide(prefix))
    
    universal_utils.print_arr(prefix_weights)

    for k in range(0, modification_num + 1):
        for i in range(peptide_weight + 1): 
            for j in range(vector_len + 1): 
                max_score = - float("inf")
                

if __name__ == "__main__": 
    # 1.1 step 2
    # spectral_vector = [int(item) for item in universal_utils.parse_arr("0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8")]
    # proteome = "XZZXZXXXZXZZXZXXZ"
    # peptide = identify_peptide_from_proteome(spectral_vector, proteome)[1]
    # print(peptide)

    # with open("./datasets/dataset_30270_2.txt", "r") as f: 
    #     spectral_vector = [int(item) for item in universal_utils.parse_arr(f.readline().strip())]
    #     proteome = f.readline().strip()
    #     peptide = identify_peptide_from_proteome(spectral_vector, proteome)[0]
    #     print(peptide)

    # 1.1 step 7
#     spectral_vector_list_str = """-1 5 -4 5 3 -1 -4 5 -1 0 0 4 -1 0 1 4 4 4
# -4 2 -2 -4 4 -5 -1 4 -1 2 5 -3 -1 3 2 -3"""
#     spectral_vector_list = []
#     for line in spectral_vector_list_str.split("\n"): 
#         spectral_vector_list.append([int(item) for item in universal_utils.parse_arr(line)])
#     proteome = "XXXZXZXXZXZXXXZXXZX"
#     score_threshold = 5
#     print(search_peptide_spectrum_matches(spectral_vector_list, proteome, score_threshold))

    # with open("./datasets/dataset_30270_7.txt", "r") as f: 
    #     spectral_vector_list_str = ""
    #     line = f.readline()
    #     while True: 
    #         if not(" ") in line: 
    #             # the current line read is the proteome
    #             break
    #         spectral_vector_list_str += line
    #         line = f.readline()

    #     proteome = line.strip()
    #     score_threshold = int(f.readline().strip())

    #     spectral_vector_list = []
    #     spectral_vector_list_str = spectral_vector_list_str.strip()
    #     for line in spectral_vector_list_str.split("\n"): 
    #         spectral_vector_list.append([int(item) for item in universal_utils.parse_arr(line)])
        
    #     PSMs = search_peptide_spectrum_matches(spectral_vector_list, proteome, score_threshold)
    #     peptides = []
    #     for PSM in PSMs: 
    #         peptides.append(PSM[0])
        
    #     universal_utils.print_arr(peptides)

    # 1.3 step 3
    # spectral_vector = [int(item) for item in universal_utils.parse_arr("4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3")]
    # # spectral_vector = [int(item) for item in universal_utils.parse_arr("4 -3 -2 3 3 -4 5 -3 -4 -4 3 4 1 6")]
    # threshold = 1
    # max_score = 8
    # print(calculate_size_of_spectrum_dictionary(spectral_vector, threshold, max_score))

    # extra dataset
    # with open("./datasets/size_spectral_dictionary.txt", "r") as f: 
    #     f.readline()
    #     spectral_vector = [int(item) for item in universal_utils.parse_arr(f.readline().strip())]
    #     threshold = int(f.readline().strip())
    #     max_score = int(f.readline().strip())

    #     print(calculate_size_of_spectrum_dictionary(spectral_vector, threshold, max_score))

    # with open("./datasets/dataset_30266_3.txt", "r") as f: 
    #     spectral_vector = [int(item) for item in universal_utils.parse_arr(f.readline().strip())]
    #     threshold = int(f.readline().strip())
    #     max_score = int(f.readline().strip())

    #     print(calculate_size_of_spectrum_dictionary(spectral_vector, threshold, max_score))

    # 1.3 step 8
    # spectral_vector = [int(item) for item in universal_utils.parse_arr("4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 3")]
    # # spectral_vector = [int(item) for item in universal_utils.parse_arr("4 -3 -2 3 3 -4 5 -3 -4 -4 3 4 1 6")]
    # threshold = 1
    # max_score = 8
    # print(calculate_probability_of_spectrum_dictionary(spectral_vector, threshold, max_score)) 

    # with open("./datasets/dataset_30266_8.txt", "r") as f: 
    #     spectral_vector = [int(item) for item in universal_utils.parse_arr(f.readline().strip())]
    #     threshold = int(f.readline().strip())
    #     max_score = int(f.readline().strip()) 

    #     print(calculate_probability_of_spectrum_dictionary(spectral_vector, threshold, max_score))

    # 1.6 step 3
    peptide = "XXZ"
    spectral_vector = universal_utils.parse_arr("4 -3 -2 3 3 -4 5 -3 -1 -1 3 4 1 -1")
    modification_num = 2
    spectral_alignment_problem(peptide, spectral_vector, modification_num)

