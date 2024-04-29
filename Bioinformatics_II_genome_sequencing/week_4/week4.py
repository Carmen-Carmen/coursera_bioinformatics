import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pyperclip
import time
from collections import Counter
import random
import math

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

# dicts

RNA_codon_to_AA = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N", "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S", "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H", "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R", "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G", "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y", "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C", "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F",
}

AA_3_letter_to_abbr = {
    "Gly": "G", "Ala": "A", "Val": "V", "Leu": "L", 
    "Gly": "G", "Ala": "A", "Val": "V", "Leu": "L", 
    "Ile": "I", "Pro": "P", "Phe": "F", "Tyr": "Y", 
    "Trp": "W", "Ser": "S", "Thr": "T", "Cys": "C", 
    "Met": "M", "Asn": "N", "Gln": "Q", "Asp": "D", 
    "Glu": "E", "Lys": "K", "Arg": "R", "His": "H"
}

AA_abbr_to_3_letter = {
    "G": "Gly", "A": "Ala", "V": "Val", "L": "Leu", 
    "I": "Ile", "P": "Pro", "F": "Phe", "Y": "Tyr", 
    "W": "Trp", "S": "Ser", "T": "Thr", "C": "Cys", 
    "M": "Met", "N": "Asn", "Q": "Gln", "D": "Asp", 
    "E": "Glu", "K": "Lys", "R": "Arg", "H": "His",
}

AA_to_RNA_codon = {
    "*": ['UAA', 'UAG', 'UGA'], 
    "A": ['GCA', 'GCC', 'GCG', 'GCU'], 
    "C": ['UGC', 'UGU'], 
    "D": ['GAC', 'GAU'], 
    "E": ['GAA', 'GAG'], 
    "F": ['UUC', 'UUU'], 
    "G": ['GGA', 'GGC', 'GGG', 'GGU'], 
    "H": ['CAC', 'CAU'], 
    "I": ['AUA', 'AUC', 'AUU'], 
    "K": ['AAA', 'AAG'], 
    "L": ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 
    "M": ['AUG'], 
    "N": ['AAC', 'AAU'], 
    "P": ['CCA', 'CCC', 'CCG', 'CCU'], 
    "Q": ['CAA', 'CAG'], 
    "R": ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'], 
    "S": ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'], 
    "T": ['ACA', 'ACC', 'ACG', 'ACU'], 
    "V": ['GUA', 'GUC', 'GUG', 'GUU'], 
    "W": ['UGG'], 
    "Y": ['UAC', 'UAU'],
}

AA_weight = {
    "G": 57, 
    "A": 71, 
    "S": 87, 
    "P": 97, 
    "V": 99, 
    "T": 101, 
    "C": 103, 
    "I": 113, 
    "L": 113, 
    "N": 114, 
    "D": 115, 
    "K": 128, 
    "Q": 128, 
    "E": 129, 
    "M": 131, 
    "H": 137, 
    "F": 147, 
    "R": 156, 
    "Y": 163, 
    "W": 186
}

weight_AA = {
    57: 'G', 71: 'A', 87: 'S', 97: 'P',
    99: 'V', 101: 'T', 103: 'C', 113:'I/L',
    114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
    131: 'M', 137: 'H', 147: 'F', 156: 'R', 
    163: 'Y', 186: 'W'
}

def translate_RNA_to_peptide(RNA_strand):
    peptide_string = ""
    i = 0
    while i < len(RNA_strand):
        codon = RNA_strand[i : i + 3]
        amino_acid = RNA_codon_to_AA[codon]
        if amino_acid == "*":
            break
        peptide_string += amino_acid

        i += 3

    return peptide_string

def get_AA_abbr(AA_3_letter):
    return AA_3_letter_to_abbr[AA_3_letter]

def get_AA_3_letter(AA_abbr):
    return AA_abbr_to_3_letter[AA_abbr]

def get_weight_by_AA(AA):
    return AA_weight[AA]

def get_theoretical_spectrum_of_cyclic_peptide(peptide_strand):
    peptide_len = len(peptide_strand)
    prefix_mass = []
    for i in range(peptide_len + 1):
        prefix_mass.append(get_weight_of_peptide(peptide_strand[0 : i]))
    peptice_weight = get_weight_of_peptide(peptide_strand)

    theoretical_spectrum = [0]
    for i in range(len(prefix_mass) - 1):
        for j in range(i + 1, len(prefix_mass)):
            # masses found by LinearSpectrum 
            theoretical_spectrum.append(prefix_mass[j] - prefix_mass[i])
            # masses corresponding to subpeptides wrapping around the end of peptide
            if i > 0 and j < len(prefix_mass) - 1:
                theoretical_spectrum.append(peptice_weight - (prefix_mass[j] - prefix_mass[i]))

    theoretical_spectrum.sort()
    return theoretical_spectrum

def get_theoretical_spectrum_of_linear_peptide(peptide_strand):
    # more elegant method
    # Our approach to generating its theoretical spectrum is based on the assumption that 
    # the mass of any subpeptide is equal to the difference between the masses of two prefixes of Peptide. 
    theoretical_spectrum = [0]
    prefix_mass = []
    peptide_len = len(peptide_strand)
    for i in range(peptide_len + 1):
        prefix_mass.append(get_weight_of_peptide(peptide_strand[0 : i]))
    for i in range(peptide_len):
        for j in range(i + 1, peptide_len + 1):
            theoretical_spectrum.append(prefix_mass[j] - prefix_mass[i])
    
    theoretical_spectrum.sort()

    return theoretical_spectrum

def get_weight_of_peptide(peptide_strand):
    if len(peptide_strand) == 0:
        return 0

    peptide_mass = 0
    for AA in [*peptide_strand]:
        peptide_mass += get_weight_by_AA(AA)
    return peptide_mass

# calculate the score of a cyclic peptide against an "experimental" spectrum
# in practice, generate the theoretical spectrum of this cyclic peptide
# then calculate the score of each mass contributes from theoretical or experimental spectrum
# and add the smaller score to the total score
def score_cyclic_peptide_against_spectrum(peptide_strand, spectrum):
    theoretical_spectrum = get_theoretical_spectrum_of_cyclic_peptide(peptide_strand)
    combined_spectrum = set(theoretical_spectrum + spectrum)
    counter_theoretical = Counter(theoretical_spectrum)
    counter_experimental = Counter(spectrum)
    score = 0
    for mass in combined_spectrum:
        score += min(counter_theoretical[mass], counter_experimental[mass])
    
    return score

def score_linear_peptide_against_spectrum(peptide_strand, spectrum):
    theoretical_spectrum = get_theoretical_spectrum_of_linear_peptide(peptide_strand)
    combined_spectrum = set(theoretical_spectrum + spectrum)
    # the use of collections.Counter significantly accelerated the speed of this method, 
    # which will iterate for thousands of times in the sequencing algorithm
    # 不再需要每次调用list.count()，相当于已经生成了一张 数组元素-count数 的map
    counter_theoretical = Counter(theoretical_spectrum)
    counter_experimental = Counter(spectrum)
    score = 0
    for mass in combined_spectrum:
        score += min(counter_theoretical[mass], counter_experimental[mass])
    
    return score

def get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, N):
    AAs = list(AA_weight.keys())
    AAs.remove("I")
    AAs.remove("K")
    # 在这里就remove掉I和K是有意义的，因为这样的话刚开始leaderboard中就只会剩下18个只有1个AA的peptide

    # leaderboard should store only the N-highest scoring linear peptides
    leaderboard = AAs.copy()
    # leader_peptide stores the linear peptide with the highest score
    leader_peptides = []
    last_score_of_leader_peptide = 0

    target_peptide_weight = spectrum[-1]
    while len(leaderboard) != 0:
        to_remove = []
        to_add = []

        for subpeptide in leaderboard:
            if get_weight_of_peptide(subpeptide) == target_peptide_weight:
                score_temp = score_cyclic_peptide_against_spectrum(subpeptide, spectrum)
                if score_temp > last_score_of_leader_peptide:
                    last_score_of_leader_peptide = score_temp
                    leader_peptides = []
                    leader_peptides.append(subpeptide)
                elif score_temp == last_score_of_leader_peptide:
                    leader_peptides.append(subpeptide)
            elif get_weight_of_peptide(subpeptide) > target_peptide_weight:
                to_remove.append(subpeptide)

        for subpeptide in to_remove:
            leaderboard.remove(subpeptide)
        
        # trim_leaderboard_by_N(leaderboard, spectrum, N)
        # trim method1: too slow to use the mapping structure
        # score_to_subpeptides = {}
        # # {
        # #     score1 - [subpeptide11, subpeptide12, ...]
        # #     score2 - [subpeptide21, subpeptide22, ...]
        # #     ...
        # # }
        # for subpeptide in leaderboard:
        #     score = score_linear_peptide_against_spectrum(subpeptide, spectrum)
        #     if not (score in score_to_subpeptides.keys()):
        #         # initialize the val to the key as an array
        #         score_to_subpeptides[score] = []
        #     score_to_subpeptides[score].append(subpeptide)
        # scores = list(score_to_subpeptides.keys())
        # scores.sort()
        # scores.reverse()
        # count = 0
        # trim_index = 0
        # for score in scores:
        #     count += len(score_to_subpeptides[score])
        #     trim_index += 1
        #     if count >= N:
        #         break
        # scores = scores[: trim_index]

        # # expand
        # for score in scores:
        #     for subpeptide in score_to_subpeptides[score]:
        #         for AA in AAs:
        #             to_add.append(subpeptide + AA)
        # leaderboard = []
        # for subpeptide in to_add:
        #     leaderboard.append(subpeptide)
        
        # trim method 2, use array to store scores
        # sort scores_arr and leaderboard simultaneously
        leaderboard = trim_leaderboard(leaderboard, spectrum, N)

        # expanding
        for subpeptide in leaderboard:
            for AA in AAs:
                to_add.append(subpeptide + AA)
        leaderboard.clear()
        for subpeptide in to_add:
            leaderboard.append(subpeptide)
        
    return leader_peptides

def trim_leaderboard(leaderboard, spectrum, N):
    # start_time = time.perf_counter()

    # 1. compute scores of all peptides in the leaderboard
    scores = []
    for subpeptide in leaderboard:
        scores.append(score_linear_peptide_against_spectrum(subpeptide, spectrum))

    # end_time = time.perf_counter()
    # print("leaderboard size: %d, generating scores used %dms" %(len(scores), (end_time - start_time) * 1000))


    # start_time = time.perf_counter()

    # 2. sort the leaderboard by scores
    # implement bubble sort in simultaneously sorting scores and leaderboard
    # for i in range(len(scores)):
    #     for j in range(len(scores) - i - 1):
    #         if scores[j] < scores[j + 1]:
    #             temp_score = scores[j]
    #             scores[j] = scores[j + 1]
    #             scores[j + 1] = temp_score
    #             temp_peptide = leaderboard[j]
    #             leaderboard[j] = leaderboard[j + 1]
    #             leaderboard[j + 1] = temp_peptide
    
    # use zip() to sort the 2 list simultaneously implementing quick sort (much faster than bubble sort...)
    # 使用zip绑定2个数组
    combined = zip(leaderboard, scores)
    # 使用sorted方法（应该是快速排序法）
    # 第二个参数是按照哪个排序，lambda表达式中的x[1]表示按照第二个数组也就是score来排序
    sorted_combination = sorted(combined, key=lambda x: x[1], reverse=True)
    # 将数据从sorted_combination中取出
    scores = [item[1] for item in sorted_combination]
    leaderboard = [item[0] for item in sorted_combination]

    # end_time = time.perf_counter()
    # print("leaderboard size: %d, sorting scores used %dms" %(len(scores), (end_time - start_time) * 1000))
    # print()

    # 3. trim the leaderboard
    trim_bound = len(scores)
    if N < len(scores):
        score_index_N = scores[N]
        for i in range(N + 1, len(scores)):
            # find the first score smaller than last recorded score, 
            # and the count of peptides bigger than N
            if scores[i] < score_index_N:
                break
            trim_bound = i

    return leaderboard[: trim_bound + 1]

def transfer_peptide_to_masses_string(peptide):
    result = ""
    for AA in peptide:
        result += str(get_weight_by_AA(AA)) + "-"
    return result[: -1]   

# spectrum text: "mass1 mass2 mass3 ..."
def transfer_spectrum_text_to_array(spectrum_text):
    return [float(mass) for mass in spectrum_text.split(" ")]    

###### below are for the cyclopeptide composed of extended AAs sequecing problem ######
extended_AA_weight = {
    chr(57): 57, 
    chr(58): 58, 
    chr(59): 59, 
    chr(60): 60, 
    chr(61): 61, 
    chr(62): 62, 
    chr(63): 63, 
    chr(64): 64, 
    chr(65): 65, 
    chr(66): 66, 
    chr(67): 67, 
    chr(68): 68, 
    chr(69): 69, 
    chr(70): 70, 
    chr(71): 71, 
    chr(72): 72, 
    chr(73): 73, 
    chr(74): 74, 
    chr(75): 75, 
    chr(76): 76, 
    chr(77): 77, 
    chr(78): 78, 
    chr(79): 79, 
    chr(80): 80, 
    chr(81): 81, 
    chr(82): 82, 
    chr(83): 83, 
    chr(84): 84, 
    chr(85): 85, 
    chr(86): 86, 
    chr(87): 87, 
    chr(88): 88, 
    chr(89): 89, 
    chr(90): 90, 
    chr(91): 91, 
    chr(92): 92, 
    chr(93): 93, 
    chr(94): 94, 
    chr(95): 95, 
    chr(96): 96, 
    chr(97): 97, 
    chr(98): 98, 
    chr(99): 99, 
    chr(100): 100, 
    chr(101): 101, 
    chr(102): 102, 
    chr(103): 103, 
    chr(104): 104, 
    chr(105): 105, 
    chr(106): 106, 
    chr(107): 107, 
    chr(108): 108, 
    chr(109): 109, 
    chr(110): 110, 
    chr(111): 111, 
    chr(112): 112, 
    chr(113): 113, 
    chr(114): 114, 
    chr(115): 115, 
    chr(116): 116, 
    chr(117): 117, 
    chr(118): 118, 
    chr(119): 119, 
    chr(120): 120, 
    chr(121): 121, 
    chr(122): 122, 
    chr(123): 123, 
    chr(124): 124, 
    chr(125): 125, 
    chr(126): 126, 
    chr(127): 127, 
    chr(128): 128, 
    chr(129): 129, 
    chr(130): 130, 
    chr(131): 131, 
    chr(132): 132, 
    chr(133): 133, 
    chr(134): 134, 
    chr(135): 135, 
    chr(136): 136, 
    chr(137): 137, 
    chr(138): 138, 
    chr(139): 139, 
    chr(140): 140, 
    chr(141): 141, 
    chr(142): 142, 
    chr(143): 143, 
    chr(144): 144, 
    chr(145): 145, 
    chr(146): 146, 
    chr(147): 147, 
    chr(148): 148, 
    chr(149): 149, 
    chr(150): 150, 
    chr(151): 151, 
    chr(152): 152, 
    chr(153): 153, 
    chr(154): 154, 
    chr(155): 155, 
    chr(156): 156, 
    chr(157): 157, 
    chr(158): 158, 
    chr(159): 159, 
    chr(160): 160, 
    chr(161): 161, 
    chr(162): 162, 
    chr(163): 163, 
    chr(164): 164, 
    chr(165): 165, 
    chr(166): 166, 
    chr(167): 167, 
    chr(168): 168, 
    chr(169): 169, 
    chr(170): 170, 
    chr(171): 171, 
    chr(172): 172, 
    chr(173): 173, 
    chr(174): 174, 
    chr(175): 175, 
    chr(176): 176, 
    chr(177): 177, 
    chr(178): 178, 
    chr(179): 179, 
    chr(180): 180, 
    chr(181): 181, 
    chr(182): 182, 
    chr(183): 183, 
    chr(184): 184, 
    chr(185): 185, 
    chr(186): 186, 
    chr(187): 187, 
    chr(188): 188, 
    chr(189): 189, 
    chr(190): 190, 
    chr(191): 191, 
    chr(192): 192, 
    chr(193): 193, 
    chr(194): 194, 
    chr(195): 195, 
    chr(196): 196, 
    chr(197): 197, 
    chr(198): 198, 
    chr(199): 199, 
    chr(200): 200,
}

def get_weight_of_peptide_extended_AAs(peptide_string):
    mass = 0
    for AA in list(peptide_string):
        mass += extended_AA_weight[AA]
    return mass

def get_theoretical_spectrum_of_cyclic_peptide_extended_AAs(peptide_strand):
    peptide_len = len(peptide_strand)
    prefix_mass = []
    for i in range(peptide_len + 1):
        prefix_mass.append(get_weight_of_peptide_extended_AAs(peptide_strand[0 : i]))
    peptide_weight = get_weight_of_peptide_extended_AAs(peptide_strand)

    theoretical_spectrum = [0]
    for i in range(len(prefix_mass) - 1):
        for j in range(i + 1, len(prefix_mass)):
            theoretical_spectrum.append(prefix_mass[j] - prefix_mass[i])
            if i > 0 and j < len(prefix_mass) - 1:
                theoretical_spectrum.append(peptide_weight - (prefix_mass[j] - prefix_mass[i]))

    theoretical_spectrum.sort()
    return theoretical_spectrum

def get_theoretical_spectrum_of_linear_peptide_extended_AAs(peptide_strand):
    theoretical_spectrum = [0]
    prefix_mass = []
    peptide_len = len(peptide_strand)
    for i in range(peptide_len + 1):
        prefix_mass.append(get_weight_of_peptide_extended_AAs(peptide_strand[0 : i]))
    for i in range(peptide_len):
        for j in range(i + 1, peptide_len + 1):
            theoretical_spectrum.append(prefix_mass[j] - prefix_mass[i])
    
    theoretical_spectrum.sort()

    return theoretical_spectrum

def score_cyclic_peptide_against_spectrum_extended_AAs(peptide_strand, spectrum):
    theoretical_spectrum = get_theoretical_spectrum_of_cyclic_peptide_extended_AAs(peptide_strand)
    combined_spectrum = set(theoretical_spectrum + spectrum)
    counter_theoretical = Counter(theoretical_spectrum)
    counter_experimental = Counter(spectrum)
    score = 0
    for mass in combined_spectrum:
        score += min(counter_theoretical[mass], counter_experimental[mass])
    
    return score

def score_linear_peptide_against_spectrum_extended_AAs(peptide_strand, spectrum):
    theoretical_spectrum = get_theoretical_spectrum_of_linear_peptide_extended_AAs(peptide_strand)
    combined_spectrum = set(theoretical_spectrum + spectrum)
    counter_theoretical = Counter(theoretical_spectrum)
    counter_experimental = Counter(spectrum)
    score = 0
    for mass in combined_spectrum:
        score += min(counter_theoretical[mass], counter_experimental[mass])
    
    return score

def trim_leaderboard_extended_AAs(leaderboard, spectrum, N):
    start_time = time.perf_counter()

    scores = []
    for subpeptide in leaderboard:
        scores.append(score_linear_peptide_against_spectrum_extended_AAs(subpeptide, spectrum))

    combined = zip(leaderboard, scores)
    sorted_combination = sorted(combined, key=lambda x: x[1], reverse=True)
    scores = [item[1] for item in sorted_combination]
    leaderboard = [item[0] for item in sorted_combination]

    end_time = time.perf_counter()
    print("leaderboard size: %d, generating and sorting scores used %dms" \
          %(len(scores), (end_time - start_time) * 1000))

    trim_bound = len(scores)
    if N < len(scores):
        score_index_N = scores[N]
        for i in range(N + 1, len(scores)):
            if scores[i] < score_index_N:
                break
            trim_bound = i

    print("trim bound:", trim_bound)
    print()

    return leaderboard[: trim_bound + 1]

# this method returns peptides in the form of masses string
def get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs(spectrum, N):
    AAs = list(extended_AA_weight.keys())

    # leaderboard should store only the N-highest scoring linear peptides
    leaderboard = AAs.copy()
    # leader_peptide stores the linear peptide with the highest score
    leader_peptides = []
    last_score_of_leader_peptide = 0

    target_peptide_weight = spectrum[-1]
    while len(leaderboard) != 0:
        to_remove = []
        to_add = []

        for subpeptide in leaderboard:
            if get_weight_of_peptide_extended_AAs(subpeptide) == target_peptide_weight:
                score_temp = score_cyclic_peptide_against_spectrum_extended_AAs(subpeptide, spectrum)
                if score_temp > last_score_of_leader_peptide:
                    last_score_of_leader_peptide = score_temp
                    leader_peptides = []
                    leader_peptides.append(subpeptide)
                elif score_temp == last_score_of_leader_peptide:
                    leader_peptides.append(subpeptide)
            elif get_weight_of_peptide_extended_AAs(subpeptide) > target_peptide_weight:
                to_remove.append(subpeptide)

        for subpeptide in to_remove:
            leaderboard.remove(subpeptide)
        
        # trim_leaderboard_by_N(leaderboard, spectrum, N)        
        leaderboard = trim_leaderboard_extended_AAs(leaderboard, spectrum, N)

        # expanding
        for subpeptide in leaderboard:
            for AA in AAs:
                to_add.append(subpeptide + AA)
        leaderboard.clear()
        for subpeptide in to_add:
            leaderboard.append(subpeptide)
    
    results = []
    for peptide in leader_peptides:
        temp = ""
        for AA in list(peptide):
            temp += str(extended_AA_weight[AA]) + "-"
        temp = temp[: -1]
        results.append(temp)

    return results

###### end of this problem ######

# spectral convultion problem
# input: a collection of integers in increasing order
# output: the list of elements in the convulsion of this spectrum. 
#       if an element has multiplicity k, it should appear exactly k times
def get_convolution_of_spectrum(spectrum):
    spectrum.sort()
    convulsion_list = []
    for row in range(len(spectrum)):
        for col in range(row):
            val = spectrum[row] - spectrum[col]
            if val != 0:
                convulsion_list.append(val)
    
    return convulsion_list

# A cyclic peptide LeaderPeptide with amino acids taken only 
# from the top M elements (and ties) of the convolution of Spectrum that fall between 57 and 200, 
# and where the size of Leaderboard is restricted to the top N (and ties).
def get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution(spectrum, N, M):
    # get the top M elements (and ties) of the convolution of spectrum
    # which fall between range(57, 201)
    convolution_list = get_convolution_of_spectrum(spectrum)
    AA_mass_list = [mass for mass in convolution_list if mass >= 57 and mass <= 200]
    AA_mass_set = list(set(AA_mass_list))
    AA_mass_count = []
    counter_item = Counter(AA_mass_list)
    for mass in AA_mass_set:
        AA_mass_count.append(counter_item[mass])
    combined = zip(AA_mass_set, AA_mass_count)
    sorted_combination = sorted(combined, key=lambda x:x[1], reverse=True)
    AA_mass_set = [item[0] for item in sorted_combination]
    AA_mass_count = [item[1] for item in sorted_combination]
    trim_bound = len(AA_mass_count)
    if M < len(AA_mass_count):
        count_index_M = AA_mass_count[M]
        for i in range(M + 1, len(AA_mass_count)):
            if AA_mass_count[i] < count_index_M:
                break
            trim_bound = i
    AA_mass_set = AA_mass_set[: trim_bound + 1]

    AAs = []
    for AA in extended_AA_weight:
        if extended_AA_weight[AA] in AA_mass_set:
            AAs.append(AA)
    
    # leaderboard should store only the N-highest scoring linear peptides
    leaderboard = AAs.copy()
    # leader_peptide stores the linear peptide with the highest score
    leader_peptides = []
    last_score_of_leader_peptide = 0

    target_peptide_weight = spectrum[-1]
    target_peptide_weight = 1322
    while len(leaderboard) != 0:
        to_remove = []
        to_add = []

        for subpeptide in leaderboard:
            if get_weight_of_peptide_extended_AAs(subpeptide) == target_peptide_weight:
                score_temp = score_cyclic_peptide_against_spectrum_extended_AAs(subpeptide, spectrum)
                if score_temp > last_score_of_leader_peptide:
                    last_score_of_leader_peptide = score_temp
                    leader_peptides = []
                    leader_peptides.append(subpeptide)
                elif score_temp == last_score_of_leader_peptide:
                    leader_peptides.append(subpeptide)
            elif get_weight_of_peptide_extended_AAs(subpeptide) > target_peptide_weight:
                to_remove.append(subpeptide)

        for subpeptide in to_remove:
            leaderboard.remove(subpeptide)
        
        # trim_leaderboard_by_N(leaderboard, spectrum, N)        
        leaderboard = trim_leaderboard_extended_AAs(leaderboard, spectrum, N)

        # expanding
        for subpeptide in leaderboard:
            for AA in AAs:
                to_add.append(subpeptide + AA)
        leaderboard.clear()
        for subpeptide in to_add:
            leaderboard.append(subpeptide)
    
    # format the peptides into masses_strings
    results = []
    for peptide in leader_peptides:
        temp = ""
        for AA in list(peptide):
            temp += str(extended_AA_weight[AA]) + "-"
        temp = temp[: -1]
        results.append(temp)

    print("highest score:", last_score_of_leader_peptide)
    return results

###### for the last problem of tyrocidine experimental spectrum sequencing ######
tyrocidine_spectrum = transfer_spectrum_text_to_array("0 97 99 113 114 128 128 147 147 163 186 227 241 242 244 260 261 262 283 291 333 340 357 388 389 390 390 405 430 430 447 485 487 503 504 518 543 544 552 575 577 584 631 632 650 651 671 672 690 691 738 745 747 770 778 779 804 818 819 835 837 875 892 892 917 932 932 933 934 965 982 989 1031 1039 1060 1061 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1225 1322")

def trim_leaderboard_extended_AAs_for_tyrocydine(leaderboard, spectrum, N):
    # spectrum = tyrocidine_spectrum

    start_time = time.perf_counter()

    scores = []
    # change the score to the sum of the score to spectrum and tyrocidine spectrum
    for subpeptide in leaderboard:
        scores.append(score_linear_peptide_against_spectrum_extended_AAs(subpeptide, spectrum) \
                      + score_linear_peptide_against_spectrum_extended_AAs(subpeptide, tyrocidine_spectrum))

    combined = zip(leaderboard, scores)
    sorted_combination = sorted(combined, key=lambda x: x[1], reverse=True)
    scores = [item[1] for item in sorted_combination]
    leaderboard = [item[0] for item in sorted_combination]

    end_time = time.perf_counter()
    print("leaderboard size: %d, generating and sorting scores used %dms" \
          %(len(scores), (end_time - start_time) * 1000))

    trim_bound = len(scores)
    if N < len(scores):
        score_index_N = scores[N]
        for i in range(N + 1, len(scores)):
            if scores[i] < score_index_N:
                break
            trim_bound = i

    print("trim bound:", trim_bound)
    print()

    return leaderboard[: trim_bound + 1]

def get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution_for_tyrocydine(spectrum, N, M):
    AAs = []
    for mass in "97 99 113 114 128 147 163 186".split(" "):
        for AA in extended_AA_weight.keys():
            if extended_AA_weight[AA] == int(mass):
                AAs.append(AA)

    # leaderboard should store only the N-highest scoring linear peptides
    leaderboard = AAs.copy()
    # leader_peptide stores the linear peptide with the highest score
    leader_peptides = []
    last_score_of_leader_peptide = 0

    target_peptide_weight = spectrum[-1]
    target_peptide_weight = 1322
    while len(leaderboard) != 0:
        to_remove = []
        to_add = []

        for subpeptide in leaderboard:
            # if get_weight_of_peptide_extended_AAs(subpeptide) == target_peptide_weight:
            if len(subpeptide) == 10:
                score_temp = score_cyclic_peptide_against_spectrum_extended_AAs(subpeptide, spectrum)
                            # score_cyclic_peptide_against_spectrum_extended_AAs(subpeptide, tyrocidine_spectrum)
                if score_temp >= 34:
                    print(score_temp)
                    leader_peptides.append(subpeptide)
            elif len(subpeptide) > 10:
                to_remove.append(subpeptide)

        for subpeptide in to_remove:
            leaderboard.remove(subpeptide)
        
        # trim_leaderboard_by_N(leaderboard, spectrum, N)        
        leaderboard = trim_leaderboard_extended_AAs_for_tyrocydine(leaderboard, spectrum, N)

        # expanding
        for subpeptide in leaderboard:
            for AA in AAs:
                to_add.append(subpeptide + AA)
        leaderboard.clear()
        for subpeptide in to_add:
            leaderboard.append(subpeptide)
    
    # format the peptides into masses_strings
    results = []
    for peptide in leader_peptides:
        temp = ""
        for AA in list(peptide):
            temp += str(extended_AA_weight[AA]) + " "
        temp = temp[: -1]
        results.append(temp)

    return results

if __name__ == "__main__":
    # 1.1 step 3
    # spectrum = [int(mass) for mass in "0 99 113 114 128 227 257 299 355 356 370 371 484".split(" ")]
    # print(score_linear_peptide_against_spectrum("NQEL", spectrum))
    # with open("dataset_102_3.txt", "r") as f:
    #     peptide_strand = f.readline().strip()
    #     spectrum = [int(mass) for mass in f.readline().strip().split(" ")]
    #     print("score:", score_cyclic_peptide_against_spectrum(peptide_strand, spectrum))

    # 1.5 step 1
    # with open("dataset_4913_1.txt", "r") as f:
    #     peptide_strand = f.readline().strip()
    #     spectrum = [int(mass) for mass in f.readline().strip().split(" ")]
    #     print("score:", score_linear_peptide_against_spectrum(peptide_strand, spectrum))

    # 1.5 step 3
    # leaderboard = "LAST ALST TLLT TQAS".split(" ")
    # spectrum = [int(mass) for mass in "0 71 87 101 113 158 184 188 259 271 372".split(" ")]
    # print_arr(trim_leaderboard(leaderboard, spectrum, 2))
    # with open("dataset_4913_3.txt", "r") as f:
    #     leaderboard = f.readline().strip().split(" ")
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())
    #     N = int(f.readline().strip())
    #     print_arr(trim_leaderboard(leaderboard, spectrum, N))

    # 1.1 step 8
    # spectrum = [int(mass) for mass in "0 71 113 129 147 200 218 260 313 331 347 389 460".split(" ")]
    # print_arr(list(set(
    #     [transfer_peptide_to_masses_string(peptide) \
    #      for peptide in get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, 10)]
    #     )))

    # extra dataset
    # with open("./leaderboard_cyclopeptide_sequencing.txt", "r") as f:
    #     f.readline()
    #     N = int(f.readline().strip())
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())

    #     print_arr(list(set(
    #         [transfer_peptide_to_masses_string(peptide) \
    #         for peptide in get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, N)]
    #         )))

    # with open("./dataset_102_8.txt", "r") as f:
    #     N = int(f.readline().strip())
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())
    #     print_arr(list(set(
    #         [transfer_peptide_to_masses_string(peptide) \
    #         for peptide in get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, N)]
    #         )))

    # 1.1 step 9
    # with open("./Tyrocidine_B1_Spectrum_10.txt", "r") as f:
    #     spectrum = transfer_spectrum_text_to_array(f.read().strip())
    #     possible_peptides = get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, 1000)
    #     scores = [score_cyclic_peptide_against_spectrum(peptide, spectrum) for peptide in possible_peptides]
    #     print("%s has highest score of %d" %(possible_peptides[0], scores[0]))

    #     print(get_weight_of_peptide("VYKNFWPFIK"))

    # 1.1 step 10
    # it is not the order that matters, it is because the starting leaderboard should
    # only contain 18 subpeptides (single AA) rather than 20, since I/L and K/Q share the same mass
    # with open("./dataset_102_10.txt", "r") as f:
    #     N = int(f.readline().strip())
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())
    #     # possible_peptides = list(set(
    #     #     get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, N)
    #     # ))
    #     possible_peptides = get_cyclopeptide_by_theoretical_spectrum_N_highest_score(spectrum, N)
    #     print_arr(possible_peptides)
    #     print("score of the leader peptides:", \
    #           score_cyclic_peptide_against_spectrum(possible_peptides[0], spectrum))
    #     # from the comments:
    #     # I ordered the output peptides in linear score order, from less to higher. It worked!!!.
    #     # scores = [score_linear_peptide_against_spectrum(peptide, spectrum) for peptide in possible_peptides]
    #     # combined = zip(possible_peptides, scores)
    #     # print(scores)
    #     # sorted_combination = sorted(combined, key=lambda x: x[1])
    #     # possible_peptides = [item[0] for item in sorted_combination]
    #     # scores = [item[1] for item in sorted_combination]
    #     # print(scores)

    #     # print(len(possible_peptides))
    #     # masses_strs = list(set(
    #     #     [transfer_peptide_to_masses_string(peptide) for peptide in possible_peptides]
    #     # ))
    #     masses_strs = [transfer_peptide_to_masses_string(peptide) for peptide in possible_peptides]

    #     # print_arr(masses_strs)
    #     for str in masses_strs:
    #         print(str)
    #     print("num of leader peptides:", len(masses_strs))

    # 1.2 step 2
    # for i in range(57, 201):
    #     print("chr(%d): %d, " %(i, i))
    # print(extended_AA_weight)
    # test_peptide_extended = ""
    # mass_of_test_peptide = 0
    # for _ in range(100):
    #     mass = int(random.random() * (201 - 57)) + 57
    #     test_peptide_extended += chr(mass)
    #     mass_of_test_peptide += mass
    # print(get_weight_of_peptide_extended_AAs(test_peptide_extended) == mass_of_test_peptide)

    # with open("./dataset_103_2.txt", "r") as f:
    #     N = int(f.readline().strip())
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())

    #     possible_peptides = get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs(spectrum, N)
    #     print("num:", len(possible_peptides))
    #     for peptide in possible_peptides:
    #         print(peptide)

    # 1.3 step 4
    # print_arr(get_convulsion_of_spectrum(transfer_spectrum_text_to_array("0 137 186 323")))
    # with open("./dataset_104_4.txt", "r") as f:
    #     spectrum = transfer_spectrum_text_to_array(f.read().strip())
    #     print_arr(get_convolution_of_spectrum(spectrum))

    # 1.3 step 7
    # M = 20
    # N = 60
    # spectrum = transfer_spectrum_text_to_array("57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493")
    # print_arr(get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution(spectrum, N, M))
    # with open("./dataset_104_7.txt", "r") as f:
    #     M = int(f.readline().strip())
    #     N = int(f.readline().strip())
    #     spectrum = transfer_spectrum_text_to_array(f.readline().strip())

    #     print_arr(get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution(spectrum, N, M))

    # 1.3 step 8
    with open("./dataset_104_8.txt", "r") as f:
        M = int(f.readline().strip())
        N = int(f.readline().strip())
        spectrum = transfer_spectrum_text_to_array(f.readline().strip())
        possible_peptides = get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution(spectrum, N, M)

        print("num:", len(possible_peptides))
        print_arr(possible_peptides)

        # peptide = possible_peptides[1]
        # peptide_in_AA_alphabet = ""
        # for mass_str in peptide.split("-"):
        #     for AA_alph in [AA for AA in AA_weight.keys() if AA != "I" and AA != "K"]:
        #         if get_weight_by_AA(AA_alph) == int(mass_str):
        #             peptide_in_AA_alphabet += AA_alph
        # print(peptide_in_AA_alphabet)

    # 1.4 step 5
    # with open("./real_spectrum.txt", "r") as f:
    #     spectrum = transfer_spectrum_text_to_array(f.read().strip())
    #     print_arr(spectrum)
    #     rounded_spectrum = []
    #     for mass in spectrum:
    #         # mass = (mass - 1.007) * (1 - 0.0004522)
    #         # rounded_spectrum.append(round(mass))

    #         # rounded_spectrum.append(round(mass - 1 - 0.3))
    #         # rounded_spectrum.append(round(mass - 1 + 0.3))

    #         rounded_spectrum.append(math.floor(mass) - 1)
    #     print_arr(rounded_spectrum)

    #     M = 20
    #     N = 1000
    #     possible_peptides = get_cyclopeptide_by_theoretical_spectrum_N_highest_score_extended_AAs_convolution_for_tyrocydine(rounded_spectrum, N, M)

    #     print("num:", len(possible_peptides))
    #     for peptide in possible_peptides:
    #         print(peptide)

    # week 4 Quiz
    # problem 3
    # print("### problem 3 ###")
    # peptide = "MAMA"
    # spectrum = transfer_spectrum_text_to_array("0 71 178 202 202 202 333 333 333 404 507 507")
    # print(score_cyclic_peptide_against_spectrum(peptide, spectrum))

    # print("\n### problem 4 ###")
    # peptide = "PEEP"
    # spectrum = transfer_spectrum_text_to_array("0 97 129 129 129 194 226 323 323 355 452")
    # print(score_linear_peptide_against_spectrum(peptide, spectrum))

    # print("\n### problem 5 ###")
    # spectrum = transfer_spectrum_text_to_array("0 86 160 234 308 320 382")
    # convolution_list = get_convolution_of_spectrum(spectrum)
    # counter = Counter(convolution_list)
    # largest_multiplicity = 0
    # elements_with_largest_multiplicity = []
    # print(counter)
    # for element in counter.keys():
    #     count = counter[element]
    #     if count > largest_multiplicity:
    #         largest_multiplicity = count
    #         elements_with_largest_multiplicity = []
    #         elements_with_largest_multiplicity.append(element)
    #     elif count == largest_multiplicity:
    #         elements_with_largest_multiplicity.append(element)
    
    # print_arr(elements_with_largest_multiplicity)