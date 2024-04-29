import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pyperclip

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

# a DNA string Pattern encodes an amino acid string Peptide 
# if the RNA string transcribed from either Pattern 
# or its reverse complement Pattern translates into Peptide
# For example, the DNA string GAAACT is transcribed into GAAACU and translated into ET. 
# The reverse complement of this DNA string, AGTTTC, is transcribed into AGUUUC and translated into SF. 
# Thus, GAAACT encodes both ET and SF.
def get_substrings_which_encode_peptide(DNA_strand, peptide_strand):
    substring_len = len(peptide_strand) * 3
    results = []
    for i in range(len(DNA_strand) - substring_len + 1):
        substring_DNA = DNA_strand[i : i + substring_len]
        RC_substring_DNA = get_reverse_complement_DNA_strand(substring_DNA)

        substring_transcribed_RNA = substring_DNA.replace("T", "U")
        RC_substring_transcribed_RNA = RC_substring_DNA.replace("T", "U")

        if translate_RNA_to_peptide(substring_transcribed_RNA) == peptide_strand or \
            translate_RNA_to_peptide(RC_substring_transcribed_RNA) == peptide_strand:
            results.append(substring_DNA)
    
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

# The theoretical spectrum of a cyclic peptide Peptide, denoted Cyclospectrum(Peptide), 
# is the collection of all of the masses of its subpeptides, 
# in addition to the mass 0 and the mass of the entire peptide, 
# with masses ordered from smallest to largest.
def get_theoretical_spectrum_of_cyclic_peptide(peptide_strand):
    # my method, simply adding mass of each subpeptide into the spectrum
    # theoretical_spectrum = []
    # peptide_len = len(peptide_strand)
    # for _ in range(peptide_len):
    #     for i in range(1, peptide_len):
    #         subpeptide = peptide_strand[0 : i]
    #         theoretical_spectrum.append(get_weight_of_peptide(subpeptide))
        
    #     # change the starting AA of the cyclic peptide
    #     first_AA = peptide_strand[0]
    #     peptide_strand = peptide_strand[1 :]
    #     peptide_strand += first_AA

    # theoretical_spectrum.append(0)
    # theoretical_spectrum.append(get_weight_of_peptide(peptide_strand))
    # theoretical_spectrum.sort()

    # return theoretical_spectrum

    # more elegant method
    # ### though fully not understood ###
    # If Peptide represents a cyclic peptide instead, 
    # then the masses in its theoretical spectrum can be divided into those found by LinearSpectrum 
    # and those corresponding to subpeptides wrapping around the end of Peptide. 
    # Furthermore, each such subpeptide has mass equal to 
    # the difference between Mass(Peptide) and a subpeptide mass identified by LinearSpectrum.
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
            # e.g. when peptide_strand = "NQEL"
            # Mass(LN)  = Mass(NQEL) - Mass(QE)
            #           = Mass(NQEL) - (Mass(NQE) - Mass(Q))
            if i > 0 and j < len(prefix_mass) - 1:
                theoretical_spectrum.append(peptice_weight - (prefix_mass[j] - prefix_mass[i]))

    theoretical_spectrum.sort()
    return theoretical_spectrum

# the theoretical spectrum of a linear peptide
def get_theoretical_spectrum_of_peptide(peptide_strand):
    # more elegant method
    # Our approach to generating its theoretical spectrum is based on the assumption that 
    # the mass of any subpeptide is equal to the difference between the masses of two prefixes of Peptide. 
    # e.g., for Peptide = NQEL, PrefixMass = (0, 114, 242, 371, 484).
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

# brute force cyclopeptide sequencing 
# the number of peptides which weigh mass is too large, so the brute force algorithm is impractical
def get_cyclopeptide_by_theoretical_spectrum_BF(spectrum):
    # the last mass in the theoretical spectrum is the largest
    # and is the mass of the entire peptide
    peptide_mass = spectrum[-1]

    results = []
    for peptide in enumerate_all_peptides_by_mass(peptide_mass):
        if get_theoretical_spectrum_of_cyclic_peptide(peptide) == spectrum:
            results.append(peptide)
    
    return results

# generate all possible peptides whose mass is equal to given mass
# the number of peptides is exponential to given mass, thus impractical
def enumerate_all_peptides_by_mass(mass):
    pass

# counting peptides with given mass, I/L and K/Q are considered the same
# implementation of recursion in dynamic programming
# mass is the goal mass to achieve in each turn of recursive call
# mass_map records the mapping of current goal mass and the count number of peptides which sum up to the mass
def count_peptides_by_mass(mass, mass_map):
    # terminus of the recursion
    if mass == 0:
        # 原来python的方法可以有2个返回值啊
        # 实际上这个方法返回了一个tuple
        # 取的时候可以用角标，也可以用ret1, ret2 = method(arg)

        # only one possibility
        return 1, mass_map
    
    if mass < 57:
        # no weight of AA is smaller than 57
        return 0, mass_map
    
    if mass in mass_map:
        return mass_map[mass], mass_map
    
    # process of recursion
    count = 0
    for i in weight_AA.keys():
        k, mass_map = count_peptides_by_mass(mass - i, mass_map)
        count += k
    mass_map[mass] = count
    return count, mass_map

# definition of an exponential function for fitting using numpy
def exponential_func(x, k, C):
    return k * np.power(C, x)

# definition of y = a*x + b
def linear_func(x, a, b):
    return a * x + b

# branch-and-bound algorigthm for cyclopeptide sequencing
# but this algorithm also has not been proven to be polynominal, 
# thus as inefficient as the brute force method from the perspective of theoretical computer science
def get_cyclopeptide_by_theoretical_spectrum_branch_and_bound(spectrum):
    AAs = list(AA_weight.keys())
    subpeptides = AAs.copy()
    results = []
    while len(subpeptides) != 0:
        to_remove = []
        to_add = []
        for subpeptide in subpeptides:
            is_subpeptide = True
            # bound
            # when comparing whether the subpeptide belongs to spectrum, use the th-spectrum of linear peptide
            for mass in get_theoretical_spectrum_of_peptide(subpeptide):
                if not (mass in spectrum):
                    is_subpeptide = False
                    break

            if is_subpeptide:
                # when comparing whether the subpeptide is identical to the spectrum, use the th-spectrum of cyclic peptide
                if get_theoretical_spectrum_of_cyclic_peptide(subpeptide) == spectrum:
                    # terminus
                    results.append(subpeptide)
                else:
                    # branch
                    for AA in AAs:
                        # less branches, because if weight(AA) is not in spectrum, 
                        # then the expanded subpeptide will never be a subpeptide of the target peptide
                        if get_weight_by_AA(AA) in spectrum:
                            to_add.append(subpeptide + AA)
            to_remove.append(subpeptide)
        
        # remove all the subpeptides in the last loop
        for subpeptide in to_remove:
            subpeptides.remove(subpeptide)
        
        # add all the expanded subpeptides
        for subpeptide in to_add:
            subpeptides.append(subpeptide)
    
    return results
# transfer the peptide strand into "w0-w1-w2-w3"-like string
def transfer_peptide_to_masses_string(peptide):
    result = ""
    for AA in peptide:
        result += str(get_weight_by_AA(AA)) + "-"
    return result[: -1]

if __name__ == "__main__":
    # for codon in RNA_codon_to_AA:
    #     AA = RNA_codon_to_AA[codon]
    #     if not (AA in AA_to_RNA_codon.keys()):
    #         AA_to_RNA_codon[AA] = []
    #     AA_to_RNA_codon[AA].append(codon)
    
    # AAs = list(AA_to_RNA_codon.keys())
    # AAs.sort()
    # for AA in AAs:
    #     print("\"%s\": %s, " %(AA, str(AA_to_RNA_codon[AA])))

    # 1.2 step 4
    # with open("dataset_96_4.txt", "r") as f:
    #     RNA_strand = f.read().strip()

    #     print(translate_RNA_to_peptide(RNA_strand))

    # 1.2 step 5
    # count = 1
    # for AA in "Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr".split("-"):
    #     AA_abbr = get_AA_abbr(AA)
    #     codon_list = AA_to_RNA_codon[AA_abbr]
    #     count *= len(codon_list)
    
    # print("num of strings: %d" %count)

    # 1.2 step 7
    # print(get_reverse_complement_DNA_strand("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"))
    # print_arr(get_substrings_which_encode_peptide("ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA", "MA"))
    # with open("dataset_96_7.txt", "r") as f:
    #     DNA_strand = f.readline().strip()
    #     peptide_strand = f.readline().strip()
    #     print_arr(get_substrings_which_encode_peptide(DNA_strand, peptide_strand))

    # 1.2 step 8
    # with open("./Bacillus_brevis.txt", "r") as f:
    #     Bacillus_brevis_genome = f.read().strip().replace("\n", "")
    #     # print(Bacillus_brevis_genome)
    #     peptide_strand = ""
    #     for AA in "Val-Lys-Leu-Phe-Pro-Trp-Phe-Asn-Gln-Tyr".split("-"): 
    #         peptide_strand += get_AA_abbr(AA)
        
    #     # tyrocidines are actually cyclic peptides, so the linear representaton of
    #     # peptide_strand should start at each of the AAs in the strand
    #     peptide_strands = []
    #     for _ in range(len(peptide_strand)):
    #         peptide_strands.append(peptide_strand)

    #         first_AA = peptide_strand[0]
    #         peptide_strand = peptide_strand[1 :]
    #         peptide_strand += first_AA
        
    #     for peptide_strand in peptide_strands:
    #         print("%s, count: %d" \
    #               %(peptide_strand, len(get_substrings_which_encode_peptide(Bacillus_brevis_genome, peptide_strand))))
        
        # results: 
        # VKLFPWFNQY, count: 0
        # KLFPWFNQYV, count: 0
        # LFPWFNQYVK, count: 0
        # FPWFNQYVKL, count: 0
        # PWFNQYVKLF, count: 0
        # WFNQYVKLFP, count: 0
        # FNQYVKLFPW, count: 0
        # NQYVKLFPWF, count: 0
        # QYVKLFPWFN, count: 0
        # YVKLFPWFNQ, count: 0

    # 1.4 step 1
    # with open("integer_mass_table.txt", "r") as f:
    #     while True:
    #         line = f.readline().strip()
    #         if line == "":
    #             break

    #         AA = line.split(" ")[0]
    #         weight = line.split(" ")[1]

    #         print("\"%s\": %s, " %(AA, weight))

    # 1.4 step 3
    # the Cyclopeptide Sequencing Problem
    # with open("dataset_98_3.txt", "r") as f:
    #     # peptide_strand = "NQEL"
    #     peptide_len = int(f.read().strip())
    #     subpeptide_num = 0
    #     for _ in range(peptide_len):
    #         # first cycle means starting from each of the AA in the peptide
    #         for i in range(1, peptide_len):
    #             # print(peptide_strand[0 : i])
    #             subpeptide_num += 1
            
    #         # first_AA = peptide_strand[0]
    #         # peptide_strand = peptide_strand[1 :]
    #         # peptide_strand += first_AA

    #     print(subpeptide_num)

    # 1.4 step 4
    # with open("dataset_98_4.txt", "r") as f:
    #     peptide_strand = f.read().strip()
    #     peptide_len = len(peptide_strand)
    #     theoretical_spectrum = get_theoretical_spectrum_of_cyclic_peptide(peptide_strand)
    #     print_arr(theoretical_spectrum)

    #     x = range(len(theoretical_spectrum))
    #     y = theoretical_spectrum
    #     bar = plt.bar(x, y, width=0.5, fc="black")
    #     plt.show()

    # 1.7 step 2
    # peptide_strand = "NQEL"
    # print_arr(get_theoretical_spectrum_of_peptide(peptide_strand))
    # with open("dataset_4912_2.txt", "r") as f:
    #     peptide_strand = f.read().strip()
    #     theoretical_spectrum = get_theoretical_spectrum_of_peptide(peptide_strand)
    #     print_arr(theoretical_spectrum)
    
    #     x = range(len(theoretical_spectrum))
    #     y = theoretical_spectrum
    #     bar = plt.bar(x, y, width=0.5, fc="black")
    #     plt.show()

    # 1.5 step 2
    # with open("dataset_99_2.txt", "r") as f:
    #     mass = int(f.read().strip())
    #     count, mass_map = count_peptides_by_mass(mass, {})
    #     print(count)

    # 1.5 step 3
    x_list = []
    y_list = []
    
    for mass in range(1, 5000):
        count = count_peptides_by_mass(mass, {})[0]
        if count != 0:
            x_list.append(mass)
            y_list.append(count)
    print("All data generated!")

    x_data = np.array(x_list)
    y_data = np.array(y_list)
    # change the type of int into float
    # so that numpy can process these real big numbers
    # https://stackoverflow.com/questions/59297543/why-do-i-get-the-loop-of-ufunc-does-not-support-argument-0-of-type-int-error-f
    y_data = y_data.astype(float)

    # # method 1 using exponential function
    # # encountered error:
    # # RuntimeWarning: overflow encountered in multiply
    # # params, _ = curve_fit(exponential_func, x_data, y_data)
    # # k_fit = params[0]
    # # C_fit = params[1]

    # # method 2 using linear function fitting
    # numpy can process all the elements in an array, 
    # simply by passing the array as a param into any calculation methods
    y_data_log = np.log(y_data)
    # y = a * x + b
    # perform curve fitting
    params, _ = curve_fit(linear_func, x_data, y_data_log)
    # extract the fitted params
    a_fit = params[0]
    b_fit = params[1]
    # y = k * C^x
    # ln(y) = ln(k * C^x) = ln(k) + x * ln(C)
    # ln(C) = a_fit, C = e^a_fit
    C_fit = np.power(np.e, a_fit)
    k_fit = np.power(np.e, b_fit)
    
    print("k:", k_fit, "C:", C_fit)
    plt.plot(x_data, y_data)
    plt.show()

    # 1.6 step 3
    # peptide_strand = ""
    # for i in range(19105):
    #     peptide_strand += "A"
    # print(len(get_theoretical_spectrum_of_peptide(peptide_strand)))
    # actually, the number should equal to n * (n + 1) / 2 + 1
    # the last 1 is the empty string

    # 1.6 step 6
    # spectrum = [int(mass) for mass in "0 113 128 186 241 299 314 427".split(" ")]
    # masses_strings = set()
    # for peptide in get_cyclopeptide_by_theoretical_spectrum_branch_and_bound(spectrum):
    #     masses_strings.add(transfer_peptide_to_masses_string(peptide))
    # masses_strings = list(masses_strings)
    
    # print_arr(masses_strings)
    # with open("./dataset_100_6.txt", "r") as f:
    #     spectrum = [int(mass) for mass in f.read().strip().split(" ")]
    #     masses_strings = set()
    #     peptides = get_cyclopeptide_by_theoretical_spectrum_branch_and_bound(spectrum)
    #     for peptide in peptides:
    #         masses_strings.add(transfer_peptide_to_masses_string(peptide))
    #     masses_strings = list(masses_strings)

    #     print_arr(masses_strings)

    # Week 3 Quiz
    # problem 2
    # print("###### problem 2 ######")
    # peptide_strand = "PRTEIN"
    # for RNA_strand in "CCUCGUACUGAUAUUAAU  CCCAGGACUGAGAUCAAU  CCCAGUACCGAGAUGAAU  CCUCGUACAGAAAUCAAC".split("  "):
    #     if translate_RNA_to_peptide(RNA_strand) == peptide_strand:
    #         print(RNA_strand)

    # print("\n###### problem 3 ######")
    # peptide_strand = "SYNGE"
    # count = 1
    # for AA in peptide_strand:
    #     codon_list = AA_to_RNA_codon[AA]
    #     count *= len(codon_list)
    # print("num of DNA:", count)

    # print("\n###### problem 4 ######")
    # print(get_weight_by_AA(get_AA_abbr("Gly")))

    # print("\n###### problem 5 ######")
    # peptide_strands = "MLAT  MAIT  IAMT  TAIM  TMLA  TMIA".split("  ")
    # spectrum = [0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416]
    # index_of_option = 1
    # for peptide_strand in peptide_strands:
    #     if get_theoretical_spectrum_of_cyclic_peptide(peptide_strand) == spectrum:
    #         print(index_of_option, peptide_strand)
    #     index_of_option += 1
    
    # print("\n###### problem 6 ######")
    # peptide_strands = "CTV  TCE  CTQ  AQV  ETC  QCV".split("  ")
    # spectrum = [int(mass) for mass in "0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333".split(" ")]
    # index_of_option = 1
    # for peptide_strand in peptide_strands:
    #     is_consistence = True
    #     for mass in get_theoretical_spectrum_of_peptide(peptide_strand):
    #         if not (mass in spectrum):
    #             is_consistence = False
    #             break
    #     if is_consistence:
    #         print(index_of_option, peptide_strand)
    #     index_of_option += 1
