import random
import math
import itertools

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

DNA_BASES = ["A", "T", "C", "G"]

# generate a random DNA strand with given length
def get_random_strand(length):
    result = ""
    for i in range(0, length):
        temp = int(random.random() * 4)
        result += DNA_BASES[temp]
    
    return result

def get_frequency_table(text, k):
    freq_map = {}
    for i in range(0, len(text) - k + 1):
        pattern = text[i : i + k]
        if pattern in freq_map.keys():
            freq_map[pattern] += 1
        else: 
            freq_map[pattern] = 1
    
    return freq_map

# calculate the mean frequency in a freq_map
def get_mean_frequency(freq_map):
    return sum(freq_map.values()) / len(freq_map)

# for a given k-mer substring pattern of text
# increase 1 to the count of every k-mer that has Hamming distance at most d from Pattern
def frequent_words_approximate(text, k, d):
    patterns = []
    freq_map = {}
    for i in range(0, len(text) - k + 1): # 假设text长度为3，要找2-mer，那么[0:1]和[1:2]都要被遍历到，因此i必须能取到1才行
        pattern = text[i : i + k]
        neighborhood = neighbors(pattern, d)
        for j in range(0, len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor not in freq_map.keys():
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] = freq_map[neighbor] + 1
    
    max_count = max_map(freq_map)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_count:
            patterns.append(pattern)
    
    return patterns

# generate all the d-neighbors
# they are defined as a set of all k-mers whose Hamming distance from Pattern does not exceed d
def neighbors(pattern, d):
    # 递归终点
    if d == 0:
        # hamming distance为0，则就是要原来的pattern
        return pattern
    if len(pattern) == 1:
        # pattern长度仅为1，且d不为0，则四种碱基都可以作为1个碱基长的pattern的d-neighbor
        return DNA_BASES
    
    neighborhood = []
    # 取pattern的后缀，即后(k - 1)个
    suffix_pattern = pattern[1:]
    # 找到后缀的所有d-neighbors
    suffix_neighbors = neighbors(suffix_pattern, d)
    for s in suffix_neighbors:
        # 遍历所有后缀的d-neighbors
        if get_hamming_distance(s, suffix_pattern) < d:
            # 如果当前遍历到的后缀的neighbor与后缀的hamming distance < d
            # 则可以将任意碱基放在当前neighbor之前，形成新的neighbor，加入neighbors数组中
            for base in DNA_BASES:
                neighborhood.append(base + s)
        else:
            # hamming distance(s, suffix_pattern) == d
            # 反之，则必须将原来pattern(传入的pattern)的第一个碱基作为新的neighbor的第一个碱基
            # 不然hamming distance就要 > d了
            neighborhood.append(pattern[:1] + s)
    
    return neighborhood

def max_map(d):
    if d.values() == None: 
        return None
    else:
        values = list(d.values())
        max_val = values[0]
        for i in range(0, len(values)):
            if values[i] > max_val:
                max_val = values[i]
        return max_val
    
# The number of mismatches between strings p and q is called the Hamming distance between these strings
def get_hamming_distance(s1, s2):
    result = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            result += 1
    
    return result

# Dna: a collection of strings
# find a (k,d)-motif which appears in every string from Dna with at most d mismatches
def motif_enumeration(Dna, k, d):
    # patterns = set([])
    # # get all k-mer pattern from strings in Dna
    # all_patterns = set([])
    # dna_patterns = []
    # for s in Dna: 
    #     dna_patterns.append(set([]))
    #     for i in range(0, len(s) - k + 1):
    #         # get a k-mer pattern
    #         temp = s[i : i + k]
    #         all_patterns.add(temp)
    #         dna_patterns[len(dna_patterns) - 1].add(temp)
    
    # # traverse each k-mer in Dna
    # for pattern in all_patterns:
    #     # traverse each k-mer pattern' differing from pattern by at most d mismatches
    #     for neighbor in neighbors(pattern, d):
    #         if neighbor != pattern:
    #             pass
    
    # return list(patterns)
    # course里给的pseudo code看不懂啊阿啊阿

    # a solution provided in comments
    # enumerate all the k-mers and their d-neighbors in each string from Dna
    # and get an intersection(交集) of these k-mer sets
    # dna_patterns = []
    # for s in Dna:
    #     dna_patterns.append(set([]))
    #     for i in range(0, len(s) - k + 1):
    #         temp = s[i : i + k]
    #         # generate all neighbors of this k-mer and add them to a set() in dna_patterns
    #         for pattern in neighbors(temp, d):
    #             dna_patterns[len(dna_patterns) - 1].add(pattern)
    # 以上写法感觉很容易角标越界，还是老老实实用2个指针吧，python的语法糖用不习惯
    dna_patterns = [set() for _ in range(len(Dna))]
    for i in range(len(Dna)):
        s = Dna[i]
        for j in range(len(s) - k + 1):
            temp = s[j : j + k]
            for pattern in neighbors(temp, d):
                dna_patterns[i].add(pattern)

        # print(s, "patterns:", dna_patterns[i])
    
    # use "and" logic to all the set() in dna_patterns 取交集！
    result = dna_patterns[0]
    for patterns in dna_patterns:
        result &= patterns

    return list(result) 

# 1.3 step 2
# defining Score(Motifs) as the number of unpopular (lower case) letters in the motif matrix Motifs
def score_motif1(motifs):
    score = 0
    for i in range(len(motifs[0])):
        temp = []
        for motif in motifs:
            temp.append(motif[i])

        score += len(motifs) - max(count_bases(temp).values())
    
    return score

def count_bases(bases_arr):
    counts = {
        "A" : 0, 
        "T" : 0, 
        "C" : 0, 
        "G" : 0
    }
    for base in bases_arr:
        counts[base] += 1
    
    return counts

# 1.3 step 8
# the entropy of a motif matrix is defined as the sum of the entropies of its columns
def score_motif2(motifs):
    profile_matrix = get_profile_matrix(motifs)
    entropy_arr = []
    for i in range(len(motifs[0])):
        for base in profile_matrix.keys():
            temp = 0
            current_P = profile_matrix[base][i]
            if current_P == 0.0:
                temp += 0
            else:
                temp += -(current_P * math.log2(current_P))
            
            entropy_arr.append(temp)
    
    return round(sum(entropy_arr), 4)


# 每个碱基在motif的该位置处出现的次数
def get_count_matrix(motifs):
    counts = {
        "A" : [0 for _ in range(len(motifs[0]))], 
        "T" : [0 for _ in range(len(motifs[0]))], 
        "C" : [0 for _ in range(len(motifs[0]))], 
        "G" : [0 for _ in range(len(motifs[0]))]
    }
    for i in range(len(motifs[0])):
        for motif in motifs:
            current_base = motif[i]
            counts[current_base][i] += 1
    
    return counts

# profile matrix 实际上是每个碱基在motif的该位置处出现的概率
def get_profile_matrix(motifs):
    counts = get_count_matrix(motifs)
    for base in counts.keys():
        for i in range(len(counts[base])):
            counts[base][i] /= len(motifs)
    
    return counts

# consensus motif defined as: 
# the nucleotide base which occurs in each position in the motifs with the highest probability
def get_consensus_motif(motifs):
    profile_matrix = get_profile_matrix(motifs)
    motif = ""
    for i in range(len(motifs[0])):
        current_base = ""
        max_P = 0
        for base in profile_matrix.keys():
            current_P = profile_matrix[base][i]
            if current_P > max_P:
                current_base = base
                max_P = current_P
        motif += current_base
    
    return motif

def get_consensus_motif_by_profile_matrix(profile_matrix, motif_length):
    motif = ""
    for i in range(motif_length):
        current_base = ""
        max_P = 0
        for base in profile_matrix.keys():
            current_P = profile_matrix[base][i]
            if current_P > max_P:
                current_base = base
                max_P = current_P
        motif += current_base
    
    return motif

# given pattern, find a collection of k-mers in Dna that minimizes d(pattern, motifs)
# d(pattern, motifs) represents the total value of the hamming distance from pattern to each motif in motifs
def get_motifs_by_pattern(pattern, Dna):
    k = len(pattern)
    motif_indices = []  # to store the starting indices of all motifs found in each strand
    for strand in Dna:
        l = len(strand)

        if l < k:
            continue

        # assuming that the first k-mer is the motif with minimized hamming distance
        min_d = get_hamming_distance(pattern, strand[: k])
        min_d_index = 0
        for i in range(1, l - k + 1):
            d = get_hamming_distance(pattern, strand[i : i + k])
            if d < min_d:
                min_d = d
                min_d_index = i

        # the index of the k-mer with the minimized hamming distance with pattern is found
        motif_indices.append(min_d_index)
    
    # now that all the indices of the motifs were found
    result = [Dna[i][motif_indices[i] : motif_indices[i] + k] for i in range(len(Dna))]

    return result

# Given a k-mer Pattern and a set of strings Dna = {Dna1, … , Dnat}, 
# we define d(Pattern, Dna) as the sum of distances between Pattern and all strings in Dna,
def get_d_from_Dna_by_pattern(pattern, Dna):
    sum_d = 0
    for motif in get_motifs_by_pattern(pattern, Dna):
        sum_d += get_hamming_distance(motif, pattern)
    
    return sum_d

# a brute force solution to the median string problem
def get_median_string_brute(Dna, k):
    d = k * len(Dna)
    median_string = ""
    for kmer in kmer_enumeration(k):
        current_d = get_d_from_Dna_by_pattern(kmer, Dna)
        if current_d < d:
            d = current_d
            median_string = kmer
    
    return median_string

# Find a Profile-most probable k-mer in a string.
# the parameter profile_matrix can either be a 2-dim array (which was provided in the course material)
# or be a map, the keys are A, T, C, G, and values are arrays
def get_profile_most_probable_kmer(strand, k, profile_matrix):
    matrix = {}
    if isinstance(profile_matrix, dict):
        matrix = profile_matrix
    else:
        # course中matrix的碱基顺序是 A C G T
        matrix = {
            "A" : profile_matrix[0], 
            "C" : profile_matrix[1], 
            "G" : profile_matrix[2], 
            "T" : profile_matrix[3], 
        }
    
    max_prob = -1   # 如果是0的话，有可能根据这个profile_matrix，在传入strand中出现概率确实为0
    most_probable_kmer = "" # 或者把这个变量设为strand[:k]也可以，不然会传回空字符串
    for i in range(len(strand) - k + 1):
        probability = 1
        kmer = strand[i : i + k]
        for j in range(k):
            prob_base = matrix[kmer[j]][j]
            probability *= prob_base

        if probability > max_prob:
            max_prob = probability
            most_probable_kmer = kmer
    
    return most_probable_kmer


# enumerate all the possible combination of 4 DNA bases to generate a k-mer
# this function can also be implemented using the neighbors function when d is set to len(pattern)
def kmer_enumeration(k):
    return [''.join(p) for p in itertools.product(DNA_BASES, repeat=k)]

# GreedyMotifSearch
def greedy_motif_search(Dna, k, t):
    # form a matrix by first k-mers in each strand from Dna
    best_motifs = [strand[:k] for strand in Dna]

    # traverse all the k-mers in the first string from Dna
    for kmer in [Dna[0][i : i + k] for i in range(len(Dna[0]) - k + 1)]:
        motifs = [kmer]
        for strand in Dna[1:]:
            # form a profile_matrix with motif_1, motif_2, ..., motif_i-1
            profile_matrix = get_profile_matrix_with_pseudocounts(motifs)
            # get the most probable motif_i (next_motif) from Dna[i] (current strand)
            next_motif = get_profile_most_probable_kmer(strand, k, profile_matrix)
            # add this most probable motif to motifs
            motifs.append(next_motif)
        
        # now that all the probable motifs based on one k-mer from Dna[0] were found
        # a scoring between motifs obtained in this loop and best_motifs should be done
        # to determine which is better
        if score_motif3(motifs) < score_motif3(best_motifs):
            best_motifs = motifs
    
    return best_motifs

# 通过motifs数组中各个motif得到consensus string后，计算每个motif相对于其的hamming distance     
def score_motif3(motifs):
    consensus_motif = get_consensus_motif(motifs)
    score = 0
    for motif in motifs:
        score += get_hamming_distance(consensus_motif, motif)
    
    return score

# Laplace's succession rule
def get_profile_matrix_with_pseudocounts(motifs):
    counts = get_count_matrix(motifs)
    for base in counts.keys():
        for i in range(len(counts[base])):
            counts[base][i] += 1
            counts[base][i] /= (len(motifs) * 2)
    
    return counts

# 根据传入的motif和profile_matrix，计算这个motif出现的概率
def get_motif_prob_by_profile_matrix(motif, profile_matrix):
    prob = 1
    for i in range(len(motif)):
        current_P = profile_matrix[motif[i]][i]
        prob *= current_P
    
    return prob

if __name__ == "__main__":
    # exer 1.2 step 2
    # 应该是指特定的9-mer，不是计算一个平均值
    # What is the expected number of occurrences of a 9-mer in 500 random DNA strings, 
    # each of length 1000? Assume that the sequences are formed by selecting each nucleotide (A, C, G, T) 
    # with the same probability (0.25).
    # for i in range(0, 10):
        # 1.0019
        # frequency = 0.0
        # for i in range(0, 500):
        #     frequency += get_mean_frequency(get_frequency_table(get_random_strand(1000), 9))
        # print(frequency / 500)
        
        # 2.2398
        # strand = ""
        # for i in range(0, 500):
        #     strand += get_random_strand(1000)
        # frequency = get_mean_frequency(get_frequency_table(strand, 9))
        # print(frequency)

        # 1.0018
        # freq_map = {}
        # for i in range(0, 500):
        #     freq_map.update(get_frequency_table(get_random_strand(1000), 9))
        
        # print(get_mean_frequency(freq_map))
    
    # The question is: we have 500 stretches of 1000 bases; 
    # how many times would we expect to see some specific 9-mer (like AAAATATCT) in there? 
    # Given that all four bases have an equal likelihood to occur (1/4), 
    # any sequence of 9 bases has the same likelihood as any other sequence of 9 bases to occur. 
    # There are 262144 possible different sequences of 9 bases (4 x ... x 4), 
    # therefore each has a likelihood of 1/262144 to occur.
    # From a stretch of 1000 bases we can draw 992 arbitrary sequences of 9 bases. 
    # We are looking for occurrences of a specific 9-mer (like AAAATATCT). 
    # At each draw the likelihood to see it is 1/262144, therefore after 992 draws it is 992 * 1/262144, 
    # and when we repeat this 500 times it is 500 * 992 * 1/262144 = 1.892 expected occurrences.
    # print(500 * (1000 - 9 + 1) / 4**9)

    # 1.2 step 6
    # try to solve the motif finding problem with frequent_words_approximate
    # strand = "ATGACCGGGATACTGATAGAAGAAAGGTTGGGGGCGTACACATTAGATAAACGTATGAAGTACGTTAGACTCGGCGCCGCCGACCCCTATTTTTTGAGCAGATTTAGTGACCTGGAAAAAAAATTTGAGTACAAAACTTTTCCGAATACAATAAAACGGCGGGATGAGTATCCCTGGGATGACTTAAAATAATGGAGTGGTGCTCTCCCGATTTTTGAATATGTAGGATCATTCGCCAGGGTCCGAGCTGAGAATTGGATGCAAAAAAAGGGATTGTCCACGCAATCGCGAACCAACGCGGACCCAAAGGCAAGACCGATAAAGGAGATCCCTTTTGCGGTAATGTGCCGGGAGGCTGGTTACGTAGGGAAGCCCTAACGGACTTAATATAATAAAGGAAGGGCTTATAGGTCAATCATGTTCTTGTGAATGGATTTAACAATAAGGGCTGGGACCGCTTGGCGCACCCAAATTCAGTGTGGGCGAGCGCAACGGTTTTGGCCCTTGTTAGAGGCCCCCGTATAAACAAGGAGGGCCAATTATGAGAGAGCTAATCTATCGCGTGCGTGTTCATAACTTGAGTTAAAAAATAGGGAGCCCTGGGGCACATACAAGAGGAGTCTTCCTTATCAGTTAATGCTGTATGACACTATGTATTGGCCCATTGGCTAAAAGCCCAACTTGACAAATGGAAGATAGAATCCTTGCATACTAAAAAGGAGCGGACCGAAAGGGAAGCTGGTGAGCAACGACAGATTCTTACGTGCATTAGCTCGCTTCCGGGGATCTAATAGCACGAAGCTTACTAAAAAGGAGCGGA"
    # Dna = [
    #     "atgaccgggatactgatAAAAAAAAGGGGGGGggcgtacacattagataaacgtatgaagtacgttagactcggcgccgccg",
    #     "acccctattttttgagcagatttagtgacctggaaaaaaaatttgagtacaaaacttttccgaataAAAAAAAAGGGGGGGa",
    #     "tgagtatccctgggatgacttAAAAAAAAGGGGGGGtgctctcccgatttttgaatatgtaggatcattcgccagggtccga",
    #     "gctgagaattggatgAAAAAAAAGGGGGGGtccacgcaatcgcgaaccaacgcggacccaaaggcaagaccgataaaggaga",
    #     "tcccttttgcggtaatgtgccgggaggctggttacgtagggaagccctaacggacttaatAAAAAAAAGGGGGGGcttatag",
    #     "gtcaatcatgttcttgtgaatggatttAAAAAAAAGGGGGGGgaccgcttggcgcacccaaattcagtgtgggcgagcgcaa",
    #     "cggttttggcccttgttagaggcccccgtAAAAAAAAGGGGGGGcaattatgagagagctaatctatcgcgtgcgtgttcat",
    #     "aacttgagttAAAAAAAAGGGGGGGctggggcacatacaagaggagtcttccttatcagttaatgctgtatgacactatgta",
    #     "ttggcccattggctaaaagcccaacttgacaaatggaagatagaatccttgcatAAAAAAAAGGGGGGGaccgaaagggaag",
    #     "ctggtgagcaacgacagattcttacgtgcattagctcgcttccggggatctaatagcacgaagcttAAAAAAAAGGGGGGGa"
    # ]
    # for i in range(len(Dna)):
    #     Dna[i] = Dna[i].upper()

    # print(Dna)

    # print_arr(motif_enumeration(Dna, 15, 4))
    
    # print(frequent_words_approximate(strand, 15, 4))

    # exer 1.2 step 8
    # with open("dataset_156_8.txt", "r") as f:
    #     f.readline()
    #     # Dna = "ATTTGGC TGCCTTA CGGTATC GAAAATT".split(" ")
    #     # k = 3
    #     # d = 1
    #     Dna = f.readline().strip().split(" ")
    #     k = 5
    #     d = 2
    #     print_arr(motif_enumeration(Dna, k, d))

    # 1.3 step 2
    # results = []
    # for _ in range(999999):
    #     motifs = [get_random_strand(15) for _ in range(10)]
    #     results.append(score_motif1(motifs))
    
    # print(max(results))

    # 1.3 step 8
    # motifs = [
    #     "TCGGGGGTTTTT",
    #     "CCGGTGACTTAC",
    #     "ACGGGGATTTTC",
    #     "TTGGGGACTTTT",
    #     "AAGGGGACTTCC",
    #     "TTGGGGACTTCC",
    #     "TCGGGGATTCAT",
    #     "TCGGGGATTCCT",
    #     "TAGGGGAACTAC",
    #     "TCGGGTATAACC"
    # ]
    # print(get_count_matrix(motifs))
    # print(get_profile_matrix(motifs))
    # print("entropy:", score_motif2(motifs))
    # print("consensus motif:", get_consensus_motif(motifs))

    # 1.4 step 5
    # Dna = [
    #     "TTACCTTAAC", 
    #     "GATATCTGTC", 
    #     "ACGGCGTTCG", 
    #     "CCCTAAAGAG", 
    #     "CGTCAGAGGT", 
    # ]
    # print_arr(get_motifs_by_pattern("AAA", Dna))

    # 1.4 step 8
    # print_arr(kmer_enumeration(3))
    # Dna = [
    #     "TTACCTTAAC", 
    #     "GATATCTGTC", 
    #     "ACGGCGTTCG", 
    #     "CCCTAAAGAG", 
    #     "CGTCAGAGGT", 
    # ]
    # print(get_median_string_brute(Dna, 3))
    # print(get_d_from_Dna_by_pattern("CCT", Dna))

    # 1.7 step 1
    # with open("dataset_5164_1.txt", "r") as f:
    #     pattern = f.readline().strip()
    #     Dna = f.readline().strip().split(" ")
        
    #     print(get_d_from_Dna_by_pattern(pattern, Dna))

    # 1.4 step 9
    # with open("dataset_158_9.txt", "r") as f:
    #     k = int(f.readline().strip())

    #     Dna = []
    #     while True:
    #         strand = f.readline().strip()
    #         if strand:
    #             Dna.append(strand)
    #         else:
    #             break
        
    #     print(get_median_string_brute(Dna, k))
    #     print(get_d_from_Dna_by_pattern("ATGCCC", Dna))

    # 1.5 step 3
    # strand = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
    # k = 5
    # profile_matrix = [
    #     [float(i) for i in "0.2 0.2 0.3 0.2 0.3".split(" ")], 
    #     [float(i) for i in "0.4 0.3 0.1 0.5 0.1".split(" ")], 
    #     [float(i) for i in "0.3 0.3 0.5 0.2 0.4".split(" ")], 
    #     [float(i) for i in "0.1 0.2 0.1 0.1 0.2".split(" ")], 
    # ]
    # print(get_profile_most_probable_kmer(strand, k, profile_matrix))
    # with open("dataset_159_3.txt", "r") as f:
    #     strand = f.readline().strip()
    #     k = int(f.readline().strip())
    #     profile_matrix = [
    #         [float(i) for i in f.readline().strip().split(" ")], 
    #         [float(i) for i in f.readline().strip().split(" ")], 
    #         [float(i) for i in f.readline().strip().split(" ")], 
    #         [float(i) for i in f.readline().strip().split(" ")], 
    #     ]
    #     print(get_profile_most_probable_kmer(strand, k, profile_matrix))

    # 1.5 step 5
    # Dna = "GGCGTTCAGGCA AAGAATCAGTCA CAAGGAGTTCGC CACGTCAATCAC CAATAATATTCG".split(" ")
    # k = 3
    # t = 5
    # print(greedy_motif_search(Dna, k, t))
    # with open("dataset_159_5.txt", "r") as f:
    #     params = f.readline().strip().split(" ")
    #     k = int(params[0])
    #     t = int(params[1])

    #     Dna = f.readline().strip().split(" ")
    #     # print(Dna)
        
    #     print_arr(greedy_motif_search(Dna, k, t))

    # 1.6 step 9
    # with open("dataset_160_9.txt", "r") as f:
    #     params = f.readline().strip().split(" ")
    #     k = int(params[0])
    #     t = int(params[1])

    #     Dna = f.readline().strip().split(" ")
    #     # print(Dna)
        
    #     print_arr(greedy_motif_search(Dna, k, t))

    # deal with the final question to find AAAAAAAATTTTTTT in 10 strands
    # with open("subtle_motif_dataset.txt", "r") as f:
    #     k = 15
    #     Dna = []
    #     while True:
    #         strand = f.readline().strip()
    #         if strand:
    #             strand = strand.replace("*", "")
    #             Dna.append(strand)
    #         else:
    #             break
        
        # too slow
        # print(get_median_string_brute(Dna, k))

        # not accurate enough
        # However, it trades speed for accuracy and returns gtAAAtAgaGatGtG (total distance: 58), 
        # which is very different from the true implanted motif AAAAAAAAGGGGGGG.
        # after pseudocount is implemented in generating profile matrix, the accuracy improved remarkably
        # print("result:")
        # print_arr(greedy_motif_search(Dna, k, 10))
        # print("score of motifs found:", score_motif3(greedy_motif_search(Dna, k, 10)))
        # print("consensus strand found:", get_consensus_motif(greedy_motif_search(Dna, k, 10)))

        # # 40
        # print("score of motifs implanted:", score_motif3([
        #     "ACGAGAAAGGGAGGG",
        #     "AAAAAATAGCAGGGT",
        #     "AAATAAACAGCGGGG",
        #     "GAAAAAAAGGGGTTT",
        #     "ATAAAGTAGAGGGGG",
        #     "CTAAAAATGGGGCGG",
        #     "AAAAAGAGAAGGGGG",
        #     "ATAGAAAAGGAAGGG",
        #     "AAAAAAGAGAGGAGT",
        #     "AAGCTAAAGGGGGGT",
        # ]))
        # print("consensus strand implanted:", get_consensus_motif([
        #     "ACGAGAAAGGGAGGG",
        #     "AAAAAATAGCAGGGT",
        #     "AAATAAACAGCGGGG",
        #     "GAAAAAAAGGGGTTT",
        #     "ATAAAGTAGAGGGGG",
        #     "CTAAAAATGGGGCGG",
        #     "AAAAAGAGAAGGGGG",
        #     "ATAGAAAAGGAAGGG",
        #     "AAAAAAGAGAGGAGT",
        #     "AAGCTAAAGGGGGGT",
        # ]))

    # Week 3 Quiz
    # temp += -(current_P * math.log2(current_P))
    # problem 3
    print("problem 3:")
    prob_A = [0.5, 0, 0, 0.5]
    prob_B = [0.25, 0.25, 0.25, 0.25]
    prob_C = [0, 0, 0, 1]
    prob_D = [0.25, 0, 0.5, 0.25]
    entropy_A = 0
    entropy_B = 0
    entropy_C = 0
    entropy_D = 0
    for p in prob_A:
        if p == 0:
            entropy_A += 0
        else:
            entropy_A += -(p * math.log2(p))
    for p in prob_B:
        if p == 0:
            entropy_B += 0
        else:
            entropy_B += -(p * math.log2(p))
    for p in prob_C:
        if p == 0:
            entropy_C += 0
        else:
            entropy_C += -(p * math.log2(p))
    for p in prob_D:
        if p == 0:
            entropy_D += 0
        else:
            entropy_D += -(p * math.log2(p))
    print("entropy A:", entropy_A, 
          "\nentropy B:", entropy_B, 
          "\nentropy C:", entropy_C, 
          "\nentropy D:", entropy_D)
    
    # problem 4
    print("\nproblem 4:")
    profile_matrix = {
        "A" : [float(p) for p in "0.4  0.3  0.0  0.1  0.0  0.9".split("  ")], 
        "C" : [float(p) for p in "0.2  0.3  0.0  0.4  0.0  0.1".split("  ")], 
        "G" : [float(p) for p in "0.1  0.3  1.0  0.1  0.5  0.0".split("  ")],
        "T" : [float(p) for p in "0.3  0.1  0.0  0.4  0.5  0.0".split("  ")]
    }
    # print(get_consensus_motif_by_profile_matrix(profile_matrix, 6))
    options = [
        "ACGTTA", 
        "AAGAGA", 
        "TCGCGA", 
        "ACGCGA", 
        "AGGTCA", 
        "AAGCTA"
    ]
    for option in options:
        print(option, get_motif_prob_by_profile_matrix(option, profile_matrix))

    # problem 5
    print("\nproblem 5:")
    Dna = [
        "CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC", 
        "GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC", 
        "GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG", 
    ]
    answer = get_median_string_brute(Dna, 7)
    options = [
        "GGTTACT", 
        "AATCCTA", 
        "ATAACGG", 
        "GAACCAC", 
        "CGTGTAA", 
        "AACGCTG"
    ]
    for option in options:
        print(option, get_d_from_Dna_by_pattern(option, Dna))
    print("\nThe answer by my algorithm:", answer, get_d_from_Dna_by_pattern(answer, Dna), )

    # problem 6
    print("problem 6:")
    motif = "TCGGTA"
    profile_matrix = {
        "A" : [float(p) for p in "0.4  0.3  0.0  0.1  0.0  0.9".split("  ")], 
        "C" : [float(p) for p in "0.2  0.3  0.0  0.4  0.0  0.1".split("  ")], 
        "G" : [float(p) for p in "0.1  0.3  1.0  0.1  0.5  0.0".split("  ")],
        "T" : [float(p) for p in "0.3  0.1  0.0  0.4  0.5  0.0".split("  ")]
    }
    profile_matrix_with_pseudocounts = {}
    prob = 1
    for i in range(len(motif)):
        current_P = profile_matrix[motif[i]][i]
        prob *= current_P
    print(prob)
    