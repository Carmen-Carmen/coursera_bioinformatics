import random
import math

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

DNA_BASES = ["A", "T", "C", "G"]

# The number of mismatches between strings p and q is called the Hamming distance between these strings
def get_hamming_distance(s1, s2):
    result = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            result += 1
    
    return result

# Find a Profile-most probable k-mer in a string.
# the parameter profile_matrix can either be a 2-dim array (which was provided in the course material)
# or be a map, the keys are A, T, C, G, and values are arrays
# it is better if the profile_matrix is with pseudocount
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

# For each k-mer Pattern in Text, compute the probability Pr(Pattern | Profile), 
# resulting in n = |Text| - k + 1 probabilities (p1, …, pn).
# GibbsSampler uses this random number generator to select a Profile-randomly generated k-mer at each step
def get_profile_randomly_generated_kmer(strand, k, profile_matrix):
    kmer_arr = [strand[i : i + k] for i in range(len(strand) - k + 1)]
    p_arr = []
    for kmer in kmer_arr:
        kmer_P = 1
        for i in range(k):
            # the possibility of this kmer equals the multiplying of 
            # all the possibilities of all DNA bases appering in the i-place
            current_base = kmer[i]
            kmer_P *= profile_matrix[current_base][i]
        
        kmer_P = int(kmer_P * math.pow(10, k))
        p_arr.append(kmer_P)
    
    return random.choices(
        kmer_arr, 
        weights=p_arr, 
        k=1
    )[0]

def get_profile_matrix_with_pseudocounts(motifs):
    counts = get_count_matrix(motifs)
    for base in counts.keys():
        for i in range(len(counts[base])):
            counts[base][i] += 1
            counts[base][i] /= (len(motifs) * 2)
    
    return counts

# 通过motifs数组中各个motif得到consensus string后，计算每个motif相对于其的hamming distance     
def score_motif3(motifs):
    consensus_motif = get_consensus_motif(motifs)
    score = 0
    for motif in motifs:
        score += get_hamming_distance(consensus_motif, motif)
    
    return score

# consensus motif defined as: 
# the nucleotide base which occurs in each position in the motifs with the highest probability
def get_consensus_motif(motifs):
    profile_matrix = get_profile_matrix_with_pseudocounts(motifs)
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

# given a collection of strings Dna and an arbitrary 4 x k matrix Profile, 
# we define Motifs(Profile, Dna) as the collection of k-mers formed by the 
# Profile-most probable k-mers in each string from Dna
# get motifs based on given profile matrix from a Dna list
def get_motifs_by_profile(profile_matrix, Dna):
    motifs = []
    
    k = 0
    if isinstance(profile_matrix, dict):
        for arr in profile_matrix.values():
            k = len(arr)
            break
    else:
        k = len(profile_matrix[0])

    for strand in Dna:
        motifs.append(get_profile_most_probable_kmer(strand, k, profile_matrix))
    return motifs

# bioinformaticians usually run this algorithm thousands of times. 
# On each run, they begin from a new randomly selected set of k-mers, 
# selecting the best set of k-mers found in all these runs.
def randomized_motif_search(Dna, k, t):
    motifs_arr = []
    for _ in range(1000):
        motifs = []
        for strand in Dna:
            index = int(random.random() * (len(strand) - k + 1))
            motifs.append(strand[index : index + k])
        
        best_motifs = motifs

        while True:
            # 每次都根据上一轮产生的motifs，产生新的profile
            profile_matrix = get_profile_matrix_with_pseudocounts(motifs)

            # 根据新的profile，从Dna的每条strand中找到最可能的motif，返回数组
            motifs = get_motifs_by_profile(profile_matrix, Dna)

            if score_motif3(motifs) < score_motif3(best_motifs):
                best_motifs = motifs
            else:
                motifs_arr.append(best_motifs)
                break
    
    return get_best_motifs(motifs_arr)

# GibbsSampler further generalizes the random number generator by 
# using the function Random(p1, …, pn) defined for any set of non-negative numbers
def random_with_P(p_arr):
    normalized_p_arr = [int(p * 1000) for p in p_arr]
    nums = [i for i in range(len(p_arr))]
    return random.choices(nums, weights=normalized_p_arr, k=1)[0]

# repeat for N times
def gibbs_sampler_motif_search(Dna, k, t, N):
    motifs_arr = []
    for __ in range(2000):
        # print("round_:", __)
        motifs = []
        for strand in Dna:
            index = int(random.random() * (len(strand) - k + 1))
            motifs.append(strand[index : index + k])
        
        best_motifs = motifs

        for _ in range(N):
            # print("round N:", _)
            i = int(random.random() * t)
            motif_i = motifs[i]
            motifs.remove(motif_i)

            profile_matrix = get_profile_matrix_with_pseudocounts(motifs)
            motif_i = get_profile_randomly_generated_kmer(Dna[i], k, profile_matrix)
            # be cautious that the index of motif_i in motifs is correlated to the index of strand which 
            # it comes from in the Dna array!!
            motifs.insert(i, motif_i)

            if score_motif3(motifs) < score_motif3(best_motifs):
                best_motifs = motifs
        
        motifs_arr.append(best_motifs)
    
    return get_best_motifs(motifs_arr)
    
def get_best_motifs(motifs_arr):
    low_score = score_motif3(motifs_arr[0])
    best_motifs = None
    for motifs in motifs_arr:
        current_score = score_motif3(motifs)
        if current_score < low_score:
            low_score = current_score
            best_motifs = motifs
    
    return best_motifs
    
    
if __name__ == "__main__":
    # 1.1 step 5
    # Dna = "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA".split(" ")
    # k = 8
    # t = 5
    # motifs = randomized_motif_search(Dna, k, t)
    # print_arr(motifs)
    # print(score_motif3(motifs))
    # print(score_motif3("TCTCGGGG CCAAGGTG TACAGGCG TTCAGGTG TCCACGTG".split(" ")))
    # with open("./dataset_161_5.txt", "r") as f:
    #     params = f.readline().strip()
    #     k = int(params.split(" ")[0])
    #     t = int(params.split(" ")[1])
    #     Dna = f.readline().strip().split(" ")

    #     motifs = randomized_motif_search(Dna, k, t)
    #     print_arr(motifs)
    #     print(score_motif3(motifs))

    # the implanted (15, 4) problem
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
        
    #     print("result:")
    #     # motifs = randomized_motif_search(Dna, k, 10)
    #     motifs = gibbs_sampler_motif_search(Dna, k ,10, 200)
    #     print_arr(motifs)
    #     print("score of motifs found:", score_motif3(motifs))
    #     print("consensus strand found:", get_consensus_motif(motifs))

    #     # 40
    #     print("score of motifs implanted:", score_motif3([
    #         "ACGAGAAAGGGAGGG",
    #         "AAAAAATAGCAGGGT",
    #         "AAATAAACAGCGGGG",
    #         "GAAAAAAAGGGGTTT",
    #         "ATAAAGTAGAGGGGG",
    #         "CTAAAAATGGGGCGG",
    #         "AAAAAGAGAAGGGGG",
    #         "ATAGAAAAGGAAGGG",
    #         "AAAAAAGAGAGGAGT",
    #         "AAGCTAAAGGGGGGT",
    #     ]))
    #     print("consensus strand implanted:", get_consensus_motif([
    #         "ACGAGAAAGGGAGGG",
    #         "AAAAAATAGCAGGGT",
    #         "AAATAAACAGCGGGG",
    #         "GAAAAAAAGGGGTTT",
    #         "ATAAAGTAGAGGGGG",
    #         "CTAAAAATGGGGCGG",
    #         "AAAAAGAGAAGGGGG",
    #         "ATAGAAAAGGAAGGG",
    #         "AAAAAAGAGAGGAGT",
    #         "AAGCTAAAGGGGGGT",
    #     ]))

    # 1.2 step 3
    # Compute the probability that ten randomly selected 15-mers from 
    # the ten 600-nucleotide long strings in the Subtle Motif Problem capture at least one implanted 15-mer.
    # possibility of not selecting any 15-mer in 10 600-nucleotide long strands: 
    # p1 = ((600 - 15) / (600 - 15 + 1)) ^ 10
    # so the possibility of selecting at least 1 15-mer: 1 - p1 = 0.016934

    # 1.2 step 4
    # Compute the probability that ten randomly selected 15-mers from 
    # ten 600-nucleotide long strings (as in the Subtle Motif Problem) capture at least two implanted 15-mers
    # the possibility of selecting exactly 1 15-mer in 10 600-nucleotide long strands: 
    # p2 = (1 / (600 - 15 + 1) * 10) * ((600 - 15) / (600 - 15 + 1)) ^ 9
    # 15-mer selected in any strand     *     no 15-mer selected in other strands
    # in step 3, the possibility of selecting at least 1 15-mer is calculated
    # so the possibility of selecting at least 2 15-mers could be calculated by 
    # minusing the possibility of selecting exactly 1 15-mer: 0.16934 - p2 = 0.000129

    # 1.3 step 2
    # for _ in range(100):
    #     print(random_with_P([0.1, 0.2, 0.9]))

    # 1.3 step 4
    # Dna = "CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG TAGTACCGAGACCGAAAGAAGTATACAGGCGT TAGATCAAGTTTCAGGTGCACGTCGGTGAACC AATCCACCAGCTCCACGTGCAATGTTGGCCTA".split(" ")
    # k = 8
    # t = 5
    # N = 100
    # print_arr(gibbs_sampler_motif_search(Dna, k, t, N))
    # with open("./dataset_163_4.txt", "r") as f:
    #     params = f.readline().strip()
    #     k = int(params.split(" ")[0])
    #     t = int(params.split(" ")[1])
    #     N = int(params.split(" ")[2])
    #     Dna = f.readline().strip().split(" ")

    #     motifs = gibbs_sampler_motif_search(Dna, k, t, N)
    #     print_arr(motifs)
    #     print(score_motif3(motifs))

    # the dormancy survival regulator DosR motif finding problem
    # with open("./DosR.txt", "r") as f:
    #     Dna = []
    #     while True:
    #         strand = f.readline().strip()
    #         if strand:
    #             Dna.append(strand)
    #         else:
    #             break
        
    #     t = 10
    #     N = 100
    #     for k in range(8, 13):
    #         print("********************")
    #         print("k =", k)
    #         randomized_motifs = randomized_motif_search(Dna, k, t)

    #         print("randomized result:")
    #         print_arr(randomized_motifs)
    #         print("consensus motif:", get_consensus_motif(randomized_motifs))
    #         print("score:", score_motif3(randomized_motifs))

    #         gibbs_motifs = gibbs_sampler_motif_search(Dna, k, t, N)
            
    #         print("gibbs sampling result:")
    #         print_arr(gibbs_motifs)
    #         print("consensus motif:", get_consensus_motif(gibbs_motifs))
    #         print("score:", score_motif3(gibbs_motifs))

    # Week 4 Quiz
    # problem 5
    Dna = [
        "AAGCCAAA", 
        "AATCCTGG",
        "GCTACTTG",
        "ATGTTTTG"
    ]

    motifs = [
        "CCA", 
        "CCT", 
        "CTT", 
        "TTG"
    ]

    profile_matrix = get_profile_matrix_with_pseudocounts(motifs)
    motifs = get_motifs_by_profile(profile_matrix, Dna)
    print_arr(motifs)
