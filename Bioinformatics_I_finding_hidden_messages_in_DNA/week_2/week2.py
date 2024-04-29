import matplotlib.pyplot as plt

def get_skew_arr(strand):
    skew_arr = [0]
    for i in range(0, len(strand)):
        temp = skew_arr[i]
        if strand[i] == "C":
            skew_arr.append(temp - 1)
        elif strand[i] == "G":
            skew_arr.append(temp + 1)
        else:
            skew_arr.append(temp)

    
    return skew_arr

def print_arr(arr):
    result = ""
    for item in arr:
        result += str(item) + " "
    
    print(result.strip())

# return the index where the skew achieve a minimum
# that is where the reverse half-strand ends and the forward half-strand begins
def get_min_skew_indices(skew_arr):
    min_skew = min(skew_arr)
    results = []
    for i in range(0, len(skew_arr)):
        if skew_arr[i] == min_skew:
            results.append(i)
    
    return results

# The number of mismatches between strings p and q is called the Hamming distance between these strings
def get_hamming_distance(s1, s2):
    result = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            result += 1
    
    return result

# from text, find the count which pattern appears
# with smaller or equal hamming distance than tolerance
def get_pattern_indices_approximate(text, pattern, tolerance):
    indices = []
    for i in range(0, len(text) - len(pattern) + 1):
        temp = text[i : i + len(pattern)]
        if get_hamming_distance(temp, pattern) <= tolerance:
            indices.append(i)
    
    return indices

def pattern_count_approximate(text, pattern, tolerance):
    count = 0
    for i in range(0, len(text) - len(pattern) + 1):
        temp = text[i:i + len(pattern)]
        if get_hamming_distance(temp, pattern) <= tolerance:
            count += 1
    return count

DNA_BASES = ["A", "T", "C", "G"]
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

# return the d-neighbors with exact d mismatches from pattern
def neighbors_exact(pattern, d):
    # 递归终点
    if d == 0:
        # hamming distance为0，则就是要原来的pattern
        return set([pattern])
    if len(pattern) == 1:
        # pattern长度仅为1，且d不为0，则四种碱基都可以作为1个碱基长的pattern的d-neighbor
        return set(DNA_BASES)
    
    neighborhood = []
    # 取pattern的后缀，即后(k - 1)个
    suffix_pattern = pattern[1:]
    # 找到后缀的所有d-neighbors
    suffix_neighbors = neighbors_exact(suffix_pattern, d) | neighbors_exact(suffix_pattern, d - 1)
    for s in suffix_neighbors:
        # 遍历所有后缀的d-neighbors
        if len([a for a, b in zip(suffix_pattern, s) if a != b]) == d - 1:
            for char in set(DNA_BASES) - set(pattern[0]):
                neighborhood.append(char + s)
        elif len([a for a, b in zip(suffix_pattern, s) if a != b]) == d:
            neighborhood.append(pattern[0] + s)
    
    return set(neighborhood)

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

# Find the most frequent k-mers (with mismatches and reverse complements) in a string.
def frequent_words_rc_approximate(text, k, d):
    patterns = []
    freq_map = {}
    for i in range(0, len(text) - k + 1): # 假设text长度为3，要找2-mer，那么[0:1]和[1:2]都要被遍历到，因此i必须能取到1才行
        pattern = text[i : i + k]
        pattern_rc = get_complement_strand(pattern)[::-1]
        neighborhood = neighbors(pattern, d) + neighbors(pattern_rc, d)
        for j in range(0, len(neighborhood)):
            neighbor = neighborhood[j]
            if neighbor not in freq_map.keys():
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] = freq_map[neighbor] + 1
    
    max_count = max_map(freq_map)
    print("max count:", max_count)
    for pattern in freq_map.keys():
        if freq_map[pattern] == max_count:
            patterns.append(pattern)
    
    return patterns

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
    
def get_complement_strand(s):
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
    
    return result
    
if __name__ == "__main__":
    # exer 1.3 step 8
    # print_arr(get_skew_arr("GAGCCACCGCGATA"))

    # exer 1.3 step 10
    # with open("E_coli.txt", "r") as f:
    #     # strand = "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
    #     strand = f.read()
    #     skew_arr = get_skew_arr(strand)
    #     print_arr(get_min_skew_indices(skew_arr))
    #     # plotting
    #     # generate the x axis
    #     x = [i for i in range(0, len(strand))]
    #     x.append(len(strand))
    #     y = skew_arr
    #     plt.plot(x, y)
    #     plt.show()

    # exer 1.4 step 3
    # with open("dataset_9_3.txt", "r") as f:
    #     s1 = f.readline()
    #     s2 = f.readline()
    #     print(get_hamming_distance(s1, s2))

    # exer 1.4 step 4
    # with open("dataset_9_4.txt", "r") as f:
    #     # text = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
    #     # pattern = "ATTCTGGA"
    #     # tolerance = 3
    #     # readline()函数会把\n也读进来，如果不strip的话实际上进行比较的是一个带了\n的(k + 1)-mer pattern，一定会有一个hamming distance
    #     pattern = f.readline().strip()
    #     text = f.readline().strip()
    #     tolerance = int(f.readline().strip())
    #     print_arr(get_pattern_indices_approximate(text, pattern, tolerance))

    # exer 1.4 step 5
    # with open("dataset_9_6.txt", "r") as f:
    #     # text = "AACAAGCTGATAAACATTTAAAGAG"
    #     # pattern = "AAAAA"
    #     # tolerance = 2
    #     pattern = f.readline().strip()
    #     text = f.readline().strip()
    #     tolerance = int(f.readline().strip())
    #     print(pattern_count_approximate(text, pattern, tolerance))

    # exer 1.7 step 4
    # with open("dataset_3014_4.txt", "r") as f:
    #     # pattern = "ACG"
    #     # d = 1
    #     pattern = f.readline().strip()
    #     d = int(f.readline().strip())
    #     # print_arr(neighbors(pattern, d))
    #     # print_arr(neighbors_exact(pattern, d))
    #     print(len(neighbors(pattern, d)))
    #     print(len(neighbors_exact(pattern, d)))

    # exer 1.4 step 9
    # with open("dataset_9_9.txt", "r") as f:
    #     text = f.readline().strip()
    #     k = 7
    #     d = 3
    #     print(frequent_words_approximate(text, k, d))

    # exer 1.4 step 10
    # with open("dataset_9_10.txt", "r") as f:
    #     # text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    #     text = f.readline().strip()
    #     k = 6
    #     d = 3
    #     print_arr(frequent_words_rc_approximate(text, k, d))

    # find dnaA box for E coli
    # with open("E_coli.txt", "r") as f:
    #     text = f.read()
    #     skew_arr = get_skew_arr(text)
    #     ori_index = get_min_skew_indices(skew_arr)[0]
    #     print("ori_index:", ori_index)
    #     text = text[ori_index : ori_index + 500]
    #     k = 9
    #     d = 1
    #     print_arr(frequent_words_rc_approximate(text, k, d))

    # final challange find a DnaA box in Salmonella enterica
    # with open("Salmonella_enterica.txt", "r") as f:
    #     f.readline()    # 第一行不是基因组
    #     # 在读基因组进来的时候要把每一行的\n去掉
    #     text = ""
    #     line = ""
    #     while True:
    #         line = f.readline().strip()
    #         if line: 
    #             text += line
    #         else:
    #             break
    #     skew_arr = get_skew_arr(text)

    #     # plot
    #     x = [i for i in range(0, len(text))]
    #     x.append(len(text))
    #     y = skew_arr
    #     plt.plot(x, y)
    #     plt.show()

    #     ori_index = get_min_skew_indices(skew_arr)[0]
    #     print("ori_index:", ori_index)
    #     text = text[ori_index : ori_index + 500]
    #     k = 9
    #     d = 1
    #     print_arr(frequent_words_rc_approximate(text, k, d))

    # Quiz
    # q2
    s1 = "TGACCCGTTATGCTCGAGTTCGGTCAGAGCGTCATTGCGAGTAGTCGTTTGCTTTCTCAAACTCC"
    s2 = "GAGCGATTAAGCGTGACAGCCCCAGGGAACCCACAAAACGTGATCGCAGTCCATCCGATCATACA"
    print("q2 Hamming distance:", get_hamming_distance(s1, s2))

    # q3
    print_arr(get_min_skew_indices(get_skew_arr("CATTCCAGTACTTCGATGATGGCGTGAAGA")))

    # q4
    print("count1('TACGCATTACAAAGCACA', 'AA')", pattern_count_approximate("TACGCATTACAAAGCACA", "AA", 1))

    # q5
    print("number of 1-neighbor of 'CCAGTCAATG':", len(neighbors("CCAGTCAATG", 1)))


