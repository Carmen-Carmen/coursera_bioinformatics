import random
import time
import pyperclip

# choose the reversal sites randomly
def random_reversal_sorting(permutation, identity):
    i, j = 0, 0
    rearrange_sites = []
    step = 0
    while True:
        i = random.choice(range(len(permutation)))
        while j <= i:
            j = random.choice(range(len(permutation) + 1))
        
        cutted = permutation[i : j]
        cutted = [item for item in reversed(cutted)]
        for index in range(len(cutted)):
            cutted[index] = -cutted[index]
        
        prefix = permutation[: i]
        suffix = permutation[j: ]
        permutation = prefix + cutted + suffix
        # time.sleep(0.1)
        # print(permutation)
        if permutation == identity:
            print("steps count: %d" % step)
            return step, rearrange_sites
        
        rearrange_sites.append([i, j])
        step += 1

def fewest_step_reversal_sorting(permutation, identity, count):
    min_rearrange_sites = []
    min_step = 100
    for _ in range(count):
        step, rearrange_sites = random_reversal_sorting(permutation, identity)
        if step < min_step:
            min_step = step
            min_rearrange_sites = rearrange_sites
    
    return min_rearrange_sites

# Input: A permutation P.
# Output: The sequence of permutations corresponding to 
# applying GreedySorting to P, ending with the identity permutation.
def greedy_reversal_sorting(permutation):
    # print_permutation(permutation)
    greedy_reversal_distance = 0
    len_p = len(permutation)
    output = ""
    for i in range(1, len_p + 1):
        element = permutation[i - 1]
        element_abs = abs(element)
        # i-th element not in correct position, do the reversal from i to index(i)
        if element_abs != i:
            # find the index of element with abs == i
            j = 0
            try:
                j = permutation.index(i)
            except ValueError:
                j = permutation.index(-i)

            cutted = permutation[i - 1 : j + 1]
            cutted = [item for item in reversed(cutted)]
            for index in range(len(cutted)):
                cutted[index] = -cutted[index]
            
            # glue permutation
            prefix = permutation[: i - 1]
            suffix = permutation[j + 1: ]
            permutation = prefix + cutted + suffix

            greedy_reversal_distance += 1
            output += print_permutation(permutation) + "\n"

        # i-th element in correct position, but not in correct direction
        if permutation[i - 1] != i:
            permutation[i - 1] = - permutation[i - 1]
            greedy_reversal_distance += 1
            output += print_permutation(permutation) + "\n"

    output = output.strip()
    pyperclip.copy(output)

    return greedy_reversal_distance

def print_permutation(permutation):
    temp = ""
    for item in permutation:
        if item > 0:
            temp += "+%d " %item
        else:
            temp += "%d " %item
    temp = temp.strip()
    print(temp)
    return temp

# Number of Breakpoints Problem: Find the number of breakpoints in a permutation.
# Input: A permutation.
# Output: The number of breakpoints in this permutation.
def count_breakpoints(permutation):
    permutation = [0] + permutation
    permutation.append(len(permutation))

    # print_permutation(permutation)
    breakpoints_count = 0
    if permutation[1] != 1:
        breakpoints_count += 1
    for i in range(1, len(permutation) - 1):
        prev = permutation[i]
        next = permutation[i + 1]
        if next != prev + 1:
            breakpoints_count += 1
    
    return breakpoints_count

if __name__ == "__main__":
    # random reverse sorting is impractical
    # permutation = [int(item) for item in "+2 -4 -3 +5 -8 -7 -6 +1".split(" ")]
    # identity = [item for item in range(1, len(permutation) + 1)]
    # for site in fewest_step_reversal_sorting(
    #     permutation, 
    #     identity, 
    #     1000
    # ):
    #     print(site)

    # 1.4 step 4
    # permutation = [int(item) for item in "-3 +4 +1 +5 -2".split(" ")]
    # print("steps count: %d" %greedy_reversal_sorting(permutation))
    # with open("./dataset_286_4.txt", "r") as f:
    #     permutation = [int(item) for item in f.read().strip().split(" ")]
    #     print("steps count: %d" %greedy_reversal_sorting(permutation))

    # 1.4 step 6
    # permutation = [int(item) for item in "+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14".split(" ")]
    # print("breakpoints count: %d" %count_breakpoints(permutation))
    # with open("./dataset_287_6.txt", "r") as f:
    #     permutation = [int(item) for item in f.read().strip().split(" ")]
    #     print("breakpoints count: %d" %count_breakpoints(permutation))

    # week 4 Quiz
    # problem 2
    print("\n### problem 2 ###")
    permutation = [int(item) for item in "+20 +7 +10 +9 +11 +13 +18 -8 -6 -14 +2 -4 -16 +15 +1 +17 +12 -5 +3 -19".split(" ")]
    print("steps count: %d" %greedy_reversal_sorting(permutation))

    print("\n### problem 3 ###")
    permutation = [int(item) for item in "+20 +8 +9 +10 +11 +12 +18 -7 -6 -14 +2 -17 -16 -15 +1 +4 +13 -5 +3 -19".split(" ")]
    print("breakpoints count: %d" %count_breakpoints(permutation))
    




        
