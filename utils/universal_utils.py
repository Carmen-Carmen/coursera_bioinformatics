import pyperclip

def print_arr(arr):
    to_print = ""
    for item in arr:
        to_print += str(item) + " "
    
    to_print = to_print.strip()
    print(to_print)
    return to_print

def parse_arr_2_dimension(arr_2_dim):
    parsed = []
    rows = arr_2_dim.split("\n")
    for row in rows:
        row = row.strip()
        splitting_char = " "
        if "\t" in row:
            splitting_char = "\t"
        elif " " in row: 
            splitting_char = " "
        
        temp = [int(item) for item in row.split(splitting_char)]
        parsed.append(temp)
    
    return parsed

def parse_arr(arr):
    if "\t" in arr: 
        return arr.split("\t")
    elif "\n" in arr: 
        return arr.split("\n")
    elif "; " in arr: 
        return arr.split("; ")
    elif " " in arr: 
        return arr.split(" ")
    
def get_key_by_distincet_val(dict, val): 
    for key in dict.keys(): 
        if dict[key] == val: 
            return key
    
    return None

def get_hamming_distance(s1, s2):
    result = 0
    for i in range(0, len(s1)):
        if s1[i] != s2[i]:
            result += 1
    
    return result