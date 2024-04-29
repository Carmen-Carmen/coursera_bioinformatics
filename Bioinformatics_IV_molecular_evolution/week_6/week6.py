from coursera_bioinformatics.utils import DNA_RNA_utils
from coursera_bioinformatics.utils import universal_utils
import pyperclip
import re
import copy
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_1.week1 import *
from coursera_bioinformatics.Bioinformatics_IV_molecular_evolution.week_2.week2 import *

if __name__ == "__main__": 
    distance_matrix_str = """0 11 2 16
11 0 13 15
2 13 0 9
16 15 9 0"""
    mapping = {
        '0': 'W',
        '1': 'X',
        '2': 'Y',
        '3': 'Z',
        '4': 'A',
        '5': 'B' 
    }
    distance_matrix = universal_utils.parse_arr_2_dimension(distance_matrix_str)
    print(distance_matrix)
    neighbor_joining_matrix = get_neighbor_joining_matrix(distance_matrix)
    # for line in neighbor_joining_matrix: 
    #     universal_utils.print_arr(line)

    tree = neighbor_joining_generate_tree(distance_matrix)
    Tree.show_tree_in_adjacent_list(tree, 'f', mapping)