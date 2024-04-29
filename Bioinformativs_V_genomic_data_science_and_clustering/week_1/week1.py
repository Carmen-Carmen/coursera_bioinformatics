import math
import numpy as np
import pandas as pd
from coursera_bioinformatics.utils.universal_utils import *
import random
import pyperclip

def get_euclidean_distance(data_point1, data_point2): 
    if len(data_point1) != len(data_point2): 
        print("Dimension inconsistency of %s and %s" %(str(data_point1), str(data_point2)))
        return None
    
    distance = 0
    for i in range(len(data_point1)): 
        distance += (data_point1[i] - data_point2[i]) ** 2
    distance = math.sqrt(distance)
    return distance

# d(DataPoint,Centers) = min all points x from Centers d(DataPoint, x).
def get_min_distance_of_data_point_to_centers(data_point, center_list): 
    min_distance = float("inf")
    for center in center_list: 
        distance = get_euclidean_distance(data_point, center)
        if distance < min_distance: 
            min_distance = distance
    
    return min_distance

# MaxDistance(Data, Centers) = max all points DataPoint from Data d(DataPoints,Centers).
def get_max_distance_of_all_data_points(data_point_list, center_list): 
    distance_list = []
    for data_point in data_point_list: 
        distance_list.append(get_min_distance_of_data_point_to_centers(data_point, center_list))
    
    return max(distance_list)

def get_data_point_with_max_distance(data_point_list, center_list): 
    max_distance = -float("inf")
    data_point_with_max_distance = None
    for data_point in data_point_list: 
        distance = get_min_distance_of_data_point_to_centers(data_point, center_list)
        if distance > max_distance: 
            max_distance = distance
            data_point_with_max_distance = data_point
    
    return data_point_with_max_distance
    
# FarthestFirstTraversal(Data, k) 
#     Centers ← the set consisting of a single randomly chosen point from Data
#     while |Centers| < k 
#         DataPoint ← the point in Data maximizing d(DataPoint, Centers) 
#         add DataPoint to Centers 
#     return Centers 
def farthest_first_traversal(data_point_list, k):
    center_list = [data_point_list[0]]
    while len(center_list) < k: 
        data_point_with_max_distance = get_data_point_with_max_distance(data_point_list, center_list)
        center_list.append(data_point_with_max_distance)
    
    return center_list

def parse_coordinates_list(coordinates_list_str): 
    coordinates_list = []
    for coordinates_str in coordinates_list_str.strip().split("\n"):
        splitted = coordinates_str.split(" ")
        coordinates = []
        for item in splitted: 
            coordinates.append(float(item))
        coordinates_list.append(coordinates)
    
    return coordinates_list

def output_coordinates(coordinates): 
    to_output = ""
    for coordinate in coordinates: 
        to_output += "%.1f " %(coordinate)
    
    return to_output.strip()

def output_coordinates_list(coordinates_list): 
    to_output = ""
    for coordinates in coordinates_list: 
        to_output += output_coordinates(coordinates) + "\n"
    return to_output.strip()

def get_squared_distortion(data_point_list, center_list): 
    total_distance = 0
    for data_point in data_point_list: 
        total_distance += get_min_distance_of_data_point_to_centers(data_point, center_list) ** 2
    
    return total_distance / len(data_point_list)

if __name__ == "__main__": 
    # d1 = DataPoint((1, 2, 3))
    # d2 = DataPoint((1, 2))
    # get_euclidean_distance(d1, d2)

    # 1.5 step 3
    # data_point_list = []
    # center_list = []
    # data_point_coordinates_list = [
    #     (1, 6), 
    #     (1, 3), 
    #     (3, 4), 
    #     (5, 6), 
    #     (5, 2), 
    #     (7, 1), 
    #     (8, 7), 
    #     (10, 3) 
    # ]
    # center_coordinates_list = [
    #     (2, 4), 
    #     (6, 7), 
    #     (7, 3)
    # ]
    # for coordinates in data_point_coordinates_list: 
    #     data_point_list.append(coordinates)
    # for coordinates in center_coordinates_list: 
    #     center_list.append(coordinates)
    # print(get_max_distance_of_all_data_points(data_point_list, center_list))

    # 1.6 step 2
#     k = 3
#     dimension = 2
#     data_point_coordinates_list = parse_coordinates_list("""0.0 0.0
# 5.0 5.0
# 0.0 5.0
# 1.0 1.0
# 2.0 2.0
# 3.0 3.0
# 1.0 2.0""")
#     # print_arr(data_point_coordinates_list)
#     data_point_list = []
#     for coordinates in data_point_coordinates_list: 
#         data_point_list.append(coordinates)
#     center_list = farthest_first_traversal(data_point_list, k)
#     output = output_coordinates_list(center_list)
#     print(output)
#     pyperclip.copy(output)
    
    # # extra dataset
    # with open("./datasets/FarthestFirstTraversal.txt") as f: 
    #     f.readline()
    #     first_line = f.readline()
    #     k = int(first_line.strip().split(" ")[0])
    #     m = int(first_line.strip().split(" ")[1])
    #     data_point_list_str = ""
    #     while True: 
    #         temp = f.readline()
    #         if "Output" in temp: 
    #             break
    #         data_point_list_str += temp
    #     data_point_list = parse_coordinates_list(data_point_list_str)
    #     # print(output_coordinates_list(data_point_list))
    #     center_list = farthest_first_traversal(data_point_list, k)
    #     output = output_coordinates_list(center_list)
    #     print(output)
    #     pyperclip.copy(output)

    # with open("./datasets/dataset_30181_2.txt") as f: 
    #     first_line = f.readline()
    #     k = int(first_line.strip().split(" ")[0])
    #     m = int(first_line.strip().split(" ")[1])
    #     data_point_list_str = f.read()
    #     data_point_list = parse_coordinates_list(data_point_list_str)
    #     center_list = farthest_first_traversal(data_point_list, k)
    #     output = output_coordinates_list(center_list)
    #     print(output)
    #     pyperclip.copy(output)

    # 1.7 step 2
    # data_point_list = [
    #     (1, 6), 
    #     (1, 3), 
    #     (3, 4), 
    #     (5, 6), 
    #     (5, 2), 
    #     (7, 1), 
    #     (8, 7), 
    #     (10, 3) 
    # ]
    # center_list1 = [
    #     (3, 4.5), 
    #     (6, 1.5), 
    #     (9, 5)
    # ]
    # center_list2 = [
    #     (5/3, 13/3), 
    #     (6.5, 6.5), 
    #     (22/3, 2)
    # ]

    # print("for 1st centers: max_distance = %.3f, distorsion = %.3f"
    #       %(
    #           get_max_distance_of_all_data_points(data_point_list, center_list1), 
    #           get_squared_distortion(data_point_list, center_list1)
    #       ))
    # print("for 2nd centers: max_distance = %.3f, distorsion = %.3f"
    #       %(
    #           get_max_distance_of_all_data_points(data_point_list, center_list2), 
    #           get_squared_distortion(data_point_list, center_list2)
    #       ))

    # 1.7 step 3
#     k = 2
#     dimension = 2
#     center_list = parse_coordinates_list("""2.31 4.55
# 5.96 9.08""")
#     data_point_list = parse_coordinates_list("""3.42 6.03
# 6.23 8.25
# 4.76 1.64
# 4.47 4.33
# 3.95 7.61
# 8.93 2.97
# 9.74 4.03
# 1.73 1.28
# 9.72 5.01
# 7.27 3.77""")
#     print("%.3f" %get_squared_distortion(data_point_list, center_list))
        
    with open("./datasets/dataset_30170_3.txt", "r") as f: 
        first_line = f.readline()
        k = int(first_line.split(" ")[0])
        dimension = int(first_line.split(" ")[1])
        centers_and_data_points_str = f.read()
        centers_str = centers_and_data_points_str.split("--------\n")[0]
        data_points_str = centers_and_data_points_str.split("--------\n")[1]
        data_point_list = parse_coordinates_list(data_points_str)
        center_list = parse_coordinates_list(centers_str)
        print("%.3f" %get_squared_distortion(data_point_list, center_list))
        print("changed on ubuntu")