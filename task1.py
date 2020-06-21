"""
    Task1
    Find the User Equilibrium traffic assignment solution for the network of Sioux-Falls, using the
    Method of Successive Averages (MSA) w/ 200 iterations

    Report:
        - Link flows and link travel times in the UE solution
        - Calculate the ratio of Volume/Capacity for all links
        - Report the total system travel time and the average travel time for all cars
"""


import pandas as pd
import numpy as np
from pathlib import Path

DATA_PATH = Path(__file__).parent / 'Data'

class Data:
    def __init__(self):
        self.link_data = pd.read_csv(DATA_PATH / 'links.csv')
        self.od_matrix = np.genfromtxt(DATA_PATH / 'od_matrix_drive.csv', delimiter=',') * 1.5


class Dijkstra:
    """
        1. Let distance of start vertex from start vertex = 0
        2. Let the distance of all other vertices from start = infinity

        Repeat until all vertices visited:
            1. Visit the unvisited vertex with the smallest known distance from the start vertex
            2. For the current vertex, examine its unvisited neighbours
            3. For the current vertex, calculate the distance of each neighbor from start vertex
            4. If the calculated distance if a vertex is less than the known distance, update the
               shortest distance
            5. Update the previous vertex for each of the updated distances
            6. Add the current vertex to the list of visited vertices
    """

    def __init__(self):
        self.vertices = []
        self.values = {}

    def calculate(self):
        visited = []


if __name__ == "__main__":
    data = Data()




