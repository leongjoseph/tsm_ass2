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
from pprint import pprint

DATA_PATH = Path(__file__).parent / 'Data'

class Data:
    def __init__(self):
        self.link_data = pd.read_csv(DATA_PATH / 'links.csv')
        self.od_matrix = np.genfromtxt(DATA_PATH / 'od_matrix_drive.csv', delimiter=',') * 1.5
    
    def generate_tt(self):
        pass
    
    def _generate_tt_equation(self, t_freeflow, v, U):
        t_actual = t_freeflow * (1 + 0.15 * ( v / U ) ** 4)
        
        return t_actual
    

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

    def __init__(self, data):
        self.data = data
        self.values = {}
        self.verticies = self._get_verticies()
        self.adjacency_list = self._get_adjacency_list()
    
    def _initialise_travel_times(self, data):
        pass
        
    def _get_adjacency_list(self):
        # 'adjacency' refers to the ability to travel TO a vertex FROM another vertex
        adjacency_list = {}
        
        for link_data in self.data.itertuples():
            starting_node = link_data[2]
            ending_node = link_data[3]
            
            current_list = adjacency_list.get(starting_node)
            
            if current_list is None:
                adjacency_list.update({
                    starting_node: [ending_node]
                    })
            else:
                current_list.append(ending_node)
        
        return adjacency_list
 
    def _get_verticies(self):
        verticies = set()
        
        for link_data in self.data.itertuples():
            starting_node = link_data[2]
            ending_node = link_data[3]
            
            verticies.add(starting_node)
            verticies.add(ending_node)
        
        return list(verticies)

    def calculate(self, origin_vertex):
        unvisited_verticies = self.verticies.copy()
        visited_verticies = []
        tt_from_origin = 0
        
        # initialise dictionary for values = {vertex: [min_travel_time, prior_vertex]}
        known_values = {vertex: [np.inf, None] for vertex in unvisited_verticies if vertex is not origin_vertex}
        known_values.update({origin_vertex: [0, None]})
        
        def _calculate_recursion(starting_vertex, tt_from_origin):
            if len(unvisited_verticies) == 0:
                return
            
            visited_verticies.append(starting_vertex)
            unvisited_verticies.remove(starting_vertex)
            
            neighbours = self.adjacency_list.get(starting_vertex)
            neighbours_values = {}
            
            for neighbour in neighbours:
                if neighbour not in visited_verticies:
                    tt = self._calculate_get_tt(starting_vertex, neighbour)
                    
                    neighbours_values.update({
                        neighbour: tt
                        })
                    
                    current_value = known_values.get(neighbour)
                    
                    if tt < current_value[0]:
                        known_values.update({
                            neighbour: [tt, starting_vertex]
                            })
            
            if len(neighbours_values) != 0:
                closest_vertex = min(neighbours_values, key=neighbours_values.get)
                closest_vertex_tt = min(neighbours_values.values())
                tt_from_origin += closest_vertex_tt
                
                return _calculate_recursion(closest_vertex, tt_from_origin)
        
        _calculate_recursion(origin_vertex, tt_from_origin)
        
        return known_values
            
    def _calculate_get_tt(self, orig_vertex, dest_vertex):
        travel_time = self.data['actual_tt'].loc[
            (self.data['init_node'] == orig_vertex) & (self.data['term_node'] == dest_vertex)].item()
        
        return travel_time
        
    def _calculate_min_tt(self, starting_vertex, neighbours):
        travel_times = {}
        
        for neighbour in neighbours:
            travel_time = self._calculate_get_tt(starting_vertex, neighbour)
            
            travel_times.update({
                neighbour: travel_time
            })
        
        return min(travel_times, key=travel_times.get), min(travel_times.values())
    
    def calculate_shortest_route(self, origin, destination):
        pass
        

    
    

if __name__ == "__main__":
    data = Data()
    dijkstra = Dijkstra(data.link_data)
    


