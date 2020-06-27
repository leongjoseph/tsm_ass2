"""
    Task1
    Find the User Equilibrium traffic assignment solution for the network of Sioux-Falls, using the
    Method of Successive Averages (MSA) w/ 200 iterations

    Report:
        - Link flows and link travel times in the UE solution
        - Calculate the ratio of Volume/Capacity for all links
        - Report the total system travel time and the average travel time for all cars

    Procedure:
        - Determine the shortest path for each OD pair using DIJKSTRA (semi-complete)
        - Assign total demand (OD matrix) to each of these paths
"""


import pandas as pd
import numpy as np
from pathlib import Path
from pprint import pprint

DATA_PATH = Path(__file__).parent / 'Data'


class Data:
    def __init__(self):
        self.link_data = self._open_and_generate_link_data()
        self.od_matrix = np.genfromtxt(DATA_PATH / 'od_matrix_drive.csv', delimiter=',')
        self.update_tt({link: 0 for link in self.link_data['link_id']})
        self.adjacency_list = self._get_adjacency_list()
        self.verticies = self._get_verticies()

    def _open_and_generate_link_data(self):
        raw_link_data = pd.read_csv(DATA_PATH / 'links.csv')
        raw_link_data['actual_tt'] = None

        return raw_link_data

    def update_tt(self, network_load: dict) -> None:
        for link_id, link_load in network_load.items():
            actual_tt = self._generate_tt(link_id, link_load)
            np.where(self.link_data['link_id'] == link_id)
            self.link_data.loc[self.link_data['link_id'] == link_id, 'actual_tt'] = actual_tt

    def _generate_tt(self, link_id, link_volume):
        link_data = self.link_data.loc[self.link_data['link_id'] == link_id]
        link_capacity = link_data['capacity'].item()
        link_fft = link_data['free_flow_time'].item()
        link_actual_tt = self._generate_actual_tt(link_fft, link_volume, link_capacity)

        return link_actual_tt

    def _generate_actual_tt(self, link_fft, link_volume, link_capacity):
        link_actual_tt = link_fft * (1 + 0.15 * (link_volume / link_capacity) ** 4)

        return link_actual_tt

    def _get_adjacency_list(self):
        # 'adjacency' refers to the ability to travel TO a vertex FROM another vertex
        adjacency_list = {}

        for link_data in self.link_data.itertuples():
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

        for link_data in self.link_data.itertuples():
            starting_node = link_data[2]
            ending_node = link_data[3]

            verticies.add(starting_node)
            verticies.add(ending_node)

        return list(verticies)


class Dijkstra(Data):
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
        Data.__init__(self)

    def calculate_tree(self, origin_vertex):
        unvisited_verticies = self.verticies.copy()
        visited_verticies = []
        tt_from_origin = 0
        current_path = []

        # initialise dictionary for values = {vertex: [min_travel_time, prior_vertex]}
        self.known_values = {vertex: {'tt': np.inf,
                                 'path': None} for vertex in unvisited_verticies if vertex is not origin_vertex}

        self._calculate_update_values(origin_vertex, tt=0, path=None)

        def _calculate_recursion(starting_vertex, tt_from_origin, current_path):
            visited_verticies.append(starting_vertex)
            unvisited_verticies.remove(starting_vertex)
            current_path.append(starting_vertex)

            # Get neighbours of the current vertex
            neighbours = self.adjacency_list.get(starting_vertex)
            neighbours_values = {}

            for neighbour in neighbours:
                if neighbour not in visited_verticies:
                    tt = self._calculate_get_tt(starting_vertex, neighbour)
                    
                    neighbours_values.update({
                        neighbour: tt
                        })

                    current_value = self.known_values.get(neighbour)

                    if tt + tt_from_origin < current_value['tt']:
                        self._calculate_update_values(neighbour,
                                                      tt + tt_from_origin,
                                                      current_path.copy())

            # Condition for no unvisited neighbours
            if len(neighbours_values) != 0:
                closest_vertex = min(neighbours_values, key=neighbours_values.get)
                closest_vertex_tt = min(neighbours_values.values())
                tt_from_origin += closest_vertex_tt
                
                return _calculate_recursion(closest_vertex, tt_from_origin, current_path)
        
        _calculate_recursion(origin_vertex, tt_from_origin, current_path)

        return self.known_values

    def _calculate_update_values(self, dest_vertex, tt, path):
        self.known_values.update({dest_vertex: {
            'tt': tt,
            'path': path
        }})

    def _calculate_get_tt(self, orig_vertex, dest_vertex):
        travel_time = self.link_data['actual_tt'].loc[
            (self.link_data['init_node'] == orig_vertex) & (self.link_data['term_node'] == dest_vertex)].item()
        
        return travel_time


class MSA(Dijkstra):
    def __init__(self):
        super().__init__()

    def _get_route(self, orig, dest):
        tree = self.calculate_tree(orig)
        dest_values = tree.get(dest)
        path_to_dest = dest_values['path']
        path_links = self._get_route_convert_to_link(path_to_dest, dest)

    def _get_route_convert_to_link(self, path, dest):
        if len(path) == 1:
            path_links = [tuple([path[0], dest])]
        else:
            path_links = []
            for idx, vertex in enumerate(path):
                if idx == 0:
                    pass
                elif idx == len(path) - 1:
                    path_links.append(tuple([path[idx - 1], path[idx]]))
                    path_links.append(tuple([path[idx], dest]))
                else:
                    path_links.append(tuple([path[idx - 1], vertex]))

        return path_links

    def _assign_demand(self):
        pass


if __name__ == "__main__":
    msa = MSA()
    msa._get_route(1, 13)



