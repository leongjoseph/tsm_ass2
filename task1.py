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
    """
    DESCRIPTION:
    Used for the access/storage/generation of data
    """
    def __init__(self):
        self.link_data = self._open_and_generate_link_data()
        self.od_matrix = np.genfromtxt(DATA_PATH / 'od_matrix_drive.csv', delimiter=',')
        self.current_demand = {link_id: 0 for link_id in self.link_data['link_id']}
        self.update_tt(self.current_demand)
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
        link_actual_tt = self._generate_tt_vdf(link_fft, link_volume, link_capacity)

        return link_actual_tt

    def _generate_tt_vdf(self, link_fft, link_volume, link_capacity):
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
        super().__init__()

    def calc_sptree(self, starting_vertex):
        # Correct method was used prior, but need to include it to calculate for 
        # EACH OD pair i.e. searches for the destination

        unvisited_verticies = self.verticies.copy()

        verticies_tts = {vertex: np.inf for vertex in unvisited_verticies} 
        verticies_paths = {vertex: [] for vertex in unvisited_verticies}

        verticies_tts.update({starting_vertex: 0})
        verticies_paths.update({starting_vertex: []})

        while len(unvisited_verticies) != 0:
            u = self._calc_sptree_find_min_unvisited_dist(unvisited_verticies, verticies_tts)
            u_tt = verticies_tts.get(u)
            u_path = verticies_paths.get(u)

            unvisited_verticies.remove(u)

            neighbours = self.adjacency_list.get(u)
            
            for neighbour in neighbours:
                neighbour_tt = self._calc_sptree_get_tt(u, neighbour) + u_tt
                neighbour_current_tt = verticies_tts.get(neighbour)

                if neighbour_current_tt > neighbour_tt:
                    temp_path = u_path.copy()
                    temp_path.append(neighbour)

                    verticies_tts.update({neighbour: neighbour_tt})
                    verticies_paths.update({neighbour: temp_path})

        return verticies_tts, verticies_paths

    def _calc_sptree_find_min_unvisited_dist(self, unvisited_verticies, verticies_tts):
        temp_tts = verticies_tts.copy()

        while True:
            temp_solution = min(temp_tts, key=temp_tts.get)

            if temp_solution in unvisited_verticies:
                return temp_solution
            else:
                del temp_tts[temp_solution]

    def _calc_sptree_get_tt(self, orig_vertex, dest_vertex):
        travel_time = self.link_data['actual_tt'].loc[
            (self.link_data['init_node'] == orig_vertex) & 
            (self.link_data['term_node'] == dest_vertex)].item()
        
        return travel_time


class MSA(Dijkstra):
    """
    Procedure here is to:
    - Assign all traffic demand to the best route (all-or-nothing)
    - Update the actual travel times
    - Correct the assignment by re-routing a portion (1/i) of the previously
      assigned traffic to the new best route
    - Continue until equilibrium is reached
    """
    def __init__(self):
        super().__init__()
        self.sp_trees = [self.calc_sptree(vertex) for vertex in self.verticies]
        pprint(self.sp_trees)

    def _get_route_info(self, orig, dest):
        sp_tree = self.sp_trees.get(orig)

        dest_values = sp_tree.get(dest)
        path_to_dest = dest_values['path']

        path_links_ids, path_links = self._get_route_convert_to_links(path_to_dest, dest)
        path_tt = dest_values['tt']

        return path_links, path_tt, path_links_ids

    def solve(self):
        pass

    def solve_single_iter(self):
        for orig in self.verticies:
            for dest in self.verticies:
                if orig == dest:
                    continue
                else:
                    path, tt, link_ids = self._get_route_info(orig, dest) # incomplete

    def _get_route_convert_to_links(self, path, dest):
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

        path_link_ids = self._get_link_ids(path_links)

        return path_link_ids, path_links

    def _assign_demand_to_best_route(self):
        pass

    def _get_link_ids(self, path_links):
        link_ids = []
        for link in path_links:
            link_id = self.link_data['link_id'].loc[
                        (self.link_data['init_node'] == link[0]) & 
                        (self.link_data['term_node'] == link[1])
                        ]

            link_ids.append(link_id.item())

        return link_ids


if __name__ == "__main__":
    dij = Dijkstra()
    dij.calc_sptree(21)

