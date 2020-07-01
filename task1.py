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
    def __init__(self, mode='ue'):
        self.mode = mode
        self.link_data = self._open_and_generate_link_data()
        self.od_matrix = np.genfromtxt(DATA_PATH / 'od_matrix_drive.csv', delimiter=',')
        self.link_demand = {link_id: 0 for link_id in self.link_data['link_id']}
        self.update_tt()
        self.adjacency_list = self._get_adjacency_list()
        self.verticies = self._get_verticies()

    def _open_and_generate_link_data(self):
        raw_link_data = pd.read_csv(DATA_PATH / 'links.csv')
        raw_link_data['actual_tt'] = None

        return raw_link_data

    def update_tt(self) -> None:
        network_load = self.link_demand

        for link_id, link_load in network_load.items():
            actual_tt = self._generate_tt(link_id, link_load)
            np.where(self.link_data['link_id'] == link_id)
            self.link_data.loc[self.link_data['link_id'] == link_id, 'actual_tt'] = actual_tt

    def get_actual_tt(self, link_id) -> int:
        actual_tt = self.link_data['actual_tt'].loc[self.link_data['link_id'] == link_id]

        return actual_tt.item()

    def get_link_capacity(self, link_id) -> int:
        capacity = self.link_data['capacity'].loc[self.link_data['link_id'] == link_id]

        return capacity.item()

    def _generate_tt(self, link_id, link_volume):
        link_data = self.link_data.loc[self.link_data['link_id'] == link_id]
        link_capacity = link_data['capacity'].item()
        link_fft = link_data['free_flow_time'].item()

        if self.mode == 'ue':
            link_actual_tt = self._generate_tt_vdf(link_fft, link_volume, link_capacity)
        elif self.mode == 'so':
            link_actual_tt = self._generate_tt_lmc(link_fft, link_volume, link_capacity)
        elif self.mode == 'capacity_increase':
            if link_id == 33 or link_id == 40:
                link_actual_tt = self._generate_tt_vdf(link_fft, link_volume, link_capacity * 1.25)
            else:
                link_actual_tt = self._generate_tt_vdf(link_fft, link_volume, link_capacity)

        return link_actual_tt

    def _generate_tt_vdf(self, link_fft, link_volume, link_capacity):
        link_actual_tt = link_fft * (1 + 0.15 * (link_volume / link_capacity) ** 4)

        return link_actual_tt
    
    def _generate_tt_lmc(self, link_fft, link_volume, link_capacity):
        link_actual_tt = link_fft * (1 + 0.15 * (link_volume / link_capacity) ** 4)
        link_lmc = link_actual_tt + link_volume * (link_fft * 0.6 * ((link_volume ** 3) / (link_capacity ** 4)))

        return link_lmc

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
        verticies_paths.update({starting_vertex: [starting_vertex]})

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


class PathHistory(Data):
    def __init__(self):
        super().__init__()
        self.path_history = self._init_path_history() # {od: {path: demand}}

    def _init_path_history(self):
        path_history = {}
        for vertex_orig in self.verticies:
            for vertex_dest in self.verticies:
                if vertex_orig != vertex_dest:
                    path_history.update({f'{vertex_orig}-{vertex_dest}': {}})
        
        return path_history

    def get_path_history(self, orig, dest):
        return self.path_history[f'{orig}-{dest}']

    def update_path_history(self, orig, dest, path, tt, demand):
        self.path_history[f'{orig}-{dest}'][tuple(path)] = [tt, demand]


class Report:
    def __init__(self):
        pass

    def report_volume_to_capacity(self):
        volume_to_capacity = {}

        for link_id, link_volume in self.link_demand.items():
            link_capacity = self.get_link_capacity(link_id)
            volume_to_capacity.update({link_id: link_volume / link_capacity})

        return volume_to_capacity

    def report_fft(self):
        fft = {}

        for od_pair, od_pair_values in self.ods_data.items():
            fft[od_pair] = od_pair_values['tt']

        return fft


class MSA(Dijkstra, PathHistory, Report):
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
        self.sp_trees = self._calc_sp_trees()
        self.iteration_scaling = [1 / i for i in range(1, 202)]
        self.iteration = 1

    def _calc_sp_trees(self):
        sp_trees = {}

        for vertex in self.verticies:
            vertex_tts, vertex_paths = self.calc_sptree(vertex)
            sp_trees.update({vertex: {'tt': vertex_tts,
                                      'paths': vertex_paths}})

        return sp_trees

    def _update_past_tts(self, orig, dest):
        for path, path_info in self.get_path_history(orig, dest).items():
            path_new_tt = self._calc_path_tt(path)
            path_old_demand = path_info[1]
            self.update_path_history(orig, dest, path, path_new_tt, path_old_demand)

    def _calc_path_tt(self, link_ids):
        current_path_tt = 0

        for link_id in link_ids:
            current_path_tt += self.get_actual_tt(link_id)

        return current_path_tt

    def _shift_demand(self, orig, dest):
        current_sp = self.ods_data[f'{orig}-{dest}']['link_ids']
        current_sp_tt = self.ods_data[f'{orig}-{dest}']['tt']

        if self.iteration == 1:
            od_demand = self.od_matrix[orig - 1, dest - 1]
            self.update_path_history(orig, dest, current_sp, current_sp_tt, od_demand)
            self._update_link_demand(current_sp, od_demand, add=True)
        else:
            demand_shift_frac = self.iteration_scaling[self.iteration - 1]
            od_path_history = self.get_path_history(orig, dest)

            total_load_shifted = 0

            # Add current sp into list of shortest paths with no tt, if it doesnt already exist
            if tuple(current_sp) not in od_path_history.keys():
                self.update_path_history(orig, dest, current_sp, current_sp_tt, 0)

            # For each path in history, if not the lowest tt, shift demand away
            for od_path, current_path_data in od_path_history.items():
                current_path_tt = current_path_data[0]
                current_path_load = current_path_data[1]

                # Shift away
                minus_load_shift = current_path_load * demand_shift_frac
                shifted_load = current_path_load - minus_load_shift
                total_load_shifted += minus_load_shift
                self.update_path_history(orig, dest, od_path, current_path_tt, shifted_load)
                self._update_link_demand(od_path, shifted_load)

            # Define path with the lowest travel time
            temp_od_path_history = {path: data[0] for path, data in od_path_history.items()}
            lowest_tt_path = min(temp_od_path_history, key=temp_od_path_history.get)

            # Shift into lowest tt path
            lowest_tt_path_data = od_path_history.get(lowest_tt_path) 
            lowest_tt_path_tt = lowest_tt_path_data[0]
            lowest_tt_path_load = lowest_tt_path_data[1]

            add_load_shift = lowest_tt_path_load + total_load_shifted
            self.update_path_history(orig, dest, lowest_tt_path, lowest_tt_path_tt, add_load_shift)
            self._update_link_demand(lowest_tt_path, add_load_shift, add=True)

    def _check_load_match(self, orig, dest):
        original_load = self.od_matrix[orig - 1, dest - 1]
        od_path_history = self.get_path_history(orig, dest)
        od_path_load = 0
        for od_path, od_path_data in od_path_history.items():
            od_path_load += od_path_data[1]
        print(f'OD: {orig}-{dest}') 
        pprint(od_path_history)
        print(f'original load: {original_load}') 
        print(f'current load: {od_path_load} \n') 


    def _update_link_demand(self, link_ids, demand, add=False):
        for link_id in link_ids:
            print('before: ', self.link_demand[link_id])
            if not add:
                demand_delta = self.link_demand[link_id] - demand
                self.link_demand[link_id] = demand_delta
            else:
                self.link_demand[link_id] = self.link_demand[link_id] + demand
            print('after: ', self.link_demand[link_id], demand, add, '\n')

    def _get_od_data(self): 
        ods_data = {}

        for orig in self.verticies:
            sp_tree = self.sp_trees.get(orig)

            for dest in self.verticies:
                if orig == dest:
                    continue
                path_tt = sp_tree['tt'].get(dest)
                path_link_ids, path_links = self._convert_verticies_to_links(orig, dest, sp_tree)
                od_data = {f'{orig}-{dest}': {'tt': path_tt,
                                              'link_ids': path_link_ids,
                                              'links': path_links}}

                ods_data.update(od_data)

        return ods_data

    def solve(self, iterations=2, task=None):
        """

        """
        if task == 1:
            self.mode = 'ue'
            i = 0

            while i < 200:
                self._solve_single()
                print(i)
                i += 1

            link_flows = pd.DataFrame.from_dict(self.link_demand, orient='index')
            self.link_data.to_csv(DATA_PATH / 'task1-link-data.csv')
            link_flows.to_csv(DATA_PATH / 'task1-link-flows.csv')

        elif task == 2:
            self.mode = 'capacity_increase'
            i = 0

            while i < 200:
                self._solve_single()
                i += 1

            link_flows = pd.DataFrame.from_dict(self.link_demand, orient='index')
            self.link_data.to_csv(DATA_PATH / 'task2-link-data.csv')
            link_flows.to_csv(DATA_PATH / 'task2-link-flows.csv')

        elif task == 3:
            self._solve_single()
            task_three = self.report_fft()
            task_three_df = pd.DataFrame.from_dict(task_three, orient='index', dtype=str)
            task_three_df.to_csv(DATA_PATH / 'task3-fft.csv')

        elif task == 4:
            self.mode = 'so'
            i = 0

            while i < 200:
                self._solve_single()
                i += 1

            link_flows = pd.DataFrame.from_dict(self.link_demand, orient='index')
            self.link_data.to_csv(DATA_PATH / 'task4-link-data.csv')
            link_flows.to_csv(DATA_PATH / 'task4-link-flows.csv')

        else:
            i = 0

            while i <= iterations:
                self._solve_single()
                pprint(self.link_demand)
                i += 1

    def _solve_single(self):
        self.ods_data = self._get_od_data()
        
        for idx_orig in range(len(self.verticies)):
            orig = idx_orig + 1
            for idx_dest in range(len(self.verticies)):
                if idx_orig == idx_dest:
                    continue

                dest = idx_dest + 1
                self._update_past_tts(orig, dest)
                self._shift_demand(orig, dest)
                #self._check_load_match(orig, dest)

        self.update_tt()
        self.sp_trees = self._calc_sp_trees()
        self.iteration += 1

    def _convert_verticies_to_links(self, orig, dest, sp_tree):
        paths_links = {vertex: [] for vertex in self.verticies}
        paths = sp_tree['paths'].get(dest)

        if len(paths) == 1:
            path_links = [tuple([paths[0], dest])]
        else:
            path_links = []

            for idx, vertex in enumerate(paths):
                if idx == 0:
                    pass
                else:
                    path_links.append(tuple([paths[idx - 1], vertex]))

        path_link_ids = self._get_link_id_from_data(path_links)

        return path_link_ids, path_links

    def _get_link_id_from_data(self, path_links):
        link_ids = []
        for link in path_links:
            link_id = self.link_data['link_id'].loc[
                        (self.link_data['init_node'] == link[0]) & 
                        (self.link_data['term_node'] == link[1])
                        ]

            link_ids.append(link_id.item())

        return link_ids


if __name__ == "__main__":
    msa = MSA()
    msa.solve(iterations=5)

    # Link flow and link travel times in UE solution
    # Ratio of volume to capacity for all links
    # Total system travel time and the average travel time for all cars

