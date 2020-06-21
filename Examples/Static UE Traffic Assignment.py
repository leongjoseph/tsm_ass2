#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import time
import pandas as pd
links = pd.read_csv('links.csv')


# In[ ]:


links.head()


# In[ ]:


nodes = [] # initiate an empty list for nodes
adjacencyLists = {} # initiate an empty dictionaries for adjencency lists
for i in range(len(links)):
    tmpFromNode = links['init_node'][i]
    tmpToNode = links['term_node'][i]
    if tmpFromNode not in nodes:
        nodes.append(tmpFromNode)
        adjacencyLists[tmpFromNode] = []
    adjacencyLists[tmpFromNode].append(i)


# In[ ]:


def myDijkstra(_nodes, _links, _flows, _adjacencyLists, _origin):
    labels = {}
    predecessors = {}
    upLinks = {}
    for i in _nodes:
        labels[i] = 99999999
    labels[_origin] = 0
    predecessors[_origin] = -1
    upLinks[_origin] = -1
    searchList = _nodes.copy()
    while len(searchList) > 0:
        smallestLabel = 99999999
        for i in searchList:
            if labels[i] < smallestLabel:
                smallestLabel = labels[i]
                smallestLabelNode = i
        searchList.remove(smallestLabelNode)
        for i in _adjacencyLists[smallestLabelNode]:
            tmpFromNode = smallestLabelNode
            tmpToNode = _links['term_node'][i]
            tmpLinkTravelTime = _links['free_flow_time'][i]*(1+0.15*(_flows[i]/_links['capacity'][i])**4)
            if labels[tmpToNode] > labels[tmpFromNode] + tmpLinkTravelTime :
                labels[tmpToNode] = labels[tmpFromNode] + tmpLinkTravelTime
                predecessors[tmpToNode] = tmpFromNode
                upLinks[tmpToNode] = i             
    output = {}
    output['SP_labels'] = labels.copy()
    output['SP_predecessors'] = predecessors.copy()
    output['SP_upLinks'] = upLinks.copy()
    return output


# In[ ]:


# read and store demand matrix in a numpy array
from numpy import genfromtxt
od_demand_matrix = genfromtxt('od_matrix_drive.csv', delimiter=',') * 1.5
od_demand_matrix[9][8]


# In[ ]:


totalDemand = 0
for i in range(len(od_demand_matrix)):
    for j in range(len(od_demand_matrix[i])):
        print ('travel demand from node (or zone)', nodes[i], 'to node (or zone)', nodes[j], 'is', od_demand_matrix[i][j])
        totalDemand += od_demand_matrix[i][j]
totalDemand


# In[ ]:


# All-or-Nothing assignment:


# In[ ]:


link_flows = {}
for i in range(len(links)):
    link_flows[i]= 0

for i in range(len(od_demand_matrix)):
    tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
    for j in range(len(od_demand_matrix[i])):
        if j != i: 
            tempNode = nodes[j] # start from destination
            while 1:
                tempPredecessor = tempShortestPathTree['SP_predecessors'][tempNode]
                tmpUpLink = tempShortestPathTree['SP_upLinks'][tempNode]
                if tempNode == nodes[i]:
                    break
                link_flows[tmpUpLink] += od_demand_matrix[i][j]
                tempNode = tempPredecessor
                
for i in range(len(link_flows)):
    print (i, links['link_id'][i], link_flows[i],  links['free_flow_time'][i],            round(links['free_flow_time'][i]*(1+0.15*(link_flows[i]/links['capacity'][i])**4),2))


# In[ ]:


total_travel_time = 0
for i in range(len(link_flows)):
    total_travel_time += link_flows[i] *  links['free_flow_time'][i] * (1+0.15*(link_flows[i]/links['capacity'][i])**4)
print(round(total_travel_time/totalDemand, 2))


# In[ ]:


# Incremental Assignment:


# In[ ]:


link_flows = {}
for i in range(len(links)):
    link_flows[i]= 0

num_increments = 100
for k in range(num_increments): # repeat by the number of increments
    for i in range(len(od_demand_matrix)):
        tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
        for j in range(len(od_demand_matrix[i])):
            if j != i: 
                tempNode = nodes[j] # start from destination
                while 1:
                    tempPredecessor = tempShortestPathTree['SP_predecessors'][tempNode]
                    tmpUpLink = tempShortestPathTree['SP_upLinks'][tempNode]
                    if tempNode == nodes[i]:
                        break
                    link_flows[tmpUpLink] += od_demand_matrix[i][j]/num_increments
                    tempNode = tempPredecessor

total_travel_time = 0
for i in range(len(link_flows)):
    total_travel_time += link_flows[i] *  links['free_flow_time'][i] * (1+0.15*(link_flows[i]/links['capacity'][i])**4)
print(round(total_travel_time/totalDemand, 2))    


# In[ ]:


# Incremental assignment with variable number of incerements:


# In[ ]:


for num_increments in range (1,20):
    time_start = time.process_time()
    link_flows = {}
    for i in range(len(links)):
        link_flows[i]= 0

    for k in range(num_increments):
        for i in range(len(od_demand_matrix)):
            tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
            for j in range(len(od_demand_matrix[i])):
                if j != i: 
                    tempNode = nodes[j] # start from destination
                    while 1:
                        tempPredecessor = tempShortestPathTree['SP_predecessors'][tempNode]
                        tmpUpLink = tempShortestPathTree['SP_upLinks'][tempNode]
                        if tempNode == nodes[i]:
                            break
                        link_flows[tmpUpLink] += od_demand_matrix[i][j]/num_increments
                        tempNode = tempPredecessor

# Compute total travel time from the solution in hand 
    total_travel_time = 0
    for i in range(len(link_flows)):
        total_travel_time += link_flows[i] *  links['free_flow_time'][i] * (1+0.15*(link_flows[i]/links['capacity'][i])**4)
    time_elapsed = (time.process_time() - time_start)
    print('num_increment=', num_increments, 'avg_travel_time= ', round(total_travel_time/totalDemand, 1),          'min comp_time= ',time_elapsed, 'sec') 


# In[ ]:


# Method of Successive Averages (MSA):


# In[ ]:


path_set = {}
path_link_set = {}
link_flows = {}
path_flows = {}
for i in range(len(links)):
    link_flows[i]= 0

for msa_iter in range (1,100):
    for i in range(len(od_demand_matrix)):
        if nodes[i] not in path_set.keys():
            path_set[nodes[i]] = {}
        tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
        for j in range(len(od_demand_matrix[i])):
            if j != i:
                
                # Take away 1/n of the existing pathflow from previous iterations
                if nodes[j] not in path_set[nodes[i]].keys():
                    path_set[nodes[i]][nodes[j]]= set()
                for previous_path in path_set[nodes[i]][nodes[j]]:
                    if previous_path in path_flows.keys():
                        path_flows[previous_path] = path_flows[previous_path] * (msa_iter - 1)/msa_iter 
                    
                # Add the new path
                newPathKey = 'dest'
                newPathLinks = []        
                tempNode = nodes[j] # start from destination and chase the predecessors to create the new path
                while tempNode != nodes[i]:
                    tempPredecessor = tempShortestPathTree['SP_predecessors'][tempNode]
                    tmpUpLink = tempShortestPathTree['SP_upLinks'][tempNode]
                    newPathKey =  str(links['link_id'][tmpUpLink]) + '-' + newPathKey
                    newPathLinks.insert(0,tmpUpLink)
                    tempNode = tempPredecessor
                path_set[nodes[i]][nodes[j]].add(newPathKey)
                path_link_set[newPathKey] = newPathLinks
                
                # Add new flow to current shortest paths
                if newPathKey not in path_flows.keys():
                    path_flows[newPathKey] = od_demand_matrix[i][j]/(msa_iter)
                else:
                    path_flows[newPathKey] += od_demand_matrix[i][j]/(msa_iter)
                    
    # compute link flows from updated path flows    
    for i in range(len(links)):
        link_flows[i]= 0
    for path in path_flows:
        for pathLink in path_link_set[path]:
            link_flows[pathLink] += path_flows[path]
    total_travel_time = 0
    for i in range(len(link_flows)):
        total_travel_time += link_flows[i] *  links['free_flow_time'][i] * (1+0.15*(link_flows[i]/links['capacity'][i])**4)
    print('iteration#=', msa_iter, 'avg_travel_time= ', round(total_travel_time/totalDemand, 1))   


# In[ ]:


# Method of Successive Averages (MSA) with convergence criteria


# In[ ]:


path_set = {}
path_link_set = {}
link_flows = {}
path_flows = {}
path_times = {}
for i in range(len(links)):
    link_flows[i]= 0
for msa_iter in range (1,100):
    for i in range(len(od_demand_matrix)):
        if nodes[i] not in path_set.keys():
            path_set[nodes[i]] = {}
        tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
        for j in range(len(od_demand_matrix[i])):
            if j != i:
                # Remove part of the existing pathflow from previous iterations
                if nodes[j] not in path_set[nodes[i]].keys():
                    path_set[nodes[i]][nodes[j]]= set()
                for previous_path in path_set[nodes[i]][nodes[j]]:
                    if previous_path in path_flows.keys():
                        path_flows[previous_path] = path_flows[previous_path] * (msa_iter - 1)/msa_iter 
                    
                # Add the new path
                newPathKey = 'dest'
                newPathLinks = []  
                tempNode = nodes[j] # start from destination and chase the predecessors to create the new path
                while tempNode != nodes[i]:
                    tempPredecessor = tempShortestPathTree['SP_predecessors'][tempNode]
                    tmpUpLink = tempShortestPathTree['SP_upLinks'][tempNode]
                    newPathKey =  str(links['link_id'][tmpUpLink]) + '-' + newPathKey
                    newPathLinks.insert(0,tmpUpLink)
                    tempNode = tempPredecessor
                path_set[nodes[i]][nodes[j]].add(newPathKey)
                path_link_set[newPathKey] = newPathLinks
                # Add new flow to current shortest paths
                if newPathKey not in path_flows.keys():
                    path_flows[newPathKey] = od_demand_matrix[i][j]/(msa_iter)
                else:
                    path_flows[newPathKey] += od_demand_matrix[i][j]/(msa_iter)
                    
    # compute link flows from updated path flows    
    for i in range(len(links)):
        link_flows[i]= 0
    for path in path_flows:
        for pathLink in path_link_set[path]:
            link_flows[pathLink] += path_flows[path]
            
    # compute convergence measures--
    total_time_difference = 0
    for i in range(len(od_demand_matrix)):
        tempShortestPathTree = myDijkstra(nodes, links, link_flows, adjacencyLists, nodes[i])
        for j in range(len(od_demand_matrix[i])):
            if j != i:
                minODTime = tempShortestPathTree['SP_labels'][nodes[j]]
                OD_time_difference = 0
                for pathkey in path_set[nodes[i]][nodes[j]]:
                    pathTravelTime = 0
                    for pathLink in path_link_set[pathkey]:
                        pathTravelTime += links['free_flow_time'][pathLink]*                        (1+0.15*(1.0*link_flows[pathLink]/links['capacity'][pathLink])**4)
                    time_difference = pathTravelTime - minODTime
                    OD_time_difference += path_flows[pathkey] * time_difference
                total_time_difference += OD_time_difference
    total_travel_time = 0
    for i in range(len(link_flows)):
        total_travel_time += link_flows[i] *  links['free_flow_time'][i] * (1+0.15*(link_flows[i]/links['capacity'][i])**4)
    print('iteration#=', msa_iter, 'avg_travel_time= ', round(total_travel_time/totalDemand, 1),           'convergence criteria =', total_time_difference/totalDemand)


# In[ ]:




