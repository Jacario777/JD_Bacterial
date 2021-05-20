# Firstname Lastname
# NetID
# COMP 182 Spring 2021 - Homework 5, Problem 3

# You may NOT import anything apart from already imported libraries.
# You can use helper functions from provided.py, but they have
# to be copied over here.

from typing import Tuple
from collections import *
from copy import *

def reverse_digraph_representation(graph: dict) -> dict:
    """
    input: a graph represented by a mapping
    output: a graph with reversed representation 
    Change outdegree to indegree for each node
    """
    #add all the nodes into in_graph
    in_graph = {}
    for keys in graph:
        in_graph[keys] = {}
    #reverse the edges for each nodes
    for keys in graph:
        for keys2 in graph[keys]:
            in_graph[keys2][keys] = graph[keys][keys2]
    return in_graph

#testing graph:
#test1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
#test2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
#print(reverse_digraph_representation(test1))

def modify_edge_weights(rgraph: dict, root: int) -> None:
    """
    input: a reversed graph represented by a mapping, and an integer
    representing the root node
    output: None
    """
    for keys in rgraph:
        nbr_list = rgraph[keys]
        if len(nbr_list) > 1:
            nbr_weight = min(nbr_list.values())  
            #modify the weight of every other edge    
            for keys2 in nbr_list:
                nbr_list[keys2] = nbr_list[keys2] - nbr_weight

"""
m1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
m2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
print(modify_edge_weights(m1, 0))
"""

def compute_rdst_candidate(rgraph: dict, root: int) -> dict:
    """
    input: a reversed graph represented by a mapping, and an integer
    representing the root node
    output: an RDST candidate
    Take in a reversed graph and root, and find a rdst candidate
    """
    candidate = {}
    explored = {}
    for keys in rgraph:
        explored[keys] = 0
    for keys in rgraph:
        candidate[keys] = {}
        #ignore edges that start from the root
        if keys == root:
            continue
        min_weight = list(rgraph[keys].values())
        min_weight.sort()
        #finding the minimum weight connected to keys
        for keys2 in rgraph[keys]:
            if rgraph[keys][keys2] == min_weight[0] and explored[keys] == 0:
                candidate[keys][keys2] = min_weight[0]
                explored[keys] += 1
    return candidate

"""
#Test:
g1 = {0:{}, 1:{0:2, 4:2}, 2:{0:2, 1:2}, 3:{0:2, 2:2}, 4:{2:2, 3:2}, 5:{1:2, 3:2}}
rdst1 = compute_rdst_candidate(g1, 0)
#print(rdst1)
g2 = {0:{}, 1:{0:16, 4:0}, 2:{0:2, 1:0}, 3:{0:12, 2:0}, 4:{2:16, 3:0}, 5:{1:8, 3:0}}
rdst2 = compute_rdst_candidate(g2, 0)
#print(rdst2)
"""
def compute_cycle(rdst_candidate: dict) -> tuple:
    """
    Input: a rdst candidate graph (reversed) represented by a dictioanry
    Output: tuple showing the circuit in the graph (if there is a circuit)
    Find a cycle in the rdst candidate and return it
    """
    
    for node in rdst_candidate:
        if len(rdst_candidate[node]) != 0:
            queue = [node]
            path = []
            times_explored = {}
            #helps mark how much time a node has been explored
            for node2 in rdst_candidate:
                times_explored[node2] = 0
            #search for a loop in the graph (similar to bfs)
            while len(queue) != 0:
                node3 = queue.pop(0)
                times_explored[node3] += 1
                for nbr in rdst_candidate[node3]:
                    queue.append(nbr)
                    if nbr in path:
                        return tuple(path)
                    if times_explored[nbr] > 3:
                        path.append(nbr)
"""
#test cases:
c1 = {1:{4:1}, 2:{1:1}, 3:{2:1}, 4:{3:1}}
print(compute_cycle(c1))
c2 = {1: {}, 2: {1: 0}, 3: {6: 0}, 4: {3: 0}, 5: {2: 0}, 6: {4: 0}}
print(compute_cycle(c2))
"""

def contract_cycle(graph: dict, cycle: tuple) -> Tuple[dict, int]:
    """
    Input: a graph in standard representation, and a cycle
    Output: a contracted graph (standard), and number of the new node added
    Given a cycle, compress the cycle into a new node
    """
    contract_graph = {}
    new_node = max(graph.keys()) + 1
    #Remove all the nodes that are in cycle
    for node in graph:
        if node not in cycle:
            contract_graph[node] = {}
    contract_graph[new_node] = {}
    #Add edges into the contract_graph based on step 4a (ii)
    min1 = float('inf')
    min2 = float('inf')
    for node in graph:
        edge1 = node
        for node2 in list(graph[node].keys()):
            edge2 = node2
            if edge1 in cycle and edge2 not in cycle:
                weight1 = graph[edge1][edge2]
                #if parallel edges exist, ensure that the edge with smallest weight is kept
                if weight1 <= min1:
                    contract_graph[new_node][edge2] = graph[edge1][edge2]
                    min1 = weight1
            elif edge2 in cycle and edge1 not in cycle:
                weight2 = graph[edge1][edge2]
                #if parallel edges exist, ensure that the edge with smallest weight is kept
                if weight2 <= min2:
                    contract_graph[edge1][new_node] = graph[edge1][edge2]
                    min2 = weight2
            elif edge1 not in cycle and edge2 not in cycle:
                contract_graph[edge1][edge2] = graph[edge1][edge2]

    return (contract_graph, new_node)
"""
c1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
c2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
c3 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0}}
print(contract_cycle(c2, (1,2)))
"""
def expand_graph(graph: dict, rdst_candidate: dict, cycle: tuple, cstar: int) -> dict:
    """
    Input:
        graph: the original (uncontracted) graph in standard representation
        rdst_candidate: RDST candidate in standard representation
        cycle: tuples of nodes that are contracted
        cstar: the node that replaces the cycle
    Output: the original (uncontracted) graph
    Expanding the graph with a previously compressed cycle
    """
    """
    #reverse the order of cycle
    new_cycle = []
    for pos in range(len(cycle) - 1, -1, -1):
        new_cycle.append(cycle[pos])
    """
    rev_graph = reverse_digraph_representation(graph)
    expand_graph = {}
    number_in = {}
    #initiate expand_graph and number_in
    for node in graph:
        expand_graph[node] = {}
        number_in[node] = 0
    #if an edge contains a node in the cycle, replace it with the origin of edge
    for node in rdst_candidate:
        for node2 in rdst_candidate[node]:
            #if the starting node is cstar
            if node == cstar:
                min_edge1 = float('inf')
                for node3 in cycle:
                    if node2 in graph[node3]:
                        #find the minimum edge whose start node is in the cycle
                        if graph[node3][node2] < min_edge1:
                            min_edge1 = graph[node3][node2]
                            min_node1 = node3
                #check in degree and add the edge into expand_graph
                if number_in[node2] == 0:
                    expand_graph[min_node1][node2] = rdst_candidate[node][node2]
                    number_in[node2] += 1
            #if the ending node is cstar
            elif node2 == cstar:
                min_edge2 = float('inf')
                for node4 in cycle:
                    if node4 in graph[node]:
                        #find the minimum edge whose end node is in the cycle
                        if graph[node][node4] < min_edge2:
                            min_edge2 = graph[node][node4]
                            min_node2 = node4
                #check in degree and add the edge into expand_graph
                if number_in[min_node2] == 0:
                    expand_graph[node][min_node2] = rdst_candidate[node][node2]
                    number_in[min_node2] += 1
            #if both nodes are not cstar
            elif number_in[node2] == 0:
                expand_graph[node][node2] = graph[node][node2]
                number_in[node2] += 1
    #reverse the order of cycle
    new_cycle = [cycle[0]]
    for pos in range(len(cycle) - 1, 0, -1):
        new_cycle.append(cycle[pos])

    #add edges within the cycle back to expand_graph
    for node in new_cycle:
        for node2 in graph[node]:
            if node2 in new_cycle and number_in[node2] == 0 and graph[node][node2] == 0 and node != node2:
                expand_graph[node][node2] = graph[node][node2]
                number_in[node2] += 1
                
    return expand_graph


#test case (standard, after edge modified):
"""
g1 = {0:{1:16, 2:2, 3:12}, 1:{2:0, 5:8}, 2:{3:0, 4:16}, 3:{4:0, 5:0}, 4:{1:0}, 5:{}}
cycle = [1, 2, 3, 4]

g1_con, new_node = contract_cycle(g1, cycle)
g1_rev = reverse_digraph_representation(g1)
#g1_mod = modify_edge_weights(g1, 0)
rdst_cand = compute_rdst_candidate(g1_rev, 0)
#print(rdst_cand)
#print(expand_graph(g1, g1_con, cycle, new_node))
"""
def compute_rdmst(graph, root):
    """
        This function checks if:
        (1) root is a node in digraph graph, and
        (2) every node, other than root, is reachable from root
        If both conditions are satisfied, it calls compute_rdmst_helper
        on (graph, root).
        
        Since compute_rdmst_helper modifies the edge weights as it computes,
        this function reassigns the original weights to the RDMST.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node id.
        
        Returns:
        An RDMST of graph rooted at r and its weight, if one exists;
        otherwise, nothing.
    """
    
    if root not in graph:
        print ("The root node does not exist")
        return
    
    distances = bfs(graph, root)
    for node in graph:
        if distances[node] == float('inf'):
            print ("The root does not reach every other node in the graph")
            return

    rdmst = compute_rdmst_helper(graph, root)
    
    # reassign the original edge weights to the RDMST and computes the total
    # weight of the RDMST
    rdmst_weight = 0
    for node in rdmst:
        for nbr in rdmst[node]:
            rdmst[node][nbr] = graph[node][nbr]
            rdmst_weight += rdmst[node][nbr]

    return (rdmst,rdmst_weight)



def bfs(graph, startnode):
    """
        Perform a breadth-first search on digraph graph starting at node startnode.
        
        Arguments:
        graph -- directed graph
        startnode - node in graph to start the search from
        
        Returns:
        The distances from startnode to each node
    """
    dist = {}
    
    # Initialize distances
    for node in graph:
        dist[node] = float('inf')
    dist[startnode] = 0
    
    # Initialize search queue
    queue = deque([startnode])
    
    # Loop until all connected nodes have been explored
    while queue:
        node = queue.popleft()
        for nbr in graph[node]:
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                queue.append(nbr)
    return dist

def compute_rdmst_helper(graph,root):
    """
        Computes the RDMST of a weighted digraph rooted at node root.
        It is assumed that:
        (1) root is a node in graph, and
        (2) every other node in graph is reachable from root.
        
        Arguments:
        graph -- a weighted digraph in standard dictionary representation.
        root -- a node in graph.
        
        Returns:
        An RDMST of graph rooted at root. The weights of the RDMST
        do not have to be the original weights.
        """
    
    # reverse the representation of graph
    rgraph = reverse_digraph_representation(graph)

    
    # Step 1 of the algorithm
    modify_edge_weights(rgraph, root)

    
    # Step 2 of the algorithm
    rdst_candidate = compute_rdst_candidate(rgraph, root)

    
    # compute a cycle in rdst_candidate
    cycle = compute_cycle(rdst_candidate)

    
    # Step 3 of the algorithm
    if not cycle:
        return reverse_digraph_representation(rdst_candidate)
    else:
        # Step 4 of the algorithm
        
        g_copy = deepcopy(rgraph)
        g_copy = reverse_digraph_representation(g_copy)
        
        # Step 4(a) of the algorithm
        (contracted_g, cstar) = contract_cycle(g_copy, cycle)
        #cstar = max(contracted_g.keys())
        
        # Step 4(b) of the algorithm
        new_rdst_candidate = compute_rdmst_helper(contracted_g, root)
        
        # Step 4(c) of the algorithm
        rdmst = expand_graph(reverse_digraph_representation(rgraph), new_rdst_candidate, cycle, cstar)
        
        return rdmst

g0 = {0: {1: 2, 2: 2, 3: 2}, 1: {2: 2, 5: 2}, 2: {3: 2, 4: 2}, 3: {4: 2, 5: 2}, 4: {1: 2}, 5: {}}
# Results for compute_rdmst(g0, 0):
# ({0: {1: 2, 2: 2, 3: 2}, 1: {5: 2}, 2: {4: 2}, 3: {}, 4: {}, 5: {}}, 10)

g1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}
# Results for compute_rdmst(g1, 0):
# ({0: {2: 4}, 1: {}, 2: {3: 8}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}, 28)

g2 = {0: {1: 5, 2: 4}, 1: {2: 2}, 2: {1: 2}}
# Results for compute_rdmst(g2, 0):
# ({0: {2: 4}, 1: {}, 2: {1: 2}}, 6)

g3 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0}}
# Results for compute_rdmst(g3, 1):
# ({1: {3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {}, 5: {}}, 11.1)

g4 = {1: {2: 2.1, 3: 1.0, 4: 9.1, 5: 1.1, 6: 10.1, 7: 10.1, 8: 6.1, 9: 11.0, 10: 10.1}, 2: {1: 2.1, 3: 1.0, 4: 17.0, 5: 1.0, 6: 18.1, 7: 18.1, 8: 14.1, 9: 19.1, 10: 18.0}, 3: {1: 1.0, 2: 1.0, 4: 16.0, 5: 0.0, 6: 17.0, 7: 17.0, 8: 13.1, 9: 18.1, 10: 17.0}, 4: {1: 9.1, 2: 17.1, 3: 16.0, 5: 16.0, 6: 5.1, 7: 5.1, 8: 15.1, 9: 6.1, 10: 5.0}, 5: {1: 1.1, 2: 1.0, 3: 0.0, 4: 16.0, 6: 17.1, 7: 17.1, 8: 13.1, 9: 18.1, 10: 17.0}, 6: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 7: 0.0, 8: 16.1, 9: 7.1, 10: 0.0}, 7: {1: 10.1, 2: 18.1, 3: 17.0, 4: 5.1, 5: 17.1, 6: 0.0, 8: 16.0, 9: 7.1, 10: 0.0}, 8: {1: 6.1, 2: 14.1, 3: 13.1, 4: 15.1, 5: 13.1, 6: 16.1, 7: 16.0, 9: 17.1, 10: 16.1}, 9: {1: 11.1, 2: 19.1, 3: 18.1, 4: 6.1, 5: 18.1, 6: 7.1, 7: 7.1, 8: 17.1, 10: 7.0}, 10: {1: 10.1, 2: 18.1, 3: 17.1, 4: 5.1, 5: 17.0, 6: 0.0, 7: 0.0, 8: 16.1, 9: 7.0}}
# Results for compute_rdmst(g4, 1):
# ({1: {8: 6.1, 3: 1.0, 4: 9.1}, 2: {}, 3: {2: 1.0, 5: 0.0}, 4: {9: 6.1, 10: 5.0}, 5: {}, 6: {7: 0.0}, 7: {}, 8: {}, 9: {}, 10: {6: 0.0}}, 28.3)

#print(compute_rdmst(g2, 0))
#print(g4)
#print(compute_rdmst(g3, 1))
#print(compute_rdmst(g4, 1))