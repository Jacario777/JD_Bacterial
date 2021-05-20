# Firstname Lastname
# NetID
# COMP 182 Spring 2021 - Homework 5, Problem 4

# You can import any standard library, as well as Numpy and Matplotlib.
# You can use helper functions from provided.py, and autograder.py,
# but they have to be copied over here.

# Your code here...

from graphviz import Digraph
from typing import Tuple
from collections import *
from copy import *

def compute_genetic_distance(seq1, seq2):
    """
    Input: seq1 and seq2, two sequences of binary integers and equal lengths
    Output: the Hamming distance between seq1 and seq2
    Calculate the hamming distance between two sequences of integers
    """
    hamming = 0
    for pos in range(len(seq1)):
        if seq1[pos] != seq2[pos]:
            hamming += 1
    return hamming

"""
#test code:
#print(compute_genetic_distance([0,1,0,1], [1,1,1,1]))
seq1 = [1,0,1,1,0,1]
seq2 = [0,1,0,1,0,0]
#print(compute_genetic_distance(seq1, seq2))
#print(compute_genetic_distance([],[]))
"""

def construct_complete_weighted_digraph(gene_data, epi_data):
    """
    Input: gene_data, epi_data, the filename for the genetic and epidemiological data
    Output: a digraph whose nodes are patients and edges are based on equation 1
    """
    seq1 = read_patient_sequences('patient_sequences.txt')
    pair_epi_dict = read_patient_traces('patient_traces.txt')
    gene_dict = compute_pairwise_gen_distances(seq1, compute_genetic_distance)
    digraph = {}
    total_dis = []
    #To find the maximum distance in pair_epi_dict and to set up digraph
    for origin in pair_epi_dict:
        digraph[origin] = {}
        for end in pair_epi_dict[origin]:
            total_dis.append(pair_epi_dict[origin][end])
    max_dis = max(total_dis)
    #do the digraph calculation
    for node1 in pair_epi_dict:
        for node2 in pair_epi_dict[node1]:
            inner_di = ((999 * (pair_epi_dict[node1][node2] / max_dis)) / 10**5)
            digraph[node1][node2] = gene_dict[node1][node2] + inner_di
    return digraph



def infer_transmap(gen_data, epi_data, patient_id):
    """
        Infers a transmission map based on genetic
        and epidemiological data rooted at patient_id
        
        Arguments:
        gen_data -- filename with genetic data for each patient
        epi_data -- filename with epidemiological data for each patient
        patient_id -- the id of the 'patient 0'
        
        Returns:
        The most likely transmission map for the given scenario as the RDMST 
        of a weighted, directed, complete digraph
        """
    
    complete_digraph = construct_complete_weighted_digraph(gen_data, epi_data)
    return compute_rdmst(complete_digraph, patient_id)

def read_patient_sequences(filename):
    """
        Turns the bacterial DNA sequences (obtained from patients) into a list containing tuples of
        (patient ID, sequence).
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A list of (patient ID, sequence) tuples.
        """
    sequences = []
    with open(filename) as f:
        line_num = 0
        for line in f:
            if len(line) > 5:
                patient_num, sequence = line.split("\t")
                sequences.append( (int(patient_num), ''.join(e for e in sequence if e.isalnum())) )
    return sequences

def read_patient_traces(filename):
    """
        Reads the epidemiological data file and computes the pairwise epidemiological distances between patients
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A dictionary of dictionaries where dict[i][j] is the
        epidemiological distance between i and j.
    """
    trace_data = []
    patient_ids = []
    first_line = True
    with open(filename) as f:
        for line in f:
            if first_line:
                patient_ids = line.split()
                patient_ids = list(map(int, patient_ids))
                first_line = False
            elif len(line) > 5:
                trace_data.append(line.rstrip('\n'))
    return compute_pairwise_epi_distances(trace_data, patient_ids)

def compute_pairwise_gen_distances(sequences, distance_function):
    """
        Computes the pairwise genetic distances between patients (patients' isolate genomes)
        
        Arguments:
        sequences -- a list of sequences that correspond with patient id's
        distance_function -- the distance function to apply to compute the weight of the 
        edges in the returned graph
        
        Returns:
        A dictionary of dictionaries where gdist[i][j] is the
        genetic distance between i and j.
        """
    gdist = {}
    cultures = {}
    
    # Count the number of differences of each sequence
    for i in range(len(sequences)):
        patient_id = sequences[i][0]
        seq = sequences[i][1]
        if patient_id in cultures:
            cultures[patient_id].append(seq)
        else:
            cultures[patient_id] = [seq]
            gdist[patient_id] = {}
    # Add the minimum sequence score to the graph
    for pat1 in range(1, max(cultures.keys()) + 1):
        for pat2 in range(pat1 + 1, max(cultures.keys()) + 1):
            min_score = float("inf")
            for seq1 in cultures[pat1]:
                for seq2 in cultures[pat2]:
                    score = distance_function(seq1, seq2)
                    if score < min_score:
                        min_score = score
            gdist[pat1][pat2] = min_score
            gdist[pat2][pat1] = min_score
    return gdist



### HELPER FUNCTIONS. ###

def find_first_positives(trace_data):
    """
        Finds the first positive test date of each patient
        in the trace data.
        Arguments:
        trace_data -- a list of data pertaining to location
        and first positive test date
        Returns:
        A dictionary with patient id's as keys and first positive
        test date as values. The date numbering starts from 0 and
        the patient numbering starts from 1.
        """
    first_pos = {}
    for pat in range(len(trace_data[0])):
        first_pos[pat + 1] = None
        for date in range(len(trace_data)):
            if trace_data[date][pat].endswith(".5"):
                first_pos[pat + 1] = date
                break
    return first_pos



def compute_epi_distance(pid1, pid2, trace_data, first_pos1, first_pos2, patient_ids):
    """
        Computes the epidemiological distance between two patients.
        
        Arguments:
        pid1 -- the assumed donor's index in trace data
        pid2 -- the assumed recipient's index in trace data
        trace_data -- data for days of overlap and first positive cultures
        first_pos1 -- the first positive test day for pid1
        first_pos2 -- the first positive test day for pid2
        patient_ids -- an ordered list of the patient IDs given in the text file
        
        Returns:
        Finds the epidemiological distance from patient 1 to
        patient 2.
        """
    first_overlap = -1
    assumed_trans_date = -1
    pid1 = patient_ids.index(pid1)
    pid2 = patient_ids.index(pid2)
    # Find the first overlap of the two patients
    for day in range(len(trace_data)):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            first_overlap = day
            break
    if (first_pos2 < first_overlap) | (first_overlap < 0):
        return len(trace_data) * 2 + 1
    # Find the assumed transmission date from patient 1 to patient 2
    for day in range(first_pos2, -1, -1):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            assumed_trans_date = day
            break
    sc_recip = first_pos2 - assumed_trans_date

    if first_pos1 < assumed_trans_date:
        sc_donor = 0
    else:
        sc_donor = first_pos1 - assumed_trans_date
    return sc_donor + sc_recip



def compute_pairwise_epi_distances(trace_data, patient_ids):
    """
        Turns the patient trace data into a dictionary of pairwise 
        epidemiological distances.
        
        Arguments:
        trace_data -- a list of strings with patient trace data
        patient_ids -- ordered list of patient IDs to expect
        
        Returns:
        A dictionary of dictionaries where edist[i][j] is the
        epidemiological distance between i and j.
        """
    edist = {}
    proc_data = []
    # Reformat the trace data
    for i in range(len(trace_data)):
        temp = trace_data[i].split()[::-1]
        proc_data.append(temp)
    # Find first positive test days and remove the indication from the data
    first_pos = find_first_positives(proc_data)
    for pid in first_pos:
        day = first_pos[pid]
        proc_data[day][pid - 1] = proc_data[day][pid - 1].replace(".5", "")
    # Find the epidemiological distance between the two patients and add it
    # to the graph
    for pid1 in patient_ids:
        edist[pid1] = {}
        for pid2 in patient_ids:
            if pid1 != pid2:
                epi_dist = compute_epi_distance(pid1, pid2, proc_data,
                                                first_pos[pid1], first_pos[pid2], patient_ids)
                edist[pid1][pid2] = epi_dist
    return edist

#For RDMST:
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
    
    #Intializes the lists to store cycle
    explored = []
    nodes_list = list(rdst_candidate.keys())
    current_node = nodes_list.pop(0)

    while current_node not in explored and len(nodes_list) != 0:
        if rdst_candidate[current_node] == {}:
            #Reset exploration because there are no more nodes to explore deadend)
            explored = []
            current_node = nodes_list.pop()
        
        else:
            #if there is a path, record every node in it
            explored.append(current_node)
            paths = list(rdst_candidate[current_node].keys())

            for node in paths:
                if node in explored:
                    current_node = node
                    break
                else:
                    current_node = paths.pop()
    
    #Make sure to remove the node that is part of the path leading into cycle (but not part of cycle)
    if len(explored) > 1:
        while current_node != explored[0]:
            explored.pop(0)
        return tuple(explored)

    #If there is no cycle, return false
    else:
        return False
    
"""
#test cases:
c1 = {1:{4:1}, 2:{1:1}, 3:{2:1}, 4:{3:1}}
#print(compute_cycle(c1))
c2 = {1:{}, 2:{1: 0}, 3:{6: 0}, 4:{3: 0}, 5:{2: 0}, 6:{4: 0}}
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
    rev_graph = reverse_digraph_representation(graph)
    expand_graph = {}
    number_in = {}
    
    #reverse the order of cycle to match standard representation
    new_cycle = []
    for pos in range(len(cycle) - 1, -1, -1):
        new_cycle.append(cycle[pos])

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
                for node3 in new_cycle:
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
                for node4 in new_cycle:
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
    """
    #reverse the order of cycle
    new_cycle = [cycle[0]]
    for pos in range(len(cycle) - 1, 0, -1):
        new_cycle.append(cycle[pos])
    """

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


#complete_graph = construct_complete_weighted_digraph("patient_sequences", "patient_traces")
rdmst_graph = infer_transmap("patient_sequences", "patient_traces", 1)
print(rdmst_graph)

#For determining if any other RDMST can exist
"""
def find_rdmst(complete, rdmst):
    min_list = {}
    for node1 in rdmst:
        for node2 in rdmst[node1]:
            edge_list = []
            for node3 in complete[node2]:
                if rdmst[node1][node2] == complete[node2][node3]:
                    edge_list.append(node3)
            min_list[node2] = edge_list
    print(min_list)
    
reverse = reverse_digraph_representation(complete_graph)
find_rdmst(reverse,rdmst_graph[0])
"""
