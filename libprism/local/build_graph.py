#########################################################################
# File Name: build_graph.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 26 Nov 2019 20:22:50 AEDT
# according to Yulin, no trick and add an optmal function
# now, only add a fuction update_weight_for_conflict_edges
# may 21
#########################################################################
#!/bin/bash

import os
import sys
import logging
import time
import networkx as nx
import tools
from networkx.drawing.nx_pydot import write_dot
import logging
import pdb

import group_function
import matrix
import myGraph

#tricks, or remove noisy

# for edges
def check_symmetrical_connected_component(g):

    #logging.debug (" check symmetrical connected component ")
    weightEdges= {}
    dealedPair = set()
    for (n1, n2) in list(g.edges):
        w = g.get_weight(n1, n2)
        if w == 'A': # remove edges between symmetrical nodes
            g.remove_edge(n1, n2)
            continue
        
        if (n1, n2) in dealedPair or (n2, n1) in dealedPair:
            continue
        symN1 = tools.get_symmetrical_node(n1)
        symN2 = tools.get_symmetrical_node(n2)
        symW = g.get_weight(symN1, symN2)
        
        #dealedPair.add( (n1,n2) )
        dealedPair.add( (symN1,symN2) )
        wKey = w + symW  
        if wKey not in weightEdges:
            weightEdges[ wKey ] = []
        #weightEdges[wKey].append( sorted((n1, n2, symN1, symN2)) )
        weightEdges[wKey].append( (n1, n2, symN1, symN2) ) # debug should not sort
             
    if nx.number_connected_components(g) == 2:
        graphs = list( nx.connected_component_subgraphs(g) ) 
        group1 = set( graphs[0].nodes )
        group2 = set( graphs[1].nodes )
        objScore = g.size(weight='weight') # weight sum in g
        if len(group1) == len(group2):
            #logging.debug ("it's an symmetrical connected component and they phasing natrually")
            print (group1)
            print (group2) 
            print ("it's an symmetrical connected component and they phasing natrually")
            return True, group1, group2, objScore

    return weightEdges  



def debug_unuse_weight(sortedWeight, weightThreshold):

    # debug, remove in public version
    remainW = 0
    print ( "weightThreshold:", weightThreshold )
    for w in sortedWeight:
        if w <= weightThreshold:
            remainW += w
    print ("Debug: sum of those weight who smaller than threshold", remainW)
    return 

def add_pair_by_pair(graph, sortedWeight, weightEdges, weightThreshold):

    group = group_function.inital_group(graph)

    confilictSet = set()
    dealedEdge = 0
    nodeNum = graph.number_of_nodes()
    preTime = time.time()

    #debug_unuse_weight(sortedWeight, weightThreshold)
    for w in sortedWeight:
        if w<=weightThreshold: 
            logging.debug ("weight smaller than threshold stop merge, threshold is: %s" % weightThreshold)
            break

        finish =  add_pair_at_same_weight(weightEdges[w], graph, 
                group, confilictSet, nodeNum)
        #conScore += temp
        if finish:
            break
        dealedEdge += len( weightEdges[w] )
        if dealedEdge % 10000 == 0:
            curTime = time.time()
            logging.debug  ("deal %s reads cost %s" % (dealedEdge ,curTime-preTime) )
            preTime = curTime
       
    return group


def generate_finalGroups(group, graph, mergeGroupThreshold, strict, mergeType, groupSum):
    #merge group and phase unphased node

    group_function.keep_group_set2(group)
    
    '''
    #pdb.set_trace()
    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes
    print( "[debug] phased Node number: ",  len(phasedNodes))
    print( "[debug] all node number in original graph: ", graph.number_of_nodes()   )
    addScore = group_function.phased_nodes_in_groups(group, sorted(unphasedNode), graph, True)
    #pdb.set_trace()
    '''

    count = 0 
    while len(group) >2:
        if mergeType == 0:
            if strict == True:
                label = group_function.merge_group_according_threshold_strict(graph, group, mergeGroupThreshold, groupSum)
            else:    
                label = group_function.merge_group_according_threshold(graph, group, mergeGroupThreshold, groupSum)
        elif mergeType == 1:  
            label = group_function.merge_two_groups_once(graph, group , mergeGroupThreshold)
        count += 1
        if label == False: # or count >= cnt:
            break
   
    group_function.check_group(group) 
        
    return 

   
def generate_finalGroups_loose(group, graph, mergeGroupThreshold, mergeType, groupSum):
    #merge group and phase unphased node
    
    group_function.keep_group_set2(group)
   
    '''
    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes
    print( "[debug] phased Node number: ",  len(phasedNodes))
    print( "[debug] all node number in original graph: ", graph.number_of_nodes()   )
    addScore = group_function.phased_nodes_in_groups(group, sorted(unphasedNode), graph, True)
    

    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes

    print( "[debug] unphased Node number: ",  len(unphasedNode))
    for n in unphasedNode:
        assert n not in group
        group[n]=set()
        group[n].add(n)
    '''
    print ("merge two groups once a time  ")
    while len(group) >2:
        label = group_function.merge_two_groups_once(graph, group , 0.01, groupSum)
        if label == False: 
            break
     
    #group_function.keep_group_set(group)
    '''
    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes
    print( "[debug] phased Node number: ",  len(phasedNodes))
    print( "[debug] all node number in original graph: ", graph.number_of_nodes()   )
    addScore = group_function.phased_nodes_in_groups(group, sorted(unphasedNode), graph, True)

    while len(group) >2:
        #label = group_function.merge_group_according_threshold_loose(graph, group, mergeGroupThreshold, groupSum)
        label = group_function.merge_two_groups_once(graph, group , 0.005, groupSum)
        if label == False: 
            break

    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes
    print( "[debug] phased Node number: ",  len(phasedNodes))
    print( "[debug] all node number in original graph: ", graph.number_of_nodes()   )
    addScore = group_function.phased_nodes_in_groups(group, unphasedNode, graph, False)
    '''
    #pdb.set_trace()

    group_function.merge_groups_only_choice(graph, group, groupSum)

    #pdb.set_trace()
    ''' 
    phasedNodes = group_function.get_phased_nodes(group) 
    unphasedNode = set(graph.nodes() ) -  phasedNodes 
    print( "[debug] unphased Node number: ",  len(unphasedNode))
    '''
    unphasedNode = set()
    g1, g2 = build_final_graph(graph, group, groupSum, unphasedNode)
    '''
    unphasedNode = set(graph.nodes() ) - g1 - g2
    print( "[debug] unphased Node number: ",  len(unphasedNode))

    #pdb.set_trace()
    g1, g2 = group_function.phased_nodes_in_two_groups(g1, g2, unphasedNode, graph)
    '''
    g1,g2=[], []
    return g1, g2


def build_final_graph(graph, group, groupSum, unphasedNode):
   

    # add unphasedNode in group
    for n in unphasedNode:
        assert n not in group
        group[n] = set()
        group[n].add(n)

    groupKeys = sorted(group.keys()) # deep copy
    groupLen = len(group)
    #print ("before merge group: %s pairs of group" %  (groupLen/2) )

    edgesList = []
    for i in range(groupLen-1):
        if groupKeys[i].endswith('_0'):
            continue
        for j in range(i+1, groupLen):
            if groupKeys[i].endswith('_1') and groupKeys[j].endswith('_1'):
                n1, n2 = groupKeys[i], groupKeys[j]   
                symN1 = tools.get_symmetrical_node( n1 )
                symN2 = tools.get_symmetrical_node( n2 ) 
                w1, w2, w3, w4, l1, l2 = group_function.get_groups_sum(n1, n2, symN1, symN2, group, graph, groupSum)
                
                if w1 > 0: 
                    edgesList.append((n1, n2, w1))
                if w2 > 0: 
                    edgesList.append((n1, symN2, w2))
                if w3 > 0: 
                    edgesList.append((symN1, n2, w3))
                if w4 > 0: 
                    edgesList.append((symN1, symN2, w4))
                
    
    finalG = myGraph.myGraph()
    finalG.add_weighted_edges_from(edgesList)
    totalScore = finalG.displayGraph()
    finalG.add_type_A_edges()
    print ("After add A edges, number of components: ", nx.number_connected_components(finalG))

    write_dot(finalG, "finalGraph" + ".dot")
    '''  
    isolatedGroup = set(groupKeys) - set( finalG.nodes() )
    print ("isolated group", len(isolatedGroup), isolatedGroup)
    
    isolatedNodeNum = 0
    for node in isolatedGroup:
        isolatedNodeNum += len(group[node])
    
    print ("node in isolated group", isolatedNodeNum)

    #pdb.set_trace()
    finalG1Group, finalG2Group = [], [] 
    graphs= list( (finalG.subgraph(c)).copy() for c in nx.connected_components(finalG))
    for g in graphs:
        #pdb.set_trace() 
        group2 = group_function.inital_group( g )
        groupSum = {}
        bestG1, bestG2 = random_merge_and_switch_node(g, group2, groupSum)
        finalG1Group.extend(bestG1)
        finalG2Group.extend(bestG2)


    for node in isolatedGroup:
        if node.endswith('_1'):
            symN = tools.get_symmetrical_node(node)
            finalG1Group.append(node)
            finalG2Group.append(symN)
    
    finalG1node, finalG2node=set(), set()
    #pdb.set_trace()
    for n in finalG1Group:
        assert isinstance(group[n], set)
        finalG1node = finalG1node.union(group[n])
    
    for n in finalG2Group:
        assert isinstance(group[n], set)
        finalG2node = finalG2node.union(group[n])
        
    pdb.set_trace()
    '''

    finalG1node, finalG2node=set(), set()
    return finalG1node, finalG2node



    
def phasing_in_connected_component(graph, weightThreshold, mergeType):
    
    finalGroups = [] 
    nodeNum = graph.number_of_nodes() 
    if nodeNum <= 3:
        return finalGroups, 0, {}
    t1 = time.time()

    #logging.debug ("remove/update conflict edges")   
    graph.update_weight_for_conflict_edges()
    res = check_symmetrical_connected_component(graph) 
    #logging.debug("after remove 'A' edges")
    graph.displayGraph()
    t2 = time.time()
    #logging.debug("Debug: check symmetrical_connected_component cost %s" % (t2-t1) )
    
    group= {} 
    if len(res) == 4 and not isinstance(res, dict) and res[0] == True:
        k1, k2 = min(res[1]), min(res[2])
        group[ k1  ] = res[1]
        group[ k2  ] = res[2]
        assert k1[:-2] ==  k2[:-2] 
        objScore = res[3]
        return objScore, group  
    
    weightEdges = res
    sortedWeight = sorted(weightEdges.keys(), reverse=True)

    group = add_pair_by_pair(graph, sortedWeight, weightEdges, weightThreshold)
   
    
    groupSum = {}
    generate_finalGroups(group, graph, 20, True, mergeType, groupSum)
    generate_finalGroups(group, graph, 20, False, mergeType, groupSum)
    
    g1, g2 = generate_finalGroups_loose(group, graph, 0, 1, groupSum)
    objScore = graph.calculate_objScore_in_groupMap(group)
    logging.info ("after phased unphased node, final objScore is %s" % objScore)
    return objScore, group

def add_pair_at_same_weight(weightEdges, graph, group, confilictSet, nodeNum):

    #es = sorted( weightEdges ) # debug must give a order, otherwise result will change
    
    #conScore = 0 
    wOrder = {}
    finish = False
    for eleList in weightEdges:
    #for eleList in es:
        n1, n2, symN1, symN2 = eleList[0], eleList[1], eleList[2], eleList[3]
        # way 2 check ../local-back
        if group_function.is_two_nodes_merged(n1, n2, symN1, symN2, group):
            continue   
        w1 = graph.get_weight(n1, n2)
        symW = graph.get_weight(symN1, symN2)

        if n1[:-2] == n2[:-2]:
            logging.debug ("%s %s cannot merge, weight is %s" % (n1, n2, w1) )
            continue
        key = w1*symW
        if key == 0:
            continue
        if key not in wOrder:
            wOrder[key] = []
        wOrder[key].append( (n1, n2, symN1, symN2) )

    es = sorted(wOrder.items(), reverse=True)
    for (w, eleList) in es:
        fourList = sorted(eleList, reverse=True)
        '''
        print ('debug')
        if len(fourList) < 100:
            print (w, fourList)
        '''
        if w < 4:
            continue
        #sys.exit()    
        for (n1, n2, symN1, symN2) in fourList: # = eleList[0], eleList[1], eleList[2], eleList[3]
            assert n1[:-2] != n2[:-2]
         
            if group_function.is_two_nodes_merged(n1, n2, symN1, symN2, group):
                continue   
            g1, g2, r1, r2 = group_function.get_group_set(group, n1, n2)
            symG1, symG2, symR1, symR2 = group_function.get_group_set(group, symN1, symN2)
            label = group_function.update_two_groups(r1, r2, symR1, symR2, group, confilictSet)
            
            if label == True:
                # debug remove late
                '''
                w1, w2 = graph.get_weight(n1,n2), graph.get_weight(n1, symN2)
                w3, w4 = graph.get_weight(symN1, n2), graph.get_weight(symN1, symN2)
                print ("merge pair", n1, n2, w, w1, w2, w3, w4) 
                ''' 
                g1, g2, r1, r2 = group_function.get_group_set(group, n1, symN1)
                phasedNum = len( g1 ) + len( g2 ) 
                if phasedNum >= 0.9999*nodeNum:
                    finish = True
                    logging.info ("phasing %s nodes in all %s nodes, break" % (phasedNum, nodeNum) )
                    break
    return finish #, conScore

def swith_pair(): 

    '''
    for ele in groupList:
        temp = graph.subgraph(ele) 
        score = temp.displayGraph() 

    n = groupList[0][0]   
    symN = tools.get_symmetrical_node(n)
    
    groupList[0].remove(n)
    groupList[0].append(symN)
    groupList[1].remove(symN)
    groupList[1].append(n)

    for ele in groupList:
        temp = graph.subgraph(ele) 
        score = temp.displayGraph() 
    '''
    return



def build_original_graph(matrixFile):
    
    t1 = time.time()
    edges  = matrix.load_all_edges_from_matrix(matrixFile)
    t2 = time.time() 
    logging.info ("load matrix cost %s" % (t2-t1) )
    edgesList = matrix.filter_edges_according_weight(edges, 0)
    originalG = myGraph.myGraph()
    originalG.add_weighted_edges_from(edgesList)
    totalScore = originalG.displayGraph()
    originalG.add_type_A_edges()
    print ("After add A edges, number of components: ", nx.number_connected_components(originalG))
    
    # debug 
    #originalG.remove_type_A_edges()
    #totalScore = originalG.displayGraph()
    #sys.exit()

    return originalG, totalScore


def build_graph_with_more_component(matrixFile, k, wType, cov):
    
    print ("weight type", wType)
    t1 = time.time()
    edges, removeNodes = matrix.load_undirected_edges_from_matrix( matrixFile, k, wType , cov)
    t2 = time.time() 
    logging.info ("load matrix (remove self loop and self conflict) cost %s" % (t2-t1) )

    #edgesList, wTemp = matrix.filter_edges_according_rate(edges, wRate)
    wTemp = 1.0/wType # wType=1 => wTypw=1; wType=2 => wTemp=0.5
    edgesList = matrix.filter_edges_according_weight(edges, wTemp)
    t3 = time.time()
    print ("filter edges cost ", t3-t2)
    G = myGraph.myGraph()
    G.add_weighted_edges_from(edgesList)
    G.displayGraph()

    print ("remove node whose read coverage higher than homo coverage") 
    print ("node number", len(removeNodes))
    G.remove_nodes(removeNodes)
    G.displayGraph()
    
    '''
    print ("remove/update confict edge") # (1_0, 2_0) (1_0, 2_1) both exisit
    G.update_weight_for_conflict_edges()
    print ("after remove/update conflict edges")   
    G.displayGraph()
    '''
    print ("check symmetrical node")    
    G.check_symmetrical_node() # symmetrical nodes must have intersection 
    print ("after remove not symmetrical node")    
    G.displayGraph()
    dT1, degreeMost = G.get_node_degree_distribution()

    logging.debug("node degree threshold %s" % dT1)
    G.remove_high_degree_node(dT1) # tricki
    
    G.displayGraph()
    G.add_type_A_edges()
    #edgesDis = G.calculate_edge_weight_distribution()
    
    t4 = time.time()
    print ("build graph and preprocess graph cost ", t4-t3)

    return G, wTemp 

def random_merge_and_switch_node(graph, group, groupSum):

    randomCnt = 0
    bestScore = 0
    groupLen = len(group)/2
    if groupLen > 100:
        sys.exit()
    cnt = 0
    while True:
        if randomCnt >= 100 or randomCnt >= 2**groupLen:
            break
        randomCnt += 1
        mergedGroupList, nodeGroup = group_function.random_merge_all_group(group) 
        objScore = graph.calculate_objScore_in_groupList( mergedGroupList )
        print ("random merge node, objScore %s" % (objScore) )
        objScore, nodeGroup = group_function.iterate_switch_group(objScore, graph, nodeGroup, group, groupSum)
        print ("random %s time, finalScore %s" % (randomCnt, objScore))
        if bestScore < objScore:
            print ("updata best nodeGroup", sorted( nodeGroup) )
            bestG1, bestG2 =set(), set()
            cnt = 1
            for node in nodeGroup:
                symNode = tools.get_symmetrical_node(node)
                bestG1 = bestG1.union(group[node])
                bestG2 = bestG2.union(group[symNode])
            bestScore = objScore
        elif bestScore == objScore:
            cnt += 1
    print ("final bestScore", bestScore)

    print ("[debug] %s times get final bestScore" % cnt)
    # for check
    objScore = graph.calculate_objScore2(bestG1, bestG2)
    print ("after merge group, bestScore", objScore)
    '''
    bestG1, bestG2 = phased_nodes_in_two_groups(bestG1, bestG2, unphasedNode, graph)

    objScore, conScore = calculate_objScore(bestG1, bestG2, graph)

    print ("after add phased node, bestScore", objScore)
    '''
    return list(bestG1), list(bestG2)


def swith_pair(): 

    '''
    for ele in groupList:
        temp = graph.subgraph(ele) 
        score = temp.displayGraph() 

    n = groupList[0][0]   
    symN = tools.get_symmetrical_node(n)
    
    groupList[0].remove(n)
    groupList[0].append(symN)
    groupList[1].remove(symN)
    groupList[1].append(n)

    for ele in groupList:
        temp = graph.subgraph(ele) 
        score = temp.displayGraph() 
    '''
    return


logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
#logging.basicConfig(stream=sys.stderr, level=logging.INFO)

if __name__ == '__main__':

    if len(sys.argv) < 6:
        print ("matrix_file kmer_size filter_edge_rate homo_cov mergeType")
        sys.exit()

    matrixFile = sys.argv[1]
    k =int(sys.argv[2])
    wRate = float(sys.argv[3])
    cov = int(sys.argv[4])
    mergeType=int(sys.argv[5]) # debug fast way
    wType = 1

    t1 = time.time()
    originalG, totalScore = build_original_graph(matrixFile)

    ''' 
    G, wTemp = build_graph_with_more_component(matrixFile, k, wType, cov)
    # remove some nodes, edges, to get more connected component
    # not necessary

    #pdb.set_trace()    
    #graphs = list( nx.connected_component_subgraphs(G) )
    #graphs= list(G.subgraph(c) for c in nx.connected_components(G))
    #graphs= list( (G.subgraph(c)).copy() for c in nx.connected_components(G))

    graphs= list( (originalG.subgraph(c)).copy() for c in nx.connected_components(G))
    '''
    wTemp = 1.0
    graphs= list( (originalG.subgraph(c)).copy() for c in nx.connected_components(originalG) )
    logging.info("number of components: %s " % ( len(graphs) ) )
    iComponent = 0
    groupNum  = 0
    finalGroupMap = {}
    for g in graphs:
        gNodesNum = g.number_of_nodes()
        if ( gNodesNum <=3 ):
            continue
        iComponent += 1  
        #logging.debug("------------------------------------------------------------------")
        print("------------------------------------------------------------------")
        #logging.debug("phasing connected component %s node number: %s" % (iComponent, gNodesNum ) )
        print("phasing connected component %s node number: %s" % (iComponent, gNodesNum ) )
        #if gNodesNum < 100: write_dot(g, "componemt" + str(iComponent) + ".dot")

        objScore, groupMap = phasing_in_connected_component(g, wTemp, mergeType)
        
        #debug
        #pdb.set_trace()
        #group_function.check_group(groupMap)
        
        # debug, only need when just use add_pair
        #group_function.keep_group_set(groupMap)
        
        finalGroupMap.update(groupMap)
        #logging.debug("objScore is %s" % objScore)
        print("objScore is %s" % objScore)
        
        groupsLen = len(groupMap)
        groupNum += groupsLen
        #logging.debug("connected component %s range from phasing group %s %s" % (iComponent, groupNum/2 - groupsLen/2+1, groupNum/2))
        print("connected component %s range from phasing group %s %s" % (iComponent, groupNum/2 - groupsLen/2+1, groupNum/2))
    originalG.remove_type_A_edges()

    # can be remove 
    objSum = originalG.calculate_objScore_in_groupMap(finalGroupMap)
    logging.debug ( "the objSocre is : %s" % objSum )
    logging.debug ( "the number of group : %s" % len(finalGroupMap) )
 
    '''
    #pdb.set_trace()
    group_function.keep_group_set(finalGroupMap)
    phasedNodes = group_function.get_phased_nodes(finalGroupMap) 
    unphasedNode = set(originalG.nodes() ) -  phasedNodes
    print( "[debug] in whole graph, phased Node number: ",  len(phasedNodes))
    print( "[debug] all node number in original graph: ", originalG.number_of_nodes()   )
 
    for n in unphasedNode:
        assert n not in finalGroupMap
        finalGroupMap[n]=set()
        finalGroupMap[n].add(n)

    #pdb.set_trace()
    groupSum ={} 
    generate_finalGroups(finalGroupMap, originalG, 10, True, mergeType, groupSum)
    generate_finalGroups(finalGroupMap, originalG, 0, False, mergeType, groupSum)
    
    g1, g2 = generate_finalGroups_loose(finalGroupMap, originalG, 0, 1, groupSum)
 
    #random_merge_and_switch(originalG, finalGroupMap, groupSum) 
    finalList = []
    finalList.append(g1)
    finalList.append(g2)

    objSum = originalG.calculate_objScore_in_groupList(finalList)
    #originalG.remove_type_A_edges()
    
    logging.info ( "The objSocre is : %s " % objSum )
    logging.info ( "Total score in original graph: %s" %  totalScore )
    logging.info ( "Choose score in original graph rate: %.6f" % (objSum/totalScore) )

    group_function.write_groupList(finalList, "phased_kmer_ID")
    '''

    group_function.keep_group_set(finalGroupMap)
    GroupList = group_function.transfer_group_map_into_list(finalGroupMap)
    objSum = originalG.calculate_objScore_in_groupList(GroupList)
    logging.info ( "The objSocre in groupsList is : %s " % objSum )
    logging.info ( "Choose score in original graph rate: %.6f" % (objSum/totalScore) )
    group_function.write_groupList(GroupList, "patition_phased_kmer_ID")
    t5 = time.time()
    logging.info  ("total cost %s" % (t5-t1) )
