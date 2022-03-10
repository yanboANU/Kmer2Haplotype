#########################################################################
#File Name: build_graph.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 26 Nov 2019 20:22:50 AEDT
# according to Yulin, no trick and add an optmal function
# now, only add a fuction update_weight_for_conflict_edges
#########################################################################
#!/bin/bash

import os
import sys
import logging
import time
import networkx as nx
#import tools
from networkx.drawing.nx_pydot import write_dot
import logging
import pdb

from libprism.local import tools
from libprism.local import group_function
from libprism.local import matrix
from libprism.local import myGraph
#from libprism.local import train


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

        if n2 < n1:
            n1, n2, symN1, symN2 = n2, n1, symN2, symN1
        if symN1 < n1:
            n1, symN1 = symN1, n1
            n2, symN2 = symN2, n2
        
            
        weightEdges[wKey].append( (n1, n2, symN1, symN2) ) # debug should not sort
            
    comCnt = nx.number_connected_components(g) #         
    if comCnt == 2:
        #graphs = list( nx.connected_component_subgraphs(g) ) 
        graphs = g.get_subgraphs()
        group1 = set( graphs[0].nodes )
        group2 = set( graphs[1].nodes )
        objScore = g.size(weight='weight') # weight sum in g
        if len(group1) == len(group2):
            #print (group1)
            #print (group2) 
            print ("it's an symmetrical connected component and they phasing natrually")
            return True, group1, group2, objScore
    elif comCnt == 3:   
        #graphs = sorted (list( nx.connected_component_subgraphs(g) ) , key=lambda x:len( set(x.nodes)  ), reverse=True )
        graphs = sorted ( g.get_subgraphs() , key=lambda x:len( set(x.nodes)  ), reverse=True )
        l1,l2,l3 = len(graphs[0].nodes), len(graphs[1].nodes), len(graphs[2].nodes)
        print ("debug: three component", l1, l2, l3)
        if l1!=l2:
            group1 = set( graphs[0].nodes )
            group2 = set( graphs[1].nodes ).union( set(graphs[2].nodes ) )
            objScore = g.size(weight='weight') # weight sum in g
            if len(group1) == len(group2):
                #print (group1)
                #print (group2) 
                print ("three connected component and they phasing natrually")
                return True, group1, group2, objScore
    elif comCnt == 4:  
        #graphs = sorted (list( nx.connected_component_subgraphs(g) ) , key=lambda x:len( set(x.nodes)  ), reverse=True )

        graphs = sorted ( g.get_subgraphs() , key=lambda x:len( set(x.nodes)  ), reverse=True )
        l1,l2,l3,l4 = len(graphs[0].nodes), len(graphs[1].nodes), len(graphs[2].nodes), len(graphs[3].nodes)
        print ("debug: four component", l1, l2, l3, l4)
        if l1!=l2 and l3!=l4:
            group1 = set( graphs[0].nodes ).union( set(graphs[3].nodes ) )
            group2 = set( graphs[1].nodes ).union( set(graphs[2].nodes ) )
            objScore = g.size(weight='weight') # weight sum in g
            if len(group1) == len(group2):
                #print (group1)
                #print (group2) 
                print ("four connected component and they phasing natrually")
                return True, group1, group2, objScore
    return weightEdges  


def add_pair_by_pair(group, graph, sortedWeight, weightEdges, weightThreshold, conflictSet):

    #group = group_function.inital_group(graph)

    dealedEdge, dealedRate = 0, 0.1
    nodeNum = graph.number_of_nodes()
    edgeNum = graph.number_of_nodes()
    #debug_unuse_weight(sortedWeight, weightThreshold)
    w = -1
    finish = False
    for w in sortedWeight:
        if w<=weightThreshold: 
            logging.debug ("weight smaller than threshold stop merge, threshold is: %s" % weightThreshold)
            break
        #print ("dealed edges number:", dealedEdge)
        
        if finish==True or (dealedEdge*2.0/edgeNum) > 5*dealedRate:
            print ("dealed %.5f edges, w is %s" % ( (dealedEdge*2.0/edgeNum), w))
            break
        elif (dealedEdge*2.0/edgeNum) > dealedRate:
            finish =  add_pair_2_nodes(weightEdges[w], graph, 
                    group, conflictSet, nodeNum)
        else:
            finish =  add_pair_at_same_weight(weightEdges[w], graph, 
                    group, conflictSet, nodeNum) 
        dealedEdge += len( weightEdges[w] )

    return group, w



def merge_group_strict(group, phasedNodes, graph, groupSum, conflictSet, l, w):
    #merge group and phase unphased node

    #group_function.keep_group_set2(group, phasedNodes)
    sLabel = False
    lenG = len(group)

    
    t1 = time.time()
    cnt = 4
    label = True
    while (cnt >= 1):
        while (label==True and len(group) >2):
            print ("merge node group strict", cnt*0.5)
            label = group_function.merge_node_group_strict(graph, 
                    group, phasedNodes, groupSum, conflictSet, cnt*0.5) 
        cnt = cnt -1
        label =True
    t2 = time.time()    
    logging.debug("Debug: merge node graph strict cost %s" % (t2-t1) )
 
    cnt = 4
    t2 = time.time()    
    label = True
    while (cnt >= 1):
        while (label==True and len(group) >2):
            print ("merge_group_according_threshold_strict", cnt*0.5)
            label = group_function.merge_group_according_threshold_strict(graph, 
                    group, phasedNodes, groupSum, conflictSet, cnt*0.5, l)
        cnt = cnt - 1    
        label =True

    t3 = time.time()    
    logging.debug("Debug: merge strict cost %s" % (t3-t2) ) 
    
    #group_function.check_group(group) 
    if len(group) < lenG:
        sLabel = True

    return sLabel

def merge_group(group, phasedNodes, graph, conflictSet, l, w=1):
    
    groupSum= {}
    sLabel, lLabel = True, True
    cnt = 1
    while sLabel==True or lLabel==True:
        print ("round ", cnt)
        sLabel = merge_group_strict(group, phasedNodes, graph, groupSum, 
                conflictSet, l, w) 
        
        lLabel = merge_group_loose(group, phasedNodes, graph, groupSum, 
                conflictSet, l) 
        cnt += 1 
    
    return 

def merge_group_loose(group, phasedNodes, graph, groupSum, conflictSet, l):
    #merge group and phase unphased node
    
    group_function.keep_group_set2(group, phasedNodes) 
    lLabel = False
    lenG = len(group)
    
    t1 = time.time()
    cnt = 5 
    label = True
    while (cnt >= 1):
        while (label==True and len(group) >2):
            print ("merge_group_according_threshold", 0.4*cnt)
            label = group_function.merge_group_according_threshold(graph, 
                    group, phasedNodes, 0.4*cnt, groupSum, conflictSet, l)
            #print ("merge two groups rate ", 0.2*cnt)
            #label = group_function.merge_two_groups_rate(graph, group, phasedNodes, 0.2*cnt, groupSum, conflictSet)
        cnt = cnt - 1  
        label =True

    t2 = time.time()
    logging.debug("Debug: merge threshold cost %s" % (t2-t1) )

    t3 = time.time() 
    while len(group) >2:
        print ("merge only choice  ")
        label= group_function.merge_groups_only_choice(graph, group, phasedNodes, groupSum, conflictSet)
        if label == False:
            break
    t4 = time.time()
    logging.debug("Debug: merge only cost %s" % (t4-t3) )
    
    if len(group) < lenG:
        lLabel = True
    return lLabel 


def build_final_graph(graph, group, groupSum, unphasedNode=[]):
   

    # add unphasedNode in group
    for n in unphasedNode:
        #assert n not in group
        group[n] = set()
        group[n].add(n)
     
    # debug 
    group_function.keep_group_set(group)
    
    groupKeys = sorted(group.keys()) # deep copy
    groupLen = len(group)
    #print ("before merge group: %s pairs of group" %  (groupLen/2) )


    print ("final group keys:", groupKeys)
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
                if w1+w2+w3+w4==0:
                    continue
                #trainningList.append( (n1, n2, w1, w2, w3, w4, l1, l2) ) 
                if w1 > 0: 
                    edgesList.append((n1, n2, w1))
                if w2 > 0: 
                    edgesList.append((n1, symN2, w2))
                if w3 > 0: 
                    edgesList.append((symN1, n2, w3))
                if w4 > 0: 
                    edgesList.append((symN1, symN2, w4))
                 
    # write trainning set
    '''
    fout = open("trainning_set", "w")
    for (n1, n2, w1, w2, w3, w4, l1, l2) in trainningList:
        fout.write("%s %s %s %s %s %s %s %s\n" % (n1, n2, w1, w2, w3, w4, l1, l2))
        print("%s %s %s %s %s %s %s %s\n" % (n1, n2, w1, w2, w3, w4, l1, l2))
    fout.close()  
    '''
    finalG = myGraph.myGraph()
    finalG.add_weighted_edges_from(edgesList)
    totalScore, num = finalG.displayGraph()
    print ("analyse final graph")
    for (n1, n2, w) in edgesList:
        print (n1, n2, w, len(group[n1]), len(group[n2]))
    finalG.add_type_A_edges()
    print ("After add A edges, number of components: ", nx.number_connected_components(finalG))
    write_dot(finalG, "finalGraph" + ".dot")
    
    #final_step(finalG, groupKeys)
    finalG1node, finalG2node=set(), set()
    return finalG1node, finalG2node

 
def phasing_in_connected_component(graph, weightThreshold, l):
    
    phasedNodes={}
    group = group_function.inital_group(graph)

    nodeNum = graph.number_of_nodes() 
    if nodeNum <= 3:
        return group, phasedNodes
    t1 = time.time()

    #logging.debug ("remove/update conflict edges")   
    graph.update_weight_for_conflict_edges()
    res = check_symmetrical_connected_component(graph) 
    #logging.debug("after remove 'A' edges")
    graph.displayGraph()
    t2 = time.time()
    #logging.debug("Debug: check symmetrical_connected_component cost %s" % (t2-t1) )
    
    if len(res) == 4 and not isinstance(res, dict) and res[0] == True:
        k1, k2 = min(res[1]), min(res[2])
        group[ k1  ] = res[1]
        group[ k2  ] = res[2]
        for ele in res[1]:
            if ele != k1:
                group[ele] = k1
        for ele in res[2]:
            if ele != k2:
                group[ele] = k2
        #assert k1[:-2] ==  k2[:-2] 
        objScore = res[3]
        return group, phasedNodes 
    
    weightEdges = res
    sortedWeight = sorted(weightEdges.keys(), reverse=True)

    t1 = time.time()
    conflictSet = set()
    w = add_pair_by_pair(group, graph, sortedWeight, weightEdges, weightThreshold, conflictSet)
    t2 = time.time()
    logging.debug("Debug: add pair cost %s" % (t2-t1) )
    print (len(conflictSet))
    print("before add pair, group size, phased size, sum", len(group), len(phasedNodes), len(group) + len(phasedNodes))
    group_function.keep_group_set2(group, phasedNodes)
    print("after add pair, group size, phased size, sum ", len(group), len(phasedNodes), len(group) + len(phasedNodes)) 
    merge_group(group, phasedNodes, graph, conflictSet, l, w)

    return group, phasedNodes

def add_pair_at_same_weight(weightEdges, graph, group, conflictSet, nodeNum):

    
    wOrder = {}
    finish = False
    for eleList in weightEdges:
    #for eleList in es:
        n1, n2, symN1, symN2 = eleList[0], eleList[1], eleList[2], eleList[3]
        # way 2 check ../local-back

        if n1[:-2] == n2[:-2]:
            logging.debug ("%s %s cannot merge" % (n1, n2) )
            continue

        if group_function.is_two_nodes_merged(n1, n2, symN1, symN2, group):
            continue   
        w1 = graph.get_weight(n1, n2)
        symW = graph.get_weight(symN1, symN2)

        if isinstance(w1, int) == False or  isinstance(symW, int) == False:
            continue
        
        if w1 <2 or symW < 2: # debug1 k=31, 2 is ok
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
        if w < 4:
            continue
        for (n1, n2, symN1, symN2) in fourList: # = eleList[0], eleList[1], eleList[2], eleList[3]
            #assert n1[:-2] != n2[:-2]
            # debug
            #g1, g2, r1, r2 = group_function.get_group_set(group, n1, n2) # can remove later  
            label = group_function.merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
            if label == True:

                #debug    
                '''
                w1 = graph.get_weight(n1, n2)
                w4 = graph.get_weight(symN1, symN2)
                w2 = graph.get_weight(n1, symN2)
                w3 = graph.get_weight(symN1, n2)
                print (n1, n2, w1, w2, w3, w4, r1, r2, len(g1), len(g2))
                '''  
                g1, g2, r1, r2 = group_function.get_group_set(group, n1, symN1)
                phasedNum = len( g1 ) + len( g2 ) 
                if phasedNum >= 0.9999*nodeNum:
                    finish = True
                    logging.info ("phasing %s nodes in all %s nodes, break" % (phasedNum, nodeNum) )
                    break
    return finish #, conScore

def add_pair_2_nodes(weightEdges, graph, group, conflictSet, nodeNum):
 
    wOrder = {}
    finish = False
    for eleList in weightEdges:
        n1, n2, symN1, symN2 = eleList[0], eleList[1], eleList[2], eleList[3]

        if n1[:-2] == n2[:-2]:
            logging.debug ("%s %s cannot merge" % (n1, n2) )
            continue

        if group_function.is_two_nodes_merged(n1, n2, symN1, symN2, group):
            continue   
        w1 = graph.get_weight(n1, n2)
        symW = graph.get_weight(symN1, symN2)

        if isinstance(w1, int) == False or  isinstance(symW, int) == False:
            continue
        
        if w1 <2 or symW < 2: # debug1 k=31, 2 is ok
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
        if w < 4:
            continue
        for (n1, n2, symN1, symN2) in fourList: 
            g1, g2, r1, r2 = group_function.get_group_set(group, n1, n2) 
            if len(g1) + len(g2) > 2:
                continue
            label = group_function.merge_two_groups(n1, n2, symN1, symN2, group, conflictSet)
            if label == True:
                g1, g2, r1, r2 = group_function.get_group_set(group, n1, symN1)
                phasedNum = len( g1 ) + len( g2 ) 
                if phasedNum >= 0.9999*nodeNum:
                    finish = True
                    logging.info ("phasing %s nodes in all %s nodes, break" % (phasedNum, nodeNum) )
                    break
    return finish #, conScore





def build_original_graph(matrixFile, mType, k, cov):
    
    t1 = time.time()

    print ("mtype", mType, "matrix file", matrixFile)
    if mType == 0:
        edges, removeNode = matrix.load_all_edges_from_matrix(matrixFile)
    elif mType == 1:
        edges, removeNode = matrix.load_undirected_edges_from_matrix(matrixFile, k, 1, cov)
    else:
        print ("mtype error", mType)
        sys.exit()
    t2 = time.time() 
    logging.info ("load matrix cost %s" % (t2-t1) )
    edgesList = matrix.transfer_edges_from_Map2List(edges)
    originalG = myGraph.myGraph()
    originalG.add_weighted_edges_from(edgesList)
    totalScore, num = originalG.displayGraph()
    if len(removeNode) >0:
        originalG.remove_nodes(removeNode)
        totalScore, num = originalG.displayGraph()
    originalG.add_type_A_edges()
    num = nx.number_connected_components(originalG) 
    print ("After add A edges, number of components: ", num)
    
    return originalG, totalScore, num


def write_phasing_results(finalGroupMap, originalG, totalScore):

    #group_function.keep_group_set2(finalGroupMap)
    print("final group size:", len(finalGroupMap))
    GroupList = group_function.transfer_group_map_into_list(finalGroupMap)
    objSum = originalG.calculate_objScore_in_groupList(GroupList)
    logging.info ( "The objSocre in groupsList is : %s " % objSum )
    logging.info ( "Choose score in original graph rate: %.6f" % (objSum/totalScore) )
    group_function.write_groupList(GroupList, "patition_phased_kmer_ID") 
    group_function.keep_group_set(finalGroupMap)
    phasedNodes = group_function.get_phased_nodes(finalGroupMap) 
    unphasedNode = set(originalG.nodes() ) -  phasedNodes
    print( "[debug] in whole graph, phased Node number: ",  len(phasedNodes))
    print( "[debug] in whole graph, unphased Node number: ",  len(unphasedNode))
    print( "[debug] all node number in original graph: ", originalG.number_of_nodes()   )
 
    print("final group size (group size >=10):", len(finalGroupMap))
    GroupList = group_function.transfer_group_map_into_list(finalGroupMap)
    group_function.write_groupList(GroupList, "phased_kmer_ID")
    return

def phasing_components(graphs, wTemp, l):
    
    iComponent = 0
    groupNum  = 0
    f_groupMap, f_phasedNodes = {}, {}
    #cfg = train.trainning() 
    for g in graphs:
        gNodesNum = g.number_of_nodes()
        if ( gNodesNum <=3 ):
            continue
        iComponent += 1  
        #logging.debug("------------------------------------------------------------------")
        print("------------------------------------------------------------------")
        #logging.debug("phasing connected component %s node number: %s" % (iComponent, gNodesNum ) )
        print("phasing connected component %s node number: %s" % (iComponent, gNodesNum ) )
        #if gNodesNum < 100: 
        #    write_dot(g, "componemt" + str(iComponent) + ".dot")

        groupMap, phasedNodes = phasing_in_connected_component(g, wTemp, l)
        f_groupMap.update(groupMap)
        f_phasedNodes.update(phasedNodes)
        #group_function.remove_small_group(groupMap) #
        
        groupsLen = len(groupMap)
        print("group number %s" % ( groupsLen))
    return f_groupMap, f_phasedNodes 

def merge_on_reliable_graph(originalG, finalGroupMap, phasedNodes, k, cov, l):

    originalG.remove_type_A_edges()
    print ("final graph")
    originalG.remove_type_A_edges()
    originalG.simple_show()

    groupSum = {} #remove edges, groupSum need recalulate
    conflictset = set()
    
    print ("debug1", len(finalGroupMap), len(phasedNodes), len(finalGroupMap)+len(phasedNodes) )
    t3 = time.time() 
    while len(finalGroupMap) >2:
        print ("merge only choice 1  ")
        label= group_function.merge_groups_only_choice(originalG, finalGroupMap, phasedNodes, groupSum, conflictset)
        if label == False:
            break
    t4 = time.time()
    logging.debug("Debug: merge only cost %s" % (t4-t3) )
     
    print ("debug2", len(finalGroupMap), len(phasedNodes), len(finalGroupMap)+len(phasedNodes) )
    label = True
    if l == True:
       cnt = 2 
       rate=0.01
    else:
       cnt = 5
       rate = 0.001
    while (cnt >= 1):
        while label and len(finalGroupMap) >2:
            print ("merge on reliable", rate*cnt)
            label = group_function.merge_group_on_reliable_graph(originalG, 
                    finalGroupMap, phasedNodes, groupSum, conflictset, cnt*rate, l)
        cnt = cnt - 1  
        label =True
   
    print ("debug3", len(finalGroupMap), len(phasedNodes), len(finalGroupMap)+len(phasedNodes) )

    t3 = time.time() 
    while len(finalGroupMap) >2:
        print ("merge only choice 2 ")
        label= group_function.merge_groups_only_choice(originalG, finalGroupMap, phasedNodes, groupSum, conflictset)
        if label == False:
            break
    t4 = time.time()
    logging.debug("Debug: merge only cost %s" % (t4-t3) )

    vipPair, temp = group_function.show_group_connection(finalGroupMap, originalG, groupSum)
    


    return vipPair 


def merge_on_original_graph(matrixFile, finalGroupMap, phasedNodes, k, cov, vipPair, l):

    originalG, totalScore, num = build_original_graph(matrixFile, 0, k, cov)
    print ("original graph")
    originalG.remove_type_A_edges()
    originalG.simple_show()
    print ("on original graph")
    groupSum = {}
    conflictset = set()
    
    group_function.merge_vip_pair(originalG, finalGroupMap, 
            phasedNodes, groupSum, conflictset, vipPair)
    '''
    if l == False:
        return

    label = True
    cnt = 2
    while (cnt >= 1):
        while label and len(finalGroupMap) >2:
            print ("merge on original", 0.01*cnt)
            label = group_function.merge_group_on_original_graph(originalG, 
                    finalGroupMap, phasedNodes, groupSum, conflictset, cnt*0.01) #, groupDegree)
        cnt = cnt - 1  
        label =True
     
    #print ("debug3", len(finalGroupMap), len(phasedNodes), len(finalGroupMap)+len(phasedNodes) )

    temp, groupDegree = group_function.show_group_connection(finalGroupMap, originalG, groupSum)
    
    t3 = time.time() 
    while len(finalGroupMap) >2:
        print ("merge only choice 2 ")
        label= group_function.merge_groups_only_choice(originalG, finalGroupMap, phasedNodes, groupSum, conflictset)
        if label == False:
            break
    t4 = time.time()
    logging.debug("Debug: merge only cost %s" % (t4-t3) )
    '''


    return 

def phasing_kmer(mType, k, cov, l=True):

    t1 = time.time()
    matrixFile = "all_matrix"
    originalG, totalScore, num = build_original_graph(matrixFile, mType, k, cov)
    
    #try
    #print ("debug l, k", l, k)
    #assert l==False
    #assert int(k)>=29
    if l == False and int(k)>=29:
        print ("for simulate data default a=1")
        a=1
    else:
        print ("for simulate data (small k) default a=2")
        a=2
    originalG.remove_edges_according_weight(a) # real data:3 ;  simulate data:2
    print( "[debug] remove edges (w<= %s and symW==0) " % (a) )
    
    th1, th2 = originalG.get_node_degree_distribution()
    originalG.remove_high_degree_node(th1)  # cov/2 is not right, we don't know node degree
    
    originalG.simple_show()

    wTemp = 1.0
    if num > 1:
        graphs = list( (originalG.subgraph(c)).copy() for c in nx.connected_components(originalG) )
    else:
        graphs = [ originalG.copy() ]
    logging.info("number of components: %s " % ( len(graphs) ) )

    groupMap, phasedNodes = phasing_components(graphs, wTemp, l)
    
    vipPair = merge_on_reliable_graph(originalG, groupMap, phasedNodes, k, cov, l)
    merge_on_original_graph(matrixFile, groupMap, phasedNodes, k, cov, vipPair, l)

    print ("before remove small group", len(groupMap))
    group_function.remove_small_group(groupMap)
    print ("after remove small group", len(groupMap))

    write_phasing_results(groupMap, originalG, totalScore)

    t5 = time.time()
    logging.info  ("total cost %s" % (t5-t1) )
    return

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
#logging.basicConfig(stream=sys.stderr, level=logging.INFO)

trainningList = []
if __name__ == '__main__':

    if len(sys.argv) < 6:
        print ("matrix_file kmer_size filter_edge_rate(0) homo_cov mType(0)")
        sys.exit()

    matrixFile = sys.argv[1]
    k =int(sys.argv[2])
    wRate = float(sys.argv[3])
    cov = int(sys.argv[4])
    mType=int(sys.argv[5]) # matrix Type
    wType = 1
   
    phasing_kmer(mType, k, cov)
    '''
    originalG, totalScore, num = build_original_graph(matrixFile, mType, k, cov)

    wTemp = 1.0
    if num > 1:
        graphs = list( (originalG.subgraph(c)).copy() for c in nx.connected_components(originalG) )
    else:
        graphs = [ originalG.copy() ]
    logging.info("number of components: %s " % ( len(graphs) ) )
 
    finalGroupMap =  phasing_components(graphs)
    originalG.remove_type_A_edges()

    write_phasing_results(finalGroupMap)
    '''
