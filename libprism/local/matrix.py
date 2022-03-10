#########################################################################
# File Name: matrix.py
# Author: Yanbo Li
# mail: liyanbotop@gmail.com
# Created Time: Tue 22 Sep 2020 12:49:43 PM AEST
#########################################################################
#!/bin/bash
import os
import sys
#import tools
from libprism.local import tools

def load_all_edges_from_matrix(Filename):

    edges = {}     
    removeNode = set()
    with open(Filename, "r") as f:
        for line in f:
            words = line.split()
            wordsLen = len( words )
            nodeNumber = int( words[0] )
            nodeID = []
            #numList = []
            for i in range(2, wordsLen-1, 2):
                #num = int(words[i+1].split('_')[0]) # support kmer number
                #numList.append(num)    
                nodeID.append(words[i])  
            
            for i in range(nodeNumber - 1):
                for j in range(i+1, nodeNumber):
                    key = (nodeID[i], nodeID[j])
                    Rkey = (nodeID[j], nodeID[i])

                    if key in edges:
                        edges[ key ] = edges[ key ] + 1
                    elif Rkey in edges:
                        edges[ Rkey ] = edges[ Rkey ] + 1
                    elif key not in edges and Rkey not in edges:
                        edges[ key ] =  1
    return edges, removeNode


def load_all_edges_from_matrix2(Filename):
    # only give edges between neighbor nodes

    edges = {}              
    with open(Filename, "r") as f:
        for line in f:
            words = line.split()
            wordsLen = len( words )
            nodeNumber = int( words[0] )
            nodeID = []
            #numList = []
            for i in range(2, wordsLen-1, 2):
                #num = int(words[i+1].split('_')[0]) # support kmer number
                #numList.append(num)    
                nodeID.append(words[i])  
            
            for i in range(nodeNumber - 1):
                j = i+1 
                key = (nodeID[i], nodeID[j])
                Rkey = (nodeID[j], nodeID[i])

                if key in edges:
                    edges[ key ] = edges[ key ] + 1
                elif Rkey in edges:
                    edges[ Rkey ] = edges[ Rkey ] + 1
                elif key not in edges and Rkey not in edges:
                    edges[ key ] =  1
    return edges


def load_edges_from_matrix(Filename): 
    # ID and attribution 

    edges = {}               
    with open(Filename, "r") as f:
        for line in f:
            words = line.split()
            wordsLen = len( words )
            nodeNumber = int( words[0] )
            nodeID = []
            nodeAtt = [] #nodeAttribute
            for i in range(2, wordsLen-1, 2):
                nodeID.append(words[i])
            for i in range(3, wordsLen-1, 2):
                nodeAtt.append(words[i])

            for i in range(nodeNumber - 1):
                for j in range(i+1, nodeNumber):
                    key = (nodeID[i], nodeID[j])
                    if key not in edges:
                        edges[ key ] = []
                    edges[ key ].append( (nodeAtt[i], nodeAtt[j]) ) 
    return edges 

def load_undirected_edges_from_matrix(Filename, k, wType, cov): 
    # ID and attribution 
    # remove self-loop and self-conflict
    #k=31.0 # maximum support
    
    edges = {}              
    nodeSupp = {} # how many reads support this k-mer

    with open(Filename, "r") as f:
        for line in f:
            #print ("debug", line)
            words = line.split()
            wordsLen = len( words )
            #nodeNumber = int( words[0] )
            nodeID = []
            numList = []
            for i in range(2, wordsLen-1, 2):
                num = int(words[i+1].split('_')[0]) # support kmer number
                if num < 2: # high coverge , high hete can use this
                    continue
                numList.append(num)    
                nodeID.append(words[i])  
                if words[i] in nodeSupp:
                    nodeSupp[ words[i] ] += 1
                else:
                    nodeSupp[ words[i] ] = 1
                
            nodeID, numList, removeNodeTemp  = check_self_loop_and_conflict(nodeID, numList)    
            #for now, we only delete conflict node in a read
            #another choice, remove those conflict node in all reads
            update_edges(nodeID, numList, wType, edges)
    print ("edges number %s" % format( len(edges) , ',' ) ) 

    removeNode = set()
    #fout = open("high_reads_suport_nodes", "w")
    #fout.write("nodeID supportNum\n")
    for node in nodeSupp:
        if nodeSupp[node] > 0.5*cov:
            removeNode.add(node)
            #fout.write("%s %s\n" % (node, nodeSupp[node] ))
    #fout.close()        
    print ("remove node whose read coverage higher than hete coverage", 0.5*cov)
    print ("heterozygous node read coverage should smaller than hete coverage (hete coverge=1/2homo)") 
    print ("high reads support node number", len(removeNode))
    #print ("[info] write high reads support nodes in high_reads_suport_nodes")
    return edges, removeNode


def transfer_edges_from_Map2List(edgesMap): 

    edgesList = []
    for key in edgesMap:
        w = edgesMap[key]
        edgesList.append( (key[0], key[1], w) )

    return edgesList


def filter_edges_according_weight(edgesMap, wTemp): # suit for wType =1

    edgesList = []
    totalEdgeNum = len(edgesMap)
    edges = sorted(edgesMap.items(), key = lambda item: item[1])
    cnt = 0
    for (key, w) in edges:
        k1 = tools.get_symmetrical_node(key[0])
        k2 = tools.get_symmetrical_node(key[1])
        if w <= wTemp and (k1, k2) not in edgesMap and (k2, k1) not in edgesMap:
            cnt += 1
            continue
        edgesList.append( (key[0], key[1], w) )

    print ("filter weight small and equal than %s and symmetrical edges not exist" % wTemp)
    print ("filter edges number %s" % format(cnt, ',') )
    print ("weight >%s or (weight <=%s and symmetrical edge exist) edges number %s" % (wTemp, wTemp, format(len(edgesList) , ',' ) ) )
    # wType =1 , use as weightThreshold
    return edgesList


def filter_edges_according_rate(edgesMap, rate): # filter edges without symmetrical edge and weigt <= wTemp

    edgesList = []
    totalEdgeNum = len(edgesMap)
    filterCnt = totalEdgeNum*rate
    edges = sorted(edgesMap.items(), key = lambda item: item[1])
    wTemp, cnt, wMin = 0, 0, 10000

    show_weight_distribution(edgesMap)
    for (key, w) in edges:
        k1 = tools.get_symmetrical_node(key[0])
        k2 = tools.get_symmetrical_node(key[1])
        if cnt < filterCnt and (k1, k2) not in edgesMap and (k2, k1) not in edgesMap:
            cnt += 1
            wTemp = w
            continue
        if w < wMin:
            wMin = w
        edgesList.append( (key[0], key[1], w) )

       
    print ("filter weight small and equal than %s and symmetrical edges not exist" % wTemp)
    print ("filter edges number %s" % cnt)
    print ("weight >%s or (weight <=%s and symmetrical edge exist) edges number %s" % (wTemp, wTemp, len(edgesList) ) )
    # wType =1 , weightThreshold=1
    wTemp = max(wMin, wTemp) # use as weightThreshold
    return edgesList, wTemp

def show_weight_distribution(edgesMap):
 
    weightCnt = {} # weight distribution
    for (k1, k2) in edgesMap:
       w = edgesMap[(k1, k2)]
       if w not in weightCnt:
           weightCnt[w] = 0
       weightCnt[w] += 1
    print ("weight distribution" , sorted(weightCnt.items(), key=lambda item:item[1], reverse=True)[1:10] )
    
    return

def check_self_loop_and_conflict(nodeID, numList):
   
    removeIndexSet =set()
    removeNodeSet = set()
    l = len(nodeID)
    for i in range(l - 1):
        for j in range(i+1, l):
            if nodeID[i] == nodeID[j]:
                removeIndexSet.add(j)
            elif nodeID[i][:-2] == nodeID[j][:-2]: 
                # can change, if coverage not enough, want to include more edges
                removeIndexSet.add(i)
                removeIndexSet.add(j)
                removeNodeSet.add(nodeID[i])
                removeNodeSet.add(nodeID[j])
    newNodeID, newNumList = [], []
    for i in range(l):
        if i not in removeIndexSet:
            newNodeID.append(nodeID[i])
            newNumList.append(numList[i])
    return newNodeID, newNumList, removeNodeSet 


def weight_type(t, newNodeNum):

    edgeNum = newNodeNum*(newNodeNum-1)/2
    if  t == 1:
        return 1
    elif t == 2:
        return round(1/(newNodeNum-1), 5)
    elif t == 3:
        return round(1/(edgeNum), 5)


def update_edges(nodeID, numList, wType, edges):
    
    newNodeNum = len(nodeID)
    for i in range(newNodeNum - 1):
        for j in range(i+1, newNodeNum):
            assert nodeID[i] != nodeID[j]
            assert nodeID[i][:-2] != nodeID[j][:-2]
            key = (nodeID[i], nodeID[j])
            Rkey = (nodeID[j], nodeID[i])
            if key in edges:
                edges[ key ] = round(edges[ key ] + 
                       # (numList[i]*numList[j]/(k-1)/(k-1))*
                       weight_type(wType, newNodeNum), 5) # +1/2
            elif Rkey in edges:
                edges[ Rkey ] = round(edges[ Rkey ] + 
                        #(numList[i]*numList[j]/(k-1)/(k-1))*
                        weight_type(wType, newNodeNum), 5) #+ 1/2
            elif key not in edges and Rkey not in edges:
                edges[ key ] = round(weight_type(wType, newNodeNum)
                        #*  (numList[i]*numList[j]/(k-1)/(k-1))
                        ,5)
    return 
